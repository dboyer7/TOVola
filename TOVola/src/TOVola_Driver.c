#include "TOVola_ODE.h" //ODE and associated functions
#include "GRHayLib.h" //Access to GRHayL library in the ET
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_odeiv2.h>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

/********************************
// TOVola is a TOV solver inspired by NRpy designed specifically for use in the Einstein Toolkit.
// It integrates the TOV equations for a solution to static, spherically symmetric stars (Neutron stars in specific), and interpollates the solution the the ET grid.
//
// MOST IMPORTANTLY, it handles multiple types of EOS in the solution: Simple Polytrope, Piecewise Polytrope, and Tabulated EOS. An upgrade from the original ET TOVSolver thorn.
//
// Simple Polytrope:
//	-- The Simple Polytrope is a special version of the Piecewise Polytrope with one EOS region relating Pressure to baryon density (P=K*rho_baryon^Gamma), where K and Gamma are parameters to be input
//      -- Only one K and Gamma are necessary to run this implementation, and a Gamma_Thermal defaulted to 2.0
//      -- Hamiltonian Constraint Violation Validated with Baikal
// Piecewise Polytrope:
//	-- The more general verion of the simple polytrope where you can choose the number of regions for your polytrope, each region denoted by (P=K_i*rho_baryon^Gamma_i)
//      -- You can get your necessary parameters by following the method outlined in Read et. al. (2008)
//      -- Since we follow Read's method, you only need to input one K value (The GRHayL library will handle the rest), and all the Gammas and rho_baryon boundary points for the number of specified regions.
//      -- Hamiltonian Constraint Violation Validated with Baikal
// Tabulated EOS:
//	-- The tabulated solver reads in a EOS table from your computer.
//      -- Once located, GRHayL slices the table for beta_equilibirum and then uses an interpollator to find the necessary densities through the calculations.
//      -- Hamiltonian Constraint Violation Validated with Baikal
//
// Users who wish to interface with this program need only create a parfile for a simulation. Output variables will be adjusted for use in ADMbase and Hydrobase as initial data.
// IllinoisGRMHD is used to set Tmunu from the data and Baikal tests the constraint Violation.
// 
// David Boyer (8/23/24)
********************************/

//Global pointers to save raw data. (Deal with it.)
CCTK_REAL *TOVola_Raw_rSchw = NULL;
CCTK_REAL *TOVola_Raw_rho_energy = NULL;
CCTK_REAL *TOVola_Raw_rho_baryon = NULL;
CCTK_REAL *TOVola_Raw_P = NULL;
CCTK_REAL *TOVola_Raw_M = NULL;
CCTK_REAL *TOVola_Raw_nu = NULL;
CCTK_REAL *TOVola_Raw_Iso_r = NULL;
int TOVola_Numpoints = 0;

//First Integration: Count how much memory needed to allocate for solution.
void TOVola_TOV_Integrator_1(CCTK_ARGUMENTS){

    //Declare value from ET
    DECLARE_CCTK_PARAMETERS
    CCTK_INFO("First Integration: Analyzing how much memory to allocate");
    CCTK_INFO("Beginning ODE Solver using GSL for TOVola...");

    //Track current position
    double TOVola_current_position = 0.0;

    //Checking and setting EOS
    if(CCTK_EQUALS("Simple",TOVola_EOS_type)){
      CCTK_INFO("Simple Polytrope");
      if(ghl_eos->neos!=1){
        CCTK_INFO("Error: Too many regions for the simple polytrope.");
        CCTK_INFO("Check your value for neos, or use a piecewise polytrope");
        CCTK_ERROR("Shutting down due to error...\n");}
        }
    else if(CCTK_EQUALS("Piecewise",TOVola_EOS_type)){
      CCTK_INFO("Piecewise Polytrope");}
    else if(CCTK_EQUALS("Tabulated",TOVola_EOS_type)){
      CCTK_INFO("Tabulated EOS");
      ghl_tabulated_compute_Ye_P_eps_of_rho_beq_constant_T(TOVola_Tin, ghl_eos);
    }
    else{
      CCTK_INFO("ERROR: Invalid EOS type. Must be either 'Simple', 'Piecewise', or 'Tabulated'");
      CCTK_ERROR("Shutting down due to error...\n");}
    
    //Create the GSL struct for the TOV ODE system
    int TOVola_number_of_equations = 4; //P,M,nu,rbar
    int TOVola_number_of_constants = 2; //rho_energy, rho_baryon
    struct TOVola_parameters TOVola_theParams = { 
      .dimension = TOVola_number_of_constants,
    };
    TOVola_theParams.rhoCentral_baryon = TOVola_central_baryon_density;
    
    //Declare the GSL ODE system and Driver.
    gsl_odeiv2_system TOVola_system = {TOVola_ODE,TOVola_jacobian_placeholder,TOVola_number_of_equations,&TOVola_theParams};
    gsl_odeiv2_driver * TOVola_d;
    //Set the ODE type (RK4(5) or DP7(8))
    if (CCTK_EQUALS("ARKF",TOVola_ODE_method)) {
       TOVola_d = gsl_odeiv2_driver_alloc_y_new(&TOVola_system, gsl_odeiv2_step_rkf45, TOVola_step, TOVola_error_limit, TOVola_error_limit);}
    else if (CCTK_EQUALS("ADP8",TOVola_ODE_method)) {
       TOVola_d = gsl_odeiv2_driver_alloc_y_new(&TOVola_system, gsl_odeiv2_step_rk8pd, TOVola_step, TOVola_error_limit, TOVola_error_limit);}
    else{
      CCTK_INFO("Invalid Step type. May only use ARKF or ADP8.");
      CCTK_ERROR("Shutting down due to error.\n");
      } 
    
    //Set maximum and minimum step size in adaptive method.
    gsl_odeiv2_driver_set_hmax(TOVola_d, TOVola_absolute_max_step);
    gsl_odeiv2_driver_set_hmin(TOVola_d, TOVola_absolute_min_step);

    //Declare ODE variables and begin SOLVING!
    double TOVola_y[TOVola_number_of_equations];
    double TOVola_c[TOVola_number_of_constants];
    
    TOVola_get_initial_condition(TOVola_y,&TOVola_theParams); //Get initial conditions
    TOVola_assign_constants(TOVola_c,&TOVola_theParams); //Set initial rho_baryon and rho_energy
    
    TOVola_Numpoints++; //First step in the solution: For first integration

    //Until the surface of the star
    for (int i = 0; i < TOVola_size; i++){
        
        TOVola_exception_handler(TOVola_current_position,TOVola_y); //Check to make sure pressure didn't go negative near the surface (numerical error check)
        gsl_odeiv2_driver_apply (TOVola_d, &TOVola_current_position, TOVola_current_position+TOVola_step, TOVola_y); //Step through the integration
        TOVola_exception_handler(TOVola_current_position,TOVola_y); //Just another check to be sure. Can't be too safe.
        TOVola_evaluate_rho_and_eps(TOVola_current_position,TOVola_y,&TOVola_theParams); //Evaluate rho_baryon and rho_energy at the new location.
        TOVola_assign_constants(TOVola_c,&TOVola_theParams); //Assign rho_baryon and rho_energy.
    	
    	TOVola_Numpoints++; //Count how much memory is needed for allocation, so we dont waste memory.

        //Are we at the surface?
        if (TOVola_do_we_terminate(TOVola_current_position, TOVola_y, &TOVola_theParams) == 1) {
            i = TOVola_size-1;             
        } 
        if (i == TOVola_size-1) { 
            CCTK_INFO("Finished Integration!"); //if so, the integration is complete, and we can exit out.
            CCTK_VINFO("FINAL: Mass :,\t%15.14e,\t",TOVola_y[2]); //See what the final mass is.
            }
    }

    gsl_odeiv2_driver_free(TOVola_d); //free the driver
    CCTK_INFO("ODE Solver using GSL for TOVola Shutting Down...\n"); //Shut down the first integration
    return;
    
}
  
//Second Integration: Save Solution to Memory
void TOVola_TOV_Integrator_2(CCTK_ARGUMENTS){

    //Declare value from ET
    DECLARE_CCTK_PARAMETERS
    CCTK_INFO("Second Integration: Writing solution to memory");
    CCTK_INFO("Beginning ODE Solver using GSL for TOVola...");

    int TOVola_this_point = 0; //Need this to track where we are in the solution. (for writing purposes)
    double TOVola_current_position = 0.0; //Track current position

    //Checking and setting EOS
    if(CCTK_EQUALS("Simple",TOVola_EOS_type)){
      CCTK_INFO("Simple Polytrope");
      if(ghl_eos->neos!=1){
        CCTK_INFO("Error: Too many regions for the simple polytrope.");
        CCTK_INFO("Check your value for neos, or use a piecewise polytrope");
        CCTK_ERROR("Shutting down due to error...\n");}
        }
    else if(CCTK_EQUALS("Piecewise",TOVola_EOS_type)){
      CCTK_INFO("Piecewise Polytrope");}
    else if(CCTK_EQUALS("Tabulated",TOVola_EOS_type)){
      CCTK_INFO("Tabulated EOS");
      ghl_tabulated_compute_Ye_P_eps_of_rho_beq_constant_T(TOVola_Tin, ghl_eos);
    }
    else{
      CCTK_INFO("ERROR: Invalid EOS type. Must be either 'Simple', 'Piecewise', or 'Tabulated'");
      CCTK_ERROR("Shutting down due to error...\n");}
    
    //Create the GSL struct for the TOV ODE system
    int TOVola_number_of_equations = 4; //P,M,nu,rbar
    int TOVola_number_of_constants = 2; //rho_energy, rho_baryon
    struct TOVola_parameters TOVola_theParams = { 
      .dimension = TOVola_number_of_constants,
    };
    TOVola_theParams.rhoCentral_baryon = TOVola_central_baryon_density;
    
    gsl_odeiv2_system TOVola_system = {TOVola_ODE,TOVola_jacobian_placeholder,TOVola_number_of_equations,&TOVola_theParams};
    gsl_odeiv2_driver * TOVola_d;
    //Set the ODE type (RK4(5) or DP7(8))
    if (CCTK_EQUALS("ARKF",TOVola_ODE_method)) {
       TOVola_d = gsl_odeiv2_driver_alloc_y_new(&TOVola_system, gsl_odeiv2_step_rkf45, TOVola_step, TOVola_error_limit, TOVola_error_limit);}
    else if (CCTK_EQUALS("ADP8",TOVola_ODE_method)) {
       TOVola_d = gsl_odeiv2_driver_alloc_y_new(&TOVola_system, gsl_odeiv2_step_rk8pd, TOVola_step, TOVola_error_limit, TOVola_error_limit);}
    else{
      CCTK_INFO("Invalid Step type. May only use ARKF or ADP8.");
      CCTK_ERROR("Shutting down due to error.\n");
      } 

    //Set maximum and minimum step size in adaptive method.
    gsl_odeiv2_driver_set_hmax(TOVola_d, TOVola_absolute_max_step);
    gsl_odeiv2_driver_set_hmin(TOVola_d, TOVola_absolute_min_step);

    //Declare ODE variables and begin SOLVING!
    double TOVola_y[TOVola_number_of_equations];
    double TOVola_c[TOVola_number_of_constants];
    
    TOVola_get_initial_condition(TOVola_y,&TOVola_theParams); //Get initial conditions
    TOVola_assign_constants(TOVola_c,&TOVola_theParams); //Set initial rho_baryon and rho_energy

    //For the Second Integration, we are actually writing the solution to memory
    TOVola_Raw_rSchw[TOVola_this_point] = TOVola_current_position;
    TOVola_Raw_rho_energy[TOVola_this_point] = TOVola_c[0];
    TOVola_Raw_rho_baryon[TOVola_this_point] = TOVola_c[1];
    TOVola_Raw_P[TOVola_this_point] = TOVola_y[0];
    TOVola_Raw_M[TOVola_this_point] = TOVola_y[2];
    TOVola_Raw_nu[TOVola_this_point] = TOVola_y[1];
    TOVola_Raw_Iso_r[TOVola_this_point] = TOVola_y[3];
    TOVola_this_point++;

    for (int i = 0; i < TOVola_size; i++){
	
	TOVola_exception_handler(TOVola_current_position,TOVola_y); //Check to make sure pressure didn't go negative near the surface (numerical error check)
        gsl_odeiv2_driver_apply (TOVola_d, &TOVola_current_position, TOVola_current_position+TOVola_step, TOVola_y); //Step through the integration

        TOVola_exception_handler(TOVola_current_position,TOVola_y); //Just another check to be sure. Can't be too safe.
        TOVola_evaluate_rho_and_eps(TOVola_current_position,TOVola_y,&TOVola_theParams); //Evaluate rho_baryon and rho_energy at the new location.
        TOVola_assign_constants(TOVola_c,&TOVola_theParams); //Assign rho_baryon and rho_energy.

        //Save the raw data to memory
	TOVola_Raw_rSchw[TOVola_this_point] = TOVola_current_position;
        TOVola_Raw_rho_energy[TOVola_this_point] = TOVola_c[0];
        TOVola_Raw_rho_baryon[TOVola_this_point] = TOVola_c[1];
        TOVola_Raw_P[TOVola_this_point] = TOVola_y[0];
        TOVola_Raw_M[TOVola_this_point] = TOVola_y[2];
        TOVola_Raw_nu[TOVola_this_point] = TOVola_y[1];
        TOVola_Raw_Iso_r[TOVola_this_point] = TOVola_y[3];
        TOVola_this_point++;

        //Are we at the surface?
        if (TOVola_do_we_terminate(TOVola_current_position, TOVola_y, &TOVola_theParams) == 1) {
            i = TOVola_size-1;             
        } 
        if (i == TOVola_size-1) { 
            CCTK_INFO("Finished Integration!"); //if so, the integration is complete, and we can exit out.
            CCTK_VINFO("FINAL: Mass :,\t%15.14e,\t",TOVola_y[2]); //See what the final mass is.
            }
    }

    gsl_odeiv2_driver_free(TOVola_d); //free the driver
    CCTK_INFO("ODE Solver using GSL for TOVola Shutting Down...\n"); //Shut down the Second integration
    return;
    
}

//Allocate space for solution  
void TOVola_TOV_Allocate(CCTK_ARGUMENTS){

	CCTK_INFO("TOVola Allocating required Memory");
	TOVola_Raw_rSchw      =    (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*TOVola_Numpoints);
	TOVola_Raw_rho_energy =    (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*TOVola_Numpoints);
	TOVola_Raw_rho_baryon =    (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*TOVola_Numpoints);
	TOVola_Raw_P 	       =    (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*TOVola_Numpoints);
	TOVola_Raw_M 	       =    (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*TOVola_Numpoints);
	TOVola_Raw_nu 	       =    (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*TOVola_Numpoints);
	TOVola_Raw_Iso_r      =    (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*TOVola_Numpoints);
	return;
}

//Free Memory
void TOVola_TOV_Cleanup(CCTK_ARGUMENTS){
	
	CCTK_INFO("TOVola's Final Cleanup: Freeing Memory");
	if (TOVola_Raw_rSchw != NULL) { free(TOVola_Raw_rSchw);}
	if (TOVola_Raw_rho_energy != NULL) { free(TOVola_Raw_rho_energy);}
	if (TOVola_Raw_rho_baryon != NULL) { free(TOVola_Raw_rho_baryon);}
	if (TOVola_Raw_P != NULL) { free(TOVola_Raw_P);}
	if (TOVola_Raw_M != NULL) { free(TOVola_Raw_M);}
	if (TOVola_Raw_nu != NULL) { free(TOVola_Raw_nu);}
	if (TOVola_Raw_Iso_r != NULL) { free(TOVola_Raw_Iso_r);
	TOVola_Numpoints = 0;}//This last point is for if you are planning to do multiple TOVs, or work on multiple processors.
	
	//Finish up!
  	CCTK_INFO("TOVola Shutting down...\n");
	return;
}

//Useful defines for later. Inspired from ET TOVSolver
#define velx (&vel[0*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])
#define vely (&vel[1*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])
#define velz (&vel[2*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])
#define velx_p (&vel_p[0*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])
#define vely_p (&vel_p[1*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])
#define velz_p (&vel_p[2*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])
#define velx_p_p (&vel_p_p[0*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])
#define vely_p_p (&vel_p_p[1*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])
#define velz_p_p (&vel_p_p[2*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])


//This function is for normalizing rbar, and calculating expnu and exp4phi for use in the toolkit.
void TOVola_Normalize_and_set_data(CCTK_REAL *restrict r_Schw_arr,
				 CCTK_REAL *restrict rho_energy_arr,
				 CCTK_REAL *restrict rho_baryon_arr,
				 CCTK_REAL *restrict P_arr,
				 CCTK_REAL *restrict M_arr,
				 CCTK_REAL *restrict expnu_arr,
				 CCTK_REAL *restrict exp4phi_arr,
				 CCTK_REAL *restrict rbar_arr){

  DECLARE_CCTK_PARAMETERS

  int this_line = 0;
  //Start with the raw data
  while (this_line < TOVola_Numpoints){
    r_Schw_arr[this_line]     = TOVola_Raw_rSchw[this_line];
    rho_energy_arr[this_line] = TOVola_Raw_rho_energy[this_line];
    rho_baryon_arr[this_line] = TOVola_Raw_rho_baryon[this_line];
    P_arr[this_line]          = TOVola_Raw_P[this_line];
    M_arr[this_line]          = TOVola_Raw_M[this_line];
    expnu_arr[this_line]      = TOVola_Raw_nu[this_line];
    rbar_arr[this_line]       = TOVola_Raw_Iso_r[this_line];	
    this_line++;
  }

  CCTK_INFO("TOVola Normalizing raw TOV data..");
  
  //Surface values are going to be important when normalizing
  double R_Schw_surface = r_Schw_arr[TOVola_Numpoints-1];
  double M_surface      = M_arr[TOVola_Numpoints-1];
  double Rbar_surface   = rbar_arr[TOVola_Numpoints-1];
  double nu_surface     = expnu_arr[TOVola_Numpoints-1]; //Remember, nu hasn't been adjusted yet, so the expnu_arr just refers to nu itself.
  double normalize      = 0.5*(sqrt(R_Schw_surface*(R_Schw_surface-2.0*M_surface))+R_Schw_surface-M_surface)/Rbar_surface;

  //Start with the first values, since rbar center is going to be problematic (r=0)
  rbar_arr[0] = rbar_arr[0]*normalize;
  rbar_arr[1] = rbar_arr[1]*normalize;
  expnu_arr[0] = exp(expnu_arr[0] - nu_surface + log(1 -2.0*M_surface/R_Schw_surface));
  expnu_arr[1] = exp(expnu_arr[1] - nu_surface + log(1 -2.0*M_surface/R_Schw_surface));
  exp4phi_arr[1] = pow((r_Schw_arr[1]/rbar_arr[1]),2.0);
  exp4phi_arr[0] = exp4phi_arr[1];
  
  //Now we can normalize the rest of rbars and get the remaining expnu and exp4phi.
  for(int i=2;i<this_line;i++){
      rbar_arr[i]    = rbar_arr[i]*normalize;
      expnu_arr[i]   = exp(expnu_arr[i] - nu_surface + log(1 -2.0*M_surface/R_Schw_surface));
      exp4phi_arr[i] = pow((r_Schw_arr[i]/rbar_arr[i]),2.0);
    }

  CCTK_INFO("Normalization of raw data complete!");
  
}

// Find interpolation index using Bisection root-finding algorithm:
// Generated in nrpytutorial, with very minor edits made
static inline int TOVola_bisection_idx_finder(const CCTK_REAL rrbar, const int numlines_in_file, const CCTK_REAL *restrict rbar_arr) {
  int x1 = 0;
  int x2 = numlines_in_file-1;
  CCTK_REAL y1 = rrbar-rbar_arr[x1];
  CCTK_REAL y2 = rrbar-rbar_arr[x2];
  if(y1*y2 > 0) {
    CCTK_VINFO("INTERPOLATION BRACKETING ERROR %e | %e %e\n",rrbar,y1,y2);
    CCTK_ERROR("Shutting down due to error");
  }
  for(int i=0;i<numlines_in_file;i++) {
    int x_midpoint = (x1+x2)/2;
    CCTK_REAL y_midpoint = rrbar-rbar_arr[x_midpoint];
    if(y_midpoint*y1 <= 0) {
      x2 = x_midpoint;
      y2 = y_midpoint;
    } else {
      x1 = x_midpoint;
      y1 = y_midpoint;
    }
    if( abs(x2-x1) == 1 ) {
      // If rbar_arr[x1] is closer to rrbar than rbar_arr[x2] then return x1:
      if(fabs(rrbar-rbar_arr[x1]) < fabs(rrbar-rbar_arr[x2])) {return x1;}
      // Otherwiser return x2:
      return x2;
    }
  }
  CCTK_ERROR("INTERPOLATION BRACKETING ERROR: DID NOT CONVERGE.\n");
}

//Interpolate to grid point
//Generated in nrpytutorial, with edits made for compatibility in TOVola
void TOVola_TOV_interpolate_1D(CCTK_REAL rrbar,const CCTK_REAL Rbar,const int Rbar_idx,const int interp_stencil_size,
                        const int numlines_in_file,const CCTK_REAL *restrict r_Schw_arr,const CCTK_REAL *restrict rho_energy_arr,const CCTK_REAL *restrict rho_baryon_arr,const CCTK_REAL *restrict P_arr,
                        const CCTK_REAL *restrict M_arr,const CCTK_REAL *restrict expnu_arr,const CCTK_REAL *restrict exp4phi_arr,const CCTK_REAL *restrict rbar_arr,
                        CCTK_REAL *restrict rho_energy,CCTK_REAL *restrict rho_baryon,CCTK_REAL *restrict P,CCTK_REAL *restrict M,CCTK_REAL *restrict expnu,CCTK_REAL *restrict exp4phi) {

  // For this case, we know that for all functions, f(r) = f(-r)
  if(rrbar < 0) rrbar = -rrbar;

  // First find the central interpolation stencil index:
  int idx = TOVola_bisection_idx_finder(rrbar,numlines_in_file,rbar_arr);


#ifdef MAX
#undef MAX
#endif
#define MAX(A, B) ( ((A) > (B)) ? (A) : (B) )

  int idxmin = MAX(0,idx-interp_stencil_size/2-1);

#ifdef MIN
#undef MIN
#endif
#define MIN(A, B) ( ((A) < (B)) ? (A) : (B) )

  // -= Do not allow the interpolation stencil to cross the star's surface =-
  // max index is when idxmin + (interp_stencil_size-1) = Rbar_idx
  //  -> idxmin at most can be Rbar_idx - interp_stencil_size + 1
  if(rrbar < Rbar) {
    idxmin = MIN(idxmin,Rbar_idx - interp_stencil_size + 1);
  } else {
    idxmin = MAX(idxmin,Rbar_idx+1);
    idxmin = MIN(idxmin,numlines_in_file - interp_stencil_size + 1);
  }
  // Now perform the Lagrange polynomial interpolation:

  // First set the interpolation coefficients:
  CCTK_REAL rbar_sample[interp_stencil_size];
  for(int i=idxmin;i<idxmin+interp_stencil_size;i++) {
    rbar_sample[i-idxmin] = rbar_arr[i];
  }
  CCTK_REAL l_i_of_r[interp_stencil_size];
  for(int i=0;i<interp_stencil_size;i++) {
    CCTK_REAL numer = 1.0;
    CCTK_REAL denom = 1.0;
    for(int j=0;j<i;j++) {
      numer *= rrbar - rbar_sample[j];
      denom *= rbar_sample[i] - rbar_sample[j];
    }
    for(int j=i+1;j<interp_stencil_size;j++) {
      numer *= rrbar - rbar_sample[j];
      denom *= rbar_sample[i] - rbar_sample[j];
    }
    l_i_of_r[i] = numer/denom;
  }

  // Then perform the interpolation:
  *rho_energy = 0.0;
  *rho_baryon = 0.0;
  *P = 0.0;
  *M = 0.0;
  *expnu = 0.0;
  *exp4phi = 0.0;

  CCTK_REAL r_Schw = 0.0;
  for(int i=idxmin;i<idxmin+interp_stencil_size;i++) {
    r_Schw      += l_i_of_r[i-idxmin] * r_Schw_arr[i];
    *rho_energy += l_i_of_r[i-idxmin] * rho_energy_arr[i];
    *rho_baryon += l_i_of_r[i-idxmin] * rho_baryon_arr[i];
    *P          += l_i_of_r[i-idxmin] * P_arr[i];
    *M          += l_i_of_r[i-idxmin] * M_arr[i];
    *expnu      += l_i_of_r[i-idxmin] * expnu_arr[i];
    *exp4phi    += l_i_of_r[i-idxmin] * exp4phi_arr[i];
  }

  //Jusst in case we are at the surface point.
  if(rrbar > Rbar) {
    *rho_energy = 0;
    *rho_baryon = 0;
    *P          = 0;
    *M          = M_arr[Rbar_idx+1];
    *expnu      = 1. - 2.*(*M) / r_Schw;
    *exp4phi    = pow(r_Schw / rrbar,2.0);
  }
}

//For populating time levels
//Ripped directly from original TOVsolver in the toolkit
void TOVola_TOV_Copy(CCTK_INT size, CCTK_REAL *var_p, CCTK_REAL *var)
{
#pragma omp parallel for
    for(int i=0; i<size; i++)
        var_p[i] = var[i];
}

//interpollate TOV data to the grid. This is the second half of TOVola, and where you get the data you actually want to use.
void TOVola_interp_Driver(CCTK_ARGUMENTS){

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  
  //Allocate memory for data storage
  CCTK_REAL *TOVola_r_Schw_arr =     (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*TOVola_Numpoints);
  CCTK_REAL *TOVola_rho_energy_arr = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*TOVola_Numpoints);
  CCTK_REAL *TOVola_rho_baryon_arr = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*TOVola_Numpoints);
  CCTK_REAL *TOVola_P_arr =          (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*TOVola_Numpoints);
  CCTK_REAL *TOVola_M_arr =          (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*TOVola_Numpoints);
  CCTK_REAL *TOVola_expnu_arr =      (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*TOVola_Numpoints);
  CCTK_REAL *TOVola_exp4phi_arr =    (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*TOVola_Numpoints);
  CCTK_REAL *TOVola_rbar_arr =       (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*TOVola_Numpoints);

  //Normalize and set arrays
  TOVola_Normalize_and_set_data(TOVola_r_Schw_arr,TOVola_rho_energy_arr,TOVola_rho_baryon_arr,TOVola_P_arr,TOVola_M_arr,TOVola_expnu_arr,TOVola_exp4phi_arr,TOVola_rbar_arr);
  
  //get surface values for later calculations.
  CCTK_REAL TOVola_Rbar = TOVola_rbar_arr[TOVola_Numpoints-1]; 
  CCTK_REAL TOVola_Mass = TOVola_M_arr[TOVola_Numpoints-1];
  
  
  //Now for the actual grid placements. Go over all grid points
  CCTK_INFO("TOVola Beginning Grid Placements...");
  for(int i=0; i<cctk_lsh[0]; i++){ 
  	for(int j=0; j<cctk_lsh[1]; j++){ 
  		for(int k=0; k<cctk_lsh[2]; k++){ 
  			int i3d=CCTK_GFINDEX3D(cctkGH,i,j,k); //3D index
  			CCTK_REAL TOVola_rGrid = sqrt((x[i3d]*x[i3d])+(y[i3d]*y[i3d])+(z[i3d]*z[i3d])); //magnitude of r on the grid
  			CCTK_REAL TOVola_rho_energy, TOVola_rho_baryon, TOVola_P, TOVola_M, TOVola_expnu, TOVola_exp4phi; //Declare TOV quantities
  			if (TOVola_rGrid < TOVola_Rbar){ //If we are INSIDE the star, we need to interpollate the data to the grid.
  				TOVola_TOV_interpolate_1D(TOVola_rGrid, TOVola_Rbar, TOVola_Numpoints-1,TOVola_Interpolation_Stencil,TOVola_Numpoints, 
  					TOVola_r_Schw_arr,TOVola_rho_energy_arr,TOVola_rho_baryon_arr,TOVola_P_arr,TOVola_M_arr,TOVola_expnu_arr,TOVola_exp4phi_arr,TOVola_rbar_arr,
  					&TOVola_rho_energy,&TOVola_rho_baryon,&TOVola_P,&TOVola_M,&TOVola_expnu,&TOVola_exp4phi);		
  				rho[i3d] = TOVola_rho_baryon;
				press[i3d] = TOVola_P;
				// tiny number prevents 0/0.
				eps[i3d] = (TOVola_rho_energy / (TOVola_rho_baryon+1e-30)) - 1.0;
				if (eps[i3d]<0){eps[i3d]=0.0;}
				alp[i3d] = pow(TOVola_expnu,0.5);//This is the lapse
				gxx[i3d] = TOVola_exp4phi;//This is the values for the metric in the coordinates we chose.
				gyy[i3d] = gxx[i3d];
				gzz[i3d] = gxx[i3d];}
			else { //If we are OUTSIDE the star, we need to calculate the grid functions directly.
				CCTK_REAL TOVola_rSchw_outside = (TOVola_rGrid+TOVola_Mass) + TOVola_Mass*TOVola_Mass/(4.0*TOVola_rGrid);//Need to know what rSchw is at our current grid location.
				rho[i3d] = 0.0;
				press[i3d] = 0.0;
				eps[i3d] = 0.0;
				alp[i3d] = pow(1-2*TOVola_Mass/TOVola_rSchw_outside,0.5); //Goes to Schwarschild
				gxx[i3d] = pow((TOVola_rSchw_outside/TOVola_rGrid),2.0);
				gyy[i3d] = gxx[i3d];
				gzz[i3d] = gxx[i3d];
				}
			
			betax[i3d] = 0.0;
			betay[i3d] = 0.0;
			betaz[i3d] = 0.0;
			gxy[i3d] = 0.0;
			gxz[i3d] = 0.0;
			gyz[i3d] = 0.0;
			//Curvature is zero for our slice.
			kxx[i3d] = 0.0;
			kyy[i3d] = 0.0;
			kzz[i3d] = 0.0;
			kxy[i3d] = 0.0;
			kxz[i3d] = 0.0;
			kyz[i3d] = 0.0;
			//velocities are zero: It is a static solution.
			velx[i3d] = 0.0;
			vely[i3d] = 0.0;
			velz[i3d] = 0.0;
	        	w_lorentz[i3d] = 1.0;
  		}
  	}
  }
  			
  CCTK_INFO("Grid Placement Successful!");
  
  //This piece of the code was VERY heavily inspired by the original ET TOVSolver, with minor edits for TOVola
  CCTK_INFO("TOVola Finalizing Grid...\n");
  int i3d = cctk_lsh[2]*cctk_lsh[1]*cctk_lsh[0];
  switch(TOVola_TOV_Populate_Timelevels)
  {
    case 3:
        TOVola_TOV_Copy(i3d, gxx_p_p,  gxx);
        TOVola_TOV_Copy(i3d, gyy_p_p,  gyy);
        TOVola_TOV_Copy(i3d, gzz_p_p,  gzz);
        TOVola_TOV_Copy(i3d, gxy_p_p,  gxy);
        TOVola_TOV_Copy(i3d, gxz_p_p,  gxz);
        TOVola_TOV_Copy(i3d, gyz_p_p,  gyz);
        TOVola_TOV_Copy(i3d, rho_p_p,  rho);
        TOVola_TOV_Copy(i3d, eps_p_p,  eps);
        TOVola_TOV_Copy(i3d, velx_p_p, velx);
        TOVola_TOV_Copy(i3d, vely_p_p, vely);
        TOVola_TOV_Copy(i3d, velz_p_p, velz);
        TOVola_TOV_Copy(i3d, w_lorentz_p_p, w_lorentz);
        // fall through
    case 2:
        TOVola_TOV_Copy(i3d, gxx_p,  gxx);
        TOVola_TOV_Copy(i3d, gyy_p,  gyy);
        TOVola_TOV_Copy(i3d, gzz_p,  gzz);
        TOVola_TOV_Copy(i3d, gxy_p,  gxy);
        TOVola_TOV_Copy(i3d, gxz_p,  gxz);
        TOVola_TOV_Copy(i3d, gyz_p,  gyz);
        TOVola_TOV_Copy(i3d, rho_p,  rho);
        TOVola_TOV_Copy(i3d, eps_p,  eps);
        TOVola_TOV_Copy(i3d, velx_p, velx);
        TOVola_TOV_Copy(i3d, vely_p, vely);
        TOVola_TOV_Copy(i3d, velz_p, velz);
        TOVola_TOV_Copy(i3d, w_lorentz_p, w_lorentz);
        // fall through
    case 1:
        break;
    default:
        CCTK_VWARN(CCTK_WARN_ABORT,
                   "Unsupported number of TOVola_TOV_Populate_TimelevelsL: %d",
                   (int)TOVola_TOV_Populate_Timelevels);
        break;
  }
  
  //Let's free up the memory
  free(TOVola_r_Schw_arr);
  free(TOVola_rho_energy_arr);
  free(TOVola_rho_baryon_arr);
  free(TOVola_P_arr);
  free(TOVola_M_arr);
  free(TOVola_expnu_arr);
  free(TOVola_exp4phi_arr);
  free(TOVola_rbar_arr);
  
  //Finish up. Goodbye!
  CCTK_INFO("Complete! Enjoy your Initial Data!\n\n");
}

//TOVola: By David Boyer
//Code made possible by: GSL (ODE solver), NRpy (Generated grid interpollator), GRHayL (EOS library), and the original ET TOVSolver (For reference and a starting point learning the Toolkit)
