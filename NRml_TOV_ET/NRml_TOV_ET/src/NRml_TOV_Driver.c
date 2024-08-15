#include "NRml_TOV_funcs.h" //nrpy_odiegm itself, for use in NRml_TOV.
#include "NRml_TOV_user_methods.h" //Edited specifically for the TOV solver.
#include "GRHayLib.h" //Access to GRHayL library in the ET
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

//IMPORTANT: NRml_TOV_ET src code does NOT need to be edited. Look at parfiles and param.ccl for use.

/********************************
// NRml_TOV is a special version of the Odie main driver. It is built to handle the TOV equations specifically, and is built to handle a variety of different types of EOS
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
// EDIT: NRml_TOV_ET is The Einstein Toolkit version of the solver, which also places the TOV data to the ET grid.
// One does not need to edit any of the source code.
// Users who wish to interface with this program need only create a parfile for a simulation. Output variables will be adjusted for use in ADMbase and Hydrobase as initial data.
// IllinoisGRMHD is used to set Tmunu from the data and Baikal tests the constraint Violation.
// 
// David Boyer (6/21/24)
********************************/

//Global pointers to save data. (Deal with it.)
CCTK_REAL *NRml_Raw_rSchw = NULL;
CCTK_REAL *NRml_Raw_rho_energy = NULL;
CCTK_REAL *NRml_Raw_rho_baryon = NULL;
CCTK_REAL *NRml_Raw_P = NULL;
CCTK_REAL *NRml_Raw_M = NULL;
CCTK_REAL *NRml_Raw_nu = NULL;
CCTK_REAL *NRml_Raw_Iso_r = NULL;
int NRml_Numpoints = 0;


void NRml_TOV_Integrator_1(CCTK_ARGUMENTS){

    //Einstein Toolkit main function

    //Declare value from ET
    DECLARE_CCTK_PARAMETERS
    CCTK_INFO("First Integration: Analyzing how much memory to allocate");
    CCTK_INFO("Beginning ODE Solver \"Odie\" for NRml_TOV_ET...");

    double rhoCentral_baryon = NRml_central_baryon_density;

    //Values for ODE steps and Error
    double step = NRml_step;
    double current_position = 0.0;
    const int size = NRml_size;
    int adams_bashforth_order = 4; //Only edit if you are using an AB method (NOT RECCOMMENDED)
    bool adaptive_step = true;
    double error_limit = NRml_error_limit;

    //Choose your desired ODE method
    //Thorn only has two methods available: (ARKF: (RK4(5))) or (ADP8: (DP7(8)))
    const nrpy_odiegm_step_type * step_type;

    if (CCTK_EQUALS("ARKF",NRml_ODE_method)) {
       step_type = nrpy_odiegm_step_ARKF;}
    else if (CCTK_EQUALS("ADP8",NRml_ODE_method)) {
       step_type = nrpy_odiegm_step_ADP8;}
    else{
      CCTK_INFO("Invalid Step type. May only use ARKF or ADP8.");
      CCTK_INFO("Shutting down due to error.\n");
      exit(1);}

    //Initializing the EOS
    if(CCTK_EQUALS("Simple",NRml_EOS_type)){
      CCTK_INFO("Simple Polytrope");
      if(ghl_eos->neos!=1){
        CCTK_INFO("Error: Too many regions for the simple polytrope.");
        CCTK_INFO("Check your value for neos, or use a piecewise polytrope (type = 'p')");
        CCTK_ERROR("Shutting down due to error...\n");}
        }
    else if(CCTK_EQUALS("Piecewise",NRml_EOS_type)){
      CCTK_INFO("Piecewise Polytrope");}
    else if(CCTK_EQUALS("Tabulated",NRml_EOS_type)){
      CCTK_INFO("Tabulated EOS");
      ghl_tabulated_compute_Ye_P_eps_of_rho_beq_constant_T(NRml_Tin, ghl_eos);
      //printf("hit\n");
    }
    else{
      CCTK_INFO("ERROR: Invalid EOS type. Must be either 'Simple', 'Piecewise', or 'Tabulated'");
      CCTK_ERROR("Shutting down due to error...\n");}


    const nrpy_odiegm_step_type * step_type_2;
    step_type_2 = step_type;

    bool no_adaptive_step;
    if(adaptive_step==true){
      no_adaptive_step=false;
      if(step_type != nrpy_odiegm_step_AHE
	 && step_type != nrpy_odiegm_step_ABS
	 && step_type != nrpy_odiegm_step_ARKF
	 && step_type != nrpy_odiegm_step_ACK
	 && step_type != nrpy_odiegm_step_ADP5
	 && step_type != nrpy_odiegm_step_ADP8){
	CCTK_INFO("ERROR: Incompatible step type with adaptive method.");
	CCTK_ERROR("Shutting down due to error.\n");}
    }
    else{
    	no_adaptive_step=true;
      if(step_type == nrpy_odiegm_step_AHE
         || step_type == nrpy_odiegm_step_ABS
         || step_type == nrpy_odiegm_step_ARKF
         || step_type == nrpy_odiegm_step_ACK
         || step_type == nrpy_odiegm_step_ADP5
         || step_type == nrpy_odiegm_step_ADP8){
        CCTK_INFO("ERROR: Incompatible step type with non-adaptive method.");
        CCTK_ERROR("Shutting down due to error.\n");}
    }
    
    //Define an ODE struct
    double absolute_error_limit = error_limit;
    double relative_error_limit = error_limit;
    int number_of_equations = 4;
    int number_of_constants = 2;
    struct constant_parameters cp = { 
      .dimension = number_of_constants,
    };
    cp.rhoCentral_baryon = rhoCentral_baryon;
    

    nrpy_odiegm_system system = {diffy_Q_eval,known_Q_eval,number_of_equations,&cp};
    

    nrpy_odiegm_driver *d;
    d = nrpy_odiegm_driver_alloc_y_new(&system, step_type, step, absolute_error_limit, relative_error_limit); 

    //setting the control struct. See OdieGM and the GSL ODE solver for more info.
    d->c->scale_factor = NRml_scale_factor;
    d->c->error_safety = NRml_error_safety;
    d->c->ay_error_scaler = NRml_ay_error_scaler;
    d->c->ady_error_scaler = NRml_ady_error_scaler;
    d->c->max_step_adjustment = NRml_max_step_adjustment;
    d->c->min_step_adjustment = NRml_min_step_adjustment;
    d->c->absolute_max_step = NRml_absolute_max_step;
    d->c->absolute_min_step = NRml_absolute_min_step;
    d->c->error_upper_tolerance = NRml_error_upper_tolerance;
    d->c->error_lower_tolerance = NRml_error_lower_tolerance;

    //Some OdieGM checks for AB methods
    int method_type = 1;
    if (step_type->rows == step_type->columns) {
        method_type = 0;  
    } 
    if (step_type->rows == 19) { 
        method_type = 2;
    } else {
        adams_bashforth_order = 0;
    }
    d->s->adams_bashforth_order = adams_bashforth_order;
    d->e->no_adaptive_step = no_adaptive_step;
    

    if (method_type == 2) {
        CCTK_VINFO("Method Order: %i.",adams_bashforth_order);
    } else {
        CCTK_VINFO("Method Order: %i.",step_type->order);            
    }

    //Declare ODEs and begin SOLVING!
    double y[number_of_equations];
    double c[number_of_constants];
    
    get_initial_condition(y,&cp); 
    assign_constants(c,&cp); 

    
    //For the first integration, we are just counting up how much memory we need to allocate.
    NRml_Numpoints++;

    for (int i = 0; i < size; i++){
        
        if (method_type == 2 && i == 0 && step_type_2 != nrpy_odiegm_step_AB) {
            d->s->type = step_type_2;
            d->s->rows = step_type_2->rows;
            d->s->columns = step_type_2->columns;
            d->s->method_type = 0;
            d->s->adams_bashforth_order = adams_bashforth_order;
            d->e->no_adaptive_step = true;
        } else if (step_type != step_type_2 && method_type == 2 && i == adams_bashforth_order) {
            d->s->type = step_type;
            d->s->rows = step_type->rows;
            d->s->columns = step_type->columns;
            d->s->method_type = 2;
            d->s->adams_bashforth_order = adams_bashforth_order;
            d->e->no_adaptive_step = true;
        }

        nrpy_odiegm_evolve_apply(d->e, d->c, d->s, &system, &current_position, current_position+step, &step, y);
        

        exception_handler(current_position,y);
        const_eval(current_position,y,&cp);
        assign_constants(c,&cp);
    	
    	NRml_Numpoints++;

        if (do_we_terminate(current_position, y, &cp) == 1) {
            i = size-1;
             
        } 
        if (i == size-1) { 
            CCTK_INFO("Finished Integration!");
            CCTK_VINFO("FINAL: Mass :,\t%15.14e,\t",y[2]);
            }
    }


    //fclose(fp2);

    nrpy_odiegm_driver_free(d);
    CCTK_INFO("ODE Solver \"Odie\" for NRml_TOV_ET Shutting Down...\n");
    return;
    
} // -GM, master of dogs (OdieGM code)
  // -David Boyer, Edited for NRml_TOV and NRml_TOV_ET
  
void NRml_TOV_Integrator_2(CCTK_ARGUMENTS){

    //Einstein Toolkit main function

    //Declare value from ET
    DECLARE_CCTK_PARAMETERS
    CCTK_INFO("Second Integration: Writing solution to Memory");
    CCTK_INFO("Beginning ODE Solver \"Odie\" for NRml_TOV_ET...");

    int NRml_this_point = 0;
    double rhoCentral_baryon = NRml_central_baryon_density;

    //Values for ODE steps and Error
    double step = NRml_step;
    double current_position = 0.0;
    const int size = NRml_size;
    int adams_bashforth_order = 4; //Only edit if you are using an AB method (NOT RECCOMMENDED)
    bool adaptive_step = true;
    double error_limit = NRml_error_limit;

    //Choose your desired ODE method
    //Thorn only has two methods available: (ARKF: (RK4(5))) or (ADP8: (DP7(8)))
    const nrpy_odiegm_step_type * step_type;

    if (CCTK_EQUALS("ARKF",NRml_ODE_method)) {
       step_type = nrpy_odiegm_step_ARKF;}
    else if (CCTK_EQUALS("ADP8",NRml_ODE_method)) {
       step_type = nrpy_odiegm_step_ADP8;}
    else{
      CCTK_INFO("Invalid Step type. May only use ARKF or ADP8.");
      CCTK_ERROR("Shutting down due to error.\n");}

    //Initializing the EOS
    if(CCTK_EQUALS("Simple",NRml_EOS_type)){
      CCTK_INFO("Simple Polytrope");
      if(ghl_eos->neos!=1){
        CCTK_INFO("Error: Too many regions for the simple polytrope.");
        CCTK_INFO("Check your value for neos, or use a piecewise polytrope (type = 'p')");
        CCTK_ERROR("Shutting down due to error...\n");}
        }
    else if(CCTK_EQUALS("Piecewise",NRml_EOS_type)){
      CCTK_INFO("Piecewise Polytrope");}
    else if(CCTK_EQUALS("Tabulated",NRml_EOS_type)){
      CCTK_INFO("Tabulated EOS");
      ghl_tabulated_compute_Ye_P_eps_of_rho_beq_constant_T(NRml_Tin, ghl_eos);
      //printf("hit\n");
    }
    else{
      CCTK_INFO("ERROR: Invalid EOS type. Must be either 'Simple', 'Piecewise', or 'Tabulated'");
      CCTK_ERROR("Shutting down due to error...\n");}


    const nrpy_odiegm_step_type * step_type_2;
    step_type_2 = step_type;

    bool no_adaptive_step;
    if(adaptive_step==true){
      no_adaptive_step=false;
      if(step_type != nrpy_odiegm_step_AHE
	 && step_type != nrpy_odiegm_step_ABS
	 && step_type != nrpy_odiegm_step_ARKF
	 && step_type != nrpy_odiegm_step_ACK
	 && step_type != nrpy_odiegm_step_ADP5
	 && step_type != nrpy_odiegm_step_ADP8){
	CCTK_INFO("ERROR: Incompatible step type with adaptive method.");
	CCTK_ERROR("Shutting down due to error.\n");}
    }

    if(adaptive_step==false){
      no_adaptive_step=true;
      if(step_type == nrpy_odiegm_step_AHE
         || step_type == nrpy_odiegm_step_ABS
         || step_type == nrpy_odiegm_step_ARKF
         || step_type == nrpy_odiegm_step_ACK
         || step_type == nrpy_odiegm_step_ADP5
         || step_type == nrpy_odiegm_step_ADP8){
        CCTK_INFO("ERROR: Incompatible step type with non-adaptive method.");
        CCTK_ERROR("Shutting down due to error.\n");}
    }

    //Define an ODE struct
    double absolute_error_limit = error_limit;
    double relative_error_limit = error_limit;
    int number_of_equations = 4;
    int number_of_constants = 2;
    struct constant_parameters cp = { 
      .dimension = number_of_constants,
    };
    cp.rhoCentral_baryon = rhoCentral_baryon;
    

    nrpy_odiegm_system system = {diffy_Q_eval,known_Q_eval,number_of_equations,&cp};
    

    nrpy_odiegm_driver *d;
    d = nrpy_odiegm_driver_alloc_y_new(&system, step_type, step, absolute_error_limit, relative_error_limit); 

    //setting the control struct. See OdieGM and the GSL ODE solver for more info.
    d->c->scale_factor = NRml_scale_factor;
    d->c->error_safety = NRml_error_safety;
    d->c->ay_error_scaler = NRml_ay_error_scaler;
    d->c->ady_error_scaler = NRml_ady_error_scaler;
    d->c->max_step_adjustment = NRml_max_step_adjustment;
    d->c->min_step_adjustment = NRml_min_step_adjustment;
    d->c->absolute_max_step = NRml_absolute_max_step;
    d->c->absolute_min_step = NRml_absolute_min_step;
    d->c->error_upper_tolerance = NRml_error_upper_tolerance;
    d->c->error_lower_tolerance = NRml_error_lower_tolerance;

    //Some OdieGM checks for AB methods
    int method_type = 1;
    if (step_type->rows == step_type->columns) {
        method_type = 0;  
    } 
    if (step_type->rows == 19) { 
        method_type = 2;
    } else {
        adams_bashforth_order = 0;
    }
    d->s->adams_bashforth_order = adams_bashforth_order;
    d->e->no_adaptive_step = no_adaptive_step;
    

    if (method_type == 2) {
        CCTK_VINFO("Method Order: %i.",adams_bashforth_order);
    } else {
        CCTK_VINFO("Method Order: %i.",step_type->order);            
    }

    //Declare ODEs and begin SOLVING!
    double y[number_of_equations];
    double c[number_of_constants];
    
    get_initial_condition(y,&cp); 
    assign_constants(c,&cp); 

    //For the Second Integration, we are actually writing the solution to memory
    NRml_Raw_rSchw[NRml_this_point] = current_position;
    NRml_Raw_rho_energy[NRml_this_point] = c[0];
    NRml_Raw_rho_baryon[NRml_this_point] = c[1];
    NRml_Raw_P[NRml_this_point] = y[0];
    NRml_Raw_M[NRml_this_point] = y[2];
    NRml_Raw_nu[NRml_this_point] = y[1];
    NRml_Raw_Iso_r[NRml_this_point] = y[3];
    NRml_this_point++;

    for (int i = 0; i < size; i++){
        
        if (method_type == 2 && i == 0 && step_type_2 != nrpy_odiegm_step_AB) {
            d->s->type = step_type_2;
            d->s->rows = step_type_2->rows;
            d->s->columns = step_type_2->columns;
            d->s->method_type = 0;
            d->s->adams_bashforth_order = adams_bashforth_order;
            d->e->no_adaptive_step = true;
        } else if (step_type != step_type_2 && method_type == 2 && i == adams_bashforth_order) {
            d->s->type = step_type;
            d->s->rows = step_type->rows;
            d->s->columns = step_type->columns;
            d->s->method_type = 2;
            d->s->adams_bashforth_order = adams_bashforth_order;
            d->e->no_adaptive_step = true;
        }

        nrpy_odiegm_evolve_apply(d->e, d->c, d->s, &system, &current_position, current_position+step, &step, y);
        

        exception_handler(current_position,y);
        const_eval(current_position,y,&cp);
        assign_constants(c,&cp);

	NRml_Raw_rSchw[NRml_this_point] = current_position;
    	NRml_Raw_rho_energy[NRml_this_point] = c[0];
    	NRml_Raw_rho_baryon[NRml_this_point] = c[1];
    	NRml_Raw_P[NRml_this_point] = y[0];
    	NRml_Raw_M[NRml_this_point] = y[2];
    	NRml_Raw_nu[NRml_this_point] = y[1];
    	NRml_Raw_Iso_r[NRml_this_point] = y[3];
    	NRml_this_point++;

        if (do_we_terminate(current_position, y, &cp) == 1) {
            i = size-1;
             
        } 
        if (i == size-1) { 
            CCTK_INFO("Finished Integration!");
            CCTK_VINFO("FINAL: Mass :,\t%15.14e,\t",y[2]);
        }
    }


    //fclose(fp2);

    nrpy_odiegm_driver_free(d);
    CCTK_INFO("ODE Solver \"Odie\" for NRml_TOV_ET Shutting Down...\n");
    return;
    
} // -GM, master of dogs (OdieGM code)
  // -David Boyer, Edited for NRml_TOV and NRml_TOV_ET

#define REAL double

//Allocate space for solution  
void NRml_TOV_Allocate(CCTK_ARGUMENTS){

	CCTK_INFO("NRml Allocating required Memory");
	NRml_Raw_rSchw      =    (REAL *)malloc(sizeof(REAL)*NRml_Numpoints);
	NRml_Raw_rho_energy =    (REAL *)malloc(sizeof(REAL)*NRml_Numpoints);
	NRml_Raw_rho_baryon =    (REAL *)malloc(sizeof(REAL)*NRml_Numpoints);
	NRml_Raw_P 	    =    (REAL *)malloc(sizeof(REAL)*NRml_Numpoints);
	NRml_Raw_M 	    =    (REAL *)malloc(sizeof(REAL)*NRml_Numpoints);
	NRml_Raw_nu 	    =    (REAL *)malloc(sizeof(REAL)*NRml_Numpoints);
	NRml_Raw_Iso_r      =    (REAL *)malloc(sizeof(REAL)*NRml_Numpoints);
	return;
}

//Free Memory
void NRml_TOV_Cleanup(CCTK_ARGUMENTS){
	
	CCTK_INFO("NRml's Final Cleanup: Freeing Memory");
	if (NRml_Raw_rSchw != NULL) { free(NRml_Raw_rSchw);}
	if (NRml_Raw_rho_energy != NULL) { free(NRml_Raw_rho_energy);}
	if (NRml_Raw_rho_baryon != NULL) { free(NRml_Raw_rho_baryon);}
	if (NRml_Raw_P != NULL) { free(NRml_Raw_P);}
	if (NRml_Raw_M != NULL) { free(NRml_Raw_M);}
	if (NRml_Raw_nu != NULL) { free(NRml_Raw_nu);}
	if (NRml_Raw_Iso_r != NULL) { free(NRml_Raw_Iso_r);
	NRml_Numpoints = 0;}
	
	//Finish up!
  	CCTK_INFO("NRml_TOV_ET Shutting down...\n");
	return;
}


//These extra functions are for use in the toolkit
//The raw TOV data is read in and converted into ADM and Hydrobase variables and placed on a grid.

//Grid interpollation functions generated by NRpy with edits made for ET compatibility
//Grid placements inspired heavily by ET's original TOVsolver in EinsteinInitialData by Hawke and Loeffler 2009, So NRml_TOV_ET can act as a drop-in replacement.

//First, defines
#define velx (&vel[0*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])
#define vely (&vel[1*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])
#define velz (&vel[2*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])
#define velx_p (&vel_p[0*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])
#define vely_p (&vel_p[1*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])
#define velz_p (&vel_p[2*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])
#define velx_p_p (&vel_p_p[0*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])
#define vely_p_p (&vel_p_p[1*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])
#define velz_p_p (&vel_p_p[2*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])

#define DIFF_X(a) (((i==0)?(a[CCTK_GFINDEX3D(cctkGH, i+1, j, k)] - \
                            a[CCTK_GFINDEX3D(cctkGH, i  , j, k)]): \
                    (i==(cctk_lsh[0]-1))?                          \
                           (a[CCTK_GFINDEX3D(cctkGH, i  , j, k)] - \
                            a[CCTK_GFINDEX3D(cctkGH, i-1, j, k)]): \
                       0.5*(a[CCTK_GFINDEX3D(cctkGH, i+1, j, k)] - \
                            a[CCTK_GFINDEX3D(cctkGH, i-1, j, k)]))/\
                   CCTK_DELTA_SPACE(0))
#define DIFF_Y(a) (((j==0)?(a[CCTK_GFINDEX3D(cctkGH, i, j+1, k)] - \
                            a[CCTK_GFINDEX3D(cctkGH, i, j  , k)]): \
                    (j==(cctk_lsh[1]-1))?                          \
                           (a[CCTK_GFINDEX3D(cctkGH, i, j  , k)] - \
                            a[CCTK_GFINDEX3D(cctkGH, i, j-1, k)]): \
                       0.5*(a[CCTK_GFINDEX3D(cctkGH, i, j+1, k)] - \
                            a[CCTK_GFINDEX3D(cctkGH, i, j-1, k)]))/\
                    CCTK_DELTA_SPACE(1))
#define DIFF_Z(a) (((k==0)?(a[CCTK_GFINDEX3D(cctkGH, i, j, k+1)] - \
                            a[CCTK_GFINDEX3D(cctkGH, i, j, k  )]): \
                    (k==(cctk_lsh[2]-1))?                          \
                           (a[CCTK_GFINDEX3D(cctkGH, i, j, k  )] - \
                            a[CCTK_GFINDEX3D(cctkGH, i, j, k-1)]): \
                       0.5*(a[CCTK_GFINDEX3D(cctkGH, i, j, k+1)] - \
                            a[CCTK_GFINDEX3D(cctkGH, i, j, k-1)]))/\
                    CCTK_DELTA_SPACE(2))

//Count lines in file
//Generated from nrpytutorial, with minor edits for NRml_TOV_ET

//Read and Normalize Raw TOV data, and then set to arrays
//Generated from nrpytutorial, with minor edits for use by NRml_TOV_ET
void NRml_Normalize_and_set_data(REAL *restrict r_Schw_arr,
				 REAL *restrict rho_energy_arr,
				 REAL *restrict rho_baryon_arr,
				 REAL *restrict P_arr,
				 REAL *restrict M_arr,
				 REAL *restrict expnu_arr,
				 REAL *restrict exp4phi_arr,
				 REAL *restrict rbar_arr){

  DECLARE_CCTK_PARAMETERS

    //CCTK_VINFO("Reading in raw TOV data...");

  int this_line = 0;
  while (this_line < NRml_Numpoints){
    r_Schw_arr[this_line]     = NRml_Raw_rSchw[this_line];
    rho_energy_arr[this_line] = NRml_Raw_rho_energy[this_line];
    rho_baryon_arr[this_line] = NRml_Raw_rho_baryon[this_line];
    P_arr[this_line]          = NRml_Raw_P[this_line];
    M_arr[this_line]          = NRml_Raw_M[this_line];
    expnu_arr[this_line]      = NRml_Raw_nu[this_line];
    rbar_arr[this_line]       = NRml_Raw_Iso_r[this_line];	
    this_line++;
  }

  //CCTK_VINFO("Raw data Successfully read in!");

  CCTK_INFO("Normalizing raw TOV data..");
  
  //Find the surface of the star
  double R_Schw_surface = -100;
  int Rbar_idx = -100;
  for(int i=1;i<this_line;i++) {
    if(rho_energy_arr[i-1]>0 && rho_energy_arr[i]==0) { R_Schw_surface = r_Schw_arr[i-1]; Rbar_idx = i-1; }
  }

  double M_surface      = M_arr[Rbar_idx];
  double Rbar_surface   = rbar_arr[Rbar_idx];
  double nu_surface     = expnu_arr[Rbar_idx]; //Remember, nu hasn't been adjusted yet, so the expnu_arr just refers to nu itself.
  double normalize      = 0.5*(sqrt(R_Schw_surface*(R_Schw_surface-2.0*M_surface))+R_Schw_surface-M_surface)/Rbar_surface;

  rbar_arr[0] = rbar_arr[0]*normalize;
  rbar_arr[1] = rbar_arr[1]*normalize;
  expnu_arr[0] = exp(expnu_arr[0] - nu_surface + log(1 -2.0*M_surface/R_Schw_surface));
  expnu_arr[1] = exp(expnu_arr[1] - nu_surface + log(1 -2.0*M_surface/R_Schw_surface));
  exp4phi_arr[1] = pow((r_Schw_arr[1]/rbar_arr[1]),2.0);
  exp4phi_arr[0] = exp4phi_arr[1];
  
  for(int i=2;i<this_line;i++){
      rbar_arr[i]    = rbar_arr[i]*normalize;
      expnu_arr[i]   = exp(expnu_arr[i] - nu_surface + log(1 -2.0*M_surface/R_Schw_surface));
      exp4phi_arr[i] = pow((r_Schw_arr[i]/rbar_arr[i]),2.0);
    }

  CCTK_INFO("Normalization of raw data complete!");
  
}

// Find interpolation index using Bisection root-finding algorithm:
// Generated in nrpytutorial, with very minor edits made
static inline int bisection_idx_finder(const REAL rrbar, const int numlines_in_file, const REAL *restrict rbar_arr) {
  int x1 = 0;
  int x2 = numlines_in_file-1;
  REAL y1 = rrbar-rbar_arr[x1];
  REAL y2 = rrbar-rbar_arr[x2];
  if(y1*y2 > 0) {
    CCTK_VINFO("INTERPOLATION BRACKETING ERROR %e | %e %e\n",rrbar,y1,y2);
    CCTK_ERROR("Shutting down due to error");
  }
  for(int i=0;i<numlines_in_file;i++) {
    int x_midpoint = (x1+x2)/2;
    REAL y_midpoint = rrbar-rbar_arr[x_midpoint];
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
//Generated in nrpytutorial, with edits made for compatibility in NRml_TOV_ET
void TOV_interpolate_1D(REAL rrbar,const REAL Rbar,const int Rbar_idx,const int interp_stencil_size,
                        const int numlines_in_file,const REAL *restrict r_Schw_arr,const REAL *restrict rho_energy_arr,const REAL *restrict rho_baryon_arr,const REAL *restrict P_arr,
                        const REAL *restrict M_arr,const REAL *restrict expnu_arr,const REAL *restrict exp4phi_arr,const REAL *restrict rbar_arr,
                        REAL *restrict rho_energy,REAL *restrict rho_baryon,REAL *restrict P,REAL *restrict M,REAL *restrict expnu,REAL *restrict exp4phi) {

  // For this case, we know that for all functions, f(r) = f(-r)
  if(rrbar < 0) rrbar = -rrbar;

  // First find the central interpolation stencil index:
  int idx = bisection_idx_finder(rrbar,numlines_in_file,rbar_arr);


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
  REAL rbar_sample[interp_stencil_size];
  for(int i=idxmin;i<idxmin+interp_stencil_size;i++) {
    rbar_sample[i-idxmin] = rbar_arr[i];
  }
  REAL l_i_of_r[interp_stencil_size];
  for(int i=0;i<interp_stencil_size;i++) {
    REAL numer = 1.0;
    REAL denom = 1.0;
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

  REAL r_Schw = 0.0;
  for(int i=idxmin;i<idxmin+interp_stencil_size;i++) {
    r_Schw      += l_i_of_r[i-idxmin] * r_Schw_arr[i];
    *rho_energy += l_i_of_r[i-idxmin] * rho_energy_arr[i];
    *rho_baryon += l_i_of_r[i-idxmin] * rho_baryon_arr[i];
    *P          += l_i_of_r[i-idxmin] * P_arr[i];
    *M          += l_i_of_r[i-idxmin] * M_arr[i];
    *expnu      += l_i_of_r[i-idxmin] * expnu_arr[i];
    *exp4phi    += l_i_of_r[i-idxmin] * exp4phi_arr[i];
  }

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
void NRml_TOV_Copy(CCTK_INT size, CCTK_REAL *var_p, CCTK_REAL *var)
{
#pragma omp parallel for
    for(int i=0; i<size; i++)
        var_p[i] = var[i];
}

//Interpollate TOV data to grid
//Parts of this function inspired by original TOVsolver in toolkit, namely the populate timelevels portion.
void NRml_interp_Driver(CCTK_ARGUMENTS){

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  
  //Allocate memory for data storage
  REAL *r_Schw_arr =     (REAL *)malloc(sizeof(REAL)*NRml_Numpoints);
  REAL *rho_energy_arr = (REAL *)malloc(sizeof(REAL)*NRml_Numpoints);
  REAL *rho_baryon_arr = (REAL *)malloc(sizeof(REAL)*NRml_Numpoints);
  REAL *P_arr =          (REAL *)malloc(sizeof(REAL)*NRml_Numpoints);
  REAL *M_arr =          (REAL *)malloc(sizeof(REAL)*NRml_Numpoints);
  REAL *expnu_arr =      (REAL *)malloc(sizeof(REAL)*NRml_Numpoints);
  REAL *exp4phi_arr =    (REAL *)malloc(sizeof(REAL)*NRml_Numpoints);
  REAL *rbar_arr =       (REAL *)malloc(sizeof(REAL)*NRml_Numpoints);

  //Normalize and set arrays
  NRml_Normalize_and_set_data(r_Schw_arr,rho_energy_arr,rho_baryon_arr,P_arr,M_arr,expnu_arr,exp4phi_arr,rbar_arr);
  
  //Find Star's Surface
  REAL Rbar = -100;
  int Rbar_idx = -100;
  for(int i=1;i<NRml_Numpoints;i++) {
    if(rho_energy_arr[i-1]>0 && rho_energy_arr[i]==0) { Rbar = rbar_arr[i-1]; Rbar_idx = i-1; }
  }
  
  
  //Now for the actual grid placements. Go over all grid points
  CCTK_INFO("Beginning Grid Placements...");
  for(int i=0; i<cctk_lsh[0]; i++){ //printf("hiti i=%d\n\n", i);
  	for(int j=0; j<cctk_lsh[1]; j++){ //printf("hitj j=%d\n\n", j);
  		for(int k=0; k<cctk_lsh[2]; k++){ //printf("hitk k=%d\n\n", k);
  			int i3d=CCTK_GFINDEX3D(cctkGH,i,j,k); //3D index
  			REAL NRml_rGrid = sqrt((x[i3d]*x[i3d])+(y[i3d]*y[i3d])+(z[i3d]*z[i3d])); //magnitude of r on the grid
  			REAL NRml_rho_energy, NRml_rho_baryon, NRml_P, NRml_M, NRml_expnu, NRml_exp4phi; //TOV quantities
  			TOV_interpolate_1D(NRml_rGrid, Rbar, Rbar_idx,NRml_Interpolation_Stencil,NRml_Numpoints, 
  				r_Schw_arr,rho_energy_arr,rho_baryon_arr,P_arr,M_arr,expnu_arr,exp4phi_arr,rbar_arr,
  				&NRml_rho_energy,&NRml_rho_baryon,&NRml_P,&NRml_M,&NRml_expnu,&NRml_exp4phi);		
  			rho[i3d] = NRml_rho_baryon;
			press[i3d] = NRml_P;
			// Validated that this is the correct assignment, 
			// as TOVSolver's values are comparable.
			// tiny number prevents 0/0.
			eps[i3d] = (NRml_rho_energy / (NRml_rho_baryon+1e-30)) - 1.0;
			if (eps[i3d]<0){eps[i3d]=0.0;}
			betax[i3d] = 0.0;
			betay[i3d] = 0.0;
			betaz[i3d] = 0.0;
			alp[i3d] = pow(NRml_expnu,0.5);
			gxx[i3d] = NRml_exp4phi;
			gyy[i3d] = gxx[i3d];
			gzz[i3d] = gxx[i3d];
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
			velx[i3d] = 0.0;
			vely[i3d] = 0.0;
			velz[i3d] = 0.0;
	        	w_lorentz[i3d] = 1.0;
  		}
  	}
  }
  			
  CCTK_INFO("Grid Placement Successful!");
  
  //This pice of the code was VERY heavily inspired by the original ET TOVSolver, with minor edits for NRml_TOV_ET
  CCTK_INFO("Finalizing Grid...\n");
  int i3d = cctk_lsh[2]*cctk_lsh[1]*cctk_lsh[0];
  switch(NRml_TOV_Populate_Timelevels)
  {
    case 3:
        NRml_TOV_Copy(i3d, gxx_p_p,  gxx);
        NRml_TOV_Copy(i3d, gyy_p_p,  gyy);
        NRml_TOV_Copy(i3d, gzz_p_p,  gzz);
        NRml_TOV_Copy(i3d, gxy_p_p,  gxy);
        NRml_TOV_Copy(i3d, gxz_p_p,  gxz);
        NRml_TOV_Copy(i3d, gyz_p_p,  gyz);
        NRml_TOV_Copy(i3d, rho_p_p,  rho);
        NRml_TOV_Copy(i3d, eps_p_p,  eps);
        NRml_TOV_Copy(i3d, velx_p_p, velx);
        NRml_TOV_Copy(i3d, vely_p_p, vely);
        NRml_TOV_Copy(i3d, velz_p_p, velz);
        NRml_TOV_Copy(i3d, w_lorentz_p_p, w_lorentz);
        // fall through
    case 2:
        NRml_TOV_Copy(i3d, gxx_p,  gxx);
        NRml_TOV_Copy(i3d, gyy_p,  gyy);
        NRml_TOV_Copy(i3d, gzz_p,  gzz);
        NRml_TOV_Copy(i3d, gxy_p,  gxy);
        NRml_TOV_Copy(i3d, gxz_p,  gxz);
        NRml_TOV_Copy(i3d, gyz_p,  gyz);
        NRml_TOV_Copy(i3d, rho_p,  rho);
        NRml_TOV_Copy(i3d, eps_p,  eps);
        NRml_TOV_Copy(i3d, velx_p, velx);
        NRml_TOV_Copy(i3d, vely_p, vely);
        NRml_TOV_Copy(i3d, velz_p, velz);
        NRml_TOV_Copy(i3d, w_lorentz_p, w_lorentz);
        // fall through
    case 1:
        break;
    default:
        CCTK_VWARN(CCTK_WARN_ABORT,
                   "Unsupported number of NRml_TOV_Populate_TimelevelsL: %d",
                   (int)NRml_TOV_Populate_Timelevels);
        break;
  }
  
  //Let's free up the memory
  free(r_Schw_arr);
  free(rho_energy_arr);
  free(rho_baryon_arr);
  free(P_arr);
  free(M_arr);
  free(expnu_arr);
  free(exp4phi_arr);
  free(rbar_arr);
  
  //Finish up. Goodbye!
  CCTK_INFO("Complete! Enjoy your Initial Data!\n\n");
}

//NRml_TOV_ET: By David Boyer
//Code made possible by: OdieGM (ODE solver), NRpy (Generated grid interpollator), GRHayL (EOS library), and the original ET TOVSolver (For reference and a starting point learning the Toolkit)
