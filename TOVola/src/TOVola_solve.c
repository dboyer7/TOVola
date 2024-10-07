#include "GRHayLib.h"
#include "TOVola_interp.c"
#include "TOVola_defines.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_errno.h>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#define ODE_SOLVER_DIM 4
#define TOVOLA_PRESSURE 0
#define TOVOLA_NU 1
#define TOVOLA_MASS 2
#define TOVOLA_R_ISO 3
#define NEGATIVE_R_INTERP_BUFFER 11

/* Prototype for TOVola_interp_function, to be used later */
/*static void TOVola_interp();*/

/* Structure to hold raw TOV data */
typedef struct {
  // EOS type

  int eos_type;

  // Current state variables
  CCTK_REAL rho_baryon;
  CCTK_REAL rho_energy;
  CCTK_REAL r_lengthscale;

  CCTK_REAL *restrict rSchw_arr;
  CCTK_REAL *restrict rho_energy_arr;
  CCTK_REAL *restrict rho_baryon_arr;
  CCTK_REAL *restrict P_arr;
  CCTK_REAL *restrict M_arr;
  CCTK_REAL *restrict nu_arr;
  CCTK_REAL *restrict Iso_r_arr;
  int numels_alloced_TOV_arr;
  ghl_eos_parameters *restrict ghl_eos;
  int numpoints_actually_saved;
  
  //Additional declarations, to pass through ETK parameters in parfile
  CCTK_REAL central_baryon_density;
  CCTK_REAL initial_ode_step_size;
  CCTK_REAL error_limit;
  CCTK_REAL absolute_min_step;
  CCTK_REAL absolute_max_step;

} TOVola_data_struct;

/* Structure to hold TOV data that will become the official ID after normalization */

/* Exception handler to prevent negative pressures */
static void TOVola_exception_handler(CCTK_REAL r, CCTK_REAL y[]) {
  // Ensure pressure does not become negative due to numerical errors
  if (y[TOVOLA_PRESSURE] < 0) {
    y[TOVOLA_PRESSURE] = 0;
  }
}

/* Termination condition for the integration */
static int TOVola_do_we_terminate(CCTK_REAL r, CCTK_REAL y[], TOVola_data_struct *TOVdata) {
  
  if (TOVdata->eos_type == 2) {
  	const CCTK_REAL PMin = exp(TOVdata->ghl_eos->lp_of_lr[0]); //PMin is not zero on the table, so we don't want to exceed table limits
        	if (y[0] <= PMin){
                	return 1;
        	}
  }
  else {
  	if (y[TOVOLA_PRESSURE] <= 0.0) { // For Simple and Piecewise Polytrope
    	return 1;
  	}
  }

  return 0; // Continue integration
}

/* Evaluate rho_baryon and rho_energy based on the EOS type */
static void TOVola_evaluate_rho_and_eps(CCTK_REAL r, const CCTK_REAL y[], TOVola_data_struct *TOVdata) {
  
  // Simple Polytrope
  if (TOVdata->eos_type == 0) {
    CCTK_REAL aK, aGamma;
    CCTK_REAL aRho_baryon = TOVdata->rho_baryon;
    CCTK_REAL eps, aPress;

    // Retrieve K and Gamma from GRHayL
    ghl_hybrid_get_K_and_Gamma(TOVdata->ghl_eos, aRho_baryon, &aK, &aGamma);
    TOVdata->rho_baryon = pow(y[TOVOLA_PRESSURE] / aK, 1.0 / aGamma);
    aRho_baryon = TOVdata->rho_baryon;
    ghl_hybrid_compute_P_cold_and_eps_cold(TOVdata->ghl_eos, aRho_baryon, &aPress, &eps);
    TOVdata->rho_energy = TOVdata->rho_baryon * (1.0 + eps);
  }
  
  // Piecewise Polytrope
  else if (TOVdata->eos_type == 1) {
    //Basically identical to Simple Polytrope, you just have more regions.
    CCTK_REAL aK, aGamma;
    CCTK_REAL aRho_baryon = TOVdata->rho_baryon;
    CCTK_REAL eps, aPress;

    ghl_hybrid_get_K_and_Gamma(TOVdata->ghl_eos,aRho_baryon,&aK,&aGamma);
    TOVdata->rho_baryon = pow(y[TOVOLA_PRESSURE]/aK, 1.0 / aGamma);
    aRho_baryon = TOVdata->rho_baryon;
    ghl_hybrid_compute_P_cold_and_eps_cold(TOVdata->ghl_eos,aRho_baryon,&aPress,&eps);
    TOVdata->rho_energy = TOVdata->rho_baryon*(1.0+eps);
  }
  
  // Tabulated EOS
  else if (TOVdata->eos_type == 2) {
    const CCTK_REAL PMin = exp(TOVdata->ghl_eos->lp_of_lr[0]);
      if(y[TOVOLA_PRESSURE] > PMin){ //Assure you are not exceeding table bounds
        //Use GRHayL function to find our current rho_baryon and rho_energy on the table.
        TOVdata->rho_baryon = ghl_tabulated_compute_rho_from_P(TOVdata->ghl_eos, y[TOVOLA_PRESSURE]);
        CCTK_REAL eps = ghl_tabulated_compute_eps_from_rho(TOVdata->ghl_eos, TOVdata->rho_baryon);
        TOVdata->rho_energy = (TOVdata->rho_baryon)*(1+eps);
        }
      else{
        //Outside the star, densities are zero.
        TOVdata->rho_baryon = 0;
        TOVdata->rho_energy = 0;}
  }
}

/* The main ODE function for GSL */
static int TOVola_ODE(CCTK_REAL r_Schw, const CCTK_REAL y[], CCTK_REAL dydr_Schw[], void *params) {
  // Cast params to TOVdata_struct
  TOVola_data_struct *TOVdata = (TOVola_data_struct *)params;

  // Evaluate rho_baryon and rho_energy based on current state
  TOVola_evaluate_rho_and_eps(r_Schw, y, TOVdata);

  // Dereference the struct to use rho_energy
  CCTK_REAL rho_energy = TOVdata->rho_energy;

  if (isnan(rho_energy)) {
    // Outside the star gives NaNs from the pow function, but we know they
    // should be zeros.
    rho_energy = 0.0;
  }

  // At the center of the star (r_Schw == 0), the TOV equations diverge, so we set reasonable values here.
  if (r_Schw == 0) {
    dydr_Schw[TOVOLA_PRESSURE] = 0.0; // dP/dr
    dydr_Schw[TOVOLA_NU] = 0.0;       // dnu/dr
    dydr_Schw[TOVOLA_MASS] = 0.0;     // dM/dr
    dydr_Schw[TOVOLA_R_ISO] = 1.0;    // dr_iso/dr
    return GSL_SUCCESS;
  }

// TOV Equations

  else {
    double mass_term = (2.0 * y[TOVOLA_MASS]) / r_Schw;
    double denominator = 1.0 - mass_term;
    dydr_Schw[TOVOLA_PRESSURE] = -((rho_energy + y[TOVOLA_PRESSURE]) * (mass_term + 8.0 * M_PI * r_Schw * r_Schw * y[TOVOLA_PRESSURE])) / (2.0*r_Schw * denominator); //dP/dr
    dydr_Schw[TOVOLA_NU] = ((mass_term) + 8.0 * M_PI * r_Schw * r_Schw * y[TOVOLA_PRESSURE]) / (r_Schw * denominator); //dnu/dr
    dydr_Schw[TOVOLA_MASS] = 4.0 * M_PI * r_Schw * r_Schw * rho_energy; // dM/dr
    dydr_Schw[TOVOLA_R_ISO] = y[TOVOLA_R_ISO] / (r_Schw * sqrt(1.0 - mass_term)); // dr_iso/dr (r_iso = isotropic radius, also known as rbar in literature)
  }
	 
  //Adjust Length Scale. 
  if (y[TOVOLA_R_ISO] > 0 && fabs(dydr_Schw[TOVOLA_R_ISO]) > 0) {
    TOVdata->r_lengthscale = fabs(y[TOVOLA_R_ISO] / dydr_Schw[TOVOLA_R_ISO]);
  }

  return GSL_SUCCESS;
}

/* Placeholder Jacobian function required by GSL */
static int TOVola_jacobian_placeholder(CCTK_REAL t, const CCTK_REAL y[], CCTK_REAL *restrict dfdy, CCTK_REAL dfdt[], void *params) {
  // Jacobian is not necessary for the TOV solution, but GSL requires some
  // function Leave it empty as it does not affect the final results
  return GSL_SUCCESS;
}

/* Initialize the ODE variables */
static void TOVola_get_initial_condition(CCTK_REAL y[], TOVola_data_struct *TOVdata) {
  
  // Simple Polytrope
  if (TOVdata->eos_type == 0) {
    CCTK_REAL aK, aGamma;
    CCTK_REAL rhoC_baryon = TOVdata->central_baryon_density;

    // Retrieve K and Gamma from GRHayL
    ghl_hybrid_get_K_and_Gamma(TOVdata->ghl_eos, rhoC_baryon, &aK, &aGamma);
    y[TOVOLA_PRESSURE] = aK * pow(rhoC_baryon, aGamma); // Pressure
    y[TOVOLA_NU] = 0.0;                                 // nu
    y[TOVOLA_MASS] = 0.0;                               // Mass
    y[TOVOLA_R_ISO] = 0.0;                              // r_iso

    // Assign initial conditions
    TOVdata->rho_baryon = rhoC_baryon;
    TOVdata->rho_energy = pow(y[TOVOLA_PRESSURE] / aK, 1.0 / aGamma) + y[TOVOLA_PRESSURE] / (aGamma - 1.0);
    // Pinitial is no longer needed as it's part of y[TOVOLA_PRESSURE]
  }
  
  // Piecewise Polytrope
  else if (TOVdata->eos_type == 1) {
    //Declare EOS info
      CCTK_REAL aK, aGamma;
      CCTK_REAL rhoC_baryon = TOVdata->central_baryon_density;
      CCTK_REAL eps;

      //Use GRHayL to find K and Gamma and calculate initial conditions
      ghl_hybrid_get_K_and_Gamma(TOVdata->ghl_eos, rhoC_baryon, &aK, &aGamma);
      ghl_hybrid_compute_P_cold_and_eps_cold(TOVdata->ghl_eos, rhoC_baryon, &y[TOVOLA_PRESSURE], &eps);
      y[TOVOLA_PRESSURE] = aK*pow(rhoC_baryon, aGamma); //Pressure
      y[TOVOLA_NU] = 0.0; // nu                                                                                                                                               
      y[TOVOLA_MASS] = 0.0; // mass                                                                                                                                                                                           
      y[TOVOLA_R_ISO] = 0.0; // r-bar

      //Assign the initial conditions
      TOVdata->rho_baryon=rhoC_baryon;
      TOVdata->rho_energy = TOVdata->rho_baryon*(1.0+eps);
  }
  
  // Tabulated EOS
  else if (TOVdata->eos_type == 2) {
  //Use GRHayL to find initial pressure on the table
    CCTK_REAL rhoC_baryon = TOVdata->central_baryon_density;
    y[TOVOLA_PRESSURE] = ghl_tabulated_compute_P_from_rho(TOVdata->ghl_eos, rhoC_baryon);
    y[TOVOLA_NU] = 0.0;
    y[TOVOLA_MASS] = 0.0;
    y[TOVOLA_R_ISO] = 0.0;

    //Assign the initial conditions
    TOVdata->rho_baryon=rhoC_baryon;
    CCTK_REAL eps = ghl_tabulated_compute_eps_from_rho(TOVdata->ghl_eos,rhoC_baryon);
    TOVdata->rho_energy = (TOVdata->rho_baryon)*(1.0+eps);
  }

  CCTK_VINFO("Initial Conditions Set: P = %.6e, nu = %.6e, M = %.6e, r_iso = %.6e", y[TOVOLA_PRESSURE], y[TOVOLA_NU], y[TOVOLA_MASS], y[TOVOLA_R_ISO]);
}

/* Assign constants after each integration step */
static void TOVola_assign_constants(CCTK_REAL c[], TOVola_data_struct *TOVdata) {
  // Assign the densities
  c[0] = TOVdata->rho_energy; // Total energy density
  c[1] = TOVdata->rho_baryon; // Baryon density

  // Handle NaN cases
  if (isnan(TOVdata->rho_energy)) {
    c[0] = 0.0;
  }
}

/* Function to set up the GSL ODE system and driver */
static int setup_ode_system(const char *ode_method, gsl_odeiv2_system *system, gsl_odeiv2_driver **driver, TOVola_data_struct *TOVdata) {
  

  system->function = TOVola_ODE;
  system->jacobian = TOVola_jacobian_placeholder;
  system->dimension = 4; // Hardcoded as per requirements
  system->params = TOVdata;

  if (CCTK_EQUALS(ode_method, "ARKF")) {
    *driver = gsl_odeiv2_driver_alloc_y_new(system, gsl_odeiv2_step_rkf45, TOVdata->initial_ode_step_size, TOVdata->error_limit,
                                            TOVdata->error_limit);
  } else if (CCTK_EQUALS(ode_method, "ADP8")) {
    *driver = gsl_odeiv2_driver_alloc_y_new(system, gsl_odeiv2_step_rk8pd, TOVdata->initial_ode_step_size, TOVdata->error_limit,
                                            TOVdata->error_limit);
  } else {
    CCTK_ERROR("Invalid ODE method. Use 'ARKF' or 'ADP8'.");
    return -1;
  }

  if (*driver == NULL) {
    CCTK_ERROR("Failed to allocate GSL ODE driver.");
    return -1;
  }

  /* Set minimum and maximum step sizes */
  gsl_odeiv2_driver_set_hmin(*driver, TOVdata->absolute_min_step);
  gsl_odeiv2_driver_set_hmax(*driver, TOVdata->absolute_max_step);

  return 0;
}

/* Initialize TOVola_data_struct structure with initial allocation */
static int initialize_tovola_data(TOVola_data_struct *TOVdata) {
  TOVdata->rSchw_arr = (CCTK_REAL *restrict)malloc(sizeof(CCTK_REAL) * TOVdata->numels_alloced_TOV_arr);
  TOVdata->rho_energy_arr = (CCTK_REAL *restrict)malloc(sizeof(CCTK_REAL) * TOVdata->numels_alloced_TOV_arr);
  TOVdata->rho_baryon_arr = (CCTK_REAL *restrict)malloc(sizeof(CCTK_REAL) * TOVdata->numels_alloced_TOV_arr);
  TOVdata->P_arr = (CCTK_REAL *restrict)malloc(sizeof(CCTK_REAL) * TOVdata->numels_alloced_TOV_arr);
  TOVdata->M_arr = (CCTK_REAL *restrict)malloc(sizeof(CCTK_REAL) * TOVdata->numels_alloced_TOV_arr);
  TOVdata->nu_arr = (CCTK_REAL *restrict)malloc(sizeof(CCTK_REAL) * TOVdata->numels_alloced_TOV_arr);
  TOVdata->Iso_r_arr = (CCTK_REAL *restrict)malloc(sizeof(CCTK_REAL) * TOVdata->numels_alloced_TOV_arr);
  if (!TOVdata->rSchw_arr || !TOVdata->rho_energy_arr || !TOVdata->rho_baryon_arr || !TOVdata->P_arr || !TOVdata->M_arr || !TOVdata->nu_arr ||
      !TOVdata->Iso_r_arr) {
    CCTK_ERROR("Memory allocation failed for TOVola_data_struct.");
    return -1;
  }
  return 0;
}

/* Free TOVola_data_struct structure */
static void free_tovola_data(TOVola_data_struct *TOVdata) {
  free(TOVdata->rSchw_arr);
  free(TOVdata->rho_energy_arr);
  free(TOVdata->rho_baryon_arr);
  free(TOVdata->P_arr);
  free(TOVdata->M_arr);
  free(TOVdata->nu_arr);
  free(TOVdata->Iso_r_arr);
  TOVdata->numels_alloced_TOV_arr = 0;
}


/* Normalize and set data */
static void TOVola_Normalize_and_set_data_integrated(TOVola_data_struct *TOVdata, CCTK_REAL *restrict r_Schw, CCTK_REAL *restrict rho_energy,
                                              CCTK_REAL *restrict rho_baryon, CCTK_REAL *restrict P, CCTK_REAL *restrict M, CCTK_REAL *restrict expnu,
                                              CCTK_REAL *restrict exp4phi, CCTK_REAL *restrict r_iso) {
  
	
  CCTK_INFO("TOVola Normalizing raw TOV data...");

  /* Check if there are enough points to normalize */
  if (TOVdata->numpoints_actually_saved < 2) {
    CCTK_ERROR("Not enough data points to normalize.");
  }

  /* Copy raw data to normalized arrays */
  for (int i = 0; i < TOVdata->numpoints_actually_saved; i++) {
    r_Schw[i] = TOVdata->rSchw_arr[i];
    rho_energy[i] = TOVdata->rho_energy_arr[i];
    rho_baryon[i] = TOVdata->rho_baryon_arr[i];
    P[i] = TOVdata->P_arr[i];
    M[i] = TOVdata->M_arr[i];
    expnu[i] = TOVdata->nu_arr[i];
    r_iso[i] = TOVdata->Iso_r_arr[i];
  }

  /* Surface values for normalization */
  const CCTK_REAL R_Schw_surface = r_Schw[TOVdata->numpoints_actually_saved - 1];
  const CCTK_REAL M_surface = M[TOVdata->numpoints_actually_saved - 1];
  const CCTK_REAL r_iso_surface = r_iso[TOVdata->numpoints_actually_saved - 1];
  const CCTK_REAL nu_surface = expnu[TOVdata->numpoints_actually_saved - 1];

  const CCTK_REAL normalize = 0.5 * (sqrt(R_Schw_surface * (R_Schw_surface - 2.0 * M_surface)) + R_Schw_surface - M_surface) / r_iso_surface;

  /* Normalize r_iso and calculate expnu and exp4phi */
  for (int i = 0; i < TOVdata->numpoints_actually_saved; i++) {
    r_iso[i] *= normalize;
    expnu[i] = exp(expnu[i] - nu_surface + log(1.0 - 2.0 * M_surface / R_Schw_surface));
    exp4phi[i] = (r_Schw[i] / r_iso[i]) * (r_Schw[i] / r_iso[i]);
  }
  CCTK_INFO("Normalization of raw data complete!");
}

/* Extend data to r<0, to ensure we can interpolate to r=0 */
void extend_to_negative_r(CCTK_REAL *restrict arr, const CCTK_REAL parity, CCTK_REAL *restrict tmp, const TOVola_data_struct *restrict TOVdata) {
  for(int i=0;i<NEGATIVE_R_INTERP_BUFFER; i++) tmp[i] = parity * arr[NEGATIVE_R_INTERP_BUFFER - i - 1];
  for(int i=0;i<TOVdata->numpoints_actually_saved; i++) tmp[i+NEGATIVE_R_INTERP_BUFFER] = arr[i];
  memcpy(arr, tmp, sizeof(CCTK_REAL) * (TOVdata->numpoints_actually_saved+NEGATIVE_R_INTERP_BUFFER));
}

//For timelevel population, from original TOVsolver in the toolkit.
void TOVola_TOV_Copy(CCTK_INT size, CCTK_REAL *var_p, CCTK_REAL *var)
{
#pragma omp parallel for
    for(int i=0; i<size; i++)
        var_p[i] = var[i];
}


//Helpful defines for later.
#define velx (&vel[0*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])
#define vely (&vel[1*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])
#define velz (&vel[2*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])
//For timelevel population
#define velx_p (&vel_p[0*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])
#define vely_p (&vel_p[1*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])
#define velz_p (&vel_p[2*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])
#define velx_p_p (&vel_p_p[0*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])
#define vely_p_p (&vel_p_p[1*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])
#define velz_p_p (&vel_p_p[2*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])



//Perform the TOV integration using GSL
void TOVola_Solve_and_Interp(CCTK_ARGUMENTS){

  DECLARE_CCTK_PARAMETERS;
  DECLARE_CCTK_ARGUMENTS;
  TOVola_ID_persist_struct TOVola_ID_persist_tmp;  // allocates memory for the pointer below.
  TOVola_ID_persist_struct *restrict TOVola_ID_persist = &TOVola_ID_persist_tmp;
  CCTK_REAL current_position = 0;

  /* Set up ODE system and driver */
  TOVola_data_struct TOVdata_tmp; // allocates memory for the pointer below.
  TOVola_data_struct *restrict TOVdata = &TOVdata_tmp;
  gsl_odeiv2_system system;
  gsl_odeiv2_driver *driver;

  //Checking and setting EOS
  if(CCTK_EQUALS("Simple",TOVola_EOS_type)){
    CCTK_INFO("Simple Polytrope");
    TOVdata->eos_type = 0;
    if(ghl_eos->neos!=1){
      CCTK_INFO("Error: Too many regions for the simple polytrope.");
      CCTK_INFO("Check your value for neos, or use a piecewise polytrope");
      CCTK_ERROR("Shutting down due to error...");}       
    }
  else if(CCTK_EQUALS("Piecewise",TOVola_EOS_type)){
    CCTK_INFO("Piecewise Polytrope");
    TOVdata->eos_type=1;}
  else if(CCTK_EQUALS("Tabulated",TOVola_EOS_type)){
    CCTK_INFO("Tabulated EOS");
    TOVdata->eos_type=2;
    ghl_tabulated_compute_Ye_P_eps_of_rho_beq_constant_T(TOVola_Tin, ghl_eos);
  }
  else{
    CCTK_INFO("ERROR: Invalid EOS type. Must be either 'Simple', 'Piecewise', or 'Tabulated'");
    CCTK_ERROR("Shutting down due to error...");}

  TOVdata->numpoints_actually_saved = 0;
  TOVdata->error_limit = TOVola_error_limit;
  TOVdata->initial_ode_step_size = TOVola_initial_ode_step_size;
  TOVdata->absolute_max_step = TOVola_absolute_max_step;
  TOVdata->absolute_min_step = TOVola_absolute_min_step;
  TOVdata->central_baryon_density = TOVola_central_baryon_density;
  TOVdata->ghl_eos = ghl_eos;
  if (setup_ode_system(TOVola_ODE_method, &system, &driver, TOVdata) != 0) {
    CCTK_ERROR("Failed to set up ODE system.");
  }
  
  CCTK_INFO("Starting TOV Integration using GSL for TOVola...");
  /* Initialize ODE variables */
  CCTK_REAL TOVola_eq[ODE_SOLVER_DIM];
  CCTK_REAL c[2];
  TOVola_get_initial_condition(TOVola_eq, TOVdata);
  TOVola_assign_constants(c, TOVdata);

  /* Initial memory allocation */
  TOVdata->numels_alloced_TOV_arr = 1024;
  if (initialize_tovola_data(TOVdata) != 0) {
    gsl_odeiv2_driver_free(driver);
    CCTK_ERROR("Failed to initialize TOVola_data_struct.");
  }

  /* Integration loop */
  TOVdata->r_lengthscale = TOVola_initial_ode_step_size; // initialize dr to a crazy small value in double precision.
  for (int i = 0; i < TOVola_size; i++) {
    CCTK_REAL dr = 0.01 * TOVdata->r_lengthscale;
    if (TOVdata->rho_baryon < 0.05 * TOVola_central_baryon_density) {
      // To get a super-accurate mass, reduce the dr sampling near the surface of the star.
      dr = 1e-6 * TOVdata->r_lengthscale;
    }
    /* Exception handling */
    TOVola_exception_handler(current_position, TOVola_eq);

    /* Apply ODE step */
    int status = gsl_odeiv2_driver_apply(driver, &current_position, current_position + dr, TOVola_eq);
    if (status != GSL_SUCCESS) {
      CCTK_VINFO("GSL ODE solver failed with status %d.", status);
      gsl_odeiv2_driver_free(driver);
      CCTK_ERROR("Shutting down due to error");
    };

    /* Post-step exception handling */
    TOVola_exception_handler(current_position, TOVola_eq);

    /* Evaluate densities */
    TOVola_evaluate_rho_and_eps(current_position, TOVola_eq, TOVdata);
    TOVola_assign_constants(c, TOVdata);

    /* Check if reallocation is needed */
    if (TOVdata->numpoints_actually_saved >= TOVdata->numels_alloced_TOV_arr) {
      // Update arr_size instead of modifying the macro
      const int new_arr_size = 1.5 * TOVdata->numels_alloced_TOV_arr;
      TOVdata->numels_alloced_TOV_arr = new_arr_size;
      TOVdata->rSchw_arr = realloc(TOVdata->rSchw_arr, sizeof(CCTK_REAL) * new_arr_size);
      TOVdata->rho_energy_arr = realloc(TOVdata->rho_energy_arr, sizeof(CCTK_REAL) * new_arr_size);
      TOVdata->rho_baryon_arr = realloc(TOVdata->rho_baryon_arr, sizeof(CCTK_REAL) * new_arr_size);
      TOVdata->P_arr = realloc(TOVdata->P_arr, sizeof(CCTK_REAL) * new_arr_size);
      TOVdata->M_arr = realloc(TOVdata->M_arr, sizeof(CCTK_REAL) * new_arr_size);
      TOVdata->nu_arr = realloc(TOVdata->nu_arr, sizeof(CCTK_REAL) * new_arr_size);
      TOVdata->Iso_r_arr = realloc(TOVdata->Iso_r_arr, sizeof(CCTK_REAL) * new_arr_size);

      if (!TOVdata->rSchw_arr || !TOVdata->rho_energy_arr || !TOVdata->rho_baryon_arr || !TOVdata->P_arr || !TOVdata->M_arr || !TOVdata->nu_arr ||
          !TOVdata->Iso_r_arr) {
        CCTK_ERROR("Memory reallocation failed during integration.\n");
        gsl_odeiv2_driver_free(driver);
      }
    }

    /* Store data */
    TOVdata->rSchw_arr[TOVdata->numpoints_actually_saved] = current_position;
    TOVdata->rho_energy_arr[TOVdata->numpoints_actually_saved] = c[0];
    TOVdata->rho_baryon_arr[TOVdata->numpoints_actually_saved] = c[1];
    TOVdata->P_arr[TOVdata->numpoints_actually_saved] = TOVola_eq[TOVOLA_PRESSURE];
    TOVdata->M_arr[TOVdata->numpoints_actually_saved] = TOVola_eq[TOVOLA_MASS];
    TOVdata->nu_arr[TOVdata->numpoints_actually_saved] = TOVola_eq[TOVOLA_NU];
    TOVdata->Iso_r_arr[TOVdata->numpoints_actually_saved] = TOVola_eq[TOVOLA_R_ISO];
    TOVdata->numpoints_actually_saved++;

    /* Termination condition */
    if (TOVola_do_we_terminate(current_position, TOVola_eq, TOVdata)) {
      CCTK_VINFO("Finished Integration at position %.6e with Mass %.14e", current_position, TOVola_eq[TOVOLA_MASS]);
      break;
    }
  }

  /* Cleanup */
  gsl_odeiv2_driver_free(driver);
  CCTK_INFO("ODE Solver using GSL for TOVola Shutting Down...");

  // Data in TOVdata->*_arr are stored at r=TOVdata->commondata->initial_ode_step_size > 0 up to the stellar surface.
  // However, we may need data at r=0, which would require extrapolation.
  // To prevent that, we copy INTERP_BUFFER data points from r>0 to r<0 so that we can always interpolate.
  CCTK_REAL *restrict tmp = malloc(sizeof(CCTK_REAL) * (TOVdata->numpoints_actually_saved + NEGATIVE_R_INTERP_BUFFER));

  //printf("Rbefor = %.15e\n", TOVdata->r_Schw_arr[TOVdata->numpoints_actually_saved-1]);

  extend_to_negative_r(TOVdata->rSchw_arr, -1.0, tmp, TOVdata);
  extend_to_negative_r(TOVdata->rho_energy_arr, +1.0, tmp, TOVdata);
  extend_to_negative_r(TOVdata->rho_baryon_arr, +1.0, tmp, TOVdata);
  extend_to_negative_r(TOVdata->P_arr, +1.0, tmp, TOVdata);
  extend_to_negative_r(TOVdata->M_arr, +1.0, tmp, TOVdata);
  extend_to_negative_r(TOVdata->nu_arr, +1.0, tmp, TOVdata);
  extend_to_negative_r(TOVdata->Iso_r_arr, -1.0, tmp, TOVdata);

  free(tmp);
  TOVdata->numpoints_actually_saved += NEGATIVE_R_INTERP_BUFFER;

  //printf("Rafter = %.15e\n", TOVdata->r_Schw_arr[TOVdata->numpoints_actually_saved-1]);
  //for(int i=0;i<TOVdata->numpoints_actually_saved; i++) printf("%e %e %e iiii\n", TOVdata->r_Schw_arr[i], TOVdata->r_iso_arr[i], TOVdata->rho_energy_arr[i]);

  /* Allocate and populate TOVola_ID_persist_struct arrays */
  TOVola_ID_persist->r_Schw_arr = (CCTK_REAL *restrict)malloc(sizeof(CCTK_REAL) * TOVdata->numpoints_actually_saved);
  TOVola_ID_persist->rho_energy_arr = (CCTK_REAL *restrict)malloc(sizeof(CCTK_REAL) * TOVdata->numpoints_actually_saved);
  TOVola_ID_persist->rho_baryon_arr = (CCTK_REAL *restrict)malloc(sizeof(CCTK_REAL) * TOVdata->numpoints_actually_saved);
  TOVola_ID_persist->P_arr = (CCTK_REAL *restrict)malloc(sizeof(CCTK_REAL) * TOVdata->numpoints_actually_saved);
  TOVola_ID_persist->M_arr = (CCTK_REAL *restrict)malloc(sizeof(CCTK_REAL) * TOVdata->numpoints_actually_saved);
  TOVola_ID_persist->expnu_arr = (CCTK_REAL *restrict)malloc(sizeof(CCTK_REAL) * TOVdata->numpoints_actually_saved);
  TOVola_ID_persist->exp4phi_arr = (CCTK_REAL *restrict)malloc(sizeof(CCTK_REAL) * TOVdata->numpoints_actually_saved);
  TOVola_ID_persist->r_iso_arr = (CCTK_REAL *restrict)malloc(sizeof(CCTK_REAL) * TOVdata->numpoints_actually_saved);

  if (!TOVola_ID_persist->r_Schw_arr || !TOVola_ID_persist->rho_energy_arr || !TOVola_ID_persist->rho_baryon_arr || !TOVola_ID_persist->P_arr || !TOVola_ID_persist->M_arr ||
      !TOVola_ID_persist->expnu_arr || !TOVola_ID_persist->exp4phi_arr || !TOVola_ID_persist->r_iso_arr) {
    free_tovola_data(TOVdata);
    CCTK_ERROR("Memory allocation failed for TOVola_ID_persist_struct arrays.\n");
  }

  TOVola_ID_persist->numpoints_arr = TOVdata->numpoints_actually_saved;

  /* Normalize and set data */
  TOVola_Normalize_and_set_data_integrated(TOVdata, TOVola_ID_persist->r_Schw_arr, TOVola_ID_persist->rho_energy_arr, TOVola_ID_persist->rho_baryon_arr, TOVola_ID_persist->P_arr,
                                           TOVola_ID_persist->M_arr, TOVola_ID_persist->expnu_arr, TOVola_ID_persist->exp4phi_arr, TOVola_ID_persist->r_iso_arr);

  /* Free raw data as it's no longer needed */
  free_tovola_data(TOVdata);

  /* Now to interp, and finalize the grid. */
  CCTK_REAL TOVola_Rbar = TOVola_ID_persist->r_iso_arr[TOVola_ID_persist->numpoints_arr-1]; 
  CCTK_REAL TOVola_Mass = TOVola_ID_persist->M_arr[TOVola_ID_persist->numpoints_arr-1];
  
  
  //Now for the actual grid placements. Go over all grid points
  CCTK_INFO("TOVola Beginning Grid Placements...");
  for(int i=0; i<cctk_lsh[0]; i++){ 
  	for(int j=0; j<cctk_lsh[1]; j++){ 
  		for(int k=0; k<cctk_lsh[2]; k++){ 
  			int i3d=CCTK_GFINDEX3D(cctkGH,i,j,k); //3D index
  			CCTK_REAL TOVola_r_iso = sqrt((x[i3d]*x[i3d])+(y[i3d]*y[i3d])+(z[i3d]*z[i3d])); //magnitude of r on the grid
  			CCTK_REAL TOVola_rho_energy, TOVola_rho_baryon, TOVola_P, TOVola_M, TOVola_expnu, TOVola_exp4phi; //Declare TOV quantities
  			if (TOVola_r_iso < TOVola_Rbar){ //If we are INSIDE the star, we need to interpollate the data to the grid.
  				TOVola_TOV_interpolate_1D(TOVola_r_iso, TOVola_Interpolation_Stencil, TOVola_Max_Interpolation_Stencil,
                                      TOVola_ID_persist->numpoints_arr, TOVola_ID_persist->r_Schw_arr,
                                      TOVola_ID_persist->rho_energy_arr, TOVola_ID_persist->rho_baryon_arr, TOVola_ID_persist->P_arr,
                                      TOVola_ID_persist->M_arr, TOVola_ID_persist->expnu_arr, TOVola_ID_persist->exp4phi_arr,
                                      TOVola_ID_persist->r_iso_arr, &TOVola_rho_energy, &TOVola_rho_baryon, &TOVola_P,
                                      &TOVola_M, &TOVola_expnu, &TOVola_exp4phi);
  				rho[i3d] = TOVola_rho_baryon;
				press[i3d] = TOVola_P;
				// tiny number prevents 0/0.
				eps[i3d] = (TOVola_rho_energy / (TOVola_rho_baryon+1e-30)) - 1.0;
				if (eps[i3d]<0){eps[i3d]=0.0;}
				alp[i3d] = pow(TOVola_expnu,0.5);//This is the lapse
				gxx[i3d] = TOVola_exp4phi;//This is the values for the metric in the coordinates we chose.
				gyy[i3d] = gxx[i3d];
				gzz[i3d] = gxx[i3d];}
			else { //If we are OUTSIDE the star, we need to calculate the grid functions directly. Thank you, Schwarzchild.
				CCTK_REAL TOVola_rSchw_outside = (TOVola_r_iso+TOVola_Mass) + TOVola_Mass*TOVola_Mass/(4.0*TOVola_r_iso);//Need to know what rSchw is at our current grid location.
				rho[i3d] = 0.0;
				press[i3d] = 0.0;
				eps[i3d] = 0.0;
				alp[i3d] = pow(1-2*TOVola_Mass/TOVola_rSchw_outside,0.5); //Goes to Schwarschild
				gxx[i3d] = pow((TOVola_rSchw_outside/TOVola_r_iso),2.0);
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
  CCTK_INFO("Populating Time Levels...");

  //This is to populate time levels.
  //Luckily, this is a static solution, so the logic isn't too complicated.
  //From the original TOVsolver in the toolkit.
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


  CCTK_INFO("Population Complete!");
  
  CCTK_INFO("Complete! Enjoy your initial data!");
  CCTK_INFO("TOVola shutting down...");
}
