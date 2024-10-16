#include "cctk.h"

#include "gsl/gsl_odeiv2.h"
#include "gsl/gsl_errno.h"

#include "GRHayLib.h"
#include "TOVola_defines.h"

#pragma once

/*********************************************************************************************************************************************************
//This is the header file that contains all the information about the TOVs and integration schemes. Used in conjunction with GSL in the driver function.
*********************************************************************************************************************************************************/

#define ODE_SOLVER_DIM 4
#define TOVOLA_PRESSURE 0
#define TOVOLA_NU 1
#define TOVOLA_MASS 2
#define TOVOLA_R_ISO 3
#define NEGATIVE_R_INTERP_BUFFER 11

/* Structure to hold raw TOV data. Only used by these functions, so is declared in solve.h */
typedef struct {
  // EOS type

  CCTK_INT eos_type;

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
  CCTK_INT numels_alloced_TOV_arr;
  ghl_eos_parameters *restrict ghl_eos;
  CCTK_INT numpoints_actually_saved;
  
  //Additional declarations, to pass through ETK parameters in parfile without causing namespace pollution.
  CCTK_REAL central_baryon_density;
  CCTK_REAL initial_ode_step_size;
  CCTK_REAL error_limit;
  CCTK_REAL absolute_min_step;
  CCTK_REAL absolute_max_step;

} TOVola_data_struct;

/* Exception handler to prevent negative pressures */
static void TOVola_exception_handler(CCTK_REAL r, CCTK_REAL y[]) {
  // Ensure pressure does not become negative due to numerical errors
  if (y[TOVOLA_PRESSURE] < 0) {
    y[TOVOLA_PRESSURE] = 0;
  }
}

/* Termination condition for the integration */
static CCTK_INT TOVola_do_we_terminate(CCTK_REAL r, CCTK_REAL y[], TOVola_data_struct *TOVdata) {
  
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
static CCTK_INT TOVola_ODE(CCTK_REAL r_Schw, const CCTK_REAL y[], CCTK_REAL dydr_Schw[], void *params) {
  // Cast params to TOVdata_struct
  TOVola_data_struct *TOVdata = (TOVola_data_struct *)params;

  // Evaluate rho_baryon and rho_energy based on current state
  TOVola_evaluate_rho_and_eps(r_Schw, y, TOVdata);

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
static CCTK_INT TOVola_jacobian_placeholder(CCTK_REAL t, const CCTK_REAL y[], CCTK_REAL *restrict dfdy, CCTK_REAL dfdt[], void *params) {
  // Jacobian is not necessary for the TOV solution, but GSL requires some
  // function. Leave it empty as it does not affect the final results
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
static CCTK_INT setup_ode_system(const char *ode_method, gsl_odeiv2_system *system, gsl_odeiv2_driver **driver, TOVola_data_struct *TOVdata) {
  

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
static CCTK_INT initialize_tovola_data(TOVola_data_struct *TOVdata) {
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

/* Initialize TOVola_ID_persist_struct with initial allocation */
static CCTK_INT initialize_ID_persist_data(TOVola_ID_persist_struct *TOVola_ID_persist, TOVola_data_struct *TOVdata) {
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
    return -1;
  }
  return 0;
}

/* Free TOVola_ID_persist_struct */
static void free_ID_persist_data(TOVola_ID_persist_struct *TOVola_ID_persist) {
  free(TOVola_ID_persist->r_Schw_arr);
  free(TOVola_ID_persist->rho_energy_arr);
  free(TOVola_ID_persist->rho_baryon_arr);
  free(TOVola_ID_persist->P_arr);
  free(TOVola_ID_persist->M_arr);
  free(TOVola_ID_persist->expnu_arr);
  free(TOVola_ID_persist->r_iso_arr);
  free(TOVola_ID_persist->exp4phi_arr);
}


/* Normalize and set data */
static void TOVola_Normalize_and_set_data_integrated(TOVola_data_struct *TOVdata, CCTK_REAL *restrict r_Schw, CCTK_REAL *restrict rho_energy,
                                              CCTK_REAL *restrict rho_baryon, CCTK_REAL *restrict P, CCTK_REAL *restrict M, CCTK_REAL *restrict expnu,
                                              CCTK_REAL *restrict exp4phi, CCTK_REAL *restrict r_iso) {
  
	
  /* Check if there are enough points to normalize */
  if (TOVdata->numpoints_actually_saved < 2) {
    CCTK_ERROR("Not enough data points to normalize.");
  }

  /* Copy raw data to normalized arrays */
  for (CCTK_INT i = 0; i < TOVdata->numpoints_actually_saved; i++) {
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
  for (CCTK_INT i = 0; i < TOVdata->numpoints_actually_saved; i++) {
    r_iso[i] *= normalize;
    expnu[i] = exp(expnu[i] - nu_surface + log(1.0 - 2.0 * M_surface / R_Schw_surface));
    exp4phi[i] = (r_Schw[i] / r_iso[i]) * (r_Schw[i] / r_iso[i]);
  }
}

/* Extend data to r<0, to ensure we can interpolate to r=0 */
static void extend_to_negative_r(CCTK_REAL *restrict arr, const CCTK_REAL parity, CCTK_REAL *restrict tmp, const TOVola_data_struct *restrict TOVdata) {
  for(CCTK_INT i=0;i<NEGATIVE_R_INTERP_BUFFER; i++) tmp[i] = parity * arr[NEGATIVE_R_INTERP_BUFFER - i - 1];
  for(CCTK_INT i=0;i<TOVdata->numpoints_actually_saved; i++) tmp[i+NEGATIVE_R_INTERP_BUFFER] = arr[i];
  memcpy(arr, tmp, sizeof(CCTK_REAL) * (TOVdata->numpoints_actually_saved+NEGATIVE_R_INTERP_BUFFER));
}
