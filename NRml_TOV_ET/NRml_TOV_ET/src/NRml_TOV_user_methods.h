#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include "GRHayLib.h"

//ETK interface. 
#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

// UPDATE: This is a special user methods file edited from Odie that defines the TOV equations themselves.
// For use in NRml_TOV.
// DO NOT TOUCH THIS FILE!!! All edits that need to be made are in the NRml_TOV_Driver.c file.


struct constant_parameters { 
    int dimension;
    double rho_energy;
    double rho_baryon;
    double Pinitial;
    int neos;
    double* Gamma;
    double* rhoBound;
    double rhoCentral_baryon;
    char type;
    ghl_eos_parameters ghl_eos;
};

// Here are the prototypes for the functions in this file, stated explicitly for the sake of clarity. 
void exception_handler (double x, double y[]); 
// Handles any exceptions the user may wish to define.

int do_we_terminate (double x, double y[], struct constant_parameters *params); 
// User-defined endpoint.
// Generally used if the code won't terminate itself from outside, or if there's a variable condition.

void const_eval (double x, const double y[], struct constant_parameters *params);
// Assign constants to the constant_parameters struct based on values in y[].

int diffy_Q_eval (double x, double y[], double dydx[], void *params);
// The definition for the system of equations itself goes here.

int known_Q_eval (double x, double y[]);
// If an exact solution is known, it goes here, otherwise leave empty.

void get_initial_condition (double y[], struct constant_parameters*params);
// Initial conditions for the system of differential equations.

void assign_constants (double c[], struct constant_parameters *params);
// Used to read values from constant_parameters into an array so they can be reported in sequence. 



void exception_handler (double x, double y[])
{
    if (y[0] < 0) {
        y[0] = 0;
    } 
}

int do_we_terminate (double x, double y[], struct constant_parameters *params)
{
  DECLARE_CCTK_PARAMETERS
  //ET Termination condition: Keep Calculating to an outer limit, for grid interpollation purposes.
  if (x>NRml_outer_radius_limit){ return 1;}
  
  //Standalone Termination Condition: Stop at star's surface
  //if (y[0] < 1e-30) {
    //Asymptotic termination condition
  //      return 1;}
  //  }
  //else if (isnan(params->rho_energy) == 1){
    //The TOV is undefined outside the surface of the star.
    //This is a kill condition, so we don't continue on.
  //  return 1;}

  else { return 0;}
    // return 1; for termination.
}

void const_eval (double x, const double y[], struct constant_parameters *params)
{

  DECLARE_CCTK_PARAMETERS
  //IF SIMPLE POLYTROPE
  if (strcmp("Simple",NRml_EOS_type) == 0){
    double aK;
    double aGamma;
    double aRho_baryon = params->rho_baryon;
    ghl_hybrid_get_K_and_Gamma(ghl_eos,aRho_baryon,&aK,&aGamma);
    params->rho_energy = pow(y[0] / aK , 1.0 / aGamma) + y[0] / (aGamma - 1.0);
    params->rho_baryon = pow(y[0]/aK, 1.0 / aGamma);}

  //IF PIECEWISE POLYTROPE
  else if (strcmp("Piecewise",NRml_EOS_type) == 0){
    double aK;
    double aGamma;
    double aRho_baryon = params->rho_baryon;
    double eps;
    double aPress;
    ghl_hybrid_get_K_and_Gamma(ghl_eos,aRho_baryon,&aK,&aGamma);
    params->rho_baryon = pow(y[0]/aK, 1.0 / aGamma);
    aRho_baryon = params->rho_baryon;
    ghl_hybrid_compute_P_cold_and_eps_cold(ghl_eos,aRho_baryon,&aPress,&eps);
    params->rho_energy = params->rho_baryon*(1.0+eps);
  }

  //IF TABULATED EOS
  else if (strcmp("Tabulated",NRml_EOS_type) == 0){
    if(y[0] > 0.0){
      params->rho_baryon = ghl_tabulated_compute_rho_from_P(ghl_eos, y[0]);
      double eps = ghl_tabulated_compute_eps_from_rho(ghl_eos,params->rho_baryon);
      params->rho_energy = (params->rho_baryon)*(1+eps);
      }
    else{
      //Outside the star
      params->rho_baryon = 0;
      params->rho_energy = 0;}
  }
    
}

int diffy_Q_eval (double x, double y[], double dydx[], void *params)
{

    exception_handler(x,y);
    const_eval(x,y,params);

    // Dereference the struct
    double rho_energy = (*(struct constant_parameters*)params).rho_energy;

    if (isnan(rho_energy)==1){
      //Outside the star gives nans from the pow function, but we know they should be zeros.
      rho_energy=0.0;}

    //AT THE CENTER
    if(x == 0.0) {
        dydx[0] = 0.0;
        dydx[1] = 0.0;
        dydx[2] = 0.0;
        dydx[3] = 1.0;
    }
    //TOV EQUATIONS
    else {
      dydx[0] = -((rho_energy+y[0])*( (2.0*y[2])/(x) + 8.0*M_PI*x*x*y[0] ))/(x*2.0*(1.0 - (2.0*y[2])/(x)));//Pressure
      dydx[1] =  ((2.0*y[2])/(x) + 8.0*M_PI*x*x*y[0])/(x*(1.0 - (2.0*y[2])/(x)));//nu
      dydx[2] = 4.0*M_PI*x*x*rho_energy;//mass
      dydx[3] = (y[3])/(x*sqrt(1.0-(2.0*y[2])/x));//Iso_r
    }

    return 0;
}

int known_Q_eval (double x, double y[])
{
  //No Analytic solution for general TOV, leave empty
    return 1;
}

void get_initial_condition (double y[], struct constant_parameters *params)
{
  // be sure to have these MATCH the equations in diffy_Q_eval
  DECLARE_CCTK_PARAMETERS
  //IF SIMPLE POLYTROPE
  if (strcmp("Simple",NRml_EOS_type) == 0){
    double aK;
    double aGamma;
    double rhoC_baryon = params->rhoCentral_baryon;
    ghl_hybrid_get_K_and_Gamma(ghl_eos,rhoC_baryon,&aK,&aGamma);
    y[0] = aK*pow(rhoC_baryon, aGamma); // Pressure, can be calcualated from central baryon density.
    y[1] = 0.0; // nu
    y[2] = 0.0; // mass
    y[3] = 0.0; // r-bar
    params->rho_baryon=rhoC_baryon;
    params->rho_energy = pow(y[0] / aK , 1.0 / aGamma) + y[0] / (aGamma - 1.0);
    params->Pinitial = y[0];
  }

  //IF PIECEWISE POLYTROPE
  else if (strcmp("Piecewise",NRml_EOS_type) == 0){
    double aK;
    double aGamma;
    double rhoC_baryon = params->rhoCentral_baryon;
    double eps;
    ghl_hybrid_get_K_and_Gamma(ghl_eos,rhoC_baryon,&aK,&aGamma);
    ghl_hybrid_compute_P_cold_and_eps_cold(ghl_eos,rhoC_baryon,&y[0],&eps);
    y[0] = aK*pow(rhoC_baryon, aGamma); //Pressure
    y[1] = 0.0; // nu                                                                                                                                               
    y[2] = 0.0; // mass                                                                                                                                                                                           
    y[3] = 0.0; // r-bar
    params->rho_baryon=rhoC_baryon;
    params->rho_energy = params->rho_baryon*(1.0+eps);
    params->Pinitial = y[0];
  }

  //IF TABULATED EOS
  else if (strcmp("Tabulated",NRml_EOS_type) == 0){
    double rhoC_baryon = params->rhoCentral_baryon;
    y[0] = ghl_tabulated_compute_P_from_rho(ghl_eos, rhoC_baryon);
    y[1] = 0.0;
    y[2] = 0.0;
    y[3] = 0.0;
    params->rho_baryon=rhoC_baryon;
    double eps = ghl_tabulated_compute_eps_from_rho(ghl_eos,rhoC_baryon);
    params->rho_energy = (params->rho_baryon)*(1.0+eps);
    params->Pinitial = y[0];
  }
  CCTK_VINFO("Got Initial Conditions!");
}

void assign_constants (double c[], struct constant_parameters *params)
{
    c[0] = params->rho_energy; // Total energy density.
    c[1] = params->rho_baryon;
    if (isnan(params->rho_energy)==1){
      c[0]=0;}
}
