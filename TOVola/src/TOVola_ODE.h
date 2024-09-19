#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <gsl/gsl_odeiv2.h>
#include "GRHayLib.h"

//ETK interface. 
#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

// This is TOVola's ODE functions for setting up the ODE system in GSL.
// There are also additional function definitions for an exception handler and termination condition
// When going through the integration.

//struct to hold information about the EOS
struct TOVola_parameters { 
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

//Prototypes

//An exception handler, in case numerical error causes pressure to go negative.
void TOVola_exception_handler(double x, double y[]); 

//Termination condition to end the integration
int TOVola_do_we_terminate(double x, double y[], struct TOVola_parameters *params); 

//Function to evaluate rho_baryon and rho_energy, as well as epsilon. Values required for the TOV.
void TOVola_evaluate_rho_and_eps(double x, const double y[], struct TOVola_parameters *params);

//The function the actually holds the TOV equations for GSL
int TOVola_ODE(double x, const double y[], double dydx[], void *params);

//Empty function. GSL requires and jacobian to setup, but it is not neccessary for the TOV solution. Can just leave empty.
int TOVola_jacobian_placeholder(double t, const double y[], double *dfdy, double dfdt[], void *params);

//Calculate initial pressure to get the solution rolling
void TOVola_get_initial_condition(double y[], struct TOVola_parameters*params);

//Save rho_baryon and rho_energy to memory for later use.
void TOVola_assign_constants(double c[], struct TOVola_parameters *params);



void TOVola_exception_handler(double x, double y[])
{
//Just a check to make sure pressure doesn't accidently goes negative do to numerical roundoff error.
//Neccessary the closer you get to the surface and pressure goes to zero.
    if (y[0] < 0) {
        y[0] = 0;
    } 
}

int TOVola_do_we_terminate(double x, double y[], struct TOVola_parameters *params)
{
  DECLARE_CCTK_PARAMETERS
  if (CCTK_EQUALS("Tabulated",TOVola_EOS_type)){
  	const CCTK_REAL PMin = exp(ghl_eos->lp_of_lr[0]); //PMin is not zero on the table, so we don't want to exceed table limits
  	if (y[0] <= PMin){
  		return 1;
  	}
  }
  else if (y[0] <= 0.0){ //PMin is not an issue when not using a table, though.
  	return 1;
  }
  
  return 0;
    // return 1; for termination.
}

void TOVola_evaluate_rho_and_eps(double x, const double y[], struct TOVola_parameters *params)
{

  DECLARE_CCTK_PARAMETERS
  //IF SIMPLE POLYTROPE
  if (CCTK_EQUALS("Simple",TOVola_EOS_type)){
    //Declare EOS info
    double aK;
    double aGamma;
    double aRho_baryon = params->rho_baryon;
    double eps;
    double aPress;
    
    //Use GRHayL function to calculate our current rho_baryon and rho_energy
    ghl_hybrid_get_K_and_Gamma(ghl_eos,aRho_baryon,&aK,&aGamma);
    params->rho_baryon = pow(y[0]/aK, 1.0 / aGamma);
    aRho_baryon = params->rho_baryon;
    ghl_hybrid_compute_P_cold_and_eps_cold(ghl_eos,aRho_baryon,&aPress,&eps);
    params->rho_energy = params->rho_baryon*(1.0+eps);}

  //IF PIECEWISE POLYTROPE
  else if (CCTK_EQUALS("Piecewise",TOVola_EOS_type)){
    //Basically identical to Simple Polytrope, you just have more regions.
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
  else if (CCTK_EQUALS("Tabulated",TOVola_EOS_type)){
    const CCTK_REAL PMin = exp(ghl_eos->lp_of_lr[0]);
    if(y[0] > PMin){ //Assure you are not exceeding table bounds
      //Use GRHayL function to find our current rho_baryon and rho_energy on the table.
      params->rho_baryon = ghl_tabulated_compute_rho_from_P(ghl_eos, y[0]);
      double eps = ghl_tabulated_compute_eps_from_rho(ghl_eos,params->rho_baryon);
      params->rho_energy = (params->rho_baryon)*(1+eps);
      }
    else{
      //Outside the star, densities are zero.
      params->rho_baryon = 0;
      params->rho_energy = 0;}
  }
    
}

int TOVola_ODE(double x, const double y[], double dydx[], void *params)
{
    //Need to have the correct value for rho_energy in the TOVs.
    //Only assigned later when the step has gone all the way through.
    TOVola_evaluate_rho_and_eps(x,y,params);

    // Dereference the struct to use rho_energy
    double rho_energy = (*(struct TOVola_parameters*)params).rho_energy;

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

int TOVola_jacobian_placeholder(double x, const double y[], double *dfdy, double dydx[], void *params)
{

  //Jacobian isn't necessary for TOVola to run properly, but GSL wants SOME function
  //Leave blank, doesn't affect final results
  
  return 0;
}

void TOVola_get_initial_condition (double y[], struct TOVola_parameters *params)
{
  DECLARE_CCTK_PARAMETERS
  //IF SIMPLE POLYTROPE
  if (CCTK_EQUALS("Simple",TOVola_EOS_type)){
    //Declare EOS info
    double aK;
    double aGamma;
    double rhoC_baryon = params->rhoCentral_baryon;
    
    //Use GRHayL to find K and Gamma and calculate initial conditions
    ghl_hybrid_get_K_and_Gamma(ghl_eos,rhoC_baryon,&aK,&aGamma);
    y[0] = aK*pow(rhoC_baryon, aGamma); // Pressure, can be calcualated from central baryon density.
    y[1] = 0.0; // nu
    y[2] = 0.0; // mass
    y[3] = 0.0; // r-bar
    
    //Assign the initial conditions
    params->rho_baryon=rhoC_baryon;
    params->rho_energy = pow(y[0] / aK , 1.0 / aGamma) + y[0] / (aGamma - 1.0);
    params->Pinitial = y[0];
  }

  //IF PIECEWISE POLYTROPE
  else if (CCTK_EQUALS("Piecewise",TOVola_EOS_type)){
    //Declare EOS info
    double aK;
    double aGamma;
    double rhoC_baryon = params->rhoCentral_baryon;
    double eps;
    
    //Use GRHayL to find K and Gamma and calculate initial conditions
    ghl_hybrid_get_K_and_Gamma(ghl_eos,rhoC_baryon,&aK,&aGamma);
    ghl_hybrid_compute_P_cold_and_eps_cold(ghl_eos,rhoC_baryon,&y[0],&eps);
    y[0] = aK*pow(rhoC_baryon, aGamma); //Pressure
    y[1] = 0.0; // nu                                                                                                                                               
    y[2] = 0.0; // mass                                                                                                                                                                                           
    y[3] = 0.0; // r-bar
    
    //Assign the initial conditions
    params->rho_baryon=rhoC_baryon;
    params->rho_energy = params->rho_baryon*(1.0+eps);
    params->Pinitial = y[0];
  }

  //IF TABULATED EOS
  else if (CCTK_EQUALS("Tabulated",TOVola_EOS_type)){
    //Use GRHayL to find initial pressure on the table
    double rhoC_baryon = params->rhoCentral_baryon;
    y[0] = ghl_tabulated_compute_P_from_rho(ghl_eos, rhoC_baryon);
    y[1] = 0.0;
    y[2] = 0.0;
    y[3] = 0.0;
    
    //Assign the initial conditions
    params->rho_baryon=rhoC_baryon;
    double eps = ghl_tabulated_compute_eps_from_rho(ghl_eos,rhoC_baryon);
    params->rho_energy = (params->rho_baryon)*(1.0+eps);
    params->Pinitial = y[0];
  }
  CCTK_VINFO("Got Initial Conditions!");
}

void TOVola_assign_constants (double c[], struct TOVola_parameters *params)
{
    //This is where we actually assign the densities
    c[0] = params->rho_energy; // Total energy density.
    c[1] = params->rho_baryon;
    if (isnan(params->rho_energy)==1){ //In case outside the surface.
      c[0]=0;}
}

//TOVola: by David Boyer
//auxillary functions for solving the ODE.
