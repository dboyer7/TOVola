#include "GRHayLib.h" //Access to GRHayL library in the ET
#include "TOVola_defines.h" //Access to the ID_persist struct
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_odeiv2.h>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
    
/* Bisection index finder using binary search */
static int TOVola_bisection_idx_finder(const CCTK_REAL rr_iso, const int numpoints_arr, const CCTK_REAL *restrict r_iso_arr) {
  int x1 = 0;
  int x2 = numpoints_arr - 1;
  CCTK_REAL y1 = rr_iso - r_iso_arr[x1];
  CCTK_REAL y2 = rr_iso - r_iso_arr[x2];
  if (y1 * y2 > 0) {
    CCTK_VINFO("INTERPOLATION BRACKETING ERROR: r_iso_min = %e ?<= r_iso = %.15e ?<= %e = r_iso_max\n", r_iso_arr[0], rr_iso,
            r_iso_arr[numpoints_arr - 1]);
    CCTK_ERROR("Shutting down due to error...");
  }
  for (int i = 0; i < numpoints_arr; i++) {
    int x_midpoint = (x1 + x2) / 2;
    CCTK_REAL y_midpoint = rr_iso - r_iso_arr[x_midpoint];
    if (y_midpoint * y1 <= 0) {
      x2 = x_midpoint;
      y2 = y_midpoint;
    } else {
      x1 = x_midpoint;
      y1 = y_midpoint;
    }
    if (abs(x2 - x1) == 1) {
      // If r_iso_arr[x1] is closer to rr_iso than r_iso_arr[x2] then return x1:
      if (fabs(rr_iso - r_iso_arr[x1]) < fabs(rr_iso - r_iso_arr[x2])) {
        return x1;
      }
      // Otherwise return x2:
      return x2;
    }
  }
  CCTK_VINFO("INTERPOLATION BRACKETING ERROR: r_iso_min = %e ?<= r_iso = %.15e ?<= %e = r_iso_max\n", r_iso_arr[0], rr_iso,
            r_iso_arr[numpoints_arr - 1]);
  CCTK_ERROR("Shutting down due to error...");
}

/* Interpolation Function using Lagrange Polynomial */
static void TOVola_TOV_interpolate_1D(CCTK_REAL rr_iso,
                                      const int Interpolation_Stencil, const int Max_Interpolation_Stencil, const int numpoints_arr, const CCTK_REAL *restrict r_Schw_arr,
                                      const CCTK_REAL *restrict rho_energy_arr, const CCTK_REAL *restrict rho_baryon_arr, const CCTK_REAL *restrict P_arr,
                                      const CCTK_REAL *restrict M_arr, const CCTK_REAL *restrict expnu_arr, const CCTK_REAL *restrict exp4phi_arr,
                                      const CCTK_REAL *restrict r_iso_arr, CCTK_REAL *restrict rho_energy, CCTK_REAL *restrict rho_baryon, CCTK_REAL *restrict P,
                                      CCTK_REAL *restrict M, CCTK_REAL *restrict expnu, CCTK_REAL *restrict exp4phi) {
  
	
  const int R_idx = numpoints_arr - 1;
  const CCTK_REAL M_star = M_arr[R_idx];
  const CCTK_REAL r_iso_max_inside_star = r_iso_arr[R_idx];
  CCTK_REAL r_Schw = 0.0;
  if (rr_iso < r_iso_max_inside_star) { // If we are INSIDE the star, we need to interpollate the data to the grid.
    // For this case, we know that for all functions, f(r) = f(-r)
    if (rr_iso < 0)
      rr_iso = -rr_iso;

    // First find the central interpolation stencil index:
    int idx_mid = TOVola_bisection_idx_finder(rr_iso, numpoints_arr, r_iso_arr);

    /* Use standard library functions instead of redefining macros */
    int idxmin = MAX(0, idx_mid - Interpolation_Stencil / 2 - 1);

    // -= Do not allow the interpolation stencil to cross the star's surface =-
    // max index is when idxmin + (TOVola_Interpolation_stencil-1) = R_idx
    //  -> idxmin at most can be R_idx - TOVola_Interpolation_stencil + 1
    idxmin = MIN(idxmin, R_idx - Interpolation_Stencil + 1);

    // Ensure that Interpolation_Stencil does not exceed the maximum
    if (Interpolation_Stencil > Max_Interpolation_Stencil) {
      CCTK_ERROR("Interpolation stencil size exceeds maximum allowed.\n");
    }

    // Now perform the Lagrange polynomial interpolation:

    // First compute the interpolation coefficients:
    CCTK_REAL r_iso_sample[Max_Interpolation_Stencil];
    for (int i = idxmin; i < idxmin + Interpolation_Stencil; i++) {
      //if(i < 0 || i >= R_idx-1) { fprintf(stderr, "ERROR!\n"); exit(1); }
      r_iso_sample[i - idxmin] = r_iso_arr[i];
    }
    CCTK_REAL l_i_of_r[Max_Interpolation_Stencil];
    for (int i = 0; i < Interpolation_Stencil; i++) {
      CCTK_REAL numer = 1.0;
      CCTK_REAL denom = 1.0;
      for (int j = 0; j < Interpolation_Stencil; j++) {
        if (j != i) {
          numer *= (rr_iso - r_iso_sample[j]);
          denom *= (r_iso_sample[i] - r_iso_sample[j]);
        }
      }
      l_i_of_r[i] = numer / denom;
    }

    // Then perform the interpolation:
    *rho_energy = 0.0;
    *rho_baryon = 0.0;
    *P = 0.0;
    *M = 0.0;
    *expnu = 0.0;
    *exp4phi = 0.0;

    for (int i = idxmin; i < idxmin + Interpolation_Stencil; i++) {
      r_Schw += l_i_of_r[i - idxmin] * r_Schw_arr[i];
      *rho_energy += l_i_of_r[i - idxmin] * rho_energy_arr[i];
      *rho_baryon += l_i_of_r[i - idxmin] * rho_baryon_arr[i];
      *P += l_i_of_r[i - idxmin] * P_arr[i];
      *M += l_i_of_r[i - idxmin] * M_arr[i];
      *expnu += l_i_of_r[i - idxmin] * expnu_arr[i];
      *exp4phi += l_i_of_r[i - idxmin] * exp4phi_arr[i];
    }

  } else {
    // If we are OUTSIDE the star, the solution is just Schwarzschild.
    r_Schw = (rr_iso + M_star) + M_star * M_star / (4.0 * rr_iso); // Need to know what r_Schw is at our current grid location.
    *rho_energy = 0;
    *rho_baryon = 0;
    *P = 0;
    *M = M_star;
    *expnu = 1. - 2.0 * (M_star) / r_Schw;
    *exp4phi = (r_Schw * r_Schw) / (rr_iso * rr_iso);
  }
  //printf("%.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e hhhh\n", rr_iso, r_Schw, *rho_energy, *rho_baryon, *P, *M, *expnu, *exp4phi);
}
/*
static void TOVola_interp(TOVola_ID_persist_struct *TOVola_ID_persist, int Interpolation_Stencil, int Max_Interpolation_Stencil){

   

  //get surface values for later calculations.
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
  				TOVola_TOV_interpolate_1D(TOVola_r_iso, Interpolation_Stencil, Max_Interpolation_Stencil,
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
  }*/
  
/*
  const CCTK_REAL x = xCart[0];
  const CCTK_REAL y = xCart[1];
  const CCTK_REAL z = xCart[2];
  const CCTK_REAL r_iso = sqrt(x * x + y * y + z * z);
  // self.xxSph[1] = sp.acos(self.xx[2] / self.xxSph[0])
  const CCTK_REAL theta = acos(z / r_iso);

  // Perform pointwise interpolation to radius r using TOVola_ID_persist data
  CCTK_REAL rho_energy_val, rho_baryon_val, P_val, M_val, expnu_val, exp4phi_val;
  TOVola_TOV_interpolate_1D(r_iso, Interpolation_Stencil, Max_Interpolation_Stencil, TOVola_ID_persist->numpoints_arr,
                            TOVola_ID_persist->r_Schw_arr, TOVola_ID_persist->rho_energy_arr, TOVola_ID_persist->rho_baryon_arr, TOVola_ID_persist->P_arr, TOVola_ID_persist->M_arr,
                            TOVola_ID_persist->expnu_arr, TOVola_ID_persist->exp4phi_arr, TOVola_ID_persist->r_iso_arr, &rho_energy_val, &rho_baryon_val, &P_val, &M_val,
                            &expnu_val, &exp4phi_val);

  // Assign interpolated values to initial_data_struct
  initial_data->alpha = sqrt(expnu_val);

  // Assuming beta and B fields are zero in this context
  initial_data->betaSphorCartU0 = 0.0;
  initial_data->betaSphorCartU1 = 0.0;
  initial_data->betaSphorCartU2 = 0.0;
  initial_data->BSphorCartU0 = 0.0;
  initial_data->BSphorCartU1 = 0.0;
  initial_data->BSphorCartU2 = 0.0;

  // Metric components (assuming diagonal for simplicity)
  initial_data->gammaSphorCartDD00 = exp4phi_val;
  initial_data->gammaSphorCartDD01 = 0.0;
  initial_data->gammaSphorCartDD02 = 0.0;
  initial_data->gammaSphorCartDD11 = exp4phi_val * r_iso * r_iso;
  initial_data->gammaSphorCartDD12 = 0.0;
  initial_data->gammaSphorCartDD22 = exp4phi_val * r_iso * r_iso * sin(theta) * sin(theta);

  // Extrinsic curvature components set to zero
  initial_data->KSphorCartDD00 = 0.0;
  initial_data->KSphorCartDD01 = 0.0;
  initial_data->KSphorCartDD02 = 0.0;
  initial_data->KSphorCartDD11 = 0.0;
  initial_data->KSphorCartDD12 = 0.0;
  initial_data->KSphorCartDD22 = 0.0;

  initial_data->T4SphorCartUU00 = rho_energy_val / expnu_val;
  initial_data->T4SphorCartUU01 = 0.0;
  initial_data->T4SphorCartUU02 = 0.0;
  initial_data->T4SphorCartUU03 = 0.0;
  initial_data->T4SphorCartUU11 = P_val / exp4phi_val;
  initial_data->T4SphorCartUU12 = 0.0;
  initial_data->T4SphorCartUU13 = 0.0;
  initial_data->T4SphorCartUU22 = P_val / (exp4phi_val * r_iso * r_iso);
  initial_data->T4SphorCartUU23 = 0.0;
  initial_data->T4SphorCartUU33 = P_val / (exp4phi_val * r_iso * r_iso * sin(theta) * sin(theta));
  
}*/
