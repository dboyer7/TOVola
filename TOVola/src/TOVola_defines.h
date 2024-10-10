#pragma once

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

/* Structure to hold TOV data that will become the official ID after normalization */
typedef struct {
        CCTK_REAL *restrict r_Schw_arr;
        CCTK_REAL *restrict rho_energy_arr;
        CCTK_REAL *restrict rho_baryon_arr;
        CCTK_REAL *restrict P_arr;
        CCTK_REAL *restrict M_arr;
        CCTK_REAL *restrict expnu_arr;
        CCTK_REAL *restrict r_iso_arr;
        CCTK_REAL *restrict exp4phi_arr;
        int numpoints_arr;
} TOVola_ID_persist_struct;
