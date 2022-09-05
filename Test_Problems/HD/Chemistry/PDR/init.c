/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief PRD Region.

  \author G. Picogna T. Grassi (picogna@usm.lmu.de)
  \date   Dec 15, 2020

  \b References
     - Grassi, T. et al. (2020). 
       "Modelling thermochemical processes in protoplanetary discs I: numerical methods".       
       MNRAS, Volume 494, Issue 3, pp.4471-4491

*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

extern void prizmo_get_rho_c(double *, double *);
extern void prizmo_n2frac_c(double *, double *);
extern void prizmo_set_d2g_c(double);

/* ********************************************************************* */
void Init (double *us, double x1, double x2, double x3)
/*
 *
 *********************************************************************** */
{
  int nv;
  double ngas=g_inputParam[N_GAS];
  double rhod;
  double x[NTRACER], n[NTRACER];
  #if EOS == IDEAL
   g_gamma = g_inputParam[GAMMA_EOS];
  #endif

  us[VX1] = 0.0;

  prizmo_set_d2g_c(1.e-2);

  NTRACER_LOOP(nv) {
    n[nv-TRC] = 0.;
    x[nv-TRC] = 0.;
  }

  n[IDX_CHEM_H2-TRC] = ngas/2.;
  n[IDX_CHEM_C-TRC] = ngas*1.e-4;
  n[IDX_CHEM_O-TRC] = ngas*3.e-4;
  n[IDX_CHEM_He-TRC] = ngas*1.e-1;

  // convert to mass fraction (and also get rho the total mass density)
  prizmo_get_rho_c(n, &rhod);
  prizmo_n2frac_c(n, x);

  NTRACER_LOOP(nv) {
    us[nv] = x[nv-TRC];
  }

  us[RHO] = rhod/(UNIT_DENSITY);
  #if EOS == IDEAL
   us[PRS] = us[RHO]*g_inputParam[TEMP]/(g_inputParam[MU]*KELVIN);
  #endif
}

/* ********************************************************************* */
void InitDomain (Data *d, Grid *grid)
/*! 
 * Assign initial condition by looping over the computational domain.
 * Called after the usual Init() function to assign initial conditions
 * on primitive variables.
 * Value assigned here will overwrite those prescribed during Init().
 *
 *
 *********************************************************************** */
{
  read_jflux();
}

/* ********************************************************************* */
void Analysis (const Data *d, Grid *grid)
/* 
 *
 *
 *********************************************************************** */
{
}
/* ********************************************************************* */
void UserDefBoundary (const Data *d, RBox *box, int side, Grid *grid) 
/*
 *
 *********************************************************************** */
{
  int k,j,i;
  int l, n;
  double density, Tgas, dr;
  double x[NTRACER];
  double jflux[NPHOTO];

  if (side == 0) 
  {    /* -- check solution inside domain -- */
    calculate_ColumnDensity(d->Vc, grid); // calculating the column density at each cell
    calculate_Attenuation(d->Vc, grid);
  }
}

