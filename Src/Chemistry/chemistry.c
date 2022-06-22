#include "pluto.h"

extern void prizmo_set_radial_ncol_h2_c(double *);
extern void prizmo_set_radial_ncol_co_c(double *);
extern void prizmo_set_vertical_ncol_co_c(double *);
extern void prizmo_evolve_rho_c(double *, double *, double *, double *, double *);
extern void prizmo_frac2n_c(double *, double *, double *);
extern void prizmo_n2frac_c(double *, double *, double *);

/* ********************************************************************* */
void Chemistry(Data_Arr v, double dt, Grid *grid)
/*!
 *  Chemistry
 *********************************************************************** */
{
  int i,j,k,n;
  double x[NTRACER];
  double Tgas, Tdust, dx;
  double dt_s;
  double rho, T;

  DOM_LOOP(k,j,i){
    rho = v[RHO][k][j][i]*UNIT_DENSITY;
    Tgas = v[PRS][k][j][i]/v[RHO][k][j][i]*(KELVIN*g_inputParam[MU]);

    NTRACER_LOOP(n) x[n-TRC] = v[n][k][j][i];
 
    dt_s = dt*UNIT_LENGTH/UNIT_VELOCITY;
    dx = grid->dx[IDIR][i]*UNIT_LENGTH;

    // set incoming column density to the cell
    prizmo_set_radial_ncol_h2_c(&irradiation.column_density[1][k][j][i]);
    prizmo_set_radial_ncol_co_c(&irradiation.column_density[2][k][j][i]);

    // set excaping column density from the cell
    prizmo_set_vertical_ncol_co_c(&irradiation.column_density[2][k][j][i]);

    // CALL PRIZMO
    prizmo_evolve_rho_c(x, &rho, &Tgas, irradiation.jflux[i], &dt_s);

    v[PRS][k][j][i] = v[RHO][k][j][i]*Tgas/(KELVIN*g_inputParam[MU]);

    NTRACER_LOOP(n) v[n][k][j][i] = x[n-TRC];
  }
}
