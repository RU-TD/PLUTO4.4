#include "pluto.h"

extern void prizmo_rho_c(double *, double *, double *, double *, double *);
extern void prizmo_init_c();
extern void prizmo_rho_c(double *, double *, double *, double *, double *);
extern void prizmo_set_cooling_mode_c(int *);
extern void prizmo_set_heating_mode_c(int *);
extern void prizmo_set_variable_crflux_c(double *);
extern void prizmo_get_draine_flux_c(double *);
extern void prizmo_attenuate_rho_c(double *, double *, double *, double *, double *);
extern void prizmo_set_variable_nh2_incoming_c(double *);
extern void prizmo_set_variable_nco_incoming_c(double *);
extern void prizmo_set_variable_nh2_escaping_c(double *);
extern void prizmo_set_variable_nco_escaping_c(double *);
extern void prizmo_set_variable_nh2o_escaping_c(double *);
extern void prizmo_set_variable_av_c(double *);
extern void prizmo_set_variable_g0_c(double *);
extern void prizmo_x2ngas_c(double *, double *, double *);
extern void prizmo_x2n_c(double *, double *, double *);
extern void prizmo_n2x_c(double *, double *, double *);

/* ********************************************************************* */
void Chemistry(Data_Arr v, double dt, Grid *grid)
/*!
 *  Chemistry
 *********************************************************************** */
{
  int i,j,k,n;
  double x[NTRACER];
  double Av;
  double Tgas, dx;
  double dt_chem, rho_chem;
  double rho, T;
  double mu=2.;
  double pi = acos(-1e0);
  int first = 1;

  // set slab model paramters
  double Gnot = 1e1 / 1.15;
  double N2Av = 6.289e-22;

  DOM_LOOP(k,j,i){
    rho = v[RHO][k][j][i];
    Tgas = v[PRS][k][j][i]/v[RHO][k][j][i]*(KELVIN*mu);

    NTRACER_LOOP(n) x[n-TRC] = v[n][k][j][i];
    dt_chem = dt*UNIT_LENGTH/UNIT_VELOCITY;
    rho_chem = rho*UNIT_DENSITY;
    Av = irradiation.column_density[k][j][i][0]*N2Av;

    // set radiation variables
    prizmo_set_variable_g0_c(&Gnot);
    prizmo_set_variable_av_c(&Av);

    // set incoming column density to the cell
    prizmo_set_variable_nh2_incoming_c(&irradiation.column_density[k][j][i][1]);
    prizmo_set_variable_nco_incoming_c(&irradiation.column_density[k][j][i][2]);

    // set excaping column density from the cell
    prizmo_set_variable_nh2_escaping_c(&irradiation.column_density[k][j][i][1]);
    prizmo_set_variable_nco_escaping_c(&irradiation.column_density[k][j][i][2]);
    prizmo_set_variable_nh2o_escaping_c(&irradiation.column_density[k][j][i][3]);

    //print("%d %d %le\n",prank,i,irradiation.jflux[k][j][i][1]);
    // CALL PRIZMO
    prizmo_rho_c(x, &rho_chem, &Tgas, irradiation.jflux[k][j][i], &dt_chem);

    v[PRS][k][j][i] = v[RHO][k][j][i]*Tgas/(KELVIN*mu);

    NTRACER_LOOP(n) v[n][k][j][i] = x[n-TRC];
  }
}
