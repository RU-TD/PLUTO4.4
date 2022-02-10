#include "pluto.h"

void lkrm_rho_c(double *, double *, double *, double *);

/* ********************************************************************* */
void Chemistry(Data_Arr v, double dt)
/*!
 *  Chemistry
 *********************************************************************** */
{
  int i,j,k,n;
//  double rate;
  double x[NTRACER];
//  double Tgas = 100.;
  double rho, T;
  double tstep;
  double variable_g0, variable_crflux, variable_av;
  //tstep = dt*3.15576e7;
  tstep = 100.*3.15576e7;
  variable_g0 = 1.;
  variable_crflux = 5.e-17;

  DOM_LOOP(k,j,i){ 
    variable_av = irradiation.S[k][j][i]/1.8e21;
    rho = v[RHO][k][j][i]*5.9738398e-7;
    T   = v[PRS][k][j][i]/v[RHO][k][j][i]*KELVIN*g_inputParam[MU];
    set_variable_g0_c(&variable_g0);
    set_variable_crflux_c(&variable_crflux);
    set_variable_av_c(&variable_av);
    NTRACER_LOOP(n) x[n-TRC] = v[n][k][j][i];
    lkrm_rho_c(x, &rho, &T, &tstep);
    NTRACER_LOOP(n) v[n][k][j][i] = x[n-TRC];
  }

}
