#include <stdio.h>
#include <math.h>

// interfaces
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

// PLUTO variables
#define NTRACER 31
#define NPHOTO 1000
#define NFLX 0
#define NIONS 0

#define IDX_CHEM_H (NFLX + NIONS + 0)
#define IDX_CHEM_CH (NFLX + NIONS + 1)
#define IDX_CHEM_C (NFLX + NIONS + 2)
#define IDX_CHEM_H2 (NFLX + NIONS + 3)
#define IDX_CHEM_CH3 (NFLX + NIONS + 4)
#define IDX_CHEM_CH2 (NFLX + NIONS + 5)
#define IDX_CHEM_CH4 (NFLX + NIONS + 6)
#define IDX_CHEM_OH (NFLX + NIONS + 7)
#define IDX_CHEM_O (NFLX + NIONS + 8)
#define IDX_CHEM_H2O (NFLX + NIONS + 9)
#define IDX_CHEM_CO (NFLX + NIONS + 10)
#define IDX_CHEM_O2 (NFLX + NIONS + 11)
#define IDX_CHEM_CH2j (NFLX + NIONS + 12)
#define IDX_CHEM_CHj (NFLX + NIONS + 13)
#define IDX_CHEM_CH3j (NFLX + NIONS + 14)
#define IDX_CHEM_Hej (NFLX + NIONS + 15)
#define IDX_CHEM_He (NFLX + NIONS + 16)
#define IDX_CHEM_Hj (NFLX + NIONS + 17)
#define IDX_CHEM_Cj (NFLX + NIONS + 18)
#define IDX_CHEM_Oj (NFLX + NIONS + 19)
#define IDX_CHEM_H2j (NFLX + NIONS + 20)
#define IDX_CHEM_COj (NFLX + NIONS + 21)
#define IDX_CHEM_E (NFLX + NIONS + 22)
#define IDX_CHEM_H3j (NFLX + NIONS + 23)
#define IDX_CHEM_CH4j (NFLX + NIONS + 24)
#define IDX_CHEM_OHj (NFLX + NIONS + 25)
#define IDX_CHEM_CH5j (NFLX + NIONS + 26)
#define IDX_CHEM_H2Oj (NFLX + NIONS + 27)
#define IDX_CHEM_H3Oj (NFLX + NIONS + 28)
#define IDX_CHEM_HCOj (NFLX + NIONS + 29)
#define IDX_CHEM_O2j (NFLX + NIONS + 30)

#define IMAX 200

// ****************************
int main (void){
  double xall[NTRACER][IMAX] = {0e0};
  double x[NTRACER] = {0e0};
  double n[NTRACER] = {0e0};
  double xold[IMAX] = {0e0};
  double jflux[NPHOTO] = {0e0};
  double tgas[IMAX] = {0e0};
  double z[IMAX] = {0e0};
  double Av[IMAX] = {0e0};
  double pi = acos(-1e0);
  double rho, dt, crflux, spy, t, Gnot, ngas, N2Av, zmin_log, zmax_log;
  double zold, NH2, NCO, NH2O, Ntot, dz, tend;
  FILE *fout;

  spy = 365. * 24. * 3600.;  // seconds per year
  dt = spy; // max integration time, s
  tend = spy * 1e7; // total integration time
  crflux = 5e-17;  // cosmic rays ionization rate, 1/s

  // initial chemical conditions (cm-3)
  ngas = 1e3;
  n[IDX_CHEM_H2] = 0.5 * ngas;
  n[IDX_CHEM_O] = 3e-4 * ngas;
  n[IDX_CHEM_C] = 1e-4 * ngas;
  n[IDX_CHEM_Hej] = 1e-1 * ngas;
  n[IDX_CHEM_E] = n[IDX_CHEM_Hej];

  // convert to mass fraction (and also get rho the total mass density)
  prizmo_n2x_c(n, x, &rho);

  // loop on grid points to copy initial conditions
  for(int i=0;i<IMAX;i++){
    // loop on species
    for(int j=0;j<NTRACER;j++){
      xall[j][i] = x[j];
    }
    // set temperature to grid points
    tgas[i] = 5e1;
  }

  // set slab model paramters
  Gnot = 1e1 / 1.15;
  N2Av = 6.289e-22;

  // loop to compute grid position and corresponding Av
  for(int i=0;i<IMAX;i++){
    zmin_log = log10(1e-6 / N2Av / ngas);
    zmax_log = log10(3e1 / N2Av / ngas);
    z[i] = pow(1e1, i * (zmax_log - zmin_log) / (IMAX - 1e0) + zmin_log);
    Av[i] = z[i] * ngas * N2Av;
  }

  // init prizmo
  prizmo_init_c();

  // set cosmic rays ionization
  prizmo_set_variable_crflux_c(&crflux);

  // open file to write
  fout = fopen("slab_final.out", "w");

  t = 0e0;  // initial time, s
  // loop on time-steps to evolve chemistry in time
  for (;;) {

    // print current and final time
    printf("%e %e\n", t / spy, tend / spy);

    // set Draine multifrequency intensity
    prizmo_get_draine_flux_c(jflux);
    // loop to multiply 2*pi*G0
    for(int j = 0; j < NPHOTO; j++){
      jflux[j] *= Gnot * 2e0 * pi;
    }

    // values at the first time of the grid
    zold = 0e0; // old grid position, cm
    NH2 = 0e0; // H2 column density, cm-2
    NCO = 0e0; // CO column density, cm-2
    NH2O = 0e0; // H2O column density, cm-2
    Ntot = 0e0; // total column density, cm-2

    // loop on grid points
    for(int i=0;i<IMAX;i++){
      dz = z[i] - zold; // compute cell size, cm

      // compute number densities from mass density
      prizmo_x2n_c(xold, n, &rho);

      // compute column densities
      NH2 += n[IDX_CHEM_H2] * dz;
      NCO += n[IDX_CHEM_CO] * dz;
      NH2O += n[IDX_CHEM_H2O] * dz;
      Ntot += (2e0 * n[IDX_CHEM_H2] + n[IDX_CHEM_H]) * dz;

      // loop to stride 2D array
      for(int j = 0; j < NTRACER; j++){
        x[j] = xall[j][i];
      }

      // attenuate radiation in the grid cell
      prizmo_attenuate_rho_c(jflux, x, &tgas[i], &dz, &rho);

      // set radiation variables
      prizmo_set_variable_g0_c(&Gnot);
      prizmo_set_variable_av_c(&Av[i]);

      // set incoming column density to the cell
      prizmo_set_variable_nh2_incoming_c(&NH2);
      prizmo_set_variable_nco_incoming_c(&NCO);

      // set excaping column density from the cell
      prizmo_set_variable_nh2_escaping_c(&NH2);
      prizmo_set_variable_nco_escaping_c(&NCO);
      prizmo_set_variable_nh2o_escaping_c(&NH2O);

      // CALL PRIZMO HERE!
      prizmo_rho_c(x, &rho, &tgas[i], jflux, &dt);

      // save to file
      fprintf(fout, "%e %e %e %e %e %e %e %e %e %e %e %e ",
        t / spy, z[i], tgas[i], 0e0, Ntot, NH2, NCO, Av[i], crflux, 0e0, 0e0, 0e0);
      for(int j = 0; j < NTRACER; j++){
        fprintf(fout, "%e ", xall[j][i]);
      }
      fprintf(fout, "\n");

      // strided variable to 2D array
      // and store current variables to compute difference with the next cell
      for(int j = 0; j < NTRACER; j++){
        xall[j][i] = x[j];
        xold[j] = x[j];
      }
      zold = z[i];

    }
    fprintf(fout, "\n");

    // increment time-step
    dt *= 1.1e0;
    // update total time
    t += dt;
    // exit when max time is reached
    if(t > tend){
      break;
    }
  }
  fclose(fout);

  printf("done!");

  return 0;
}
