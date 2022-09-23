/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief XEUV Irradiated Disk.

  \author G. Picogna T. Grassi (picogna@usm.lmu.de)
  \date   Sep 06, 2022

  \b References
     - Grassi, T. et al. (2020). 
       "Modelling thermochemical processes in protoplanetary discs I: numerical methods".       
       MNRAS, Volume 494, Issue 3, pp.4471-4491
     - Ercolano, B. et al. (2009).
       "X-Ray Irradiated Protoplanetary Disk Atmospheres. II. Predictions from Models in Hydrostatic Equilibrium"
       MNRAS, Volume 699, Issue 2, pp.1639-1649

*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

extern void prizmo_get_rho_c(double *, double *);
extern void prizmo_n2frac_c(double *, double *);
extern void prizmo_set_d2g_c(double);
void read_jflux();

/* ********************************************************************* */
void Init (double *us, double x1, double x2, double x3)
/*
 *
 *********************************************************************** */
{
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
  int i,j,k,m,nv,id,size[3];
  long int offset;
  double *x1 = grid->x[IDIR];
  double *x2 = grid->x[JDIR];
  double *x3 = grid->x[KDIR];
  char grid_file[32];
  double x[NTRACER], n[NTRACER];
  double rhod, ms1, R, r, Omega, OmegaK, ngas;
  double unit_mass = UNIT_DENSITY*UNIT_LENGTH*UNIT_LENGTH*UNIT_LENGTH;
  double mpart = g_inputParam[GAMMA_EOS]*CONST_amu;

  #if EOS == IDEAL
   g_gamma = g_inputParam[GAMMA_EOS];
  #else
   print ("! Only ideal gas implemented for this test case.\n");
   QUIT_PLUTO(1);
  #endif

  read_jflux();

  prizmo_set_d2g_c(1.e-2);

  ms1       = g_inputParam[MSTAR]*CONST_Msun/unit_mass;
  g_mu      = CONST_G/POW2(UNIT_VELOCITY)/UNIT_LENGTH*unit_mass*ms1;

  sprintf(grid_file,"grid0.out");

  id = InputDataOpen("Temp0.dbl", grid_file, " ", 0, CENTER);
 
  InputDataGridSize(id, size); // Get grid data size 
  offset = (long)size[0]*(long)size[1]*(long)size[2]; // Offset of a single data block 

  // Interpolate density
  TOT_LOOP(k,j,i) d->Vc[RHO][k][j][i] = InputDataInterpolate(id,x1[i],x2[j],x3[k]);
  InputDataClose(id);

  // Interpolate dust temperature
  id = InputDataOpen("Temp0.dbl", grid_file, " ", offset, CENTER);
  TOT_LOOP(k,j,i) d->Vc[PRS][k][j][i] = InputDataInterpolate(id,x1[i],x2[j],x3[k]);
  InputDataClose(id);

  TOT_LOOP(k,j,i) {

    // init chemical condition
    NTRACER_LOOP(nv) {
      n[nv-TRC] = 0.;
      x[nv-TRC] = 0.;
    }

    ngas = d->Vc[RHO][k][j][i]*UNIT_DENSITY/mpart;

    n[IDX_CHEM_H2-TRC] = ngas*5.e-1;
    n[IDX_CHEM_C-TRC] = ngas*1.e-4;
    n[IDX_CHEM_O-TRC] = ngas*3.e-4;
    n[IDX_CHEM_He-TRC] = ngas*1.e-1;

    // convert to mass fraction (and also get rho the total mass density)
    prizmo_get_rho_c(n, &rhod);
    prizmo_n2frac_c(n, x);

    NTRACER_LOOP(nv) {
      d->Vc[nv][k][j][i] = x[nv-TRC];
    }

    d->Vc[RHO][k][j][i] = rhod/(UNIT_DENSITY);

    if (d->Vc[RHO][k][j][i] < g_smallDensity) {
      d->Vc[RHO][k][j][i] = g_smallDensity;
    }

    //d->Vc[PRS][k][j][i] *= d->Vc[RHO][k][j][i]/(KELVIN*g_inputParam[MU]);

   #if GEOMETRY == SPHERICAL
    r  = x1[i];
    R  = r*sin(x2[j]);
   #elif GEOMETRY == CYLINDRICAL
    R = x1[i];
   #endif

   #if ROTATING_FRAME == YES
    g_OmegaZ  = sqrt(g_mu);
   #endif
    d->Vc[VX1][k][j][i] = 0.0;
    d->Vc[VX2][k][j][i] = 0.0;
    OmegaK   = sqrt(g_mu/R)/R;
    Omega    = OmegaK * sqrt((2.2+0.5)*POW2(0.05) + 1.0 + 0.5 - 0.5*R/r);
    d->Vc[VX3][k][j][i] = R*(Omega);
  }

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
  int   i, j, k, nv;
  double *x1,*x2,*x3;
  double R, z, r;
  double minP;
  double damping, OmegaK;
  double damping_outer = 38.0;

  x1 = grid->x[IDIR];
  x2 = grid->x[JDIR];
  x3 = grid->x[KDIR];

  if (side == 0) {    /* -- check solution inside domain -- */

    calculate_ColumnDensity(d->Vc, grid); // calculating the column density at each cell
    calculate_Attenuation(d->Vc, grid);

    TOT_LOOP(k,j,i){

      #if GEOMETRY == SPHERICAL
        r  = x1[i];
        R = r * sin(x2[j]);
      #else
        QUIT_PLUTO(1);
      #endif

      if (r > damping_outer) {
        d->flag[k][j][i]   |= FLAG_INTERNAL_BOUNDARY;
      }

      if (d->Vc[RHO][k][j][i]  < g_smallDensity) {
        d->Vc[RHO][k][j][i] = g_smallDensity;
      }

    }
  }
}

/*********************************************************************** */
void read_jflux()
/*
 * Import radiation flux at 1 au from input file
 *
 *********************************************************************** */
{
    //TODO: verify that the file exist and print en error statement otherwise
    FILE *fout;
    double Jscale = 1.e1;
    int n;

    // load radiation at 1 AU
    fout = fopen("runtime_data/radiation_field.dat", "r");
    NPHOTO_LOOP(n) {
        fscanf(fout, "%le", &irradiation.jflux0[n]);
        irradiation.jflux0[n] *= 4.*CONST_PI*Jscale;
    }
    fclose(fout);
}

#if (BODY_FORCE & VECTOR)
/* ************************************************************************ */
void BodyForceVector(double *v, double *g, double x1, double x2, double x3)
/*
 *
 *
 *
 *************************************************************************** */
{
}
#endif

#if (BODY_FORCE & POTENTIAL)
/* ************************************************************************ */
double BodyForcePotential(double x1, double x2, double x3)
/*
 *
 *
 *
 *************************************************************************** */
{
  double phi  = - g_mu/x1;
  return phi;
}
#endif

