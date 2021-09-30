/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Photoevaporating PPD as described
  in section 2 of Owen et al., MNRAS (2010) 401, 1415
  and in Picogna et al., MNRAS (2019), 487, 691

  This version has been modified to include an embedded planet. It 
  relies heavily on Giovanni Picogna's original implementation 
  from 2016.

  \author G. Picogna (picogna@usm.lmu.de)
  \date   Dec 17, 2020
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "init.h"

#if GEOMETRY != SPHERICAL
print("! Only spherical geometry is implemented for this test case.\n");
QUIT_PLUTO(1);
// All functions defined in this file only give correct results with a
// spherical geometry. If you ever get rid of this condition, please check
// everything carfully.
#endif

#if EOS != IDEAL
print("! Only ideal gas implemented for this test case.\n");
QUIT_PLUTO(1);
#endif


IrradiationData irradiation;
BoundaryData innerBoundaryData;
InitialData initialData;
PlanetData planet;


/* ********************************************************************* */
void Init (double *us, double x1, double x2, double x3)
/*
 *
 *
 *
 *********************************************************************** */
{
  double R;
  double gm = UNIT_G * g_inputParam[MS];
  static int first_call = 1;

  if (first_call) {
    if (abs(sqrt(gm) - 2*CONST_PI) > 1e-6)
    {
      print("! sqrt(G*M_star) is not equal to 2pi (%.6f != %.6f). Please check your units!\n", sqrt(gm), 2*CONST_PI);
      QUIT_PLUTO(1);
    }

    #if ROTATING_FRAME == YES
      g_OmegaZ = OmegaP;
    #endif

    g_gamma = g_inputParam[GAMMA];
    g_smallDensity = RHO_MIN;
    g_smallPressure = g_smallDensity * T_MIN / (KELVIN * g_inputParam[MU]);

    first_call = 0;
  }

  #ifdef FARGO
    R = x1 * sin(x2);
    us[FARGO_W] = 2.*CONST_PI / sqrt(R) - R * g_OmegaZ;
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
  int i, j, k, l, m, id, size[3];
  long int offset;
  double R, OmegaK;
  double *x1 = grid->x[IDIR];
  double *x2 = grid->x[JDIR];
  double *x3 = grid->x[KDIR];
  char grid_file[32];

  #if PLANET_HEATING == SML16
    double x, y, z;
    double dpx, dpy, dp;
    double Tp;
    static double a_sml16 = -0.616428949637868; // fit to linear part of Szulágyi et al. 2016
    static double b_sml16 = 231.7971169551819;
  #endif


  if (grid->xl_glob[JDIR][0] < 0) {
    print("! Ghost cells extend into negative x2-direction. Please choose a higher inner boundary.\n");
    QUIT_PLUTO(1);
  }
  
  CDInit(grid);  

  // initialize planet properties
  planet.M = g_inputParam[MP] * CONST_Mearth / UNIT_MASS;
  planet.rhill = pow(g_inputParam[MP] / (g_inputParam[MP] + g_inputParam[MS] * CONST_Msun / CONST_Mearth) / 3.0, 
    1.0 / 3.0);
  planet.x = 1.0;
  planet.y = 0.0;
  planet.M_accreted = 0; // The cumulative mass that has been removed in the vicinity of the planet in the
  // local domain. Only used if REMOVE_ACCRETING_MASS == YES. Note that the mass is not added to the planet

  sprintf (grid_file, INPUT_GRID_FILE);

  // Interpolate density (1st record in the file)
  id = InputDataOpen(INPUT_DATA_FILE,grid_file," ", 0, CENTER);

  InputDataGridSize (id, size); /* Get grid data size */
  offset = (long)size[0]*(long)size[1]*(long)size[2]; // Offset of a single data block

  TOT_LOOP(k,j,i) {
    d->Vc[RHO][k][j][i] = InputDataInterpolate (id,x1[i],x2[j],x3[k]);
    d->Vc[RHO][k][j][i] = MAX(d->Vc[RHO][k][j][i], g_smallDensity);
  }
  InputDataClose(id);

  // Interpolate r-velocity (2nd record in the file)
  id = InputDataOpen(INPUT_DATA_FILE,grid_file," ",1*offset, CENTER);
  TOT_LOOP(k,j,i) {
    d->Vc[VX1][k][j][i] = InputDataInterpolate (id,x1[i],x2[j],x3[k]);
  }
  InputDataClose(id);

  // Interpolate theta-velocity (3rd record in the file)
  id = InputDataOpen(INPUT_DATA_FILE,grid_file," ",2*offset, CENTER);
  TOT_LOOP(k,j,i) {
    d->Vc[VX2][k][j][i] = InputDataInterpolate (id,x1[i],x2[j],x3[k]);
  }
  InputDataClose(id);

  /* Interpolate phi-velocity (4th record in the file) and subtract rotational velocity. 
     Make sure the input data provides the absolute velocity in a static frame */
  id = InputDataOpen(INPUT_DATA_FILE,grid_file," ",3*offset, CENTER);
  TOT_LOOP(k,j,i) {
    R = x1[i] * sin(x2[j]);
    d->Vc[VX3][k][j][i] = InputDataInterpolate (id,x1[i],x2[j],x3[k]);
    #ifdef FARGO
      OmegaK = 2.*CONST_PI / (R*sqrt(R));
      d->Vc[VX3][k][j][i] -= R*(OmegaK - g_OmegaZ);
    #else
      d->Vc[VX3][k][j][i] -= R * g_OmegaZ;
    #endif  
  }
  InputDataClose(id);

  // Interpolate Tdust (5th record in the file)
  id = InputDataOpen(INPUT_DATA_FILE,grid_file," ",4*offset, CENTER);
  TOT_LOOP(k,j,i) {
    irradiation.Tdust[k][j][i] = InputDataInterpolate (id,x1[i],x2[j],x3[k]);
    irradiation.Tdust[k][j][i] = MAX(irradiation.Tdust[k][j][i], T_MIN);

    #if PLANET_HEATING == SML16
      // calculate the temperature based on a fit to figure 3 in Szulágyi et al. 2016.
      // the dust temperature is only used to limit the minimum gas temperature, so we can safely increase it here.
      // note that it is only valid for a fixed planet 

      x = x1[i] * sin(x2[j]) * cos(x3[k]);
      y = x1[i] * sin(x2[j]) * sin(x3[k]);
      z = x1[i] * cos(x2[j]);

      dpx = x - planet.x;
      dpy = y - planet.y;
      dp = (sqrt(dpx * dpx + dpy * dpy + z * z) + 1. / UNIT_LENGTH) / planet.rhill; // add 1cm to prevent (unlikely) division by 0

      Tp = b_sml16 * pow(dp, a_sml16) * ( 1 - 1 / (exp(-(dp-3)/0.5) + 1) );
      if (Tp > irradiation.Tdust[k][j][i]) {
        irradiation.Tdust[k][j][i] = Tp;
      }
    #endif

    irradiation.Tdust[k][j][i] = MIN(irradiation.Tdust[k][j][i], T_MAX);
  }
  InputDataClose(id);

  // Interpolate Tgas (6th record in the file)
  id = InputDataOpen(INPUT_DATA_FILE,grid_file," ",5*offset, CENTER);
  TOT_LOOP(k,j,i) {
    irradiation.Tgas[k][j][i] = InputDataInterpolate (id,x1[i],x2[j],x3[k]);
    irradiation.Tgas[k][j][i] = MAX(irradiation.Tgas[k][j][i], irradiation.Tdust[k][j][i]);
    irradiation.Tgas[k][j][i] = MIN(irradiation.Tgas[k][j][i], T_MAX);
  }
  InputDataClose(id);

  // Interpolate the column density at the inner boundary (7th record in the file)
  if (innerBoundaryData.cd == NULL) {
    // First allocate the memory
    innerBoundaryData.cd = (double **)Array2D(NX3_TOT, NX2_TOT, sizeof(double));
  }
  id = InputDataOpen(INPUT_DATA_FILE,grid_file," ",6*offset, CENTER);
  KTOT_LOOP(k) { 
    JTOT_LOOP(j) {
      innerBoundaryData.cd[k][j] = InputDataInterpolate (id,grid->x_glob[IDIR][grid->gbeg[IDIR]],x2[j],x3[k]);
    }
  }
  InputDataClose(id);

  // Initialize the initial conditions
  if (initialData.rho == NULL) {
    // First allocate the memory
    initialData.rho = (double ***)Array3D(NX3_TOT, NX2_TOT, NX1_TOT, sizeof(double));
    initialData.vx1 = (double ***)Array3D(NX3_TOT, NX2_TOT, NX1_TOT, sizeof(double));
    initialData.vx2 = (double ***)Array3D(NX3_TOT, NX2_TOT, NX1_TOT, sizeof(double));
  }
  TOT_LOOP(k,j,i) {
    initialData.rho[k][j][i] = d->Vc[RHO][k][j][i];
    initialData.vx1[k][j][i] = d->Vc[VX1][k][j][i];
    initialData.vx2[k][j][i] = d->Vc[VX2][k][j][i];
  }

  // CalculateS(grid, d); // calculate the column density
  // UpdateGasTemp(grid, d); // update the gas temperature based on the column density

  // Set the pressure
  TOT_LOOP(k,j,i) {
    d->Vc[PRS][k][j][i] = irradiation.Tgas[k][j][i] * d->Vc[RHO][k][j][i] / (KELVIN * g_inputParam[MU]);
  }
  
  // now add passive scalar tracers
#if NTRACER == 2
  // we will add two separate tracers inside and outside the gap
  double tr_Rin = 4.4*CONST_au / UNIT_LENGTH;
  double tr_Rout = 6.0*CONST_au / UNIT_LENGTH;
  DOM_LOOP(k,j,i) {
    R = x1[i] * sin(x2[j]);
    if (R < tr_Rin) {
      d->Vc[TRC][k][j][i] = 1.0;
    } else if (R > tr_Rout) {
      d->Vc[TRC + 1][k][j][i] = 1.0;
    }
  }

#elif NTRACER == 39
  // we will add three vertical layers of tracers, limited by the following meridional angles
  double tr_thetabounds[] = {1.5, 1.571, 1.34, 1.41, 1.19, 1.25};
  // in each layer we will add 13 sections distributed along the outer gap edge
  double tr_rmin = 6*CONST_au/UNIT_LENGTH;
  double tr_rmax = 6.5*CONST_au/UNIT_LENGTH;
  double tr_dphi = 6./360. * 2*CONST_PI; // 6° azimuthal extent of a tracer in radians
  double tr_phiint; // azimuthal spacing between tracers in radians
  double tr_phibounds[26];
  double tr_phitmp = 0;
  for(m = 0; m<13; m++) {
    tr_phiint = 39./360. * 2*CONST_PI;
    if (tr_phitmp + tr_dphi <= CONST_PI/4. || tr_phitmp + tr_dphi >= 7*CONST_PI/4.) {
      // we want more tracers close to the planet location, so we use smaller spacing here
      tr_phiint = 8./360. * 2*CONST_PI;
    }
    tr_phibounds[2*m] = tr_phitmp;
    tr_phibounds[2*m+1] = tr_phitmp + tr_dphi;
    tr_phitmp += tr_dphi + tr_phiint;
  }
  
  int tr_id = 0;
  DOM_LOOP(k,j,i) {
    if (x1[i] >= tr_rmin && x1[i] <= tr_rmax) {
      for (l=0; l<3; l++) {
        if (x2[j] >= tr_thetabounds[2*l] && x2[j] <= tr_thetabounds[2*l+1]) {
          for (m=0; m<13; m++) {
            if (x3[k] >= tr_phibounds[2*m] && x3[k] <= tr_phibounds[2*m+1]) {
              tr_id = 13*l + m;
              if (tr_id < NTRACER) {
                d->Vc[TRC + tr_id][k][j][i] = 1.0;
              }
            }
          }
        }
      }
    }
  }
#endif
}


/* ********************************************************************* */
void UserDefBoundary(const Data *d, RBox *box, int side, Grid *grid)
/*
 * 
 *
 *********************************************************************** */
{
  int i, j, k, nv;
  double r, th, R;
  static int previous_stepNum;

  if (side == 0) { /* -- check solution inside domain -- */

    double rin, rout, f_damp, tau_damp;

    CalculateS(grid, d); // calculate the column density

    UpdateGasTemp(grid, d); // update the gas temperature based on the column density

    DOM_LOOP(k, j, i)
    {
      r = grid->x[IDIR][i];
      th = grid->x[JDIR][j];
      R = r*sin(th);

      if (r < damping_rin || r > damping_rout) {

        rin = grid->x_glob[IDIR][grid->gbeg[IDIR]];
        rout = grid->x_glob[IDIR][grid->gend[IDIR]];

        if (r < damping_rin) {
          f_damp = (damping_rin - r) / (damping_rin - rin);
        } else {
          f_damp = (r - damping_rout) / (rout - damping_rout);
        }
        f_damp = f_damp*f_damp;
        if (f_damp > 1.) { 
          f_damp = 1.;
        }

        tau_damp = f_damp * g_dt / (g_inputParam[TAU_DAMP] * R * sqrt(R)); // R * sqrt(R) = period in c.u.

        d->Vc[RHO][k][j][i] -= (d->Vc[RHO][k][j][i] - initialData.rho[k][j][i]) * tau_damp;
        d->Vc[VX1][k][j][i] -= (d->Vc[VX1][k][j][i] - initialData.vx1[k][j][i]) * tau_damp;
        d->Vc[VX2][k][j][i] -= (d->Vc[VX2][k][j][i] - initialData.vx2[k][j][i]) * tau_damp;
      }

#if REMOVE_ACCRETING_MASS == YES
      double phi, x, y, z;
      double dp2;

      phi = grid->x[KDIR][k];
      x = r * sin(th) * cos(phi);
      y = r * sin(th) * sin(phi);
      z = r * cos(th);

      // remove a fraction of the mass within 1/2 R_Hill similar to D'Angelo & Lubow (2008)
      dp2 = (x - planet.x) * (x - planet.x) + (y - planet.y) * (y - planet.y) + z * z;
      if (dp2 < 0.25 * planet.rhill * planet.rhill)
      {
        d->Vc[RHO][k][j][i] *= (1. - g_dt * OmegaP / (2*CONST_PI * g_inputParam[TAU_ACC]));
        if (g_stepNumber > previous_stepNum)
        {
          // this condition ensures that we only add the removed mass during the prediction
          // step, and not the other integration steps. this is not very accurate,
          // but it is the simplest way to remove mass in a way that is  independent of the
          // chosen time-stepping and our treatment of accretion is very unrealistic anyway
          // and it only serves to avoid overdensities in the vicinity of the planet
          planet.M_accreted += d->Vc[RHO][k][j][i] * grid->dV[k][j][i] * g_dt * OmegaP / (2*CONST_PI * g_inputParam[TAU_ACC]);
        }
      }
#endif

      if (d->Vc[RHO][k][j][i] < g_smallDensity) {
        d->Vc[RHO][k][j][i] = g_smallDensity;
      }

      // update pressure
      d->Vc[PRS][k][j][i] = (d->Vc[RHO][k][j][i] * irradiation.Tgas[k][j][i]) / (KELVIN * g_inputParam[MU]);
    }

    RBox dom_box;
    RBoxDefine (IBEG, IEND, JBEG, JEND, KBEG, KEND, CENTER, &dom_box);
    PrimToCons3D (d->Vc, d->Uc, &dom_box); 

    previous_stepNum = g_stepNumber;
  }

  if (side == X1_BEG){
    X1_BEG_LOOP(k,j,i){
        // NVAR_LOOP(nv) d->Vc[nv][k][j][i] = d->Vc[nv][k][j][IBEG];
        NVAR_LOOP(nv) d->Vc[nv][k][j][i] = d->Vc[nv][k][j][2*IBEG - i - 1];
        d->Vc[VX1][k][j][i] *= -1.0;
        // d->Vc[VX1][k][j][i] = 0.0;
        // d->Vc[VX2][k][j][i] = 0.0;
        #ifdef FARGO
          d->Vc[VX3][k][j][i] = 0.0;
        #else
          r = grid->x[IDIR][i];
          th = grid->x[JDIR][j]; 
          R = r * sin(th);
          d->Vc[VX3][k][j][i] = 2*CONST_PI/sqrt(R) - R * g_OmegaZ;
        #endif
    }
  }

  if (side == X1_END){
    X1_END_LOOP(k,j,i){
        // NVAR_LOOP(nv)  d->Vc[nv][k][j][i] = d->Vc[nv][k][j][IEND];
        NVAR_LOOP(nv)  d->Vc[nv][k][j][i] = d->Vc[nv][k][j][2*IEND - i + 1];
        d->Vc[VX1][k][j][i] *= -1.0;
        // d->Vc[VX1][k][j][i] = 0.0;
        // d->Vc[VX2][k][j][i] = 0.0;
        #ifdef FARGO
          d->Vc[VX3][k][j][i] = 0.0;
        #else
          r = grid->x[IDIR][i];
          th = grid->x[JDIR][j]; 
          R = r * sin(th);
          d->Vc[VX3][k][j][i] = 2*CONST_PI/sqrt(R) - R * g_OmegaZ;
        #endif
    }
  }
}


#if (BODY_FORCE & POTENTIAL)
/* ************************************************************************ */
double BodyForcePotential(double x1, double x2, double x3)
/*
 *
 * 
 *
 *************************************************************************** */
{
  double x, y, z;
  double d, rsm, phiPlanet, muPlanet_eff;
  static double gm = 4. * CONST_PI * CONST_PI; // G*M_star

  x = x1 * sin(x2) * cos(x3);
  y = x1 * sin(x2) * sin(x3);
  z = x1 * cos(x2);

#if ROTATING_FRAME == NO
  planet.x = cos(OmegaP * g_time);
  planet.y = sin(OmegaP * g_time);
#endif

  if (g_time < g_inputParam[PLANETRAMPUPTIME])
  {
    // slowly ramp up the planet mass in order to avoid strong shocks
    muPlanet_eff = UNIT_G * (PlanetM0 + 
      (planet.M - PlanetM0) * sin(CONST_PI / 2. / g_inputParam[PLANETRAMPUPTIME] * g_time));
  }
  else
  {
    muPlanet_eff = UNIT_G * planet.M;
  }

  d = sqrt((x - planet.x) * (x - planet.x) + (y - planet.y) * (y - planet.y) + z * z);
  rsm = g_inputParam[RSM] * planet.rhill;
  phiPlanet = -muPlanet_eff / d;
  if (d <= rsm)
  {
    phiPlanet *= (pow(d / rsm, 4) - 2 * pow(d / rsm, 3) + 2 * d / rsm);
  }

  return - gm / x1 + phiPlanet;
}
#endif


/* ********************************************************************* */
void Analysis(const Data *d, Grid *grid)
/*
 *
 *
 *********************************************************************** */
{
  double sendbuf[2], recvbuf[2];

  sendbuf[0] = planet.M_accreted;
  sendbuf[1] = calculateTorque(d, grid, NULL);

#ifdef PARALLEL
  MPI_Reduce(sendbuf, recvbuf, 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
#else
  recvbuf[0] = planet.M_accreted;
  recvbuf[1] = sendbuf[1];
#endif

  if (prank == 0)
  {
    char fname[512];
    static double tpos = -1.0;
    FILE *fp;

    sprintf(fname, "%s/planet.dat", RuntimeGet()->output_dir);
    if (g_stepNumber == 0)
    {
      fp = fopen(fname, "w");
      fprintf(fp, "# %7s  %12s  %14s\n", "t", "M_acc", "torque_z");
    }
    else
    {
      if (tpos < 0.)
      {
        char sline[512];
        fp = fopen(fname, "r");
        while (fgets(sline, 512, fp))
        {
        }
        sscanf(sline, "%lf\n", &tpos);
        fclose(fp);
      }
      fp = fopen(fname, "a");
    }

    if (g_time > tpos)
    {
      fprintf(fp, "%12.6e  %12.6e  %12.6e \n", g_time, recvbuf[0], recvbuf[1]);
    }
    fclose(fp);
  }
}


/* ********************************************************************* */
double calculateTorque(const Data *d, Grid *grid, double ***elemental_torques)
{
  /*
 * Calculates the z-component of the torque acting on the planet. 
 * note that this is includes only one hemisphere and the total torque 
 * would be twice this value.
 * if elemental_torques is given, the torques excerted by each zone are
 * calculated and stored in it.
 * 
 *********************************************************************** */
  int i, j, k;
  double r, th, phi, x, y, z, V;
  double rp, dp, dpx, dpy, fb;
  double apx, apy;
  double sum_apx = 0, sum_apy = 0;

  DOM_LOOP(k, j, i)
  {
    r = grid->x[IDIR][i];
    th = grid->x[JDIR][j];
    phi = grid->x[KDIR][k];
    x = r * sin(th) * cos(phi);
    y = r * sin(th) * sin(phi);
    z = r * cos(th);
    V = grid->dV[k][j][i];

    // gravitational acceleration acting on the planet due to the mass in this zone
    dpx = x - planet.x;
    dpy = y - planet.y;
    dp = sqrt(dpx * dpx + dpy * dpy + z * z) + 1. / UNIT_LENGTH; // add 1cm to prevent (unlikely) division by 0
    fb = 1. / (exp(-(dp / planet.rhill - 0.8) / 0.08) + 1);      // tapering function of Kley et al. (2009)
    apx = (fb * UNIT_G * V * d->Vc[RHO][k][j][i] / (dp * dp * dp) * dpx);
    apy = (fb * UNIT_G * V * d->Vc[RHO][k][j][i] / (dp * dp * dp) * dpy);
    sum_apx += apx;
    sum_apy += apy;
    if (elemental_torques != NULL)
    {
      elemental_torques[k][j][i] = planet.M * (planet.x * apy - planet.y * apx);
    }
  }

  return planet.M * (planet.x * sum_apy - planet.y * sum_apx);
}


void UpdateGasTemp(Grid *grid, const Data *d) {
  int i, j, k;
  double r;

  DOM_LOOP(k, j, i) {
    r = grid->x[IDIR][i];

    if (irradiation.S[k][j][i] < 2.5e22) {
      if (g_inputParam[LX] > 0) {
        irradiation.Tgas[k][j][i] = CalculateT(irradiation.S[k][j][i], d->Vc[RHO][k][j][i], irradiation.Tdust[k][j][i], r);
      } else {
        irradiation.Tgas[k][j][i] = d->Vc[PRS][k][j][i] / d->Vc[RHO][k][j][i] * KELVIN * g_inputParam[MU];
      }

      if (irradiation.Tgas[k][j][i] > T_MAX) {
        irradiation.Tgas[k][j][i] = T_MAX;
      } else if (irradiation.Tgas[k][j][i] < irradiation.Tdust[k][j][i]) {
        irradiation.Tgas[k][j][i] = irradiation.Tdust[k][j][i];
      }
    } else {
      irradiation.Tgas[k][j][i] = irradiation.Tdust[k][j][i];
    }
  }
}

double CalculateT(double column, double rho, double Tdust, double dist)
{
  double Temp;
  double en, xi, zeta2;
  double b, c, d, m;

  double mpart = g_inputParam[MU] * CONST_amu;
  double Lx = g_inputParam[LX];

  en = rho * UNIT_DENSITY / mpart;
  xi = Lx / (en * dist * dist * UNIT_LENGTH * UNIT_LENGTH);

  zeta2 = log10(xi);

  if (zeta2 < -8.0)
  {
    Temp = Tdust;
  }
  else if (zeta2 > -2.0)
  {
    Temp = T_MAX;
  }
  else
  {
    Temp = read_table_temperature(column, zeta2);
  }
  return Temp;
}

double read_table_temperature(double column, double xi)
{
  int jlo, jhi, jmid, klo, khi, kmid, i;
  double lgT, column_mid;
  double ylo, yhi, t;
  static double *b_fit, *c_fit, *d_fit, *m_fit;
  static int n_fits = 10;
  static double fit_interval = 2.5e21;

  /* -------------------------------------------
   *    Read tabulated cooling function
   * ------------------------------------------- */

  if (b_fit == NULL)
  {
    print(" > Reading table from disk...\n");
    FILE *ftemp = fopen("fitParameter.dat", "r");
    if (ftemp == NULL)
    {
      print("! fitParameter.dat could not be found.\n");
      QUIT_PLUTO(1);
    }
    b_fit = ARRAY_1D(n_fits, double);
    c_fit = ARRAY_1D(n_fits, double);
    d_fit = ARRAY_1D(n_fits, double);
    m_fit = ARRAY_1D(n_fits, double);

    i = 0;
    while (fscanf(ftemp, "%lf %lf %lf %lf \n", &b_fit[i], &c_fit[i], &d_fit[i], &m_fit[i]) != EOF)
    {
      i++;
    }

    print(" > Finished reading table from disk...\n");
  }

  /* ----------------------------------------------
   *    Table lookup
   * ----------------------------------------------*/
  if (column <= fit_interval)
  {
    jlo = 0;
    lgT = d_fit[jlo] + (1.5 - d_fit[jlo]) / pow((1.0 + pow((xi / c_fit[jlo]), b_fit[jlo])), m_fit[jlo]);
  }
  else if (column >= 2.5e22)
  {
    jhi = n_fits - 1;
    lgT = d_fit[jhi] + (1.5 - d_fit[jhi]) / pow((1.0 + pow((xi / c_fit[jhi]), b_fit[jhi])), m_fit[jhi]);
  }
  else
  {
    jlo = column / fit_interval;
    if (jlo >= n_fits - 1)
    {
      jlo = n_fits - 2;
    }
    jhi = jlo + 1;

    ylo = d_fit[jlo] + (1.5 - d_fit[jlo]) / pow((1.0 + pow((xi / c_fit[jlo]), b_fit[jlo])), m_fit[jlo]);
    yhi = d_fit[jhi] + (1.5 - d_fit[jhi]) / pow((1.0 + pow((xi / c_fit[jhi]), b_fit[jhi])), m_fit[jhi]);

    //Linear interpolation between the fits
    t = (column - fit_interval * jhi) / fit_interval;
    lgT = ylo + (yhi - ylo) * t;
  }

  return pow(10., lgT);
}


void CalculateS(Grid *grid, const Data *d)
{
  int k, j, i, n;
  double column_density = 0.;
  double density = 0.0, dr = 0.0;
  double mpart = g_inputParam[MU] * CONST_amu;
  MPI_Status status;

  if (irradiation.neighbour.receive_rank != irradiation.neighbour.send_rank)
  { //check if more than one process is in the communicator (nproc > 1)
    if (irradiation.neighbour.receive_rank == -1 && irradiation.neighbour.send_rank != -1)
    { //no neighbor in between the star and the current domain
      n = 0;
      for (k = KBEG; k <= KEND; k++)
      {
        for (j = JBEG; j <= JEND; j++)
        {
          column_density = innerBoundaryData.cd[k][j];
          irradiation.S[k][j][IBEG] = column_density;
          for (i = IBEG; i <= IEND; i++)
          {
            density = d->Vc[RHO][k][j][i] * UNIT_DENSITY;
            dr = grid->dx[IDIR][i] * UNIT_LENGTH;
            column_density += density / mpart * dr;
            irradiation.S[k][j][i + 1] = column_density;
          }
          irradiation.data_buffer[n] = irradiation.S[k][j][IEND + 1];
          irradiation.column_density_offset[n] = 0.;
          n++;
        }
      }
      MPI_Send(irradiation.data_buffer, sizeof(double) * (NX2 * NX3), MPI_CHAR, irradiation.neighbour.send_rank, 1001, MPI_COMM_WORLD);
    }
    else
    {
      MPI_Recv(irradiation.column_density_offset, sizeof(double) * (NX2 * NX3), MPI_CHAR, irradiation.neighbour.receive_rank, 1001, MPI_COMM_WORLD, &status);

      n = 0;
      for (k = KBEG; k <= KEND; k++)
      {
        for (j = JBEG; j <= JEND; j++)
        {
          column_density = irradiation.column_density_offset[n];
          irradiation.S[k][j][IBEG] = column_density;
          for (i = IBEG; i <= IEND; i++)
          {
            density = d->Vc[RHO][k][j][i] * UNIT_DENSITY;
            dr = grid->dx[IDIR][i] * UNIT_LENGTH;
            column_density += density / mpart * dr;
            irradiation.S[k][j][i + 1] = column_density;
          }
          irradiation.data_buffer[n] = irradiation.S[k][j][IEND + 1];
          n++;
        }
      }

      if (irradiation.neighbour.send_rank != -1)
      {
        MPI_Send(irradiation.data_buffer, sizeof(double) * (NX2 * NX3), MPI_CHAR, irradiation.neighbour.send_rank, 1001, MPI_COMM_WORLD);
      }
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
}


void CDInit(Grid *grid)
{
  int n = 0,k,j,i;
  irradiation.neighbour.receive_rank = -1;
  irradiation.neighbour.send_rank = -1;

  irradiation.S = (double***)Array3D(NX3_TOT, NX2_TOT, NX1_TOT, sizeof(double));
  irradiation.Tdust = (double***)Array3D(NX3_TOT, NX2_TOT, NX1_TOT, sizeof(double));
  irradiation.Tgas = (double***)Array3D(NX3_TOT, NX2_TOT, NX1_TOT, sizeof(double));
  irradiation.data_buffer = (double *)Array1D(NX2 * NX3, sizeof(double));
  irradiation.column_density_offset = (double *)Array1D(NX2 * NX3, sizeof(double));

  for(n = 0; n < NX2 * NX3; ++n) {
    irradiation.column_density_offset[n] = 0.;
  }

  FindCommunicationNeighbours(grid);
}

void FindCommunicationNeighbour(int current_rank, LocalDomainInfo *domain_info_array, int nproc, CommunicationNeighbour *cn)
{
  int n = 0;
  LocalDomainInfo *cdi = &domain_info_array[current_rank];

  cn->receive_rank = -1;
  cn->send_rank = -1;

  for (n = 0; n < nproc; ++n)
  {
    if (n != current_rank)
    {
      if (cdi->x1_begin - 1 == domain_info_array[n].x1_end && cdi->x2_begin == domain_info_array[n].x2_begin && cdi->x2_end == domain_info_array[n].x2_end && cdi->x3_begin == domain_info_array[n].x3_begin && cdi->x3_end == domain_info_array[n].x3_end)
      {
        cn->receive_rank = n;
      }
      else if (cdi->x1_end + 1 == domain_info_array[n].x1_begin && cdi->x2_begin == domain_info_array[n].x2_begin && cdi->x2_end == domain_info_array[n].x2_end && cdi->x3_begin == domain_info_array[n].x3_begin && cdi->x3_end == domain_info_array[n].x3_end)
      {
        cn->send_rank = n;
      }
    }
  }
}

void FindCommunicationNeighbours(Grid *grid)
{
  int nproc = 0, n = 0;
  LocalDomainInfo ldi;
  LocalDomainInfo *domain_info_array = NULL;
  LocalDomainInfo *item = NULL;
  CommunicationNeighbour *neighbour_info_array = NULL;
  MPI_Status status;

  MPI_Comm_size(MPI_COMM_WORLD, &nproc);

  ldi.x1_begin = grid->beg[IDIR];
  ldi.x1_end = grid->end[IDIR];

  ldi.x2_begin = grid->beg[JDIR];
  ldi.x2_end = grid->end[JDIR];

  ldi.x3_begin = grid->beg[KDIR];
  ldi.x3_end = grid->end[KDIR];

  if (prank == 0)
  {
    domain_info_array = malloc(sizeof(LocalDomainInfo) * nproc);
    neighbour_info_array = malloc(sizeof(CommunicationNeighbour) * nproc);
  }

  MPI_Gather(&ldi, sizeof(LocalDomainInfo), MPI_CHAR, domain_info_array, sizeof(LocalDomainInfo), MPI_CHAR, 0, MPI_COMM_WORLD);

  if (prank == 0)
  {
    for (n = 0; n < nproc; ++n)
    {
      FindCommunicationNeighbour(n, domain_info_array, nproc, &neighbour_info_array[n]);
    }

    for (n = 0; n < nproc; ++n)
    {
      if (n != 0)
      {
        MPI_Send(&neighbour_info_array[n], sizeof(CommunicationNeighbour), MPI_CHAR, n, 0, MPI_COMM_WORLD);
      }
    }

    irradiation.neighbour.receive_rank = neighbour_info_array[0].receive_rank;
    irradiation.neighbour.send_rank = neighbour_info_array[0].send_rank;

    free(domain_info_array);
    free(neighbour_info_array);
  }
  else
  {
    MPI_Recv(&irradiation.neighbour, sizeof(CommunicationNeighbour), MPI_CHAR, 0, 0, MPI_COMM_WORLD, &status);
  }
}
