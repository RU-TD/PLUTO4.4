#include "pluto.h"
#include "particles.h"


#define init_part_thetamin 1.08 // limits the particles roughly to within the bound disc


double randn()
{
  double U1, U2, W, mult;
  static double X1, X2;
  static int call = 0;

  if (call == 1)
    {
      call = !call;
      return ((double) X2);
    }

  do
    {
      U1 = -1 + ((double) rand() / RAND_MAX) * 2;
      U2 = -1 + ((double) rand() / RAND_MAX) * 2;
      W = pow (U1, 2) + pow (U2, 2);
    }
  while (W >= 1 || W == 0);

  mult = sqrt ((-2 * log (W)) / W);
  X1 = U1 * mult;
  X2 = U2 * mult;

  call = !call;

  return ((double) X1);
}

/* ************************************************************** */
void ParticlesInit(Particle * pl, Grid *grid, const Data *d, int j)
/*
 *
 * PURPOSE
 *
 *  Initialize the particles
 *
 * ARGUMENTS
 *
 *  PART_HEAD  *list : pointer to the list of particles. Needed to allocate
 *                     each new particle
 *  Grid       *grid : The pluto grid, found in other functions in the code.
 *  const Data *d    : The full dataset containing the physical variables.
 *                     Needed for interpolation of pyhisical variables in the
 *                     particle position.
 *  int         j    : integer related to the number of the particle.
 *                     This function is called from inside a loop over all
 *                     particles, where j is the running index.
 *
 * USAGE
 *
 **************************************************************** */
{
  int i;
  double R;
  double part_domain[DIMENSIONS][2] = {
    {grid->xbeg_glob[0], grid->xend_glob[0]},
    #if DIMENSIONS > 1
      {init_part_thetamin, grid->xend_glob[1]},
    #endif
    #if DIMENSIONS > 2
      {grid->xbeg_glob[2], grid->xend_glob[2]}
    #endif
  };
  static double micron = 1e-4 / UNIT_LENGTH;

  static int first_call = 1;

  if (first_call == 1){
    srand(0);
    first_call = 0;
  }
  pl->identity=j;

  pl->rank = 0;

  pl->flag = 'I';

  // Intialize particles positions randomly within the particle domain
  for (i=0;i<DIMENSIONS;++i){
    pl->coor[i] = part_domain[i][0] + (part_domain[i][1] - part_domain[i][0]) * (rand() / (1.+RAND_MAX));
    pl->mom[i]=0.0; // Intialize particles momentum to zero
  }

  // Define particles radius in cgs units
  if (j%6 == 0) pl->radius = 1e-4*micron;
  else if (j%6 == 1) pl->radius = 1e-2*micron;
  else if (j%6 == 2) pl->radius = micron;
  else if (j%6 == 3) pl->radius = 1e2*micron;
  else if (j%6 == 4) pl->radius = 1e3*micron;
  else pl->radius = 1e4*micron;

#if GEOMETRY == SPHERICAL
  R = pl->coor[0]*sin(pl->coor[1]);
#elif GEOMETRY == POLAR
  R = pl->coor[0];
#elif GEOMETRY == CARTESIAN
  R = sqrt(POW2(pl->coor[0])+POW2(pl->coor[1]));
#else
  printLog("GEOMETRY NOT SUPPORTED: CHOOSE BETWEEN POLAR AND SPHERICAL");
  QUIT_PLUTO(1);
#endif

  // Set Keplerian speed for the Lagrangian particles
  double OmegaK = 2.0 * CONST_PI / (R * sqrt(R));
#if GEOMETRY == POLAR
  //pl->coor[2] = 0.0;
  pl->mom[1] += R*R*OmegaK;
#elif GEOMETRY == SPHERICAL
  //pl->coor[1] = CONST_PI/2.0;
  pl->mom[2] += R*R*OmegaK;
#else
  printLog("GEOMETRY NOT SUPPORTED: CHOOSE BETWEEN POLAR AND SPHERICAL\n");
  QUIT_PLUTO(1);
#endif

}

/* ******************************************************************************* */
void ParticlesInitVelocity(Particle *pl, double ***uu[], Grid *grid)
/*
 *
 * PURPOSE
 *
 *  Initialize the particles
 *
 * ARGUMENTS
 *
 *  PART_HEAD  *list : pointer to the list of particles. Needed to allocate
 *                     each new particle
 *  Grid       *grid : The pluto grid, found in other functions in the code.
 *  const Data *d    : The full dataset containing the physical variables.
 *                     Needed for interpolation of pyhisical variables in the
 *                     particle position.
 *  int         j    : integer related to the number of the particle.
 *                     This function is called from inside a loop over all
 *                     particles, where j is the running index.
 *
 * USAGE
 *
 **************************************************************** */

{
  double cs;
  calculateGravitation(uu, grid, pl);
  pl->tstop = calculateDrag(uu, grid, pl);

#if PARTICLES_TURBULENCE == YES
  if (pl->flag == 'R') exit;
/*
  double R = pl->coor[0]*sin(pl->coor[1]);
  double OmegaK = sqrt(g_mu/R)/R;

  double Hg = 0.05*R;
  double tau = pl->tstop*OmegaK;

 #if EOS == IDEAL
  cs = OmegaK*Hg;
 #elif EOS == ISOTHERMAL
  cs = g_isoSoundSpeed/sqrt(R);
 #endif

  //Diffusion parameters for dust and gas
  double Sc = POW2(1.0 + POW2(tau))/(1.0 + 4.0*tau*tau);
  double Dg = g_inputParam[ALPHA]*cs*Hg;
  double Dd = Dg/Sc;

  double Hp = Hg/sqrt(1.0 + tau/g_inputParam[ALPHA]*(1.0 + 2.0*tau)/(1.0 + tau));

  //Setting initial vertical position and speed of dust particles
  double z = Hp*randn();
#if GEOMETRY == SPHERICAL
  pl->coor[1] = acos(z/pl->coor[0]);
#elif GEOMETRY == POLAR
  pl->coor[2] = z;
#endif
*/
  pl->flag = 'R';
#endif
}

void ParticlesAnalysis(Particle * particles, int nop, int nstep, double glob_time)
/*
 *
 * PURPOSE
 *
 *   Analyse the particles (and write on disk?)
 *
 *
 *
 * NOTE
 *
 *   This is actually just a lopp over the particles, with the option
 *   to destroy them if needed (whenever check==0 the particle is deleted)
 *
 *
 *
 **************************************************************** */
{
	int check=1;
        int i = 0;
	Particle *pl;

	if(particles != NULL ){
	  for(i = 0; i < nop; i++) {

            pl=&(particles[i]);

			/* BEGINNING OF USERDEFINED SECTION */



			/* END OF USERDEFINED SECTION       */

			if(check==0){
		    /* if a particle is destroied, update nop, and destroy particle*/

				--nop;
			}

        } /* endwhile*/
    } /* endif*/
}
