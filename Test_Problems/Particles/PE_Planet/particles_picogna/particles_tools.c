#include "pluto.h"
#include "particles.h"


/* ******************************************************************************* */
void ParticlesInterpol(double *a_interpol, Particle *pl,
		   double ***uu[], Grid *grid, int VAR)
/*
 *
 * Perform interpolation of the variable VAR on the coordinates of the particle pl
 * The result is saved in the a_interpol variable.
 *
 * The interpolation is performed on 4 grid cells (2D) or 8 grid cells (3D)

 ********************************************************************************* */
{
  /* variables for interpolation*/
  int i, il, ir, jr, jl, kr, kl, ind0, ind1, ind2;
  double ix, iy, iz;

  DIM_EXPAND(ind0 = pl->cell[0];, ind1 = pl->cell[1];, ind2 = pl->cell[2];);

  /*-------------------------------------------------------*/
  /*                                                       */
  /*   This Perform Lagrangian interpolation of velocity   */
  /*            base of polynoms: < 1, x, y, xy >          */
  /*                                                       */
  /*-------------------------------------------------------*/

  if (grid->x[0][ind0] > pl->coor[0]) {
    il = ind0 - 1;
    ir = ind0;
  } else {
    il = ind0;
    ir = ind0 + 1;
  }
  ix = (pl->coor[0] - grid->x[0][il]) / (grid->x[0][ir] - grid->x[0][il]);

#if ( DIMENSIONS == 2 ||DIMENSIONS ==3 )

  if (grid->x[1][ind1] > pl->coor[1]) {
    jl = ind1 - 1;
    jr = ind1;
  } else {
    jl = ind1;
    jr = ind1 + 1;
  }
  iy = (pl->coor[1] - grid->x[1][jl]) / (grid->x[1][jr] - grid->x[1][jl]);

#endif

#if DIMENSIONS ==3
  if (grid->x[2][ind2] > pl->coor[2]) {
    kl = ind2 - 1;
    kr = ind2;
  } else {
    kl = ind2;
    kr = ind2 + 1;
  }
  iz = (pl->coor[2] - grid->x[2][kl]) / (grid->x[2][kr] - grid->x[2][kl]);
#endif
  //print("ix %e iy %e iz %e number %d\n",ix,iy,iz,pl->identity);


#if PHYSICS == MHD || PHYSICS == RMHD
  if (VAR == VX1 || VAR == BX1) {
#else
  if (VAR == VX1) {
#endif


#if DIMENSIONS == 1
    a_interpol[0] = uu[var][0][0][il]
      + ix * (uu[VAR][0][0][ir] - uu[VAR][0][0][il]);
#endif

#if DIMENSIONS == 2

    for (i = 0; i < DIMENSIONS; ++i) {
      a_interpol[i] =
        (1. - ix)*(uu[VAR + i][0][jl][il] + iy * (uu[VAR + i][0][jr][il] - uu[VAR + i][0][jl][il]))
        + ix * (uu[VAR + i][0][jl][ir] + iy * (uu[VAR + i][0][jr][ir] - uu[VAR + i][0][jl][ir]));
    }

#endif

#if DIMENSIONS == 3

    /* We use a formulation which asks for 18 calls to variables and 8 products*/
    for (i = 0; i < DIMENSIONS; ++i) {
      a_interpol[i] =
        (1. - ix)*(uu[VAR + i][kl][jl][il]
        + iz * (uu[VAR + i][kr][jl][il] - uu[VAR + i][kl][jl][il])
        + iy * ((uu[VAR + i][kl][jr][il] - uu[VAR + i][kl][jl][il])
        + iz * (uu[VAR + i][kr][jr][il] - uu[VAR + i][kl][jr][il]
        - uu[VAR + i][kr][jl][il] + uu[VAR + i][kl][jl][il])
        )
        )
        + ix * (uu[VAR + i][kl][jl][ir]
        + iz * (uu[VAR + i][kr][jl][ir] - uu[VAR + i][kl][jl][ir])
        + iy * ((uu[VAR + i][kl][jr][ir] - uu[VAR + i][kl][jl][ir])
        + iz * (uu[VAR + i][kr][jr][ir] - uu[VAR + i][kl][jr][ir]
        - uu[VAR + i][kr][jl][ir] + uu[VAR + i][kl][jl][ir])
        )
        );

    }
#endif
  } else {


#if DIMENSIONS == 1
    *a_interpol = uu[VAR][0][0][il]
      + ix * (uu[VAR][0][0][ir] - uu[VAR][0][0][il]);
#endif

#if DIMENSIONS == 2
    *a_interpol =
      (1. - ix)*(uu[VAR][0][jl][il] + iy * (uu[VAR][0][jr][il] - uu[VAR][0][jl][il]))
      + ix * (uu[VAR][0][jl][ir] + iy * (uu[VAR][0][jr][ir] - uu[VAR][0][jl][ir]));


#endif

#if DIMENSIONS == 3

    /* We use a formulation which asks for 18 calls to variables and 8 products*/

    *a_interpol =
      (1. - ix)*(uu[VAR][kl][jl][il]
      + iz * (uu[VAR][kr][jl][il] - uu[VAR][kl][jl][il])
      + iy * ((uu[VAR][kl][jr][il] - uu[VAR][kl][jl][il])
      + iz * (uu[VAR][kr][jr][il] - uu[VAR][kl][jr][il]
      - uu[VAR][kr][jl][il] + uu[VAR][kl][jl][il])
      )
      )
      + ix * (uu[VAR][kl][jl][ir]
      + iz * (uu[VAR][kr][jl][ir] - uu[VAR][kl][jl][ir])
      + iy * ((uu[VAR][kl][jr][ir] - uu[VAR][kl][jl][ir])
      + iz * (uu[VAR][kr][jr][ir] - uu[VAR][kl][jr][ir]
      - uu[VAR][kr][jl][ir] + uu[VAR][kl][jl][ir])
      )
      );
#endif

  }

}


/* ********************************************************** */
/*void ParticlesLocate(Particle *pl, Grid *grid) */
 void ParticlesLocate(double *CoordArr, int *CellArr, Grid *grid)
/*
 * determine the cell of the particle, even in stretched grids
 * works by dividing the interval in half, and checking in which
 * half the particle is, then repeating untill the interval is
 * basically only 1 cell
 * It converges fast, although apparently there could be better methods
 *
 * *************************************************************** */
{
    int l_ind, r_ind;
    int m_ind, i, j;
    int true = 1;
#ifdef PARALLEL
    extern int prank;
#endif
    for (i = 0; i < DIMENSIONS; ++i) {

#ifdef PARALLEL
        l_ind = grid->nghost[i];
        r_ind = grid->np_int[i] + grid->nghost[i] - 1;
#else
        l_ind = grid->nghost[i];
        r_ind = grid->np_int_glob[i] + grid->nghost[i] - 1;
#endif
        while (true) {

            m_ind = l_ind + (r_ind - l_ind) / 2;

            if (CoordArr[i] <= grid->xr[i][m_ind]) {
                r_ind = m_ind;
            } else {
                l_ind = m_ind + 1;
            }

            if (l_ind == r_ind) break;

        }
        CellArr[i] = r_ind;
    }
}
void transformSphericalCartesian(double * cartesian,double * spherical,double r, double theta, double phi)
{
  cartesian[0]  = spherical[0]*sin(theta)*cos(phi);
  cartesian[0] +=-spherical[2]*sin(phi);
  cartesian[0] += spherical[1]*cos(theta)*cos(phi);
  cartesian[1]  = spherical[0]*sin(theta)*sin(phi);
  cartesian[1] += spherical[2]*cos(phi);
  cartesian[1] += spherical[1]*cos(theta)*sin(phi);
  cartesian[2]  = spherical[0]*cos(theta);
  cartesian[2] += 0.0;
  cartesian[2] +=-spherical[1]*sin(theta);
}
void transformCartesianSpherical(double * cartesian,double * spherical,double r, double theta, double phi)
{
  spherical[0]  = cartesian[0]*sin(theta)*cos(phi);
  spherical[0] += cartesian[1]*sin(theta)*sin(phi);
  spherical[0] += cartesian[2]*cos(theta);
  spherical[1]  = cartesian[0]*cos(theta)*cos(phi);
  spherical[1] += cartesian[1]*cos(theta)*sin(phi);
  spherical[1] +=-cartesian[2]*sin(theta);
  spherical[2]  =-cartesian[0]*sin(phi);
  spherical[2] += cartesian[1]*cos(phi);
  spherical[2] += 0.0;
}

void transformParticleSpherical(Particle * pl)
{
    /* Only for fully implicit solver, translates velocity and acceleration
     *
     */
    int i;
    double temp[DIMENSIONS];
    transformCartesianSpherical(pl->mom,temp,pl->coor[0],pl->coor[1],pl->coor[2]);
    for(i=0;i<DIMENSIONS;i++)
        pl->mom[i] = temp[i];
#if 0
    transformCartesianSpherical(pl->acceleration,temp,pl->coor[0],pl->coor[1],pl->coor[2]);
    for(i=0;i<DIMENSIONS;i++)
        pl->acceleration[i] = temp[i];
#endif
}
void transformParticleCartesian(Particle * pl)
{
    double temp[DIMENSIONS];
    int i;
    transformSphericalCartesian(temp,pl->mom,pl->coor[0],pl->coor[1],pl->coor[2]);
    for(i=0;i<DIMENSIONS;i++)
        pl->mom[i] = temp[i];
#if 0
    transformSphericalCartesian(temp,pl->acceleration,pl->coor[0],pl->coor[1],pl->coor[2]);
    for(i=0;i<DIMENSIONS;i++)
        pl->acceleration[i] = temp[i];
#endif
}
