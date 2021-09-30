#include "pluto.h"
#include "particles.h"


int ParticleOutflowBoundary(Particle *pl, double ***uu[], Grid *grid, int dim, int side);
void ParticleReflectiveBoundary(Particle *pl, double ***uu[], Grid *grid, int dim, int side);
void ParticlePeriodicBoundary(Particle *pl, double ***uu[], Grid *grid, int dim, int side);
int ParticleUserdefBoundary(Particle *pl, double ***uu[], Grid *grid, int dim, int side);


/* **************************************************** */
int ParticlesBoundary(double ***uu[], Grid *grid, Particle * pl)
/*
 * Checks for boundary conditions and calls the relevant
 * routines in a way similar to the rest of the code.
 *
 * It returns 1 if the particle is still inside the domain
 * or 0 if the particle has to be destroyed (done with
 * specific instructions at the end of the loop on particles)
 *
 ********************************************************/
{

    int dim, side;
    int boundaryType;
    int check = 1;

    for (dim = 0; dim < DIMENSIONS; ++dim) {       
        boundaryType = particlesConf.left_bound[dim]; // treat the left boundary first
        for (side = 0; side < 2; ++side) {

            #if DEBUG == YES
                print("PRANK: %d, dim: %d, side: %d - this_side_type: %d\n", prank, dim, side, boundaryType);
            #endif

            switch(boundaryType) {
                case OUTFLOW:
                    check = ParticleOutflowBoundary(pl, uu, grid, dim, side);
                    break;

                case REFLECTIVE:
                case AXISYMMETRIC:
                case EQTSYMMETRIC:
                    ParticleReflectiveBoundary(pl, uu, grid, dim, side);
                    break;

                case PERIODIC:
                    ParticlePeriodicBoundary(pl, uu, grid, dim, side);
                    break;

                case USERDEF:
                    check = ParticleUserdefBoundary(pl, uu, grid, dim, side);
                    break;

                default:
                    print("ERROR: Boundary condition %i for dimension %i side %i is undefined.\n", 
                        boundaryType, dim, side);
                    QUIT_PLUTO(1);
            }
            if (check == 0) return 0;

            boundaryType = particlesConf.right_bound[dim]; // in the next loop treat the right boundary
        }
    }
    return 1;
}


/* **************************************************** */
int ParticleUserdefBoundary(Particle *pl, double ***uu[], Grid *grid,
        int dim, int side)

    /* By default userdefined boundaries are set to discriminate between
     * particles flowing from inside towards OUTSIDE of the domain,
     * and viceversa.
     *
     * TO SET INFLOW CONDITIONS, LEAVE THE INFLOWING PARTICLES ALONE HERE
     * AND SET CONDITIONS IN PART_INFLOW, BELOW.
     * side is either 0 or 1, so k=2*side-1 is either -1 or 1
     * In this way, k*pl->mom[dim] is positive for outflow,
     * and negative for inflow, regardless of the side value.
     * This allows to use just one if statement and not one per side
     * returns 1 if the particle is still inside the domain
     * and returns 0 if the particle must be destroyed,
     * using the integer 'check'
     *
     * **************************************************** */
{
    int check = 1;

    if ((2 * side - 1) * pl->mom[dim] > 0.0) { /* outflowing particles */

        /* customize or leave blank */

    } else { /* inflowing particles */

        /* customize or leave blank */

    }

    return check;

}


/* **************************************************** */
int ParticleOutflowBoundary(Particle *pl, double ***uu[], Grid *grid,
        int dim, int side)
    /*
     * returns 1 if the particle is still inside the domain
     * and returns 0 if the particle must be destroyed
     *
     * **************************************************** */
{
    if (side == 0) {
        return pl->coor[dim] >= grid->xbeg_glob[dim] ? 1 : 0;
    } else {
        return pl->coor[dim] <= grid->xend_glob[dim] ? 1 : 0;
    }
}

   
/* **************************************************** */
void ParticleReflectiveBoundary(Particle *pl, double ***uu[], Grid *grid,
        int dim, int side)
    /*
     * Refelective Boundary conditions for Particles.
     *
     * **************************************************** */
{
    pl->mom[dim] *= -1.0;
    if (side == 0) {
        if (pl->coor[dim] < grid->xbeg_glob[dim]) {
            pl->coor[dim] = 2 * grid->xbeg_glob[dim] - pl->coor[dim];
            pl->coor_old[dim] = 2 * grid->xbeg_glob[dim] - pl->coor_old[dim];
        }
    } else {
        if (pl->coor[dim] > grid->xend_glob[dim]) {
            pl->coor[dim] = 2 * grid->xend_glob[dim] - pl->coor[dim];
            pl->coor_old[dim] = 2 * grid->xend_glob[dim] - pl->coor_old[dim];
        }  
    }
}


/* **************************************************** */
void ParticlePeriodicBoundary(Particle *pl, double ***uu[], Grid *grid,
        int dim, int side)
    /*
     * Periodic Boundary conditions for Particles.
     *
     * **************************************************** */
{    
    if (side == 0) {
        while (pl->coor[dim] < grid->xbeg_glob[dim]) {
            pl->coor[dim] += grid->xend_glob[dim] - grid->xbeg_glob[dim];
            pl->coor_old[dim] += grid->xend_glob[dim] - grid->xbeg_glob[dim];
        }
    } else {
        while (pl->coor[dim] > grid->xend_glob[dim]) {
            pl->coor[dim] += grid->xbeg_glob[dim] - grid->xend_glob[dim];
            pl->coor_old[dim] = grid->xbeg_glob[dim] + pl->coor_old[dim] - grid->xend_glob[dim];
        }
    }
}
