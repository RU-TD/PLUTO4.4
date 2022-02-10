/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Restricted three-body problem

  \author D. Thun (daniel.thun@uni-tuebingen.de)
  \date   Jun 28, 2017

  \b References:
     - "Forward Symplectic Integrators for Solving Gravitational Few-Body 
        Problems"
        Chin & Chen, Celestial Mechanics and Dynamical Astronomy (2005)
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
void Init (double *v, double x1, double x2, double x3)
/*! 
 * The Init() function can be used to assign initial conditions as
 * as a function of spatial position.
 *
 * \param [out] v   a pointer to a vector of primitive variables
 * \param [in] x1   coordinate point in the 1st dimension
 * \param [in] x2   coordinate point in the 2nd dimension
 * \param [in] x3   coordinate point in the 3rdt dimension
 *
 * The meaning of x1, x2 and x3 depends on the geometry:
 * \f[ \begin{array}{cccl}
 *    x_1  & x_2    & x_3  & \mathrm{Geometry}    \\ \noalign{\medskip}
 *     \hline
 *    x    &   y    &  z   & \mathrm{Cartesian}   \\ \noalign{\medskip}
 *    R    &   z    &  -   & \mathrm{cylindrical} \\ \noalign{\medskip}
 *    R    & \phi   &  z   & \mathrm{polar}       \\ \noalign{\medskip}
 *    r    & \theta & \phi & \mathrm{spherical} 
 *    \end{array}
 *  \f]
 *
 * Variable names are accessed by means of an index v[nv], where
 * nv = RHO is density, nv = PRS is pressure, nv = (VX1, VX2, VX3) are
 * the three components of velocity, and so forth.
 *
 *********************************************************************** */
{
    g_isoSoundSpeed = g_inputParam[ASPECT_RATIO];
    g_smallDensity = g_inputParam[SIGMA_FLOOR]*g_inputParam[SIGMA_REF];

    v[RHO] = g_inputParam[SIGMA_REF] * pow(x1, -g_inputParam[ALPHA_SIGMA]);
    v[VX1] = 0.0;
    v[VX2] = sqrt(CONST_G_CODE_UNITS*g_inputParam[M_CO]/x1);
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
}


/* ********************************************************************* */
void Analysis (const Data *d, Grid *grid)
/*! 
 *  Perform runtime data analysis.
 *
 * \param [in] d the PLUTO Data structure
 * \param [in] grid   pointer to array of Grid structures  
 *
 *********************************************************************** */
{
}

/* ********************************************************************* */
void UserDefBoundary (const Data *d, RBox *box, int side, Grid *grid) 
/*! 
 *  Assign user-defined boundary conditions.
 *
 * \param [in,out] d  pointer to the PLUTO data structure containing
 *                    cell-centered primitive quantities (d->Vc) and 
 *                    staggered magnetic fields (d->Vs, when used) to 
 *                    be filled.
 * \param [in] box    pointer to a RBox structure containing the lower
 *                    and upper indices of the ghost zone-centers/nodes
 *                    or edges at which data values should be assigned.
 * \param [in] side   specifies the boundary side where ghost zones need
 *                    to be filled. It can assume the following 
 *                    pre-definite values: X1_BEG, X1_END,
 *                                         X2_BEG, X2_END, 
 *                                         X3_BEG, X3_END.
 *                    The special value side == 0 is used to control
 *                    a region inside the computational domain.
 * \param [in] grid  pointer to an array of Grid structures.
 *
 *********************************************************************** */
{
    int i, j, k;

    if (side == 0)
    {
        /* Density floor */
        TOT_LOOP(k,j,i) 
        {
            if (d->Vc[RHO][k][j][i] < g_smallDensity)
            {
                d->Vc[RHO][k][j][i] = g_smallDensity;
            }
        }
    } 
    else if (side == X1_BEG)
    {
        X1_BEG_LOOP(k, j, i)
        {
            /* drho /dr = 0 */
            d->Vc[RHO][k][j][i] = d->Vc[RHO][k][j][2*IBEG-i-1];

            if (d->Vc[VX1][k][j][IBEG] > 0.0)
                d->Vc[VX1][k][j][i] = -d->Vc[VX1][k][j][2*IBEG-i-1];
            else
                d->Vc[VX1][k][j][i] =  d->Vc[VX1][k][j][2*IBEG-i-1];

            /* domega / dr = 0 */
            d->Vc[VX2][k][j][i] = grid->x[IDIR][i] / grid->x[IDIR][2*IBEG-i-1]
                                 *d->Vc[VX2][k][j][2*IBEG-i-1];
        }
    }
    else if (side == X1_END)
    {
        X1_END_LOOP(k, j, i)
        {
            d->Vc[RHO][k][j][i] =  d->Vc[RHO][k][j][2*IEND-i+1];
            d->Vc[VX1][k][j][i] = -d->Vc[VX1][k][j][2*IEND-i+1];
            d->Vc[VX2][k][j][i] =  sqrt( CONST_G_CODE_UNITS*g_inputParam[M_CO]
                                        /grid->x[IDIR][i]);
        }
    }
}

#if (BODY_FORCE & VECTOR)
/* ********************************************************************* */
void BodyForceVector(double *v, double *g, double x1, double x2, double x3)
/*!
 * Prescribe the acceleration vector as a function of the coordinates
 * and the vector of primitive variables *v.
 *
 * \param [in] v  pointer to a cell-centered vector of primitive 
 *                variables
 * \param [out] g acceleration vector
 * \param [in] x1  position in the 1st coordinate direction \f$x_1\f$
 * \param [in] x2  position in the 2nd coordinate direction \f$x_2\f$
 * \param [in] x3  position in the 3rd coordinate direction \f$x_3\f$
 *
 *********************************************************************** */
{
    g[IDIR] = 0.0;
    g[JDIR] = 0.0;
    g[KDIR] = 0.0;
}
#endif

#if (BODY_FORCE & POTENTIAL)
/* ********************************************************************* */
double BodyForcePotential(double x1, double x2, double x3)
/*!
 * Return the gravitational potential as function of the coordinates.
 *
 * \param [in] x1  position in the 1st coordinate direction \f$x_1\f$
 * \param [in] x2  position in the 2nd coordinate direction \f$x_2\f$
 * \param [in] x3  position in the 3rd coordinate direction \f$x_3\f$
 * 
 * \return The body force potential \f$ \Phi(x_1,x_2,x_3) \f$.
 *
 *********************************************************************** */
{
    return 0.0;
}
#endif
