/* /////////////////////////////////////////////////////////////////// */
/*! \file
 *  \brief Specification of explicit first and second viscosity coefficients*/
/* /////////////////////////////////////////////////////////////////// */
#include "pluto.h"
/* ************************************************************************** */

void Visc_nu(double *v, double x1, double x2, double x3,
             double *nu1, double *nu2)
/*!
 *
 *  \param [in]      v  pointer to data array containing cell-centered quantities
 *  \param [in]      x1 real, coordinate value
 *  \param [in]      x2 real, coordinate value
 *  \param [in]      x3 real, coordinate value
 *  \param [in, out] nu1  pointer to first viscous coefficient
 *  \param [in, out] nu2  pointer to second viscous coefficient
 *
 *  \return This function has no return value.
 * ************************************************************************** */
{
  double R;

#if GEOMETRY == SPHERICAL
  R = x1 * sin(x2);
#elif GEOMETRY == CYLINDRICAL || (GEOMETRY == POLAR) 
  R = x1;
#endif

  // OmegaK = 2.0 * CONST_PI / (R * sqrt(R)); // sqrt(G * M_star) = 2*pi in code units
  // *nu1 = g_inputParam[ALPHA] * g_gamma * v[PRS] / OmegaK; // cs^2 = gamma * prs/rho
  *nu1 = g_inputParam[CONST_VISCOSITY] * 2.0 * CONST_PI * sqrt(R) * v[RHO]; // visc * R**2 * OmegaK
  *nu2 = 0.0;
}
