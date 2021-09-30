#ifndef MLW_H_INIT
#define MLW_H_INIT

#include "pluto.h"
#include "units.h"
#include "planet.h"
#include "irradiation.h"

#define T_MIN 10.0
#define T_MAX 1.e4
#define RHO_MIN ( 1e-21 / UNIT_DENSITY )

#define damping_rin (g_inputParam[DAMPING_RIN] * CONST_au / UNIT_LENGTH)
#define damping_rout (g_inputParam[DAMPING_ROUT] * CONST_au / UNIT_LENGTH)

#if ROTATING_FRAME == NO
  #define g_OmegaZ 0.0
#endif

/* structs */
typedef struct BOUNDARY_DATA
{
  double x1;
  double ** rho;
  double ** vx1;
  double ** vx2;
  double ** cd;
} BoundaryData;

typedef struct INITIAL_DATA
{
  double *** rho;
  double *** vx1;
  double *** vx2;
} InitialData;

extern BoundaryData innerBoundaryData;
extern InitialData initialData;

/* prototypes */
double calculateTorque(const Data *d, Grid *grid, double ***elemental_torques);

#endif