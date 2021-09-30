#ifndef MLW_H_PLANET
#define MLW_H_PLANET

#define SML16 1

#ifndef PLANET_HEATING
  #define PLANET_HEATING NO
#endif

#define OmegaP ( 2*CONST_PI * sqrt(1.0 + g_inputParam[MP] / g_inputParam[MS] * CONST_Mearth / CONST_Msun))
#define PlanetM0  ( g_inputParam[PLANETRAMPUPM0] * CONST_Mearth / UNIT_MASS )

/* structs */
typedef struct PLANET_DATA
{
  double x;
  double y;
  double rhill;
  double M;
  double M_accreted;
} PlanetData;

extern PlanetData planet;

#endif