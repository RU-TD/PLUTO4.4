#ifndef MLW_H_INIT
#define MLW_H_INIT

#include "pluto.h"
#include "irradiation.h"

#define T_MIN 10.0
#define T_MAX 1.e4
#define RHO_MIN ( 1e-21 / UNIT_DENSITY )

#define UNIT_MASS (UNIT_DENSITY * UNIT_LENGTH * UNIT_LENGTH * UNIT_LENGTH)
#define UNIT_G (CONST_G / POW2(UNIT_VELOCITY) / UNIT_LENGTH * UNIT_MASS) // gravitational constant in code units

#define damping_rout (g_inputParam[DAMPING_ROUT] * CONST_au / UNIT_LENGTH)

#if ROTATING_FRAME == NO
  #define g_OmegaZ 0.0
#endif

#endif