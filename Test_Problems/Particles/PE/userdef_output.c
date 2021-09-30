#include "pluto.h"
#include "irradiation.h"

#define UNIT_MASS (UNIT_DENSITY * UNIT_LENGTH * UNIT_LENGTH * UNIT_LENGTH)
#define UNIT_G (CONST_G / (UNIT_VELOCITY*UNIT_VELOCITY) / UNIT_LENGTH * UNIT_MASS) // gravitational constant in code units


/* *************************************************************** */
void ComputeUserVar (const Data *d, Grid *grid)
/*
 *
 *  PURPOSE
 *
 *    Define user-defined output variables
 *
 *
 *
 ***************************************************************** */
{
  int i, j, k;
  double ***cd, ***Td, ***Tg;

  cd = GetUserVar("cd");
  Td = GetUserVar("Tdust");
  Tg = GetUserVar("Tgas");

  DOM_LOOP(k,j,i){
    cd[k][j][i] = irradiation.S[k][j][i];
    Td[k][j][i] = irradiation.Tdust[k][j][i];
    Tg[k][j][i] = irradiation.Tgas[k][j][i];
  }
}


/* ************************************************************* */
void ChangeOutputVar ()
/* 
 *
 * 
 *************************************************************** */
{ 
  Image *image;

#if PARTICLES
  //SetOutputVar ("energy",PARTICLES_FLT_OUTPUT, NO);
//  SetOutputVar ("x1",    PARTICLES_FLT_OUTPUT, NO);
  //SetOutputVar ("vx1",   PARTICLES_FLT_OUTPUT, NO);
#endif

}
