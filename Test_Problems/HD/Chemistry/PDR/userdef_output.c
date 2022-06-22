#include "pluto.h"

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
    double ***cd, ***cdh2, ***cdtot;
    cd = GetUserVar("cd");
    cdh2 = GetUserVar("cdh2");
    cdtot = GetUserVar("cdtot");

    DOM_LOOP(k,j,i){
        cd[k][j][i] = irradiation.column_density[0][k][j][i];
	cdh2[k][j][i] = irradiation.column_density[1][k][j][i];
	cdtot[k][j][i] = irradiation.column_density[2][k][j][i];
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

}





