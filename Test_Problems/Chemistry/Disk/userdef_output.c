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
    double ***cd, ***cdH2, ***cdCO;
    cd = GetUserVar("cd");
    cdH2 = GetUserVar("cdH2");
    cdCO = GetUserVar("cdCO");

    DOM_LOOP(k,j,i){
        cd[k][j][i] = irradiation.column_density[0][k][j][i];
	cdH2[k][j][i] = irradiation.column_density[1][k][j][i];
	cdCO[k][j][i] = irradiation.column_density[2][k][j][i];
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





