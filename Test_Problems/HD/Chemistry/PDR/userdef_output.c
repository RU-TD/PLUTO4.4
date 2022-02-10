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
    //double ***cd, ***cdh2, ***cdco, ***cdh2o, ***cdtot;
    double ***cd, ***cdh2, ***cdtot;//, ***cdh2o;// ***cdtot;
    cd = GetUserVar("cd");
    cdh2 = GetUserVar("cdh2");
    //cdco = GetUserVar("cdco");
    //cdh2o = GetUserVar("cdh2o");
    cdtot = GetUserVar("cdtot");

    DOM_LOOP(k,j,i){
        cd[k][j][i] = irradiation.column_density[k][j][i][0];
	cdh2[k][j][i] = irradiation.column_density[k][j][i][1];
	//cdco[k][j][i] = irradiation.column_density[k][j][i][2];
	//cdh2o[k][j][i] = irradiation.column_density[k][j][i][3];
	cdtot[k][j][i] = irradiation.column_density[k][j][i][4];
    }
}
/* ************************************************************* */
void ChangeDumpVar ()
/* 
 *
 * 
 *************************************************************** */
{ 
  Image *image;

}





