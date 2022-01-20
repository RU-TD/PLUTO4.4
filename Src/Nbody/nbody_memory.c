#include "pluto.h"

double *k1x, *k1y, *k1z;
double *k1vx, *k1vy, *k1vz;

double *k2x, *k2y, *k2z;
double *k2vx, *k2vy, *k2vz;

double *k3x, *k3y, *k3z;
double *k3vx, *k3vy, *k3vz;

double *k4x, *k4y, *k4z;
double *k4vx, *k4vy, *k4vz;

void nbodyAllocateMemory()
{
    /* memory allocation */
    g_nb.m = ARRAY_1D(NB_N, double);
    g_nb.x = ARRAY_1D(NB_N, double);
    g_nb.y = ARRAY_1D(NB_N, double);
    g_nb.z = ARRAY_1D(NB_N, double);

    g_nb.vx = ARRAY_1D(NB_N, double);
    g_nb.vy = ARRAY_1D(NB_N, double);
    g_nb.vz = ARRAY_1D(NB_N, double);

    g_nb.ax = ARRAY_1D(NB_N, double);
    g_nb.ay = ARRAY_1D(NB_N, double);
    g_nb.az = ARRAY_1D(NB_N, double);

    g_nb.xold = ARRAY_1D(NB_N, double);
    g_nb.yold = ARRAY_1D(NB_N, double);
    g_nb.zold = ARRAY_1D(NB_N, double);

    g_nb.vxold = ARRAY_1D(NB_N, double);
    g_nb.vyold = ARRAY_1D(NB_N, double);
    g_nb.vzold = ARRAY_1D(NB_N, double);

    g_nb.axdisk = ARRAY_1D(NB_N, double);
    g_nb.aydisk = ARRAY_1D(NB_N, double);
    g_nb.azdisk = ARRAY_1D(NB_N, double);

    g_nb.feelsDisk = ARRAY_1D(NB_N, int);

    g_nb.a = ARRAY_1D(NB_N, double);
    g_nb.e = ARRAY_1D(NB_N, double);
    g_nb.inc = ARRAY_1D(NB_N, double);
    g_nb.Omega = ARRAY_1D(NB_N, double);
    g_nb.omega = ARRAY_1D(NB_N, double);
    g_nb.f = ARRAY_1D(NB_N, double);
    g_nb.P = ARRAY_1D(NB_N, double);
    g_nb.ea = ARRAY_1D(NB_N, double);
    g_nb.ma = ARRAY_1D(NB_N, double);
}

void nbodyFreeMemory()
{
    FreeArray1D(g_nb.x);
    FreeArray1D(g_nb.y);
    FreeArray1D(g_nb.z);

    FreeArray1D(g_nb.vx);
    FreeArray1D(g_nb.vy);
    FreeArray1D(g_nb.vz);

    FreeArray1D(g_nb.ax);
    FreeArray1D(g_nb.ay);
    FreeArray1D(g_nb.az);

    FreeArray1D(g_nb.xold);
    FreeArray1D(g_nb.yold);
    FreeArray1D(g_nb.zold);

    FreeArray1D(g_nb.vxold);
    FreeArray1D(g_nb.vyold);
    FreeArray1D(g_nb.vzold);

    FreeArray1D(g_nb.axdisk);
    FreeArray1D(g_nb.aydisk);
    FreeArray1D(g_nb.azdisk);

    FreeArray1D(g_nb.feelsDisk);

    FreeArray1D(g_nb.a);
    FreeArray1D(g_nb.e);
    FreeArray1D(g_nb.inc);
    FreeArray1D(g_nb.Omega);
    FreeArray1D(g_nb.omega);
    FreeArray1D(g_nb.f);
    FreeArray1D(g_nb.P);
    FreeArray1D(g_nb.ea);
    FreeArray1D(g_nb.ma);
}

void nbodyAllocateRK4Memory()
{
    k1x = ARRAY_1D(NB_N, double);
    k1y = ARRAY_1D(NB_N, double);
    k1z = ARRAY_1D(NB_N, double);
    k1vx = ARRAY_1D(NB_N, double);
    k1vy = ARRAY_1D(NB_N, double);
    k1vz = ARRAY_1D(NB_N, double);

    k2x = ARRAY_1D(NB_N, double);
    k2y = ARRAY_1D(NB_N, double);
    k2z = ARRAY_1D(NB_N, double);
    k2vx = ARRAY_1D(NB_N, double);
    k2vy = ARRAY_1D(NB_N, double);
    k2vz = ARRAY_1D(NB_N, double);

    k3x = ARRAY_1D(NB_N, double);
    k3y = ARRAY_1D(NB_N, double);
    k3z = ARRAY_1D(NB_N, double);
    k3vx = ARRAY_1D(NB_N, double);
    k3vy = ARRAY_1D(NB_N, double);
    k3vz = ARRAY_1D(NB_N, double);

    k4x = ARRAY_1D(NB_N, double);
    k4y = ARRAY_1D(NB_N, double);
    k4z = ARRAY_1D(NB_N, double);
    k4vx = ARRAY_1D(NB_N, double);
    k4vy = ARRAY_1D(NB_N, double);
    k4vz = ARRAY_1D(NB_N, double);
}

void nbodyFreeRK4Memory()
{
    FreeArray1D(k1x);
    FreeArray1D(k1y);
    FreeArray1D(k1z);
    FreeArray1D(k1vx);
    FreeArray1D(k1vy);
    FreeArray1D(k1vz);

    FreeArray1D(k2x);
    FreeArray1D(k2y);
    FreeArray1D(k2z);
    FreeArray1D(k2vx);
    FreeArray1D(k2vy);
    FreeArray1D(k2vz);
    
    FreeArray1D(k3x);
    FreeArray1D(k3y);
    FreeArray1D(k3z);
    FreeArray1D(k3vx);
    FreeArray1D(k3vy);
    FreeArray1D(k3vz);
    
    FreeArray1D(k4x);
    FreeArray1D(k4y);
    FreeArray1D(k4z);
    FreeArray1D(k4vx);
    FreeArray1D(k4vy);
    FreeArray1D(k4vz);
}
