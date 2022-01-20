#include "pluto.h"

/* Advance nbody system in time from t to t+dt using a RK4 integrator */
void nbodyAdvanceSystem(double dt)
{
    double one_third = 1.0/3.0;
    double one_sixth = 1.0/6.0;

    nbodySaveOldState();

    /* calculation of k1 */
    nbodyCalcAccelerations();
    int l;
    for (l = 0; l < NB_N; l++)
    {
        DIM_EXPAND(k1x[l] = g_nb.vx[l]*dt;,
               k1y[l] = g_nb.vy[l]*dt;,
               k1z[l] = g_nb.vz[l]*dt;)

        DIM_EXPAND(k1vx[l] = g_nb.ax[l]*dt;,
               k1vy[l] = g_nb.ay[l]*dt;,
               k1vz[l] = g_nb.az[l]*dt;)
    }

    /* calculation of k2 */
    for (l = 0; l < NB_N; l++)
    {
        DIM_EXPAND(g_nb.x[l] = g_nb.xold[l] + 0.5*k1x[l];,
               g_nb.y[l] = g_nb.yold[l] + 0.5*k1y[l];, 
               g_nb.z[l] = g_nb.zold[l] + 0.5*k1z[l];)

        DIM_EXPAND(g_nb.vx[l] = g_nb.vxold[l] + 0.5*k1vx[l];, 
               g_nb.vy[l] = g_nb.vyold[l] + 0.5*k1vy[l];, 
               g_nb.vz[l] = g_nb.vzold[l] + 0.5*k1vz[l];)
    }
    nbodyCalcAccelerations();
    for (l = 0; l < NB_N; l++)
    {
        DIM_EXPAND(k2x[l] = g_nb.vx[l]*dt;,
               k2y[l] = g_nb.vy[l]*dt;,
               k2z[l] = g_nb.vz[l]*dt;)

        DIM_EXPAND(k2vx[l] = g_nb.ax[l]*dt;,
               k2vy[l] = g_nb.ay[l]*dt;,
               k2vz[l] = g_nb.az[l]*dt;)
    }

    /* calculation of k3 */
    for (l = 0; l < NB_N; l++)
    {
        DIM_EXPAND(g_nb.x[l] = g_nb.xold[l] + 0.5*k2x[l];,
               g_nb.y[l] = g_nb.yold[l] + 0.5*k2y[l];, 
               g_nb.z[l] = g_nb.zold[l] + 0.5*k2z[l];)

        DIM_EXPAND(g_nb.vx[l] = g_nb.vxold[l] + 0.5*k2vx[l];,
               g_nb.vy[l] = g_nb.vyold[l] + 0.5*k2vy[l];, 
               g_nb.vz[l] = g_nb.vzold[l] + 0.5*k2vz[l];)
    }
    nbodyCalcAccelerations();
    for (l = 0; l < NB_N; l++)
    {
        DIM_EXPAND(k3x[l] = g_nb.vx[l]*dt;,
               k3y[l] = g_nb.vy[l]*dt;,
               k3z[l] = g_nb.vz[l]*dt;)

        DIM_EXPAND(k3vx[l] = g_nb.ax[l]*dt;,
               k3vy[l] = g_nb.ay[l]*dt;,
               k3vz[l] = g_nb.az[l]*dt;)
    }

    /* calculation of k4 */
    for (l = 0; l < NB_N; l++)
    {
        DIM_EXPAND(g_nb.x[l] = g_nb.xold[l] + k3x[l];, 
               g_nb.y[l] = g_nb.yold[l] + k3y[l];, 
               g_nb.z[l] = g_nb.zold[l] + k3z[l];) 

        DIM_EXPAND(g_nb.vx[l] = g_nb.vxold[l] + k3vx[l];, 
               g_nb.vy[l] = g_nb.vyold[l] + k3vy[l];, 
               g_nb.vz[l] = g_nb.vzold[l] + k3vz[l];) 
    }
    nbodyCalcAccelerations();
    for (l = 0; l < NB_N; l++)
    {
        DIM_EXPAND(k4x[l] = g_nb.vx[l]*dt;,
               k4y[l] = g_nb.vy[l]*dt;,
               k4z[l] = g_nb.vz[l]*dt;)

        DIM_EXPAND(k4vx[l] = g_nb.ax[l]*dt;,
               k4vy[l] = g_nb.ay[l]*dt;,
               k4vz[l] = g_nb.az[l]*dt;)
    }

    /* calculate solution at t+dt */
    for (l = 0; l < NB_N; l++)
    {
        DIM_EXPAND(g_nb.x[l] = g_nb.xold[l] + (k1x[l] + k4x[l])*one_sixth 
                                        + (k2x[l] + k3x[l])*one_third;,
               g_nb.y[l] = g_nb.yold[l] + (k1y[l] + k4y[l])*one_sixth
                                        + (k2y[l] + k3y[l])*one_third;,
               g_nb.z[l] = g_nb.zold[l] + (k1z[l] + k4z[l])*one_sixth 
                                        + (k2z[l] + k3z[l])*one_third;)

        DIM_EXPAND(g_nb.vx[l] = g_nb.vxold[l] + (k1vx[l] + k4vx[l])*one_sixth 
                                          + (k2vx[l] + k3vx[l])*one_third;,
               g_nb.vy[l] = g_nb.vyold[l] + (k1vy[l] + k4vy[l])*one_sixth
                                          + (k2vy[l] + k3vy[l])*one_third;,
               g_nb.vz[l] = g_nb.vzold[l] + (k1vz[l] + k4vz[l])*one_sixth 
                                          + (k2vz[l] + k3vz[l])*one_third;)
    }

    /* Move back to centre of mass of central object */
    nbodyCenterOfMassCoordinatesAndVelocity(CENTRAL_OBJECT);
    for (l = 0; l < NB_N; l++)
    {
        g_nb.x[l] -= g_nb.comX;
        g_nb.y[l] -= g_nb.comY;
        g_nb.z[l] -= g_nb.comZ;
        g_nb.vx[l] -= g_nb.comVX;
        g_nb.vy[l] -= g_nb.comVY;
        g_nb.vz[l] -= g_nb.comVZ;
    }
}
