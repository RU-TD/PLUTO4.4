#include "pluto.h"

/* Indirect term: acceleration of the centre of mass of the central object */
void nbodyCalcIndirectTerm()
{
    double inv_m = 1.0/g_inputParam[M_CO];
    g_nb.ax_indirect = 0.0;
    g_nb.ay_indirect = 0.0;
    g_nb.az_indirect = 0.0;

    int l;
    for (l = 0; l < CENTRAL_OBJECT; l++)
    {
        DIM_EXPAND(g_nb.ax_indirect += g_nb.m[l]*g_nb.ax[l];,
                   g_nb.ay_indirect += g_nb.m[l]*g_nb.ay[l];,
                   g_nb.az_indirect += g_nb.m[l]*g_nb.az[l];)
    }

    DIM_EXPAND(g_nb.ax_indirect *= inv_m;,
               g_nb.ay_indirect *= inv_m;,
               g_nb.az_indirect *= inv_m;)
}

/* Calculate acceleration of each body due to all other bodies */
void nbodyCalcAccelerations()
{
    /* Resets all accelerations to zero */
    DIM_EXPAND(g_nb.ax = memset(g_nb.ax, 0.0, NB_N*sizeof(double));,
               g_nb.ay = memset(g_nb.ay, 0.0, NB_N*sizeof(double));,
               g_nb.az = memset(g_nb.az, 0.0, NB_N*sizeof(double));)

    int i,j;
	for (i = 0; i < NB_N; i++)
	{
		for (j = i+1; j < NB_N; j++)
		{
            DIM_EXPAND(double dx = g_nb.x[i] - g_nb.x[j];,
			           double dy = g_nb.y[i] - g_nb.y[j];,
			           double dz = g_nb.z[i] - g_nb.z[j];)

            double r3 = DIM_EXPAND(dx*dx, +dy*dy, +dz*dz);
            r3 *= sqrt(r3);

			double acc_mag_i = CONST_G_CODE_UNITS*g_nb.m[j] / r3;
			double acc_mag_j = CONST_G_CODE_UNITS*g_nb.m[i] / r3;

            DIM_EXPAND(g_nb.ax[i] += -acc_mag_i * dx;
                       g_nb.ax[j] +=  acc_mag_j * dx;, 
                       g_nb.ay[i] += -acc_mag_i * dy;
                       g_nb.ay[j] +=  acc_mag_j * dy;,
                       g_nb.az[i] += -acc_mag_i * dz;
                       g_nb.az[j] +=  acc_mag_j * dz;)
		}
	}

    int l;
    for (l = 0; l < NB_N; l++)
    {
        /* Disk feedback */
        DIM_EXPAND(g_nb.ax[l] += g_nb.axdisk[l];,
                   g_nb.ay[l] += g_nb.aydisk[l];,
                   g_nb.az[l] += g_nb.azdisk[l];)
    }
}

void nbodyCalcDiskFeedback(const Data *d, Grid *grid)
{
    int nv;
    double vc[NVAR];
    int l;
    double *x1 = grid->x[IDIR];
    double *x2 = grid->x[JDIR];
    double *x3 = grid->x[KDIR];

    for (l = 0; l < NB_N; l++)
    {
        /* reseting accelerations */
        g_nb.axdisk[l] = 0.0;
        g_nb.aydisk[l] = 0.0;
        g_nb.azdisk[l] = 0.0;

        /* check if body l feels the disk */
        if (g_nb.feelsDisk[l])
        {
            /* calculation of disk feedback on body l */
            int i, j, k;
            KDOM_LOOP(k) 
            {
                #if GEOMETRY == SPHERICAL
                double sin_phi = sin(x3[k]);
                double cos_phi = cos(x3[k]);
                #endif
                JDOM_LOOP(j)
                {
                    #if GEOMETRY == POLAR
                    double sin_phi = sin(x2[j]);
                    double cos_phi = cos(x2[j]);
                    #elif GEOMETRY == SPHERICAL
                    double sin_theta = sin(x2[j]);
                    double cos_theta = cos(x2[j]);
                    #endif
                    
                    IDOM_LOOP(i)
                    {
                        #if GEOMETRY == POLAR
                        double R = x1[i];
                        #elif GEOMETRY == SPHERICAL
                        double r = x1[i];
                        #endif

			            double dV = grid->dV[k][j][i];

                        NFLX_LOOP(nv) vc[nv] = d->Vc[nv][k][j][i];
                        double smoothing_squared = nbodySmoothingSquared(
                                                      vc, 
                                                      x1[i],
                                                      x2[j],
                                                      x3[k]);

                        #if GEOMETRY == POLAR
                        DIM_EXPAND(double xc = R * cos_phi;,
                                   double yc = R * sin_phi;,
                                   double zc = x3[k];)
                        #elif GEOMETRY == SPHERICAL
                        DIM_EXPAND(double xc = r*sin_theta*cos_phi;,
                                   double yc = r*sin_theta*sin_phi;,
                                   double zc = r*cos_theta;)
                        #endif

                        DIM_EXPAND(double dx = xc - g_nb.x[l];,
                                   double dy = yc - g_nb.y[l];,
                                   double dz = zc - g_nb.z[l];)

                        double r3 = DIM_EXPAND(dx*dx, +dy*dy, +dz*dz) + smoothing_squared;
                        r3 *= sqrt(r3);

                        double magnitude =  CONST_G_CODE_UNITS 
                                           *d->Vc[RHO][k][j][i] * dV / r3;

                        DIM_EXPAND(g_nb.axdisk[l] += magnitude*dx;,
                                   g_nb.aydisk[l] += magnitude*dy;,
                                   g_nb.azdisk[l] += magnitude*dz;)
                    }
                }
            }

            /* MPI communication */
            #ifdef PARALLEL
            double tmp;
            DIM_EXPAND(
            MPI_Allreduce(g_nb.axdisk+l, &tmp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            g_nb.axdisk[l] = tmp;,
            MPI_Allreduce(g_nb.aydisk+l, &tmp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            g_nb.aydisk[l] = tmp;,
            MPI_Allreduce(g_nb.azdisk+l, &tmp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            g_nb.azdisk[l] = tmp;)
            #endif
        }
    }
}

/* calculates acceleration of cell (x1, x2, x3) due to all bodies */
double nbodyCalcCellAcceleration(double *v,
                                 double x1,
                                 double x2,
                                 double x3,
                                 int dir)
{
    #if GEOMETRY == POLAR
    double sin_phi = sin(x2);
    double cos_phi = cos(x2);
    #elif GEOMETRY == SPHERICAL
    double sin_theta = sin(x2);
    double cos_theta = cos(x2);
    double sin_phi = sin(x3);
    double cos_phi = cos(x3);
    #endif

    #if GEOMETRY == POLAR
    DIM_EXPAND(double x = x1*cos_phi;,
           double y = x1*sin_phi;,
           double z = x3;)
    #elif GEOMETRY == SPHERICAL
    DIM_EXPAND(double x = x1*sin_theta*cos_phi;,
           double y = x1*sin_theta*sin_phi;,
           double z = x1*cos_theta;)
    #endif

    double smoothing_squared = nbodySmoothingSquared(v, x1, x2, x3);

    double ax = 0.0;
    double ay = 0.0;
    double az = 0.0;

    int l;
    for (l = 0; l < NB_N; l++)
    {
        DIM_EXPAND(double dx = g_nb.x[l] - x;,
                   double dy = g_nb.y[l] - y;,
                   double dz = g_nb.z[l] - z;)

        double r3 = DIM_EXPAND(dx*dx, +dy*dy, +dz*dz) + smoothing_squared;
        r3 *= sqrt(r3);
        double magnitude = CONST_G_CODE_UNITS*g_nb.m[l]/r3;

        DIM_EXPAND(ax += magnitude*dx;,
                   ay += magnitude*dy;,
                   az += magnitude*dz;)
    }
    
    #if INDIRECT_TERM == YES
    DIM_EXPAND(ax -= g_nb.ax_indirect;,
               ay -= g_nb.ay_indirect;,
               az -= g_nb.az_indirect;)
    #endif

    DIM_EXPAND(
    if (dir == IDIR)
    {
        #if GEOMETRY == POLAR
        return ax*cos_phi + ay*sin_phi;
        #elif GEOMETRY == SPHERICAL
        return ax*sin_theta*cos_phi + ay*sin_theta*sin_phi + az*cos_theta;
        #endif
    },
    else if (dir == JDIR)
    {
        #if GEOMETRY == POLAR
        return -ax*sin_phi + ay*cos_phi;
        #elif GEOMETRY == SPHERICAL
        return ax*cos_theta*cos_phi/x1 + ay*cos_theta*sin_phi/x1 - az*sin_theta/x1; 
        #endif
    },
    else
    {
        #if GEOMETRY == POLAR
        return a_z;
        #elif GEOMETRY == SPHERICAL
        return -ax*sin_phi/(x1*sin_theta) + ay*cos_phi/(x1*sin_theta);
        #endif
    })
}
