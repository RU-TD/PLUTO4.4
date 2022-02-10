#include "pluto.h"

void nbodySaveOldState()
{
    int l;
    for (l = 0; l < NB_N; l++)
    {
        g_nb.xold[l] = g_nb.x[l];
        g_nb.yold[l] = g_nb.y[l];
        g_nb.zold[l] = g_nb.z[l];

        g_nb.vxold[l] = g_nb.vx[l];
        g_nb.vyold[l] = g_nb.vy[l];
        g_nb.vzold[l] = g_nb.vz[l];
    }
}

/* Mass of the first n particles */
double nbodyMass(int n)
{
    double m = 0.0;
    int l;
    for (l = 0; l < n; l++)
    {
        m += g_nb.m[l];
    }

	return m;
}

/* Coordinates and velocity of the  center of mass of the first n particles */
void nbodyCenterOfMassCoordinatesAndVelocity(int n)
{
    double total_mass = 0.0;
    g_nb.comX = 0.0;
    g_nb.comY = 0.0;
    g_nb.comZ = 0.0;
    
    g_nb.comVX = 0.0;
    g_nb.comVY = 0.0;
    g_nb.comVZ = 0.0;

    int l;
    for (l = 0; l < n; l++)
    {
        total_mass += g_nb.m[l];
        g_nb.comX += g_nb.m[l] * g_nb.x[l];
        g_nb.comY += g_nb.m[l] * g_nb.y[l];
        g_nb.comZ += g_nb.m[l] * g_nb.z[l];

        g_nb.comVX += g_nb.m[l] * g_nb.vx[l];
        g_nb.comVY += g_nb.m[l] * g_nb.vy[l];
        g_nb.comVZ += g_nb.m[l] * g_nb.vz[l];
    }

    g_nb.comX /= total_mass;
    g_nb.comY /= total_mass;
    g_nb.comZ /= total_mass;

    g_nb.comVX /= total_mass;
    g_nb.comVY /= total_mass;
    g_nb.comVZ /= total_mass;
}

/* Orbital elements of particles 2..N in Jacobian coordinates:
   This means the orbital elements of particle n are calculated with 
   respect to the center of mass of the first 1..n-1 particles. */
void nbodyCalcOrbitalElements()
{
    int l;
    for (l = 1; l < NB_N; l++)
    {
        double M = nbodyMass(l);
        nbodyCenterOfMassCoordinatesAndVelocity(l);

        double inv_mu = 1.0/(CONST_G_CODE_UNITS*(M+g_nb.m[l]));

        /* Transform to Jacobian coordinates */
        double x = g_nb.x[l] - g_nb.comX;
        double y = g_nb.y[l] - g_nb.comY;
        double z = g_nb.z[l] - g_nb.comZ;

        double vx = g_nb.vx[l] - g_nb.comVX;
        double vy = g_nb.vy[l] - g_nb.comVY;
        double vz = g_nb.vz[l] - g_nb.comVZ;

        double r = sqrt(x*x + y*y + z*z);
        double inv_r = 1.0/r;
        double v_squared_mu = (vx*vx + vy*vy + vz*vz)*inv_mu;
        double rv = x*vx + y*vy + z*vz;

        /* Eccentricity vector */
        double ex = (v_squared_mu - 1.0*inv_r) * x - rv*inv_mu * vx;
        double ey = (v_squared_mu - 1.0*inv_r) * y - rv*inv_mu * vy;
        double ez = (v_squared_mu - 1.0*inv_r) * z - rv*inv_mu * vz;

        /* Specific angular momentum vector */
        double hx = y*vz - z*vy;
        double hy = z*vx - x*vz;
        double hz = x*vy - y*vx;
        double h = sqrt(hx*hx + hy*hy + hz*hz);

        /* This vector points to the ascending node */
        double nx = -hy;
        double ny = hx;
        double n = sqrt(nx*nx + ny*ny);

        /* Semi-major axis */
        g_nb.a[l] = 1.0/(2.0*inv_r - v_squared_mu);

        /* Period */
        g_nb.P[l] = 2.0*CONST_PI*sqrt(g_nb.a[l]*g_nb.a[l]*g_nb.a[l]*inv_mu);

        /* Eccentricity */
        g_nb.e[l] = sqrt(ex*ex + ey*ey + ez*ez);

        /* Inclination */
        g_nb.inc[l] = acos(hz/h);

        /* Planar orbits (i = 0 or i = pi) */
        if (g_nb.inc[l] < ORBITAL_EPS || g_nb.inc[l] > CONST_PI - ORBITAL_EPS)
        {
            /* Longitude of ascending node */
            g_nb.Omega[l] = 0.0;

            /* Circular orbits (e = 0) */
            if (g_nb.e[l] < ORBITAL_EPS)
            {
                /* Argument of periapsis */
                g_nb.omega[l] = 0.0;

                /* True anomaly */
                g_nb.f[l] = acos(x/r);
                if (vx > 0.0)
                {
                    g_nb.f[l] = 2.0*CONST_PI - g_nb.f[l];
                }
            }
            /* Elliptic orbits (e > 0) */
            else 
            {
                /* Argument of periapsis */
                g_nb.omega[l] = atan2(ey, ex);
                if (g_nb.omega[l] < 0.0)
                {
                    g_nb.omega[l] += 2.0*CONST_PI;
                }

                /* True anomaly */
                g_nb.f[l] = acos((ex*x + ey*y + ez*z)/(g_nb.e[l]*r));
                if (rv < 0.0)
                {
                    g_nb.f[l] = 2.0*CONST_PI - g_nb.f[l];
                }
            }
        }
        /* Inclined orbits (0 < i < pi) */
        else
        {
            /* Longitude of ascending node */
            g_nb.Omega[l] = acos(nx/n);
            if (ny < 0.0)
            {
                g_nb.Omega[l] = 2.0*CONST_PI - g_nb.Omega[l];
            }

            /* Circular orbits (e = 0) */
            if (g_nb.e[l] < ORBITAL_EPS)
            {
                /* Argument of periapsis */
                g_nb.omega[l-1] = 0.0;
                
                /* True anomaly */
                g_nb.f[l] = acos((nx*x + ny*y)/(n*r));
                if (z < 0.0)
                {
                    g_nb.f[l] = 2.0*CONST_PI - g_nb.f[l];
                }
            }
            /* Elliptic orbits (e > 0) */
            else
            {
                /* Argument of periapsis */
                g_nb.omega[l] = acos((nx*ex + ny*ey)/(n*g_nb.e[l]));
                if (ez < 0.0)
                {
                    g_nb.omega[l] = 2.0*CONST_PI - g_nb.omega[l];
                }

                /* True anomaly */
                g_nb.f[l] = acos((ex*x + ey*y + ez*z)/(g_nb.e[l]*r));
                if (rv < 0.0)
                {
                    g_nb.f[l] = 2.0*CONST_PI - g_nb.f[l];
                }
            }
        }
        
        /* Eccentric anomaly */
        g_nb.ea[l] = acos((g_nb.e[l] + cos(g_nb.f[l]))/(1.0 + g_nb.e[l]*cos(g_nb.f[l])));
        if (g_nb.f[l] > CONST_PI)
        {
            g_nb.ea[l] =  2.0*CONST_PI -g_nb.ea[l];
        }

        /* Mean anomaly */
        g_nb.ma[l] = g_nb.ea[l] - g_nb.e[l]*sin(g_nb.ea[l]);
    }
}
