#include "pluto.h"

void nbodyInitializeCentralObject()
{
    double m = g_inputParam[M_CO];
    #if CENTRAL_OBJECT == STAR
    g_nb.m[0] = m;
    g_nb.x[0] = 0.0;
    g_nb.y[0] = 0.0;
    g_nb.z[0] = 0.0;

    g_nb.vx[0] = 0.0;
    g_nb.vy[0] = 0.0;
    g_nb.vz[0] = 0.0;

    g_nb.feelsDisk[0] = CO_FEELS_DISK;
    #elif CENTRAL_OBJECT == BINARY
    /*
     * This function calculates the position and the velocity
     * of both binary components, assuming they orbit around their centre
     * of mass. The reference plane is the x-y-plane and the reference
     * direction points along the x-axis. The orbital elements are the one
     * of the secondary. The following user-defined parameters are needed:
     *
     * M_CO:      Total mass of the binary (M_CO = M_primary+M_secondary)
     *            in code units
     * Q_BIN:     Mass ratio of the binary (q = M_secondary/M_primary)
     *
     * A_BIN:     Semimajor axis of the binary in code units
     * E_BIN:     Eccentricity of the binary
     * I_BIN:     Inclination, tilt of the binary with respect to the 
     *            x-y-plane, measured at the ascending node.
     *            For two-dimensional simulations the following values are
     *            possible:
     *            I_BIN = 0.0 for prograde orbits
     *            I_BIN = pi  for retrograde orbits
     * OMEGA_BIN: Longitude of the ascending node (where the orbit passes
     *            upward through the x-y-plane), measured from the
     *            x-axis.
     *            For two-dimensional simulations this value should be zero
     * PERI_BIN:  Argument of periapsis, angle between the ascending node
     *            and the pericentre of the binary
     * F_BIN:     True anomaly. Position of the secondary on its 
     *            orbit with respect to the pericentre at time t=0.
     *
     */
    double q = g_inputParam[Q_BIN];
    double a = g_inputParam[A_BIN];
    double e = g_inputParam[E_BIN];
    double inc = g_inputParam[I_BIN];
    double Omega = g_inputParam[OMEGA_BIN];
    double omega = g_inputParam[PERI_BIN];
    double f = g_inputParam[F_BIN];

    double k1 = q/(1.+q);
    double k2 = 1./(1.+q);

    double cO = cos(Omega);
    double sO = sin(Omega);
    double co = cos(omega);
    double so = sin(omega);
    double cf = cos(f);
    double sf = sin(f);
    double ci = cos(inc);
    double si = sin(inc);

    double r = a*(1.-e*e)/(1.+e*cf);
    double v0 = sqrt(CONST_G_CODE_UNITS*m/a/(1.-e*e));

    double x = r*(cO*(co*cf-so*sf) - sO*(so*cf+co*sf)*ci);
    double y = r*(sO*(co*cf-so*sf) + cO*(so*cf+co*sf)*ci);
    double z = r*(so*cf+co*sf)*si;
    
    double vx = v0*((e+cf)*(-ci*co*sO - cO*so) - sf*(co*cO - ci*so*sO));
    double vy = v0*((e+cf)*(ci*co*cO - sO*so) - sf*(co*sO + ci*so*cO));
    double vz = v0*((e+cf)*co*si - sf*si*so);

    /* primary */
    g_nb.m[0] = m/(1.+q);

    g_nb.x[0] = -k1*x;
    g_nb.y[0] = -k1*y;
    g_nb.z[0] = -k1*z;

    g_nb.vx[0] = -k1*vx;
    g_nb.vy[0] = -k1*vy;
    g_nb.vz[0] = -k1*vz;

    g_nb.feelsDisk[0] = CO_FEELS_DISK;

    /* secondary */
    g_nb.m[1] = m*q/(1.+q);

    g_nb.x[1] = k2*x;
    g_nb.y[1] = k2*y;
    g_nb.z[1] = k2*z;

    g_nb.vx[1] = k2*vx;
    g_nb.vy[1] = k2*vy;
    g_nb.vz[1] = k2*vz;

    g_nb.feelsDisk[1] = CO_FEELS_DISK;
    #endif
}

void nbodyInitialize()
{
    nbodyAllocateMemory();
    nbodyAllocateRK4Memory();

    if (prank == 0)
    {
        nbodyInitializeCentralObject();

        if (NO_OF_PLANETS > 0)
        {
            nbodyReadPlanetFile("planet.ini");
        }
    }
    #ifdef PARALLEL
    MPI_Bcast(g_nb.m, NB_N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(g_nb.x, NB_N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(g_nb.y, NB_N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(g_nb.z, NB_N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(g_nb.vx, NB_N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(g_nb.vy, NB_N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(g_nb.vz, NB_N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(g_nb.feelsDisk, NB_N, MPI_INT, 0, MPI_COMM_WORLD);
    #endif
}
