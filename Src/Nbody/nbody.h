#ifndef NBODY_H
#define NBODY_H

#define BUFFER_LENGTH 500
#define ORBITAL_EPS 1.e-8

typedef struct NBODY_SYSTEM {
    double *m;
    double *x, *y, *z;
    double *vx, *vy, *vz;
    double *ax, *ay, *az;
    double *axdisk, *aydisk, *azdisk;
    double ax_indirect, ay_indirect, az_indirect;
    double *xold, *yold, *zold;
    double *vxold, *vyold, *vzold;
    int *feelsDisk;
    
    /* orbital elements */
    double *a;           /* semi-major axis */
    double *e;           /* eccentricity */
    double *inc;         /* inclination */
    double *Omega;       /* longitude of ascending node (zero for inc = 0 or PI) */
    double *omega;       /* argument of pericentre */
    double *f;           /* true anomaly */
    double *P;           /* period */
    double *ea;          /* eccentric anomaly */
    double *ma;          /* mean anomaly */

    /* centre of mass */
    double comX, comVX;
    double comY, comVY;
    double comZ, comVZ;
} Nbody_System;

extern double *k1x; 
extern double *k1y; 
extern double *k1z;

extern double *k1vx;
extern double *k1vy;
extern double *k1vz;

extern double *k2x, *k2y, *k2z;
extern double *k2vx, *k2vy, *k2vz;

extern double *k3x, *k3y, *k3z;
extern double *k3vx, *k3vy, *k3vz;

extern double *k4x, *k4y, *k4z;
extern double *k4vx, *k4vy, *k4vz;

/* nbody.c prototypes */
double nbodySmoothingSquared(double *v, double x1, double x2, double x3);

/* nbody_gravity.c prototypes */
void nbodyCalcIndirectTerm();
void nbodyCalcAccelerations();
void nbodyCalcDiskFeedback(const Data *d, Grid *grid);
double nbodyCalcCellAcceleration(double *v, double x1, double x2, double x3, int dir);

/* nbody_init.c prototypes */
void nbodyInitialize();

/* nbody_io.c prototypes */
int nbodyReadPlanetFile(const char *filename);
void nbodyWriteCoordinates();
void nbodyWriteOrbitalElements();
void nbodyWriteCOM();
void nbodyWriteRestartFile(Output *output);
void nbodyWriteOrbitalElementsAtOutputTime(Output *output);
void nbodyReadRestartFile(Runtime *ini, int nrestart);
void nbodyPrintCoordinates();

/* nbody_memory.c prototypes */
void nbodyAllocateMemory();
void nbodyFreeMemory();
void nbodyAllocateRK4Memory();
void nbodyFreeRK4Memory();

/* nbody_rk4.c prototypes */
void nbodyAdvanceSystem(double dt);

/* nbody_tools.c prototypes */
void nbodySaveOldState();
double nbodyMass(int n);
void nbodyCenterOfMassCoordinatesAndVelocity(int n);
void nbodyCalcOrbitalElements();

#endif /* NBODY_H */
