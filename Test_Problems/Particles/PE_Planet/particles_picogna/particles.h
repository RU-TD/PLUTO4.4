/*  **************************************************************
            Header file for particle module.
 **************************************************************  */
#ifndef MLW_H_PARTICLES
#define MLW_H_PARTICLES

//#define DEBUG
//#define PROC_DEBUG
#define PI_NONE 0
#define PI_TEST 1
#define PI_ZHU2014 2


typedef struct PARTICLES_CONFIG {
    int num_particles; /* total number of particles to add to the global domain */

    int left_bound[3]; /* Array of left boundary types for particles */
    int right_bound[3];  /* Array of right boundary types for particles */

    int output_dbl_dn; /* double output frequency in number of steps */
    int output_tab_dn; /* ascii  output frequency in number of steps */
    int output_h5_dn; /* hdf5  binary output frequency in number of steps */
    int output_anl_dn; /* analysis      frequency in number of steps */

    double output_dbl_dt; /* double output frequency in time units      */
    double output_tab_dt; /* ascii  output frequency in time units      */
    double output_h5_dt; /* hdf5 binary output frequency in time units      */
    double output_anl_dt; /* analysis      frequency in time units      */

} ParticlesConfig;


typedef struct PARTICLE {
    int identity; /*   identifiant of  particle  */
    int cell [DIMENSIONS]; /*   particle'   cell index    */
    int rank;
    double coor[DIMENSIONS]; /*   coordinates of particle   */
    double mom [DIMENSIONS]; /*   momenta     of particle   */
    double radius;
    double coor_old[DIMENSIONS];
    double tstop;
    double g_grav[DIMENSIONS];
    double g_drag[DIMENSIONS];
    char   flag;
} Particle;


extern ParticlesConfig particlesConf;
extern int NPARTICLES;
extern int NOP;

extern int * particles;
extern int * particleCapacity;
extern Particle ** particlesPerProc;

extern int nprocs;


#ifdef PARALLEL
 extern MPI_Datatype ParticleType;
 extern MPI_Datatype PartConfigType;
 extern Particle * partsArray;
 extern int partsArrayLength;
#endif

//static int nprocs;
#define PROC_MAXSIZE 8192
static double *x0_proc, *x1_proc;
static double *y0_proc, *y1_proc;
static double *z0_proc, *z1_proc;

/*  ----  Prototyping goes here  ----  */
void ParticlesSet(Grid *grid, const Data *d, Runtime *runtime);
void ParticlesInit(Particle *pl, Grid *grid, const Data *d, int j);
void ParticlesPredictor(double ***uu[], Grid *grid);
void ParticlesCorrector(double ***uu[], Grid *grid);
int ParticlesBoundary(double ***uu[], Grid *grid, Particle * particles);
void ParticlesAnalysis(Particle * particles, int nop, int nstep, double glob_time);
void ParticlesSave(int nop, Runtime *);
void PartSaveDbl(Particle * particles, int nop, int nstep, double glob_time, char *output_dir);
void PartLoadDbl(Particle * particles, int nop, int nstep, char *data_dir);
void PartSaveTab(Particle * particles, int nop, int nstep, double glob_time, char *output_dir);
void PartSaveVTK(Particle * particles, int nop, int nstep, double glob_time, char *output_dir);
void WriteParticleData(Particle * particles, int);
void calculateGravitation(double ****v, Grid *grid, Particle * particles);
double calculateDrag(double ****v, Grid *grid, Particle *pl);
void CorrectPositions(double ***uu[], Grid *grid, Particle * particles);
void PredictPositions(double ***uu[], Grid *grid, Particle * particles);

/* Useful routines for initialization */
void ParticlesInterpol(double *a_interpol, Particle *pl, double ***uu[], Grid *grid, int VAR);
void ParticlesInitVelocity(Particle * pl, double ***uu[], Grid *grid);
void ParticlesChkProcessor(Particle * pl);
void ParticlesAllocateArrays();
void ParticlesReadConfig();
void PARTICLE_CONFIG_DATATYPE();
void PARTICLE_STRUCT_DATATYPE();
void CreateDomainDecomposition(Grid *grid);
int ParticlesChkInflow(Particle *pl, Grid *grid, double ***uu[], int nop);
void PartInflowBound(Particle *pl, double ***uu[], Grid *grid, int nop);
void ParticlesLocate(double *, int *, Grid *grid);
void AddTurbulence(double ***uu[], Grid *grid, Particle * pl);
void transformSphericalCartesian(double * cartesian,double * spherical,double r, double theta, double phi);
void transformCartesianSpherical(double * cartesian,double * spherical,double r, double theta, double phi);
void transformParticleSpherical(Particle *pl);
void transformParticleCartesian(Particle *pl);
double randn ();


#endif