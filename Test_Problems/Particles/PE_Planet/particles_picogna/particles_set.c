#include "pluto.h"
#include "particles.h"
#include "units.h"
#include "planet.h"


ParticlesConfig particlesConf;
int NPARTICLES = 0;
int NOP = 0;
int *particles = NULL;
int *particleCapacity = NULL;
int nprocs = 0;

Particle ** particlesPerProc;

#ifdef PARALLEL
    MPI_Datatype ParticleType;
    MPI_Datatype PartConfigType;
    Particle * partsArray;
    int partsArrayLength = 100;
#endif

// taken from runtime_setup.c
#define COMPARE(s1,s2,ii) \
        for ( (ii) = 1 ; (ii) < NOPT && !(strcmp ( (s1), (s2)) == 0); (ii)++);
#define NOPT 32

inline double wrapAngle( double angle )
{
    double twoPi = 2.0 * CONST_PI;
    return angle - twoPi * floor( angle / twoPi );
}

double randn1()
{
  double U1, U2, W, mult;
  static double X1, X2;
  static int call = 0;

  if (call == 1)
    {
      call = !call;
      return ((double) X2);
    }

  do
    {
      U1 = -1 + ((double) rand () / RAND_MAX) * 2;
      U2 = -1 + ((double) rand () / RAND_MAX) * 2;
      W = pow (U1, 2) + pow (U2, 2);
    }
  while (W >= 1 || W == 0);

  mult = sqrt ((-2 * log (W)) / W);
  X1 = U1 * mult;
  X2 = U2 * mult;

  call = !call;

  return ((double) X1);
}

void AddTurbulence(double ***uu[], Grid *grid, Particle * pl)
/****************************************************************
 *
 * Add kicks in particle positions due to turbulence motions
 * (see Charnoz et al. 2011)
 *
 **************************************************************** */
{
    double dt = g_dt;
    double density, pressure, Sc, Dd, mean, vari, c_s, H, R, tau, OmegaK;
    double hscale[DIMENSIONS], dx[DIMENSIONS], drho[DIMENSIONS];
    int i, j, k, l, ip, im, jp, jm, kp, km;
    static int first_call = 1;

    if (first_call == 1){
        srand(0);
        first_call = 0;
    }
    ParticlesLocate(pl->coor, pl->cell, grid);
    ParticlesInterpol(&density, pl, uu, grid, RHO);

#if GEOMETRY == SPHERICAL
    R = pl->coor[0]*sin(pl->coor[1]);
#elif GEOMETRY == POLAR
    R = pl->coor[0];
#else
    printLog("Geometry not implemented for the particle integrator/n");
    QUIT_PLUTO(1);
#endif

#if EOS == IDEAL
    ParticlesInterpol(&pressure, pl, uu, grid, PRS);
    c_s = sqrt(g_gamma*pressure/density);
#elif EOS == ISOTHERMAL
    c_s = g_isoSoundSpeed/sqrt(R);
#else
    printLog("EOS not implemented for the particle integrator/n")
    QUIT_PLUTO(1);
#endif

    i = pl->cell[0];
    j = pl->cell[1];
    k = pl->cell[2];
    ip = pl->cell[0]+1;
    im = pl->cell[0]-1;
    jp = pl->cell[1]+1;
    jm = pl->cell[1]-1;
    kp = pl->cell[2]+1;
    km = pl->cell[2]-1;

    hscale[0] = 1.0;
   #if (GEOMETRY == POLAR) || (GEOMETRY == SPHERICAL)
    hscale[1] = pl->coor[0];
   #else
    hscale[1] = 1.0;
   #endif
   #if GEOMETRY == SPHERICAL
    hscale[2] = R;
   #else
    hscale[2] = 1.0;
   #endif

    OmegaK = 2 * CONST_PI / (sqrt(R)*R);
    H = c_s/OmegaK;
    tau = pl->tstop*OmegaK;
    Sc = POW2(1.0 + tau*tau)/(1.0+4.0*tau*tau);
    Dd = g_inputParam[ALPHA]*c_s*H/Sc;

    dx[0] = 2.*grid->dx[IDIR][i];
    dx[1] = 2.*grid->dx[JDIR][j];
    dx[2] = 2.*grid->dx[KDIR][k];

    drho[0] = (uu[RHO][k][j][ip]-uu[RHO][k][j][im])/hscale[0]/dx[0];
    drho[1] = (uu[RHO][k][jp][i]-uu[RHO][k][jm][i])/hscale[1]/dx[1];
    drho[2] = (uu[RHO][kp][j][i]-uu[RHO][km][j][i])/hscale[2]/dx[2];

    vari = sqrt(2.0*Dd*dt);

    for (l=0; l<DIMENSIONS; l++) {
        mean = Dd/density*drho[l]*dt;
        if (l != iVTH) {
            pl->coor[l] += mean+randn1()*vari;
            pl->coor_old[l] += mean+randn1()*vari;
            if (l == iVPHI) {
                wrapAngle( pl->coor[l] );
                wrapAngle( pl->coor_old[l] );
            }
        }
    }
}

/* ********************************************************************* */
void CreateDomainDecomposition(Grid *grid)
    /*!
     * Show the parallel domain decomposition by having each processor print
     * its own computational domain.
     * This is activated with the -show-dec command line argument.
     * It may be long for several thousand processors.
     *
     * \param [in] grid     a pointer to the grid structure
     *********************************************************************** */
{
    int i, j, k, ngh;
    int i0, i1, j0, j1, k0, k1;
    int nxp, nyp, nzp;
    int nghx, nghy, nghz;
    double x0, x1, y0, y1, z0, z1;

#ifdef PARALLEL
    extern int prank;
    MPI_Status status;
#endif

    static int *i0_proc, *i1_proc;
    static int *j0_proc, *j1_proc;
    static int *k0_proc, *k1_proc;

    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    /* ---- Allocate memory ---- */

    i0_proc = ARRAY_1D(nprocs, int); i1_proc = ARRAY_1D(nprocs, int); // these will hold the global start and end
    j0_proc = ARRAY_1D(nprocs, int); j1_proc = ARRAY_1D(nprocs, int); // indices of each processor's local domains
    k0_proc = ARRAY_1D(nprocs, int); k1_proc = ARRAY_1D(nprocs, int);
    x0_proc = ARRAY_1D(nprocs, double); x1_proc = ARRAY_1D(nprocs, double); // these will hold the global start and end
    y0_proc = ARRAY_1D(nprocs, double); y1_proc = ARRAY_1D(nprocs, double); // indices of each processor's local domains
    z0_proc = ARRAY_1D(nprocs, double); z1_proc = ARRAY_1D(nprocs, double);

#ifdef PARALLEL
    nxp = grid->np_tot[0]; // number of points in the local domain (boundaried included)
    nyp = grid->np_tot[1];
    nzp = grid->np_tot[2];

    i0 = nghx = grid->nghost[0]; // local start index of the local array without boundary
    j0 = nghy = grid->nghost[1];
    k0 = nghz = grid->nghost[2];

    i1 = i0 + grid->np_int[0] - 1; // local end index of the local array without boundary
    j1 = j0 + grid->np_int[1] - 1;
    k1 = k0 + grid->np_int[2] - 1;

    x0 = grid->xl[0][i0]; x1 = grid->xr[0][i1]; // local start and end point without boundary
    y0 = grid->xl[1][j0]; y1 = grid->xr[1][j1];
    z0 = grid->xl[2][k0]; z1 = grid->xr[2][k1];

    DIM_EXPAND(
        MPI_Gather (&x0, 1, MPI_DOUBLE, x0_proc, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Gather (&x1, 1, MPI_DOUBLE, x1_proc, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);  ,

        MPI_Gather (&y0, 1, MPI_DOUBLE, y0_proc, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Gather (&y1, 1, MPI_DOUBLE, y1_proc, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);  ,

        MPI_Gather (&z0, 1, MPI_DOUBLE, z0_proc, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Gather (&z1, 1, MPI_DOUBLE, z1_proc, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        )

    printLog ("> Domain Decomposition (%d procs):\n\n", nprocs);

    if (nprocs  > 1){
        DIM_EXPAND(
            MPI_Bcast(x0_proc,nprocs,MPI_DOUBLE,0,MPI_COMM_WORLD);
            MPI_Bcast(x1_proc,nprocs,MPI_DOUBLE,0,MPI_COMM_WORLD); ,

            MPI_Bcast(y0_proc,nprocs,MPI_DOUBLE,0,MPI_COMM_WORLD);
            MPI_Bcast(y1_proc,nprocs,MPI_DOUBLE,0,MPI_COMM_WORLD); ,

            MPI_Bcast(z0_proc,nprocs,MPI_DOUBLE,0,MPI_COMM_WORLD);
            MPI_Bcast(z1_proc,nprocs,MPI_DOUBLE,0,MPI_COMM_WORLD);
        )

    }
    if (prank == 0) {
        for (k = 0; k < nprocs; k++){
            printLog ("  - Proc # %d, X1: [%f, %f], %d pt\n",k, x0_proc[k], x1_proc[k], grid->np_int[0]);
            #if DIMENSIONS > 1
                printLog ("              X2: [%f, %f], %d pt\n",    y0_proc[k], y1_proc[k], grid->np_int[1]);
            #endif
            #if DIMENSIONS > 2
                printLog ("              X3: [%f, %f], %d pt\n\n",  z0_proc[k], z1_proc[k], grid->np_int[2]);
            #endif
        }
    }
    printLog ("\n");
#endif
}

void ParticlesAllocateArrays() {
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    PARTICLE_STRUCT_DATATYPE();

    particles = ARRAY_1D(nprocs, int);
    particlesPerProc = ARRAY_1D(nprocs, Particle *);
    particleCapacity = ARRAY_1D(nprocs, int);
#ifdef PARALLEL
    partsArray = ARRAY_1D(partsArrayLength, Particle);
#endif
}

void ParticlesSet(Grid *grid, const Data *d, Runtime *runtime) {
    Particle *pl;
    int i = 0, j = 0, t = 0;
#ifdef PARALLEL
    int partPerProc = 0;
    int resto = 0;
    int tragetProc = 0;
    int index = 0;
    int h = 0;
    int partsToBeReceived = 0;
    int partsToBeSent = 0;
    int * counters = ARRAY_1D(nprocs, int);
    int * partAdded = ARRAY_1D(nprocs, int);
    int ttl = 0;
    Particle * allPl = NULL;
    Particle ** appoggio = ARRAY_1D(nprocs, Particle *);
    MPI_Status recv;
    MPI_Status Bstatus;
    PARTICLE_CONFIG_DATATYPE();
    Particle * appoggio1D = NULL;
    Particle * appoggio1D_2 = NULL;
    int diff = 0;

    //hardcoded restarting:
    //TODO better restarting interface
    int restart = 0;
    int restart_file = 0;
    if (restart)
        NOP = 1000000;
#endif

    printLog("> Particles initialization \n");
    ParticlesReadConfig();
#ifdef PARALLEL
    partPerProc = NOP / nprocs;
    resto = NOP % nprocs;

    for (i = 0; i < nprocs; i++) {
        particles[i] = 0; // number of particles for each process i
        partAdded[i] = 0; // 
        particleCapacity[i] = MIN(partPerProc+resto+50, PROC_MAXSIZE); // number of particles each process can hold
        particlesPerProc[i] = ARRAY_1D(particleCapacity[i], Particle); // the actual array of particles for process i 
    }
    i = 0;

    if (prank == 0) {
        allPl = ARRAY_1D(NOP, Particle); // array that will hold all particles before distributing them among processes
        if (restart){
            PartLoadDbl(allPl, NOP, restart_file, runtime->output_dir);
        }
        for (i = 0; i < NOP; i++) {
            if (restart){
                ParticlesChkProcessor(&allPl[i]);
            }
            else
            {
                ParticlesInit(&allPl[i], grid, d, i); // creates all particles (assigns them to process 0 (p->rank=0) )
                ParticlesChkProcessor(&allPl[i]); // reassign the particles to the correct processes
            }
            pl = &(allPl[i]);

            if (particles[pl->rank] < particleCapacity[pl->rank]) {
                particlesPerProc[pl->rank][particles[pl->rank]] = *pl;
                if (pl->rank == 0) {
                    if (!(restart))
                        ParticlesInitVelocity(&particlesPerProc[pl->rank][particles[pl->rank]], d->Vc, grid);
                    partAdded[pl->rank]++;
                }
                particles[pl->rank]++;

            } else {
                particleCapacity[pl->rank] += particleCapacity[pl->rank] / 2;
                appoggio1D = ARRAY_1D(particleCapacity[pl->rank], Particle);
                for (h = 0; h < particles[pl->rank]; h++) {
                    appoggio1D[h] = particlesPerProc[pl->rank][h];
                }
                FreeArray1D((void *) particlesPerProc[pl->rank]);
                particlesPerProc[pl->rank] = appoggio1D;
                particlesPerProc[pl->rank][particles[pl->rank]] = *pl;
                if (pl->rank == 0) {
                    if (!(restart))
                        ParticlesInitVelocity(&particlesPerProc[pl->rank][particles[pl->rank]], d->Vc, grid);
                    partAdded[pl->rank]++;
                }
                particles[pl->rank]++;
            }

        } // ENDS THE LOOP ON NOP. PUTS ALL PARTICLES INTO RANK = 0 proc. and reamaining are stored in
        //PARTICLESPERPROC array corresponding to each processor.
        for (i = 1; i < nprocs; i++) {
            MPI_Send(&particles[i], 1, MPI_INT, i, 123, MPI_COMM_WORLD);
            if (particles[i] > 0) {
                MPI_Send(particlesPerProc[i], particles[i], ParticleType, i, 123, MPI_COMM_WORLD);
                {
                    partAdded[i] += particles[i];
                    particles[i] = 0;
                }
            }
        } // END LOOP FOR EACH PROCESSOR. Here the particles are send to the respective processor untill they
        // fill the counter limit. Remaining ones are send back. Update the Paricle Added counter and particle count.

    } else {
        MPI_Recv(&particles[prank], 1, MPI_INT, 0, 123, MPI_COMM_WORLD, &recv);
        if (particles[prank] > 0){
            if (particles[prank] < particleCapacity[prank]) {
                MPI_Recv(particlesPerProc[prank], particles[prank], ParticleType, 0, 123, MPI_COMM_WORLD, &recv);
            } else{
                FreeArray1D((void *) particlesPerProc[prank]);
                particleCapacity[prank] = particles[prank];
                particlesPerProc[prank] = ARRAY_1D(particleCapacity[prank], Particle);
                MPI_Recv(particlesPerProc[prank], particles[prank], ParticleType, 0, 123, MPI_COMM_WORLD, &recv);
            }
        }

        for (i = 0; i < particles[prank]; i++){
            if (!(restart))
                ParticlesInitVelocity(&particlesPerProc[prank][i], d->Vc, grid);
        }
    }      

    FreeArray1D((void *) partAdded);
#else
    particlesPerProc[0] = ARRAY_1D(NOP, Particle);
    for (i = 0; i < NOP; i++) {
        ParticlesInit(&particlesPerProc[0][i], grid, d, i);
        ParticlesInitVelocity(&particlesPerProc[0][i], d->Vc, grid);
    }
#endif
}


/* ************************************************************************* */
void ParticlesPredictor(double ***uu[], Grid *grid)
{
    /*
 *
 * Performs a predictor step for the particle pushing algorithm.
 *
 * At its core, this routine is a simple loop over the particles.
 *
 * Note that before anything is compute, new particles are injected (if needed)
 * with the CHECK_PART_INFLOW function.
 *
     *************************************************************************** */

    extern int prank;
    int i, j, k;

#ifdef PARALLEL
    int prevRank = -1;
    int h = 0;
    Particle * appoggio = NULL;
    int check = 1;
#endif

    /* Variables to handle data structure*/
    Particle *pl;
    for (i = 0; i < nprocs; i++) {
        for (j = 0; j < particles[i]; j++) {
            pl = &(particlesPerProc[i][j]);
            PredictPositions(uu, grid, pl);
        }
    }
#ifdef PARALLEL
    for (i = 0; i < nprocs; i++) {
        for (j = 0; j < particles[i]; j++) {
            pl = &(particlesPerProc[i][j]);
            prevRank = pl->rank;
            check = ParticlesBoundary(uu, grid, pl);
            if (!(check))
            {
                for (h = j; h < particles[i] - 1; h++) {
                    particlesPerProc[i][h] = particlesPerProc[i][h + 1];
                }
                particles[i]--;
                continue;
            }

            ParticlesChkProcessor(pl);
            if (prevRank != pl->rank) {
                if (particles[pl->rank] < particleCapacity[pl->rank]) {
                    // The new processor has sufficient space to receive the particle
                    particlesPerProc[pl->rank][particles[pl->rank]] = particlesPerProc[i][j];
                    //TODO remove this line? it should do nothing
                    particlesPerProc[pl->rank][particles[pl->rank]].rank = pl->rank;
                    particles[pl->rank]++;
                } else {
                    particleCapacity[pl->rank] += particleCapacity[pl->rank] / 2;
                    appoggio = ARRAY_1D(particleCapacity[pl->rank], Particle);
                    for (h = 0; h < particles[pl->rank]; h++) {
                        appoggio[h] = particlesPerProc[pl->rank][h];
                    }
                    FreeArray1D((void *) particlesPerProc[pl->rank]);
                    particlesPerProc[pl->rank] = appoggio;
                    particlesPerProc[pl->rank][particles[pl->rank]++] = particlesPerProc[i][j];
                }
                // Remove the particle from the previous processor
                for (h = j; h < particles[prevRank] - 1; h++) {
                    particlesPerProc[prevRank][h] = particlesPerProc[prevRank][h + 1];
                }
                particles[prevRank]--;
                j--;
            }
        }
    }
#endif
    for (i = 0; i < nprocs; i++) {
        for (j = 0; j < particles[i]; j++) {
            pl = &(particlesPerProc[i][j]);
            if (pl->flag == 'I') ParticlesInitVelocity(pl, uu, grid);
        }
    }
}


/* ************************************************************** */
void ParticlesCorrector(double ***uu[], Grid *grid)
    /*
     * Corrector part of RK2 Intergration for particle velocity.
     *
     **************************************************************** */
{
    int i = 0;
    int j = 0;
    int h = 0;
    int sum;
    int check;
    double *v_interpol[COMPONENTS];
    extern int prank;
    Particle *pl;
    int prevRank = -1;

    Particle * appoggio = NULL;
#ifdef PARALLEL
    MPI_Request request;
    MPI_Request requestRecv;
    MPI_Status status;
    MPI_Status statusRecv;
    int partsToBeReceived = 0;
    Particle * appoggio1D = NULL;

    for (i = 0; i < nprocs; i++) {
        //print("proc %d size %d prank %d\n",i,particles[i],prank);
        if (!(i == prank)) {
            partsToBeReceived = 0;
            MPI_Irecv(&partsToBeReceived,1,MPI_INT, i, 123, MPI_COMM_WORLD, &requestRecv);

            /*MPI_Isend(&particles[i], 1, MPI_INT, i, 123, MPI_COMM_WORLD, &request);*/
            MPI_Send(&particles[i], 1, MPI_INT, i, 123, MPI_COMM_WORLD);

            /*MPI_Recv(&partsToBeReceived, 1, MPI_INT, i, 123, MPI_COMM_WORLD, &status);*/
            MPI_Wait(&requestRecv,&status);

            if(partsToBeReceived > 0){
                if (partsToBeReceived > partsArrayLength) {
                    FreeArray1D((void *) partsArray);
                    partsArray = ARRAY_1D(partsToBeReceived, Particle);
                    partsArrayLength = partsToBeReceived;
                }
                MPI_Irecv(partsArray, partsToBeReceived, ParticleType, i , 124, MPI_COMM_WORLD, &requestRecv);
            }
            if (particles[i] > 0) {
                /*MPI_Isend(particlesPerProc[i], particles[i], ParticleType, i, 124, MPI_COMM_WORLD, &request);*/
                MPI_Send(particlesPerProc[i], particles[i], ParticleType, i, 124, MPI_COMM_WORLD);
                particles[i] = 0;
            }
            if (partsToBeReceived > 0) {
                MPI_Wait(&requestRecv,&status);
                /* MPI_Recv(partsArray, partsToBeReceived, ParticleType, i, 124, MPI_COMM_WORLD, &status);*/
                for (j = 0; j < partsToBeReceived; j++) {

                    if (particles[prank] < particleCapacity[prank]) {
                        particlesPerProc[prank][particles[prank]] = partsArray[j];
                        particles[prank]++;
                    } else {
                        particleCapacity[prank] += particleCapacity[prank] / 2;
                        appoggio1D = ARRAY_1D(particleCapacity[prank], Particle);
                        for (h = 0; h < particles[prank]; h++) {
                            appoggio1D[h] = particlesPerProc[prank][h];
                        }
                        FreeArray1D((void *) particlesPerProc[prank]);
                        particlesPerProc[prank] = appoggio1D;
                        particlesPerProc[prank][particles[prank]] = partsArray[j];
                        particles[prank]++;
                    }
                }
            }
            if (particles[i] > 0) {
                /*MPI_Irecv(particlesPerProc[i], particles[i], ParticleType, i, 125, MPI_COMM_WORLD, &request);*/
                MPI_Wait(&request, &status);
            }
        }
    }
    i = prank;
    for (j = 0; j < particles[i]; j++) {
        pl = &particlesPerProc[i][j];
        //        if (pl->coor[0]<0.04)
        //            print("ID: %d prank %d\n",pl->identity,prank);
        // check = ParticlesBoundary(uu, grid, pl);
        check = 1;
        if (check == 0) {
            //TODO changed h = 0 to h = j, is this correct?
            for (h = j; h < particles[i] - 1; h++) {
                particlesPerProc[i][h] = particlesPerProc[i][h + 1];
            }
            particles[i]--;
            j--;
            //TODO does prank 0 need to know that?
            NOP--;
        } else {
#if PARTICLES_TURBULENCE
            AddTurbulence(uu, grid, pl);
#endif
            CorrectPositions(uu, grid, pl);
        }
    }
    for (i = 0; i < nprocs; i++) {
        for (j = 0; j < particles[i]; j++) {
            pl = &(particlesPerProc[i][j]);
            prevRank = pl->rank;
            check = ParticlesBoundary(uu, grid, pl);
            if (!(check))
            {
                for (h = j; h < particles[i] - 1; h++) {
                    particlesPerProc[i][h] = particlesPerProc[i][h + 1];
                }
                particles[i]--;
                continue;
            }

            ParticlesChkProcessor(pl);
            if (prevRank != pl->rank) {
                if (particles[pl->rank] < particleCapacity[pl->rank]) {
                    // The new processor has sufficient space to receive the particle
                    particlesPerProc[pl->rank][particles[pl->rank]] = particlesPerProc[i][j];
                    //TODO remove this line? it should do nothing
                    particlesPerProc[pl->rank][particles[pl->rank]].rank = pl->rank;
                    particles[pl->rank]++;
                } else {
                    particleCapacity[pl->rank] += particleCapacity[pl->rank] / 2;
                    appoggio = ARRAY_1D(particleCapacity[pl->rank], Particle);
                    for (h = 0; h < particles[pl->rank]; h++) {
                        appoggio[h] = particlesPerProc[pl->rank][h];
                    }
                    FreeArray1D((void *) particlesPerProc[pl->rank]);
                    particlesPerProc[pl->rank] = appoggio;
                    particlesPerProc[pl->rank][particles[pl->rank]++] = particlesPerProc[i][j];
                }
                // Remove the particle from the previous processor
                for (h = j; h < particles[prevRank] - 1; h++) {
                    particlesPerProc[prevRank][h] = particlesPerProc[prevRank][h + 1];
                }
                particles[prevRank]--;
                j--;
            }
        }
    }
    for (i = 0; i < nprocs; i++) {
        //print("proc %d size %d prank %d\n",i,particles[i],prank);
        if (!(i == prank)) {
            partsToBeReceived = 0;
            MPI_Irecv(&partsToBeReceived,1,MPI_INT, i, 123, MPI_COMM_WORLD, &requestRecv);

            /*MPI_Isend(&particles[i], 1, MPI_INT, i, 123, MPI_COMM_WORLD, &request);*/
            MPI_Send(&particles[i], 1, MPI_INT, i, 123, MPI_COMM_WORLD);

            /*MPI_Recv(&partsToBeReceived, 1, MPI_INT, i, 123, MPI_COMM_WORLD, &status);*/
            MPI_Wait(&requestRecv,&status);

            if(partsToBeReceived > 0){
                if (partsToBeReceived > partsArrayLength) {
                    FreeArray1D((void *) partsArray);
                    partsArray = ARRAY_1D(partsToBeReceived, Particle);
                    partsArrayLength = partsToBeReceived;
                }
                MPI_Irecv(partsArray, partsToBeReceived, ParticleType, i , 124, MPI_COMM_WORLD, &requestRecv);
            }
            if (particles[i] > 0) {
                /*MPI_Isend(particlesPerProc[i], particles[i], ParticleType, i, 124, MPI_COMM_WORLD, &request);*/
                MPI_Send(particlesPerProc[i], particles[i], ParticleType, i, 124, MPI_COMM_WORLD);

                particles[i] = 0;
            }
            if (partsToBeReceived > 0) {
                MPI_Wait(&requestRecv,&status);
                /* MPI_Recv(partsArray, partsToBeReceived, ParticleType, i, 124, MPI_COMM_WORLD, &status);*/
                for (j = 0; j < partsToBeReceived; j++) {

                    if (particles[prank] < particleCapacity[prank]) {
                        particlesPerProc[prank][particles[prank]] = partsArray[j];
                        particles[prank]++;
                    } else {
                        particleCapacity[prank] += particleCapacity[prank] / 2;
                        appoggio1D = ARRAY_1D(particleCapacity[prank], Particle);
                        for (h = 0; h < particles[prank]; h++) {
                            appoggio1D[h] = particlesPerProc[prank][h];
                        }
                        FreeArray1D((void *) particlesPerProc[prank]);
                        particlesPerProc[prank] = appoggio1D;
                        particlesPerProc[prank][particles[prank]] = partsArray[j];
                        particles[prank]++;
                    }
                }
            }
            if (particles[i] > 0) {
                /*MPI_Irecv(particlesPerProc[i], particles[i], ParticleType, i, 125, MPI_COMM_WORLD, &request);*/
                MPI_Wait(&request, &status);
            }
        }
    }
#else
    for (i = 0; i < nprocs; i++) {
        if (i == prank) {
            for (j = 0; j < particles[i]; j++) {
                pl = &particlesPerProc[i][j];
                check = ParticlesBoundary(uu, grid, pl);
                if (check == 0) {
                    for (h = 0; h < particles[i] - 1; h++) {
                        particlesPerProc[i][h] = particlesPerProc[i][h + 1];
                        particles[i]--;
                    }
                } else {
                    ParticlesLocate(pl->coor, pl->cell, grid);
                    ParticlesInterpol(v_interpol, pl, uu, grid, VX1);
                    for (h = 0; h < DIMENSIONS; ++h) {
                        pl->coor[h] = pl->coor_old[h] + g_dt * 0.5 * (pl->mom[h] + v_interpol[h]);
                    }
                    ParticlesInitVelocity(pl, uu, grid);
                }
            }
        }
    }
#endif

    for (i = 0; i < nprocs; i++) {
        for (j = 0; j < particles[i]; j++) {
            pl = &(particlesPerProc[i][j]);
            if (pl->flag == 'I') ParticlesInitVelocity(pl, uu, grid);
        }
    }

#ifdef PARALLEL
    MPI_Barrier(MPI_COMM_WORLD);
#endif

}

/* ************************************************************** */
void CorrectPositions(double ***uu[], Grid *grid, Particle * pl)
    /*
     * Corrector positions of particles.
     *
     **************************************************************** */
{
    double R;

    int k;
    double hdt = 0.5*g_dt;
    double dtk, R0, cent, r, th;

    double g_gravold[DIMENSIONS];
    double g_dragold[DIMENSIONS];
    double tstopold;
    double dtp, Rn;
    double mom_old[DIMENSIONS];

    tstopold = pl->tstop;
    if(tstopold <= g_dt) {
        for (k = 0; k < DIMENSIONS; k++) {
            g_gravold[k] = pl->g_grav[k];
            g_dragold[k] = pl->g_drag[k];
        }
    }

    calculateGravitation(uu, grid, pl);
    pl->tstop = calculateDrag(uu, grid, pl);
    r  = pl->coor[0];
    for (k = 0; k < DIMENSIONS; k++) {
        mom_old[k] = pl->mom[k];
    }

    if(tstopold > g_dt) {
        if (pl->tstop == 0.0) {
            dtk = g_dt;
        } else {
            dtk = g_dt/(1.0+hdt/pl->tstop);
        }
    } else {
        if ((pl->tstop == 0.0) || (tstopold == 0.0)) {
            dtk = hdt;
            dtp = 1.0;
        } else {
            dtk = hdt/(1.0 + hdt/tstopold + hdt/pl->tstop + hdt*g_dt/tstopold/pl->tstop);
            dtp = (1.0+g_dt/tstopold);
        }
    }

#if GEOMETRY == SPHERICAL

    th = pl->coor[1];
    R0 = r*sin(th);
    if(tstopold > g_dt) {
        pl->mom[2] += (pl->g_grav[2] + pl->g_drag[2]) * dtk;
        cent          = 0.5*cos(th)/sin(th)*(mom_old[2]*mom_old[2] + pl->mom[2]*pl->mom[2])/R0/R0;
        pl->mom[1] += (pl->g_grav[1] + pl->g_drag[1] + cent) * dtk;
        cent          = 0.5/r*((mom_old[1]*mom_old[1] + pl->mom[1]*pl->mom[1])/r/r +
                        (mom_old[2]*mom_old[2] + pl->mom[2]*pl->mom[2])/R0/R0);
        pl->mom[0] += (pl->g_grav[0] + pl->g_drag[0] + cent) * dtk;
    } else {
        Rn            = pl->coor_old[0]*sin(pl->coor_old[1]);
        pl->mom[2] += (pl->g_grav[2] + pl->g_drag[2]) * dtk * dtp;
        pl->mom[2] += (g_gravold[2] + g_dragold[2]) * dtk;
        cent          = cos(th)/sin(th)*pl->mom[2]*pl->mom[2]/R0/R0;
        pl->mom[1] += ((pl->g_grav[1] + pl->g_drag[1]) + cent) * dtk * dtp;
        cent          = cos(pl->coor_old[1])/sin(pl->coor_old[1])*mom_old[2]*mom_old[2]/Rn/Rn;
        pl->mom[1] += (g_gravold[1] + g_dragold[1] + cent) * dtk;
        cent          = 1.0/r*(pl->mom[1]*pl->mom[1]/r/r + pl->mom[2]*pl->mom[2]/R0/R0);
        pl->mom[0] += (pl->g_grav[0] + pl->g_drag[0] + cent) * dtk * dtp;
        cent          = 1.0/pl->coor_old[0]*(mom_old[1]*mom_old[1]/pl->coor_old[0]/pl->coor_old[0] +
                        mom_old[2]*mom_old[2]/Rn/Rn);
        pl->mom[0] += (g_gravold[0] + g_dragold[0] + cent) * dtk;
    }

#elif GEOMETRY == POLAR

    if(tstopold > g_dt) {
        pl->mom[2] += (pl->g_grav[2] + pl->g_drag[2]) * dtk;
        pl->mom[1] += (pl->g_grav[1] + pl->g_drag[1]) * dtk;
        cent          = 0.5/r/r/r*(mom_old[1]*mom_old[1] + pl->mom[1]*pl->mom[1]);
        pl->mom[0] += (pl->g_grav[0] + pl->g_drag[0] + cent) * dtk;
    } else {
        pl->mom[2] += (pl->g_grav[2] + pl->g_drag[2]) * dtk * dtp;
        pl->mom[2] += (g_gravold[2] + g_dragold[2]) * dtk;
        pl->mom[1] += (pl->g_grav[1] + pl->g_drag[1]) * dtk * dtp;
        pl->mom[1] += (g_gravold[1] + g_dragold[1]) * dtk;
        cent          = 0.5*pl->mom[1]*pl->mom[1]/r/r/r;
        pl->mom[0] += (pl->g_grav[0] + pl->g_drag[0] + cent) * dtk * dtp;
        cent          = 0.5*mom_old[1]*mom_old[1]/pl->coor_old[0]/pl->coor_old[0]/pl->coor_old[0];
        pl->mom[0] += (pl->g_grav[0] + pl->g_drag[0] + cent) * dtk;
    }

#else
    printLog("Only spherical and polar geometries implemented for particles.\n");
    QUIT_PLUTO(1);
#endif
#if GEOMETRY == SPHERICAL

    if(tstopold > g_dt) {
        R0           = pl->coor[0]*sin(pl->coor[1]);
        pl->coor[0] += pl->mom[0] * hdt;
        pl->coor[1] += 0.5*(pl->mom[1]/pl->coor_old[0]/pl->coor_old[0] +
                       pl->mom[1]/pl->coor[0]/pl->coor[0]) * hdt;
        R            = pl->coor[0]*sin(pl->coor[1]);
        pl->coor[2] += 0.5*(pl->mom[2]/R0/R0 + pl->mom[2]/R/R) * hdt;
#if ROTATING_FRAME == YES
        pl->coor[2] -= g_OmegaZ * hdt;
#endif
    }
    else
    {
        R0          = pl->coor_old[0]*sin(pl->coor_old[1]);
        pl->coor[0] = pl->coor_old[0] + (pl->mom[0] + mom_old[0]) * hdt;
        pl->coor[1] = pl->coor_old[1] + (pl->mom[1]/pl->coor[0]/pl->coor[0] +
                mom_old[1]/pl->coor_old[0]/pl->coor_old[0]) * hdt;
        R           = pl->coor[0]*sin(pl->coor[1]);
        pl->coor[2] = pl->coor_old[2] + (pl->mom[2]/R/R + mom_old[2]/R0/R0) * hdt;

#if ROTATING_FRAME == YES
        pl->coor[2] -= g_OmegaZ * g_dt;
#endif
    }

#elif GEOMETRY == POLAR

    if(tstopold > g_dt) {
        pl->coor[0] += pl->mom[0] * hdt;
        pl->coor[1] += 0.5*(pl->mom[1]/r/r +
                pl->mom[1]/pl->coor[0]/pl->coor[0]) * hdt;
#if ROTATING_FRAME == YES
     pl->coor[1] -= g_OmegaZ * hdt;
#endif
        pl->coor[2] += pl->mom[2] * hdt;
    }
    else
    {
        pl->coor[0] = pl->coor_old[0] + (pl->mom[0] + mom_old[0]) * hdt;
        pl->coor[1] = pl->coor_old[1] + (pl->mom[1]/pl->coor[0]/pl->coor[0] +
                mom_old[1]/pl->coor_old[0]/pl->coor_old[0]) * hdt;
#if ROTATING_FRAME == YES
        pl->coor[1] -= g_OmegaZ * g_dt;
#endif
        pl->coor[2] = pl->coor_old[2] + (pl->mom[2] + mom_old[2]) * hdt;
    }

#else
    printLog("Only spherical and polar geometries implemented for particles.\n");
    QUIT_PLUTO(1);
#endif
}

/* ************************************************************** */
void PredictPositions(double ***uu[], Grid *grid, Particle * pl)
/*
 * Predict positions of particles.
 *
 **************************************************************** */
{
    double dtp, R, R0;
    int k;

    if (pl->tstop > g_dt) {
        dtp = 0.5*g_dt;
    } else {
        dtp = g_dt;
    }

    for (k = 0; k < DIMENSIONS; k++) {
        pl->coor_old[k] = pl->coor[k];
    }
#if GEOMETRY == SPHERICAL
    R0 = pl->coor_old[0]*sin(pl->coor_old[1]);
    pl->coor[0] += pl->mom[0] * dtp;
    pl->coor[1] += 0.5*(pl->mom[1]/POW2(pl->coor_old[0]) + pl->mom[1]/POW2(pl->coor[0])) * dtp;
    R = pl->coor[0]*sin(pl->coor[1]);
    pl->coor[2] += 0.5*(pl->mom[2]/POW2(R0) + pl->mom[2]/POW2(R)) * dtp;
  #if ROTATING_FRAME == YES
    pl->coor[2] -= g_OmegaZ * dtp;
  #endif
#elif ((GEOMETRY == POLAR) || (GEOMETRY == CYLINDRICAL))
    pl->coor[0] += pl->mom[0] * dtp;
    pl->coor[1] += 0.5*(pl->mom[1]/POW2(pl->coor_old[0]) + pl->mom[1]/POW2(pl->coor[0])) * dtp;
  #if ROTATING_FRAME == YES
    pl->coor[1] -= g_OmegaZ * dtp;
  #endif
  #if (GEOMETRY == POLAR)
    pl->coor[2] += pl->mom[2] * dtp;
  #endif
#else
    printLog("Only spherical and polar/cylindrical geometries implemented for particles.\n");
    QUIT_PLUTO(1);
#endif
}

void ParticlesChkProcessor(Particle * pl) {
    int procIndex = 0;
    int x1Bool = 0;
    int x2Bool = 0;
    int x3Bool = 0;
    int i = 0;

    procIndex = pl->rank;
    for (i = 0; i < DIMENSIONS; i++) {
        if (i == 0) {
            if (pl->coor[i] >= x0_proc[procIndex] && pl->coor[i] <= x1_proc[procIndex]) {
                x1Bool = 1;
            }
        } else if (i == 1) {
            if (pl->coor[i] >= y0_proc[procIndex] && pl->coor[i] <= y1_proc[procIndex]) {
                x2Bool = 1;
            }
        } else if (i == 2) {
            if (pl->coor[i] >= z0_proc[procIndex] && pl->coor[i] <= z1_proc[procIndex]) {
                x3Bool = 1;
            }
        }
        if (DIMENSIONS == 3 && x1Bool && x2Bool && x3Bool) {
            pl->rank = procIndex;
            return;
        } else if (DIMENSIONS == 2 && x1Bool && x2Bool) {
            pl->rank = procIndex;
            return;
        } else if (DIMENSIONS == 1 && x1Bool) {
            pl->rank = procIndex;
            return;
        }
    }

    for (procIndex = 0; procIndex < nprocs; procIndex++) {
        x1Bool = 0;
        x2Bool = 0;
        x3Bool = 0;
        for (i = 0; i < DIMENSIONS; i++) {

            if (i == 0) {
                if (pl->coor[i] >= x0_proc[procIndex] && pl->coor[i] <= x1_proc[procIndex]) {
                    x1Bool = 1;
                }
            } else if (i == 1) {
                if (pl->coor[i] >= y0_proc[procIndex] && pl->coor[i] <= y1_proc[procIndex]) {
                    x2Bool = 1;
                }
            } else if (i == 2) {
                if (pl->coor[i] >= z0_proc[procIndex] && pl->coor[i] <= z1_proc[procIndex]) {
                    x3Bool = 1;
                }
            }
        }
        if (DIMENSIONS == 3 && x1Bool && x2Bool && x3Bool) {
            pl->rank = procIndex;
            break;
        } else if (DIMENSIONS == 2 && x1Bool && x2Bool) {
            pl->rank = procIndex;
            break;
        } else if (DIMENSIONS == 1 && x1Bool) {
            pl->rank = procIndex;
            break;
        }
    }
}

/* ****************************************** */
void ParticlesReadConfig()
    /*
     * Reads the file particles.ini and set the
     * frequency of output and analysis
     *
     ******************************************** */ 
{
    int nlines, idim, itype;
    char ini_file[128];
    char *bound_opt[NOPT], *str;
    char *bbeg_label[] = {"part-X1-beg", "part-X2-beg","part-X3-beg"};
    char *bend_label[] = {"part-X1-end", "part-X2-end","part-X3-end"};

    for (itype = 0; itype < NOPT; itype++) {
        bound_opt[itype] = "0000";
    }
    bound_opt[OUTFLOW]      = "outflow";
    bound_opt[REFLECTIVE]   = "reflective";
    bound_opt[AXISYMMETRIC] = "axisymmetric";
    bound_opt[EQTSYMMETRIC] = "eqtsymmetric";
    bound_opt[PERIODIC]     = "periodic";
    bound_opt[SHEARING]     = "shearingbox";
    bound_opt[USERDEF]      = "userdef";
    bound_opt[POLARAXIS]    = "polaraxis";

#ifdef PARALLEL
    int j = 0;
    MPI_Status status;
    if (prank == 0) {
#endif

        sprintf(ini_file, "pluto.ini");
        printLog("> Reading %s (SETUP) ...\n\n", ini_file);

        nlines = ParamFileRead(ini_file);

        // number of particles
        particlesConf.num_particles = atoi(ParamFileGet("NPARTICLES", 1));

        /* ----------------------------- boundaries------------------------------ */
        // adapted from runtime_setup.c
        for (idim = 0; idim < 3; idim++){

            str = ParamFileGet(bbeg_label[idim], 1);
            COMPARE (str, bound_opt[itype], itype);
            if (itype == NOPT) {
                printf ("! ParticleReadConfig(): don't know how to put left particle boundary '%s'  \n", str);
                QUIT_PLUTO(1);
            }
            particlesConf.left_bound[idim] = itype;
        }

        for (idim = 0; idim < 3; idim++){

            str = ParamFileGet(bend_label[idim], 1);
            COMPARE (str, bound_opt[itype], itype);
            if (itype == NOPT) {
                printf ("! ParticleReadConfig(): don't know how to put right particle boundary '%s'  \n", str);
                QUIT_PLUTO(1);
            }
            particlesConf.right_bound[idim] = itype;
        }

        /* -------------------------- output --------------------------------- */
        particlesConf.output_dbl_dt = atof(ParamFileGet("dbl.pa", 1));
        particlesConf.output_dbl_dn = atoi(ParamFileGet("dbl.pa", 2));

        particlesConf.output_tab_dt = atof(ParamFileGet("tab.pa", 1));
        particlesConf.output_tab_dn = atoi(ParamFileGet("tab.pa", 2));

        particlesConf.output_h5_dt = atof(ParamFileGet("h5.pa", 1));
        particlesConf.output_h5_dn = atoi(ParamFileGet("h5.pa", 2));

        particlesConf.output_anl_dt = atof(ParamFileGet("analysis.pa", 1));
        particlesConf.output_anl_dn = atoi(ParamFileGet("analysis.pa", 2));

#ifdef PARALLEL
        for (j = 1; j < nprocs; j++) {
            MPI_Send(&particlesConf, 1, PartConfigType, j, 123, MPI_COMM_WORLD);
        }
    } else {
        MPI_Recv(&particlesConf, 1, PartConfigType, 0, 123, MPI_COMM_WORLD, &status);
    }
#endif
    NPARTICLES = particlesConf.num_particles;
    NOP = NPARTICLES;
}

/* ************************************************** */
//void ParticlesCheck(int nop, Runtime *ini)
/*
 * Decides if output or analysis are needed and call
 * the relevant routine
 *
 **************************************************** */
/*{
    for (i = 0; i < nprocs; i++) {
        for (j = 0; j < particles[i]; j++) {
            pl = &(particlesPerProc[i][j]);
            prevRank = pl->rank;
            check = ParticlesBoundary(uu, grid, pl);
            if (!(check))
            {
                for (h = j; h < particles[i] - 1; h++) {
                    particlesPerProc[i][h] = particlesPerProc[i][h + 1];
                }
                particles[i]--;
                continue;
            }

            ParticlesChkProcessor(pl);
            if (prevRank != pl->rank) {
                if (particles[pl->rank] < particleCapacity[pl->rank]) {
                    // The new processor has sufficient space to receive the particle
                    particlesPerProc[pl->rank][particles[pl->rank]] = particlesPerProc[i][j];
                    //TODO remove this line? it should do nothing
                    particlesPerProc[pl->rank][particles[pl->rank]].rank = pl->rank;
                    particles[pl->rank]++;
                } else {
                    particleCapacity[pl->rank] += particleCapacity[pl->rank] / 2;
                    appoggio = ARRAY_1D(particleCapacity[pl->rank], Particle);
                    for (h = 0; h < particles[pl->rank]; h++) {
                        appoggio[h] = particlesPerProc[pl->rank][h];
                    }
                    FreeArray1D((void *) particlesPerProc[pl->rank]);
                    particlesPerProc[pl->rank] = appoggio;
                    particlesPerProc[pl->rank][particles[pl->rank]++] = particlesPerProc[i][j];
                }
                // Remove the particle from the previous processor
                for (h = j; h < particles[prevRank] - 1; h++) {
                    particlesPerProc[prevRank][h] = particlesPerProc[prevRank][h + 1];
                }
                particles[prevRank]--;
                j--;
            }
        }
    }
}
*/
/* ************************************************** */
void ParticlesSave(int nop, Runtime *ini)
    /*
     * Decides if output or analysis are needed and call
     * the relevant routine
     *
     **************************************************** */
{
    char single_file[] = "single_file";
    static double dbl_out_time = 0.0, tab_out_time = 0.0, anl_out_time = 0.0;
    extern int prank;
    int i = 0, j = 0;
    int last_step;
#ifdef PARALLEL
    int totalParts = 0;
    int partsToBeReceived = 0;
    MPI_Status status;
    Particle * particlesArray = NULL;
    static int first_savecall = 1;
    if (first_savecall)
    {
        dbl_out_time = g_time;
        first_savecall = 0;
    }
    last_step = (fabs(g_time-ini->tstop) < 1.e-12 ? 1:0);
    if (dbl_out_time - g_time < g_dt || g_stepNumber == 0 || last_step) {

        MPI_Barrier(MPI_COMM_WORLD);

        /********* Memory Consumption Strategy ************/
        if (nop > 0) {
            if (prank == 0) {
                particlesArray = ARRAY_1D(nop, Particle);
                for (i = 0; i < nprocs; i++) {
                    partsToBeReceived = 0;
                    if (i != 0) {
                        MPI_Recv(&partsToBeReceived, 1, MPI_INT, i, 123, MPI_COMM_WORLD, &status);
                        if (partsToBeReceived > 0) {
                            MPI_Recv(&particlesArray[totalParts], partsToBeReceived, ParticleType, i, 123, MPI_COMM_WORLD, &status);
                            totalParts += partsToBeReceived;
                        }
                        for (j = 0; j < particles[i]; j++) {
                            particlesArray[totalParts++] = particlesPerProc[i][j];
                        }
                    } else {
                        for (j = 0; j < particles[i]; j++) {
                            particlesArray[totalParts++] = particlesPerProc[i][j];
                        }
                    }
                }
                if (particlesConf.output_dbl_dt > 0.0) {
                    if (dbl_out_time - g_time < g_dt || g_stepNumber == 0 || last_step) {
                        PartSaveDbl(particlesArray, totalParts, g_stepNumber, g_time, ini->output_dir);
                    }
                }

                if (particlesConf.output_dbl_dn > 0) {
                    if (g_stepNumber % particlesConf.output_dbl_dn == 0 || g_stepNumber == 0 || last_step, ini->output_dir) {
                        PartSaveDbl(particlesArray, totalParts, g_stepNumber, g_time, ini->output_dir);
                    }
                }

                if (particlesConf.output_tab_dt > 0.0) {
                    if (tab_out_time - g_time < g_dt || g_stepNumber == 0 || last_step) {
                        PartSaveVTK(particlesArray, totalParts, g_stepNumber, g_time, ini->output_dir);
                        tab_out_time += particlesConf.output_tab_dt;
                    }
                }

                FreeArray1D((void *)particlesArray);
            } else {
                for (i = 0; i < nprocs; i++) {
                    totalParts += particles[i];
                }
                particlesArray = ARRAY_1D(totalParts, Particle);
                totalParts = 0;
                for (i = 0; i < nprocs; i++) {
                    for (j = 0; j < particles[i]; j++) {
                        particlesArray[totalParts] = particlesPerProc[i][j];
                        totalParts++;
                    }
                }
                MPI_Send(&totalParts, 1, MPI_INT, 0, 123, MPI_COMM_WORLD);
                if (totalParts > 0) {
                    MPI_Send(particlesArray, totalParts, ParticleType, 0, 123, MPI_COMM_WORLD);
                }
                FreeArray1D((void *)particlesArray);
            }
        }
        dbl_out_time += particlesConf.output_dbl_dt;
    }


    /*************************************************/


#else
    if (particlesConf.output_dbl_dt > 0.0) {
        if (dbl_out_time - g_time < g_dt || g_stepNumber == 0 || last_step) {
            PartSaveDbl(particlesPerProc[0], nop, g_stepNumber, g_time, ini->output_dir);
            dbl_out_time += particlesConf.output_dbl_dt;
        }
    }

    if (particlesConf.output_dbl_dn > 0) {
        if (g_stepNumber % particlesConf.output_dbl_dn == 0 || g_stepNumber == 0 || last_step) {
            PartSaveDbl(particlesPerProc[0], nop, g_stepNumber, g_time, ini->output_dir);
        }
    }

    /* Tabulated Output */
    if (particlesConf.output_tab_dt > 0.0) {
        if (tab_out_time - g_time < g_dt || g_stepNumber == 0 || last_step) {
            PartSaveTab(particlesPerProc[0], nop, g_stepNumber, g_time, ini->output_dir);
            tab_out_time += particlesConf.output_tab_dt;
        }
    }

    if (particlesConf.output_tab_dn > 0) {
        if (g_stepNumber % particlesConf.output_tab_dn == 0 || g_stepNumber == 0 || last_step) {
            PartSaveTab(particlesPerProc[0], nop, g_stepNumber, g_time, ini->output_dir);
        }
    }

    if (particlesConf.output_anl_dt > 0.0 || g_stepNumber == 0 || last_step) {
        if (anl_out_time - g_time < g_dt) {
            ParticlesAnalysis(particlesPerProc[0], nop, g_stepNumber, g_time);
            anl_out_time += particlesConf.output_anl_dt;
        }
    }

    if (particlesConf.output_anl_dn > 0 || g_stepNumber == 0 || last_step) {
        if (g_stepNumber % particlesConf.output_anl_dn == 0) {
            ParticlesAnalysis(particlesPerProc[0], nop, g_stepNumber, g_time);
        }
    }
#endif
}

void calculateGravitation(double ****v, Grid *grid, Particle *pl)
{
    double r, r3, theta, phi, x, y, z;
    double d, d3, phiPlanet, muPlanet_eff;
    static double gm = 4. * CONST_PI * CONST_PI; // G*M_star
    double g_cart[DIMENSIONS], g_spher[DIMENSIONS];
    int i;

    r = pl->coor[0];
    theta = pl->coor[1];
    phi = pl->coor[2];
    r3 = r*r*r;

    x = r * sin(theta) * cos(phi);
    y = r * sin(theta) * sin(phi);
    z = r * cos(theta);

    if (g_time < g_inputParam[PLANETRAMPUPTIME]) {
        // slowly ramp up the planet mass in order to avoid strong shocks
        muPlanet_eff = UNIT_G * planet.M * sin(CONST_PI / 2. / g_inputParam[PLANETRAMPUPTIME] * g_time);
    } else {
        muPlanet_eff = UNIT_G * planet.M;
    }

    g_cart[0] = -gm * x / r3;
    g_cart[1] = -gm * y / r3;
    g_cart[2] = -gm * z / r3;

    d = sqrt((x - planet.x) * (x - planet.x) + (y - planet.y) * (y - planet.y) + z * z);
    d3 = d*d*d;

    g_cart[0] += -muPlanet_eff * (x - planet.x) / d3;
    g_cart[1] += -muPlanet_eff * (y - planet.y) / d3;
    g_cart[2] += -muPlanet_eff * z / d3;

    transformCartesianSpherical(g_cart, g_spher, pl->coor[0], pl->coor[1], pl->coor[2]);
    for(i=0;i<DIMENSIONS;i++) {
        pl->g_grav[i] = g_spher[i];
    }
}

/* ******************************************************************************* */
double calculateDrag(double ****v, Grid *grid, Particle *pl)
    /*
     * Calculate the drag acceleration onto the dust particles(g) and returns the
     * stopping time tstop.
     * TODO: - pure Stokes regime for reference
     *       - extend to other geometries
     *       - read values of radius a and mass m0 of the molecules from input file
     *       - calculate the geometric cross section sigma from a instead of using
     *         the fixed value
     *       - find a way to generalize the definition of the sound speed based on
     *         the law defined in eos.c
     ********************************************************************************* */
{
    static double velocity[DIMENSIONS];
    double density;
    double pressure;
    double tstop;
    double density_particle = 1.0/UNIT_DENSITY;
    double vtherm, vrel[DIMENSIONS], vrelabs;
    double Kn, Ma, Re, sigma, nu, Cd, CdE, CdS;
    double l, f, term, temper, m0, c_s, R;

    ParticlesLocate(pl->coor, pl->cell, grid);
    ParticlesInterpol(velocity, pl, v, grid, VX1);
    ParticlesInterpol(&density, pl, v, grid, RHO);
#if EOS == IDEAL
    ParticlesInterpol(&pressure, pl, v, grid, PRS);
#endif

#if GEOMETRY == POLAR
#if ROTATING_FRAME == YES
    velocity[1] += g_OmegaZ*pl->coor[0];
#endif
    DIM_EXPAND(vrel[0] = pl->mom[0]             - velocity[0]; ,
             vrel[1] = pl->mom[1]/pl->coor[0] - velocity[1]; ,
             vrel[2] = pl->mom[2]             - velocity[2];
             )
#elif GEOMETRY == SPHERICAL
        R = pl->coor[0]*sin(pl->coor[1]);
#if ROTATING_FRAME == YES
    velocity[2] += g_OmegaZ*R;
#endif
    DIM_EXPAND(vrel[0] = pl->mom[0]             - velocity[0]; ,
             vrel[1] = pl->mom[1]/pl->coor[0] - velocity[1]; ,
             vrel[2] = pl->mom[2]/R           - velocity[2];
             )
#endif
    vrelabs = DIM_EXPAND(vrel[0]*vrel[0], + vrel[1]*vrel[1], + vrel[2]*vrel[2]);
    vrelabs = sqrt(vrelabs);

#if EOS == IDEAL
    c_s = sqrt(g_gamma*pressure/density);
#elif EOS == ISOTHERMAL
    c_s = g_isoSoundSpeed/sqrt(R);
#endif
    vtherm = sqrt(8.0/CONST_PI)*c_s;
    // mean free path for molecular hydrogen
    // (see Haghighipour & Boss, 2003 eq. 20)
    l = (4.72e-9/(density*UNIT_DENSITY))/UNIT_LENGTH;
    f = pl->radius/(pl->radius+l);

    // a0 = 1.5e-8 cm for molecular hydrogen
    //sigma = CONST_PI*pow(1.5e-8/UNIT_LENGTH,2.0);
    sigma = 2.e-15/UNIT_LENGTH/UNIT_LENGTH;
    m0 = CONST_mH / (UNIT_DENSITY * UNIT_LENGTH * UNIT_LENGTH * UNIT_LENGTH);
    nu = 1.0/3.0*m0*vtherm/sigma;

    Kn = 0.5*l/pl->radius;
    Ma = vrelabs/c_s;
    Re = 2.0*pl->radius*density*vrelabs/nu;

    /************************************************
      Lyra et al., 2009 approach
pro: - valid also for high Mach numbers
- retain the exact value of the Epstein
stopping time for small size particles
cons:- too many operations
     ************************************************/
    CdE = 2.0*sqrt(Ma*Ma+128.0/9.0/CONST_PI);
    if (Re <= 1.e-3) {
        CdS = 24.0*nu/(2.0*pl->radius*density*c_s) +
            3.6/c_s*pow(vrelabs,0.687)*pow(2.0*pl->radius*density/nu,-0.313);
    } else if (Re <= 500.0) {
        CdS = 24.0*Ma/Re + 3.6*Ma*pow(Re,-0.313);
    } else if (Re <= 1500.0) {
        CdS = Ma*9.5e-5*pow(Re,1.397);
    } else {
        CdS = Ma*2.61;
    }

    Cd = (9.0*Kn*Kn*CdE + CdS)/(3.0*Kn+1.0)/(3.0*Kn+1.0);
    tstop = 4.0*l*density_particle/(3.0*density*Cd*c_s*Kn);

    /************************************************
      Nader & Boss, 2003
cons: - not valid for high Mach numbers
- not converging exactly to the Epstein
stopping time for small sizes.
pro:  - already tested on FARGO.
     ************************************************
     if (Re <= 1.0) {
     Cd = 24.0/Re;
     } else if (Re <= 800.0) {
     Cd = 24.0*pow(Re,-0.6);
     } else {
     Cd = 0.44;
     }

     term = density*((1.0-f)*vtherm + 3.0/8.0*f*Cd*vrelabs);
     tstop  = pl->radius*density_particle/term;
     */
    /************************************************
      Epstein regime
cons: - valid only for small size particles
pro:  - very fast
- good for comparison
     *************************************************
     tstop  = (pl->radius*density_particle)/(density*vtherm);
     */
    DIM_EXPAND(pl->g_drag[0] = 0.0; ,
            pl->g_drag[1] = 0.0; ,
            pl->g_drag[2] = 0.0;);

#if GEOMETRY == POLAR
    pl->g_drag[0] = -(pl->mom[0] - velocity[0]            )/tstop;
    pl->g_drag[1] = -(pl->mom[1] - velocity[1]*pl->coor[0])/tstop;
    pl->g_drag[2] = -(pl->mom[2] - velocity[2]            )/tstop;
#elif GEOMETRY == SPHERICAL
    pl->g_drag[0] = -(pl->mom[0] - velocity[0]            )/tstop;
    pl->g_drag[1] = -(pl->mom[1] - velocity[1]*pl->coor[0])/tstop;
    pl->g_drag[2] = -(pl->mom[2] - velocity[2]*R          )/tstop;
#endif

    return tstop;
}
