/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief PRD Region.

  \author G. Picogna T. Grassi (picogna@usm.lmu.de)
  \date   Dec 15, 2020

  \b References
     - Grassi, T. et al. (2020). 
       "Modelling thermochemical processes in protoplanetary discs I: numerical methods".       
       MNRAS, Volume 494, Issue 3, pp.4471-4491

*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

void calculate_ColumnDensity(Grid* grid, const Data* data);
void attenuate_jflux(const Data* data, Grid* grid);
void read_jflux();
extern void prizmo_get_rho_c(double *, double *);
extern void prizmo_n2frac_c(double *, double *);
extern void prizmo_frac2n_c(double *, double *, double *);
extern void prizmo_rt_rho_c(double *, double *, double *, double *, double *);

/* ********************************************************************* */
void Init (double *us, double x1, double x2, double x3)
/*
 *
 *********************************************************************** */
{
  int nv;
  double ngas=g_inputParam[N_GAS];
  double rhod;
  double x[NTRACER], n[NTRACER];
  #if EOS == IDEAL
   g_gamma = g_inputParam[GAMMA_EOS];
  #endif

  us[VX1] = 0.0;

  NTRACER_LOOP(nv) {
    n[nv-TRC] = 0.;
    x[nv-TRC] = 0.;
  }

  n[IDX_CHEM_H2-TRC] = ngas/2.;
  n[IDX_CHEM_C-TRC] = ngas*1.e-4;
  n[IDX_CHEM_O-TRC] = ngas*3.e-4;
  n[IDX_CHEM_He-TRC] = ngas*1.e-1;

  // convert to mass fraction (and also get rho the total mass density)
  prizmo_get_rho_c(n, &rhod);
  prizmo_n2frac_c(n, x);

  NTRACER_LOOP(nv) {
    us[nv] = x[nv-TRC];
  }

  us[RHO] = rhod/(UNIT_DENSITY);
  #if EOS == IDEAL
   us[PRS] = us[RHO]*g_inputParam[TEMP]/(g_inputParam[MU]*KELVIN);
  #endif
}

/* ********************************************************************* */
void InitDomain (Data *d, Grid *grid)
/*! 
 * Assign initial condition by looping over the computational domain.
 * Called after the usual Init() function to assign initial conditions
 * on primitive variables.
 * Value assigned here will overwrite those prescribed during Init().
 *
 *
 *********************************************************************** */
{
  int i,j,k,n;
  read_jflux();
  for (int n=0; n<3; n++){
    DOM_LOOP(k,j,i){
      irradiation.column_density[n][k][j][i] = 0.0;
    }
  }
}

/* ********************************************************************* */
void Analysis (const Data *d, Grid *grid)
/* 
 *
 *
 *********************************************************************** */
{
}
/* ********************************************************************* */
void UserDefBoundary (const Data *d, RBox *box, int side, Grid *grid) 
/*
 *
 *********************************************************************** */
{
  if (side == 0) {    /* -- check solution inside domain -- */
    calculate_ColumnDensity(grid,d); // calculating the column density at each cell
    attenuate_jflux(d,grid);
  }
}

void read_jflux()
{
  FILE *fout;
  double Jscale = 1.e1;

  // load radiation at 1 AU
  fout = fopen("runtime_data/radiation_field.dat", "r");
  for (int n=0; n<NPHOTO; n++){
    fscanf(fout, "%le", &irradiation.jflux0[n]);
    irradiation.jflux0[n] *= 2.*CONST_PI*Jscale;
  }
  fclose(fout);
}

void attenuate_jflux(const Data *d, Grid *grid)
{
    int i,j,k,n;
    double x[NTRACER];
    double dx,Tgas,rho;
    double jflux[NPHOTO];
   
    for (int n=0; n<NPHOTO; n++){
        jflux[n] = irradiation.jflux0[n];
        irradiation.jflux[IBEG][n] = irradiation.jflux0[n];
    }

    DOM_LOOP(k,j,i){
        dx = grid->dx[IDIR][i]*UNIT_LENGTH;
        Tgas = d->Vc[PRS][k][j][i]/d->Vc[RHO][k][j][i]*(KELVIN*g_inputParam[MU]);
        rho = d->Vc[RHO][k][j][i]*UNIT_DENSITY;;
        NTRACER_LOOP(n) x[n-TRC] = d->Vc[n][k][j][i];
        prizmo_rt_rho_c(x,&rho,&Tgas,jflux,&dx);
        for (int n=0; n<NPHOTO; n++){
            irradiation.jflux[i+1][n] = jflux[n];
        }
    }
}

void find_CommunicationNeighbour(int current_rank, LocalDomainInfo *domain_info_array, int nproc, CommunicationNeighbour* cn)
{
    int n = 0;
    LocalDomainInfo* cdi = &domain_info_array[current_rank];

    cn->receive_rank = -1;
    cn->send_rank = -1;

    for(n = 0; n < nproc; ++n)
    {
        if(n != current_rank)
        {
            if(cdi->x1_begin - 1 == domain_info_array[n].x1_end && cdi->x2_begin == domain_info_array[n].x2_begin && cdi->x2_end == domain_info_array[n].x2_end && cdi->x3_begin == domain_info_array[n].x3_begin && cdi->x3_end == domain_info_array[n].x3_end)
            {
                cn->receive_rank = n;
            }
            else if(cdi->x1_end + 1 == domain_info_array[n].x1_begin && cdi->x2_begin == domain_info_array[n].x2_begin && cdi->x2_end == domain_info_array[n].x2_end && cdi->x3_begin == domain_info_array[n].x3_begin && cdi->x3_end == domain_info_array[n].x3_end)
            {
                cn->send_rank = n;
            }
        }
    }
}

void find_CommunicationNeighbours(Grid *grid)
{
    int nproc = 0, n = 0;
    LocalDomainInfo ldi;
    LocalDomainInfo *domain_info_array = NULL;
    LocalDomainInfo *item = NULL;
    CommunicationNeighbour* neighbour_info_array = NULL;
    MPI_Status status;

    MPI_Comm_size(MPI_COMM_WORLD, &nproc);

    ldi.x1_begin = grid->beg[IDIR];
    ldi.x1_end = grid->end[IDIR];

    ldi.x2_begin = grid->beg[JDIR];
    ldi.x2_end = grid->end[JDIR];

    ldi.x3_begin = grid->beg[KDIR];
    ldi.x3_end = grid->end[KDIR];

    if(prank == 0)
    {
        domain_info_array = malloc(sizeof(LocalDomainInfo) * nproc);
        neighbour_info_array = malloc(sizeof(CommunicationNeighbour) * nproc);
    }

    MPI_Gather(&ldi, sizeof(LocalDomainInfo), MPI_CHAR, domain_info_array, sizeof(LocalDomainInfo), MPI_CHAR, 0, MPI_COMM_WORLD);

    if(prank == 0)
    {
        for(n = 0; n < nproc; ++n)
        {
            find_CommunicationNeighbour(n,domain_info_array,nproc,&neighbour_info_array[n]);
        }

        for(n = 0; n < nproc; ++n)
        {
            if(n != 0)
            {
                MPI_Send(&neighbour_info_array[n],sizeof(CommunicationNeighbour),MPI_CHAR,n,0,MPI_COMM_WORLD);
            }
        }

        irradiation.neighbour.receive_rank = neighbour_info_array[0].receive_rank;
        irradiation.neighbour.send_rank = neighbour_info_array[0].send_rank;

        free(domain_info_array);
        free(neighbour_info_array);
    }
    else
    {
        MPI_Recv(&irradiation.neighbour, sizeof(CommunicationNeighbour), MPI_CHAR, 0, 0, MPI_COMM_WORLD, &status);
    }
}

void calculate_ColumnDensity_perDomain(Grid* grid, const Data* d, int val)
{
    int k, j, i, l;
    double column_density = 0.;
    double density = 0.0, dr = 0.0;
    double mpart = g_inputParam[GAMMA_EOS]*CONST_amu;
    double x[NTRACER], n[NTRACER];
    MPI_Status status;

    NTRACER_LOOP(l){
        x[l-TRC] = 0.;
        n[l-TRC] = 0.;
    }

    KDOM_LOOP(k){
        JDOM_LOOP(j){
            column_density = 0.0;
            irradiation.column_density[val][k][j][IBEG] = column_density;
            IDOM_LOOP(i){
                density = d->Vc[RHO][k][j][i] * UNIT_DENSITY;
                NTRACER_LOOP(l) x[l-TRC] = d->Vc[l][k][j][i];
                prizmo_frac2n_c(x, &density, n);
                dr = grid->dx[IDIR][i]*UNIT_LENGTH;
                if(val == 0) {
                    column_density += density/mpart*dr;
                } else if (val == 1) {
                    column_density += n[IDX_CHEM_H2-TRC] * dr;
                } else {
                    column_density += n[IDX_CHEM_CO-TRC] * dr;
                }
                irradiation.column_density[val][k][j][i+1] = column_density;
            }
        }
    }

    if(irradiation.neighbour.receive_rank != irradiation.neighbour.send_rank)//check if more than one process is in the communicator (nproc > 1)
    {
        if(irradiation.neighbour.receive_rank == -1 && irradiation.neighbour.send_rank != -1) //no neighbor in between the star and the current domain
        {
            l = 0;
            i = IEND+1;
            KDOM_LOOP(k){
                JDOM_LOOP(j){
                    irradiation.data_buffer[l] = irradiation.column_density[val][k][j][i];
                    irradiation.column_density_offset[l] = 0.;
                    l++;
                }
            }
            MPI_Send(irradiation.data_buffer,NX2*NX3,MPI_DOUBLE,irradiation.neighbour.send_rank,0,MPI_COMM_WORLD);
        }
        else
        {
            MPI_Recv(irradiation.column_density_offset ,NX2*NX3, MPI_DOUBLE, irradiation.neighbour.receive_rank, 0, MPI_COMM_WORLD, &status);

            l = 0;
            i = IEND+1;
            KDOM_LOOP(k){
		JDOM_LOOP(j){
                    irradiation.data_buffer[l] = irradiation.column_density_offset[l] + irradiation.column_density[val][k][j][i];
                    l++;
                }
            }

            if(irradiation.neighbour.send_rank != -1)
            {
                MPI_Send(irradiation.data_buffer,NX2*NX3,MPI_DOUBLE,irradiation.neighbour.send_rank,0,MPI_COMM_WORLD);
            }

        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
}

void calculate_ColumnDensity(Grid* grid, const Data* d)
{
        int k, j, i, l, n;
	
	for (l=0; l<3; l++)
	{
        	calculate_ColumnDensity_perDomain(grid, d, l);
/*
        	n = 0;
        	for (k = KBEG; k <= KEND; k++)
        	{
                	for (j = JBEG; j <= JEND; j++)
                	{
                        	for (i = IBEG; i <= IEND; i++)
                        	{
                                	irradiation.column_density[l][k][j][i] += irradiation.column_density_offset[n];
                        	}
                        	n++;
                	}

		}
*/
	}
}
