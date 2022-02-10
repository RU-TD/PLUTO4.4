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
void calculate_JFlux(Grid* grid, const Data* data);

/* ********************************************************************* */
void Init (double *us, double x1, double x2, double x3)
/*
 *
 *********************************************************************** */
{
  int nv;
  double ngas;
  double mu = 2.;
  #if EOS == IDEAL
   g_gamma = g_inputParam[GAMMA_EOS];
  #endif

  us[RHO] = mu*g_inputParam[RHO];

  us[VX1] = 0.0;

  #if EOS == IDEAL
   us[PRS] = us[RHO]*g_inputParam[TEMP]/(mu*KELVIN);
  #endif

  // initial chemical conditions
  NTRACER_LOOP(nv) us[nv] = 0.0;
  us[IDX_CHEM_H2] = 0.5;
  us[IDX_CHEM_O] = 3e-4;
  us[IDX_CHEM_C] = 1e-4;
  us[IDX_CHEM_Hej] = 1e-1;
  us[IDX_CHEM_E] = us[IDX_CHEM_Hej];
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
   #if CHEMISTRY != NO
    double Gnot = 1e1 / 1.15;
    double pi = acos(-1e0);
    for (int k = KBEG; k <= KEND; k++)
    {
      for (int j = JBEG; j <= JEND; j++)
      {
        // set Draine multifrequency intensity
        prizmo_get_draine_flux_c(irradiation.jflux[KBEG][JBEG][IBEG]);
        // loop to multiply 2*pi*G0
        for(int l = 0; l < NPHOTO; l++){
          irradiation.jflux[KBEG][JBEG][IBEG][l] *= Gnot * 2e0 * pi;
        }
      }
    }
    calculate_JFlux(grid,d); // calculating JFlux at each cell
   #endif
  }
  QUIT_PLUTO(0);
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

        ldi.x1_begin = grid[IDIR].beg;
        ldi.x1_end = grid[IDIR].end;

        ldi.x2_begin = grid[JDIR].beg;
        ldi.x2_end = grid[JDIR].end;

        ldi.x3_begin = grid[KDIR].beg;
        ldi.x3_end = grid[KDIR].end;


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
        int k, j, i, l, n;
        double column_density = 0.;
        double density = 0.0, dr = 0.0;
        double mpart = g_inputParam[GAMMA_EOS]*CONST_amu;
        MPI_Status status;

        for (k = KBEG; k <= KEND; k++)
        {
                for (j = JBEG; j <= JEND; j++)
                {
                        column_density = 0.0;
                        irradiation.column_density[k][j][IBEG][val] = column_density;
                        for (i = IBEG; i <= IEND; i++)
                        {
                                density = d->Vc[RHO][k][j][i] * UNIT_DENSITY;
                                dr = grid[IDIR].dx[i]*UNIT_LENGTH;
				if(val == 0) {
                                	column_density += density/mpart*dr;
				} else if (val == 1) {
					column_density += d->Vc[IDX_CHEM_H2][k][j][i]*density/mpart*dr;
				} else if (val == 2) {
					column_density += d->Vc[IDX_CHEM_CO][k][j][i]*density/mpart*dr;
				} else if (val == 3) {
					column_density += d->Vc[IDX_CHEM_H2O][k][j][i]*density/mpart*dr;
				} else {
					column_density += (2.*d->Vc[IDX_CHEM_H2][k][j][i] + d->Vc[IDX_CHEM_H][k][j][i])*density/mpart*dr;
				}
                                irradiation.column_density[k][j][i+1][val] = column_density;
                        }
                }
        }

        if(irradiation.neighbour.receive_rank != irradiation.neighbour.send_rank)//check if more than one process is in the communicator (nproc > 1)
        {
                if(irradiation.neighbour.receive_rank == -1 && irradiation.neighbour.send_rank != -1) //no neighbor in between the star and the current domain
                {
                        n = 0;
                        i = IEND+1;
                        for (k = KBEG; k <= KEND; k++)
                        {
                                for (j = JBEG; j <= JEND; j++)
                                {
                                        irradiation.data_buffer[n] = irradiation.column_density[k][j][i][val];
                                        irradiation.column_density_offset[n] = 0.;
                                        n++;
                                }
                        }
                        MPI_Send(irradiation.data_buffer,sizeof(double) * (NX2*NX3),MPI_CHAR,irradiation.neighbour.send_rank,0,MPI_COMM_WORLD);
                }
                else
                {
                        MPI_Recv(irradiation.column_density_offset ,sizeof(double) * (NX2*NX3), MPI_CHAR, irradiation.neighbour.receive_rank, 0, MPI_COMM_WORLD, &status);

                        n = 0;
                        i = IEND+1;
                        for (k = KBEG; k <= KEND; k++)
                        {
                                for (j = JBEG; j <= JEND; j++)
                                {
                                        irradiation.data_buffer[n] = irradiation.column_density_offset[n] + irradiation.column_density[k][j][i][val];
                                        n++;
                                }
                        }

                        if(irradiation.neighbour.send_rank != -1)
                        {
                                MPI_Send(irradiation.data_buffer,sizeof(double) * (NX2*NX3),MPI_CHAR,irradiation.neighbour.send_rank,0,MPI_COMM_WORLD);
                        }

                }
                MPI_Barrier(MPI_COMM_WORLD);
        }

}

void calculate_ColumnDensity(Grid* grid, const Data* d)
{
        int k, j, i, l, n;
	
	for (l=0; l<5; l++)
	{
        	calculate_ColumnDensity_perDomain(grid, d, l);

        	n = 0;
        	for (k = KBEG; k <= KEND; k++)
        	{
                	for (j = JBEG; j <= JEND; j++)
                	{
                        	for (i = IBEG; i <= IEND; i++)
                        	{
                                	irradiation.column_density[k][j][i][l] += irradiation.column_density_offset[n];
                        	}
                        	n++;
                	}
        	}
	}
}

void calculate_JFlux_perDomain(Grid* grid, const Data* d, int val)
{
        int k, j, i, l, n;
        double jfluxes[NPHOTO] = {0e0};
	double x[NTRACER];
        double density = 0.0, dr = 0.0, Tgas = 0.0;
	double mu=2.;
	double Gnot = 1e1 / 1.15;
	double pi = acos(-1e0);
        MPI_Status status;
       #ifdef PARALLEL
	MPI_Comm_rank (MPI_COMM_WORLD, &prank);
       #endif

        for (k = KBEG; k <= KEND; k++)
        {
                for (j = JBEG; j <= JEND; j++)
                {
                        for (i = IBEG; i <= IEND; i++)
                        {
			      NTRACER_LOOP(n) x[n-TRC] = d->Vc[n][k][j][i];
			      density = d->Vc[RHO][k][j][i] * UNIT_DENSITY;
                              dr = grid[IDIR].dx[i]*UNIT_LENGTH;
			      Tgas = d->Vc[PRS][k][j][i]/d->Vc[RHO][k][j][i]*(KELVIN*mu);
			      // attenuate radiation in the grid cell
			      for (n = 0; n<NPHOTO; n++) {
			            jfluxes[n] = irradiation.jflux[k][j][i][n];
			      }
			      prizmo_attenuate_rho_c(jfluxes, x, &Tgas, &dr, &density);
                              irradiation.jflux[k][j][i+1][val] = jfluxes[val];
                        }
                }
        }
/*
        if(irradiation.neighbour.receive_rank != irradiation.neighbour.send_rank)//check if more than one process is in the communicator (nproc > 1)
        {
                if(irradiation.neighbour.receive_rank == -1 && irradiation.neighbour.send_rank != -1) //no neighbor in between the star and the current domain
                {
                        n = 0;
                        i = IEND+1;
                        for (k = KBEG; k <= KEND; k++)
                        {
                                for (j = JBEG; j <= JEND; j++)
                                {
                                        irradiation.jflux_buffer[n] = irradiation.jflux[k][j][i][val];
                                        irradiation.jflux_offset[n] = 0.;
					n++;
                                }
                        }
                        MPI_Send(irradiation.jflux_buffer,sizeof(double) * (NX2*NX3),MPI_CHAR,irradiation.neighbour.send_rank,0,MPI_COMM_WORLD);
                }
                else
                {
                        MPI_Recv(irradiation.jflux_offset ,sizeof(double) * (NX2*NX3), MPI_CHAR, irradiation.neighbour.receive_rank, 0, MPI_COMM_WORLD, &status);

                        n = 0;
                        i = IEND+1;
                        for (k = KBEG; k <= KEND; k++)
                        {
                                for (j = JBEG; j <= JEND; j++)
                                {
                                        irradiation.jflux_buffer[n] = irradiation.jflux_offset[n] + irradiation.jflux[k][j][i][val];
					n++;
                                }
                        }

                        if(irradiation.neighbour.send_rank != -1)
                        {
                                MPI_Send(irradiation.jflux_buffer,sizeof(double) * (NX2*NX3),MPI_CHAR,irradiation.neighbour.send_rank,0,MPI_COMM_WORLD);
                        }

                }
                MPI_Barrier(MPI_COMM_WORLD);
        }
*/
}


void calculate_JFlux(Grid* grid, const Data* d)
{
        int k, j, i, l, n;

	for (l = 0; l <= NPHOTO; l++)
        {
              calculate_JFlux_perDomain(grid, d, l);
	      
              n = 0;
              for (k = KBEG; k <= KEND; k++)
              {
                    for (j = JBEG; j <= JEND; j++)
                    {
                          for (i = IBEG; i <= IEND; i++)
                          {
				if (l == 0) print("PRIMA %d %d %le %le\n",prank,i,irradiation.jflux[k][j][i][0],irradiation.jflux_offset[n]);
                                //irradiation.jflux[k][j][i][l] += irradiation.jflux_offset[n];
				if (l == 0) print("DOPO %d %d %le\n",prank,i,irradiation.jflux[k][j][i][l]);
                          }
                          n++;
                    }
              }
	      
	}
}

