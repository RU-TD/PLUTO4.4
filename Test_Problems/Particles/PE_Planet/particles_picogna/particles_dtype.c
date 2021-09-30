#include "pluto.h"
#include "particles.h"


#ifdef PARALLEL
void PARTICLE_CONFIG_DATATYPE() {
  int i = 0;
  struct PARTICLES_CONFIG pc[1];
  MPI_Datatype type[11] = {MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_DOUBLE, MPI_DOUBLE, 
    MPI_DOUBLE, MPI_DOUBLE};
  int blocklengths[11] = {1, DIMENSIONS, DIMENSIONS, 1, 1, 1, 1, 1, 1, 1, 1};
  MPI_Aint displ[11];
  MPI_Aint base;

  MPI_Get_address(&pc[0], &base);
  MPI_Get_address(&pc[0].num_particles, &displ[0]);
  MPI_Get_address(&pc[0].left_bound, &displ[1]);
  MPI_Get_address(&pc[0].right_bound, &displ[2]);
  MPI_Get_address(&pc[0].output_dbl_dn, &displ[3]);
  MPI_Get_address(&pc[0].output_tab_dn, &displ[4]);
  MPI_Get_address(&pc[0].output_h5_dn, &displ[5]);
  MPI_Get_address(&pc[0].output_anl_dn, &displ[6]);
  MPI_Get_address(&pc[0].output_dbl_dt, &displ[7]);
  MPI_Get_address(&pc[0].output_tab_dt, &displ[8]);
  MPI_Get_address(&pc[0].output_h5_dt, &displ[9]);
  MPI_Get_address(&pc[0].output_anl_dt, &displ[10]);
  for (i = 0; i < 11; i++) {
    displ[i] -= base;
  }
  MPI_Type_create_struct(11, blocklengths, displ, type, &PartConfigType);
  MPI_Type_commit(&PartConfigType);
}


void PARTICLE_STRUCT_DATATYPE() {
  int i = 0;
  Particle pl[1];

  MPI_Datatype type[11] = {MPI_INT, MPI_INT, MPI_INT,
      MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE,
      MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_CHAR};
  int blocklengths[11] = {1, 1, DIMENSIONS, 1, DIMENSIONS,
      DIMENSIONS, DIMENSIONS, 1, DIMENSIONS, DIMENSIONS, 1};

  MPI_Aint displ[11];

  MPI_Aint base;

  MPI_Get_address(&pl[0], &base);
  MPI_Get_address(&pl[0].identity, &displ[0]);
  MPI_Get_address(&pl[0].rank, &displ[1]);
  MPI_Get_address(&pl[0].cell, &displ[2]);
  MPI_Get_address(&pl[0].radius, &displ[3]);
  MPI_Get_address(&pl[0].coor, &displ[4]);
  MPI_Get_address(&pl[0].coor_old, &displ[5]);
  MPI_Get_address(&pl[0].mom, &displ[6]);
  MPI_Get_address(&pl[0].tstop, &displ[7]);
  MPI_Get_address(&pl[0].g_grav, &displ[8]);
  MPI_Get_address(&pl[0].g_drag, &displ[9]);
  MPI_Get_address(&pl[0].flag, &displ[10]);

  for (i = 0; i < 11; i++) {
      displ[i] -= base;
  }

  MPI_Type_create_struct(11, blocklengths, displ, type, &ParticleType);
  MPI_Type_commit(&ParticleType);
}
#endif
