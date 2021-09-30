#include "pluto.h"
#include "particles.h"


void PartSaveVTK(Particle * myParticles, int nop, int nstep, double glob_time, char *output_dir)
{
  /* save particles in ascii format with number of particles,step and output time
    as an unstrcutured grid in a .vtu file*/

  static int called = 0;
  int i = 0, j = 0;
  FILE *stream;
  Particle *pl;
  char filename[128];

  printLog("> Writing Particles file #%d (vtu) to disk...\n", called);

  sprintf(filename, "%s/part_data.%04d.vtu", output_dir, called);
  stream = fopen(filename, "w");
  if(myParticles != NULL) {
   fprintf(stream, "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
   fprintf(stream, "<UnstructuredGrid>\n");
   fprintf(stream, "<Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n",nop, nop);
   fprintf(stream, "<Points>\n");
   fprintf(stream, "<DataArray name=\"Position\" type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n");

  for(j = 0; j < nop; j++){
    pl = &myParticles[j];
    fprintf(stream, " %f %f %f", pl->coor[0], pl->coor[1], 0.0);
   }
  fprintf(stream, "\n");
  fprintf(stream, "</DataArray>\n");
  fprintf(stream, "</Points>\n");
  fprintf(stream, "<PointData>\n");
  fprintf(stream, "<DataArray type=\"Float32\" Name=\"Velocity\" NumberOfComponents=\"3\" format=\"ascii\">\n");

  for(j = 0; j < nop; j++) {
       pl = &myParticles[j];
       fprintf(stream, " %f %f %f", pl->mom[0], pl->mom[1], 0.0);
  }
  fprintf(stream, "\n");
  fprintf(stream, "</DataArray>\n");
  fprintf(stream, "<DataArray type=\"Int32\" Name=\"Identity\" format=\"ascii\">\n");

  for(j = 0; j < nop; j++) {
       pl = &myParticles[j];
       fprintf(stream, " %d", pl->identity);
     }
  fprintf(stream, "\n");
  fprintf(stream, "</DataArray>\n");
  }
  fprintf(stream, "</PointData>\n");
  fprintf(stream, "<Cells>\n");
  fprintf(stream, "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n");
  for(j = 0; j < nop+1; j++) {
     fprintf(stream, " %d", j );
  }
  fprintf(stream, "</DataArray>\n");
  fprintf(stream, "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n");
  for(j = 1; j < nop+1; j++) {
     fprintf(stream, " %d",j);
  }
  fprintf(stream, "</DataArray>\n");
  fprintf(stream, "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n");
  for(j = 0; j < nop; j++) {
     fprintf(stream, " %u",1);
  }
  fprintf(stream, "</DataArray>\n");
  fprintf(stream, "</Cells>\n");
  fprintf(stream, "</Piece>\n");
  fprintf(stream, "</UnstructuredGrid>\n");
  fprintf(stream, "</VTKFile>\n");
  fclose(stream);
  called += 1;
}




void PartSaveDbl(Particle * myParticles, int nop, int nstep, double glob_time, char *output_dir) {

  /* save structure particle, with time, step and number of particles */
  /* PREFERRED ROUTINE TO SAVE PARTICLES */

  static int called = 0;
  int i = 0, j = 0;
  FILE *stream;
  Particle *pl;
  char filename[128];
#ifdef PARALLEL
  extern int prank;
#endif
  printLog("> Writing Particles file #%d (dbl) to disk...\n", called);

  sprintf(filename, "%s/part_data.%04d.dbl", output_dir, called);
/*
  print("BEFORE FILE_OPEN EXECUTION\n");
*/
  /*stream = FILE_OPEN (filename, sizeof(double), "w");*/
  stream = fopen(filename,"w");

  fwrite(&nop, sizeof (int), 1, stream);
  fwrite(&nstep, sizeof (int), 1, stream);
  fwrite(&glob_time, sizeof (double), 1, stream);

  if(myParticles != NULL) {
    /*
    print("PRANK: %d - PartSaveDbl - myParticles != NULL NOP: %d\n",prank,nop);
    */

    for(j = 0; j < nop; j++) {

      pl = &myParticles[j];

      fwrite(&(pl->identity), sizeof (int), 1, stream);

      for (i = 0; i < DIMENSIONS; ++i) {
        fwrite(&(pl->cell[i]), sizeof (int), 1, stream);
      }

      for (i = 0; i < DIMENSIONS; ++i) {
        fwrite(&(pl->coor[i]), sizeof (double), 1, stream);
      }

      for (i = 0; i < DIMENSIONS; ++i) {
        fwrite(&(pl->mom[i]), sizeof (double), 1, stream);
      }

      fwrite(&(pl->tstop), sizeof (double), 1, stream);

      fwrite(&(pl->radius), sizeof (double), 1, stream);
    }
  }

  /*FILE_CLOSE (stream, sizeof(double));*/
  fclose(stream);

  called += 1;

}

void PartLoadDbl(Particle * myParticles, int nop, int nstep, char *data_dir) {

  /* save structure particle, with time, step and number of particles */
  /* PREFERRED ROUTINE TO SAVE PARTICLES */

  static int called = 0;
  int i = 0, j = 0;
  FILE *stream;
  Particle *pl;
  char filename[128];
  double dump;
#ifdef PARALLEL
  extern int prank;
#endif
  called = nstep;
  printLog("> Reading Particles file #%d (dbl) to disk...\n", called);

  sprintf(filename, "%s/part_data.%04d.dbl", data_dir, called);
/*
  print("BEFORE FILE_OPEN EXECUTION\n");
*/
  /*stream = FILE_OPEN (filename, sizeof(double), "w");*/
  stream = fopen(filename,"r");

  fread(&nop, sizeof (int), 1, stream);
  fread(&nstep, sizeof (int), 1, stream);
  fread(&dump, sizeof (double), 1, stream);

  if(myParticles != NULL) {
    /*
    print("PRANK: %d - PartSaveDbl - myParticles != NULL NOP: %d\n",prank,nop);
    */

    for(j = 0; j < nop; j++) {

      pl = &myParticles[j];

      fread(&(pl->identity), sizeof (int), 1, stream);

      for (i = 0; i < DIMENSIONS; ++i) {
        fread(&(pl->cell[i]), sizeof (int), 1, stream);
      }

      for (i = 0; i < DIMENSIONS; ++i) {
        fread(&(pl->coor[i]), sizeof (double), 1, stream);
      }

      for (i = 0; i < DIMENSIONS; ++i) {
        fread(&(pl->mom[i]), sizeof (double), 1, stream);
      }
      
      fread(&(pl->tstop), sizeof (double), 1, stream);
      
      fread(&(pl->radius), sizeof (double), 1, stream);
    }
  }

  /*FILE_CLOSE (stream, sizeof(double));*/
  fclose(stream);

  called += 1;

}


void PartSaveTab(Particle * particles, int nop, int nstep, double glob_time, char *output_dir) {

  /* save particles in ascii format with number of particles,step and output time */

  static int called = 0;
  int i = 0, j = 0;
  FILE *stream;
  Particle *pl;
  char filename[128];

  printLog("> Writing Particles file #%d (tab) to disk...\n", called);

  sprintf(filename, "%s/part_data.%04d.tab", output_dir, called);
  stream = fopen(filename, "w");

  fprintf(stream, "%d  %d  %lf", nop, nstep, glob_time);
  fprintf(stream, "\n");

  if(particles != NULL) {
    for(j = 0; j < nop; j++) {
      pl = &particles[j];

      fprintf(stream, " %d", pl->identity);

      for (i = 0; i < DIMENSIONS; ++i) {
        fprintf(stream, " %lf", pl->coor[i]);
      }
      for (i = 0; i < DIMENSIONS; ++i) {
        fprintf(stream, " %lf", pl->mom[i]);
      }

      fprintf(stream, "\n");
    }
  }

  fclose(stream);
  called += 1;
}
