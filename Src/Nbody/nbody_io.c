#include "pluto.h"

int nbodyReadPlanetFile(const char *filename)
{
    int i;
    char *line_buffer = NULL;
    FILE *f;

    /* open file */
    if ((f = fopen(filename, "r")) == NULL)
    {
        fprintf(stderr, "Error in read_planet_file: Could not open %s\n", filename);
        goto read_planet_file_error;
    }

    /* allocating memory for line_buffer */
    if ((line_buffer = malloc(BUFFER_LENGTH * sizeof(line_buffer))) == NULL)
    {
        fprintf(stderr, "Error in read_planet_file: Could not allocate memory for line buffer\n");
        goto read_planet_file_error;
    }
    
    i = CENTRAL_OBJECT;
    while (fgets(line_buffer, BUFFER_LENGTH, f) != NULL)
    {
        if (line_buffer[strlen(line_buffer)-1] != '\n')
        {
            fprintf(stderr, "Error in read_planet_file: Line buffer to short\n");
            goto read_planet_file_error;
        }

        /* skip comments beginning with # */
        if (line_buffer[0] == '#')
        {
            continue;
        }

        #if PLANET_FORMAT == CART
        if (sscanf(line_buffer, "%lf %lf %lf %lf %lf %lf %lf %d %lf", 
                   &(g_nb.m[i]), &(g_nb.x[i]), &(g_nb.y[i]), &(g_nb.z[i]), 
                   &(g_nb.vx[i]), &(g_nb.vy[i]), &(g_nb.vz[i]), 
                   &(g_nb.feelsDisk[i]), &(g_nb.rampupTime[i])) != 9)
        {
            fprintf(stderr, "Error in read_planet_file: Not enough values in line %d\n", i);
            goto read_planet_file_error;
        }
        #elif PLANET_FORMAT == ORBIT
        double a, e, inc, Omega, omega, f;
        sscanf(line_buffer, "%lf %lf %lf %lf %lf %lf %lf %d %lf", 
               &(g_nb.m[i]), &a, &e, &inc, &Omega, &omega, 
               &f, &(g_nb.feelsDisk[i]), &(g_nb.rampupTime[i]));
        //{
        //    fprintf(stderr, "Error in read_planet_file: Not enough values in line %d\n", i);
	    fprintf(stderr, "%lf %lf %lf %lf %lf %lf %lf %d %lf\n", g_nb.m[i], a, e, inc, Omega, omega, 
                f, g_nb.feelsDisk[i], g_nb.rampupTime[i]);
        //    goto read_planet_file_error;
        //}
        double cO = cos(Omega);
        double sO = sin(Omega);
        double co = cos(omega);
        double so = sin(omega);
        double cf = cos(f);
        double sf = sin(f);
        double ci = cos(inc);
        double si = sin(inc);

        double mu = CONST_G_CODE_UNITS * (nbodyMass(i) + g_nb.m[i]);
        nbodyCenterOfMassCoordinatesAndVelocity(i);

        double r = a*(1-e*e)/(1+e*cos(f));
        double v0 = sqrt(mu/(a*(1.0-e*e))); 

        g_nb.x[i] = g_nb.comX + r*(cO*(co*cf-so*sf) - sO*(so*cf+co*sf)*ci);
        g_nb.y[i] = g_nb.comY + r*(sO*(co*cf-so*sf) + cO*(so*cf+co*sf)*ci);
        g_nb.z[i] = g_nb.comZ + r*(so*cf+co*sf)*si;
        
        g_nb.vx[i] = g_nb.comVX + v0*((e+cf)*(-ci*co*sO - cO*so) - sf*(co*cO - ci*so*sO));
        g_nb.vy[i] = g_nb.comVY + v0*((e+cf)*(ci*co*cO - sO*so) - sf*(co*sO + ci*so*cO));
        g_nb.vz[i] = g_nb.comVZ + v0*((e+cf)*co*si - sf*si*so);
        #endif
        i++;
    }
    fclose(f);
    free(line_buffer);
    return i;

read_planet_file_error:
    free(line_buffer);
    fclose(f);
    return -1;
}

void nbodyWriteCoordinates()
{
    if (prank == 0)
    {
        char fname[512];
        static double tpos = -1.0;
        int dummy;
        FILE *fp;

        sprintf(fname, "%s/nbody_coordinates.dat", RuntimeGet()->output_dir);

        /* Write header */
        if (g_stepNumber == 0)
        {
            fp = fopen(fname, "w");
            fprintf(fp, "# Particle id\n");
            fprintf(fp, "# Time in code units\n");
            fprintf(fp, "# x in code units\n");
            fprintf(fp, "# y in code units\n");
            fprintf(fp, "# z in code units\n");
            fprintf(fp, "# vx in code units\n");
            fprintf(fp, "# vy in code units\n");
            fprintf(fp, "# vz in code units\n");
            fprintf(fp, "#%-2s   %-12s   %-12s   %-12s   %-12s   %-12s   %-12s   %-12s\n",
                    "id", "t", "x", "y", "z", "vx", "vy", "vz");
        }
        else
        {
            if (tpos < 0.0)
            {
                char sline[512];
                fp = fopen(fname, "r");
                while (fgets(sline, 512, fp))
                    ;
                sscanf(sline, "%d %lf\n", &dummy, &tpos);
                fclose(fp);
            }
            fp = fopen(fname, "a");
        }

        if (g_time > tpos)
        {
            int l;
            for (l = 0; l < NB_N; l++)
            {
                fprintf(fp, "%-3d  %- 12.6e  %- 12.6e  %- 12.6e  %- 12.6e  %- 12.6e  %- 12.6e  %- 12.6e\n",
                        l,
                        g_time,
                        g_nb.x[l], 
                        g_nb.y[l], 
                        g_nb.z[l], 
                        g_nb.vx[l], 
                        g_nb.vy[l], 
                        g_nb.vz[l]);
            }
        }
        fclose(fp);
    }
}

void nbodyWriteOrbitalElements()
{
    if (prank == 0)
    {
        char fname[512];
        static double tpos = -1.0;
        int dummy;
        FILE *fp;

        sprintf(fname, "%s/nbody_orbital_elements.dat", RuntimeGet()->output_dir);

        /* Write header */
        if (g_stepNumber == 0)
        {
            fp = fopen(fname, "w");
            fprintf(fp, "# id: Particle id\n");
            fprintf(fp, "# t: time in code units\n");
            fprintf(fp, "# a: semi-major axis in code units\n");
            fprintf(fp, "# e: eccentricity \n");
            fprintf(fp, "# i: inclination in rad [0, pi]\n");
            fprintf(fp, "# Omega: longitude of ascending node in rad [0, 2pi]\n");
            fprintf(fp, "# omega: argument of pericentre in rad [0, 2pi]\n");
            fprintf(fp, "# f: true anomaly in rad [0, 2pi]\n");
            fprintf(fp, "# P: period in code units\n");
            fprintf(fp, "# E: eccentric anomaly in rad [0, 2pi]\n");
            fprintf(fp, "# M: mean anomaly in rad [0, 2pi]\n");
            fprintf(fp, "#%-2s   %-12s   %-12s   %-12s   %-12s   %-12s   %-12s   %-12s   %-12s   %-12s   %-12s\n",
                    "id", "t", "a", "e", "i", "Omega", "omega", "f", "P", "E", "M");
        }
        else
        {
            if (tpos < 0.0)
            {
                char sline[512];
                fp = fopen(fname, "r");
                while (fgets(sline, 512, fp))
                    ;
                sscanf(sline, "%d %lf\n", &dummy, &tpos);
                fclose(fp);
            }
            fp = fopen(fname, "a");
        }

        if (g_time > tpos)
        {
            nbodyCalcOrbitalElements();
            int l;
            for (l = 1; l < NB_N; l++)
            {
                fprintf(fp, "%-3d  %- 12.6e  %- 12.6e  %- 12.6e  %- 12.6e  %- 12.6e  %- 12.6e  %- 12.6e  %- 12.6e  %- 12.6e  %- 12.6e\n",
                        l,
                        g_time,
                        g_nb.a[l], 
                        g_nb.e[l], 
                        g_nb.inc[l], 
                        g_nb.Omega[l], 
                        g_nb.omega[l], 
                        g_nb.f[l],
                        g_nb.P[l],
                        g_nb.ea[l],
                        g_nb.ma[l]);
            }
        }
        fclose(fp);
    }
}

void nbodyWriteCOM()
{
    if (prank == 0)
    {
        char fname[512];
        static double tpos = -1.0;
        int dummy;
        FILE *fp;

        sprintf(fname, "%s/nbody_com.dat", RuntimeGet()->output_dir);

        /* Write header */
        if (g_stepNumber == 0)
        {
            fp = fopen(fname, "w");
            fprintf(fp, "# id = 0: COM of central object\n");
            fprintf(fp, "# id = 1: COM of all objects\n");
            fprintf(fp, "# Time in code units\n");
            fprintf(fp, "# comX in code units\n");
            fprintf(fp, "# comY in code units\n");
            fprintf(fp, "# comZ in code units\n");
            fprintf(fp, "# comVelX in code units\n");
            fprintf(fp, "# comVelY in code units\n");
            fprintf(fp, "# comVelZ in code units\n");
            fprintf(fp, "#%-2s   %-12s   %-12s   %-12s   %-12s   %-12s   %-12s   %-12s\n",
                    "id", "t", "x", "y", "z", "vx", "vy", "vz");
        }
        else
        {
            if (tpos < 0.0)
            {
                char sline[512];
                fp = fopen(fname, "r");
                while (fgets(sline, 512, fp))
                    ;
                sscanf(sline, "%d %lf\n", &dummy, &tpos);
                fclose(fp);
            }
            fp = fopen(fname, "a");
        }

        if (g_time > tpos)
        {
            /* Center of mass of central object (id = 0) */
            nbodyCenterOfMassCoordinatesAndVelocity(CENTRAL_OBJECT);
            fprintf(fp, "%-3d  %- 12.6e  %- 12.6e  %- 12.6e  %- 12.6e  %- 12.6e  %- 12.6e  %- 12.6e\n",
                    0, g_time, g_nb.comX, g_nb.comY, g_nb.comZ, g_nb.comVX, g_nb.comVY, g_nb.comVZ);

            /* Center of mass of all bodies (id = 1) */
            nbodyCenterOfMassCoordinatesAndVelocity(NB_N);
            fprintf(fp, "%-3d  %- 12.6e  %- 12.6e  %- 12.6e  %- 12.6e  %- 12.6e  %- 12.6e  %- 12.6e\n",
                    1, g_time, g_nb.comX, g_nb.comY, g_nb.comZ, g_nb.comVX, g_nb.comVY, g_nb.comVZ);
        }
        fclose(fp);
    }
}

void nbodyWriteRestartFile(Output *output)
{
    FILE *fout;
    char filename[512], sline[512];
    sprintf(filename, "%s/nbody.out", output->dir);

    if (prank == 0)
    {
        if (output->nfile == 0)
        {
            fout = fopen(filename, "w");
        }
        else
        {
            fout = fopen(filename, "r+");
            int nv,l;
            for (nv = 0; nv < output->nfile; nv++) 
            {
                for (l = 0; l < NB_N; l++)
                {
                    fgets(sline, 512, fout);
                }
            }
            fseek(fout, ftell(fout), SEEK_SET);
        }

        int l;
        for (l = 0; l < NB_N; l++)
        {
            fprintf(fout, "%d %d % 22.15e % 22.15e % 22.15e % 22.15e % 22.15e % 22.15e % 22.15e %d\n",
                    output->nfile, l, g_nb.m[l], g_nb.x[l], g_nb.y[l], g_nb.z[l],
                    g_nb.vx[l], g_nb.vy[l], g_nb.vz[l], g_nb.feelsDisk[l]);
        }
        fclose(fout);
    }
}

void nbodyWriteOrbitalElementsAtOutputTime(Output *output)
{
    FILE *fout;
    char filename[512], sline[512];
    sprintf(filename, "%s/nbody_orbital_elements.out", output->dir);

    if (prank == 0)
    {
        if (output->nfile == 0)
        {
            fout = fopen(filename, "w");
        }
        else
        {
            fout = fopen(filename, "r+");
            int nv,l;
            for (nv = 0; nv < output->nfile; nv++) 
            {
                for (l = 0; l < NB_N; l++)
                {
                    fgets(sline, 512, fout);
                }
            }
            fseek(fout, ftell(fout), SEEK_SET);
        }

        nbodyCalcOrbitalElements();
        int l;
        for (l = 1; l < NB_N; l++)
        {
            fprintf(fout, "%d %d % 22.15e % 22.15e % 22.15e % 22.15e % 22.15e % 22.15e % 22.15e % 22.15e % 22.15e\n",
                    output->nfile, l, g_nb.a[l], g_nb.e[l], g_nb.inc[l], g_nb.Omega[l],
                    g_nb.omega[l], g_nb.f[l], g_nb.P[l], g_nb.ea[l], g_nb.ma[l]);
        }
        fclose(fout);
    }
}

void nbodyReadRestartFile(Runtime *ini, int nrestart)
{
    char fname[512], sline[512];
    FILE *f;

    if (prank == 0)
    {
        sprintf(fname,"%s/nbody.out", ini->output_dir);
        f = fopen(fname, "r");

        int nv,l;
        for (nv = 0; nv < nrestart; nv++)
        {
            for (l = 0; l < NB_N; l++)
            {
                fgets(sline, 512, f);
            }
        }

        for (l = 0; l < NB_N; l++)
        {
            fscanf(f, "%*d %*d %lf %lf %lf %lf %lf %lf %lf %*d",
                   &g_nb.m[l], &g_nb.x[l], &g_nb.y[l], &g_nb.z[l],
                   &g_nb.vx[l], &g_nb.vy[l], &g_nb.vz[l]);
        }

        fclose(f);
    }
    
    #ifdef PARALLEL
    MPI_Bcast(g_nb.m, NB_N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(g_nb.x, NB_N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(g_nb.y, NB_N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(g_nb.z, NB_N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(g_nb.vx, NB_N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(g_nb.vy, NB_N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(g_nb.vz, NB_N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    #endif
}

void nbodyPrintCoordinates()
{
    int l;
    for (l = 0; l < NB_N; l++)
    {
        printf("%d %d % 12.6e % 12.6e % 12.6e % 12.6e % 12.6e % 12.6e % 12.6e\n",
                        prank, l, g_nb.m[l], g_nb.x[l], g_nb.y[l], g_nb.z[l],
                        g_nb.vx[l], g_nb.vy[l], g_nb.vz[l]);
    }
}

