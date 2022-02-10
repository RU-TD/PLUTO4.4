/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Set labels, indexes and prototypes for the MHD module.

  Contains basic macro definitions, structure definitions and global
  variable declarations used by the MHD module.

  \author A. Mignone (mignone@ph.unito.it)
  \date   April, 2, 2015
*/
/* ///////////////////////////////////////////////////////////////////// */

/* *********************************************************
    Set flow variable indices.
    Extra vector components, when not needed, point to the
    last element (255) of the array stored by startup.c.  
   ********************************************************* */

#define  RHO 0
#define  MX1 1
#define  MX2 (COMPONENTS >= 2 ? 2: 255)
#define  MX3 (COMPONENTS == 3 ? 3: 255) 
#if HAVE_ENERGY
  #define ENG  (2*COMPONENTS + 1)
  #define PRS  ENG
#endif

#define VX1   MX1
#define VX2   MX2
#define VX3   MX3

#define NFLX (1 + 2*COMPONENTS + HAVE_ENERGY)

#define NTRACER 31
#define NPHOTO 1000
#define NIONS  0

#define IDX_CHEM_H (NFLX + NIONS + 0)
#define IDX_CHEM_CH (NFLX + NIONS + 1)
#define IDX_CHEM_C (NFLX + NIONS + 2)
#define IDX_CHEM_H2 (NFLX + NIONS + 3)
#define IDX_CHEM_CH3 (NFLX + NIONS + 4)
#define IDX_CHEM_CH2 (NFLX + NIONS + 5)
#define IDX_CHEM_CH4 (NFLX + NIONS + 6)
#define IDX_CHEM_OH (NFLX + NIONS + 7)
#define IDX_CHEM_O (NFLX + NIONS + 8)
#define IDX_CHEM_H2O (NFLX + NIONS + 9)
#define IDX_CHEM_CO (NFLX + NIONS + 10)
#define IDX_CHEM_O2 (NFLX + NIONS + 11)
#define IDX_CHEM_CH2j (NFLX + NIONS + 12)
#define IDX_CHEM_CHj (NFLX + NIONS + 13)
#define IDX_CHEM_CH3j (NFLX + NIONS + 14)
#define IDX_CHEM_Hej (NFLX + NIONS + 15)
#define IDX_CHEM_He (NFLX + NIONS + 16)
#define IDX_CHEM_Hj (NFLX + NIONS + 17)
#define IDX_CHEM_Cj (NFLX + NIONS + 18)
#define IDX_CHEM_Oj (NFLX + NIONS + 19)
#define IDX_CHEM_H2j (NFLX + NIONS + 20)
#define IDX_CHEM_COj (NFLX + NIONS + 21)
#define IDX_CHEM_E (NFLX + NIONS + 22)
#define IDX_CHEM_H3j (NFLX + NIONS + 23)
#define IDX_CHEM_CH4j (NFLX + NIONS + 24)
#define IDX_CHEM_OHj (NFLX + NIONS + 25)
#define IDX_CHEM_CH5j (NFLX + NIONS + 26)
#define IDX_CHEM_H2Oj (NFLX + NIONS + 27)
#define IDX_CHEM_H3Oj (NFLX + NIONS + 28)
#define IDX_CHEM_HCOj (NFLX + NIONS + 29)
#define IDX_CHEM_O2j (NFLX + NIONS + 30)

/* *************************************************
     Now define more convenient and user-friendly 
     pointer labels for geometry setting      
   ************************************************* */

#if GEOMETRY == CARTESIAN

  #define VX    VX1
  #define VY    VX2
  #define VZ    VX3

  #define MX    MX1
  #define MY    MX2
  #define MZ    MX3

#endif

#if GEOMETRY == CYLINDRICAL 

 #define iVR    VX1
 #define iVZ    VX2
 #define iVPHI  VX3

 #define iMR    MX1
 #define iMZ    MX2
 #define iMPHI  MX3

#endif

#if GEOMETRY == POLAR 
 #define iVR    VX1
 #define iVPHI  VX2
 #define iVZ    VX3

 #define iMR    MX1
 #define iMPHI  MX2
 #define iMZ    MX3
#endif

#if GEOMETRY == SPHERICAL 
 #define iVR     VX1
 #define iVTH    VX2
 #define iVPHI   VX3

 #define iMR    MX1
 #define iMTH   MX2
 #define iMPHI  MX3
#endif

 /* *************************************************
      Label the different waves in increasing order 
      following the number of vector components.
    ************************************************* */
 
 enum KWAVES {
  KSOUNDM, KSOUNDP
  #if HAVE_ENERGY
   , KENTRP
  #endif
 };

 /* ***********************************************************
                    Prototyping goes here          
    *********************************************************** */
 
 int  ConsToPrim   (double **, double **, int, int, unsigned char *);
 void Eigenvalues (double **, double *, double **, int, int);
 void PrimEigenvectors (double *, double, double, double *, double **, double **);
 void ConsEigenvectors (double *, double *, double,
                        double **, double **, double *);
 
 void Flux      (double **, double **, double *, double **, double *, int, int);
 void HLL_Speed (double **, double **, double *, double *,
                 double *, double *, int, int);
 void MaxSignalSpeed (double **, double *, double *, double *, int, int);
 void PrimToCons   (double **, double **, int, int);
 void PrimRHS    (double *, double *, double, double, double *);
 void PrimSource (const State_1D *, int, int,
                  double *, double *, double **, Grid *);
 
 Riemann_Solver TwoShock_Solver, LF_Solver, Roe_Solver, HLL_Solver,
                HLLC_Solver, RusanovDW_Solver;
 Riemann_Solver AUSMp_Solver;


typedef struct LOCAL_DOMAIN_INFO
{
        int x1_begin;
        int x1_end;

        int x2_begin;
        int x2_end;

        int x3_begin;
        int x3_end;

} LocalDomainInfo;

typedef struct COMMUNICATION_NEIGHBOUR
{
        int receive_rank;
        int send_rank;
} CommunicationNeighbour;

typedef struct IRRADIATION_DATA
{
        CommunicationNeighbour neighbour;
        double****      column_density;
	double****      jflux;
        double* data_buffer;
	double* jflux_buffer;
        double* column_density_offset;
	double* jflux_offset;

} IrradiationData;

void initialize_ColumnDensity(Grid *grid);

void find_CommunicationNeighbours(Grid *grid);
void calculate_ColumnDensity(Grid *grid, const Data* data);
void calculate_JFlux(Grid *grid, const Data* data);

void find_CommunicationNeighbour(int current_rank, LocalDomainInfo *domain_info_array, int nproc, CommunicationNeighbour* cn);
void calculate_ColumnDensity_perDomain(Grid* grid, const Data* data, int val);
void calculate_JFlux_perDomain(Grid* grid, const Data* data, int val);

Riemann_Solver TwoShock_Solver, LF_Solver, Roe_Solver, HLL_Solver,
               HLLC_Solver, RusanovDW_Solver;
Riemann_Solver AUSMp_Solver;

IrradiationData irradiation;
