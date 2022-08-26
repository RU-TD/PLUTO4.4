/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Set labels, indexes and prototypes for the HD module.

  Contains variable names and prototypes for the HD module

  \author A. Mignone (mignone@to.infn.it)
  \date   Dec 2, 2020
*/
/* ///////////////////////////////////////////////////////////////////// */

/* *********************************************************
    Set flow variable indices.
    Extra vector components, when not needed, point to the
    last element (255) of the array stored by startup.c.  
   ********************************************************* */

#define  RHO 0
#define  MX1 1
#define  MX2 2
#define  MX3 3
#if HAVE_ENERGY
  #define ENG  4
  #define PRS  ENG
#endif

#define VX1   MX1
#define VX2   MX2
#define VX3   MX3

#define NFLX (4 + HAVE_ENERGY)

#define NTRACER 40
#define NPHOTO 1000
#define NIONS  0

#define IDX_CHEM_C (NFLX + NIONS + 0)
#define IDX_CHEM_CH (NFLX + NIONS + 1)
#define IDX_CHEM_CH2 (NFLX + NIONS + 2)
#define IDX_CHEM_CH2j (NFLX + NIONS + 3)
#define IDX_CHEM_CH3 (NFLX + NIONS + 4)
#define IDX_CHEM_CH3j (NFLX + NIONS + 5)
#define IDX_CHEM_CH4 (NFLX + NIONS + 6)
#define IDX_CHEM_CH4j (NFLX + NIONS + 7)
#define IDX_CHEM_CH5j (NFLX + NIONS + 8)
#define IDX_CHEM_CHj (NFLX + NIONS + 9)
#define IDX_CHEM_CO (NFLX + NIONS + 10)
#define IDX_CHEM_CO_DUST (NFLX + NIONS + 11)
#define IDX_CHEM_COj (NFLX + NIONS + 12)
#define IDX_CHEM_Cj (NFLX + NIONS + 13)
#define IDX_CHEM_Cjj (NFLX + NIONS + 14)
#define IDX_CHEM_Cjjj (NFLX + NIONS + 15)
#define IDX_CHEM_Cjjjj (NFLX + NIONS + 16)
#define IDX_CHEM_E (NFLX + NIONS + 17)
#define IDX_CHEM_H (NFLX + NIONS + 18)
#define IDX_CHEM_H2 (NFLX + NIONS + 19)
#define IDX_CHEM_H2O (NFLX + NIONS + 20)
#define IDX_CHEM_H2O_DUST (NFLX + NIONS + 21)
#define IDX_CHEM_H2Oj (NFLX + NIONS + 22)
#define IDX_CHEM_H2j (NFLX + NIONS + 23)
#define IDX_CHEM_H3Oj (NFLX + NIONS + 24)
#define IDX_CHEM_H3j (NFLX + NIONS + 25)
#define IDX_CHEM_HCOj (NFLX + NIONS + 26)
#define IDX_CHEM_He (NFLX + NIONS + 27)
#define IDX_CHEM_Hej (NFLX + NIONS + 28)
#define IDX_CHEM_Hejj (NFLX + NIONS + 29)
#define IDX_CHEM_Hj (NFLX + NIONS + 30)
#define IDX_CHEM_O (NFLX + NIONS + 31)
#define IDX_CHEM_O2 (NFLX + NIONS + 32)
#define IDX_CHEM_O2j (NFLX + NIONS + 33)
#define IDX_CHEM_OH (NFLX + NIONS + 34)
#define IDX_CHEM_OHj (NFLX + NIONS + 35)
#define IDX_CHEM_Oj (NFLX + NIONS + 36)
#define IDX_CHEM_Ojj (NFLX + NIONS + 37)
#define IDX_CHEM_Ojjj (NFLX + NIONS + 38)
#define IDX_CHEM_Ojjjj (NFLX + NIONS + 39)

/* *************************************************
     Now define more convenient and user-friendly 
     pointer labels for geometry setting      
   ************************************************* */

#if GEOMETRY == CYLINDRICAL 
  #define iVR    VX1
  #define iMR    MX1

  #define iVZ    VX2
  #define iMZ    MX2

  #define iVPHI  VX3
  #define iMPHI  MX3
#endif

#if GEOMETRY == POLAR 
  #define iVR    VX1
  #define iMR    MX1

  #define iVPHI  VX2
  #define iMPHI  MX2

  #define iVZ    VX3
  #define iMZ    MX3
#endif

#if GEOMETRY == SPHERICAL 
  #define iVR    VX1
  #define iMR    MX1

  #define iVTH   VX2
  #define iMTH   MX2

  #define iVPHI  VX3
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

int  ConsToPrim   (double **, double **, int, int, uint16_t *);
void Eigenvalues (double **, double *, double **, int, int);
void PrimEigenvectors (const State *, int, int);
void ConsEigenvectors (double *, double *, double, 
                       double **, double **, double *);

void Flux      (const State *, int, int);
void HLL_Speed (const State *, const State *, double *, double *, int, int);
void MaxSignalSpeed (const State *, double *, double *, int, int);
void PrimToCons   (double **, double **, int, int);
void PrimRHS    (double *, double *, double, double, double *);
void PrimSource (const State *, double **, int, int, Grid *);

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
        double**** column_density;
        double**** jflux;
	double*    jflux0;
        double*    data_buffer;
	double*    jflux_buffer;
        double*    column_density_offset;

} IrradiationData;

void initialize_ColumnDensity(Grid *grid);

Riemann_Solver TwoShock_Solver, LF_Solver, Roe_Solver, HLL_Solver,
               HLLC_Solver, RusanovDW_Solver;
Riemann_Solver AUSMp_Solver;

IrradiationData irradiation;

