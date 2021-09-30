#define  PHYSICS                        HD
#define  DIMENSIONS                     3
#define  COMPONENTS                     3
#define  GEOMETRY                       SPHERICAL
#define  BODY_FORCE                     POTENTIAL
#define  COOLING                        NO
#define  RECONSTRUCTION                 LINEAR
#define  TIME_STEPPING                  RK2
#define  NTRACER                        0
#define  PARTICLES                      NO
#define  USER_DEF_PARAMETERS            13
#define  FARGO_OUTPUT_VTOT              YES
#define  FARGO_NSTEP_AVERAGE            -1

/* -- physics dependent declarations -- */

#define  DUST_FLUID                     NO
#define  EOS                            IDEAL
#define  ENTROPY_SWITCH                 NO
#define  THERMAL_CONDUCTION             NO
#define  VISCOSITY                      RK_LEGENDRE
#define  ROTATING_FRAME                 YES

/* -- user-defined parameters (labels) -- */

// #define  ALPHA                          0
#define  CONST_VISCOSITY                0
#define  GAMMA                          1
#define  MU                             2
#define  LX                             3
#define  MS                             4
#define  MP                             5
#define  RSM                            6
#define  PLANETRAMPUPTIME               7
#define  PLANETRAMPUPM0                 8
#define  TAU_ACC                        9
#define  TAU_DAMP                       10
#define  DAMPING_RIN                    11
#define  DAMPING_ROUT                   12

/* [Beg] user-defined constants (do not change this line) */

#define  WARNING_MESSAGES               NO
#define  INTERNAL_BOUNDARY              YES
#define  SHOCK_FLATTENING               MULTID
#define  CHAR_LIMITING                  YES
#define  LIMITER                        MINMOD_LIM
#define  UNIT_LENGTH                    (5.2*CONST_au)
#define  UNIT_DENSITY                   (CONST_Msun/(UNIT_LENGTH*UNIT_LENGTH*UNIT_LENGTH))
#define  UNIT_VELOCITY                  (sqrt(CONST_G*g_inputParam[MS]*CONST_Msun/UNIT_LENGTH)/(2.*CONST_PI))

/* [End] user-defined constants (do not change this line) */

#define INPUT_GRID_FILE                 "input/grid.0500.out"
#define INPUT_DATA_FILE                 "input/data.0500.dbl"
#define PLANET_HEATING                  SML16
#define REMOVE_ACCRETING_MASS           NO
#define INCLUDE_PARTICLES               NO // Giovanni Picogna's particle implementation
#define PARTICLES_TURBULENCE            NO
// #define DISABLE_HYDRODYNAMICS           YES