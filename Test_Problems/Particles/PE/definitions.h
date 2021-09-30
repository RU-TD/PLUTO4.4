#define  PHYSICS                        HD
#define  DIMENSIONS                     2
#define  COMPONENTS                     3
#define  GEOMETRY                       SPHERICAL
#define  BODY_FORCE                     POTENTIAL
#define  COOLING                        NO
#define  RECONSTRUCTION                 LINEAR
#define  TIME_STEPPING                  RK2
#define  NTRACER                        0
#define  PARTICLES                      NO
#define  USER_DEF_PARAMETERS            6

/* -- physics dependent declarations -- */

#define  DUST_FLUID                     NO
#define  EOS                            IDEAL
#define  ENTROPY_SWITCH                 NO
#define  THERMAL_CONDUCTION             NO
#define  VISCOSITY                      RK_LEGENDRE
#define  ROTATING_FRAME                 NO

/* -- user-defined parameters (labels) -- */

#define  ALPHA                          0
#define  GAMMA                          1
#define  MU                             2
#define  LX                             3
#define  MS                             4
#define  DAMPING_ROUT                   5

/* [Beg] user-defined constants (do not change this line) */

#define  INTERNAL_BOUNDARY              YES
#define  INITIAL_SMOOTHING              NO
#define  LIMITER                        MINMOD_LIM
#define  SHOCK_FLATTENING               MULTID
#define  CHAR_LIMITING                  YES
#define  UNIT_LENGTH                    (5.2*CONST_au)
#define  UNIT_DENSITY                   (CONST_Msun/(UNIT_LENGTH*UNIT_LENGTH*UNIT_LENGTH))
#define  UNIT_VELOCITY                  (sqrt(CONST_G*g_inputParam[MS]*CONST_Msun/UNIT_LENGTH)/(2.*CONST_PI))

/* [End] user-defined constants (do not change this line) */
#define INPUT_GRID_FILE                 "input/grid1_jo_converted.out"
#define INPUT_DATA_FILE                 "input/data1_jo_converted.dbl"
