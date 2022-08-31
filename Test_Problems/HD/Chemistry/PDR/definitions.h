#define  PHYSICS                        HD
#define  DIMENSIONS                     1
#define  GEOMETRY                       CARTESIAN
#define  BODY_FORCE                     NO
#define  COOLING                        NO
#define  RECONSTRUCTION                 WENO3
#define  TIME_STEPPING                  RK3
#define  PARTICLES                      NO
#define  USER_DEF_PARAMETERS            4

/* -- physics dependent declarations -- */

#define  DUST_FLUID                     NO
#define  EOS                            IDEAL
#define  ENTROPY_SWITCH                 NO
#define  THERMAL_CONDUCTION             NO
#define  VISCOSITY                      NO
#define  ROTATING_FRAME                 NO
#define  CHEMISTRY                      YES
#define  DISABLE_HYDRO                  YES

/* -- user-defined parameters (labels) -- */

#define  N_GAS                          0
#define  MU                             1
#define  TEMP                           2
#define  GAMMA_EOS                      3

/* [Beg] user-defined constants (do not change this line) */

#define  INTERNAL_BOUNDARY              YES
#define  LIMITER                        MC_LIM
#define  SHOCK FLATTENING               ONED
#define  FAILSAFE                       TRUE
#define  UNIT_DENSITY                   (1.e5*CONST_mp)
#define  UNIT_LENGTH                    1.
#define  UNIT_VELOCITY                  (1./3.154e7)
#define  MULTIPLE_LOG_FILES             NO

/* [End] user-defined constants (do not change this line) */
