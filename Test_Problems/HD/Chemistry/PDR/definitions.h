#define  PHYSICS                        HD
#define  DIMENSIONS                     1
#define  GEOMETRY                       CARTESIAN
#define  BODY_FORCE                     NO
#define  COOLING                        NO
#define  RECONSTRUCTION                 WENO3
#define  TIME_STEPPING                  RK3
#define  PARTICLES                      NO
#define  USER_DEF_PARAMETERS            3

/* -- physics dependent declarations -- */

#define  DUST_FLUID                     NO
#define  EOS                            IDEAL
#define  ENTROPY_SWITCH                 NO
#define  THERMAL_CONDUCTION             NO
#define  VISCOSITY                      NO
#define  ROTATING_FRAME                 NO
#define  CHEMISTRY                      YES

/* -- user-defined parameters (labels) -- */

#define  RHO                            0
#define  TEMP                           1
#define  GAMMA_EOS                      2

/* [Beg] user-defined constants (do not change this line) */

#define  INTERNAL_BOUNDARY              YES
#define  LIMITER                        MC_LIM
#define  UNIT_DENSITY                   (1.e3*CONST_mp)
#define  UNIT_LENGTH                    1.
#define  UNIT_VELOCITY                  1.

/* [End] user-defined constants (do not change this line) */
