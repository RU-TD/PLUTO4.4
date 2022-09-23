#define  PHYSICS                        HD
#define  DIMENSIONS                     2
#define  GEOMETRY                       SPHERICAL
#define  BODY_FORCE                     POTENTIAL
#define  COOLING                        NO
#define  RECONSTRUCTION                 LINEAR
#define  TIME_STEPPING                  RK2
#define  PARTICLES                      NO
#define  USER_DEF_PARAMETERS            5

/* -- physics dependent declarations -- */

#define  DUST_FLUID                     NO
#define  EOS                            IDEAL
#define  ENTROPY_SWITCH                 NO
#define  THERMAL_CONDUCTION             NO
#define  VISCOSITY                      EXPLICIT
#define  ROTATING_FRAME                 NO
#define  CHEMISTRY                      YES
#define  DISABLE_HYDRO                  YES

/* -- user-defined parameters (labels) -- */

#define  MSTAR                          0
#define  G_LX                           1
#define  GAMMA_EOS                      2
#define  MU                             3
#define  TEMP                           4

/* [Beg] user-defined constants (do not change this line) */

#define  INTERNAL_BOUNDARY              YES
#define  LIMITER                        MINMOD_LIM
#define  SHOCK_FLATTENING               MULTID
#define  FAILSAFE                       TRUE
#define  UNIT_LENGTH                    (10.0*CONST_au)
#define  UNIT_DENSITY                   (CONST_Msun/(UNIT_LENGTH*UNIT_LENGTH*UNIT_LENGTH))
#define  UNIT_VELOCITY                  (sqrt(CONST_G*CONST_Msun/UNIT_LENGTH)/(2.*CONST_PI))

/* [End] user-defined constants (do not change this line) */
