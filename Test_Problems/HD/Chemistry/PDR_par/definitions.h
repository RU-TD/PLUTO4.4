#define  PHYSICS                 HD
#define  DIMENSIONS              1
#define  COMPONENTS              1
#define  GEOMETRY                CARTESIAN
#define  BODY_FORCE              NO
#define  COOLING                 NO
#define  RECONSTRUCTION          WENO3
#define  TIME_STEPPING           RK3
#define  DIMENSIONAL_SPLITTING   YES
#define  USER_DEF_PARAMETERS     3
#define  DISABLE_HYDRODYNAMICS   YES

/* -- physics dependent declarations -- */

#define  EOS                     IDEAL
#define  ENTROPY_SWITCH          NO
#define  THERMAL_CONDUCTION      NO
#define  VISCOSITY               NO
#define  ROTATING_FRAME          NO
#define  CHEMISTRY               YES

/* -- user-defined parameters (labels) -- */

#define  RHO                     0
#define  TEMP                    1
#define  GAMMA_EOS               2

/* [Beg] user-defined constants (do not change this line) */

#define  UNIT_DENSITY            (1.e3*CONST_mp)
#define  UNIT_LENGTH             1.
#define  UNIT_VELOCITY           1.

/* [End] user-defined constants (do not change this line) */

/* -- supplementary constants (user editable) -- */ 

#define  INITIAL_SMOOTHING   NO
#define  WARNING_MESSAGES    YES
#define  PRINT_TO_FILE       NO
#define  INTERNAL_BOUNDARY   YES
#define  SHOCK_FLATTENING    NO
#define  CHAR_LIMITING       NO
#define  LIMITER             MC_LIM
