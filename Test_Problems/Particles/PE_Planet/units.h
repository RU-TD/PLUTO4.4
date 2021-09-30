#ifndef MLW_H_UNITS
#define MLW_H_UNITS

#define UNIT_MASS (UNIT_DENSITY * UNIT_LENGTH * UNIT_LENGTH * UNIT_LENGTH)
#define UNIT_G (CONST_G / POW2(UNIT_VELOCITY) / UNIT_LENGTH * UNIT_MASS) // gravitational constant in code units

#endif