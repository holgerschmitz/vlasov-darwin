// -*- C++ -*-
// $Id$

 
#include <valarray>
#include <stdlib.h>
 
#ifndef NUMERIC_H
#define NUMERIC_H

/** @file numeric.h
 *  Numeric helper functions and classes
 */

/**
 * @brief Random double numbers.
 *
 * This function transforms integer random numbers to double random
 * numbers. For most applications this is ok. However care has to be applied
 * because of the finite graininess of the resulting random numbers.
 *
 * Especially: There is a smallest nonzero frand which is 1/RAND_MAX
 */
inline double frand() { return rand()/double(RAND_MAX); }

double erf(double);



#endif /* NUMERIC_H */



