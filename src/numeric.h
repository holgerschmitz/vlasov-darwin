// -*- C++ -*-
// $Id$

 
#include <valarray>
#include <stdlib.h>
 
#ifndef NUMERIC_H
#define NUMERIC_H

#define PIl 3.1415926535897932384626433832795029L
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

/** @brief Returns the error function
 *  @todo Do we really need the error function. Isn't it implemented in math.h
 */
//double erf(double);

template<typename T>
T tanh(T x) {
  T p = exp(x);
  T m = exp(-x);
  return (p-m)/(p+m);
}

template<typename T>
T sinh(T x) {
  return 0.5*(exp(x)-exp(-x));
}

template<typename T>
T cosh(T x) {
  return 0.5*(exp(x)+exp(-x));
}

template<typename T>
T sech(T x) {
  return 2./(exp(x)+exp(-x));
}


#endif /* NUMERIC_H */



