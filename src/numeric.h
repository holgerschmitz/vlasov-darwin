/***************************************************************************
                          numeric.h  -  description
                             -------------------
    begin                : Fri Jun 2 2000
    copyright            : (C) 2000 by Holger Schmitz
    email                : 
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/*
 *	    OPAR - Open Architecture Particle Simulation
 *
 * Version 0.99
 *
 * (c) 1998 by Th. Daube and H. Schmitz
 */
 
#include <valarray>
#include <stdlib.h>
#include "config.h"
 
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



