/***************************************************************************
                          opar.h  -  description
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
/** @file opar.h
 *  The main header file of the simulation. Edit this file to change
 *  the configuration to one, two or three dimensional.
 *  @todo Parallelization with MPI
 *  @todo A phase space output that will draw the data into an image
 */

#include "Vlas2d.h"

// The dimensionality of the simulation
#define ONE_DIMENSIONAL
//#define TWO_DIMENSIONAL
//#define THREE_DIMENSIONAL

// When this is defined, we are compiling for the next version

#define EL 1.602e-19
#define ELMASS  9.108E-31
#define EPS0 8.8542e-12
#define CLIGHT  3.0e+8

#define PI 3.1415926535897932384626433832795029
