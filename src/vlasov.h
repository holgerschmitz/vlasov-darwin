#include "opar.h"
#include "index.h"
#include <matrix.H>

/** @file vlasov.h
 *  Defines the classes for solving Vlasov's equation.
 *  The goal is to implement a number of different schemes for advancing 
 *  Vlasov's equation to be able to compare the different advantages and
 *  disadvantages
 *
 *  @todo Test gyro-frequency (how? waves?)
 *  @todo Test magnetosonic modes
 *  @todo Rethink class structure
 *  @todo Update Documentation
 *  @todo Make parallel
 */

#ifndef VLASOV_H
#define VLASOV_H

static const int DIMENSION=2;

typedef FixedArray<int,DIMENSION> PositionI;
typedef FixedArray<double,DIMENSION> PositionD;

typedef FixedArray<int,DIMENSION+3> PhasePositionI;
typedef FixedArray<double,DIMENSION+3> PhasePositionD;

typedef FixedArray<int,3> VelocityI;
typedef FixedArray<double,3> VelocityD;

typedef NumMatrix<double, DIMENSION+3> VlasovDist;
typedef NumMatrix<double, DIMENSION>  ScalarField;



#endif
