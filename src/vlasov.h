// -*- C++ -*-
// $Id$

#include "index.h"
#include "matrix.h"

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
 *
 *  
 */


#ifndef VLASOV_H
#define VLASOV_H

#define PI 3.1415926535897932384626433832795029

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

/** @page indices Array indices
 *
 *  @section general_ind General
 *
 *  We divide the whole simulation area into a number of grid cells. The
 *  interior of the simulation is indexed from \f$1\ldots N\f$. If we have 
 *  periodic boundary conditions, this means that \f$A(0) = A(N)\f$, of in 
 *  general \f$A(i) = A(N+i)\f$.
 *
 *  The physical extent of the simulation box is \f$N\Delta x\f$.
 *
 *  @section disttrib_ind Distribution Function
 *
 *  The distribution function may be distributed among a number
 *  of processors. For the third order advancer we need at least two
 *  cells to every side of the. The whole simulation region in memory thus 
 *  has the extent \f$ -1 \ldots N+2 \f$.
 *  
 *  @section fields_ind Fields and Potentials solver
 *
 *  Internally the Poisson solver uses the Fortran hwscrt routine. 
 *  In the case of periodic boundaries this routine needs one cell to be 
 *  duplicated. This means the region spans \f$1\ldots N+1\f$, where
 *  \f$A(1) = A(N+1)\f$. If the boundaries are fixed, however we still want the 
 *  interior to span \f$1\ldots N\f$, so we need an array \f$0\ldots N+1\f$.
 *  In general potentials will thus span \f$0\ldots N+1\f$ and in the case of
 *  periodic boundary conditions only the sub-array \f$1\ldots N+1\f$ will
 *  be passed to hwscrt. 
 *
 *  The bondary points can be exploited when spatial derivatives of the 
 *  potentials have to be taken.
 *
 */
 
/** @mainpage 2d3v Vlasov Darwin Code
 *
 *  @section faq Questions and Answers
 *
 *  For the use of indices in the arrays see the page @ref indices
 */
 
