// -*- C++ -*-
// $Id$
#ifndef NUMBOUNDARY_H
#define NUMBOUNDARY_H

#include "matrix.h"

  /** @file  numboundary.h
    * @brief Implements the numeric boundary types.
    *
    * Here various boundary types are implemented, for example:
    * Neumann, Dirichlet, symmetrical and periodical conditions.
    */


/// abstract base class of the NumBoundary classes
class NumBoundary {
  public:
    ///purely virtual function
    virtual void apply(NumMatrix<double,2> &u) const = 0;
    ///purely virtual function
    virtual bool normalize() const = 0;
};

/// Symmetric boundary which is build from 2 base types
template<class BoundX, class BoundY>
class SymmetricBoundary : public NumBoundary {
  public:
	///defualt constructor
      SymmetricBoundary() {}
	/// copy contructor
      SymmetricBoundary(const SymmetricBoundary&) {}
	///apply a matrix
      void apply(NumMatrix<double,2> &u) const;
	///call the normalize members of the base boundaries
      bool normalize() const { return BoundX::normalize() && BoundY::normalize(); }
};

/// Boundary which is build from 4 mixed base types 
template<class BoundLeft, class BoundRight, class BoundBottom, class BoundTop>
class MixedBoundary : public NumBoundary {
  public:
	///default constructor
      MixedBoundary() {}
	///copy constructor
      MixedBoundary(const MixedBoundary&) {}
	///apply a matrix
      void apply(NumMatrix<double,2> &u) const;
	///call the normalize members of the base boundaries
      bool normalize() const { 
         return BoundLeft::normalize() && BoundRight::normalize()
          && BoundBottom::normalize() && BoundTop::normalize(); 
      }
};

/// Periodic Boundary
class BoundPeriodic {
  private:
	  const int low0,   ///< lower cell boundary of the processed area
	  		low1,   ///< lower cell boundary of the processed area, shifted by +1
			high0,  ///< upper cell boundary of the processed area
			high1;  ///< upper cell boundary of the processed area, shifted by -1
  public:
	///constructor
      BoundPeriodic(int low, int high) 
        : low0(low), low1(low+1), high0(high), high1(high-1) {}
	///accessor method, return low0 member
      int low() { return low0; }
	///accessor method, return high1 member
      int lowSrc() { return high1; }
	///accessor method, return high0 member
      int high() { return high0; }
	///accessor method, return low1 member
      int highSrc() { return low1; }
	///return 1
      int factor() { return 1; }
	///returns "TRUE", since normalized
      static bool normalize() { return true; }
};

/// Dirichlet Boundary
class BoundDirichlet {
    private:
	    const int low0,  ///< lower cell boundary of the processed area
	              low1,  ///< lower cell boundary of the processed area, shifted by +1
		        high0, ///< upper cell boundary of the processed area
		        high1; ///< upper cell boundary of the processed area, shifted by -1
    public:
	///constructor
      BoundDirichlet(int low, int high) 
        : low0(low), low1(low+1), high0(high), high1(high-1) {}
	///accessor method, return low0 member
      int low() { return low0; }
	///accessor method, return low1 member
      int lowSrc() { return low1; }
	///accessor method, return high0 member
      int high() { return high0; }
	///accessor method, return high1 member
      int highSrc() { return high1; }
	///return -1
      int factor() { return -1; }
	///return "FALSE", since not normalized
      static bool normalize() { return false; }
};

/// Neumann Boundary
class BoundNeumann {
    private:
	    const int low0,    ///< lower cell boundary of the processed area
	              low1,    ///< lower cell boundary of the processed area, shifted by +1
		        high0,   ///< upper cell boundary of the processed area
			  high1;   ///< upper cell boundary of the processed area, shifted by -1
    public:
	///constructor
      BoundNeumann(int low, int high) 
        : low0(low), low1(low+1), high0(high), high1(high-1) {}
	///accessor method, return low0 member
      int low() { return low0; }
	///accessor method, return low1 member
      int lowSrc() { return low1; }
	///accessor method, return high0 member
      int high() { return high0; }
	///accessor method, return high1 member
      int highSrc() { return high1; }
	///return 1
      int factor() { return 1; }
	///return "TRUE", since normalized 
      static bool normalize() { return true; }
};

#include "numboundary.t"

/// a periodic and symmetric boundary
typedef SymmetricBoundary<BoundPeriodic,BoundPeriodic>   PPBoundary;
///mixed symmetric boundary from Dirichlet and periodic type
typedef SymmetricBoundary<BoundPeriodic,BoundDirichlet>  PDBoundary;
///mixed symmetric boundary from Neumann and periodic type
typedef SymmetricBoundary<BoundPeriodic,BoundNeumann>    PNBoundary;

///mixed symmetric boundary from Dirichlet and periodic type
typedef SymmetricBoundary<BoundDirichlet,BoundPeriodic>  DPBoundary;
///symmetric boundary of Dirichlet type
typedef SymmetricBoundary<BoundDirichlet,BoundDirichlet> DDBoundary;
/// symmetric mixed boundary from Dirichlet and Neumann type
typedef SymmetricBoundary<BoundDirichlet,BoundNeumann>   DNBoundary;
///mixed symmetric boundary from Neumann and periodic type
typedef SymmetricBoundary<BoundNeumann,BoundPeriodic>    NPBoundary;
/// symmetric mixed boundary from Dirichlet and Neumann type
typedef SymmetricBoundary<BoundNeumann,BoundDirichlet>   NDBoundary;
///symmetric boundary of Neumann type
typedef SymmetricBoundary<BoundNeumann,BoundNeumann>     NNBoundary;

#endif

