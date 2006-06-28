// -*- C++ -*-
// $Id$
#ifndef NUMBOUNDARY_H
#define NUMBOUNDARY_H

#include "matrix.h"

class NumBoundary {
  public:
    virtual void apply(NumMatrix<double,2> &u) const = 0;
    virtual bool normalize() const = 0;
};

template<class BoundX, class BoundY>
class SymmetricBoundary : public NumBoundary {
  public:
      SymmetricBoundary() {}
      SymmetricBoundary(const SymmetricBoundary&) {}
      void apply(NumMatrix<double,2> &u) const;
      bool normalize() const { return BoundX::normalize() && BoundY::normalize(); }
};

template<class BoundLeft, class BoundRight, class BoundBottom, class BoundTop>
class MixedBoundary : public NumBoundary {
  public:
      MixedBoundary() {}
      MixedBoundary(const MixedBoundary&) {}
      void apply(NumMatrix<double,2> &u) const;
      bool normalize() const { 
         return BoundLeft::normalize() && BoundRight::normalize()
          && BoundBottom::normalize() && BoundTop::normalize(); 
      }
};

template<class BoundLeft, class BoundRight, class BoundBottom, class BoundTop>
class MixedBoundaryWithOffset : public NumBoundary {
  private:
    double ol;
    double or;
    double ob;
    double ot;
    
  public:
    MixedBoundaryWithOffset(double ol_=0.0, double or_=0.0, double ob_=0.0, double ot_=0.0)
     : ol(2.0*ol_), or(2.0*or_), ot(2.0*ot_), ob(2.0*ob_)
    {}
    MixedBoundaryWithOffset(const MixedBoundaryWithOffset& bound)
     : ol(bound.ol), or(bound.or), ob(bound.ob), ot(bound.ot)
    {}
    
    double setOffsetLeft(double ol_)
    {
      ol = 2.0*ol_;
    }
    
    double setOffsetRight(double or_)
    {
      or = 2.0*or_;
    }
    
    double setOffsetBottom(double ob_)
    {
      ob = 2.0*ob_;
    }
    
    double setOffsetTop(double ot_)
    {
      ot = 2.0*ot_;
    }
    
    void apply(NumMatrix<double,2> &u) const;
    bool normalize() const { 
       return false; 
    }
};



class BoundPeriodic {
  private:
      const int low0, low1, high0, high1;
  public:
      BoundPeriodic(int low, int high) 
        : low0(low), low1(low+1), high0(high), high1(high-1) {}
      int low() { return low0; }
      int lowSrc() { return high1; }
      int high() { return high0; }
      int highSrc() { return low1; }
      int factor() { return 1; }
      static bool normalize() { return true; }
};

class BoundDirichlet {
    private:
      const int low0, low1, high0, high1;
    public:
      BoundDirichlet(int low, int high) 
        : low0(low), low1(low+1), high0(high), high1(high-1) {}
      int low() { return low0; }
      int lowSrc() { return low1; }
      int high() { return high0; }
      int highSrc() { return high1; }
      int factor() { return -1; }
      static bool normalize() { return false; }
};

class BoundNeumann {
    private:
      const int low0, low1, high0, high1;
    public:
      BoundNeumann(int low, int high) 
        : low0(low), low1(low+1), high0(high), high1(high-1) {}
      int low() { return low0; }
      int lowSrc() { return low1; }
      int high() { return high0; }
      int highSrc() { return high1; }
      int factor() { return 1; }
      static bool normalize() { return true; }
};

#include "numboundary.t"

typedef SymmetricBoundary<BoundPeriodic,BoundPeriodic>   PPBoundary;
typedef SymmetricBoundary<BoundPeriodic,BoundDirichlet>  PDBoundary;
typedef SymmetricBoundary<BoundPeriodic,BoundNeumann>    PNBoundary;

typedef SymmetricBoundary<BoundDirichlet,BoundPeriodic>  DPBoundary;
typedef SymmetricBoundary<BoundDirichlet,BoundDirichlet> DDBoundary;
typedef SymmetricBoundary<BoundDirichlet,BoundNeumann>   DNBoundary;

typedef SymmetricBoundary<BoundNeumann,BoundPeriodic>    NPBoundary;
typedef SymmetricBoundary<BoundNeumann,BoundDirichlet>   NDBoundary;
typedef SymmetricBoundary<BoundNeumann,BoundNeumann>     NNBoundary;

#endif

