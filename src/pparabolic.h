// -*- C++ -*-
// $Id$

#include "vlasovbase.h"

#ifndef PPARABOLIC_H
#define PPARABOLIC_H

template<class ForceField>
class PParabolicScheme : public ForceField {
  protected:      
      /// A maximum value of the distribution function
      double f_infty;
  protected:
      /// Advance the x--direction of the distribution function
      void advanceSpace_x(double timestep);
      /// Advance the y--direction of the distribution function
      void advanceSpace_y(double timestep);
      /// Accelerate the vx--direction of the distribution function
      void advanceVel_x(double timestep);
      /// Accelerate the vy--direction of the distribution function
      void advanceVel_y(double timestep);
      /// Accelerate the vz--direction of the distribution function
      void advanceVel_z(double timestep);
  private:
      class BoundaryBase {
        public:
          virtual void apply(NumMatrix<double, 1>&) = 0;
      };
      class BoundX : public BoundaryBase {
        private:
          Boundary *boundary;
        public:
          BoundX(Boundary *boundary_) : boundary(boundary_) {}
          void apply(NumMatrix<double, 1> &f) 
          { boundary->exchangeXLine(f); }
      };
      class BoundY : public BoundaryBase {
        private:
          Boundary *boundary;
        public:
          BoundY(Boundary *boundary_) : boundary(boundary_) {}
          void apply(NumMatrix<double, 1> &f) 
          { boundary->exchangeYLine(f); }
      };
      class BoundV : public BoundaryBase {
        public:
          void apply(NumMatrix<double, 1> &f) 
          {
            const int *lo = f.getLow();
            const int *hi = f.getHigh();
            f(lo[0]) = 0;
            f(lo[0]+1) = 0;
            f(hi[0]) = 0;
            f(hi[0]-1) = 0;
          }
      };
      
      NumMatrix<double, 1> du;
      NumMatrix<double, 1> uh;
      NumMatrix<double, 1> uL;
      NumMatrix<double, 1> uR;
      NumMatrix<double, 1> u6;
      NumMatrix<double, 1> Phi;
      
      void oneDimPpm(NumMatrix<double, 1> &f,
                     NumMatrix<double, 1> &v,
                     int lx, int hx,
                     BoundaryBase *bound);
  public:
      PParabolicScheme(SpeciesData &data) : ForceField(data) {}
};

#include "pparabolic.t"

#endif
