// -*- C++ -*-
// $Id$
//-----------------------------------------------------------------------------
#include "vlasovbase.h"
//-----------------------------------------------------------------------------
/** @file
 * @brief Piecewiese parabolic integration scheme (not operational)
 *
 *  Implements a integration scheme for spatial integration, for plugging in the VlasovSpecies.
 */
//-----------------------------------------------------------------------------
#ifndef PPARABOLIC_H
#define PPARABOLIC_H
//-----------------------------------------------------------------------------
//PParabolicScheme
/** @brief Integration scheme class
  *
  * Implements a parabolic spatial inegration scheme
  */
template<class ForceField>
class PParabolicScheme : public ForceField {
  protected:      
      /// A maximum value of the distribution function
      double f_infty;
  protected:
      /// Advance the x-direction of the distribution function
      void advanceSpace_x(double timestep);
      /// Advance the y-direction of the distribution function
      void advanceSpace_y(double timestep);
      /// Accelerate the vx-direction of the distribution function
      void advanceVel_x(double timestep);
      /// Accelerate the vy-direction of the distribution function
      void advanceVel_y(double timestep);
      /// Accelerate the vz-direction of the distribution function
      void advanceVel_z(double timestep);
  private:
	//some local boundary classes
	  /// local nested boundary class, abstract base class
      class BoundaryBase {
        public:
	    ///apply a matrix
          virtual void apply(NumMatrix<double, 1>&) = 0;
      };
	///local nested boundary class, x-direction
      class BoundX : public BoundaryBase {
        private:
	    /// pointer to the boundary
          Boundary *boundary;
        public:
	    ///constructor
          BoundX(Boundary *boundary_) : boundary(boundary_) {}
	    ///apply a matrix
          void apply(NumMatrix<double, 1> &f) 
          { boundary->exchangeXLine(f); }
      };
	///local nested boundary class, y-direction
      class BoundY : public BoundaryBase {
        private:
	    ///pointer to the boundary
          Boundary *boundary;
        public:
	    ///constructor
          BoundY(Boundary *boundary_) : boundary(boundary_) {}
	    ///apply a matrix
          void apply(NumMatrix<double, 1> &f) 
          { boundary->exchangeYLine(f); }
      };
	/// local nested boundary class
      class BoundV : public BoundaryBase {
        public:
	    ///apply a matrix
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
      
	///piecewise parabolic in one dimension
      void oneDimPpm(NumMatrix<double, 1> &f,
                     NumMatrix<double, 1> &v,
                     int lx, int hx,
                     BoundaryBase *bound);
  public:
	///constructor
      PParabolicScheme(SpeciesData &data) : ForceField(data) {}
	//no detructor?
};
//PParabolicScheme
//-----------------------------------------------------------------------------
#include "pparabolic.t"
//-----------------------------------------------------------------------------
#endif //PPARABOLIC_H
