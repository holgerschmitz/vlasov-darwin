// -*- C++ -*-
// $Id$

#include "vlasovbase.h"

#ifndef SCHEME_H
#define SCHEME_H

template<class ForceField>
class PosFluxCons3rdOrder : public ForceField {
  protected:      
      /// A maximum value of the distribution function
      double f_infty;
      
      bool errmsg;
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

      /// The slope limiter towards the left
      double epsilonLeft(double fj, double fjp);
      /// The slope limiter towards the right
      double epsilonRight(double fj, double fjm);
      
      /// Third order flux-cunserving interpolation in the x-direction
      double interpolateX    (const PositionI &Xi,
                              const VelocityI &Vi,
                              int j,
                              double alpha);

      /// Third order flux-cunserving interpolation in the y-direction
      double interpolateY    (const PositionI &Xi,
                              const VelocityI &Vi,
                              int j,
                              double alpha);

      /// Third order flux-cunserving interpolation in the vx-direction
      double interpolateVx   (const PositionI &Xi,
                              const VelocityI &Vi,
                              int j,
                              double alpha);

      /// Third order flux-cunserving interpolation in the vy-direction
      double interpolateVy   (const PositionI &Xi,
                              const VelocityI &Vi,
                              int j,
                              double alpha);

      /// Third order flux-cunserving interpolation in the vz-direction
      double interpolateVz   (const PositionI &Xi,
                              const VelocityI &Vi,
                              int j,
                              double alpha);

      /// Third order flux-cunserving interpolation in the vx-direction
      double interpolateVxLow(const PositionI &Xi,
                              const VelocityI &Vi,
                              int j,
                              double alpha);

      /// Third order flux-cunserving interpolation in the vy-direction
      double interpolateVyLow(const PositionI &Xi,
                              const VelocityI &Vi,
                              int j,
                              double alpha);

      /// Third order flux-cunserving interpolation in the vz-direction
      double interpolateVzLow(const PositionI &Xi,
                              const VelocityI &Vi,
                              int j,
                              double alpha);

      /// Third order flux-cunserving interpolation in the vx-direction
      double interpolateVxHigh(const PositionI &Xi,
                              const VelocityI &Vi,
                              int j,
                              double alpha);

      /// Third order flux-cunserving interpolation in the vy-direction
      double interpolateVyHigh(const PositionI &Xi,
                              const VelocityI &Vi,
                              int j,
                              double alpha);

      /// Third order flux-cunserving interpolation in the vz-direction
      double interpolateVzHigh(const PositionI &Xi,
                              const VelocityI &Vi,
                              int j,
                              double alpha);
      /// The number of boundary cells needed to the right
      int bound_plus()   { return 1;}
      /// The number of boundary cells needed to the left
      int bound_minus()  { return 1;}
  public:
      PosFluxCons3rdOrder(SpeciesData &data) : ForceField(data), errmsg(false) {}
};

#include "scheme_temp.cpp"

#endif
