// -*- C++ -*-
// $Id$

#include "vlasovbase.h"

#ifndef FINDIFF_H
#define FINDIFF_H


template<class ForceField>
class FiniteDiffScheme : public ForceField {
  protected:
      /// A maximum value of the distribution function
      double f_infty;
      VlasovDist TempDist;
      VlasovDist CurrentDist;
      /** Advance the all directions of the distribution function
       * using the explicit finite difference scheme
       */
      void advanceAll(double timestep);
      
      void initDist();
  public:
      FiniteDiffScheme(SpeciesData &data) : ForceField(data) {}
};

#include "findiff.t"

#endif
