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

template<class ForceField>
class FiniteDiffSchemeRK2 : public ForceField {
  private:
      void forwardDist(VlasovDist &rhs, VlasovDist &lhs, double timestep);
      int RKState;
  protected:
      /// A maximum value of the distribution function
      double f_infty;
      VlasovDist D1;
      /** Advance the all directions of the distribution function
       * using the explicit finite difference scheme
       */
      void advanceAll(double timestep);
      
      void initDist();
  public:
      FiniteDiffSchemeRK2(SpeciesData &data)
        : ForceField(data), RKState(0) {}
};

template<class ForceField>
class FiniteDiffSchemeRK3 : public ForceField {
  private:
      void forwardDist(VlasovDist &rhs, VlasovDist &lhs, double timestep);
      int RKState;
      double numDensity;
  protected:
      /// A maximum value of the distribution function
      double f_infty;
      
      VlasovDist T1Dist, T2Dist;
      
      /** Advance the all directions of the distribution function
       * using the explicit finite difference scheme
       */
      void advanceAll(double timestep);
      
      void initDist();
  public:
      FiniteDiffSchemeRK3(SpeciesData &data)
        : ForceField(data), RKState(0) {}
};

template<class ForceField>
class FiniteDiff4SchemeRK3 : public ForceField {
  private:
      void forwardDist(VlasovDist &rhs, VlasovDist &lhs, double timestep);
      int RKState;
      double numDensity;
  protected:
      /// A maximum value of the distribution function
      double f_infty;
      
      VlasovDist T1Dist, T2Dist;
      
      /** Advance the all directions of the distribution function
       * using the explicit finite difference scheme
       */
      void advanceAll(double timestep);
      
      void initDist();
  public:
      FiniteDiff4SchemeRK3(SpeciesData &data)
        : ForceField(data), RKState(0) {}
};

#include "findiff.t"

#endif
