// -*- C++ -*-
// $Id$

#ifndef GLOBALS_H
#define GLOBALS_H

#include "task.h"
#include "vlasov.h"
#include "simtype.h"
#include "forcefield.h"

class Parameters : public Rebuildable {
    private:
      int GridX;
      int GridY;
      int VelRes;
      int VelResX, VelResY, VelResZ;
      double GridSpace_x;
      double GridSpace_y;
      double Dt;
      double BgDensity;
      
      int TotalTime;

      PositionI GridLow;
      PositionI GridHigh;
      PhasePositionI DistLow;
      PhasePositionI DistHigh;

      /** This velocity ratio is the ratio of the normalization velocity
       *  (normally the thermal velocity of the electrons) to the speed of
       *  light
       */
      double VelocityRatio;
      
      VelocityD BField;
      VelocityD EField;
      
      int MGNu1;
      int MGNu2;
      int MGGamma;
      double MGEpsilon;
      
      static Parameters *globals;
      
      typedef ForceField::FieldType FieldType;
      typedef ForceField::DiagnosticType FieldDiagnosticType;
      FieldType *Field;
      FieldDiagnosticType *FieldDiag;
      
      bool initialized;
      static int Argc;
      static char **Argv;
    protected:
      virtual PARAMETERMAP* MakeParamMap (PARAMETERMAP* pm = NULL);
    public:  
      Parameters();
      void init();
      int gridX() { return GridX; }
      int gridY() { return GridY; }
//      int velRes() { return VelRes; }
      double gridSpace_x() { return GridSpace_x; }
      double gridSpace_y() { return GridSpace_y; }
      double volumeQuant() { return GridSpace_x*GridSpace_y; }
      double dt() { return Dt; }
      double bgDensity() { return BgDensity; }
      int totalTime() { return TotalTime; }
      double velocityRatio() { return VelocityRatio; }
      
      const VelocityD& bField() { return BField; }
      const VelocityD& eField() { return EField; }
      
      int mgNu1() { return MGNu1; }
      int mgNu2() { return MGNu2; }
      int mgGamma() { return MGGamma; }
      double mgEpsilon() { return MGEpsilon; }

      const PositionI& gridLow() { return GridLow; }
      const PositionI& gridHigh() { return GridHigh; }
      const PhasePositionI& distLow() { return DistLow; }
      const PhasePositionI& distHigh() { return DistHigh; }

      FieldType *getField() { return Field; }

      static int getArgc() { return Argc; }
      static void setArgc(int Argc_) { Argc = Argc_; }
      static char **getArgv() { return Argv; }
      static void setArgv(char** Argv_) { Argv = Argv_; }

      static Parameters &instance() { 
       if ( !(globals->initialized) ) globals->init();
        return *globals; 
      }
      
      std::string Rebuild(std::istream& in);
};

class Boundary;

class BoundaryKeeper {
  private:
    static Boundary **ppboundary;
  public:
    static void setBoundary(Boundary **ppboundary_) { ppboundary = ppboundary_; }
    static Boundary* getBoundary() { return *ppboundary; }
};

#endif
