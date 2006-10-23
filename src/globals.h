// -*- C++ -*-
// $Id$

#ifndef GLOBALS_H
#define GLOBALS_H
//-----------------------------------------------------------------------------
#include "task.h"
#include "vlasov.h"
#include "simtype.h"
#include "forcefield.h"
//-----------------------------------------------------------------------------
/** @file globals.h
	* @brief global parameter class
	*
	* Defines a class for storing global Parameters
	*/
//-----------------------------------------------------------------------------
//Parameters, singleton	
/// a class for storing parameters, inherits Rebuildable
class Parameters : public Rebuildable {
    private:
	/// cell count x
      int GridX;
	/// cell count y
      int GridY;
	///velocity resolution of the grid, for all directions
	int VelRes; 
	int VelResX, ///< velocity resolution of the grid, x-direction only
	    VelResY, ///< velocity resolution of the grid, y-direction only
	    VelResZ; ///< velocity resolution of the grid, z-direction only
	/// horizontal extension size of a single grid cell
      double GridSpace_x;
	/// vertical extension of a single grid cell
      double GridSpace_y;
	/// time step
      double Dt;
	/// initial density
      double BgDensity;
      /// total time steps
      int TotalTime;
      /// used for initilaizing the Distribution in vlas2d.cpp
      PositionI GridLow;
	/// used for initilaizing the Distribution in vlas2d.cpp
      PositionI GridHigh;
      /// used for initilaizing the Distribution in vlasovinit.cpp
      PhasePositionI DistLow;
	/// used for initilaizing the Distribution in vlasovinit.cpp
      PhasePositionI DistHigh;

      /** @brief velocity ratio over speed of light
	 *  
	 *  This velocity ratio is the ratio of the normalization velocity
       *  (normally the thermal velocity of the electrons) to the speed of
       *  light
       */
      double VelocityRatio;
      ///magnetic field (FixedArray<double, 3>)
      VelocityD BField;
	///electrostatic field (FixedArray<double, 3>)
      VelocityD EField;
      ///will be set "true" if its a restart
      bool IsRestart;
	//variables for the poisson multigrid solver, see poisson.h
	/** @brief For the Poisson Multigrid solver
	  *
	  * Do MGNu1 Gauss steps ...
	  */
      int MGNu1;
	/** @brief For the Poisson Multigrid solver
	 *
	 * ...then repeat multigrid refinement MGGamma times...
	 */
      int MGGamma;
	/** @brief For the Poisson Multigrid solver
	 *
	 * ...then do nu2 Gauss steps.
	 */
	int MGNu2;
	/** @brief For the Poisson Multigrid solver
	 *
	 * Error to reach in the Multigrid solver.
	 */
      double MGEpsilon;
      ///pointer to global parameters, will be set to "this" in constructor
      static Parameters *globals;
      //------------------------------------------------------------
      //typedefs
	///local typedef for the Field
      typedef ForceField::FieldType FieldType;
	///local typedef for the diagnostics
      typedef ForceField::DiagnosticType FieldDiagnosticType;
      //------------------------------------------------------------
      ///pointer to field type, see typedefs
      FieldType *Field;
      ///pointer to diagnostic type, see typedefs
      FieldDiagnosticType *FieldDiag;
      /// set to "true" if parameters class is initialized
      bool initialized;
      ///static member, number of the parameters main was called with
      static int Argc;
      /// static member, first argument main was called with
      static char **Argv;
    protected:
      ///Create the parameter map
      virtual PARAMETERMAP* MakeParamMap (PARAMETERMAP* pm = NULL);
    public:  
	//----------------------------------------------------------
	/// constructor, sets up a not initialized parameters object, restart is "false", sets globals to "this"
      Parameters();
	//----------------------------------------------------------
	///initialize and set initialized to "true"
      void init();
      //--------------------------------------------------------
      //accessor methods 
      ///accessor method, returns GridX
      int gridX() { return GridX; }
      ///accessor method, returns GridY
      int gridY() { return GridY; }
	///accessor method, returns grid spacing 
      double gridSpace_x() { return GridSpace_x; }
	///accessor method, returns grid spacing 
      double gridSpace_y() { return GridSpace_y; }
	/// accessor method, returns cell volume
      double volumeQuant() { return GridSpace_x*GridSpace_y; }
	/// accessor method, returns time step
      double dt() { return Dt; }
	/// accessor method, returns starting density
      double bgDensity() { return BgDensity; }
	/// accessor method, returns total time steps
      int totalTime() { return TotalTime; }
	/// accessor method, returns velocity over speed of light, normalized velocity
      double velocityRatio() { return VelocityRatio; }
	///accessor method, returns magnetic field
      const VelocityD& bField() { return BField; }
	///accessor method, returns electrostatic field 
      const VelocityD& eField() { return EField; }
      ///accessor method, returns isRestart
      bool isRestart() { return IsRestart; }
	///sets isRestart to its argument
      void setRestart(bool IsRestart_) { IsRestart = IsRestart_; }
      ///accessor method, MGNu1
      int mgNu1() { return MGNu1; }
      ///accessor method,  MGNu2
      int mgNu2() { return MGNu2; }
      ///accessor method, MGGamma
      int mgGamma() { return MGGamma; }
      ///accessor method, returns MGEpsilon
      double mgEpsilon() { return MGEpsilon; }
	///accessor method, returns GridLow
      const PositionI& gridLow() { return GridLow; }
      ///accessor method, returns GridHigh
      const PositionI& gridHigh() { return GridHigh; }
      ///accessor method, returns DistLow
      const PhasePositionI& distLow() { return DistLow; }
      ///raccessor method, eturns DistHigh
      const PhasePositionI& distHigh() { return DistHigh; }
      ///accessor method, returns pointer to the field
      FieldType *getField() { return Field; }
	///accessor method, eturns argc
      static int getArgc() { return Argc; }
      ///sets argc to its argument
      static void setArgc(int Argc_) { Argc = Argc_; }
      ///accessor method, returns argv
      static char **getArgv() { return Argv; }
      ///stes argv to its argument
      static void setArgv(char** Argv_) { Argv = Argv_; }
	
	///if parameters is not initialized, initialize. in any case: return parameters
      static Parameters &instance() { 
       if ( !(globals->initialized) ) globals->init();
        return *globals; 
      }
	//----------------------------------------------------------------------
      ///rebuild
      std::string Rebuild(std::istream& in);
};
//Parameters
//-----------------------------------------------------------------------------
//BoundaryKeeper
class Boundary;

/// a handle class for boundary objects
class BoundaryKeeper {
  private:
    ///pointer to a pointer to boundary
    static Boundary** ppboundary;
  public:
    /// set ppboundary to argument
    static void setBoundary(Boundary **ppboundary_) { ppboundary = ppboundary_; }
    ///return pointer to boundary
    static Boundary* getBoundary() { return *ppboundary; }
};
//BOundaryKeeper
//-----------------------------------------------------------------------------

#endif //GLOBALS_H
