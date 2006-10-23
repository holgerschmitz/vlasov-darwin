#include <string>
#include <fstream>
#include "numeric.h"
#include "scalarfield.h"
#include "boundary.h"
#include "derivedfields.h"
#include "parameter.h"
#include "diagnostic.h"
#include "hdfstream.h"
//----------------------------------------------------------------------------
#ifndef VLASOVBASE_H
#define VLASOVBASE_H
//-----------------------------------------------------------------------------
    /** @file vlasovbase.h
	* @brief base class for forcefields
	*
	*  All the ForceField-classes that are used as template parameters to the
 	*  Vlasov-Species should be derived from this class. It defines 
 	*  fundamental constants (charge, mass) needed by the force fields 
 	*  and the vlasov species
 	*/
//-----------------------------------------------------------------------------
//typedefs
///phase space diagnostics, type def for a SimpleDiagnostic object
typedef SimpleDiagnostic<VlasovDist,HDFostream> PhaseDiag;
///wrapped pointer to ForceFieldBase
typedef PtrWrapper<ForceFieldBase> pForceFieldBase;
//-----------------------------------------------------------------------------
//forward declarations
class VlasovInitialiser;
class VlasovDerivedDiagnostic;
struct SpeciesData;
//-----------------------------------------------------------------------------
//ForceFieldBase
/** @brief The (abstract) Base class for the Force Field. 
  *
  *  All the ForceField-classes that are used as template parameters to the
  *  Vlasov-Species should be derived from this class. It defines 
  *  fundamental constants (charge, mass) needed by the force fields 
  *  and the vlasov species
  */
class ForceFieldBase {
  protected:
	//characteristic data of the grid
      /// the size of the numerical grid in the space direction
      PositionD GridSize;
      /// the size of the numerical grid in the space direction
      PositionI GridSizeI;
      /// the size of the numerical grid in the velocity directions
      VelocityD VelSize; 
      /// half of VelSize 
      VelocityD VelSizeH;
      ///center of the velocity grid
      VelocityD VelGridCenter;
      /// the physical extent of the spatial grid
      PositionD BoxRange;
      /// extent of the velocity grid
      VelocityD VRange;   
      /// the physical space between grid points
      PositionD dx;
	
      /** @brief The distribution function in phase space
	  *
	  * This is a five dimensional matrix of double values
	  */
      VlasovDist Distribution;
	
	//particle data
      /// The Charge of the particle species
      double Charge;
      /// The Mass of the particles
      double Mass;
      ///The boundary class, used for exchanging values across processes and for wrapping periodic systems
	
	//boundary and initialiser
      Boundary *boundary;
      /// A pointer to the initialiser
      VlasovInitialiser *init;
	
      /// The timestep
      int tstep;
      /// The diagnostic for the distribution function
      PhaseDiag *phasediag;
      ///derived fields
      DerivedFieldsContainer derivedFields;
  public:
	///local typedef for a list of diagnostic interfaces (VlasovDerivedDiagnostics)
      typedef std::list<VlasovDerivedDiagnostic*> DerivedDiagList;
  protected:
	///list for diagnostics
      DerivedDiagList diaglist;
  public:
	//Constructor and Destructor
  	///constructor, initialize the data members with the according values of the SpeciesDate
      ForceFieldBase(SpeciesData &data);
      ///virtual Destructor, deletes boundary
      virtual ~ForceFieldBase ();
      /** Initialises the distribution function with a given initialiser (init)
	  *
	  * Calls the initialize() method of the initializer for this object.
	  * The initializers are found in the vlasovinit.h file
	  */
	
	//Helper functions
      void initialise();
      /** Resizes the grid of the distribution function
        * (and the distribution function temporaries)
        */
      void resize(PhasePositionI low, PhasePositionI high);
      /** Returns the density error.
       *  The value returned is the total particle number density.
       *  So in many cases the value of 1 means no error at all
       */
      double densityError(VlasovDist &distrib);
     /** @brief Corrects any errors in the density.
	 *
       *  The FACTOR by which the density is perturbed is given as argument.
       *  The distribution will be divided by that factor
       */
      void correctDensityError(double err, VlasovDist &distrib);
	
	// ACCESSOR METHODS
      /// Returns the Charge
      double getCharge () const { return Charge; };
      /// Returns the Mass
      double getMass () const { return Mass; };
	/// Returns the Distribution function
      VlasovDist &getDistribution() { return Distribution; }
      
      /// Returns the physical position of the grid cellspecified by Xi
      PositionD position(const PositionI &Xi) {
        return BoxRange*Xi/GridSize;
      }
      ///Returns the size of the grid in the i-th  dimension in grid cells
      int gridSize(int i) {
        return GridSizeI[i];
      }
      ///The grid spacing in x-direction in physical dimensions
      double deltaX() {
        return BoxRange[0]/GridSize[0];
      }
      /// The grid spacing in y-direction in physical dimensions
      double deltaY() {
        return BoxRange[1]/GridSize[1];
      }
      /// The grid spacing in v_x--direction in physical dimensions
      double deltaVx() {
        return VRange[0]/VelSizeH[0];
      }
      /// The grid spacing in v_y-direction in physical  dimensions
      double deltaVy() {
        return VRange[1]/VelSizeH[1];
      }
      /// The grid spacing in v_z-direction in physical dimensions
      double deltaVz() {
        return VRange[2]/VelSizeH[2];
      }
     /** @brief Returns the velocity in the physically
       * normalized units.
       *
       * NOT in grid points per timestep!!
       */
      VelocityD velocity(const VelocityI &Vi) {
        return VelocityD(
          VRange[0]*(Vi[0]-VelGridCenter[0])/VelSizeH[0],
          VRange[1]*(Vi[1]-VelGridCenter[1])/VelSizeH[1],
          VRange[2]*(Vi[2]-VelGridCenter[2])/VelSizeH[2]
        );
//        return (VRange*(VelocityD(-0.5,-0.5,-0.5)+Vi))/VelSizeH - VRange;
      }
     /** @brief Returns the velocity in the physically
       * normalized units.
       *
       * NOT in grid points per timestep!!
       */
      VelocityD velocity(const VelocityD &Vi) {
        return VelocityD(
          VRange[0]*(Vi[0]-VelGridCenter[0])/VelSizeH[0],
          VRange[1]*(Vi[1]-VelGridCenter[1])/VelSizeH[1],
          VRange[2]*(Vi[2]-VelGridCenter[2])/VelSizeH[2]
        );
      }
      /// Returns true if this is the master process
      bool master() {
        return boundary->master();
      } 
      /// Returns the pocess number
      int procnum() {
        return boundary->procnum();
      }
      /// Returns DerivedField by name
      pDistributionDerivedField getDerivedField(std::string);
      /// Adds a diagnostic object to the list
      void addDerivedDiagnostic(VlasovDerivedDiagnostic*);
      ///purely virtual method
      virtual DistMomentRho *getDerivedRho() = 0;
};
//ForceFieldBase
//-----------------------------------------------------------------------------
//Species Data

/** @brief characteristics of a particle species
  *
  * This class is used to initialize the corresponding variables in the
  * classes derived from ForceFieldBase
  */
struct SpeciesData {
  /// single particle charge
  double charge;
  ///mass
  double mass;
  ///extent of the velocity grid in the x-direction, corresponds to VRange[0]
  double GridRange_vx;
  ///extent of the velocity grid in the y-direction, corresponds to VRange[1]
  double GridRange_vy;
  ///extent of the velocity grid in the z-direction, corresponds to VRange[2]
  double GridRange_vz;
  ///total density, used to enforce density conservation
  double densityGoal;
  ///pointer to the initializer
  VlasovInitialiser *init; 
  ///The phase space diagnostics
  PhaseDiag *phasediag;
};
//SpeciesData
//-----------------------------------------------------------------------------
//VlasovDerived Diagnostics

/** @brief diagnostic class for the derived fields
  * 
  * This is the unique object which performs diagnostics for 
  * Inherits the SimpleDiagnostic class.
  */
class VlasovDerivedDiagnostic : public SimpleDiagnostic<ScalarField,std::ofstream> 
{
  private:
	/// name of the field class
      std::string classname;
	/// name of the field
      std::string fieldname;
  protected:
	/// since this should be performed globaly, returns always "true"
      bool singleOut() { return true; }
  public:
	/** @brief local typedef for diagnostic interfaces list
	  *
	  * from ForceFieldBase:
        * typedef std::list<VlasovDerivedDiagnostic*> DerivedDiagList
	  */
      typedef ForceFieldBase::DerivedDiagList DerivedDiagList;
	/// list of all diagnostic objects
      static DerivedDiagList diaglist;
	/// static pointer to this object
      static VlasovDerivedDiagnostic *fielddiag;
  public:
	///default constructor, registers this object with the (static) DiagList
      VlasovDerivedDiagnostic();
	/** @brief retrieves the field specified by field- and classname, 
	  *
	  * sets the field in the inherited object accordingly
	  */
      void retrieveField(ForceFieldBase*);
	///create the parameter map
      virtual PARAMETERMAP* MakeParamMap (PARAMETERMAP* pm = NULL);
};
//VlasovDerivedDiagnostic
//-----------------------------------------------------------------------------
//CheckDensity

    /** @brief adds up the density on the inner cells and displays it
	*
	* The region of the electromagnetic fields is on
  	* grid cell smaller in every direction.
  	*
   	* It's important that the density is only calculated on the
  	* inner cells, since ScalarFieldReduce simply adds all the densities of
  	* the processes.
	*/
void CheckDensity(VlasovDist &dist, const char *Msg);
//-----------------------------------------------------------------------------

#endif
//VLASOVBASE_H
