#include <string>
#include "numeric.h"
#include "boundary.h"
#include "derivedfields.h"

#ifndef VLASOVBASE_H
#define VLASOVBASE_H

class VlasovInitialiser;

struct SpeciesData {
  Boundary *bound;
  double charge;
  double mass;
  double GridRange_vx;
  double GridRange_vy;
  double GridRange_vz;
};


/** The Base class for the Force Field. All the 
 *  ForceField-classes that are used as template parameters to the
 *  Vlasov-Species should be derived from this class. It defines 
 *  fundamental constants (charge, mass) needed by the force fields 
 *  and the vlasov species
 */
class ForceFieldBase {
  protected:
      /// the size of the numerical grid in the space direction
      PositionD GridSize;
      /// the size of the numerical grid in the space direction
      PositionI GridSizeI;
      
      /// the size of the numerical grid in the velocity directions
      VelocityD VelSize; 
      /// half of VelSize 
      VelocityD VelSizeH;

      /// the physical extent of the spatial grid
      PositionD BoxRange;
      /// range in velocty space
      VelocityD VRange;   

      /// the physical space between grid points
      PositionD dx;       

      /// The distribution function
      VlasovDist Distribution;
      /// The Charge of the particle species
      double Charge;      
      /// The Mass of the particles
      double Mass;
      /** The boundary class, used for exchanging values across
       *  processes and for wrapping periodic systems
       */
      Boundary *boundary;
      
      /// The timestep
      int tstep;
      
      DerivedFieldsContainer derivedFields;
  public:
      ForceFieldBase(SpeciesData &data);
      /// Destructor
      virtual ~ForceFieldBase ();

      /// Initialises the distribution function with a given initialiser
      void initialise(VlasovInitialiser *init);
                        
      /** Resizes the grid of the distribution function
       *  (and the distribution function temporaries)
       */
      void resize(PhasePositionI &low, PhasePositionI &high);

// ACCESSOR METHODS
      
      /// Return the Charge
      double getCharge () const { return Charge; };

      VlasovDist &getDistribution() { return Distribution; }
      
      /** Returns the physical position of the grid cell
       *  specified by Xi
       */
      PositionD position(const PositionI &Xi) {
        return BoxRange*Xi/GridSize;
      }

      /** @brief Returns the size of the grid in the i-th 
       *  dimension in grid cells
       */
      int gridSize(int i) {
        return GridSizeI[i];
      }

      /** @brief The grid spacing in x--direction in physical
       *  dimensions
       */
      double deltaX() {
        return BoxRange[0]/GridSize[0];
      }

      /** @brief The grid spacing in y--direction in physical
       *  dimensions
       */
      double deltaY() {
        return BoxRange[1]/GridSize[1];
      }

      /** @brief The grid spacing in v_x--direction in physical
       *  dimensions
       */
      double deltaVx() {
        return VRange[0]/VelSizeH[0];
      }

      /** @brief The grid spacing in v_y--direction in physical
       *  dimensions
       */
      double deltaVy() {
        return VRange[1]/VelSizeH[1];
      }

      /** @brief The grid spacing in v_z--direction in physical
       *  dimensions
       */
      double deltaVz() {
        return VRange[2]/VelSizeH[2];
      }

      /**@brief Returns the velocity in the physically
       * normalized units.
       * (NOT in grid points per timestep!!)
       */
      VelocityD velocity(const VelocityI &Vi) {
        return (VRange*(VelocityD(-0.5,-0.5,-0.5)+Vi))/VelSizeH - VRange;
      }
      
      /// Returns true if this is the master process
      bool master() {
        return boundary->master();
      }
      
      /// Returns the pocess number
      int procnum() {
        return boundary->procnum();
      }

};

#endif
