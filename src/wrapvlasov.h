// -*- C++ -*-
// $Id$


#include <string>
#include "numeric.h"
#include "potential.h"
#include "boundary.h"

#ifndef WRAPVLASOV_H
#define WRAPVLASOV_H

class Potential;
class VlasovInitialiser;

struct SpeciesData {
  Boundary *bound;
  double charge;
  double mass;
  double GridRange_vx;
  double GridRange_vy;
  double GridRange_vz;
};


/** @brief This class wraps the Vlasov advancers in a SpeciesBase class
 *  so that Opar can use it.
 *
 *  The abstract methods MakeRho and Rho are redefined to integrate over
 *  the velocity space and return the resulting density.
 *
 *  The template class ForceField must declare a method "Force" that
 *  takes a PositionI as argument an returns a VelocityD
 *
 *  The ForceField template class additionally has to supply the Init, 
 *  and MakeParamMap method. The Init mehod has to take a scaling factor
 *  of type double whith which to multiply the forcefield.
 *
 */
template<class ForceField>
class VlasovSpecies : public ForceField {
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
      /// the physical space between grid points
      PositionD dx;       
      
      /// A maximum value of the distribution function
      double f_infty;
  protected:
      /// The distribution function
      VlasovDist Distribution;
      /// The Runge Kutta Temporaries
      VlasovDist T1, T2;
      /// Another temporary needed for the Runge Kutta integration
      VlasovDist TempDist;
      /** The state of the Runge Kutta integrator.
       *  Also used by the leap-frog integrator for the start-up phase
       */
      int RKState;
  protected:
      /// contains the particle density
      ScalarField gRho;
      /// contains the kinetic energy when calculated
      ScalarField EKin;   
      /// The current densities in x,y and z-direction
      ScalarField Jx, Jy, Jz;
      /// The mixed second order velocity moments of the distribution
      ScalarField Vxx, Vxy, Vxz;   
      ScalarField Vyy, Vyz, Vzz;
      /// The Charge of the particle species
      double Charge;      
      /// The Mass of the particles
      double Mass;
      /// The timestep (in s)
      double dt;
      /// range in velocty space
      VelocityD VRange;   

      /** The boundary class, used for exchanging values across
       *  processes and for wrapping periodic systems
       */
      Boundary *boundary;

      /// The timestep
      int t;
  public:
      /// Default constructor
      VlasovSpecies (SpeciesData &data);
      /// Destructor
      virtual ~VlasovSpecies ();

      /// Sets the force field and registers the species with it
      void setForceField(typename ForceField::FieldType *pField_);

      /// Initialises the distribution function with a given initialiser
      void initialise(VlasovInitialiser *init);
                        
      /** @brief Initializes the VlasovSpecies
       *
       *  Additionally calls the Init of the ForceField class template.
       *
       *  @todo Test the normalisation.
       */
      virtual void Init ();
      
      /** Perform one timestep.
       *  Calculates the potential and then advances the distribution function.
       *  @todo The potential solver has to be drawn out of this method, otherwise it will
       *  not work for multi-species simulations
       */
      virtual void Execute();

      void MakeRho();
      void MakeJs();

      /// Return a reference to the density grid
      ScalarField& Rho () { return gRho; };

      /** Return a reference to the grid containing the 
       *  x--component of the current density
       */
      ScalarField &getJx() { return Jx; }

      /** Return a reference to the grid containing the 
       *  y--component of the current density
       */
      ScalarField &getJy() { return Jy; }

      /** Return a reference to the grid containing the 
       *  z--component of the current density
       */
      ScalarField &getJz() { return Jz; }

      /// Return the Charge
      double getCharge () const { return Charge; };

      /// Return the time step
      double getdt () const { return dt; };

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
      
      /** Resizes the grid of the distribution function
       *  (and the distribution function temporaries)
       */
      void resize(PhasePositionI &low, PhasePositionI &high);

      /// Writes the whole distribution function as text to a stream
      void write(ostream &);
      
      /** @brief Writes a y-vy slice of the distribution function.
       *  The x, vx and vz values are passed as parameters
       */
      void writeYVySlice(int x,int vx, int vz, string fname);
      /** @brief Writes a vx-vy slice of the distribution function.
       *  The x, y and vz values are passed as parameters.
       *  The value t is prepended to each line
       */
      void writeVxVySlice(int t, int x, int y, int vz, ostream &Slice);

      /// Returns the vector of the current density at a grid position
      VelocityD getJ(int i, int j);
      
      /** @brief Returns the values of the second order velocity moments
       *  in a six-component vector
       */
      FixedArray<double,6> getVVTens(int i, int j);
      
      /** @brief Integrates the kinetic energy depending on the 
       *  position and returns it in a scalar field
       */
      ScalarField &KineticEnergy();
      
      /** @brief Adds the kinetic energy and the field energy
       *  (from the ForceField) and integrates over the whole 
       *  simulation box.
       */
      double TotalEnergy();

      /// Returns true if this is the master process
      bool master() {
        return boundary->master();
      }
      
      /// Returns the pocess number
      int procnum() {
        return boundary->procnum();
      }

  private:
      /** Returns the density error.
       *  The value returned is actually the density itself.
       *  So in normal cases the value of 1 means no error at all
       */
      double densityError();
  
      /** @brief Corrects any errors in the density.
       *  The FACTOR by which the density is perturbed is given as argument.
       *  The distribution will be divided by that factor
       */
      void correctDensityError(double err);
      
      /// Advance the distribution function one timestep
      void advance(double timestep);
      /** @brief Half of the distribution advance.
       *  Only advances the spatial direction for half a timestep
       */
      void advanceStepA(double timestep);

      /** @brief Half of the distribution advance.
       *  Only advances the velocity directions one timestep
       *  and then the spatial direction for half a timestep.
       */
      void advanceStepB(double timestep);

      /** @brief The full distribution advance.
       *  Only advances the velocity directions one timestep
       *  and then the spatial direction for a timestep.
       */
      void advanceStepFull(double timestep);

      /// Initialisation before a timestep
      void InterpolationInitStep(const VlasovDist &Dist);



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
                              double alpha,
                              const VlasovDist &Dist);

      /// Third order flux-cunserving interpolation in the y-direction
      double interpolateY    (const PositionI &Xi,
                              const VelocityI &Vi,
                              int j,
                              double alpha,
                              const VlasovDist &Dist);

      /// Third order flux-cunserving interpolation in the vx-direction
      double interpolateVx   (const PositionI &Xi,
                              const VelocityI &Vi,
                              int j,
                              double alpha,
                              const VlasovDist &Dist);

      /// Third order flux-cunserving interpolation in the vy-direction
      double interpolateVy   (const PositionI &Xi,
                              const VelocityI &Vi,
                              int j,
                              double alpha,
                              const VlasovDist &Dist);

      /// Third order flux-cunserving interpolation in the vz-direction
      double interpolateVz   (const PositionI &Xi,
                              const VelocityI &Vi,
                              int j,
                              double alpha,
                              const VlasovDist &Dist);
      /// The number of boundary cells needed to the right
      int bound_plus()   { return 1;}
      /// The number of boundary cells needed to the left
      int bound_minus()  { return 1;}
}; // VlasovSpecies


class Potential;

/** @brief Implements the electrostatic force for plugging into the 
 *  VlasovSpecies
 */
class EFieldForce {
  protected:
      /// Pointer to the potential
      Potential* pPot;
  private:
      /// Scaling constant
      double dttx;
      /// A field containing the field energy
      ScalarField FEngy;
  public:
      /// The default constructor
      EFieldForce() {};

      /// Gets the x-component of the electric field from the Potential
      ScalarField &GetEx();
      /// Gets the y-component of the electric field from the Potential
      ScalarField &GetEy();

      /// Returns the reference to the electric field energy
      ScalarField &FieldEnergy();

      /** The force at one position given the velocity.
       *  Actually the displacement in the velocity space is returned
       */
      VelocityD Force(const PositionI &Pos, 
                      const VelocityD &Vel,
                      double dt);
      
      /// Initialises the force field
      void Init(double dttx_);
      
      typedef Potential FieldType;
};

/** @brief Implements the electrostatic force with a constant magnetic
 *  field for plugging into the VlasovSpecies
 */
class EBFieldForce {
  protected:
      /// Pointer to the potential
      Potential* pPot;
      /** @brief Magnetic Field normalized by 
       *  \f$\frac{v_A}{c}\frac{1}{B_0}\f$.
       */
  private:
      VelocityD B;
      /// Scaling constant
      double dttx;
  public:
      /// The default constructor
      EBFieldForce() {};

      /// Gets the x-component of the electric field from the Potential
      ScalarField &GetEx();
      /// Gets the y-component of the electric field from the Potential
      ScalarField &GetEy();

      /** The force at one position given the velocity.
       *  Actually the displacement in the velocity space is returned.
       *  The scheme of Boris is used for calculating the acceleration
       */
      VelocityD Force(const PositionI &Pos, 
                      const VelocityD &Vel,
                      double dt);

       /// Initialises the force field
      void Init(double dttx_);
      
      /// Sets the magnetic field 
      void setBField(VelocityD B_);

      typedef Potential FieldType;
};

class Darwin;

/** @brief Implements a force field that plugs into the concrete advancers
 *  taking the force from the Opar electric field.
 *
 *  Implements the Force method for the Advancer and the Init and MakeParamMap
 *  method to work together with the VlasovSpecies.
 */
class EMDarwinForce {
  protected:
      /// Pointer to the darwin field solver
      Darwin* pPot;
  private:
      /// Scaling constant
      double dttx;
      /// A field containing the field energy
      ScalarField FEngy;
  public:
      /// The default constructor
      EMDarwinForce() {};

      /// Gets the x-component of the electric field from the field solver
      ScalarField &GetEx();
      /// Gets the y-component of the electric field from the field solver
      ScalarField &GetEy();
      /// Gets the z-component of the electric field from the field solver
      ScalarField &GetEz();
      /// Gets the x-component of the magnetic field from the field solver
      ScalarField &GetBx();
      /// Gets the y-component of the magnetic field from the field solver
      ScalarField &GetBy();
      /// Gets the z-component of the magnetic field from the field solver
      ScalarField &GetBz();

      /// Returns the reference to the electric field energy
      ScalarField &FieldEnergy();

      /** The force at one position given the velocity.
       *  Actually the displacement in the velocity space is returned.
       *  The scheme of Boris is used for calculating the acceleration
       */
      VelocityD Force(const PositionI &Pos, 
                      const VelocityD &Vel,
                      double dt);

      /// Initialises the force field
      void Init(double dttx_);

      typedef Darwin FieldType;
};

void write_Distribution(VlasovDist &, const char *);

#include "wrapvlasov_temp.cpp"

#endif
