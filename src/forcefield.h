#include "vlasovbase.h"

#ifndef FORCEFIELD_H
#define FORCEFIELD_H


class Potential;


/** @brief Implements the electrostatic force for plugging into the 
 *  VlasovSpecies
 */
class EFieldForce : public ForceFieldBase {
  public:
      typedef Potential FieldType;      
  protected:
      /// Pointer to the potential
      Potential* pPot;
  private:
      /// Scaling constant
      double dttx;
      DistMomentRho *Rho;
      ScalarField FEngy;
  public:
      
      EFieldForce(SpeciesData &data) : ForceFieldBase(data) {}
      
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
      
      DistMomentRho *getDerivedRho();
};

typedef PtrWrapper<EFieldForce> pEFieldForce;

/** @brief Implements the electrostatic force with a constant magnetic
 *  field for plugging into the VlasovSpecies
 */
class EBFieldForce : public ForceFieldBase {
  public:
      typedef Potential FieldType;      
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
      DistMomentRho *Rho;
  public:
      
      EBFieldForce(SpeciesData &data) : ForceFieldBase(data) {}
      
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
           
      DistMomentRho *getDerivedRho();
};

/** @brief Implements a forcefield with a constant electric and magnetic
 *  field for plugging into the VlasovSpecies
 */
class ConstEBFieldForce : public ForceFieldBase {
  public:
      class VoidPotential {
        public:
            VoidPotential() {}
            //VoidPotential(Boundary*) {}
            void Init() {}
            void Execute() {}
            void AddSpecies(ConstEBFieldForce*) {}
      };
      typedef VoidPotential FieldType;      
  protected:
      /// Pointer to the potential
      VoidPotential* pPot;
  private:
      VelocityD B;
      VelocityD E;
      DistMomentVelocities *Veloc;
      /// Scaling constant
      double dttx;
  public:
      
      ConstEBFieldForce(SpeciesData &data) : ForceFieldBase(data) {}
      
      /** The force at one position given the velocity.
       *  Actually the displacement in the velocity space is returned.
       *  The scheme of Boris is used for calculating the acceleration
       */
      VelocityD Force(const PositionI &Pos, 
                      const VelocityD &Vel,
                      double dt);

       /// Initialises the force field
      void Init(double dttx_);
            
      DistMomentVelocities *getDerivedVelocities();
};

class Darwin;

/** @brief Implements a force field that plugs into the concrete advancers
 *  taking the force from the Opar electric field.
 *
 *  Implements the Force method for the Advancer and the Init and MakeParamMap
 *  method to work together with the VlasovSpecies.
 */
class EMDarwinForce : public ForceFieldBase {
  public:
      typedef Darwin FieldType;
  
  protected:
      /// Pointer to the darwin field solver
      Darwin* pPot;
  
  private:
      /// Scaling constant
      double dttx;
      /// 
      DistMomentRho *Rho;
      DistMomentVelocities *Veloc;
      ScalarField FEngy;
  public:
      
      EMDarwinForce(SpeciesData &data) : ForceFieldBase(data) {}
      
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

      DistMomentRho *getDerivedRho();
      DistMomentVelocities *getDerivedVelocities();
};

typedef PtrWrapper<EMDarwinForce> pEMDarwinForce;

#endif
