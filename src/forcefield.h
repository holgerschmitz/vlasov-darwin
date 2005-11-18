#include "scalarfield.h"
#include "vlasovbase.h"

#ifndef FORCEFIELD_H
#define FORCEFIELD_H


class Potential;

class EFieldDiagnostic;
class VoidFieldDiagnostic;
class EBFieldDiagnostic;


/** @brief Implements the electrostatic force for plugging into the 
 *  VlasovSpecies
 */
class EFieldForce : public ForceFieldBase {
  public:
      typedef Potential FieldType;
      typedef EFieldDiagnostic DiagnosticType;      
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
//      ScalarField &FieldEnergy();

      /** The force at one position given the velocity.
       *  Actually the displacement in the velocity space is returned
       */
      VelocityD Force(const PositionI &Pos, 
                      const VelocityD &Vel,
                      double dt);
      
      /** The forces individually by component. Calculation is the same as in 
       *  the Force method but all operations not needed by the x--component
       *  are omitted. This method is therefore slightly faster.
       */
      double ForceX(const PositionI &Pos, 
                    const VelocityD &Vel,
                    double dt);

      /** The forces individually by component. Calculation is the same as in 
       *  the Force method but all operations not needed by the y--component
       *  are omitted. This method is therefore slightly faster.
       */
      double ForceY(const PositionI &Pos, 
                    const VelocityD &Vel,
                    double dt);

      /** The forces individually by component. Calculation is the same as in 
       *  the Force method but all operations not needed by the z--component
       *  are omitted. This method is therefore slightly faster.
       */
      double ForceZ(const PositionI &Pos, 
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
class GenericEMForceBase_Electrostatic : public ForceFieldBase {
  public:
      typedef Potential FieldType;
      typedef EFieldDiagnostic DiagnosticType; 
  protected:
      /// Pointer to the potential
      Potential* pPot;
  private:
      /** @brief Magnetic Field normalized by 
       *  \f$\frac{v_A}{c}\frac{1}{B_0}\f$.
       */
      VelocityD B;
      /// Scaling constant
      double dttx;
      DistMomentRho *Rho;
  public:
      
      GenericEMForceBase_Electrostatic(SpeciesData &data) : ForceFieldBase(data) {}
      
      /// Gets the x-component of the electric field from the field solver
      double GetEx(const PositionI &Pos) ;
      
      /// Gets the y-component of the electric field from the field solver
      double GetEy(const PositionI &Pos);
      
      /// Gets the z-component of the electric field from the field solver
      double GetEz(const PositionI &Pos) { return 0.0; }
      /// Gets the x-component of the magnetic field from the field solver
      double GetBx(const PositionI &Pos) { return B[0]; }
      /// Gets the y-component of the magnetic field from the field solver
      double GetBy(const PositionI &Pos) { return B[1]; }
      /// Gets the z-component of the magnetic field from the field solver
      double GetBz(const PositionI &Pos) { return B[2]; }

       /// Initialises the force field
      void Init(double dttx_);
           
      DistMomentRho *getDerivedRho();
};

class GenericEMForceBase_ConstEB;
class VoidPotential {
  public:
      VoidPotential() {}
      //VoidPotential(Boundary*) {}
      void Init() {}
      void Execute() {}
      void AddSpecies(GenericEMForceBase_ConstEB*) {}
      ScalarField *GetByName(const std::string&) { return NULL;}
};

/** @brief Implements a forcefield with a constant electric and magnetic
 *  field for plugging into the VlasovSpecies
 */
class GenericEMForceBase_ConstEB : public ForceFieldBase {
  public:
      typedef VoidPotential FieldType;
      typedef VoidFieldDiagnostic DiagnosticType; 
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
      GenericEMForceBase_ConstEB(SpeciesData &data) : ForceFieldBase(data) {}

       /// Initialises the force field
      void Init(double dttx_);

      /// Gets the x-component of the electric field from the field solver
      double GetEx(const PositionI &Pos) { return E[0]; }
      /// Gets the y-component of the electric field from the field solver
      double GetEy(const PositionI &Pos) { return E[1]; }
      /// Gets the z-component of the electric field from the field solver
      double GetEz(const PositionI &Pos) { return E[2]; }
      /// Gets the x-component of the magnetic field from the field solver
      double GetBx(const PositionI &Pos) { return B[0]; }
      /// Gets the y-component of the magnetic field from the field solver
      double GetBy(const PositionI &Pos) { return B[1]; }
      /// Gets the z-component of the magnetic field from the field solver
      double GetBz(const PositionI &Pos) { return B[2]; }
            
      DistMomentVelocities *getDerivedVelocities();
};

template<class FType>
class GenericEMForceBase_FullEM : public ForceFieldBase
{
  public:
      typedef FType FieldType;
      typedef EBFieldDiagnostic DiagnosticType; 
  protected:
      /// Pointer to the darwin field solver
      FieldType* pPot;
  private:
      /// Scaling constant
      double dttx;
      DistMomentRho *Rho;
      DistMomentVelocities *Veloc;
  public:
      GenericEMForceBase_FullEM(SpeciesData &data) : ForceFieldBase(data) {}

      /// Initialises the force field
      void Init(double dttx_);

      /// Gets the x-component of the electric field from the field solver
      double GetEx(const PositionI &Pos);
      /// Gets the y-component of the electric field from the field solver
      double GetEy(const PositionI &Pos);
      /// Gets the z-component of the electric field from the field solver
      double GetEz(const PositionI &Pos);
      /// Gets the x-component of the magnetic field from the field solver
      double GetBx(const PositionI &Pos);
      /// Gets the y-component of the magnetic field from the field solver
      double GetBy(const PositionI &Pos);
      /// Gets the z-component of the magnetic field from the field solver
      double GetBz(const PositionI &Pos);

      DistMomentRho *getDerivedRho();
      DistMomentVelocities *getDerivedVelocities();
};



/** @brief Implements a force field that plugs into the concrete advancers
 *  taking the force from a electromagetic field specified as template
 *  parameter.
 *
 *  Implements the Force method for the Advancer and the Init and MakeParamMap
 *  method to work together with the VlasovSpecies.
 */
template<class ForceBaseType>
class GenericEMForce : public ForceBaseType {
  public:
      
      GenericEMForce(SpeciesData &data) : ForceBaseType(data) {}
      
      /// Returns the reference to the electric field energy
//      ScalarField &FieldEnergy();

      /** The force at one position given the velocity.
       *  Actually the displacement in the velocity space is returned.
       *  The scheme of Boris is used for calculating the acceleration
       *  The Boris scheme has been checked on simple trajectory integration.
       */
      VelocityD Force(const PositionI &Pos, 
                      const VelocityD &Vel,
                      double dt);
                      
      /** The forces individually by component. Calculation is the same as in 
       *  the Force method but all operations not needed by the x--component
       *  are omitted. This method is therefore slightly faster.
       */
      double ForceX(const PositionI &Pos, 
                    const VelocityD &Vel,
                    double dt);

      /** The forces individually by component. Calculation is the same as in 
       *  the Force method but all operations not needed by the y--component
       *  are omitted. This method is therefore slightly faster.
       */
      double ForceY(const PositionI &Pos, 
                    const VelocityD &Vel,
                    double dt);

      /** The forces individually by component. Calculation is the same as in 
       *  the Force method but all operations not needed by the z--component
       *  are omitted. This method is therefore slightly faster.
       */
      double ForceZ(const PositionI &Pos, 
                    const VelocityD &Vel,
                    double dt);
};

/** @brief Implements a force field that plugs into the concrete advancers
 *  taking the force from a electromagetic field specified as template
 *  parameter.
 *
 *  Implements the Force method for the Advancer and the Init and MakeParamMap
 *  method to work together with the VlasovSpecies.
 */
template<class ForceBaseType>
class GenericEMForceBoris2 : public ForceBaseType {
  public:
      
      GenericEMForceBoris2(SpeciesData &data) : ForceBaseType(data) {}
      
      /// Returns the reference to the electric field energy
//      ScalarField &FieldEnergy();

      /** The force at one position given the velocity.
       *  Actually the displacement in the velocity space is returned.
       *  The scheme of Boris is used for calculating the acceleration
       *  The Boris scheme has been checked on simple trajectory integration.
       */
      VelocityD Force(const PositionI &Pos, 
                      const VelocityD &Vel,
                      double dt);
                      
      /** The forces individually by component. Calculation is the same as in 
       *  the Force method but all operations not needed by the x--component
       *  are omitted. This method is therefore slightly faster.
       */
      double ForceX(const PositionI &Pos, 
                    const VelocityD &Vel,
                    double dt);

      /** The forces individually by component. Calculation is the same as in 
       *  the Force method but all operations not needed by the y--component
       *  are omitted. This method is therefore slightly faster.
       */
      double ForceY(const PositionI &Pos, 
                    const VelocityD &Vel,
                    double dt);

      /** The forces individually by component. Calculation is the same as in 
       *  the Force method but all operations not needed by the z--component
       *  are omitted. This method is therefore slightly faster.
       */
      double ForceZ(const PositionI &Pos, 
                    const VelocityD &Vel,
                    double dt);
};

/** @brief Implements a force field that plugs into the concrete advancers
 *  taking the force from a electromagetic field specified as template
 *  parameter.
 *
 *  Implements the Force method for the Advancer and the Init and MakeParamMap
 *  method to work together with the VlasovSpecies.
 */
template<class ForceBaseType>
class GenericEMForceExact : public ForceBaseType {
  public:
      
      GenericEMForceExact(SpeciesData &data) : ForceBaseType(data) {}
      
      /// Returns the reference to the electric field energy
//      ScalarField &FieldEnergy();

      /** The force at one position given the velocity.
       *  Actually the displacement in the velocity space is returned.
       *  The scheme of Boris is used for calculating the acceleration
       *  The Boris scheme has been checked on simple trajectory integration.
       */
      VelocityD Force(const PositionI &Pos, 
                      const VelocityD &Vel,
                      double dt);
                      
      /** The forces individually by component. Calculation is the same as in 
       *  the Force method but all operations not needed by the x--component
       *  are omitted. This method is therefore slightly faster.
       */
      double ForceX(const PositionI &Pos, 
                    const VelocityD &Vel,
                    double dt);

      /** The forces individually by component. Calculation is the same as in 
       *  the Force method but all operations not needed by the y--component
       *  are omitted. This method is therefore slightly faster.
       */
      double ForceY(const PositionI &Pos, 
                    const VelocityD &Vel,
                    double dt);

      /** The forces individually by component. Calculation is the same as in 
       *  the Force method but all operations not needed by the z--component
       *  are omitted. This method is therefore slightly faster.
       */
      double ForceZ(const PositionI &Pos, 
                    const VelocityD &Vel,
                    double dt);
};

/** @brief Implements a force field that plugs into the concrete advancers
 *  taking the force from a electromagetic field specified as template
 *  parameter.
 *
 *  Implements the Force method for the Advancer and the Init and MakeParamMap
 *  method to work together with the VlasovSpecies.
 */
template<class ForceBaseType>
class GenericEMForceExactBackSubs : public ForceBaseType {
  public:
      
      GenericEMForceExactBackSubs(SpeciesData &data) : ForceBaseType(data) {}
      
      /// Returns the reference to the electric field energy
//      ScalarField &FieldEnergy();

      /** The force at one position given the velocity.
       *  Actually the displacement in the velocity space is returned.
       *  The scheme of Boris is used for calculating the acceleration
       *  The Boris scheme has been checked on simple trajectory integration.
       */
      VelocityD Force(const PositionI &Pos, 
                      const VelocityD &Vel,
                      double dt);
                      
      /** The forces individually by component. Calculation is the same as in 
       *  the Force method but all operations not needed by the x--component
       *  are omitted. This method is therefore slightly faster.
       */
      double ForceX(const PositionI &Pos, 
                    const VelocityD &Vel,
                    double dt);

      /** The forces individually by component. Calculation is the same as in 
       *  the Force method but all operations not needed by the y--component
       *  are omitted. This method is therefore slightly faster.
       */
      double ForceY(const PositionI &Pos, 
                    const VelocityD &Vel,
                    double dt);

      /** The forces individually by component. Calculation is the same as in 
       *  the Force method but all operations not needed by the z--component
       *  are omitted. This method is therefore slightly faster.
       */
      double ForceZ(const PositionI &Pos, 
                    const VelocityD &Vel,
                    double dt);
};

/** @brief Implements a force field that plugs into the concrete advancers
 *  taking the force from a electromagetic field specified as template
 *  parameter.
 *
 *  Implements the Force method for the Advancer and the Init and MakeParamMap
 *  method to work together with the VlasovSpecies.
 */
template<class ForceBaseType>
class GenericEMForceDirect : public ForceBaseType {
  public:
      
      GenericEMForceDirect(SpeciesData &data) : ForceBaseType(data) {}
      
      /// Returns the reference to the electric field energy
//      ScalarField &FieldEnergy();

      /** The force at one position given the velocity.
       *  Actually the displacement in the velocity space is returned.
       *  The scheme of Boris is used for calculating the acceleration
       *  The Boris scheme has been checked on simple trajectory integration.
       */
      VelocityD Force(const PositionI &Pos, 
                      const VelocityD &Vel,
                      double dt);
                      
      /** The forces individually by component. Calculation is the same as in 
       *  the Force method but all operations not needed by the x--component
       *  are omitted. This method is therefore slightly faster.
       */
      double ForceX(const PositionI &Pos, 
                    const VelocityD &Vel,
                    double dt);

      /** The forces individually by component. Calculation is the same as in 
       *  the Force method but all operations not needed by the y--component
       *  are omitted. This method is therefore slightly faster.
       */
      double ForceY(const PositionI &Pos, 
                    const VelocityD &Vel,
                    double dt);

      /** The forces individually by component. Calculation is the same as in 
       *  the Force method but all operations not needed by the z--component
       *  are omitted. This method is therefore slightly faster.
       */
      double ForceZ(const PositionI &Pos, 
                    const VelocityD &Vel,
                    double dt);
};



#include "forcefield.t"

#endif
