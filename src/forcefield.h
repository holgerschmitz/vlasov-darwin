#include "scalarfield.h"
#include "vlasovbase.h"
//-----------------------------------------------------------------------------
#ifndef FORCEFIELD_H
#define FORCEFIELD_H
//-----------------------------------------------------------------------------
  /** @file forcefield.h
	* @brief Implements various types of force fields.
	*
	* Implements the electrostatic force for plugging into the 
	* VlasovSpecies.
	*/
//-----------------------------------------------------------------------------
//class declarations
class Potential;

class EFieldDiagnostic;
class VoidFieldDiagnostic;
class EBFieldDiagnostic;
//-----------------------------------------------------------------------------
// electrostatic field
/** @brief Implements the electrostatic force for plugging into the 
 *  VlasovSpecies
 */
class EFieldForce : public ForceFieldBase {
  public:
	///local typedef for the field type
      typedef Potential FieldType;
	/** @brief local typedef  for the diagnostic, 
	  *
	  * this will cause the higher level classes to treat the 
	  * diagnostics correctly. (See globals.h)
	  */
      typedef EFieldDiagnostic DiagnosticType;
  protected:
      /// Pointer to the potential
      Potential* pPot;
  private:
      /// Scaling constant
      double dttx;
	/// density distribution
      DistMomentRho *Rho;
      /// scalarfield for eletric field energy
      ScalarField FEngy;
  public:
      ///constructor, initializes the base object
      EFieldForce(SpeciesData &data) : ForceFieldBase(data) {}
      /// Gets the x-component of the electric field from the Potential
      ScalarField &GetEx();
      /// Gets the y-component of the electric field from the Potential
      ScalarField &GetEy();

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
      
      ///returns Rho
      DistMomentRho *getDerivedRho();
};
//EFieldForce
//-----------------------------------------------------------------------------
/// wrapped pointer to a electrostatic force field object
typedef PtrWrapper<EFieldForce> pEFieldForce;
//-----------------------------------------------------------------------------
//GenericEMForceBase

/** @brief Implements the electrostatic force with a constant magnetic
 *  field for plugging into the VlasovSpecies
 */
class GenericEMForceBase_Electrostatic : public ForceFieldBase {
  public:
	///local typedef for the field type
      typedef Potential FieldType;
	///local typedef  for the diagnostic
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
      ///pointer to the potential
      DistMomentRho *Rho;
  public:
      ///constructor
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
      ///returns Rho
      DistMomentRho *getDerivedRho();
};
//GenericEMForceBase
//-----------------------------------------------------------------------------
//GenericEMForceBase_ConstEB, forward declaration
class GenericEMForceBase_ConstEB;
//-----------------------------------------------------------------------------
//VoidPotential
/** @brief Dummy potential field class
  *
  * This version is used, when a Potential field class is formally needed
  * but not do to anything.
  * 
  * All members of this class do nothing.
  */
class VoidPotential {
  public:
  	///Constructor
      VoidPotential() {}
      //VoidPotential(Boundary*) {}
      /// This version of init does nothing
      void Init() {}
      /// This version of execute does nothing 
      void Execute() {}
      /// This version of addSpecies does nothing 
      void AddSpecies(GenericEMForceBase_ConstEB*) {}
      ///get scalarfield by name, returns always "NULL"
      ScalarField *GetByName(const std::string&) { return NULL;}
};
//VoidPotential
//-----------------------------------------------------------------------------
//GenericEMForceBase_ConstEB
/** @brief Implements a force field with a constant electric and magnetic
 *  field for plugging into the VlasovSpecies
 */
class GenericEMForceBase_ConstEB : public ForceFieldBase {
  public:
	///local typedef for the field type
      typedef VoidPotential FieldType;
	///local typedef  for the diagnostic
      typedef VoidFieldDiagnostic DiagnosticType; 
  protected:
      /// Pointer to the potential
      VoidPotential* pPot;
  private:
  	  ///magnetic field
      VelocityD B;
      ///electric field
      VelocityD E;
      ///pointer to velocities
      DistMomentVelocities *Veloc;
      /// Scaling constant
      double dttx;
      
  public:
  	  ///constructor
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
      ///returns Veloc
      DistMomentVelocities *getDerivedVelocities();
};
//GenericEMForceBase_ConstEB
//-----------------------------------------------------------------------------
//GenericEMForceBase_FullEM
/** @brief Implements a full electromagnetic force field 
 *  field for plugging into the VlasovSpecies
 */
template<class FType>
class GenericEMForceBase_FullEM : public ForceFieldBase
{
  public:
	///local typedef for the field type
      typedef FType FieldType;
	///local typedef  for the diagnostic
      typedef EBFieldDiagnostic DiagnosticType; 
  protected:
      /// Pointer to the darwin field solver
      FieldType* pPot;
  private:
      /// Scaling constant
      double dttx;
      ///pointer to the potential
      DistMomentRho *Rho;
      ///pointer to the velocities
      DistMomentVelocities *Veloc;
  public:
  	  ///constructor
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
	///Returns density distribution
      DistMomentRho *getDerivedRho();
      ///Returns the derived velocities
      DistMomentVelocities *getDerivedVelocities();
};
//GenericEMForceBase_FullEM
//--------------------------------------------------------------
//GenericEMForce
/** @brief Implements a force field that plugs into the concrete advancers
 *
 *  taking the force from a electromagetic field specified as template
 *  parameter.
 *
 *  Implements the Force method for the Advancer and the Init and MakeParamMap
 *  method to work together with the VlasovSpecies.
 */
template<class ForceBaseType>
class GenericEMForce : public ForceBaseType {
  public:
      ///constructor
      GenericEMForce(SpeciesData &data) : ForceBaseType(data) {}
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
//GenericEMForce
//-----------------------------------------------------------------------------
//GenericEMForceBoris2
/** @brief Implements a force field that plugs into the concrete advancers
 *  taking the force from a electromagetic field specified as template
 *  parameter.
 *
 *  This class implements the Boris scheme approximation in  second order for 
 *  the electromagnetic force.
 * 
 *  Implements the Force method for the Advancer and the Init and MakeParamMap
 *  method to work together with the VlasovSpecies.
 */
template<class ForceBaseType>
class GenericEMForceBoris2 : public ForceBaseType {
  public:
      ///Constructor
      GenericEMForceBoris2(SpeciesData &data) : ForceBaseType(data) {}

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
//GenericEMForceBoris2
//-----------------------------------------------------------------------------
//GenericEMForceExact
/** @brief Implements a force field that plugs into the concrete advancers
 *  taking the force from a electromagetic field specified as template
 *  parameter.
 *
 *  
 *
 *  Implements the Force method for the Advancer and the Init and MakeParamMap
 *  method to work together with the VlasovSpecies.
 */
template<class ForceBaseType>
class GenericEMForceExact : public ForceBaseType {
  public:
      ///constructor
      GenericEMForceExact(SpeciesData &data) : ForceBaseType(data) {}
	
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
//GenericEMForceExact
//-----------------------------------------------------------------------------
//GenericEMForceExactBackSubs
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
      ///constructor
      GenericEMForceExactBackSubs(SpeciesData &data) : ForceBaseType(data) {}
      
      // Returns the reference to the electric field energy
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
//GenericEMForceExactBackSubs
//-----------------------------------------------------------------------------
//GenericEMForceDirect
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
      ///constructor
      GenericEMForceDirect(SpeciesData &data) : ForceBaseType(data) {}
      
      // Returns the reference to the electric field energy
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
//GenericEMforceDirect
//-----------------------------------------------------------------------------
//including template member functions from forcefield.t

#include "forcefield.t"
//-----------------------------------------------------------------------------
#endif //FORCEFIELD_H
