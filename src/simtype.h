//#include "advancer.h"
#include "scheme.h"
#include "pparabolic.h"
#include "findiff.h"
//-----------------------------------------------------------------------------
/** @file simtype.h
  * @brief Defines the simulation method.
  *
  * Here the explicit simulation type is set up. Via preprocessor macro a type
  * is choosen from a prefined list of combinations of Scheme, 
  * Advancer and ForceFieldType.
  * Other combinations can be added by adding new macros.
  */
//-----------------------------------------------------------------------------
//forward declarations
template<class, template<class> class>
class SimpleAdvance;
 
template<class, template<class> class>
class SimpleLeapFrogAdvance;

template<class, template<class> class>
class RungeKuttaAdvance;

template<class, template<class> class>
class LeapFrogAdvance;
//=============================================================================
//-----------------------------------------------------------------------------
//macros selecting the scheme and advancer
// Advancers can be one of the following:
//    SimpleAdvance (works only with FiniteDiffScheme)
//    LeapFrogAdvance (deprecated)
//    SimpleLeapFrogAdvance
//    RungeKuttaAdvance
//
// Scheme can be one of the following:
//    FiniteDiffScheme (works only with SimpleAdvance)
//    PosFluxCons3rdOrder
//    PParabolicScheme
//-----------------------------------------------------------------------------
//=============================================================================
//===================== THE SELECTION OF THE SIMTYPE IS HERE===================
//=============================================================================

//############################CURRENTLY SELECTED###############################
// positive fluxconservative scheme of 3rd order with back substitution
///this is currently selected simulation type, the default is: POS_FLUX_3_BACKSUBS
#define POS_FLUX_3_BACKSUBS
//#############################################################################

//=============================================================================
//===================== THE SELECTION OF THE SIMTYPE IS HERE===================
//=============================================================================

//-----------------------------------------------------------------------------
// other setups, which can be selected
// each setup consists of a Scheme, an Advancer and a forcefield Type

// positive fluxconservative scheme of 3rd order with back substitution
#ifdef POS_FLUX_3_BACKSUBS
  ///sets the scheme to positive flux conservative in 3rd order
  #define SCHEME_TYPE PosFluxCons3rdOrder
  ///sets the advancer to the simple leap frog advancer
  #define ADVANCE_TYPE SimpleLeapFrogAdvance
  ///sets the force type to electromagnetic with boris approximation
  #define FORCE_TYPE GenericEMForceBoris2
#endif

#ifdef POS_FLUX_3_TESTING
  ///sets the scheme to positive flux conservative in 3rd order
  #define SCHEME_TYPE PosFluxCons3rdOrder
  ///sets the advancer to the simple leap frog advancer
  #define ADVANCE_TYPE SimpleLeapFrogAdvance
  ///sets the force type to electromagnetic with backsubstitution
  #define FORCE_TYPE GenericEMForceExactBackSubs
#endif

// back substitution with Runge Kutta
#ifdef POS_FLUX_3_BACKSUBS_RK
  ///sets the scheme to positive flux conservative in 3rd order
  #define SCHEME_TYPE PosFluxCons3rdOrder
  ///sets the advancer to runge kutta
  #define ADVANCE_TYPE RungeKuttaAdvance
  ///sets the force type to electromagnetic with boris approximation
  #define FORCE_TYPE GenericEMForceBoris2
#endif

//timesplitting scheme
#ifdef POS_FLUX_3_TIMESPLIT
  ///sets the scheme to positive flux conservative in 3rd order
  #define SCHEME_TYPE PosFluxCons3rdOrder
  ///sets the advancer to leap frog
  #define ADVANCE_TYPE LeapFrogAdvance
  ///sets the force type to electromagnetic 
  #define FORCE_TYPE GenericEMForce
#endif

//timesplitting in 1st order
#ifdef POS_FLUX_3_TIMESPLIT_FIRST_ORDER
  ///sets the scheme to positive flux conservative in 3rd order
  #define SCHEME_TYPE PosFluxCons3rdOrder
  ///sets the advancer to leap frog
  #define ADVANCE_TYPE LeapFrogAdvance
  ///sets the force type to electromagnetic with boris approximation
  #define FORCE_TYPE GenericEMForceDirect
#endif

//finite differences 
#ifdef FINITE_DIFF
  ///sets the scheme to finite differences
  #define SCHEME_TYPE FiniteDiffSchemeRK2
  ///sets the advancer to the simple advancer
  #define ADVANCE_TYPE SimpleAdvance
  ///sets the force type to electromagnetic with boris approximation
  #define FORCE_TYPE GenericEMForceDirect
#endif

//finite difference with three state Runge Kutta
#ifdef FINITE_DIFF_RK3
  ///sets the scheme to finite differences
  #define SCHEME_TYPE FiniteDiffSchemeRK3
  ///sets the advancer to the simple advancer
  #define ADVANCE_TYPE SimpleAdvance
  ///sets the force type to electromagnetic with boris approximation
  #define FORCE_TYPE GenericEMForceDirect
#endif

//finite difference of 4th order with three state Runge Kutta
#ifdef FINITE_DIFF4_RK3
  ///sets the scheme to finite differences
  #define SCHEME_TYPE FiniteDiff4SchemeRK3
  ///sets the advancer to the simple advancer
  #define ADVANCE_TYPE SimpleAdvance
  ///sets the force type to electromagnetic with boris approximation
  #define FORCE_TYPE GenericEMForceDirect
#endif
//------------------------------------------------------------------------------
//forward declarations

class Darwin;
class GenericEMForceBase_ConstEB;
class GenericEMForceBase_Electrostatic;
template<class> class GenericEMForceBase_FullEM;

template<class> class GenericEMForce;
template<class> class GenericEMForceBoris2;
template<class> class GenericEMForceExact;
template<class> class GenericEMForceExactBackSubs;
template<class> class GenericEMForceDirect;
//------------------------------------------------------------------------------
//typedefs for field types 

/// constant EB Field
typedef FORCE_TYPE<
  GenericEMForceBase_ConstEB
> ConstEBFieldForce;

///electrostatic force with B field
typedef FORCE_TYPE<
  GenericEMForceBase_Electrostatic
> ESWithBFieldForce;

///wrapped pointer to a ConstEBFieldForce object
typedef PtrWrapper<ConstEBFieldForce> pConstEBFieldForce;

///EM force with Darwin approximation
typedef FORCE_TYPE<
  GenericEMForceBase_FullEM<Darwin>
> EMDarwinForce;

///wrapped pointer to a EMDarwinForce object (electromagnetic Force with Darwin approximation)
typedef PtrWrapper<EMDarwinForce> pEMDarwinForce;

//forward declaration
class Magnetostatic;

///Magnetostatic Force
typedef FORCE_TYPE<
  GenericEMForceBase_FullEM<Magnetostatic>
> MagnetostaticForce;

///wrapped pointer to a MagnetostaticForce object
typedef PtrWrapper<MagnetostaticForce> pMagnetostaticForce;

//-----------------------------------------------------------------------------
#include "forcefield.h"
//-----------------------------------------------------------------------------

//=============================================================================
//=========================HERE THE FORCEFIELD IS SET==========================
//=============================================================================
//typedefs for typeswitching the ForceFieldTemplateparameter
//each sets the ForceField template parameter which plugs into the VlasovSpecies
//depending on the defined preprocessor variable from the simulation type
#ifdef SIM_ELECTROSTATIC
///local typedef for this Force Field type
typedef EFieldForce ForceField;
#endif

#ifdef SIM_ELECTROSTATIC_B
///local typedef for this Force Field type
typedef ESWithBFieldForce ForceField;
#endif

#ifdef SIM_DARWIN
///local typedef for this Force Field type
typedef EMDarwinForce ForceField;
#endif

#ifdef SIM_MAGNETOSTATIC
///local typedef for this Force Field type
typedef MagnetostaticForce ForceField;
#endif

#ifdef SIM_CONSTFIELD
///local typedef for this Force Field type
typedef ConstEBFieldForce ForceField;
#endif
//=============================================================================
//-----------------------------------------------------------------------------
//=============================================================================
//this is the base object all the policies plug into
/** @brief Here the simulation type is build from the defined "policies".
  * 
  * Forcefield, Advancer and Scheme are defined in this file (simtype.h).
  * These are plugged into the concrete ForceField object, which will be used
  * to specify the VlasovSpecies from this template.
  */
template<
  class ForceField, 
  template<class,template<class> class> class Advancer,
  template<class> class Scheme
>
class VlasovSpecies;
//--------------------------------------------------------------------------- 
//===========================================================================
//==========================VLASOVSPECIES TYPEDEF============================
//===========================================================================
///plug scheme and advancer from the macro definitions into VlasovSpecies
typedef VlasovSpecies<ForceField, ADVANCE_TYPE, SCHEME_TYPE> Vlasov;
/// global typedef for a pointer to a Vlasov object
typedef Vlasov* pVlasov;
//-----------------------------------------------------------------------------

