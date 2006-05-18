//#include "advancer.h"
#include "scheme.h"
#include "pparabolic.h"
#include "findiff.h"

template<class, template<class> class>
class SimpleAdvance;
 
template<class, template<class> class>
class SimpleLeapFrogAdvance;

template<class, template<class> class>
class RungeKuttaAdvance;

template<class, template<class> class>
class LeapFrogAdvance;

#define POS_FLUX_3_BACKSUBS

// This is the standard simulation method
#ifdef POS_FLUX_3_BACKSUBS
  #define SCHEME_TYPE PosFluxCons3rdOrder
  #define ADVANCE_TYPE SimpleLeapFrogAdvance
  #define FORCE_TYPE GenericEMForceBoris2
#endif

#ifdef POS_FLUX_3_BACKSUBS_EXACT
  #define SCHEME_TYPE PosFluxCons3rdOrder
  #define ADVANCE_TYPE SimpleLeapFrogAdvance
  #define FORCE_TYPE GenericEMForceExactBackSubs
#endif

#ifdef POS_FLUX_3_BACKSUBS_RK
  #define SCHEME_TYPE PosFluxCons3rdOrder
  #define ADVANCE_TYPE RungeKuttaAdvance
  #define FORCE_TYPE GenericEMForceBoris2
#endif

#ifdef POS_FLUX_3_TIMESPLIT
  #define SCHEME_TYPE PosFluxCons3rdOrder
  #define ADVANCE_TYPE LeapFrogAdvance
  #define FORCE_TYPE GenericEMForce
#endif

#ifdef POS_FLUX_3_TIMESPLIT_FIRST_ORDER
  #define SCHEME_TYPE PosFluxCons3rdOrder
  #define ADVANCE_TYPE LeapFrogAdvance
  #define FORCE_TYPE GenericEMForceDirect
#endif

#ifdef FINITE_DIFF
  #define SCHEME_TYPE FiniteDiffSchemeRK2
  #define ADVANCE_TYPE SimpleAdvance
  #define FORCE_TYPE GenericEMForceDirect
#endif

#ifdef FINITE_DIFF_RK3
  #define SCHEME_TYPE FiniteDiffSchemeRK3
  #define ADVANCE_TYPE SimpleAdvance
  #define FORCE_TYPE GenericEMForceDirect
#endif

#ifdef FINITE_DIFF4_RK3
  #define SCHEME_TYPE FiniteDiff4SchemeRK3
  #define ADVANCE_TYPE SimpleAdvance
  #define FORCE_TYPE GenericEMForceDirect
#endif


// Advancers can be one of the following:
//    SimpleAdvance (funktioniert nur mit FiniteDiffScheme)
//    LeapFrogAdvance (deprecated)
//    SimpleLeapFrogAdvance
//    RungeKuttaAdvance

// Scheme can be one of the following:
//    FiniteDiffScheme (funktioniert nur mit SimpleAdvance)
//    PosFluxCons3rdOrder
//    PParabolicScheme

class Darwin;
class GenericEMForceBase_ConstEB;
class GenericEMForceBase_Electrostatic;
template<class> class GenericEMForceBase_FullEM;

template<class> class GenericEMForce;
template<class> class GenericEMForceBoris2;
template<class> class GenericEMForceExact;
template<class> class GenericEMForceExactBackSubs;
template<class> class GenericEMForceDirect;

typedef FORCE_TYPE<
  GenericEMForceBase_ConstEB
> ConstEBFieldForce;

typedef FORCE_TYPE<
  GenericEMForceBase_Electrostatic
> ESWithBFieldForce;

typedef PtrWrapper<ConstEBFieldForce> pConstEBFieldForce;

typedef FORCE_TYPE<
  GenericEMForceBase_FullEM<Darwin>
> EMDarwinForce;

typedef PtrWrapper<EMDarwinForce> pEMDarwinForce;

class Magnetostatic;
typedef FORCE_TYPE<
  GenericEMForceBase_FullEM<Magnetostatic>
> MagnetostaticForce;

typedef PtrWrapper<MagnetostaticForce> pMagnetostaticForce;

#include "forcefield.h"

#ifdef SIM_ELECTROSTATIC
typedef EFieldForce ForceField;
#endif

#ifdef SIM_ELECTROSTATIC_B
typedef ESWithBFieldForce ForceField;
#endif

#ifdef SIM_DARWIN
typedef EMDarwinForce ForceField;
#endif

#ifdef SIM_MAGNETOSTATIC
typedef MagnetostaticForce ForceField;
#endif

#ifdef SIM_CONSTFIELD
typedef ConstEBFieldForce ForceField;
#endif


template<
  class ForceField, 
  template<class,template<class> class> class Advancer,
  template<class> class Scheme
>
class VlasovSpecies;


typedef VlasovSpecies<ForceField, ADVANCE_TYPE, SCHEME_TYPE> Vlasov;
typedef Vlasov* pVlasov;


