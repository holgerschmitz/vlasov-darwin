#include "advancer.h"
#include "scheme.h"
#include "pparabolic.h"
#include "findiff.h"
#include "forcefield.h"

#ifdef SIM_ELECTROSTATIC
typedef EFieldForce ForceField;
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

// Advancers can be one of the following:
//    SimpleAdvance (funktioniert nur mit FiniteDiffScheme)
//    LeapFrogAdvance (deprecated)
//    SimpleLeapFrogAdvance
//    RungeKuttaAdvance

// Scheme can be one of the following:
//    FiniteDiffScheme (funktioniert nur mit SimpleAdvance)
//    PosFluxCons3rdOrder
//    PParabolicScheme

typedef VlasovSpecies<ForceField,SimpleAdvance,FiniteDiffScheme> Vlasov;
typedef Vlasov* pVlasov;


