#include "advancer.h"
#include "scheme.h"
#include "forcefield.h"
#include "testcharacteristics.h"

#ifdef SIM_ELECTROSTATIC
typedef EFieldForce ForceField;
#endif

#ifdef SIM_DARWIN
typedef EMDarwinForce ForceField;
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

//typedef VlasovSpecies<ForceField,RungeKuttaAdvance,PosFluxCons3rdOrder> Vlasov;
typedef TestCharacteristics<ForceField> Vlasov;
typedef Vlasov* pVlasov;


