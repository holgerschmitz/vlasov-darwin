#include "advancer.h"
#include "scheme.h"
#include "forcefield.h"


//typedef EFieldForce ForceField;
//typedef EMDarwinForce ForceField;

typedef ConstEBFieldForce ForceField;

template<
  class ForceField, 
  template<class,template<class> class> class Advancer,
  template<class> class Scheme
>
class VlasovSpecies;

typedef VlasovSpecies<ForceField,LeapFrogAdvance,PosFluxCons3rdOrder> Vlasov;
typedef Vlasov* pVlasov;


