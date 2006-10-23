// -*- C++ -*-
// $Id$
//-----------------------------------------------------------------------------
    /** @file wrapvlasov.h
	* @brief wraps advancers in a SpeciesBase class
	*
	*  The template class ForceField must declare a method "Force" that
 	*  takes a PositionI as argument an returns a VelocityD
 	*
 	*  The ForceField template class additionally has to supply the Init, 
 	*  and MakeParamMap method. The Init mehod has to take a scaling factor
 	*  of type double whith which to multiply the forcefield.
	*/
//-----------------------------------------------------------------------------
#include <string>
#include "numeric.h"
#include "scalarfield.h"
#include "potential.h"
#include "boundary.h"

#include "vlasovbase.h"
#include "scheme.h"
#include "pparabolic.h"
#include "advancer.h"

#include "parameter.h"
#include "task.h"
//-----------------------------------------------------------------------------
#ifndef WRAPVLASOV_H
#define WRAPVLASOV_H
//-----------------------------------------------------------------------------
//forward declarations
class Potential;
class VlasovInitialiser;
//-----------------------------------------------------------------------------
//Vlasovspecies
/** @brief This class wraps the Vlasov advancers in a SpeciesBase class
 *
 *  The template class ForceField must declare a method "Force" that
 *  takes a PositionI as argument an returns a VelocityD
 *
 *  The ForceField template class additionally has to supply the Init, 
 *  and MakeParamMap method. The Init mehod has to take a scaling factor
 *  of type double whith which to multiply the forcefield.
 *	
 */
template<
  class ForceField, 
  template<class,template<class> class> class Advancer = LeapFrogAdvance,
  template<class> class Scheme = PosFluxCons3rdOrder
>
class VlasovSpecies 
  : public Advancer<ForceField,Scheme>
{
  protected:      
      /// The timestep (in s)
      double dt;
      /// Goal for density
      double densityGoal;
      /// Contains the kinetic energy, if calculated
      ScalarField EKin;

  public:
      /** @brief Constructor,  
	  *
        * sets dt to the global timestep from parameters
	  * and densityGoal to the appropiate value from SpeciesData
	  */
      VlasovSpecies (SpeciesData &data);
      /// Virtual destructor
      virtual ~VlasovSpecies() {}

      /// Sets the force field and registers the species with it
      void setForceField(typename ForceField::FieldType *pField_);

      /** @brief Initializes the VlasovSpecies
       *
       *  Additionally calls the Init of the ForceField class template.
       *
	 *  @todo Test the normalisation.
       */
      virtual void Init ();

      /** @brief Perform one timestep.
	 *
       *  Calculates the potential and then advances the distribution function.
	 *  Calls the advance method of the Advancer base class,
	 *  then Task::Execute to execute sub-Tasks.
       */
      virtual void Execute();

      /// Writes the whole distribution function as text to a stream
      void write(ostream &);

  private:
      /** @brief Initialisation before first timestep, 
	  *
	  * sets f_infty to the \em global maximum of the distribution function
	  */
      void InterpolationInitStep(const VlasovDist &Dist);
};
//Vlasovspecies
//-----------------------------------------------------------------------------
//including implementations
#include "wrapvlasov_temp.cpp"
#endif //WRAPVLASOV_H
