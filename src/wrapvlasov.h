// -*- C++ -*-
// $Id$


#include <string>
#include "numeric.h"
#include "scalarfield.h"
#include "potential.h"
#include "boundary.h"

#include "vlasovbase.h"
#include "scheme.h"
#include "advancer.h"

#include "parameter.h"
#include "task.h"

#ifndef WRAPVLASOV_H
#define WRAPVLASOV_H

class Potential;
class VlasovInitialiser;


/** @brief This class wraps the Vlasov advancers in a SpeciesBase class
 *  so that Opar can use it.
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
  protected:
  
      /// contains the kinetic energy when calculated
      ScalarField EKin;   

  public:
      /// Default constructor
      VlasovSpecies (SpeciesData &data);
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
      
      /** Perform one timestep.
       *  Calculates the potential and then advances the distribution function.
       */
      virtual void Execute();

      /// Writes the whole distribution function as text to a stream
      void write(ostream &);
            
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
      
      /// Initialisation before a timestep
      void InterpolationInitStep(const VlasovDist &Dist);
}; // VlasovSpecies


void write_Distribution(VlasovDist &, const char *);


#include "wrapvlasov_temp.cpp"

#endif
