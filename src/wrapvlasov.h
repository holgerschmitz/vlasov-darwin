// -*- C++ -*-
// $Id$


#include <string>
#include "numeric.h"
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
 *  The abstract methods MakeRho and Rho are redefined to integrate over
 *  the velocity space and return the resulting density.
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
      
      /// A maximum value of the distribution function
      double f_infty;
      
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
       *  @todo The potential solver has to be drawn out of this method, otherwise it will
       *  not work for multi-species simulations
       */
      virtual void Execute();

      void MakeRho();
      void MakeJs();

      /// Writes the whole distribution function as text to a stream
      void write(ostream &);
      
      /** @brief Writes a y-vy slice of the distribution function.
       *  The x, vx and vz values are passed as parameters
       */
      void writeYVySlice(int x,int vx, int vz, string fname);
      /** @brief Writes a vx-vy slice of the distribution function.
       *  The x, y and vz values are passed as parameters.
       *  The value t is prepended to each line
       */
      void writeVxVySlice(int t, int x, int y, int vz, ostream &Slice);      
      
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
