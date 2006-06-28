// -*- C++ -*-
// $Id$

#ifndef SHOCK_H
#define SHOCK_H


#include "numeric.h"
#include "vlasov.h"
#include "parameter.h"
#include "task.h"
#include "boundary.h"
#include "numboundary.h"
#include "wrapvlasov.h"

#ifndef SINGLE_PROCESSOR

#include <mpi.h>

/** Boundary class implementing a simple shockwave setup.
 * The boundary conditions are periodic in the y-direction. There is a 
 * conducting wall at x=L and an constant inflow at x=0
 * 
 * A lot of methods are inherited by the MPIPeriodicSplitXYBoundary.
 * exchangeY and ScalarFieldReduce contain the major changes
 */
class MPIShockBoundary : public MPIPeriodicSplitXYBoundary {
  private:
      MixedBoundaryWithOffset<
        BoundNeumann, BoundNeumann,
        BoundPeriodic, BoundPeriodic
      > PotentialBound;
      
      MixedBoundaryWithOffset<
        BoundDirichlet, BoundNeumann, 
        BoundPeriodic, BoundPeriodic
      > BZBound;

      MixedBoundary<
        BoundDirichlet, BoundDirichlet, 
        BoundPeriodic, BoundPeriodic
      > ExBound;

      MixedBoundaryWithOffset<
        BoundNeumann, BoundDirichlet, 
        BoundPeriodic, BoundPeriodic
      > AyBound;
      
      DPBoundary evenXBound;           

      MixedBoundary<
        BoundDirichlet, BoundNeumann, 
        BoundPeriodic, BoundPeriodic
      > evenZBound;
      
      
      MixedBoundary<
        BoundNeumann, BoundDirichlet, 
        BoundPeriodic, BoundPeriodic
      > oddYBound;
      
      NPBoundary ScalarBound;
      
      double vx;
      VelocityD v_th;
      double density;
      
  public:
      /**Default constructor used by the parsers
       */
      MPIShockBoundary();
      /** @brief Constructor using the command line arguments.
       *  Sets up all the local variables and allocates memory
       *  for the buffers
       */
      MPIShockBoundary(int argc, char **argv);
  
      /// Virtual destructor deleting all the allocated arrays
      ~MPIShockBoundary();
      
      /** Wraps the boundaries in x-direction.
       *  Only exchangeX is needed here, exchangeY is periodic.
       */
      void exchangeX(VlasovDist &field);
      
      /// Adds the scalar fields and wraps them
      void ScalarFieldCombine(ScalarField &field) const;
      /// Wraps the scalar fields
      void ScalarFieldReduce(ScalarField &field) const;

      /// Returns shock boundary conditions
      const NumBoundary& getNumBoundary(ScalarField &field) const;
      
      bool periodicX() { return false; }
      bool periodicY() { return true; }
      
    protected:
      PARAMETERMAP* MakeParamMap (PARAMETERMAP* pm = NULL);
};

#endif // SINGLE_PROCESSOR

/** @brief Initialises the Vlasov distribution function for a shock.
 *  The thermal velocity and a streaming velocity
 *  of the distribution can be supplied.
 */
class VlasovShockInit : public VlasovInitialiser {
  protected:
      /** density in particles\f$m^{-{\rm dim}}\f$ at infinity and 
       *  the perturbation in the current sheet
       */
      double Ninf, N0;
      /// The thickness of the current sheet
      double lambda;
      
      /// The velocity carrying the sheet current
      double vz0;
      /// The initial perturbation of the current to start the reconnection
      double vz1;
      
      /// thermal velocity
      VelocityD v_th;
      
      /// pointer to the owning VlasovSpecies class
      ForceFieldBase *pVlasov; 
  public:
    /// Default constructor
      VlasovShockInit();
   /// Destructor
      virtual ~VlasovShockInit();

    /** @brief Do the initialisation
     *    
     *  Iterate through the whole distribution function and assign the appropriate
     *  phase space density to every point in phase space. The Phase space density 
     *  is calculated as a Maxwellian distribution.
     */
      void initialise(ForceFieldBase *pVlasov);
    protected:
      PARAMETERMAP* MakeParamMap (PARAMETERMAP* pm = NULL);
};

#endif // SHOCK_H

