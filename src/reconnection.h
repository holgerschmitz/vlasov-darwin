// -*- C++ -*-
// $Id$

#ifndef RECONNECTION_H
#define RECONNECTION_H

#ifndef SINGLE_PROCESSOR

#include "numeric.h"
#include "vlasov.h"
#include "parameter.h"
#include "task.h"
#include "boundary.h"
#include "numboundary.h"
#include "wrapvlasov.h"

#include <mpi.h>

/** Boundary class implementing a simple reconnection setup.
 * The boundary conditions are periodic in the x-direction and
 * mirrored in the y-direction.
 * 
 * Mirrored here means, that for odd fields \f$F\f$ the condition is
 * \f$F = 0\f$, for even fields it reads \f$\partial_y F = 0\f$
 * 
 * A lot of methods are inherited by the MPIPeriodicSplitXBoundary.
 * exchangeY and ScalarFieldReduce contain the major changes
 */
class SimpleReconnectionBoundary : public MPIPeriodicSplitXBoundary {
  private:
      PDBoundary oddBound;
      PNBoundary evenBound;
//     PPBoundary oddBound;
//     PPBoundary evenBound;
  public:
      /**Default constructor used by the parsers
       */
      SimpleReconnectionBoundary();
      /** @brief Constructor using the command line arguments.
       *  Sets up all the local variables and allocates memory
       *  for the buffers
       */
      SimpleReconnectionBoundary(int argc, char **argv);
  
      /// Virtual destructor deleting all the allocated arrays
      ~SimpleReconnectionBoundary();
      
      /** Wraps the boundaries in y-direction.
       * 
       * Since the current is an odd field, we have to ensure that
       * \f$\partial_y j_z = 0\f$ and \f$j_x=j_y=0\f$. To do this we
       * set \f$f(y=y_0+\Delta, v_x, v_y, v_z) = f(y=y_0-\Delta, -v_x, v_y, -v_z)\f$
       * where \f$y_0\f$ is the location of the boundary.
       */
      void exchangeY(VlasovDist &field);
      
      /// Adds the scalar fields and wraps them
      void ScalarFieldReduce(ScalarField &field) const;

      /// Returns periodic boundary conditions
      const NumBoundary& getNumBoundary(ScalarField &field) const;
};

/** @brief Initialises the Vlasov distribution function with a
 *  current sheet. The thermal velocity and a streaming velocity
 *  of the distribution can be supplied.
 *
 *  The density is given by
 *  \f$$n(y) = n_0\text{sech}^2(y/\lambda) + n_{\infty}$\f$
 *  and the velocity by
 * \f$$v_z(y) = v_{z0}\frac{1}{n_0 + n_{\infty}\cosh^2(y/\lambda)}$\f$
 *  where \f$v_{z0}\f$ is specified and should be 
 *  \f$$v_{z0} = -\frac{m}{\sum m} \frac{B_0}{q\mu_0\lambda}$\f$
 */
class VlasovReconnectionInit : public VlasovInitialiser {
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
      VlasovReconnectionInit();
   /// Destructor
      virtual ~VlasovReconnectionInit();

      /** @brief Perform distribution initialisation.
       *  The Maxwellian distribution will be written into the
       *  dist parameter.
       */
      void initialise(ForceFieldBase *pVlasov);
    protected:
      PARAMETERMAP* MakeParamMap (PARAMETERMAP* pm = NULL);
};


#endif // SINGLE_PROCESSOR

#endif // RECONNECTION_H

