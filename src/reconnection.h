// -*- C++ -*-
// $Id$

//-----------------------------------------------------------------------------

#ifndef RECONNECTION_H
#define RECONNECTION_H

//-----------------------------------------------------------------------------

#include "numeric.h"
#include "vlasov.h"
#include "parameter.h"
#include "task.h"
#include "boundary.h"
#include "numboundary.h"
#include "wrapvlasov.h"

#ifndef SINGLE_PROCESSOR
#include <mpi.h>

//-----------------------------------------------------------------------------

/** @file reconnection.h
 *  @brief boundary classes implementing various reconnection setups.
 *
 */

//-----------------------------------------------------------------------------

/** @brief boundary class implementing a simple reconnection setup.
 *
 * The boundary conditions are periodic in the x-direction and
 * mirrored in the y-direction.
 * 
 * Mirrored here means, that for odd fields \f$F\f$ the condition is
 * \f$F = 0\f$, for even fields it reads \f$\partial_y F = 0\f$
 * 
 * A lot of methods are inherited by the MPIPeriodicSplitXBoundary.
 * exchangeY and ScalarFieldReduce contain the major changes
 */
class SimpleReconnectionBoundary : public MPIPeriodicSplitXYBoundary {
  private:
	/// symmetric mixed boundary from Dirichlet and Neumann type
      NDBoundary evenXBound;
	/// symmetric mixed boundary from Dirichlet and Neumann type
      DNBoundary evenYBound;
	///symmetric boundary of Dirichlet type
      DDBoundary evenZBound;
	///symmetric boundary of Neumann type
      NNBoundary ScalarBound;
  public:
      /**Default constructor used by the parsers
       */
      SimpleReconnectionBoundary();
      /** @brief Constructor using the command line arguments.
	     *
       *  Sets up all the local variables and allocates memory
       *  for the buffers
       */
      SimpleReconnectionBoundary(int argc, char **argv);
  
      /// Virtual destructor deleting all the allocated arrays
      ~SimpleReconnectionBoundary();
      
      /** @brief Wraps the boundaries in x-direction.
       * 
       * Since the current is an odd field, we have to ensure that
       * \f$\partial_x j_y \partial_x j_z= 0\f$ and \f$j_x=0\f$. To do this we
       * set \f$f(x=x_0+\Delta, v_x, v_y, v_z) = f(x=x_0-\Delta, -v_x, v_y, v_z)\f$
       * where \f$x_0\f$ is the location of the boundary.
       */
      void exchangeX(VlasovDist &field);

      /** @brief Wraps the boundaries in y-direction.
       * 
       * Since the current is an odd field, we have to ensure that
       * \f$\partial_y j_x = \partial_y j_z = 0\f$ and \f$j_y=0\f$. To do this we
       * set \f$f(y=y_0+\Delta, v_x, v_y, v_z) = f(y=y_0-\Delta, v_x, -v_y, v_z)\f$
       * where \f$y_0\f$ is the location of the boundary.
       */
      void exchangeY(VlasovDist &field);

      /// Adds the scalar fields and wraps them
      void ScalarFieldCombine(ScalarField &field) const;
      /// Wraps the scalar fields
      void ScalarFieldReduce(ScalarField &field) const;

      /// Returns periodic boundary conditions
      const NumBoundary& getNumBoundary(ScalarField &field) const;
      /// always return "false"
      bool periodicX() { return false; }
      /// always return "false"
      bool periodicY() { return false; }
};

//-----------------------------------------------------------------------------

/** @brief Boundary class implementing a reconnection setup.
 *
 */
class GEMReconnectionBoundary : public MPIPeriodicSplitXYBoundary {
  private:
	///a mixed boundary
      MixedBoundary<BoundNeumann,BoundNeumann,BoundDirichlet,BoundNeumann> evenXBound;
	///a mixed boundary
      MixedBoundary<BoundDirichlet,BoundDirichlet,BoundNeumann,BoundDirichlet> evenYBound;
	///a mixed boundary
      MixedBoundary<BoundDirichlet,BoundDirichlet,BoundDirichlet,BoundNeumann> evenZBound;
	///odd parity boundary
      NDBoundary oddYBound;
	///a mixed boundary
      MixedBoundary<BoundNeumann,BoundNeumann,BoundNeumann,BoundDirichlet> oddZBound;
	///even parity boundary
      DNBoundary evenXZBound;
	///scalar boundary
      NNBoundary ScalarBound;
  public:
      ///Default constructor used by the parsers
      GEMReconnectionBoundary();
      /** @brief Constructor using the command line arguments.
	 *
       *  Sets up all the local variables and allocates memory
       *  for the buffers
       */
      GEMReconnectionBoundary(int argc, char **argv);

      /// Virtual destructor deleting all the allocated arrays
      ~GEMReconnectionBoundary();

      /** @brief Wraps the boundaries in x-direction.
       * 
       * Since the current is an odd field, we have to ensure that
       * \f$\partial_x j_y \partial_x j_z= 0\f$ and \f$j_x=0\f$. To do this we
       * set \f$f(x=x_0+\Delta, v_x, v_y, v_z) = f(x=x_0-\Delta, -v_x, v_y, v_z)\f$
       * where \f$x_0\f$ is the location of the boundary.
       */
      void exchangeX(VlasovDist &field);

     /** @brief Wraps the boundaries in y-direction.
       * 
       * Since the current is an odd field, we have to ensure that
       * \f$\partial_y j_x = \partial_y j_z = 0\f$ and \f$j_y=0\f$. To do this we
       * set \f$f(y=y_0+\Delta, v_x, v_y, v_z) = f(y=y_0-\Delta, v_x, -v_y, v_z)\f$
       * where \f$y_0\f$ is the location of the boundary.
       */
      void exchangeY(VlasovDist &field);

      /// Adds the scalar fields and wraps them
      void ScalarFieldCombine(ScalarField &field) const;
      /// Wraps the scalar fields
      void ScalarFieldReduce(ScalarField &field) const;

      /// Returns periodic boundary conditions
      const NumBoundary& getNumBoundary(ScalarField &field) const;
      /// always return "false"
      bool periodicX() { return false; }
      /// always return "false"
      bool periodicY() { return false; }
};

//-----------------------------------------------------------------------------

/** Like GEMReconnectionBoundary but simulates the complete system.
 *  No symmetry requirements.
 */
class FullGEMReconnectionBoundary : public MPIPeriodicSplitXYBoundary {
  private:
      PNBoundary evenXBound;
      PDBoundary evenYBound;
      PNBoundary evenZBound;
      PDBoundary oddYBound;
      PDBoundary oddZBound;
      PNBoundary evenXZBound;
      PNBoundary ScalarBound;
  public:
      /**Default constructor used by the parsers
       */
      FullGEMReconnectionBoundary();
      /** @brief Constructor using the command line arguments.
       *  Sets up all the local variables and allocates memory
       *  for the buffers
       */
      FullGEMReconnectionBoundary(int argc, char **argv);
  
      /// Virtual destructor deleting all the allocated arrays
      ~FullGEMReconnectionBoundary();
      
      /** Wraps the boundaries in y-direction.
       * 
       * Since the current is an odd field, we have to ensure that
       * \f$\partial_y j_x = \partial_y j_z = 0\f$ and \f$j_y=0\f$. To do this we
       * set \f$f(y=y_0+\Delta, v_x, v_y, v_z) = f(y=y_0-\Delta, v_x, -v_y, v_z)\f$
       * where \f$y_0\f$ is the location of the boundary.
       */
      void exchangeY(VlasovDist &field);
      
      /// Adds the scalar fields and wraps them
      void ScalarFieldCombine(ScalarField &field) const;
      /// Wraps the scalar fields
      void ScalarFieldReduce(ScalarField &field) const;

      /// Returns periodic boundary conditions
      const NumBoundary& getNumBoundary(ScalarField &field) const;
      
      bool periodicX() { return false; }
      bool periodicY() { return false; }
      
};

//-----------------------------------------------------------------------------

/** GEM reconnection boundary with inflow at the top and outflow at the
 *  sides
 */
class OpenGEMReconnectionBoundary : public MPIPeriodicSplitXYBoundary {
  private:
      class EzNumBoundary : public NumBoundary {
        private:
          /// (inflow velocity) * (lobe magnetic field)
          mutable double inflowfactor;
          /// The inflow half width
          mutable double inflowwidth;

        public:
          EzNumBoundary(double inflowfactor_=-1.0, double inflowwidth_=1.0)
           : inflowfactor(inflowfactor_), inflowwidth(inflowwidth_)
          {}

          EzNumBoundary(const EzNumBoundary& bound)
           : inflowfactor(bound.inflowfactor), inflowwidth(bound.inflowwidth)
          {}

          double setParameters(double inflowfactor_, double inflowwidth_) const
          {
            inflowfactor = inflowfactor_;
            inflowwidth = inflowwidth_;
          }

          void apply(NumMatrix<double,2> &u) const;
          
          bool normalize() const { 
             return false; 
          }
      };


      MixedBoundary<BoundNeumann,BoundNeumann,BoundDirichlet,BoundNeumann> evenXBound;
      MixedBoundary<BoundNeumann,BoundDirichlet,BoundNeumann,BoundNeumann> evenYBound;
      MixedBoundary<BoundNeumann,BoundDirichlet,BoundDirichlet,BoundNeumann> evenZBound;
      NDBoundary oddYBound;
      NNBoundary oddZBound;
      MixedBoundary<BoundNeumann,BoundDirichlet,BoundNeumann,BoundNeumann> evenXZBound;
      NNBoundary ScalarBound;
      
      MixedBoundaryWithOffset<BoundNeumann,BoundNeumann,BoundDirichlet,BoundDirichlet> bxBound;
      EzNumBoundary ezBound;

      double vyinflow;
      double inflowwidth;
      double BxLobe;
  public:
      /**Default constructor used by the parsers
       */
      OpenGEMReconnectionBoundary();
      /** @brief Constructor using the command line arguments.
       *  Sets up all the local variables and allocates memory
       *  for the buffers
       */
      OpenGEMReconnectionBoundary(int argc, char **argv);
  
      /// Virtual destructor deleting all the allocated arrays
      ~OpenGEMReconnectionBoundary();
      
      /** Wraps the boundaries in x-direction.
       * 
       * Since the current is an odd field, we have to ensure that
       * \f$\partial_x j_y \partial_x j_z= 0\f$ and \f$j_x=0\f$. To do this we
       * set \f$f(x=x_0+\Delta, v_x, v_y, v_z) = f(x=x_0-\Delta, -v_x, v_y, v_z)\f$
       * where \f$x_0\f$ is the location of the boundary.
       */
      void exchangeX(VlasovDist &field);

      /** Wraps the boundaries in y-direction.
       * 
       * Since the current is an odd field, we have to ensure that
       * \f$\partial_y j_x = \partial_y j_z = 0\f$ and \f$j_y=0\f$. To do this we
       * set \f$f(y=y_0+\Delta, v_x, v_y, v_z) = f(y=y_0-\Delta, v_x, -v_y, v_z)\f$
       * where \f$y_0\f$ is the location of the boundary.
       */
      void exchangeY(VlasovDist &field);
      
      /// Adds the scalar fields and wraps them
      void ScalarFieldCombine(ScalarField &field) const;
      /// Wraps the scalar fields
      void ScalarFieldReduce(ScalarField &field) const;

      /// Returns periodic boundary conditions
      const NumBoundary& getNumBoundary(ScalarField &field) const;
      
      bool periodicX() { return false; }
      bool periodicY() { return false; }
  protected:
      PARAMETERMAP* MakeParamMap(PARAMETERMAP*);
      
};

#endif // SINGLE_PROCESSOR

//-----------------------------------------------------------------------------

/** @brief Initialises the Vlasov distribution function with a
 *  current sheet. 
 * 
 *  The thermal velocity and a streaming velocity
 *  of the distribution can be supplied.
 *
 *  The density is given by
 *  \f$$n(y) = n_0\text{sech}^2(y/\lambda) + n_{\infty}$\f$
 *  and the velocity by
 * \f$$v_z(y) = v_{z0}$\f$
 *  where \f$v_{z0}\f$ is specified and should be 
 *  \f$$v_{z0} = -\frac{m}{\sum m} \frac{B_0}{q\mu_0\lambda}$\f$
 * 
 * The distribution is made up of two Maxwellians, one shifted with the
 * velocity \f$v_{z0}\f$, creating the current sheet. The other is not shifted
 * and has the constant background density.
 */
class VlasovReconnectionInit : public VlasovInitialiser {
  protected:
      /** density in particles\f$m^{-{\rm dim}}\f$ at infinity ...
       */
      double Ninf, N0;
	///< ...and the perturbation in the current sheet
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

    /** @brief Do the initialisation
     *    
     *  Iterate through the whole distribution function and assign the appropriate
     *  phase space density to every point in phase space. The Phase space density 
     *  is calculated as a Maxwellian distribution.
     */
      void initialise(ForceFieldBase *pVlasov);
    protected:
	///build the parameter map
      PARAMETERMAP* MakeParamMap (PARAMETERMAP* pm = NULL);
};

//-----------------------------------------------------------------------------

/** @brief Initialises the Vlasov distribution function with a
 *  current sheet. 
 * 
 *  The thermal velocity and a streaming velocity
 *  of the distribution can be supplied.
 *
 *  The density is given by
 *  \f$$n(y) = n_0\text{sech}^2(y/\lambda) + n_{\infty}$\f$
 *  and the velocity by
 * \f$$v_z(y) = v_{z0}$\f$
 *  where \f$v_{z0}\f$ is specified and should be 
 *  \f$$v_{z0} = -\frac{m}{\sum m} \frac{B_0}{q\mu_0\lambda}$\f$
 * 
 * The distribution is made up of two Maxwellians, one shifted with the
 * velocity \f$v_{z0}\f$, creating the current sheet. The other is not shifted
 * and has the constant background density.
 */
class GEMReconnectionInit : public VlasovInitialiser {
  protected:
      /// density in particles\f$m^{-{\rm dim}}\f$ at infinity...
      double Ninf, 
	N0;  ///<...and the perturbation in the current sheet
      /// The thickness of the current sheet
      double lambda;
      
      /// The velocity carrying the sheet current
      double vz0;

      /// The initial perturbation of the current to start the reconnection
      double pert;
      
      /// thermal velocity
      VelocityD v_th;
      
      /// pointer to the owning VlasovSpecies class
      ForceFieldBase *pVlasov; 
  public:
      /// Default constructor
      GEMReconnectionInit();
      /// Destructor
      virtual ~GEMReconnectionInit();

    /** @brief Do the initialisation
     *    
     *  Iterate through the whole distribution function and assign the appropriate
     *  phase space density to every point in phase space. The Phase space density 
     *  is calculated as a Maxwellian distribution.
     */
      void initialise(ForceFieldBase *pVlasov);
    protected:
	    ///build the parameter map
      PARAMETERMAP* MakeParamMap (PARAMETERMAP* pm = NULL);
};

//-----------------------------------------------------------------------------

/** @brief Initialises the Vlasov distribution function with two opposing
 *  current sheets. 
 *
 *  The thermal velocity and a streaming velocity
 *  of the distribution can be supplied.
 *
 *  The density is given by
 *  \f$$n(y) = n_0\left( 
 *          \text{sech}^2((y-L/4)/\lambda) + 
 *          \text{sech}^2((y+L/4)/\lambda)
 *      \right) + n_{\infty}$\f$
 *  and the velocity by
 * \f$$v_z(y) = v_{z0}\frac{\left( 
 *          \text{sech}^2((y-L/4)/\lambda) - 
 *          \text{sech}^2((y+L/4)/\lambda)
 *      \right)}{n(y)}$\f$
 *  where \f$v_{z0}\f$ is specified and should be 
 *  \f$$v_{z0} = -\frac{m}{\sum m} \frac{B_0}{q\mu_0\lambda}$\f$
 */
class VlasovPeriodicReconnectionInit : public VlasovInitialiser {
  protected:
      /// density in particles\f$m^{-{\rm dim}}\f$ at infinity...
	  double Ninf, 
	         N0;///<...and the perturbation in the current sheet
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
      VlasovPeriodicReconnectionInit();
      /// Destructor
      virtual ~VlasovPeriodicReconnectionInit();

     /** @brief Perform distribution initialisation.
	 *
       *  The Maxwellian distribution will be written into the
       *  dist parameter.
       */
      void initialise(ForceFieldBase *pVlasov);
    protected:
	    ///build the parameter map
      PARAMETERMAP* MakeParamMap (PARAMETERMAP* pm = NULL);
};

//-----------------------------------------------------------------------------

/** @brief Like GEMReconnectionInit only for a full system. 
 *  Current sheet in the y-center.
 */
class FullGEMReconnectionInit : public VlasovInitialiser {
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
      double pert;
      
      /// thermal velocity
      VelocityD v_th;
      
      /// pointer to the owning VlasovSpecies class
      ForceFieldBase *pVlasov; 
  public:
    /// Default constructor
      FullGEMReconnectionInit();
   /// Destructor
      virtual ~FullGEMReconnectionInit();

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

//-----------------------------------------------------------------------------

/** @brief Initialises the Vlasov distribution function with a
 *  current sheet and inflow at the top. The thermal velocity and a 
 *  streaming velocity of the distribution can be supplied.
 *
 *  The density is given by
 *  \f$$n(y) = n_0\text{sech}^2(y/\lambda) + n_{\infty}$\f$
 *  and the velocity by
 * \f$$v_z(y) = v_{z0}$\f$
 *  where \f$v_{z0}\f$ is specified and should be 
 *  \f$$v_{z0} = -\frac{m}{\sum m} \frac{B_0}{q\mu_0\lambda}$\f$
 * 
 * The distribution is made up of two Maxwellians, one shifted with the
 * velocity \f$v_{z0}\f$, creating the current sheet. The other is not shifted
 * and has the constant background density.
 */
class OpenGEMReconnectionInit : public VlasovInitialiser {
  protected:
      /** density in particles\f$m^{-{\rm dim}}\f$ at infinity and 
       *  the perturbation in the current sheet
       */
      double Ninf, N0;
      /// The thickness of the current sheet
      double lambda;
      
      /// The velocity carrying the sheet current
      double vz0;
      
      /// The inflow velocity
      double vyinflow;
      
      /// The inflow half width
      double inflowwidth;

      /// The initial perturbation of the current to start the reconnection
      double pert;
      
      /// thermal velocity
      VelocityD v_th;
      
      /// pointer to the owning VlasovSpecies class
      ForceFieldBase *pVlasov; 
  public:
    /// Default constructor
      OpenGEMReconnectionInit();
   /// Destructor
      virtual ~OpenGEMReconnectionInit();

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

//-----------------------------------------------------------------------------

#endif // RECONNECTION_H
