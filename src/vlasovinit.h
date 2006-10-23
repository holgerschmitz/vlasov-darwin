// -*- C++ -*-
// $Id$
//-----------------------------------------------------------------------------
/** @file vlasovinit.h
  * @brief Implements different initializers for the distribution function.
  *
  * Possible initializers are: maxwellian distribution, two maxwellian distributions, 
  * maxwellian + pertubation, inhomogenous temperature and from a saved HDF file.
  */
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
#ifndef VLASOVINIT_H
#define VLASOVINIT_H
#include "wrapvlasov.h"
//-----------------------------------------------------------------------------
//VlasovInitialiser
/** @brief Base class for initialising the Vlasov Distribution function.
  *  
  *  Defines only one abstract method initialise(), that is called by the 
  *  VlasovAdvancerBase class in the derived classes.
  */
class VlasovInitialiser : public Rebuildable {
  public:
   /** @brief  Has to be overridden by the subclass to initialise the distribution.
     *
     *  The distribution is passed by reference to allow changing.
     *  The VelRange parameter speciefies in which velocity range
     *  the distribution is defined: The reqion is given by 
     *  [ -VelRange : +VelRange ]
     */
    virtual void initialise(ForceFieldBase *pVlasov)=0;
    ///virtual method, returns always "false"
    virtual bool restart() { return false; }    
};
//VlasovInitialiser
//-----------------------------------------------------------------------------
//VlasovMaxwellInit
/** @brief Initialises the Vlasov distribution function with a Maxwellian distribution. 
 *
 *  The thermal velocity and a streaming velocity
 *  of the distribution can be supplied by global variables.
 *  Additionally the density can be perturbed by some value n_perturb.
 *  The wave vector of the perturbation is k_perturb.
 */
class VlasovMaxwellInit : public VlasovInitialiser {
  protected:
      /// density in particles\f$m^{-{\rm dim}}\f$
      double N;
      /// streaming vlocity  
      VelocityD u_stream;
      /// thermal velocity
      VelocityD v_th;
      /// The perturbation wave vector
      PositionI k_perturb;
      /// The perturbation density
      double n_perturb; 

      /// pointer to the owning VlasovSpecies class
      ForceFieldBase *pVlasov; 
  public:
	 /// Default constructor
      VlasovMaxwellInit();
	 /// Destructor
      virtual ~VlasovMaxwellInit();

      /** @brief Perform distribution initialisation.
	 *
       *  The Maxwellian distribution will be written into the
       *  dist parameter.
	 *  Iterates through the whole distribution function and assigns the appropriate
	 *  phase space density to every point in phase space. The Phase space density 
	 *  is calculated as a Maxwellian distribution.
       */
      virtual void initialise(ForceFieldBase *pVlasov);
  protected:
      /// Create the parameter map for the readable object
      virtual PARAMETERMAP* MakeParamMap (PARAMETERMAP* pm = NULL);
};
//VlasovMaxwellInit
//-----------------------------------------------------------------------------
//VlasovTwoMaxwellInit
/** @brief Initialises the Vlasov distribution function with two Maxwellian distributions. 
 * 
 *  The thermal velocities and a streaming velocities
 *  of the distributions can be supplied by global variables.
 *  Additionally the densities can be perturbed by some value n_perturb.
 *  The wave vector of the perturbations are k_perturb.
 */
class VlasovTwoMaxwellInit : public VlasovInitialiser {
  protected:
      /// density of first population in particles\f$m^{-\text{dim}}\f$
      double N1;   
      /// density of second population in particles\f$m^{-\text{dim}}\f$        
      double N2;
      /// streaming vlocity of first population  
      VelocityD u_stream1;
      /// streaming vlocity of second population
      VelocityD u_stream2;
      /// thermal velocity of first population  
      VelocityD v_th1;
      /// thermal velocity of second population
      VelocityD v_th2;
      /// The perturbation wave vector of first population  
      PositionI k_perturb1;
      /// The perturbation wave vector of second population
      PositionI k_perturb2;
      /// The perturbation density of first population  
      double n_perturb1;
      /// The perturbation density of second population
      double n_perturb2;

      /// pointer to the owning VlasovSpecies class
      ForceFieldBase *pVlasov; 
  public:
	  /// Default constructor
      VlasovTwoMaxwellInit();
	  /// Destructor
      virtual ~VlasovTwoMaxwellInit();

      /** @brief Perform distribution initialisation.
	 * 
       *  The Maxwellian distribution will be written into the
       *  dist parameter.
	 *  Iterates through the whole distribution function and assigns the appropriate
	 *  phase space density to every point in phase space. The Phase space density 
	 *  is calculated as a Maxwellian distribution.
       */
      virtual void initialise(ForceFieldBase *pVlasov);
    protected:
    /// Create the parameter map for the readable object
      virtual PARAMETERMAP* MakeParamMap (PARAMETERMAP* pm = NULL);
};
//VlasovTwoMaxwellInit
//-----------------------------------------------------------------------------
//VlasovWaveGenInit
/** @brief Initialises the Vlasov distribution function with a
 *  Maxwellian distribution and a perturbation containing a wide wave
 *  spectrum. 
 *  
 *  The thermal velocity and a streaming velocity
 *  of the distributions can be supplied by global variables.
 */
class VlasovWaveGenInit : public VlasovInitialiser {
  protected:
      /// density in particles\f$m^{-{\rm dim}}\f$
      double N;
      /// thermal velocity
      VelocityD v_th;
      /// streaming velocity
      VelocityD u_stream;
      /// The perturbation density  
      double n_perturb;

      /// pointer to the owning VlasovSpecies class
      ForceFieldBase *pVlasov; 
  public:
   /// Default constructor
      VlasovWaveGenInit();
      /// Destructor
      virtual ~VlasovWaveGenInit();

      /** @brief Perform distribution initialisation.
	 *
       *  The Maxwellian distribution will be written into the
       *  dist parameter.
       */
      virtual void initialise(ForceFieldBase *pVlasov);
  protected:
  /// Create the parameter map for the readable object
      virtual PARAMETERMAP* MakeParamMap (PARAMETERMAP* pm = NULL);
};
//VlasovWaveGenInit
//-----------------------------------------------------------------------------
//VlasovGaussInit

/** @brief Initialises the Vlasov distribution function with a
 *  inhomogeneuos temperature. 
 *
 *  The temperature distribution is
 *  initialised as a Gaussian distribution in space
 *  \f$T = a\left(1 + A\exp\left[\left(-(x_0-x)/b\right)^2\right]\right)\f$
 */
class VlasovGaussTempInit : public VlasovInitialiser {
  protected:
      /// density in particles\f$m^{-{\rm dim}}\f$
      double N;
      /// thermal velocity
      VelocityD v_th;
      /// Amplitude of temperature disturbance
      double t_amp;
      /// width of Gaussian
      double t_width;      
      /// streaming vlocity
      VelocityD u_stream;

      /// pointer to the owning VlasovSpecies class
      ForceFieldBase *pVlasov; 
  public:
      /// Default constructor
      VlasovGaussTempInit();
      /// Destructor
      virtual ~VlasovGaussTempInit();

      /** @brief Perform distribution initialisation.
	 *
       *  The Maxwellian distribution will be written into the
       *  dist parameter.
	 *  Iterates through the whole distribution function and assigns the appropriate
	 *  phase space density to every point in phase space. The Phase space density 
	 *  is calculated as a Maxwellian distribution.
       */
      virtual void initialise(ForceFieldBase *pVlasov);
  protected:
  /// Create the parameter map for the readable object
      virtual PARAMETERMAP* MakeParamMap (PARAMETERMAP* pm = NULL);
};
//VlasovGaussInit
//-----------------------------------------------------------------------------
//VlaosvHDFInit

/** @brief Initialises the Vlasov distribution function 
  *  from a HDF file
  */
class VlasovHDFInit : public VlasovInitialiser {
  protected:
      /// The filename
      std::string fname;
      /// pointer to the owning VlasovSpecies class
      ForceFieldBase *pVlasov; 
  public:
      /// Default constructor
      VlasovHDFInit();
      /// Destructor
      virtual ~VlasovHDFInit();

      /** @brief Perform distribution initialisation.
	 *
       *  Read from the HDF file
       */
      void initialise(ForceFieldBase *pVlasov);
      /** @brief return "true"
        *
        * since the distribution function is reinitialized from a file
        */
      bool restart() { return true; }  
    protected:
    /// Create the parameter map for the readable object
      PARAMETERMAP* MakeParamMap (PARAMETERMAP* pm = NULL);
};
//VlasovHDFInit
//-----------------------------------------------------------------------------

//  brief Initialises the Vlasov distribution function with 
//    two different regions containing a flow velocity u_stream and
//    -u_stream respectively. In each region the distributions are
//   Maxwellian distributions.
//   The thermal velocity and a streaming velocity
//   of the distributions can be supplied by global variables.
//  
// class VlasovCurrentSheetInit : public VlasovInitialiser {
//   protected:
//        density in particles\f$m^{-{\rm dim}}\f$
//       double N;
//        thermal velocity
//       VelocityD v_th;
//        streaming vlocity
//       VelocityD u_stream;  
// 
//        pointer to the owning VlasovSpecies class
//       VlasovSpecies<ForceField> *pVlasov; 
//   public:
// 	  Default constructor
//       VlasovCurrentSheetInit(VlasovSpecies *);
// 	   Destructor
//       virtual ~VlasovCurrentSheetInit();
// 
//        brief Perform distribution initialisation.
//          The Maxwellian distribution will be written into the
//         dist parameter.
//        
//       virtual void initialise(VlasovDist &dist,
//                               const VelocityD &VelRange_);
// };


#endif //VLASOVINIT_H
