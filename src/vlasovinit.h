// -*- C++ -*-
// $Id$

#include "wrapvlasov.h"

/** @brief Implements the VlasovInitialiser interface to create a two stream
 *  distribution.
 *
 *  This class also inherits from the Task class, so it can work together with
 *  the VlasovSpecies class.
 */

#ifndef VLASOVINIT_H
#define VLASOVINIT_H

/** @brief Base class for initialising the Vlasov Distribution function
 *  
 *  Defines only an abstract method initialise, that is called by the 
 *  ValsovAdvancerBase class.
 */
class VlasovInitialiser {
    public:
        /** @brief  Should be overridden by the subclass to initialise
         *  the distribution.
         *
         *  The distribution is passed by reference to allow changing.
         *  The VelRange parameter speciefies in which velocity range
         *  the distrbution is defined: The reqion is given by 
         *  [ -VelRange : +VelRange ]
         */
        virtual void initialise(VlasovDist &dist,
                                const VelocityD &VelRange_)=0;
};

template<class ForceField>
struct VlasovMaxwellInitData {
  VlasovSpecies<ForceField> *vlasov;
  double InitStream_vx;
  double InitStream_vy;
  double InitStream_vz;
  double InitTherm_vx;
  double InitTherm_vy;
  double InitTherm_vz;
  double Init_N0;
  double Init_N1;
  int Init_kx;
  int Init_ky;
};

/** @brief Initialises the Vlasov distribution function with a
 *  Maxwellian distribution. The thermal velocity and a streaming velocity
 *  of the distribution can be supplied by global variables.
 *  Additionally the density can be perturbed by some value n_perturb.
 *  The wave vector of the perturbation is k_perturb.
 */
template<class ForceField>
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
      VlasovSpecies<ForceField> *pVlasov; 
  public:
	  /// Default constructor
      VlasovMaxwellInit(VlasovMaxwellInitData<ForceField> &data);
	  /// Destructor
      virtual ~VlasovMaxwellInit();

      /** @brief Perform distribution initialisation.
       *  The Maxwellian distribution will be written into the
       *  dist parameter.
       */
      virtual void initialise(VlasovDist &dist,
                              const VelocityD &VelRange_);
};

template<class ForceField>
struct VlasovTwoMaxwellInitData {
  VlasovSpecies<ForceField> *vlasov;
  double InitStream_vx;
  double InitStream_vy;
  double InitStream_vz;
  double InitTherm_vx;
  double InitTherm_vy;
  double InitTherm_vz;
  double Init_N0;
  double Init_N1;
  int Init_kx;
  int Init_ky;
  double InitStream2_vx;
  double InitStream2_vy;
  double InitStream2_vz;
  double InitTherm2_vx;
  double InitTherm2_vy;
  double InitTherm2_vz;
  double Init2_N0;
  double Init2_N1;
  int Init2_kx;
  int Init2_ky;
};


/** @brief Initialises the Vlasov distribution function with two
 *  Maxwellian distributions. The thermal velocities and a streaming velocities
 *  of the distributions can be supplied by global variables.
 *  Additionally the densities can be perturbed by some value n_perturb.
 *  The wave vector of the perturbations are k_perturb.
 */
template<class ForceField>
class VlasovTwoMaxwellInit : public VlasovInitialiser {
  protected:
      /// density of first population in particles\f$m^{-{\rm dim}}\f$
      double N1;   
      /// density of second population in particles\f$m^{-{\rm dim}}\f$        
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
      VlasovSpecies<ForceField> *pVlasov; 
  public:
	  /// Default constructor
      VlasovTwoMaxwellInit(VlasovTwoMaxwellInitData<ForceField> &data);
	  /// Destructor
      virtual ~VlasovTwoMaxwellInit();

      /** @brief Perform distribution initialisation.
       *  The Maxwellian distribution will be written into the
       *  dist parameter.
       */
      virtual void initialise(VlasovDist &dist,
                              const VelocityD &VelRange_);
};

/** @brief Initialises the Vlasov distribution function with a
 *  Maxwellian distribution and a perturbation containing a wide wave
 *  spectrum. The thermal velocity and a streaming velocity
 *  of the distributions can be supplied by global variables.
 */
template<class ForceField>
class VlasovWaveGenInit : public VlasovInitialiser {
  protected:
      /// density in particles\f$m^{-{\rm dim}}\f$
      double N;
      /// thermal velocity
      VelocityD v_th;
      /// streaming vlocity
      VelocityD u_stream;
      /// The perturbation density  
      double n_perturb;

      /// pointer to the owning VlasovSpecies class
      VlasovSpecies<ForceField> *pVlasov; 
  public:
   /// Default constructor
      VlasovWaveGenInit(VlasovMaxwellInitData<ForceField> &data);
   /// Destructor
      virtual ~VlasovWaveGenInit();

      /** @brief Perform distribution initialisation.
       *  The Maxwellian distribution will be written into the
       *  dist parameter.
       */
      virtual void initialise(VlasovDist &dist,
                              const VelocityD &VelRange_);
};

// /** @brief Initialises the Vlasov distribution function with 
//  *  two different regions containing a flow velocity u_stream and
//  *  -u_stream respectively. In each region the distributions are
//  *  Maxwellian distributions.
//  *  The thermal velocity and a streaming velocity
//  *  of the distributions can be supplied by global variables.
//  */
// template<class ForceField>
// class VlasovCurrentSheetInit : public VlasovInitialiser {
//   protected:
//       /// density in particles\f$m^{-{\rm dim}}\f$
//       double N;
//       /// thermal velocity
//       VelocityD v_th;
//       /// streaming vlocity
//       VelocityD u_stream;  
// 
//       /// pointer to the owning VlasovSpecies class
//       VlasovSpecies<ForceField> *pVlasov; 
//   public:
// 	  /// Default constructor
//       VlasovCurrentSheetInit(VlasovSpecies<ForceField> *);
// 	  /// Destructor
//       virtual ~VlasovCurrentSheetInit();
// 
//       /** @brief Perform distribution initialisation.
//        *  The Maxwellian distribution will be written into the
//        *  dist parameter.
//        */
//       virtual void initialise(VlasovDist &dist,
//                               const VelocityD &VelRange_);
// };


#include "vlasovinit_temp.cpp"

#endif
