// -*- C++ -*-
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
class VlasovMaxwellInit : public VlasovInitialiser {
 protected:
  double N;           ///< density in particles\f$m^{-{\rm dim}}\f$
    VelocityD u_stream;  ///< streaming vlocity
      VelocityD v_th;      ///< thermal velocity
        PositionI k_perturb;
        double n_perturb;
        
        // pointer to the owning VlasovSpecies class
        VlasovSpecies<ForceField> *pVlasov; 
    public:
		/// Default constructor
        VlasovMaxwellInit(VlasovSpecies<ForceField> *);
		/// Destructor
        virtual ~VlasovMaxwellInit();

        virtual void initialise(VlasovDist &dist,
                                const VelocityD &VelRange_);
};

template<class ForceField>
class VlasovTwoMaxwellInit : public VlasovInitialiser {
    protected:
        double N1;           ///< density in particles\f$m^{-{\rm dim}}\f$
        double N2;           ///< density in particles\f$m^{-{\rm dim}}\f$
        VelocityD u_stream1;  ///< streaming vlocity
        VelocityD u_stream2;  ///< streaming vlocity
        VelocityD v_th1;      ///< thermal velocity
        VelocityD v_th2;      ///< thermal velocity
        PositionI k_perturb1;
        PositionI k_perturb2;
        double n_perturb1;
        double n_perturb2;
        
        // pointer to the owning VlasovSpecies class
        VlasovSpecies<ForceField> *pVlasov; 
    public:
		/// Default constructor
        VlasovTwoMaxwellInit(VlasovSpecies<ForceField> *);
		/// Destructor
        virtual ~VlasovTwoMaxwellInit();

        virtual void initialise(VlasovDist &dist,
                                const VelocityD &VelRange_);
};

template<class ForceField>
class VlasovWaveGenInit : public VlasovInitialiser {
    protected:
        double N;           ///< density in particles\f$m^{-{\rm dim}}\f$
        VelocityD v_th;      ///< thermal velocity
        VelocityD u_stream;  ///< streaming vlocity
        double n_perturb;
        
        // pointer to the owning VlasovSpecies class
        VlasovSpecies<ForceField> *pVlasov; 
    public:
		/// Default constructor
        VlasovWaveGenInit(VlasovSpecies<ForceField> *);
		/// Destructor
        virtual ~VlasovWaveGenInit();

        virtual void initialise(VlasovDist &dist,
                                const VelocityD &VelRange_);
};

#include "vlasovinit_temp.cpp"

#endif
