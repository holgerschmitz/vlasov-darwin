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
        PositionI Nx;       ///< Size of the grid
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
class VlasovWaveGenInit : public VlasovInitialiser {
    protected:
        double N;           ///< density in particles\f$m^{-{\rm dim}}\f$
        VelocityD v_th;      ///< thermal velocity
        PositionI Nx;       ///< Size of the grid
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
