// -*- C++ -*-
// $Id$

/** @file potential.h
 *  Classes implementing potentials
 */


#ifndef POTENTIAL_H
#define POTENTIAL_H

#define ID_POTENTIAL "Potential"

#include "numeric.h"
#include "stlpwrapper.h"
#include "poisson.h"
#include "wrapvlasov.h"
#include "diaghelper.h"
#include <vector>

typedef VlasovSpecies<EFieldForce> ESVlasovSpecies;

typedef PtrWrapper<ESVlasovSpecies> WESVlasovSpecies;


/** @brief The class for calculating the electrostatic potential.
 *  At the moment only completely periodic boundary conditions 
 *  have been implemented.
 */
class Potential {
	protected:
        /// Size of the grid
		PositionI LBound,HBound;
        /// grid spacing
		PositionD dx;
        
        /// The grid containing the potential values
		ScalarField Pot;   

        /// Additional density
        double n0;          

        /// Contains the charge density values \f$\rho({\bf x})\f$
		ScalarField den;
        
        /// Container of all the Species that contribute to the charge density   
		vector<WESVlasovSpecies> species;
           
        /// A temporary field used for the poisson solver              
        ScalarField In;
        
        /// The Poisson solver
        Poisson *pois;

        /** @brief The electrostatic fields calculated from the derivative 
         *  of the potential
         */
        ScalarField Ex,Ey;
        
        /// A diagnostic for calculating the electrostatic field energy
        ES_EFieldEnergy DiagField;
        
        /** @brief Boolean indicating if this is the main process
         *  Used for diagnostic only
         */
        bool mainproc;
	public:  
		/// Default constructor initializing the main process variable
		Potential (bool mainproc_) : mainproc(mainproc_) {}
        /// Destructor
		virtual ~Potential () {}

        /// Returns a reference to the electric x-field
        ScalarField &GetEx() { return Ex; }
        /// Returns a reference to the electric y-field
        ScalarField &GetEy() { return Ey; }

        /// Returns the grid size (lower bound)
		virtual const PositionI &GetLBound () const { return LBound; };
        /// Returns the grid size (upper bound)
		virtual const PositionI &GetHBound () const { return HBound; };

        /** @brief Perform initialization. Setting all the field sizes
         *  and clearing the values
         */
		virtual void Init ();
        /** @brief Solve for the potential.
         *  Sums up all the charge densities from the species and
         *  solves Poissons equation to get the potential.
         *  Then it uses central differences to calculate the
         *  electric fields
         */
        virtual void Execute (double timestep);

        /** @brief Adds a species to the list of species that need 
         *  contribute to the charge density
         */
        void AddSpecies(ESVlasovSpecies* pS);

}; // Potential

/// Helper function to write a scalar potential into a file
void write_Scalar(ScalarField &, const char*, double offset=0.0);
/// Helper function to write a scalar potential into a stream
void write_Scalar(ScalarField &, ostream&, double offset=0.0);

#endif // POTENTIAL_H




