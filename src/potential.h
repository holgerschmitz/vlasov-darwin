/***************************************************************************
                          potential.h  -  description
                             -------------------
    begin                : Fri Jun 2 2000
    copyright            : (C) 2000 by Holger Schmitz
    email                : 
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
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


/** @brief A base class for potential classes
 *
 *  This base class contains the data structures for managing a potential
 *  however no algorithms for its calculation.
 *  @param geometry Name of the reference Geometry
 *  @param diagnostic Allows the output of the potantial grid
 */
class Potential {
	protected:
		PositionI LBound,HBound;    ///< Size of the grid
		PositionD dx;        ///< grid spacing

		ScalarField Pot;   ///< The grid containing the potential values

        double n0;          ///< Additional density

		ScalarField den;   ///< Contains the charge density values \f$\rho({\bf x})\f$
		vector<WESVlasovSpecies> species;
                         ///< Container of all the Species that contribute to the charge density
        ScalarField In;            
        Poisson *pois;

        ScalarField Ex,Ey;
        
        ES_EFieldEnergy DiagField;
        bool mainproc;
	public:  
		/// Default constructor
		Potential (bool mainproc_) : mainproc(mainproc_) {}
        /// Destructor
		virtual ~Potential () {}


        ScalarField &GetEx() { return Ex; }
        ScalarField &GetEy() { return Ey; }

        /// Returns the grid size
		virtual const PositionI &GetLBound () const { return LBound; };
		virtual const PositionI &GetHBound () const { return HBound; };

		virtual void Init ();
        virtual void Execute (double timestep);

        void AddSpecies(ESVlasovSpecies* pS);

}; // Potential

void write_Scalar(ScalarField &, const char*, double offset=0.0);
void write_Scalar(ScalarField &, ostream&, double offset=0.0);

#endif // POTENTIAL_H




