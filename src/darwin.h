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


#ifndef DARWIN_H
#define DARWIN_H

#include "numeric.h"
#include "stlpwrapper.h"
#include "poisson.h"
#include "helmholtz.h"
#include "wrapvlasov.h"
#include "diaghelper.h"
#include <vector>

typedef VlasovSpecies<EMDarwinForce> DarwinVlasovSpecies;

typedef PtrWrapper<DarwinVlasovSpecies> WDarwinVlasovSpecies;


/** @brief Class that implements the Darwin approximation for 
 *  Maxwell's equation.
 */
class Darwin {
	protected:
		PositionI Nx;       ///< Size of the grid
		PositionD dx;       ///< grid spacing
		double dV;		    ///< Volume of a grid cell in \f${\rm m}^3\f$

        double n0;          ///< Additional density
        double csc;         ///< the ratio (c_s / c)^2

		vector<WDarwinVlasovSpecies> species;
        ///< Container of all the Species that contribute to the charge density
        
        ScalarField In, Lambda, Out;            
        Poisson *pois;
        Helmholtz *helmh;
        
        /// The grid containing the potential values
		ScalarField Pot, Az;
        ScalarField Ex, Ey, Ez,  Bx, By, Bz;
        /** @brief Contains the charge and current densities 
         * \f$\rho({\bf x})\f$ and \f${\bf j}({\bf x})\f$
         */
		ScalarField den, om2, jx, jy, jz, sx, sy, sz;
        ScalarField vxx, vxy, vxz, vyy, vyz, vzz;
        
	public:  
		/// Default constructor
		Darwin () {};
        /// Destructor
		virtual ~Darwin () {};


        ScalarField &GetEx() { return Ex; }
        ScalarField &GetEy() { return Ey; }
        ScalarField &GetEz() { return Ez; }

        ScalarField &GetBx() { return Bx; }
        ScalarField &GetBy() { return By; }
        ScalarField &GetBz() { return Bz; }

        /// Returns the grid size
		virtual PositionI GetNx () const { return Nx; };

		virtual void Init ();
        virtual bool Execute (double timestep);

        void AddSpecies(DarwinVlasovSpecies* pS);

}; // Darwin

void write_Scalar(ScalarField &, const char*);

#endif // POTENTIAL_H




