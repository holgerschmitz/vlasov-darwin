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
#include "potential.h"
#include <vector>

typedef VlasovSpecies<EMDarwinForce> DarwinVlasovSpecies;

typedef PtrWrapper<DarwinVlasovSpecies> WDarwinVlasovSpecies;


/** @brief Class that implements the Darwin approximation for 
 *  Maxwell's equation.
 */
class Darwin {
 protected:
  /// Lower and upper bound of the numerical grid
  PositionI LBound,HBound;
  /** @brief grid spacing */
  PositionD dx;
  /** @brief A timestep counter */
  int tstep;
  /** @brief Volume of a grid cell in \f${\rm m}^3\f$ */
  double dV;		    ///< Volume of a grid cell in \f${\rm m}^3\

    double n0;          ///< Additional density
        double csc;         ///< the ratio (c_s / c)^2

        /// Container of all the Species that contribute to the charge density
		vector<WDarwinVlasovSpecies> species;
        
        ScalarField In, Lambda, Out;            
        Poisson *pois;
        Helmholtz *helmh;
        
        /// The grid containing the potential values
		ScalarField Pot, Az;
        ScalarField Ex, Ey, Ez,  Bx, By, Bz;
        ScalarField Etx, Ety, Theta, DivEt;
        /** @brief Contains the charge and current densities 
         * \f$\rho({\bf x})\f$ and \f${\bf j}({\bf x})\f$
         */
		ScalarField den, om2, jx, jy, jz, sx, sy, sz;
        ScalarField vxx, vxy, vxz, vyy, vyz, vzz;
        
        bool mainproc;
	public:  
		/// Default constructor
		Darwin (bool mainproc_) : mainproc(mainproc_) {}
        /// Destructor
		virtual ~Darwin () {}

        /// Return the scalar field that stores \f$E_x\f$
        ScalarField &GetEx() { return Ex; }
        /// Return the scalar field that stores \f$E_y\f$
        ScalarField &GetEy() { return Ey; }
        /// Return the scalar field that stores \f$E_z\f$
        ScalarField &GetEz() { return Ez; }

        /// Return the scalar field that stores \f$B_x\f$
        ScalarField &GetBx() { return Bx; }
        /// Return the scalar field that stores \f$B_y\f$
        ScalarField &GetBy() { return By; }
        /// Return the scalar field that stores \f$B_z\f$
        ScalarField &GetBz() { return Bz; }

        /// Returns the grid size
		virtual const PositionI &GetLBound () const { return LBound; };
		virtual const PositionI &GetHBound () const { return HBound; };

        /// Get parameters Nx and dx and initialize the potential grid
		virtual void Init ();
                
        virtual bool Execute (double timestep);

        void AddSpecies(DarwinVlasovSpecies* pS);
    private:
        Darwin() {}
        Darwin(Darwin &) {}
        void clearDiv(ScalarField &Fx, ScalarField &Fy); 

}; // Darwin

#endif // POTENTIAL_H




