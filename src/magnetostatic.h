// -*- C++ -*-
// $Id$

/** @file magnetostatic.h
 *  Classes implementing potentials
 */


#ifndef MAGNETOSTATIC_H
#define MAGNETOSTATIC_H

#include "scalarfield.h"
#include "numeric.h"
#include "stlpwrapper.h"
#include "poisson.h"
#include "helmholtz.h"
#include "forcefield.h"
#include "diaghelper.h"
#include "forcefield.h"
#include <vector>
#include <string>


/** @brief Class that implements the Magnetostatic approximation for 
 *  Maxwell's equation.
 */
class Magnetostatic {
  protected:
      bool firststep;
    
      /// Lower and upper bound of the numerical grid
      PositionI LBound,HBound;
  
      /// grid spacing
      PositionD dx;
      
      /// time step
      double dt;
            
      /** @brief Volume of a grid cell in \f${\rm m}^3\f$ */
      double dV;		    

      /// The background charge density
      double n0;
      
      /// Container of all the Species that contribute to the charge density
      vector<pMagnetostaticForce> species;

      /// Temporary fields needed for the Poisson and Helmholtz solvers
      ScalarField In, Lambda, Out;
      
      /// The poisson solver object
      Poisson *pois;

      /** @brief The grid containing the potential values.
       *  Scalar potential and z-component of the vector potential
       */
      ScalarField Pot, Ax, Ay, Az;
      ScalarField oldAx, oldAy, oldAz;
      
      /// All the electric and magnetic field components
      ScalarField Ex, Ey, Ez,  Bx, By, Bz;
      
      /** @brief Helper fields for clearing the divergence of the 
       *  transverse electric field.
       */
      ScalarField Theta, DivF;
      
      /** @brief Contains the charge and current densities 
       * \f$\rho({\bf x})\f$ and \f${\bf j}({\bf x})\f$
       */
      ScalarField den, jx, jy, jz;
      
  public:  
	    /// Construct passing the value of mainproc
	    Magnetostatic () {}

      /// Return the scalar field that stores the field component \f$E_x\f$
      ScalarField &GetEx() { return Ex; }
      /// Return the scalar field that stores the field component \f$E_y\f$
      ScalarField &GetEy() { return Ey; }
      /// Return the scalar field that stores the field component \f$E_z\f$
      ScalarField &GetEz() { return Ez; }

      /// Return the scalar field that stores the field component \f$B_x\f$
      ScalarField &GetBx() { return Bx; }
      /// Return the scalar field that stores the field component \f$B_y\f$
      ScalarField &GetBy() { return By; }
      /// Return the scalar field that stores the field component \f$B_z\f$
      ScalarField &GetBz() { return Bz; }

      ScalarField *GetByName(const std::string&);
      /** @brief Initialize all the physical quantities and allocate
       *  the scalar fields
       */
	    void Init ();

      /** @brief Calculates the electric and magnetic field in the 
       *  Darwin approximation. These fields can then be retrieved by
       *  the GetE* and the GetB* methods
       */
      bool Execute ();

      /** @brief Add a species to the darwin solver.
       *  The charge and current densities have to be summed over
       *  all species
       */
      void AddSpecies(MagnetostaticForce* pS);
  private:
      /// Make copy constructor private
      Magnetostatic(Magnetostatic &) {}
      
      /// Clear the divergence of a dimensional field in two dimensions
      void clearDiv(ScalarField &Fx, ScalarField &Fy); 

}; // Magnetostatic

#endif // MAGNETOSTATIC_H




