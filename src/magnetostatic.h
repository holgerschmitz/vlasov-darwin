// -*- C++ -*-
// $Id$

/** @file magnetostatic.h
 *  @brief class for magnetostatic approximation
 *
 * Implements the magnetostatic approximation of the Maxwell equations.
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
  	  /// true if its the first step
      bool firststep;
    
      /// Lower bound of the numerical grid
      PositionI LBound,HBound;
  	  ///< upper bound of the numerical grid
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

      
      ScalarField In, 	///< Temporary field needed for the Poisson and Helmholtz solvers
      			  Lambda, ///< Temporary field needed for the Poisson and Helmholtz solvers
      			  Out; ///< Temporary field needed for the Poisson and Helmholtz solvers
      
      /// The poisson solver object
      Poisson *pois;

      /// The Helmholtz solver object
      Helmholtz *helmh;
      

       
      ScalarField Pot, ///< scalar potential
      			  Ax, ///< x-component of the vector potential
      			  Ay, ///< y-component of the vector potential
      			  Az; ///< z-component of the vector potential
      ScalarField oldjx, ///< x-component of the old current density
      			  oldjy, ///< y-component of the old current density
      			  oldjz;///< z-component of the old current density
      
 
      ScalarField Ex, ///< electric field, x-component
      			  Ey, ///< electric field, x-component
      			  Ez,  ///< electric field, z-component
      			  Bx, ///< magnetic field, x-component
      			  By, ///< magnetic field, y-component
      			  Bz;///< magnetic field, z-component
      
      /** @brief Helper field for clearing the divergence of the 
       *  transverse electric field.
       */
      ScalarField Theta, DivF;
      /**< @brief Helper field for clearing the divergence of the 
       *  transverse electric field.
       */
       
      ScalarField den, ///< charge density
      			  jx,  ///< x-component of the current density
      			  jy, ///< y-component of the current density
      			  jz;///< z-component of the current density

      ScalarField om2, ///< charge density \f$\omega^2({\bf x})=\sum_i q_i\rho_i({\bf x})\f$ 
      			  sx,///< x-component of the current density \f${\bf s}({\bf x})=\sum_i q_i{\bf j}_i({\bf x})\f$
      			  sy, ///< y-component of the current density \f${\bf s}({\bf x})=\sum_i q_i{\bf j}_i({\bf x})\f$
      			  sz;///< z-component of the current density \f${\bf s}({\bf x})=\sum_i q_i{\bf j}_i({\bf x})\f$
      
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
	  ///get field by name (Ex, ..., Bx,...) default is Bz
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




