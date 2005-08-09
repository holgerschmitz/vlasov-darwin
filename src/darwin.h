// -*- C++ -*-
// $Id$

/** @file darwin.h
 *  Classes implementing potentials
 */


#ifndef DARWIN_H
#define DARWIN_H

#include "scalarfield.h"
#include "numeric.h"
#include "stlpwrapper.h"
#include "poisson.h"
#include "helmholtz.h"
#include "forcefield.h"
#include "diaghelper.h"
#include "simtype.h"
#include "forcefield.h"
#include <vector>
#include <string>


/** @brief Class that implements the Darwin approximation for 
 *  Maxwell's equation.
 */
class Darwin {
  protected:
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
      
      /** @brief The ratio of thermal speed to speed of ligh 
       *  \f$(v_th / c)^2\f$
       */
      double csc; 

      /// Container of all the Species that contribute to the charge density
      vector<pEMDarwinForce> species;

      /// Temporary fields needed for the Poisson and Helmholtz solvers
      ScalarField In, Lambda, Out;
      
      /// The poisson solver object
      Poisson *pois;
      
      /// The Helmholtz solver object
      Helmholtz *helmh;

      /** @brief The grid containing the potential values.
       *  Scalar potential and z-component of the vector potential
       */
      ScalarField Pot, Az;
      
      /// All the electric and magnetic field components
      ScalarField Ex, Ey, Ez,  Bx, By, Bz;
      
      /// The transverse electric field
      ScalarField Etx, Ety;
      
      /** @brief Helper fields for clearing the divergence of the 
       *  transverse electric field.
       */
      ScalarField Theta, DivEt;
      
      /** @brief Contains the charge and current densities 
       * \f$\rho({\bf x})\f$ and \f${\bf j}({\bf x})\f$
       */
      ScalarField den, jx, jy, jz, jxold, jyold, jzold;
      
      /** @brief Contains the charge and current densities 
       * \f$\omega^2({\bf x})=\sum_i q_i\rho_i({\bf x})\f$ 
       *  and \f${\bf s}({\bf x})=\sum_i q_i{\bf j}_i({\bf x})\f$
       */
      ScalarField om2, sx, sy, sz;
      /** @brief The components of the \f${\bf vv}\f$ tensor
       *
       */
      ScalarField vxx, vxy, vxz, vyy, vyz, vzz;
      
      bool initOldFields;
  public:  
	    /// Construct passing the value of mainproc
	    Darwin () {}

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
      void AddSpecies(EMDarwinForce* pS);
  private:
      /// Make copy constructor private
      Darwin(Darwin &) {}
      
      /// Clear the divergence of a dimensional field in two dimensions
      void clearDiv(ScalarField &Fx, ScalarField &Fy); 

}; // Darwin

#endif // POTENTIAL_H




