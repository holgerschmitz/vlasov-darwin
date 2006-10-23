// -*- C++ -*-
// $Id$
//-----------------------------------------------------------------------------
/** @file darwin.h
 *  @brief class implementing the Darwin approximation
 *
 * Implements the Darwin approximation of the Maxwell equations, in order to suppress purely electromagnetic modes,
 * which would introduce a unwanted CFL condition.
 */
//-----------------------------------------------------------------------------
#ifndef DARWIN_H
#define DARWIN_H
//-----------------------------------------------------------------------------
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
//-----------------------------------------------------------------------------
/** @brief Class that implements the Darwin approximation for 
 *  Maxwell's equation.
 */
class Darwin {
  protected:
      
	PositionI LBound, ///< Lower bound of the numerical grid
		    HBound; ///< Upper bound of the numerical grid
	 
      /// grid spacing
      PositionD dx;
      
      /// time step
      double dt;
            
      /** @brief Volume of a grid cell in \f${\rm m}^3\f$ */
      double dV;		    

      /// The background charge density
      double n0;
      
      /** @brief The ratio of thermal speed to speed of light 
       *  \f$(v_th / c)^2\f$
       */
      double csc; 

      /// Container of all the Species that contribute to the charge density
      vector<pEMDarwinForce> species;

      ScalarField In,     ///< Temporary field needed for the Poisson and Helmholtz solvers
			Lambda, ///< Temporary field needed for the Poisson and Helmholtz solvers
			Out;    ///< Temporary field needed for the Poisson and Helmholtz solvers
	  
      /// The poisson solver object
      Poisson *pois;
      
      /// The Helmholtz solver object
      Helmholtz *helmh;

	ScalarField Pot, ///< The grid containing the Scalar potential
			Az; ///< The grid containing the z-component of the vector potential
	  
      ScalarField Ex,///<  x-component of the electric field
			Ey,///<  y-component of the electric field
			Ez,///<  z-component of the electric field
			Bx,///<  x-component of the magnetic field
			By,///<  y-component of the magnetic field
			Bz; ///<  z-component of the magnetic field

      ScalarField Etx,  ///< The transverse electric field
			Ety;  ///< The transverse electric field
      
      ScalarField Theta, ///< Helper field for clearing the divergence of the transverse electric field.
			DivEt; ///<  Helper field for clearing the divergence of the transverse electric field.

      /// @brief Contains the charge  density \f$\rho({\bf x})\f$
      ScalarField den, 
			jx,    ///< x-component of the current density
			jy,    ///< y-component of the current density
			jz,    ///< z-component of the current density
			jxold, ///< x-component of the old current density
			jyold, ///< y-component of the old current density
		jzold;       ///< z-component of the old current density
      
      /// @brief Contains the charge and current densities \f$\omega^2({\bf x})=\sum_i q_i\rho_i({\bf x})\f$ 
      ScalarField om2, 
			sx, ///<x-component of the current density \f${\bf s}({\bf x})=\sum_i q_i{\bf j}_i({\bf x})\f$
			sy, ///<y-component of  the current density \f${\bf s}({\bf x})=\sum_i q_i{\bf j}_i({\bf x})\f$
			sz; ///<z-component of the current density \f${\bf s}({\bf x})=\sum_i q_i{\bf j}_i({\bf x})\f$


      ScalarField vxx, ///< xx component of the \f${\bf vv}\f$ tensor
			vxy, ///< xy component of the \f${\bf vv}\f$ tensor
			vxz, ///< xz component of the \f${\bf vv}\f$ tensor
			vyy, ///< yy component of the \f${\bf vv}\f$ tensor
			vyz, ///< yz component of the \f${\bf vv}\f$ tensor
			vzz; ///< zz component of the \f${\bf vv}\f$ tensor
 
      ///if true, init() initializes old current density with current values
      bool initOldFields;
  public:  
	/// Constructor, passing the value of mainproc
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
	  
	/// Return field component specified by string (Ex, ..., Bx, ...); Bz is returned by default
      ScalarField *GetByName(const std::string&);
	  
      /** @brief Initialize all the physical quantities and allocate
       *  the scalar fields
       */
	void Init ();

      /** @brief Calculates the electric and magnetic field in the 
       *  Darwin approximation. 
	 *
	 *  These fields can then be retrieved by
       *  the GetE* and the GetB* methods
       */
      bool Execute ();

     /** @brief Add a species to the darwin solver.
	 *
       *  The charge and current densities have to be summed over
       *  all species
       */
      void AddSpecies(EMDarwinForce* pS);
  private:
      /// Make copy constructor private
      Darwin(Darwin &) {}
      
      /// Clear the divergence of a dimensional field in two dimensions
      void clearDiv(ScalarField &Fx, ScalarField &Fy); 

}; 
// Darwin
//-----------------------------------------------------------------------------
#endif // POTENTIAL_H
