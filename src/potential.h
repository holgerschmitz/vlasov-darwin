// -*- C++ -*-
// $Id$

/** @file potential.h
 *  Classes implementing potentials
 */


#ifndef POTENTIAL_H
#define POTENTIAL_H

#define ID_POTENTIAL "Potential"

#include "numeric.h"
#include "scalarfield.h"
#include "stlpwrapper.h"
#include "poisson.h"
#include "wrapvlasov.h"
#include "diaghelper.h"
#include "forcefield.h"
#include <vector>
#include <string>

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
        
      /** @brief Volume of a grid cell in \f${\rm m}^3\f$ */
      double dV;       

      /// The grid containing the potential values
  		ScalarField Pot;   

      /// Additional density
      double n0;          

      /// Contains the charge density values \f$\rho({\bf x})\f$
		  ScalarField den;
        
      /// Container of all the Species that contribute to the charge density   
  		vector<pEFieldForce> species;
           
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
        
	public:  
		  /// Default constructor initializing the main process variable
  		Potential () {}
      /// Destructor
		  virtual ~Potential () {}

      /// Returns a reference to the electric x-field
      ScalarField &GetEx() { return Ex; }
      /// Returns a reference to the electric y-field
      ScalarField &GetEy() { return Ey; }

      ScalarField *GetByName(const std::string &);

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
      virtual void Execute ();

      /** @brief Adds a species to the list of species that need 
       *  contribute to the charge density
       */
      void AddSpecies(EFieldForce* pS);
}; // Potential

/// Helper function to write a scalar potential into a file
void write_Scalar(ScalarField &, const char*, double offset=0.0);
/// Helper function to write a scalar potential into a stream
void write_Scalar(ScalarField &, ostream&, double offset=0.0);

#endif // POTENTIAL_H




