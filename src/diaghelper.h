// -*- C++ -*-
// $Id$

/** @file diaghelper.h 
	* @brief helper classes for diagnostics 
	*
	* This file contains small helper classes to provide
 	* easy diagnostic of physical values.
 	*/
 
#include <fstream>

#ifndef DIAGHELPER_H
#define DIAGHELPER_H

//using namespace std;


class Potential;

/** @brief A class that provides diagnostic of the total electrostatic
 *  field energy
 *  @todo use a global timestep
 */
class ES_EFieldEnergy {
  private:
      /// The timestep.
      int t;
      /// pointer to the potential to be diagnosed
      Potential *pPot;
      /// The value of the total ES energy
      double Energy;
      /// The output file stream to diagnose into
      std::ofstream Output;
  public:
        /// Initialize the local values
	  void Init (Potential *pPot_);
        /// Perform diagnostic and write result to the output stream
	  void Execute (); 
};

#endif
