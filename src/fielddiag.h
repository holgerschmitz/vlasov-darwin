// -*- C++ -*-
// $Id$

#ifndef FIELDDIAG_H
#define FIELDDIAG_H

#include "scalarfield.h"
#include "forcefield.h"
#include "diagnostic.h"
/** @file
  * @brief Implements diagnostic classes for field diagnostics
  *
  * In detail these are: electrical field , electromagnetic fields 
  * and voidField,a dummy field class
  */

/** @brief simple diagnostic for electric fields, 
  *
  * inherits Rebuildable and uses the SimpleDiagnostic class
  */
class EFieldDiagnostic : public Rebuildable
{
    /** @brief local typedef 
	 *
	 *for a SimpleDiagnostic object, which writes scalar fields to the standard output
     */
    typedef SimpleDiagnostic<ScalarField,std::ofstream> DiagType;
    ///use SimpleDiagnostic for  the x-component of \f$ \vec E \f$
    DiagType *Ex_diag;
    ///use SimpleDiagnostic for the y-component of \f$ \vec E \f$
    DiagType *Ey_diag;
    ///use SimpleDiagnostic for the potential field
    DiagType *Pot_diag;
  public:
   /** @brief default constructor
     *
     * note that all diagnostic pointers (Ex_diag etc) are set to NULL.
     * They are to be assigned by the rebuild method.
     */
    EFieldDiagnostic();
    ///create parametermap, consisting of the field components
    PARAMETERMAP* MakeParamMap(PARAMETERMAP* pm);
   /** @brief initialize method
     *
     * which sets all fields in the diagnostic objects to the appropiate component
     */
    void init(ForceField::FieldType *);
};

///diagnostic class for voidField, a dummy class
class VoidFieldDiagnostic : public Rebuildable
{
  public:
    ///does nothing, but returning its parameter
    PARAMETERMAP* MakeParamMap (PARAMETERMAP* pm) { return pm; }
    ///does nothing
    void init(ForceField::FieldType *) {}
};

/** @brief diagnostic class for electromagnetic fields
  *
  * uses the SimpleDiagnostic class, inherits Rebuildable.
  */
class EBFieldDiagnostic : public Rebuildable
{
	/** @brief local typedef 
	  *
	  *for a SimpleDiagnostic object, which writes scalar fields to the standard output
	  */
    typedef SimpleDiagnostic<ScalarField,std::ofstream> DiagType;
    DiagType *Ex_diag; ///< diagnostic for the x-component of the electric field
    DiagType *Ey_diag; ///< diagnostic for the y-component of the electric field
    DiagType *Ez_diag; ///< diagnostic for the z-component of the electric field
    DiagType *Bx_diag; ///< diagnostic for the x-component of the magnetic field
    DiagType *By_diag; ///< diagnostic for the y-component of the magnetic field
    DiagType *Bz_diag; ///< diagnostic for the z-component of the magnetic field
  public:
    /** @brief default constructor
      *
      * note that all diagnostic pointers (Ex_diag etc) are set to NULL.
      * They are to be assigned by the rebuild method.
      */
    EBFieldDiagnostic();
    ///create parameter map, consisting of the field components
    PARAMETERMAP* MakeParamMap (PARAMETERMAP* pm);
    /** @brief initialize method
      *
      * which sets all fields in the diagnostic objects to the appropiate component
      */
    void init(ForceField::FieldType *);
};

#endif

