// -*- C++ -*-
// $Id$

#ifndef FIELDDIAG_H
#define FIELDDIAG_H

#include "forcefield.h"
#include "diagnostic.h"

template<class Field> 
class FieldDiagnostic { 
 public:
    PARAMETERMAP* MakeParamMap (PARAMETERMAP* pm);
    void init();
};

template<>
class FieldDiagnostic<Potential> : public Rebuildable
{
    typedef SimpleDiagnostic<ScalarField,std::ofstream> DiagType;
    DiagType *Ex_diag;
    DiagType *Ey_diag;
  public:
    FieldDiagnostic();
    PARAMETERMAP* MakeParamMap(PARAMETERMAP* pm);
    void init();
};

template<>
class FieldDiagnostic<VoidPotential> : public Rebuildable
{
  public:
    PARAMETERMAP* MakeParamMap (PARAMETERMAP* pm) { return pm; }
    void init() {}
};

template<>
class FieldDiagnostic<Darwin> : public Rebuildable
{
    typedef SimpleDiagnostic<ScalarField,std::ofstream> DiagType;
    DiagType *Ex_diag;
    DiagType *Ey_diag;
    DiagType *Ez_diag;
    DiagType *Bx_diag;
    DiagType *By_diag;
    DiagType *Bz_diag;
  public:
    FieldDiagnostic();
    PARAMETERMAP* MakeParamMap (PARAMETERMAP* pm);
    void init();
};

#include "fielddiag.t"

#endif

