// -*- C++ -*-
// $Id$

#ifndef FIELDDIAG_H
#define FIELDDIAG_H

#include "scalarfield.h"
#include "forcefield.h"
#include "diagnostic.h"


class EFieldDiagnostic : public Rebuildable
{
    typedef SimpleDiagnostic<ScalarField,std::ofstream> DiagType;
    DiagType *Ex_diag;
    DiagType *Ey_diag;
    DiagType *Pot_diag;
  public:
    EFieldDiagnostic();
    PARAMETERMAP* MakeParamMap(PARAMETERMAP* pm);
    void init(ForceField::FieldType *);
};

class VoidFieldDiagnostic : public Rebuildable
{
  public:
    PARAMETERMAP* MakeParamMap (PARAMETERMAP* pm) { return pm; }
    void init(ForceField::FieldType *) {}
};

class EBFieldDiagnostic : public Rebuildable
{
    typedef SimpleDiagnostic<ScalarField,std::ofstream> DiagType;
    DiagType *Ex_diag;
    DiagType *Ey_diag;
    DiagType *Ez_diag;
    DiagType *Bx_diag;
    DiagType *By_diag;
    DiagType *Bz_diag;
  public:
    EBFieldDiagnostic();
    PARAMETERMAP* MakeParamMap (PARAMETERMAP* pm);
    void init(ForceField::FieldType *);
};

#endif

