#include "globals.h"
#include "fielddiag.h"
#include "potential.h"
#include "darwin.h"

EFieldDiagnostic::EFieldDiagnostic()
{
  Ex_diag = NULL;
  Ey_diag = NULL;
}

PARAMETERMAP* EFieldDiagnostic::MakeParamMap (PARAMETERMAP* pm)
{
  pm = Rebuildable::MakeParamMap(pm);
  (*pm)["ex"] = WParameter(
      new ParameterRebuild<DiagType, DiagType>(&Ex_diag)
  );
  (*pm)["ey"] = WParameter(
      new ParameterRebuild<DiagType, DiagType>(&Ey_diag)
  );
  (*pm)["pot"] = WParameter(
      new ParameterRebuild<DiagType, DiagType>(&Pot_diag)
  );
  return pm;
}

void EFieldDiagnostic::init(ForceField::FieldType *field)
{
  if (Pot_diag)
  {
    Pot_diag->setField(field->GetByName("Pot"));
    Pot_diag->setSingleOut(true);
  }
  if (Ex_diag)
  {
    Ex_diag->setField(field->GetByName("Ex"));
    Ex_diag->setSingleOut(true);
  }
  if (Ey_diag) 
  {
    Ey_diag->setField(field->GetByName("Ey"));
    Ey_diag->setSingleOut(true);
  }
}

EBFieldDiagnostic::EBFieldDiagnostic()
{
  Ex_diag = NULL;
  Ey_diag = NULL;
  Ez_diag = NULL;
  Bx_diag = NULL;
  By_diag = NULL;
  Bz_diag = NULL;
}

PARAMETERMAP* EBFieldDiagnostic::MakeParamMap (PARAMETERMAP* pm)
{
  pm = Rebuildable::MakeParamMap(pm);
  (*pm)["ex"] = WParameter(
      new ParameterRebuild<DiagType, DiagType>(&Ex_diag)
  );
  (*pm)["ey"] = WParameter(
      new ParameterRebuild<DiagType, DiagType>(&Ey_diag)
  );
  (*pm)["ez"] = WParameter(
      new ParameterRebuild<DiagType, DiagType>(&Ez_diag)
  );
  (*pm)["bx"] = WParameter(
      new ParameterRebuild<DiagType, DiagType>(&Bx_diag)
  );
  (*pm)["by"] = WParameter(
      new ParameterRebuild<DiagType, DiagType>(&By_diag)
  );
  (*pm)["bz"] = WParameter(
      new ParameterRebuild<DiagType, DiagType>(&Bz_diag)
  );
  return pm;
}

void EBFieldDiagnostic::init(ForceField::FieldType *field)
{
  if (Ex_diag)
  {
    Ex_diag->setField(field->GetByName("Ex"));
    Ex_diag->setSingleOut(true);
  }
  if (Ey_diag) 
  {
    Ey_diag->setField(field->GetByName("Ey"));
    Ey_diag->setSingleOut(true);
  }
  if (Ez_diag) 
  {
    Ez_diag->setField(field->GetByName("Ez"));
    Ez_diag->setSingleOut(true);
  }
  if (Bx_diag)
  {
    Bx_diag->setField(field->GetByName("Bx"));
    Bx_diag->setSingleOut(true);
  }
  if (By_diag) 
  {
    By_diag->setField(field->GetByName("By"));
    By_diag->setSingleOut(true);
  }
  if (Bz_diag) 
  {
    Bz_diag->setField(field->GetByName("Bz"));
    Bz_diag->setSingleOut(true);
  }
}

