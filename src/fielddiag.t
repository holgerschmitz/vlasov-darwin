#include "globals.h"

FieldDiagnostic<Potential>::FieldDiagnostic()
{
  Ex_diag = NULL;
  Ey_diag = NULL;
}

PARAMETERMAP* FieldDiagnostic<Potential>::MakeParamMap (PARAMETERMAP* pm)
{
  pm = Rebuildable::MakeParamMap(pm);
  (*pm)["ex"] = WParameter(
      new ParameterRebuild<DiagType, DiagType>(&Ex_diag)
  );
  (*pm)["ey"] = WParameter(
      new ParameterRebuild<DiagType, DiagType>(&Ey_diag)
  );
}

void FieldDiagnostic<Potential>::init()
{
  if (Ex_diag)
  {
    Ex_diag->setField(
      & (Parameters::instance().getField()->GetEx())
    );
    Ex_diag->setSingleOut(true);
  }
  if (Ey_diag) 
  {
    Ey_diag->setField(
      & (Parameters::instance().getField()->GetEy())
    );
    Ey_diag->setSingleOut(true);
  }
}

FieldDiagnostic<Darwin>::FieldDiagnostic()
{
  Ex_diag = NULL;
  Ey_diag = NULL;
  Ez_diag = NULL;
  Bx_diag = NULL;
  By_diag = NULL;
  Bz_diag = NULL;
}

PARAMETERMAP* FieldDiagnostic<Darwin>::MakeParamMap (PARAMETERMAP* pm)
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
      new ParameterRebuild<DiagType, DiagType>(&Ex_diag)
  );
  (*pm)["by"] = WParameter(
      new ParameterRebuild<DiagType, DiagType>(&Ey_diag)
  );
  (*pm)["bz"] = WParameter(
      new ParameterRebuild<DiagType, DiagType>(&Ez_diag)
  );
}

void FieldDiagnostic<Darwin>::init()
{
  if (Ex_diag)
  {
    Ex_diag->setField(
      & (Parameters::instance().getField()->GetEx())
    );
    Ex_diag->setSingleOut(true);
  }
  if (Ey_diag) 
  {
    Ey_diag->setField(
      & (Parameters::instance().getField()->GetEy())
    );
    Ey_diag->setSingleOut(true);
  }
  if (Ez_diag) 
  {
    Ez_diag->setField(
      & (Parameters::instance().getField()->GetEz())
    );
    Ez_diag->setSingleOut(true);
  }
  if (Bx_diag)
  {
    Bx_diag->setField(
      & (Parameters::instance().getField()->GetBx())
    );
    Bx_diag->setSingleOut(true);
  }
  if (By_diag) 
  {
    By_diag->setField(
      & (Parameters::instance().getField()->GetBy())
    );
    By_diag->setSingleOut(true);
  }
  if (Bz_diag) 
  {
    Bz_diag->setField(
      & (Parameters::instance().getField()->GetBz())
    );
    Bz_diag->setSingleOut(true);
  }
}

