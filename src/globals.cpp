// -*- C++ -*-
// $Id$
#include "globals.h"
#include "darwin.h"
#include "potential.h"

Parameters *Parameters::globals;
int Parameters::Argc;
char **Parameters::Argv;

Parameters::Parameters() {
  globals = this;
  Field = new FieldType();
}

PARAMETERMAP* Parameters::MakeParamMap (PARAMETERMAP* pm) {
  pm = Rebuildable::MakeParamMap(pm);
  (*pm)["x-size"] = WParameter(new ParameterValue<int>(&GridX,256));
  (*pm)["y-size"] = WParameter(new ParameterValue<int>(&GridY,32));
  (*pm)["v-size"] = WParameter(new ParameterValue<int>(&VelRes,32));
  (*pm)["v-size-x"] = WParameter(new ParameterValue<int>(&VelResX,-1));
  (*pm)["v-size-y"] = WParameter(new ParameterValue<int>(&VelResY,-1));
  (*pm)["v-size-z"] = WParameter(new ParameterValue<int>(&VelResZ,-1));
  (*pm)["x-space"] = WParameter(new ParameterValue<double>(&GridSpace_x,0.39269908));
  (*pm)["y-space"] = WParameter(new ParameterValue<double>(&GridSpace_y,0.39269908));
  (*pm)["dt"] = WParameter(new ParameterValue<double>(&Dt,5e-3));
  (*pm)["vel-ratio"] = WParameter(new ParameterValue<double>(&VelocityRatio,0.01));

  (*pm)["T"] = WParameter(new ParameterValue<int>(&TotalTime,100000));

  (*pm)["Ex"] = WParameter(new ParameterValue<double>(&EField[0],0));
  (*pm)["Ey"] = WParameter(new ParameterValue<double>(&EField[1],0));
  (*pm)["Ez"] = WParameter(new ParameterValue<double>(&EField[2],0));
  (*pm)["Bx"] = WParameter(new ParameterValue<double>(&BField[0],0));
  (*pm)["By"] = WParameter(new ParameterValue<double>(&BField[1],0));
  (*pm)["Bz"] = WParameter(new ParameterValue<double>(&BField[2],0));
  
  (*pm)["bg-density"] = WParameter(new ParameterValue<double>(&BgDensity,0));

  (*pm)["mg-nu1"] = WParameter(new ParameterValue<int>(&MGNu1,1));
  (*pm)["mg-nu2"] = WParameter(new ParameterValue<int>(&MGNu2,1));
  (*pm)["mg-gamma"] = WParameter(new ParameterValue<int>(&MGGamma,1));
  (*pm)["mg-epsilon"] = WParameter(new ParameterValue<double>(&MGEpsilon,1e-6));
  return pm;
}   

std::string Parameters::Rebuild(std::istream& in) 
{
  std::string strToken = Rebuildable::Rebuild(in);
  GridLow = PositionI(0,0);
  
  if (-1==VelResX) VelResX = VelRes;
  if (-1==VelResY) VelResY = VelRes;
  if (-1==VelResZ) VelResZ = VelRes;
  
  GridHigh = PositionI(GridX+1,GridY+1);
  DistLow = PhasePositionI(-1,-1,0,0,0);
  DistHigh = PhasePositionI(GridX+2,GridY+2,VelResX+1,VelResY+1,VelResZ+1);
  return strToken;
}

Boundary **BoundaryKeeper::ppboundary;
