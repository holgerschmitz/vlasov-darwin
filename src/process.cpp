// -*- C++ -*-
// $Id$
#include "process.h"
#include "wrapvlasov.h"
#include "vlasovinit.h"

Process *Process::process;

void Process::init() {
  Parameters::instance().getField()->Init();
  for (
    SpeciesList::iterator it=species.begin(); 
    it!=species.end(); 
    ++it
  )
  {
    (*it)->Init();
  }
}

void Process::run() {
  int T = Parameters::instance().totalTime(); 
  ForceField::FieldType *field = Parameters::instance().getField();
  for (time=0; time<T; ++time) 
  {
    std::cout << "Cycle " << time << std::endl << flush;
    for (
      SpeciesList::iterator it=species.begin(); 
      it!=species.end(); 
      ++it
    )
    {
      (*it)->Execute();
    }
    field->Execute();
    DiagnosticManager::instance().execute();
  }
}

PARAMETERMAP* ProcessRebuild::MakeParamMap (PARAMETERMAP* pm) {
  pm = Rebuildable::MakeParamMap(pm);
  (*pm)["boundary"] 
      = WParameter(new ParameterRebuild<BoundaryRebuild, BoundaryRebuild>(&boundary));
  (*pm)["species"] 
      = WParameter(new ParameterRebuild<VlasovRebuild, VlasovRebuild>(&vlasov));
  (*pm)["general"] 
      = WParameter(new ParameterRebuild<Parameters, Parameters>(&globals));

  return pm;
}

std::string ProcessRebuild::Rebuild(std::istream& in)
{
  std::cerr << "Rebuilding Process!\n";
  Rebuildable::Rebuild(in);
  if (NULL == boundary) {
    std::cerr << "No boundary specified! Must exit!\n";
    exit(-1);
  }
  if (NULL == globals) {
    std::cerr << "No globals specified! Must exit!\n";
    exit(-1);
  }
  std::cerr << "Making Process!\n";
  makeProcess();
  std::cerr << "Making Process! done\n";
  return "";
}

void ProcessRebuild::makeProcess()
{
  Process *p = new Process(boundary->getBoundary(), VlasovRebuild::getSpeciesList());
}

BoundaryRebuild::BoundaryRebuild() {
  BoundaryKeeper::setBoundary(&boundary);
}

Boundary *BoundaryRebuild::getBoundary() { 
  if (NULL == boundary) {
    std::cerr << "No boundary type specified inside boundary! Must exit!\n";
    exit(-1);
  }
  return boundary; 
}

PARAMETERMAP* BoundaryRebuild::MakeParamMap (PARAMETERMAP* pm) {
  pm = Rebuildable::MakeParamMap(pm);
  (*pm)["single-periodic"] 
      = WParameter(new ParameterRebuild<SinglePeriodicBoundary, Boundary>(&boundary));
#ifndef SINGLE_PROCESSOR
  (*pm)["mpi-x-periodic"] 
      = WParameter(new ParameterRebuild<MPIPeriodicSplitXBoundary, Boundary>(&boundary));
  (*pm)["mpi-xy-periodic"] 
      = WParameter(new ParameterRebuild<MPIPeriodicSplitXYBoundary, Boundary>(&boundary));
#endif // single processor
  return pm;
}

VlasovRebuild::RebuildList VlasovRebuild::slist;

VlasovRebuild::VlasovRebuild() {
  slist.push_back(this);
}

pVlasov VlasovRebuild::getVlasovInstance() {
  vlasovData.init = initRebuild->getInitInstance();
  pVlasov vl = new Vlasov(vlasovData);
  for (
    VlasovDerivedDiagnostic::DiagList::iterator it=VlasovDerivedDiagnostic::diaglist.begin();
    it !=VlasovDerivedDiagnostic::diaglist.end();
    ++it
  )
  {
    (*it)->retrieveField(vl);
  }
  
  VlasovDerivedDiagnostic::diaglist.clear();
  return vl;
}

SpeciesList VlasovRebuild::getSpeciesList() {
  SpeciesList list;
  for (
    RebuildList::iterator it=slist.begin();
    it !=slist.end();
    ++it
  )
  {
    list.push_back((*it)->getVlasovInstance());
  }
  return list; 
}

PARAMETERMAP* VlasovRebuild::MakeParamMap (PARAMETERMAP* pm) {
  pm = Rebuildable::MakeParamMap(pm);
  (*pm)["mass"] = WParameter(new ParameterValue<double>(&vlasovData.mass, 1));
  (*pm)["charge"] = WParameter(new ParameterValue<double>(&vlasovData.charge, 1));
  (*pm)["vx-range"] = WParameter(new ParameterValue<double>(&vlasovData.GridRange_vx, 10));
  (*pm)["vy-range"] = WParameter(new ParameterValue<double>(&vlasovData.GridRange_vy, 10));
  (*pm)["vz-range"] = WParameter(new ParameterValue<double>(&vlasovData.GridRange_vz, 10));
  (*pm)["init"] = WParameter(
      new ParameterRebuild<VlasovInitRebuild, VlasovInitRebuild>(&initRebuild)
  );
  (*pm)["phase-space"] = WParameter(
      new ParameterRebuild<PhaseDiag, PhaseDiag>(&vlasovData.phasediag)
  );
  (*pm)["fielddiag"] = WParameter(
      new ParameterRebuild<
        VlasovDerivedDiagnostic, 
        VlasovDerivedDiagnostic
      >(&VlasovDerivedDiagnostic::fielddiag)
  );
  return pm;
}

PARAMETERMAP* VlasovInitRebuild::MakeParamMap (PARAMETERMAP* pm) {
  pm = Rebuildable::MakeParamMap(pm);
  (*pm)["maxwell"] 
      = WParameter(new ParameterRebuild<VlasovMaxwellInit, VlasovInitialiser>(&initialiser));
  (*pm)["two-maxwell"] 
      = WParameter(new ParameterRebuild<VlasovTwoMaxwellInit, VlasovInitialiser>(&initialiser));
  (*pm)["wavegen"] 
      = WParameter(new ParameterRebuild<VlasovWaveGenInit, VlasovInitialiser>(&initialiser));
  (*pm)["gausstemp"] 
      = WParameter(new ParameterRebuild<VlasovGaussTempInit, VlasovInitialiser>(&initialiser));
  (*pm)["hdfrestart"] 
      = WParameter(new ParameterRebuild<VlasovHDFInit, VlasovInitialiser>(&initialiser));
  return pm;
}
