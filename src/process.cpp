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
  for (int t=0; t<T; ++t) 
  {
    std::cout << "Cycle " << t << std::endl << flush;
    for (
      SpeciesList::iterator it=species.begin(); 
      it!=species.end(); 
      ++it
    )
    {
      (*it)->Execute();
    }
    field->Execute();
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
  Rebuildable::Rebuild(in);
  makeProcess();
}

void ProcessRebuild::makeProcess()
{
  Process *p = new Process(boundary->getBoundary(), VlasovRebuild::getSpeciesList());
}

BoundaryRebuild::BoundaryRebuild() {
  BoundaryKeeper::setBoundary(&boundary);
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
  return new Vlasov(vlasovData);
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
      new ParameterRebuild<PhaseDiag, PhaseDiag>(&phasediag)
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
  return pm;
}
