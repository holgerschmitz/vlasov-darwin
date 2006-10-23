// -*- C++ -*-
// $Id$
#include "process.h"
#include "wrapvlasov.h"
#include "vlasovinit.h"
#include "reconnection.h"
#include "bernstein.h"
#include "magnetostatic.h"
//-----------------------------------------------------------------------------
Process *Process::process;
//-----------------------------------------------------------------------------
Process::~Process() {
  delete boundary;
}
//-----------------------------------------------------------------------------
void Process::init() {
  //init field of the parameter instance
  Parameters::instance().getField()->Init();
  //init all species
  for ( SpeciesList::iterator it=species.begin();  it!=species.end();  ++it )
  {
    (*it)->Init();
  }
}
//-----------------------------------------------------------------------------
void Process::run() {
  int T = Parameters::instance().totalTime(); 
  ForceField::FieldType *field = Parameters::instance().getField();
  //display time step and call execute for all species +field 
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
//-----------------------------------------------------------------------------
//ProcessRebuild
//build parametermap, containing boundary, species and general parameters
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
//-----------------------------------------------------------------------------
std::string ProcessRebuild::Rebuild(std::istream& in)
{
  std::cerr << "Rebuilding Process ...\n";
  Rebuildable::Rebuild(in);
  if (NULL == boundary) {
    std::cerr << "No boundary specified! Must exit!\n";
    exit(-1);
  }
  if (NULL == globals) {
    std::cerr << "No globals specified! Must exit!\n";
    exit(-1);
  }
  std::cerr << "Making Process ...\n";
  makeProcess();
  std::cerr << "Making Process: done\n";
  return "";
}
//-----------------------------------------------------------------------------
//allocate new process, get species and boundary
// p has no use other that to allocate process, which will the constructor
// p's lifetime is limited to this method, but, since in Process' constructor a static pointer is assigned,
// this does not extend to the process object itself.
void ProcessRebuild::makeProcess()
{
  Process *p = new Process(boundary->getBoundary(), VlasovRebuild::getSpeciesList());
}
//-----------------------------------------------------------------------------
BoundaryRebuild::BoundaryRebuild() {
  BoundaryKeeper::setBoundary(&boundary);
}
//-----------------------------------------------------------------------------
Boundary *BoundaryRebuild::getBoundary() { 
  if (NULL == boundary) {
    std::cerr << "No boundary type specified inside boundary! Must exit!\n";
    exit(-1);
  }
  return boundary; 
}
//-----------------------------------------------------------------------------
PARAMETERMAP* BoundaryRebuild::MakeParamMap (PARAMETERMAP* pm) {
  pm = Rebuildable::MakeParamMap(pm);
  (*pm)["single-periodic"] 
      = WParameter(new ParameterRebuild<SinglePeriodicBoundary, Boundary>(&boundary));
#ifndef SINGLE_PROCESSOR
  (*pm)["mpi-x-periodic"] 
      = WParameter(new ParameterRebuild<MPIPeriodicSplitXBoundary, Boundary>(&boundary));
  (*pm)["mpi-xy-periodic"] 
      = WParameter(new ParameterRebuild<MPIPeriodicSplitXYBoundary, Boundary>(&boundary));
  (*pm)["simple-reconnection"] 
      = WParameter(new ParameterRebuild<SimpleReconnectionBoundary, Boundary>(&boundary));
  (*pm)["gem-reconnection"] 
      = WParameter(new ParameterRebuild<GEMReconnectionBoundary, Boundary>(&boundary));
#endif // single processor
  return pm;
}
//-----------------------------------------------------------------------------
//VlasovRebuild::RebuildList VlasovRebuild::slist;
SpeciesList VlasovRebuild::spList;
//-----------------------------------------------------------------------------
VlasovRebuild::VlasovRebuild() {
//  slist.push_back(this);
}
//-----------------------------------------------------------------------------
// allocate a new vlasov object ( typedef for:VlasovSpecies<ForceField,ADVANCE_TYPE,SCHEME_TYPE>)
//add all derived diagnostics to it.
//return pointer to this vlasovSpecies object
pVlasov VlasovRebuild::getVlasovInstance() {
  //set initializer pointer in the species data object to the instance pointer of initRebuild
  vlasovData.init = initRebuild->getInitInstance();
  //allocate new VlasovSpecies 
  pVlasov vl = new Vlasov(vlasovData);
  //for all items of the derived diagnostic list:
  // add them to the new object
  //
  for (
    VlasovDerivedDiagnostic::DerivedDiagList::iterator it=VlasovDerivedDiagnostic::diaglist.begin();
    it !=VlasovDerivedDiagnostic::diaglist.end();
    ++it
  )
  { 
    //add iterator contents to the new object
    vl->addDerivedDiagnostic(*it);
    // call the retrieve method for iterator content
    //this will set the field in the derived diagnotics 
    (*it)->retrieveField(vl);
  }
  
  VlasovDerivedDiagnostic::diaglist.clear();
  return vl;
}
//-----------------------------------------------------------------------------
//add instance to species list
void VlasovRebuild::finalize() {
  spList.push_back(getVlasovInstance());
}
//-----------------------------------------------------------------------------
SpeciesList VlasovRebuild::getSpeciesList() {
  return spList;
//  SpeciesList list;
//  for (
//    RebuildList::iterator it=slist.begin();
//    it !=slist.end();
//    ++it
//  )
//  {
//    list.push_back((*it)->getVlasovInstance());
//  }
//  return list; 
}
//-----------------------------------------------------------------------------
PARAMETERMAP* VlasovRebuild::MakeParamMap (PARAMETERMAP* pm) {
  pm = Rebuildable::MakeParamMap(pm);
  (*pm)["mass"] = WParameter(new ParameterValue<double>(&vlasovData.mass, 1));
  (*pm)["charge"] = WParameter(new ParameterValue<double>(&vlasovData.charge, 1));
  (*pm)["vx-range"] = WParameter(new ParameterValue<double>(&vlasovData.GridRange_vx, 10));
  (*pm)["vy-range"] = WParameter(new ParameterValue<double>(&vlasovData.GridRange_vy, 10));
  (*pm)["vz-range"] = WParameter(new ParameterValue<double>(&vlasovData.GridRange_vz, 10));
  (*pm)["density"] = WParameter(new ParameterValue<double>(&vlasovData.densityGoal, 0));
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
//-----------------------------------------------------------------------------
PARAMETERMAP* VlasovInitRebuild::MakeParamMap (PARAMETERMAP* pm) {
  pm = Rebuildable::MakeParamMap(pm);
  (*pm)["maxwell"] 
      = WParameter(new ParameterRebuild<VlasovMaxwellInit, VlasovInitialiser>(&initialiser));
  (*pm)["two-maxwell"] 
      = WParameter(new ParameterRebuild<VlasovTwoMaxwellInit, VlasovInitialiser>(&initialiser));
  (*pm)["wavegen"] 
      = WParameter(new ParameterRebuild<VlasovWaveGenInit, VlasovInitialiser>(&initialiser));
  (*pm)["bernstein"] 
      = WParameter(new ParameterRebuild<VlasovBernsteinInit, VlasovInitialiser>(&initialiser));
  (*pm)["gausstemp"] 
      = WParameter(new ParameterRebuild<VlasovGaussTempInit, VlasovInitialiser>(&initialiser));
  (*pm)["hdfrestart"] 
      = WParameter(new ParameterRebuild<VlasovHDFInit, VlasovInitialiser>(&initialiser));
  (*pm)["reconnection"] 
      = WParameter(new ParameterRebuild<VlasovReconnectionInit, VlasovInitialiser>(&initialiser));
 (*pm)["gem-reconnection"] 
      = WParameter(new ParameterRebuild<GEMReconnectionInit, VlasovInitialiser>(&initialiser));
  (*pm)["periodic-reconnection"] 
      = WParameter(new ParameterRebuild<VlasovPeriodicReconnectionInit, VlasovInitialiser>(&initialiser));
  return pm;
}
//-----------------------------------------------------------------------------
