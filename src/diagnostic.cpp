// -*- C++ -*-
// $Id$

#include "diagnostic.h"

DiagnosticInterface::DiagnosticInterface() {
  DiagnosticManager::instance().addDiagnostic(this);
  t=0;
}

void DiagnosticInterface::execute() {
  if ((0==t) && appending()) open(fname);
  ++t;
  if ( (t % interval) != 0 ) return;
  
  if (!appending()) open(fname);
  write();
  if (!appending()) close();
}

PARAMETERMAP* DiagnosticInterface::MakeParamMap(PARAMETERMAP* pm)
{
  pm = Rebuildable::MakeParamMap(pm);
  (*pm)["file"] = WParameter(
    new ParameterValue<std::string>(&fname, "")
  );
  (*pm)["append"] = WParameter(
    new ParameterValue<std::string>(&append, "n")
  );
  (*pm)["interval"] = WParameter(
    new ParameterValue<int>(&interval, 100)
  );
  return pm;
}

bool DiagnosticInterface::appending() {
  return 'y' == append[0];
}


DiagnosticManager *DiagnosticManager::theManager = NULL;

DiagnosticManager::DiagnosticManager() {}

void DiagnosticManager::addDiagnostic(DiagnosticInterface *diag)
{
  diags.push_back(pDiagnosticInterface(diag));
}

void DiagnosticManager::execute()
{
  for (
    DiagList::iterator it=diags.begin();
    it != diags.end();
    ++it
  )
  {
    (*it)->execute();
  }
}

DiagnosticManager& DiagnosticManager::instance()
{
  if (theManager==NULL) theManager = new DiagnosticManager;
  return *theManager;
}


