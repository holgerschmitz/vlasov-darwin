// -*- C++ -*-
// $Id$


#include "wrapvlasov.h"
#include "boundary.h"
#include <fstream>
#include <string>
#include <signal.h>

bool VlasovDiagnostic = false;

list<pVlasov> species;

ForceField::FieldType *field;

void FPE_Exception(int)
{
  std::cerr << "Floating point exception\n";
  assert(0);
}

/** @brief Print some initial diagnostics and start the simulation.
 */
int main (int argc, char** argv) {

  Parameters::setArgc(argc);
  Parameters::setArgv(argv);

  signal(SIGFPE, FPE_Exception);
  
  // If argument is specified read setup from file
  std::string setupfilename = "setup.dat";
  if (argc > 1) setupfilename = argv[1];
  
  ifstream setupfile(setupfilename.c_str());
  cerr << "MAIN: Reading Inputfile " << setupfilename << " ...\n";
  
  ProcessRebuild().Rebuild(setupfile);
  
  cerr << "MAIN: initializing\n";
  Process::instance().init();
  
  cout << "MAIN: running\n";
  Process::instance().run();
      
  Process::finalize();
  return 0;
}




