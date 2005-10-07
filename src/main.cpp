// -*- C++ -*-
// $Id$


#include "wrapvlasov.h"
#include "boundary.h" 
#include <fstream>
#include <string>
#include <unistd.h>

int myLocalFiniteDiff4SchemeRK3Counter=0;
bool VlasovDiagnostic = false;

list<pVlasov> species;

ForceField::FieldType *field;

/** @brief Print some initial diagnostics and start the simulation.
 */
int main (int argc, char** argv) {

#ifndef SINGLE_PROCESSOR
      MPI_Init(&argc, &argv);
#endif

  Parameters::setArgc(argc);
  Parameters::setArgv(argv);
 
  // If argument is specified read setup from file
  
  std::string setupfilename = "setup.dat";
//  if (argc > 1) setupfilename = argv[1];

  char buf[200];
  getcwd(buf,200);
  cerr << "CWD is " << buf << "\n";
  
  ifstream setupfile(setupfilename.c_str());
  cerr << "MAIN: Reading Inputfile " << setupfilename.c_str() << " ...\n";
  
  if (!setupfile)
  {
    cerr << "MAIN: Could not open " << setupfilename.c_str() << " ...\n"; 
    exit(-1);
  }
  ProcessRebuild().Rebuild(setupfile);
  
  cerr << "MAIN: initializing\n";
  Process::instance().init();
  
  cout << "MAIN: running\n";
  Process::instance().run();
      
  Process::finalize();
  return 0;
}




