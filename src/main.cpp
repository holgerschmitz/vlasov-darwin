// -*- C++ -*-
// $Id$


#include "wrapvlasov.h"
#include "boundary.h"
#include <fstream>
#include <string>


int argc;
char **argv;

list<pVlasov> species;

ForceField::FieldType *field;



/** @brief Print some initial diagnostics and start the simulation.
 */
int main (int argc_, char** argv_) {

  argc = argc_;
  argv = argv_;

  // If argument is specified read setup from file
  std::string setupfilename = "setup.dat";
  if (argc > 1) setupfilename = argv[1];
  
  ifstream setupfile(setupfilename.c_str());
  cout << "MAIN: Reading Inputfile " << setupfilename << " ...\n";
  
  ProcessRebuild().Rebuild(setupfile);
  
  cout << "MAIN: initializing\n";
  Process::instance().init();
  
  cout << "MAIN: running\n";
  Process::instance().run();
      
  return 0;
}




