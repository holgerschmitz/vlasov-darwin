// -*- C++ -*-
// $Id$

#include "wrapvlasov.h"
#include "boundary.h" 
#include <fstream>
#include <string>
//equivalent of windows' stdlib.h 
#include <unistd.h>

//use of this is what exactly?
int myLocalFiniteDiff4SchemeRK3Counter=0;
//bool VlasovDiagnostic = false;
//pVlasov: pointer to the owning VlasovSpecies class, member of init class
list<pVlasov> species;

ForceField::FieldType *field;

// Print some initial diagnostics and start the simulation.
//---------------------------------------------------------------
//begin of main
int main (int argc, char** argv) {
	// if SINGLE_PROCESSOR is not defined, the MPI protocol is used 
	#ifndef SINGLE_PROCESSOR
    	  MPI_Init(&argc, &argv);
	#endif
	// parameters is found in globals.h, not parameter.h!
	// set  the static member of parameters
  	Parameters::setArgc(argc);
  	// set  the static member of parameters
  	Parameters::setArgv(argv);
 
  	//--------------------------------------------------------------- 
	// if specified get setup parameters from file
      // default setup file name, to be parsed if specified differently in the
	// argument
  	std::string setupfilename = "setup.dat";
	// if (argc > 1) setupfilename = argv[1];
	//buffer for current working directory
  	char buf[200];
	//get current working directory
  	getcwd(buf,200);
	cerr << "*******************************************************************************\n";
	cerr << "****************************** This is Vlasov2d3v *****************************\n";
	cerr << "*******************************************************************************\n";
  	cerr << "CWD is " << buf << "\n";
    //open setup file
  	ifstream setupfile(setupfilename.c_str());
	cerr << "MAIN: Reading Inputfile: " << setupfilename.c_str() << " ...\n";
    //if "setup.dat" does not exist, display error message
  	if (!setupfile)
  	{
		cerr << "MAIN: Could not open: " << setupfilename.c_str() << " ...\n"; 
    	exit(-1);
  	}
  	//--------------------------------------------------------------
      /*  Construct a local ProcessRebuild object, call its Rebuild method.
       *
  	 * This object is created and destroyed immediately after its method was called.
  	 * Rebuild calls makeProcess(), which dynamically allocates a new Process,
  	 * The pointer to this object ceases to exist after the makeProcess memeber is done.
  	 * But in Process' constructor a static pointer is set to the new process.
  	 * In result we get a static pointer to a dynamical allocated object.
  	 */
  	ProcessRebuild().Rebuild(setupfile);
	//---------------------------------------------------------------
	// instance returns a pointer to owning process,static data member
  	cerr << "                          *** MAIN:  INITIALIZING ***\n";
  	//turn over control to process class, initialize...
  	Process::instance().init();
  	cout << "                             *** MAIN: RUNNING ***\n";
  	//.. and run the process
  	Process::instance().run();
    //delete the allocated process
  	Process::finalize();
	cerr << "*************************** Vlasov2d3v has finished ***************************\n";
  return 0;
}

// end of main
