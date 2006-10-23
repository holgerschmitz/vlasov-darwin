// -*- C++ -*-
// $Id$
//-----------------------------------------------------------------------------
#ifndef PROCESS_H
#define PROCESS_H
//-----------------------------------------------------------------------------
#include "task.h"
#include "globals.h"
#include "simtype.h"
#include "diagnostic.h"
#include "hdfstream.h"
//-----------------------------------------------------------------------------
/** @file process.h
  * @brief contains the process classes
  *
  * Here the process classes are defined, which are the global objects of this code.
  */
//-----------------------------------------------------------------------------
//pVlasov: pointer to the owning VlasovSpecies class, member of init class, 
//VlasovSpecies is defined in wrapvlasov.h
///list of pointer to their owning VlasovSpecies objects
typedef std::list<pVlasov> SpeciesList;
//-----------------------------------------------------------------------------
// singleton implemetation seems useful here...
// since there can be only one process, and it would make things easier regarding makeprocess etc
// and process IS already almost a singleton
/** @brief The process class
  *
  * This class has only one instance and controls setup, rebuilding, initializing and
  * executing of all objects of the programm. First it calls all the Rebuild methods, 
  * builds the correct objects and then distributes the execution to the execxte() members
  * of the appropiate classes. It also prints global (error) message on the screen.
  */
class Process {
  private:
    /// boundary conditions of the process
    Boundary *boundary;
    /// processed species
    SpeciesList species;
    /// pointer to process
    static Process *process;
    /// time step
    int time;
  public:
	//-----------------------------------------------------------------------
	//constructor & destructor
    /** @brief constructor, 
      *
      * sets static member process to "this" and time step to 0. Since process is
      * static, the lifetime of the dynamical allocated Process object is prolonged 
      * indefinetly. Which means, it has to be deleted in finalize()
      */
    Process(Boundary *boundary_, const SpeciesList &species_)
        : boundary(boundary_), species(species_) 
    {
      process = this;
      time = 0;
    }
	///destructor, deletes boundary
    ~Process();
	//-----------------------------------------------------------
	/// returns Boundary
    const Boundary& getBoundary() const { return *boundary; }
	/// returns time
    int getTime() { return time; }
	/// returns reference to current instance of the process class
    static Process& instance() { return *process; }
	/** @brief Initialize the fields and species
	  *
        * Call the initialize members of each items in the SpeciesList species and 
	  * the field from parameters.
	  */
    void init();
	/** @brief Run the process.
	  *
	  * Which means: For the desired number of timesteps call the execute() members of
	  * the objects stored in the speciesList
	  */
    void run();
	/// deletes pointer to process
    static void finalize() { delete process; }
};
//
//-----------------------------------------------------------------------------
//forward declarations
class BoundaryRebuild;
class VlasovRebuild;
//-----------------------------------------------------------------------------
//ProcessRebuild
// uses the other rebuildX classes to rebuild the process
///class for rebuilding the processs, inherits Rebuidable
class ProcessRebuild : public Rebuildable {
    private:
      ///boundary
      BoundaryRebuild *boundary;
      ///pointer to vlasovRebuild
      VlasovRebuild *vlasov;
      ///parameters
      Parameters *globals;
    protected:
      ///build parametermap, containing boundary, species and general parameters
      virtual PARAMETERMAP* MakeParamMap (PARAMETERMAP* pm = NULL);
      /** @brief allocate new process, get species and boundary.
        * 
        * Set pointer to process, which is destroyed after this member is done
        *.But when allocating the process constructor is called,
        * which will set a static pointer to the allocated object.
        * Thus we get a static pointer to a dynamically allocated object.
        * Will be used in main() to get a new process started.
        */
      void makeProcess();
    public:
      ///rebuild from file, call Rebuildable::Rebuild and makeProcess
      virtual std::string Rebuild (std::istream& in);
};
//ProcessRebuild
//------------------------------------------------------------------------------
///class for rebuilding the boundary, inherits Rebuidable
class BoundaryRebuild : public Rebuildable {
    private:
	    ///pointer to boundary
      Boundary *boundary;
    protected:
	    ///buidl the parameter map
      virtual PARAMETERMAP* MakeParamMap (PARAMETERMAP* pm = NULL);
    public:
      ///constructor
      BoundaryRebuild();
	///accessor method, return boundary
      Boundary *getBoundary();
};
//BoundaryRebuild
//-----------------------------------------------------------------------------
//forward declaration
class VlasovInitRebuild;
//----------------------------------------------------------------------------
//VlasovRebuild
///class for rebuilding the vlasov object, inherits Rebuidable
class VlasovRebuild : public Rebuildable {
    private:
//typedef std::list<VlasovRebuild*> RebuildList;
	    ///Species List
      static SpeciesList spList;
      ///characteristic data of species
      SpeciesData vlasovData;
      ///pointer to  VlasovInitRebuild
      VlasovInitRebuild *initRebuild;
//static RebuildList slist;
    protected:
      ///build parameter map
      virtual PARAMETERMAP* MakeParamMap (PARAMETERMAP* pm = NULL);
    public:
      ///constructor, does nothing
      VlasovRebuild();
      /// return pointer to VlasovSpecies
      pVlasov getVlasovInstance();
      ///return static member SpeciesList
      static SpeciesList getSpeciesList();
      //add instance to global specieslist sp of this class 
      /// Add this instance to species list.
        void finalize();
};

///class for rebuilding the initialiser, inherits Rebuidable
class VlasovInitRebuild : public Rebuildable {
    private:
	///pointer to the initialiser
      VlasovInitialiser *initialiser;
    protected:
      ///Create the parameter map
      virtual PARAMETERMAP* MakeParamMap (PARAMETERMAP* pm = NULL);
    public:
	///returns the pointer to the initialisier
      VlasovInitialiser *getInitInstance() { return initialiser; }
};
//VlasovRebuild
//-----------------------------------------------------------------------------
#endif //PROCESS_H
