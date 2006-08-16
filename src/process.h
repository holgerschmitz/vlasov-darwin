// -*- C++ -*-
// $Id$

#ifndef PROCESS_H
#define PROCESS_H

#include "task.h"
#include "globals.h"
#include "simtype.h"
#include "diagnostic.h"
#include "hdfstream.h"

typedef std::list<pVlasov> SpeciesList;

class Process {
  private:
    Boundary *boundary;
    SpeciesList species;
    static Process *process;
    int time;
  public:
    Process(Boundary *boundary_, const SpeciesList &species_)
        : boundary(boundary_), species(species_) 
    {
      process = this;
      time = 0;
    }
    ~Process();
    Boundary& getBoundary() const { return *boundary; }
    int getTime() { return time; }
    static Process& instance() { return *process; }
    void init();
    void run();
    static void finalize() { delete process; }
};

class BoundaryRebuild;
class VlasovRebuild;

class ProcessRebuild : public Rebuildable {
    private:
      BoundaryRebuild *boundary;
      VlasovRebuild *vlasov;
      Parameters *globals;
    protected:
      virtual PARAMETERMAP* MakeParamMap (PARAMETERMAP* pm = NULL);
      void makeProcess();
    public:
      virtual std::string Rebuild (std::istream& in);
};

class BoundaryRebuild : public Rebuildable {
    private:
      Boundary *boundary;
    protected:
      virtual PARAMETERMAP* MakeParamMap (PARAMETERMAP* pm = NULL);
    public:
      BoundaryRebuild();
      Boundary *getBoundary();
};

class VlasovInitRebuild;
class VlasovRebuild : public Rebuildable {
    private:
//      typedef std::list<VlasovRebuild*> RebuildList;
      static SpeciesList spList;
      SpeciesData vlasovData;
      VlasovInitRebuild *initRebuild;
//      static RebuildList slist;
    protected:
      virtual PARAMETERMAP* MakeParamMap (PARAMETERMAP* pm = NULL);
    public:
      VlasovRebuild();
      pVlasov getVlasovInstance();
      static SpeciesList getSpeciesList();

      void finalize();
};

class VlasovInitRebuild : public Rebuildable {
    private:
      VlasovInitialiser *initialiser;
    protected:
      virtual PARAMETERMAP* MakeParamMap (PARAMETERMAP* pm = NULL);
    public:
      VlasovInitialiser *getInitInstance() { return initialiser; }
};

#endif
