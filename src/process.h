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
  public:
    Process(Boundary *boundary_, const SpeciesList &species_)
        : boundary(boundary_), species(species_) 
    {
      process = this;
    }
    const Boundary& getBoundary() const { return *boundary; }
    static Process& instance() { return *process; }
    void init();
    void run();
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
      Boundary *getBoundary() { return boundary; }
};

class VlasovInitRebuild;
class VlasovRebuild : public Rebuildable {
    private:
      typedef std::list<VlasovRebuild*> RebuildList;
      SpeciesData vlasovData;
      VlasovInitRebuild *initRebuild;
      static RebuildList slist;
      typedef SimpleDiagnostic<VlasovDist,HDFostream> PhaseDiag;
      PhaseDiag *phasediag;
    protected:
      virtual PARAMETERMAP* MakeParamMap (PARAMETERMAP* pm = NULL);
    public:
      VlasovRebuild();
      pVlasov getVlasovInstance();
      static SpeciesList getSpeciesList();
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
