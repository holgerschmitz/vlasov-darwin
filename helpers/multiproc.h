#include "vlasov.h"

#ifndef MULTIPROC_H
#define MULTIPROC_H
class PhaseSplit
{
  public:
    virtual void setNProc(int n);
    int getNProc() { return nproc; }
    virtual void setSize(int x, int y, int vx, int vy, int vz);
    /// Returns the global lower bound of the distribution function
    virtual PhasePositionI DistLow(int procnum) = 0;
    /// Returns the global upper bound of the distribution function
    virtual PhasePositionI DistHigh(int procnum) = 0;
    const PhasePositionI &getGlobalLow() { return globalLow; }
    const PhasePositionI &getGlobalHigh() { return globalHigh; }
  protected:
    int nproc;
    PhasePositionI globalLow;
    PhasePositionI globalHigh;
};

class PhaseSplitX : public PhaseSplit
{
  public:
    PhaseSplitX(int nproc);
    /// Returns the global lower bound of the distribution function
    PhasePositionI DistLow(int procnum);
    /// Returns the global upper bound of the distribution function
    PhasePositionI DistHigh(int procnum);
};

class PhaseSplitXY : public PhaseSplit
{
  public:
    PhaseSplitXY(int nproc);
    /// Returns the global lower bound of the distribution function
    PhasePositionI DistLow(int procnum);
    /// Returns the global upper bound of the distribution function
    PhasePositionI DistHigh(int procnum);
  private:
    void makeDims();
    void makeMyCoord(int procnum);
    int dims[2];
    int mycoord[2];
};


#endif
