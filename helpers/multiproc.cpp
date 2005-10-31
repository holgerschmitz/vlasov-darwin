#include "multiproc.h"



void PhaseSplit::setNProc(int n)
{
  nproc = n;
}

void PhaseSplit::setSize(int x, int y, int vx, int vy, int vz)
{
  globalLow = PhasePositionI(-1,-1,0,0,0);
  globalHigh = PhasePositionI(x+2,y+2,vx+1,vy+1,vz+1);
}

PhaseSplitX::PhaseSplitX(int nproc_)
{
  this->setNProc(nproc_);
}

PhasePositionI PhaseSplitX::DistLow(int procnum)
{
  PhasePositionI Low = globalLow;
  double width = (globalHigh[0]-2.)/double(nproc);
  if (procnum>0) 
    Low[0] = int(width*procnum)-1;
    
  return Low;
}

PhasePositionI PhaseSplitX::DistHigh(int procnum)
{
  PhasePositionI High = globalHigh;
  double width = (globalHigh[0]-2.)/double(nproc);
  if (procnum<nproc-1) 
    High[0] = int(width*(procnum+1))+2;
  return High;
}

PhaseSplitXY::PhaseSplitXY(int nproc_)
{
  this->setNProc(nproc_);
}

PhasePositionI PhaseSplitXY::DistLow(int procnum)
{
  makeDims();
  makeMyCoord(procnum);

  PhasePositionI Low = globalLow;
  double width[2];
  width[0] = (globalHigh[0]-2.)/double(dims[0]);
  width[1] = (globalHigh[1]-2.)/double(dims[1]);

  if (mycoord[0]>0) 
      Low[0] = int(width[0]*mycoord[0])-1;

  if (mycoord[1]>0) 
      Low[1] = int(width[1]*mycoord[1])-1;

  return Low;
}

PhasePositionI PhaseSplitXY::DistHigh(int procnum)
{
  makeDims();
  makeMyCoord(procnum);

  PhasePositionI High = globalHigh;
  double width[2];
  width[0] = (globalHigh[0]-2.)/double(dims[0]);
  width[1] = (globalHigh[1]-2.)/double(dims[1]);

  if (mycoord[0]<(dims[0]-1))
      High[0] = int(width[0]*(mycoord[0]+1))+2;

  if (mycoord[1]<(dims[1]-1))
      High[1] = int(width[1]*(mycoord[1]+1))+2;

  return High;
}

void PhaseSplitXY::makeDims()
{
  double Lx = globalHigh[0]-globalLow[0]-3;
  double Ly = globalHigh[1]-globalLow[1]-3;
  dims[0] = int(sqrt( Lx*(nproc+0.001) / Ly ));
  if (dims[0]<1) dims[0]=1;

  dims[1] = int(nproc/dims[0]);
  
  if (dims[1]<1)
  {
    dims[1]=1;
    dims[0]=nproc;
  }
}


void PhaseSplitXY::makeMyCoord(int procnum)
{
  mycoord[0] = procnum/dims[1];
  mycoord[1] = procnum - dims[1]*mycoord[0];
}
