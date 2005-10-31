#include "filename.h"
#include "multiproc.h"
#include "hdfstream.h"

#include "matrix.h"
#include <list>

void redist
  (
    SimFileName srcname, 
    SimFileName destname, 
    int NProcSrc, 
    PhaseSplit &phsplit
  )
{
  Matrix<double,5> srcphase;
  Matrix<double,5> destphase;
  
  for (int k=0; k<phsplit.getNProc(); ++k)
  {
    destname.setReplace("#p", k);
    PhasePositionI DestLow = phsplit.DistLow(k);
    PhasePositionI DestHigh = phsplit.DistHigh(k);
    
    destphase.resize(DestLow, DestHigh);
    
    std::cout << "Process " << k << ":  ("
      << DestLow[0] << ", " << DestLow[1] << ") -- ("
      << DestHigh[0] << ", " << DestHigh[1] << ")" << std::endl;
    
    // Setting flags for checking
    for (int x = DestLow[0]; x<=DestHigh[0]; ++x)
      for (int y = DestLow[1]; y<=DestHigh[1]; ++y)
        destphase(x,y,DestHigh[2],DestHigh[3],DestHigh[4]) = -1;
        
    HDFistream hdfin;
    for (int i=0; i<NProcSrc; ++i)
    {
      srcname.setReplace("#p",i);
      hdfin.open(srcname.filename().c_str());
      std::cout << "Reading Phasespace " << i << std::endl;
      hdfin >> srcphase;
      hdfin.close();
      const int *lowS = srcphase.getLow();
      const int *highS = srcphase.getHigh();

      int low[5], high[5];
      for (int s=0; s<5; ++s)
      {
        low[s] = lowS[s];
        high[s] = highS[s];
      }
      
      if (low[0]>-1) low[0] += 2;
      if (low[1]>-1) low[1] += 2;

      if (high[0] < phsplit.getGlobalHigh()[0] ) high[0] -= 2;
      if (high[1] < phsplit.getGlobalHigh()[0] ) high[0] -= 2;
      
      for (int x = low[0]; x<=high[0]; ++x)
      {
        for (int y = low[1]; y<=high[1]; ++y)
        {
        
          if ( (x>=DestLow[0]) && (x<=DestHigh[0]) && 
               (y>=DestLow[1]) && (y<=DestHigh[1]) )
          {
            for (int vx=low[2]; vx<=high[2]; vx++)
              for (int vy=low[3]; vy<=high[3]; vy++)
                for (int vz=low[4]; vz<=high[4]; vz++)
                  destphase(x,y,vx, vy, vz) = srcphase(x, y, vx, vy, vz);
          }

        }
      }
    }
    
    // Checking
    for (int x = DestLow[0]; x<=DestHigh[0]; ++x)
      for (int y = DestLow[1]; y<=DestHigh[1]; ++y)
      {
        if (destphase(x,y,DestHigh[2],DestHigh[3],DestHigh[4]) == -1)
          std::cerr << "Error at " << x << ", " << y << std::endl;
      }

    HDFostream hdfout(destname.filename().c_str());
    hdfout << destphase;
    hdfout.close();
  }
}

int main()
{
  SimFileName srcname("electron#t-#p.out");
  SimFileName destname("redist/electron#t-#p.out");
  
  int NProcSrc = 32;
  int NProcDest = 18;
  
  PhaseSplitXY phsplit(18);
  phsplit.setSize(512, 256, 20, 20, 20);
  
  int Time = 2200;  
  srcname.setReplace("#t", Time);
  destname.setReplace("#t", Time);
  redist(srcname, destname, NProcSrc, phsplit);

}
