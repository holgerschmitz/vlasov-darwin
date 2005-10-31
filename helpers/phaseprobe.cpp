
#include "filename.h"
#include "hdfstream.h"

#include "matrix.h"
#include <list>

typedef std::pair<int,int> Point;
typedef std::list<Point> PointList;

void phaseprobe(SimFileName srcname, SimFileName destname, const PointList &vec)
{
  Matrix<double,3> probe;
  Matrix<double,5> phase;
  HDFistream hdfin;
  for (int i=0; i<32; ++i)
  {
    srcname.setReplace("#p",i);
    hdfin.open(srcname.filename().c_str());
    std::cout << "Reading Phasespace " << i << std::endl;
    hdfin >> phase;
    hdfin.close();
    const int *low = phase.getLow();
    const int *high = phase.getHigh();
    
    for (PointList::const_iterator it=vec.begin(); it!=vec.end(); ++it)
    {
      int posx = it->first;
      int posy = it->second;
      if ( (posx>low[0]+1) && (posx<high[0]-1) && 
        (posy>low[1]+1) && (posy<high[1]-1) )
      {
        destname.setReplace("#x", posx);
        destname.setReplace("#y", posy);
        
        probe.resize( &(low[2]), &(high[2]) );
        
        for (int vx=low[2]; vx<=high[2]; vx++)
          for (int vy=low[3]; vy<=high[3]; vy++)
            for (int vz=low[4]; vz<=high[4]; vz++)
              probe(vx, vy, vz) = phase(posx, posy, vx, vy, vz);

        HDFostream hdfout(destname.filename().c_str());
        hdfout << probe;
        hdfout.close();
      }
    }
  }
}

int main()
{
  SimFileName srcname("electron#t-#p.out");
  SimFileName destname("probe/electron_#t_#xx#y.out");
  PointList points;
  for (int i=-20; i<=20; ++i)
  {
    int posx = int( 7*exp( 0.3*double(i) ))+1;
    int posy = int( 7*exp(-0.3*double(i) ))+1;
    points.push_back(Point(posx,posy));
  }
  
  for (int time=0; time<=1600; time+=200)
  {
    srcname.setReplace("#t", time);
    destname.setReplace("#t", time);
    phaseprobe(srcname, destname, points);
  }
}
