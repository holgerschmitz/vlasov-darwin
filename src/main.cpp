
#include "Vlas2d.h"
#include "wrapvlasov.h"
#include "boundary.h"
#include <fstream>
#include <strstream>


/** @brief Print some initial diagnostics and start the simulation.
 *
 *  If an argument is specified at the command line then it
 *  is interpreted as a file name and the setup is read from the
 *  corresponding file. If no argument is specified OPAR tries to read
 *  the setup from the standard input.
 */
int main (int argc, char** argv) {

//    VlasovSpecies<EFieldForce> vlasov(new SinglePeriodicBoundary);
//    VlasovSpecies<EFieldForce> vlasov(new MPIPeriodicSplitXBoundary(argc,argv));
    VlasovSpecies<EMDarwinForce> vlasov(new MPIPeriodicSplitXBoundary(argc,argv));
        
    vlasov.Init();
    
    ofstream Exfile;   
    ofstream Eyfile;   
    ofstream Ezfile;   
    ofstream Bxfile;   
    ofstream Byfile;   
    ofstream Bzfile;   
    ofstream Jxfile;   
    ofstream Jyfile;   
    ofstream Jzfile;   
    ofstream Rhofile;
    ofstream TotalEngy;
    ofstream Phase;
    
    if (vlasov.master()) {
        Exfile.open("ExWaves.out");   
        Eyfile.open("EyWaves.out");   
        Ezfile.open("EzWaves.out");   
        Bxfile.open("BxWaves.out");   
        Byfile.open("ByWaves.out");   
        Bzfile.open("BzWaves.out");   
        Jxfile.open("jxWaves.out");   
        Jyfile.open("jyWaves.out");   
        Jzfile.open("jzWaves.out");   
        Rhofile.open("RhoWaves.out");
        TotalEngy.open("TotalEngy.out");
        Phase.open("PhaseVxVySlice.out");
    }

/*       
    PositionI pos(1,1);
    VelocityD v(1,0,0);
    VelocityD F;
    PositionD x(0,0);
    
    double dt = 1;
    ofstream path("path.out");
    
    vlasov.dttx = dt;
    vlasov.vnorm = 5;
    
    for (int t=0; t<2000; ++t) {
      F = vlasov.Force(pos,v,0.01);
      v = v+F;
      x[0] += dt*v[0];
      x[1] += dt*v[1];
      path << dt*t << " " << x[0] << " " << x[1]
           << " " << v[0] << " " << v[1] << " " << v[1]
           << " " << F[0] << " " << F[1] << " " << F[1] << "\n";
    }
    
*/    
    for (int t=0; t<300000; ++t) {
        cerr << "Cycle " << t << endl;
        vlasov.Execute();
        if ( (t%10)==0 ) {
            double En = vlasov.TotalEnergy();
            if (vlasov.master())
                TotalEngy << t << " " << En << endl;
            //vlasov.writeVxVySlice(3,3,0,"PhaseSlice.out");
            vlasov.writeVxVySlice(t,2,3,16,Phase);        
        }
        if ( (t>0) && vlasov.master()  && ((t%5) == 0)) {
//        if (false) {
             write_Scalar(vlasov.GetEx(),Exfile);
             write_Scalar(vlasov.GetEy(),Eyfile);
             write_Scalar(vlasov.GetEz(),Ezfile);
             write_Scalar(vlasov.GetBx(),Bxfile, Gl_B0x);
             write_Scalar(vlasov.GetBy(),Byfile, Gl_B0y);
             write_Scalar(vlasov.GetBz(),Bzfile, Gl_B0z);
             write_Scalar(vlasov.getJx(),Jxfile);
             write_Scalar(vlasov.getJy(),Jyfile);
             write_Scalar(vlasov.getJz(),Jzfile);
             write_Scalar(vlasov.Rho(),Rhofile);
        }
    }


//    vlasov.writeVxVySlice(2,3,16,"PhaseVxVySlice.out");
    
/*    char fname[20];
    ostrstream ost(fname,20);
    ost << "PhaseSpace" << vlasov.procnum() << ".out";
    ofstream PhaseOut(fname);
    vlasov.write(PhaseOut);
    PhaseOut.close();
*/    
    return 0;
}




