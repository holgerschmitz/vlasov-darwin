
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
    ofstream Jyfile;   
    ofstream Jzfile;   
    ofstream Rhofile;
    ofstream TotalEngy;
    
    if (vlasov.master()) {
        Exfile.open("ExWaves.out");   
        Eyfile.open("EyWaves.out");   
        Ezfile.open("EzWaves.out");   
        Bxfile.open("BxWaves.out");   
        Byfile.open("ByWaves.out");   
        Bzfile.open("BzWaves.out");   
        Jyfile.open("jyWaves.out");   
        Jzfile.open("jzWaves.out");   
        Rhofile.open("RhoWaves.out");
        TotalEngy.open("TotalEngy.out");
    }
    
    
    for (int t=0; t<4100; ++t) {
        cerr << "Cycle " << t << endl;
        vlasov.Execute();
        if ( (t%10)==0 ) {
            double En = vlasov.TotalEnergy();
            if (vlasov.master())
                TotalEngy << t << " " << En << endl;
            vlasov.writeYVySlice(10,8,8,"PhaseSlice.out");
        }
        if ( (t>=0) && vlasov.master() /*&& ((t%5) == 0) */ ) {
//        if (false) {
             write_Scalar(vlasov.GetEx(),Exfile);
             write_Scalar(vlasov.GetEy(),Eyfile);
             write_Scalar(vlasov.GetEz(),Ezfile);
             write_Scalar(vlasov.GetBx(),Bxfile, Gl_B0x);
             write_Scalar(vlasov.GetBy(),Byfile, Gl_B0y);
             write_Scalar(vlasov.GetBz(),Bzfile, Gl_B0z);
             write_Scalar(vlasov.getJy(),Jyfile);
             write_Scalar(vlasov.getJz(),Jzfile);
             write_Scalar(vlasov.Rho(),Rhofile);
        }
    }
    
    char fname[20];
    ostrstream ost(fname,20);
    ost << "PhaseSpace" << vlasov.procnum() << ".out";
    ofstream PhaseOut(fname);
    vlasov.write(PhaseOut);
    PhaseOut.close();
    
    return 0;
}




