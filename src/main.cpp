
#include "Vlas2d.h"
#include "wrapvlasov.h"
#include <fstream>
#define X 0
#define Y 1
#define Z 2

#define R     0 // Radius
#define PHI   1 // gemessen von der X-Achse aus (positiv) [0..360]
#define THETA 2 // gemessen von der Z-Achse aus (positiv) [0..180]



/** @brief Print some initial diagnostics and start the simulation.
 *
 *  If an argument is specified at the command line then it
 *  is interpreted as a file name and the setup is read from the
 *  corresponding file. If no argument is specified OPAR tries to read
 *  the setup from the standard input.
 */
int main (int argc, char** argv) {

    init_VlasovParams();
    
    VlasovSpecies<EMDarwinForce> vlasov;
        
    vlasov.Init();
    
    ofstream Exfile("ExWaves.out");
    ofstream Eyfile("EyWaves.out");
    ofstream Ezfile("EzWaves.out");
    ofstream Bxfile("BxWaves.out");
    ofstream Byfile("ByWaves.out");
    ofstream Bzfile("BzWaves.out");
    ofstream Rhofile("RhoWaves.out");
    
    for (int t=0; t<11264; ++t) {
        cerr << "Cycle " << t << endl;
        vlasov.Execute();
        if ( (t>=3072) && (3*(t/3)==t) ) {
            write_Scalar(vlasov.GetEx(),Exfile);
            write_Scalar(vlasov.GetEy(),Eyfile);
            write_Scalar(vlasov.GetEz(),Ezfile);
            write_Scalar(vlasov.GetBx(),Bxfile);
            write_Scalar(vlasov.GetBy(),Byfile);
            write_Scalar(vlasov.GetBz(),Bzfile);
            write_Scalar(vlasov.Rho(),Rhofile);
        }
    }
     
    ofstream PhaseOut("PhaseSpace.out");
    vlasov.write(PhaseOut);
    PhaseOut.close();
    
    return 0;
}




