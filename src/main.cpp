// -*- C++ -*-
// $Id$


#include "Vlas2d.h"
#include "wrapvlasov.h"
#include "boundary.h"
#include <fstream>
#include <sstream>
#include <string>
#include <list>

//typedef EFieldForce ForceField;
typedef EMDarwinForce ForceField;

typedef VlasovSpecies<ForceField,LeapFrogAdvance,PosFluxCons3rdOrder> *pVlasov;

int argc;
char **argv;

list<pVlasov> species;
Boundary *bound;
ForceField::FieldType *field;

void readGeneral(istream &In) {
  string token;
  while (true) {
    In >> token;
    if (token=="end") break;
    else if (token=="x-size") In >> GlGridX;
    else if (token=="y-size") In >> GlGridY;
    else if (token=="v-size") In >> GlVelRes;
    else if (token=="x-space") In >> GlGridSpace_x;
    else if (token=="y-space") In >> GlGridSpace_y;
//     else if (token=="vx-range") In >> GlGridRange_vx;
//     else if (token=="vy-range") In >> GlGridRange_vy;
//     else if (token=="vz-range") In >> GlGridRange_vz;
    else if (token=="dt") In >> Gl_dt;
    else if (token=="vel-ratio") In >> Gl_VelocityRatio;
    else
      cerr << "Unrecognised token " << token << endl;
  }
  GlGridHigh = PositionI(GlGridX+1,GlGridY+1);
  GlDistHigh = PhasePositionI(GlGridX+2,GlGridY+2,GlVelRes+1,GlVelRes+1,GlVelRes+1);
}

void readBoundary(istream &In) {
  string token;
  
  In >> token;
  if (token=="single-periodic") 
    bound = new SinglePeriodicBoundary;
#ifndef SINGLE_PROCESSOR
  else if (token=="mpi-x-periodic") 
    bound = new MPIPeriodicSplitXBoundary(argc,argv);
  else if (token=="mpi-xy-periodic") 
    bound = new MPIPeriodicSplitXYBoundary(argc,argv);
#endif // single processor
  else
    cerr << "Unrecognised token " << token << endl;

  field = new (ForceField::FieldType)(bound->master());
  field->Init();
}

VlasovInitialiser *readMaxwellInit(istream &In, pVlasov spec) {
  string token;
  VlasovMaxwellInitData data;
  
  while (true) {
    In >> token;
    if (token=="end") break;
    else if (token=="InitStream_vx") In >> data.InitStream_vx;
    else if (token=="InitStream_vy") In >> data.InitStream_vy;
    else if (token=="InitStream_vz") In >> data.InitStream_vz;
    else if (token=="InitTherm_vx") In >> data.InitTherm_vx;
    else if (token=="InitTherm_vy") In >> data.InitTherm_vy;
    else if (token=="InitTherm_vz") In >> data.InitTherm_vz;
    else if (token=="Init_N0") In >> data.Init_N0;
    else if (token=="Init_N1") In >> data.Init_N1;
    else if (token=="Init_kx") In >> data.Init_kx;
    else if (token=="Init_ky") In >> data.Init_ky;
    else
      cerr << "Unrecognised token " << token << endl;
  }
  data.vlasov = spec;
  return new VlasovMaxwellInit(data);
}

VlasovInitialiser *readTwoMaxwellInit(istream &In, pVlasov spec) {
  string token;
  VlasovTwoMaxwellInitData data;

    while (true) {
    In >> token;
    if (token=="end") break;
    else if (token=="InitStream_vx") In >> data.InitStream_vx;
    else if (token=="InitStream_vy") In >> data.InitStream_vy;
    else if (token=="InitStream_vz") In >> data.InitStream_vz;
    else if (token=="InitTherm_vx") In >> data.InitTherm_vx;
    else if (token=="InitTherm_vy") In >> data.InitTherm_vy;
    else if (token=="InitTherm_vz") In >> data.InitTherm_vz;
    else if (token=="Init_N0") In >> data.Init_N0;
    else if (token=="Init_N1") In >> data.Init_N1;
    else if (token=="Init_kx") In >> data.Init_kx;
    else if (token=="Init_ky") In >> data.Init_ky;
    else if (token=="InitStream2_vx") In >> data.InitStream2_vx;
    else if (token=="InitStream2_vy") In >> data.InitStream2_vy;
    else if (token=="InitStream2_vz") In >> data.InitStream2_vz;
    else if (token=="InitTherm2_vx") In >> data.InitTherm2_vx;
    else if (token=="InitTherm2_vy") In >> data.InitTherm2_vy;
    else if (token=="InitTherm2_vz") In >> data.InitTherm2_vz;
    else if (token=="Init2_N0") In >> data.Init2_N0;
    else if (token=="Init2_N1") In >> data.Init2_N1;
    else if (token=="Init2_kx") In >> data.Init2_kx;
    else if (token=="Init2_ky") In >> data.Init2_ky;
    else
      cerr << "Unrecognised token " << token << endl;
  }
  data.vlasov = spec;
  return new VlasovTwoMaxwellInit(data);
}

VlasovInitialiser *readWaveGenInit(istream &In, pVlasov spec) {
  string token;
  VlasovMaxwellInitData data;
  
  while (true) {
    In >> token;
    if (token=="end") break;
    else if (token=="InitStream_vx") In >> data.InitStream_vx;
    else if (token=="InitStream_vy") In >> data.InitStream_vy;
    else if (token=="InitStream_vz") In >> data.InitStream_vz;
    else if (token=="InitTherm_vx") In >> data.InitTherm_vx;
    else if (token=="InitTherm_vy") In >> data.InitTherm_vy;
    else if (token=="InitTherm_vz") In >> data.InitTherm_vz;
    else if (token=="Init_N0") In >> data.Init_N0;
    else if (token=="Init_N1") In >> data.Init_N1;
    else
      cerr << "Unrecognised token " << token << endl;
  }
  data.vlasov = spec;
  return new VlasovWaveGenInit(data);
}


void readInitialise(istream &In, SpeciesData &data) {
  string token;
  pVlasov spec 
      = new VlasovSpecies<ForceField,LeapFrogAdvance,PosFluxCons3rdOrder>(data);
  spec->setForceField(field);
  spec->Init();
  
  In >> token;
  if (token=="maxwell") spec->initialise(readMaxwellInit(In,spec));
  else if (token=="two-maxwell") spec->initialise(readTwoMaxwellInit(In,spec));
  else if (token=="wavegen") spec->initialise(readWaveGenInit(In,spec));
  else
    cerr << "Unrecognised token " << token << endl;
 
  species.push_back(spec);
}

void readSpecies(istream &In) {
  string token;
  SpeciesData data;

  data.bound = bound;
  
  while (true) {
    In >> token;
    if (token=="end") break;
    else if (token=="charge") In >> data.charge;
    else if (token=="mass") In >> data.mass;
    else if (token=="vx-range") In >> data.GridRange_vx;
    else if (token=="vy-range") In >> data.GridRange_vy;
    else if (token=="vz-range") In >> data.GridRange_vz;
    else if (token=="init") readInitialise(In,data);
    else
      cerr << "Unrecognised token " << token << endl;
  }
  
  
}

void readParameters() {
  ifstream In("setup.dat");
  string token;
      
  while (In >> token) {
    if (token=="general") 
      readGeneral(In);
    else if (token=="boundary") 
      readBoundary(In);
    else if (token=="species") 
      readSpecies(In);
    else
      cerr << "Unrecognised token " << token << endl;
  }

}

/** @brief Print some initial diagnostics and start the simulation.
 */
int main (int argc_, char** argv_) {

//    VlasovSpecies<EFieldForce> vlasov(new SinglePeriodicBoundary);
//    VlasovSpecies<EFieldForce> vlasov(new MPIPeriodicSplitXBoundary(argc,argv));
  argc = argc_;
  argv = argv_;
  
  readParameters();

  ofstream Exfile;   
  ofstream Eyfile;   
//   ofstream Ezfile;   
//   ofstream Bxfile;   
//   ofstream Byfile;   
//   ofstream Bzfile;   
//   ofstream Jxfile;   
//   ofstream Jyfile;   
//   ofstream Jzfile;   
//   ofstream Rhofile;
//   ofstream TotalEngy;
//   ofstream Phase;
    
  if (bound->master()) {
      Exfile.open("ExWaves.out");   
      Eyfile.open("EyWaves.out");   
//       Ezfile.open("EzWaves.out");   
//       Bxfile.open("BxWaves.out");   
//       Byfile.open("ByWaves.out");   
//       Bzfile.open("BzWaves.out");   
//       Jxfile.open("jxWaves.out");   
//       Jyfile.open("jyWaves.out");   
//       Jzfile.open("jzWaves.out");   
//       Rhofile.open("RhoWaves.out");
//       TotalEngy.open("TotalEngy.out");
//       Phase.open("PhaseVxVySlice.out");
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
  list<pVlasov>::iterator it;
  for (int t=0; t<300000; ++t) {
      cerr << "Cycle " << t << endl;
      for (it=species.begin(); it!=species.end(); ++it) {
        (*it)->Execute();
      }
      field->Execute(Gl_dt);
      
//       if ( (t%10)==0 ) {
//           double En = vlasov.TotalEnergy();
//           if (vlasov.master())
//               TotalEngy << t << " " << En << endl;
//           //vlasov.writeVxVySlice(3,3,0,"PhaseSlice.out");
//           vlasov.writeVxVySlice(t,2,3,16,Phase);        
//       }
      if ( (t>0) && bound->master() /* && ((t%5) == 0) */ ) {
//        if (false) {
           write_Scalar(field->GetEx(),Exfile);
           write_Scalar(field->GetEy(),Eyfile);
//            write_Scalar(field.GetEz(),Ezfile);
//            write_Scalar(field.GetBx(),Bxfile, Gl_B0x);
//            write_Scalar(field.GetBy(),Byfile, Gl_B0y);
//            write_Scalar(field.GetBz(),Bzfile, Gl_B0z);
//            write_Scalar(field.getJx(),Jxfile);
//            write_Scalar(field.getJy(),Jyfile);
//            write_Scalar(field.getJz(),Jzfile);
//            write_Scalar(field.Rho(),Rhofile);
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




