// -*- C++ -*-
// $Id$

#include "Vlas2d.h"

int GlGridX = 256;
int GlGridY = 4;
//const int GlGridX = 64;
//const int GlGridY = 64;
int GlVelRes = 20;

PositionI GlGridLow(0,0);
PositionI GlGridHigh(GlGridX+1,GlGridY+1);

PhasePositionI GlDistLow(-1,-1,0,0,0);
PhasePositionI GlDistHigh(GlGridX+2,GlGridY+2,GlVelRes+1,GlVelRes+1,GlVelRes+1);

//double GlGridRange_vx = 5.5;
//double GlGridRange_vy = 5.5;
//double GlGridRange_vz = 5.5;

const double InitStream_vx = 0.0;
const double InitStream_vy = 0.2;
const double InitStream_vz = 0.2;

const double InitStream2_vx = 1;
const double InitStream2_vy = 1;
const double InitStream2_vz = 1;

// const double InitTherm_vx = 1.4142136;
// const double InitTherm_vy = 1.4142136;
// const double InitTherm_vz = 1.4142136;
const double InitTherm_vx = 0.3;
const double InitTherm_vy = 0.3;
const double InitTherm_vz = 0.3;

const double InitTherm2_vx = 0.5;
const double InitTherm2_vy = 0.5;
const double InitTherm2_vz = 0.5;

const double Init_N0 = 1.0;
const double Init_N1 = 0.00;

const double Init2_N0 = 0.5;
const double Init2_N1 = 0.1;

const int MG_nu1 = 1;
const int MG_nu2 = 1;
const int MG_gama = 1;
const double MG_epsilon = 1e-6;

double GlGridSpace_x = 0.39269908; //4*M_PI/double(GlGridX);
double GlGridSpace_y = 0.39269908; //4*M_PI/double(GlGridY);

//const double GlGridSpace_x = 2*M_PI/double(GlGridX);
//const double GlGridSpace_y = 2*M_PI/double(GlGridY);

const double GlVolumeQuant = 1;

const int Init_kx = 1;
const int Init_ky = 0;
const int Init_kz = 0;

const int Init2_kx = 1;
const int Init2_ky = 2;
const int Init2_kz = 0;

const int Gl_B0x = 3;
const int Gl_B0y = 0;
const int Gl_B0z = 0;

double Gl_dt = 5e-3;

double Gl_VelocityRatio = 0.1;
