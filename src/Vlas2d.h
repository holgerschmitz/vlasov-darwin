#ifndef VLAS2D_H
#define VLAS2D_H

#include "vlasov.h"

//  Are being replaced by the following:
//  No! They aren't

extern const int GlGridX;
extern const int GlGridY;
extern const int GlVelRes;

extern PositionI GlGridLow; // The full field size for scalar field
extern PositionI GlGridHigh; 

extern PhasePositionI GlDistLow; // The processes part of the distribution
extern PhasePositionI GlDistHigh;

const extern double GlGridRange_vx;
const extern double GlGridRange_vy;
const extern double GlGridRange_vz;

const extern double InitStream_vx;
const extern double InitStream_vy;
const extern double InitStream_vz;

const extern double InitStream2_vx;
const extern double InitStream2_vy;
const extern double InitStream2_vz;

const extern double InitTherm_vx;
const extern double InitTherm_vy;
const extern double InitTherm_vz;

const extern double InitTherm2_vx;
const extern double InitTherm2_vy;
const extern double InitTherm2_vz;

const extern double Init_N0;
const extern double Init_N1;

const extern double Init2_N0;
const extern double Init2_N1;

const extern int Init_kx;
const extern int Init_ky;
const extern int Init_kz;

const extern int Init2_kx;
const extern int Init2_ky;
const extern int Init2_kz;

const extern int MG_nu1;
const extern int MG_nu2;
const extern int MG_gama;
const extern double MG_epsilon;

const extern double GlGridSpace_x;
const extern double GlGridSpace_y;

const extern double GlVolumeQuant;

const extern int Gl_B0x;
const extern int Gl_B0y;
const extern int Gl_B0z;

const extern double Gl_dt;
const extern double Gl_Charge;
const extern double Gl_Mass;

// This velocity ratio is the ratio of the normalization velocity
// (normally the thermal velocity of the electrons) to the speed of
// light
const extern double Gl_VelocityRatio;

#endif
