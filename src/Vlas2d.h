// -*- C++ -*-
// $Id$

#ifndef VLAS2D_H
#define VLAS2D_H

#include "vlasov.h"

extern int GlGridX;
extern int GlGridY;
extern int GlVelRes;

extern PositionI GlGridLow; // The full field size for scalar field
extern PositionI GlGridHigh; 

extern PhasePositionI GlDistLow; // The processes part of the distribution
extern PhasePositionI GlDistHigh;

//extern double GlGridRange_vx;
//extern double GlGridRange_vy;
//extern double GlGridRange_vz;

const extern int MG_nu1;
const extern int MG_nu2;
const extern int MG_gama;
const extern double MG_epsilon;

extern double GlGridSpace_x;
extern double GlGridSpace_y;

const extern double GlVolumeQuant;

const extern int Gl_B0x;
const extern int Gl_B0y;
const extern int Gl_B0z;

extern double Gl_dt;

// This velocity ratio is the ratio of the normalization velocity
// (normally the thermal velocity of the electrons) to the speed of
// light
extern double Gl_VelocityRatio;

#endif
