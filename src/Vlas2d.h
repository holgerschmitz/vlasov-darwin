#ifndef VLAS2D_H
#define VLAS2D_H

extern int GlGridSize_x;
extern int GlGridSize_y;

extern int GlGridSize_vx;
extern int GlGridSize_vy;
extern int GlGridSize_vz;

extern double GlGridRange_vx;
extern double GlGridRange_vy;
extern double GlGridRange_vz;

extern double InitStream_vx;
extern double InitStream_vy;
extern double InitStream_vz;

extern double InitTherm_vx;
extern double InitTherm_vy;
extern double InitTherm_vz;

extern double Init_N0;
extern double Init_N1;

extern int Init_kx;
extern int Init_ky;
extern int Init_kz;

extern int MG_nu1;
extern int MG_nu2;
extern int MG_gama;
extern double MG_epsilon;

extern double GlGridSpace_x;
extern double GlGridSpace_y;

extern double GlVolumeQuant;

void init_VlasovParams();

#endif
