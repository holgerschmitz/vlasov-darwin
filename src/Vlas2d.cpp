#include "Vlas2d.h"

int GlGridSize_x;
int GlGridSize_y;

int GlGridSize_vx;
int GlGridSize_vy;
int GlGridSize_vz;

double GlGridRange_vx;
double GlGridRange_vy;
double GlGridRange_vz;

double InitStream_vx;
double InitStream_vy;
double InitStream_vz;

double InitTherm_vx;
double InitTherm_vy;
double InitTherm_vz;

double Init_N0;
double Init_N1;

int MG_nu1;
int MG_nu2;
int MG_gama;
double MG_epsilon;

double GlGridSpace_x;
double GlGridSpace_y;

double GlVolumeQuant;

int Init_kx;
int Init_ky;
int Init_kz;

void init_VlasovParams() {
    GlGridSize_x = 32;
    GlGridSize_y = 32;

    GlGridSize_vx = 22;
    GlGridSize_vy = 22;
    GlGridSize_vz = 22;

    GlGridRange_vx = 9;
    GlGridRange_vy = 9;
    GlGridRange_vz = 9;
    
    GlGridSpace_x = 0.392699;
    GlGridSpace_y = 0.392699;
    
    InitStream_vx = 0;
    InitStream_vy = 0;
    InitStream_vz = 0;

    InitTherm_vx = 1.4142136;
    InitTherm_vy = 1.4142136;
    InitTherm_vz = 1.4142136;

    Init_N0 = 1;
    Init_N1 = 0.1;

    Init_kx = 1;
    Init_ky = 0;
    Init_kz = 0;

    GlVolumeQuant = 1;
    
    MG_nu1 = 1;
    MG_nu2 = 1;
    MG_gama = 1;
    MG_epsilon = 1e-6;
}
