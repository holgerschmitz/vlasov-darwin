#include "poisson.h"

extern "C" {
  void hwscrt_(
    double *a, double *b, int *m, int *mbdcnd, 
    double *bda, double *bdb,
	double *c, double *d, int *n, int *nbdcnd, 
    double *bdc, double *bdd,
    double *elmbda, double *f, int *idimf,
    double *pertrb, int *ierror, double *w);
}

Poisson::Poisson() {
    WorkArray = NULL;
}

Poisson::Poisson(const PositionD &pMiN_, const PositionD &pMaX_, 
                 const PositionI &steps_,
                 bcond boundx, bcond boundy) {
    resize(pMiN_,pMaX_,steps_,boundx,boundy);
}

void Poisson::do_resize(const PositionD &pMiN_, const PositionD &pMaX_, 
                     const PositionI &steps_) {
    for (int i=0; i<2; ++i) {
        pMiN[i] = pMiN_[i];
        pMaX[i] = pMaX_[i];
        steps[i] = steps_[i]-1;
    }
    
    if (WorkArray) delete[] WorkArray;
    
    WorkArray = new double[ 4*(steps[0]+1) + 
                            (13 + int(log2(double(steps[0]+1))))
                                *(steps[1]+1)
                          ];
}

void Poisson::resize(const PositionD &pMiN_, const PositionD &pMaX_, 
                     const PositionI &steps_,
                     bcond boundx, bcond boundy) {
    do_resize(pMiN_, pMaX_, steps_);
    boundary[0] = boundx;
    boundary[1] = boundy;

}

void Poisson::resize(const PositionD &pMiN_, const PositionD &pMaX_, 
                     const PositionI &steps_) {
    do_resize(pMiN_, pMaX_, steps_);
}

void Poisson::BCond(bcond boundx, bcond boundy) {
    boundary[0] = boundx;
    boundary[1] = boundy;
}

Poisson::~Poisson() {
    if (WorkArray) delete[] WorkArray;
}

void Poisson::solve(ScalarField &In, ScalarField &Out) {
    PositionI i;
    int lowx, highx, lowy, highy;

    int imin[2], imax[2], nstep[2];

    assert( (boundary[0]==periodic) && (boundary[1]==periodic) );

    if(boundary[0]==periodic) {
        lowx = 0;
        highx = Out.getHigh(0);
        imin[0] = 1;
        imax[0] = highx;
    }
    
    if(boundary[1]==periodic) {
        lowy = 0;
        highy = Out.getHigh(1);
        imin[1] = 1;
        imax[1] = highy;
    }
    
    ScalarField Temp(imin,imax);
    
    for (i[0]=imin[0]; i[0]<=imax[0]; ++i[0])
        for (i[1]=imin[1]; i[1]<=imax[1]; ++i[1]) {
            Temp(i[0],i[1]) = In(i[0],i[1]);
//            std::cerr << "temp " << i << " " << Temp(i[0],i[1]) << std::endl;
        }
           
    for (int j=0; j<2; ++j) nstep[j] = imax[j]-imin[j];
    
    static double lambda=0;
    double perturb;
    int ierror;
    int idimf = nstep[0]+1;
    
    hwscrt_(&pMiN[0], &pMaX[0], &nstep[0], &boundary[0], NULL, NULL,
            &pMiN[1], &pMaX[1], &nstep[1], &boundary[1], NULL, NULL,
            &lambda, Temp.data(), &idimf,
            &perturb, &ierror, WorkArray);
            
    if (ierror != 0) {
        cerr << "Error "<<ierror<<" in Poisson" << endl;
        cerr << "Params are:"<< endl;
        cerr << pMiN[0] <<" "<< pMaX[0] <<" "<< nstep[0] <<" "<< boundary[0]<<endl;
        cerr << pMiN[1] <<" "<< pMaX[1] <<" "<< nstep[1] <<" "<< boundary[1]<<endl;
    }
    assert(ierror == 0);
    
    if(boundary[1]==periodic) {
        for (i[0]=imin[0]; i[0]<=imax[0]; ++i[0])
            for (i[1]=imin[1]; i[1]<=imax[1]; ++i[1]) {
                Out(i[0],i[1]) = Temp(i[0],i[1]);
//                std::cerr << "temp " << i << " " << Out(i[0],i[1]) << std::endl;
            }
        
        for (i[0]=0; i[0]<=imax[0]; ++i[0])
            Out(i[0],0) = Out(i[0],highy-1);
            
        for (i[1]=0; i[1]<=imax[1]; ++i[1])
            Out(0,i[1]) = Out(highx-1,i[1]);
    }
    
}   
