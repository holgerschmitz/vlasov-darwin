// -*- C++ -*-
// $Id$

#include "poisson.h"

#include "helmholtz.h"
#include "globals.h"

#define ISODD(i) ((i)&1)     

Poisson::Poisson() {

    epsilon = Parameters::instance().mgEpsilon();
    gama    = Parameters::instance().mgGamma();
    nu1     = Parameters::instance().mgNu1();
    nu2     = Parameters::instance().mgNu2();

    // dx
    dx[0] = Parameters::instance().gridSpace_x();
    dx[1] = Parameters::instance().gridSpace_y();
    
    boundary = NULL;
}

// solves -laplace u = f
void Poisson::solve( NumMatrix<double,2> &u,
                     NumMatrix<double,2> &f,
                     const NumBoundary &boundary_) {

    boundary = &boundary_;

    normalize(f);
    normalize(u);
    
    NumMatrix<double,2> uold(u);
    double error;
    do {
//    	cerr << "Poisson step\n";
        uold=u;
        mgi(u,f);
        normalize(u);
        error=distance(uold,u);
//        cerr << "Poisson Iterating Error = " << error << " with epsilon = " << epsilon << endl;
    } while (error > epsilon);
    
    boundary = NULL;
}

void Poisson::mgi( NumMatrix<double,2> &u,
                   NumMatrix<double,2> &f,
		   bool norm) 
{
    int  i;
    int mx=u.getHigh(0)-1;
    int my=u.getHigh(1)-1;

    if (ISODD(mx) || ISODD(my)) {
        double error;

        NumMatrix<double,2> uold(u);

        do {
            uold=u;
            gauss(u,f);
            error=distance(uold,u);
//            cout << "Poisson Iterating Error = " << error << " with epsilon = " << epsilon 
//                 <<"("<<mx<<","<<my<<")"<< endl;
        } while(error > epsilon);
        
    } else {
        NumMatrix<double,2> un(Index::set(0,0),Index::set(mx/2+1,my/2+1));
        NumMatrix<double,2> fn(Index::set(0,0),Index::set(mx/2+1,my/2+1));

        for(i = 1; i <= nu1; i++) {
            gauss(u,f);
	    if (norm || boundary->normalize()) normalize(u);
	}

        fn.clear();
        defect(u,f,fn);

        un.clear();
        for(i = 1; i <= gama; i++)
            mgi(un,fn);

        prolongate(u,un);

        for(i=1; i <= nu2; i++) {
            gauss(u,f);
	    if (norm || boundary->normalize()) normalize(u);
	}

    }
}



void Poisson::gauss( NumMatrix<double,2> &u,
                     NumMatrix<double,2> &f) {
    int  i,j;
    int mx=u.getHigh(0)-1;
    int my=u.getHigh(1)-1;
    double invdx2=sqr(1./dx[0]);
    double invdy2=sqr(1./dx[1]);

    //   even block

    for(j = 1; j <= my; j += 2) {
        for(i = 1; i <= mx; i += 2) {
            u(i,j)=0.5/(invdx2+invdy2)
                *(invdy2*(u(i,j+1)+u(i,j-1))
              +invdx2*(u(i+1,j)+u(i-1,j))+f(i,j));
        }
    }

    for(j = 2; j <= my; j += 2) {
        for(i = 2; i <= mx; i += 2) {
            u(i,j)=0.5/(invdx2+invdy2)
                *(invdy2*(u(i,j+1)+u(i,j-1))
              +invdx2*(u(i+1,j)+u(i-1,j))+f(i,j));
        }
    }

    boundary->apply(u);

    //   odd block

    for(j = 2; j <= my; j += 2) {
        for(i = 1; i <= mx; i += 2) {
            u(i,j)=0.5/(invdx2+invdy2)
                *(invdy2*(u(i,j+1)+u(i,j-1))
               +invdx2*(u(i+1,j)+u(i-1,j))+f(i,j));
        }
    }

    for(j = 1; j <= my; j += 2) {
        for(i = 2; i <= mx; i += 2) {
            u(i,j)=0.5/(invdx2+invdy2)
                *(invdy2*(u(i,j+1)+u(i,j-1))
              +invdx2*(u(i+1,j)+u(i-1,j))+f(i,j));
        }
    }

    boundary->apply(u);
}

void Poisson::defect( NumMatrix<double,2> &u,
                        NumMatrix<double,2> &f,
                        NumMatrix<double,2> &fn) {
    int i,j;
    int mx=u.getHigh(0)-1;
    int my=u.getHigh(1)-1;
    double invdx2=sqr(1./dx[0]);
    double invdy2=sqr(1./dx[1]);

    NumMatrix<double,2> d(f);
    d.clear();
    
    //   defect

    for(j = 1; j <= my; j++) {
        for(i = 1; i <= mx; i++) {
            d(i,j)=2. * (invdx2+invdy2) * u(i,j)
                        -invdy2*(u(i,j+1)+u(i,j-1))
                        -invdx2*(u(i+1,j)+u(i-1,j))
                        -f(i,j);
        }
    }

    boundary->apply(d);

    //   restriction
    int ih,jh;

    for(j = 1; j <= my; j += 2) {
        for(i = 1; i <= mx; i += 2) {
            ih = (i-1)/2+1;
            jh = (j-1)/2+1;
            fn(ih,jh)=0.25*( d(i  ,j)+d(i  ,j+1)
                   +d(i+1,j)+d(i+1,j+1));
        }
    }

}


double Poisson::distance( NumMatrix<double,2> &uold,
                          NumMatrix<double,2> &u) 
{
    int mx=u.getHigh(0)-1;
    int my=u.getHigh(1)-1;

    double diff=0.;
    
    for(int j = 1; j <= my; j++) {
        for(int i = 1; i <= mx; i++) {
            diff = max(diff,abs(uold(i,j)-u(i,j)));
        }
    }
    
    return diff;
}


void Poisson::normalize(NumMatrix<double,2> &u) {
    
    if (!boundary->normalize()) return;

    int mx=u.getHigh(0)-1;
    int my=u.getHigh(1)-1;

    int i,j;

    double h=1./(mx*my);

    double sum = 0.;

    for(j = 1; j <= my; j++) {
        for(i = 1; i <= mx; i++) {
            sum += u(i,j);
        }
    }

    sum *= h;

    for (j = 0; j <= my+1; j++) { 
        for (i = 0; i <= mx+1; i++) {
            u(i,j) -= sum;
        }
    }
}

void Poisson::prolongate( NumMatrix<double,2> &u,
                            NumMatrix<double,2> &un) {
    int mx=u.getHigh(0)-1;
    int my=u.getHigh(1)-1;

    int i,j;


    NumMatrix<double,2> w(u);

    for (j = 1; j <= my; j++) { 
      for (i = 1; i <= mx; i++) {
        w(i,j) = un((i-1)/2+1,(j-1)/2+1);
      }
    }

    u -= w;

    boundary->apply(u);
}

