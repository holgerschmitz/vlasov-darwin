// -*- C++ -*-
// $Id$

#include "helmholtz.h"
#include "Vlas2d.h"

#define ISODD(i) ((i)&1)     

Helmholtz::Helmholtz() {

    epsilon = MG_epsilon;
    gama    = MG_gama;
    nu1     = MG_nu1;
    nu2     = MG_nu2;

    // dx
    dx[0] = GlGridSpace_x;
    dx[1] = GlGridSpace_y;


}

// solves -laplace u = f
void Helmholtz::solve( NumMatrix<double,2> &u,
		               NumMatrix<double,2> &f,
                       NumMatrix<double,2> &lambda) {

    NumMatrix<double,2> uold(u);
    double error;
    do {
        uold=u;
        mgi(u,f,lambda);
        error=distance(uold,u);
//        cout << "Iterating Error = " << error << " with epsilon = " << epsilon << endl;
    } while (error > epsilon);
}

void Helmholtz::mgi( NumMatrix<double,2> &u,
		             NumMatrix<double,2> &f,
                     NumMatrix<double,2> &lambda) {
    int	i;
    int mx=u.getHigh(0)-1;
    int my=u.getHigh(1)-1;

    if (ISODD(mx) || ISODD(my)) {
        double error;

        NumMatrix<double,2> uold(u);

        do {
            uold=u;
            gauss(u,f,lambda);
            error=distance(uold,u);
        } while(error > epsilon);
        
    } else {
        NumMatrix<double,2> un(Index::set(0,0),Index::set(mx/2+1,my/2+1));
        NumMatrix<double,2> fn(Index::set(0,0),Index::set(mx/2+1,my/2+1));
        NumMatrix<double,2> lambdan(Index::set(0,0),Index::set(mx/2+1,my/2+1));

        for(i = 1; i <= nu1; i++)
            gauss(u,f,lambda);

        fn.clear();
        defect(u,f,lambda,fn,lambdan);

        un.clear();
        for(i = 1; i <= gama; i++)
            mgi(un,fn,lambdan);

        prolongate(u,un);

        for(i=1; i <= nu2; i++)
            gauss(u,f,lambda);
    }
}



void Helmholtz::gauss( NumMatrix<double,2> &u,
		               NumMatrix<double,2> &f,
                       NumMatrix<double,2> &lambda) {
    int	i,j;
    int mx=u.getHigh(0)-1;
    int my=u.getHigh(1)-1;
    double invdx2=sqr(1./dx[0]);
    double invdy2=sqr(1./dx[1]);

    //   even block

    for(j = 1; j <= my; j += 2) {
        for(i = 1; i <= mx; i += 2) {
            u(i,j)=0.5/(invdx2+invdy2+0.5*lambda(i,j))
                *(invdy2*(u(i,j+1)+u(i,j-1))
	             +invdx2*(u(i+1,j)+u(i-1,j))+f(i,j));
        }
    }

    for(j = 2; j <= my; j += 2) {
        for(i = 2; i <= mx; i += 2) {
            u(i,j)=0.5/(invdx2+invdy2+0.5*lambda(i,j))
                *(invdy2*(u(i,j+1)+u(i,j-1))
	             +invdx2*(u(i+1,j)+u(i-1,j))+f(i,j));
        }
    }

    boundary(u);

    //   odd block

    for(j = 2; j <= my; j += 2) {
        for(i = 1; i <= mx; i += 2) {
            u(i,j)=0.5/(invdx2+invdy2+0.5*lambda(i,j))
                *(invdy2*(u(i,j+1)+u(i,j-1))
	             +invdx2*(u(i+1,j)+u(i-1,j))+f(i,j));
        }
    }

    for(j = 1; j <= my; j += 2) {
        for(i = 2; i <= mx; i += 2) {
            u(i,j)=0.5/(invdx2+invdy2+0.5*lambda(i,j))
                *(invdy2*(u(i,j+1)+u(i,j-1))
	             +invdx2*(u(i+1,j)+u(i-1,j))+f(i,j));
        }
    }

    boundary(u);
}

void Helmholtz::defect( NumMatrix<double,2> &u,
                        NumMatrix<double,2> &f,
                        NumMatrix<double,2> &lambda,
                        NumMatrix<double,2> &fn,
                        NumMatrix<double,2> &lambdan) {
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
            d(i,j)=2. * (invdx2+invdy2+0.5*lambda(i,j)) * u(i,j)
                        -invdy2*(u(i,j+1)+u(i,j-1))
	                    -invdx2*(u(i+1,j)+u(i-1,j))
                        -f(i,j);
        }
    }

    boundary(d);

    //   restriction
    int ih,jh;

    for(j = 1; j <= my; j += 2) {
        for(i = 1; i <= mx; i += 2) {
            ih = (i-1)/2+1;
            jh = (j-1)/2+1;
            fn(ih,jh)=0.25*( d(i  ,j)+d(i  ,j+1)
				            +d(i+1,j)+d(i+1,j+1));
                                     
            lambdan(ih,jh)=0.25*( lambda(i  ,j)+lambda(i  ,j+1)
				                 +lambda(i+1,j)+lambda(i+1,j+1));
        }
    }
}


double Helmholtz::distance( NumMatrix<double,2> &uold,
			                NumMatrix<double,2> &u) {
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


void Helmholtz::normalize(NumMatrix<double,2> &u) {
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

void Helmholtz::prolongate( NumMatrix<double,2> &u,
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

    boundary(u);
}

void Helmholtz::boundary(NumMatrix<double,2> &u) {
    int mx=u.getHigh(0)-1;
    int my=u.getHigh(1)-1;

    int i,j;

    for(i = 1; i <= mx; i++) {
        u(i,0   ) = u(i,my);
        u(i,my+1) = u(i,1 );
    }
    for(j = 0; j <= my+1; j++) {
        u(0,j   ) = u(mx,j);
        u(mx+1,j) = u(1,j);
    }

}
