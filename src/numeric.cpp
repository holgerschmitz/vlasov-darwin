/***************************************************************************
                          numeric.cpp  -  description
                             -------------------
    begin                : Fri Jun 2 2000
    copyright            : (C) 2000 by Holger Schmitz
    email                : 
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include "numeric.h"
#include <iostream>


#ifdef TWO_DIMENSIONAL

#endif

#ifdef THREE_DIMENSIONAL

#endif

//using namespace std;

double gammln(double xx)
{
	double x,tmp,ser;
	static double cof[6]={76.18009173,-86.50532033,24.01409822,
		-1.231739516,0.120858003e-2,-0.536382e-5};
	int j;

	x=xx-1.0;
	tmp=x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser=1.0;
	for (j=0;j<=5;j++) {
		x += 1.0;
		ser += cof[j]/x;
	}
	return -tmp+log(2.50662827465*ser);
}


#define ITMAX 100
#define EPS 3.0e-7

void gser(double &gamser, double a, double x, double &gln)
{
	int n;
	double sum,del,ap;

	gln=gammln(a);
	if (x <= 0.0) {
		if (x < 0.0) {
            cerr << "x less than 0 in routine GSER" << endl;
            exit(-1);
        }
		gamser=0.0;
		return;
	} else {
		ap=a;
		del=sum=1.0/a;
		for (n=1;n<=ITMAX;n++) {
			ap += 1.0;
			del *= x/ap;
			sum += del;
			if (fabs(del) < fabs(sum)*EPS) {
				gamser=sum*exp(-x+a*log(x)-gln);
				return;
			}
		}
        cerr << "a too large, ITMAX too small in routine GSER" << endl;
        exit(-1);
		return;
	}
}

void gcf(double &gammcf, double a, double x, double &gln)
{
	int n;
	double gold=0.0,g,fac=1.0,b1=1.0;
	double b0=0.0,anf,ana,an,a1,a0=1.0;

	gln=gammln(a);
	a1=x;
	for (n=1;n<=ITMAX;n++) {
		an=(double) n;
		ana=an-a;
		a0=(a1+a0*ana)*fac;
		b0=(b1+b0*ana)*fac;
		anf=an*fac;
		a1=x*a0+anf*a1;
		b1=x*b0+anf*b1;
		if (a1) {
			fac=1.0/a1;
			g=b1*fac;
			if (fabs((g-gold)/g) < EPS) {
				gammcf=exp(-x+a*log(x)-gln)*g;
				return;
			}
			gold=g;
		}
	}
    cerr << "a too large, ITMAX too small in routine GCF" << endl;
    exit(-1);
}

#undef ITMAX
#undef EPS

double gammp(double a, double x)
{
	double gamser,gammcf,gln;

	if (x < 0.0 || a <= 0.0) {
        cerr << "Invalid arguments in routine GAMMP" << endl;
        exit(-1);
    }
	if (x < (a+1.0)) {
		gser(gamser,a,x,gln);
		return gamser;
	} else {
		gcf(gammcf,a,x,gln);
		return 1.0-gammcf;
	}
}


double erf(double x)
{
	return x < 0.0 ? -gammp(0.5,x*x) : gammp(0.5,x*x);
}




