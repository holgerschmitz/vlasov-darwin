// -*- C++ -*-
// $Id$

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
    PositionI steps;
    for (int i=0; i<2; ++i) {
        pMiN[i] = pMiN_[i];
        pMaX[i] = pMaX_[i];
        steps[i] = 2*steps_[i];
    }
    
    if (WorkArray) delete[] WorkArray;
    
    WorkArray = new double[ (4*(steps[1]+2) + 
                            (13 + int(log2(double(steps[1]+2))))
                                *(steps[0]+2))
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
    // only periodic boundariy conditions are implemented
    assert( (boundary[0]==periodic) && (boundary[1]==periodic) );

    const int *low = In.getLow();
    const int *high = In.getHigh();
    
    int im, iM, ip, iP;
    int jm, jM, jp, jP;
    int lx = high[0]-low[0]-2;
    int ly = high[1]-low[1]-2;
    double dx = (pMaX[0]-pMiN[0])/double(lx+1);
    double dy = (pMaX[1]-pMiN[1])/double(ly+1);
    
    ScalarField ResIn(low,high);
    ScalarField ResOut(low,high);
    
//    cerr << "Poisson solving step 1\n";
    
//     for (int i=low[0]; i<=high[0]; ++i)
//         for (int j=low[1]; j<=high[1]; ++j) {
//             cerr << i << " " << j << " " << In(i,j) << " " << Out(i,j) << "\n";
//             TmpIn(i,j) = In(i,j);
// 
//         }

    solve_step(In,Out);
    
//     for (int i=low[0]; i<=high[0]; ++i)
//         for (int j=low[1]; j<=high[1]; ++j) {
//             cerr << i << " " << j << " " << TmpIn(i,j) << " " << In(i,j) << " " << Out(i,j) << "\n";
// 
//         }

//    cerr << "Poisson solving step 2\n";
    
    for (int i=low[0]; i<=high[0]; ++i)
        for (int j=low[1]; j<=high[1]; ++j) {
            im = i  > low[0]  ? i  - 1: i +lx;
            iM = im > low[0]  ? im - 1: im+lx;
            ip = i  < high[0] ? i  + 1: i -lx;
            iP = ip < high[0] ? ip + 1: ip-lx;
            jm = j  > low[1]  ? j  - 1: j +ly;
            jM = jm > low[1]  ? jm - 1: jm+ly;
            jp = j  < high[1] ? j  + 1: j -ly;
            jP = jp < high[1] ? jp + 1: jp-ly;
            ResIn(i,j) = In(i,j) - ( -60*Out(i,j)
                    + 16*(
                        Out(ip,j) + Out(i,jp) + Out(im,j) + Out(i,jm)
                    )
                    - (
                        Out(iP,j) + Out(i,jP) + Out(iM,j) + Out(i,jM)
                    )
                )/(12*dx*dy);
        }
    
//    cerr << "Poisson solving step 3\n";
    solve_step(ResIn,ResOut); 
    
//    cerr << "Poisson solving step 4\n";
    
    for (int i=low[0]; i<=high[0]; ++i) {
        for (int j=low[1]; j<=high[1]; ++j) {
            Out(i,j) += ResOut(i,j);
            // cout << i << " " << j << " " << Out(i,j) << "\n" << flush;
        }
        
    }
    
}

void Poisson::solve_step(ScalarField &In, ScalarField &Out) {
    PositionI i;
    int lowx, highx, lowy, highy;

    int imin[2], imax[2], nstep[2];


    // for periodic boundaries we only take part of the array
    if(boundary[0]==periodic) {
        lowx = In.getLow(0);
        highx = In.getHigh(0);
        imin[0] = lowx+1;
        imax[0] = highx;
    }
    
    if(boundary[1]==periodic) {
        lowy = In.getLow(1);
        highy = In.getHigh(1);
        imin[1] = lowy+1;
        imax[1] = highy;
    }
    
//    cerr << "Poisson\n";
//      cerr << "    " << lowx << " " << highx << "\n";
//      cerr << "    " << lowy << " " << highy << "\n";
//      cerr << "    " << imin[0] << " " << imax[0] << "\n";
//      cerr << "    " << imin[1] << " " << imax[1] << "\n";
    
    ScalarField Temp(imin,imax);
    
    for (i[0]=imin[0]; i[0]<=imax[0]; ++i[0])
        for (i[1]=imin[1]; i[1]<=imax[1]; ++i[1])
            Temp(i[0],i[1]) = In(i[0],i[1]);
           
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
            for (i[1]=imin[1]; i[1]<=imax[1]; ++i[1])
                Out(i[0],i[1]) = Temp(i[0],i[1]);
        
        for (i[0]=lowx; i[0]<=highx; ++i[0])
            Out(i[0],lowy) = Out(i[0],highy-1);
            
        for (i[1]=lowy; i[1]<=highy; ++i[1])
            Out(lowx,i[1]) = Out(highx-1,i[1]);
    }
    
} 

/** @page hwscrt Fortran hwscrt routine
 *
 *
 *  @section purpose_hwscrt PURPOSE
 *
 *          SUBROUTINE HWSCRT SOLVES THE STANDARD FIVE-POINT FINITE
 *     DIFFERENCE APPROXIMATION TO THE HELMHOLTZ EQUATION IN CARTESIAN
 *     COORDINATES:
 *
 *          (D/DX)(DU/DX) + (D/DY)(DU/DY) + LAMBDA*U = F(X,Y).
 *
 *
 *
 *     * * * * * * * *    PARAMETER DESCRIPTION     * * * * * * * * * *
 *
 *             * * * * * *   ON INPUT    * * * * * *
 *
 *     A,B
 *       THE RANGE OF X, I.E., A .LE. X .LE. B.  A MUST BE LESS THAN B.
 *
 *     M
 *       THE NUMBER OF PANELS INTO WHICH THE INTERVAL (A,B) IS
 *       SUBDIVIDED.  HENCE, THERE WILL BE M+1 GRID POINTS IN THE
 *       X-DIRECTION GIVEN BY X(I) = A+(I-1)DX FOR I = 1,2,...,M+1,
 *       WHERE DX = (B-A)/M IS THE PANEL WIDTH. M MUST BE GREATER THAN 3.
 *
 *     MBDCND
 *       INDICATES THE TYPE OF BOUNDARY CONDITIONS AT X = A AND X = B.
 *
 *       = 0  IF THE SOLUTION IS PERIODIC IN X, I.E., U(I,J) = U(M+I,J).
 *       = 1  IF THE SOLUTION IS SPECIFIED AT X = A AND X = B.
 *       = 2  IF THE SOLUTION IS SPECIFIED AT X = A AND THE DERIVATIVE OF
 *            THE SOLUTION WITH RESPECT TO X IS SPECIFIED AT X = B.
 *      = 3  IF THE DERIVATIVE OF THE SOLUTION WITH RESPECT TO X IS
 *           SPECIFIED AT X = A AND X = B.
 *      = 4  IF THE DERIVATIVE OF THE SOLUTION WITH RESPECT TO X IS
 *           SPECIFIED AT X = A AND THE SOLUTION IS SPECIFIED AT X = B.
 *
 *    BDA
 *      A ONE-DIMENSIONAL ARRAY OF LENGTH N+1 THAT SPECIFIES THE VALUES
 *      OF THE DERIVATIVE OF THE SOLUTION WITH RESPECT TO X AT X = A.
 *      WHEN MBDCND = 3 OR 4,
 *
 *           BDA(J) = (D/DX)U(A,Y(J)), J = 1,2,...,N+1  .
 *
 *      WHEN MBDCND HAS ANY OTHER VALUE, BDA IS A DUMMY VARIABLE.
 *
 *    BDB
 *      A ONE-DIMENSIONAL ARRAY OF LENGTH N+1 THAT SPECIFIES THE VALUES
 *      OF THE DERIVATIVE OF THE SOLUTION WITH RESPECT TO X AT X = B.
 *      WHEN MBDCND = 2 OR 3,
 *
 *           BDB(J) = (D/DX)U(B,Y(J)), J = 1,2,...,N+1  .
 *
 *      WHEN MBDCND HAS ANY OTHER VALUE BDB IS A DUMMY VARIABLE.
 *
 *    C,D
 *      THE RANGE OF Y, I.E., C .LE. Y .LE. D.  C MUST BE LESS THAN D.
 *
 *    N
 *      THE NUMBER OF PANELS INTO WHICH THE INTERVAL (C,D) IS
 *      SUBDIVIDED.  HENCE, THERE WILL BE N+1 GRID POINTS IN THE
 *      Y-DIRECTION GIVEN BY Y(J) = C+(J-1)DY FOR J = 1,2,...,N+1, WHERE
 *      DY = (D-C)/N IS THE PANEL WIDTH.  N MUST BE GREATER THAN 3.
 *
 *    NBDCND
 *      INDICATES THE TYPE OF BOUNDARY CONDITIONS AT Y = C AND Y = D.
 *
 *      = 0  IF THE SOLUTION IS PERIODIC IN Y, I.E., U(I,J) = U(I,N+J).
 *      = 1  IF THE SOLUTION IS SPECIFIED AT Y = C AND Y = D.
 *      = 2  IF THE SOLUTION IS SPECIFIED AT Y = C AND THE DERIVATIVE OF
 *           THE SOLUTION WITH RESPECT TO Y IS SPECIFIED AT Y = D.
 *      = 3  IF THE DERIVATIVE OF THE SOLUTION WITH RESPECT TO Y IS
 *           SPECIFIED AT Y = C AND Y = D.
 *      = 4  IF THE DERIVATIVE OF THE SOLUTION WITH RESPECT TO Y IS
 *           SPECIFIED AT Y = C AND THE SOLUTION IS SPECIFIED AT Y = D.
 *
 *    BDC
 *      A ONE-DIMENSIONAL ARRAY OF LENGTH M+1 THAT SPECIFIES THE VALUES
 *      OF THE DERIVATIVE OF THE SOLUTION WITH RESPECT TO Y AT Y = C.
 *      WHEN NBDCND = 3 OR 4,
 *
 *           BDC(I) = (D/DY)U(X(I),C), I = 1,2,...,M+1  .
 *
 *      WHEN NBDCND HAS ANY OTHER VALUE, BDC IS A DUMMY VARIABLE.
 *
 *    BDD
 *      A ONE-DIMENSIONAL ARRAY OF LENGTH M+1 THAT SPECIFIES THE VALUES
 *      OF THE DERIVATIVE OF THE SOLUTION WITH RESPECT TO Y AT Y = D.
 *      WHEN NBDCND = 2 OR 3,
 *
 *           BDD(I) = (D/DY)U(X(I),D), I = 1,2,...,M+1  .
 *
 *      WHEN NBDCND HAS ANY OTHER VALUE, BDD IS A DUMMY VARIABLE.
 *
 *    ELMBDA
 *      THE CONSTANT LAMBDA IN THE HELMHOLTZ EQUATION.  IF
 *      LAMBDA .GT. 0, A SOLUTION MAY NOT EXIST.  HOWEVER, HWSCRT WILL
 *      ATTEMPT TO FIND A SOLUTION.
 *
 *    F
 *      A TWO-DIMENSIONAL ARRAY WHICH SPECIFIES THE VALUES OF THE RIGHT
 *      SIDE OF THE HELMHOLTZ EQUATION AND BOUNDARY VALUES (IF ANY).
 *      FOR I = 2,3,...,M AND J = 2,3,...,N
 *
 *           F(I,J) = F(X(I),Y(J)).
 *
 *      ON THE BOUNDARIES F IS DEFINED BY
 *
 *           MBDCND     F(1,J)        F(M+1,J)
 *           ------     ---------     --------
 *
 *             0        F(A,Y(J))     F(A,Y(J))
 *             1        U(A,Y(J))     U(B,Y(J))
 *             2        U(A,Y(J))     F(B,Y(J))     J = 1,2,...,N+1
 *             3        F(A,Y(J))     F(B,Y(J))
 *             4        F(A,Y(J))     U(B,Y(J))
 *
 *
 *           NBDCND     F(I,1)        F(I,N+1)
 *           ------     ---------     --------
 *
 *             0        F(X(I),C)     F(X(I),C)
 *             1        U(X(I),C)     U(X(I),D)
 *             2        U(X(I),C)     F(X(I),D)     I = 1,2,...,M+1
 *             3        F(X(I),C)     F(X(I),D)
 *             4        F(X(I),C)     U(X(I),D)
 *
 *      F MUST BE DIMENSIONED AT LEAST (M+1)*(N+1).
 *
 *      NOTE
 *
 *      IF THE TABLE CALLS FOR BOTH THE SOLUTION U AND THE RIGHT SIDE F
 *      AT  A CORNER THEN THE SOLUTION MUST BE SPECIFIED.
 *
 *    IDIMF
 *      THE ROW (OR FIRST) DIMENSION OF THE ARRAY F AS IT APPEARS IN THE
 *      PROGRAM CALLING HWSCRT.  THIS PARAMETER IS USED TO SPECIFY THE
 *      VARIABLE DIMENSION OF F.  IDIMF MUST BE AT LEAST M+1  .
 *
 *    W
 *      A ONE-DIMENSIONAL ARRAY THAT MUST BE PROVIDED BY THE USER FOR
 *      WORK SPACE.  W MAY REQUIRE UP TO 4*(N+1) +
 *      (13 + INT(LOG2(N+1)))*(M+1) LOCATIONS.  THE ACTUAL NUMBER OF
 *      LOCATIONS USED IS COMPUTED BY HWSCRT AND IS RETURNED IN LOCATION
 *      W(1).
 *
 *
 *            * * * * * *   ON OUTPUT     * * * * * *
 *
 *    F
 *      CONTAINS THE SOLUTION U(I,J) OF THE FINITE DIFFERENCE
 *      APPROXIMATION FOR THE GRID POINT (X(I),Y(J)), I = 1,2,...,M+1,
 *      J = 1,2,...,N+1  .
 *
 *    PERTRB
 *      IF A COMBINATION OF PERIODIC OR DERIVATIVE BOUNDARY CONDITIONS
 *      IS SPECIFIED FOR A POISSON EQUATION (LAMBDA = 0), A SOLUTION MAY
 *      NOT EXIST.  PERTRB IS A CONSTANT, CALCULATED AND SUBTRACTED FROM
 *      F, WHICH ENSURES THAT A SOLUTION EXISTS.  HWSCRT THEN COMPUTES
 *      THIS SOLUTION, WHICH IS A LEAST SQUARES SOLUTION TO THE ORIGINAL
 *      APPROXIMATION.  THIS SOLUTION PLUS ANY CONSTANT IS ALSO A
 *      SOLUTION.  HENCE, THE SOLUTION IS NOT UNIQUE.  THE VALUE OF
 *      PERTRB SHOULD BE SMALL COMPARED TO THE RIGHT SIDE F.  OTHERWISE,
 *      A SOLUTION IS OBTAINED TO AN ESSENTIALLY DIFFERENT PROBLEM.
 *      THIS COMPARISON SHOULD ALWAYS BE MADE TO INSURE THAT A
 *      MEANINGFUL SOLUTION HAS BEEN OBTAINED.
 *
 *    IERROR
 *      AN ERROR FLAG THAT INDICATES INVALID INPUT PARAMETERS.  EXCEPT
 *      FOR NUMBERS 0 AND 6, A SOLUTION IS NOT ATTEMPTED.
 *
 *      = 0  NO ERROR.
 *      = 1  A .GE. B.
 *      = 2  MBDCND .LT. 0 OR MBDCND .GT. 4  .
 *      = 3  C .GE. D.
 *      = 4  N .LE. 3
 *      = 5  NBDCND .LT. 0 OR NBDCND .GT. 4  .
 *      = 6  LAMBDA .GT. 0  .
 *      = 7  IDIMF .LT. M+1  .
 *      = 8  M .LE. 3
 *
 *      SINCE THIS IS THE ONLY MEANS OF INDICATING A POSSIBLY INCORRECT
 *      CALL TO HWSCRT, THE USER SHOULD TEST IERROR AFTER THE CALL.
 *
 *    W
 *      W(1) CONTAINS THE REQUIRED LENGTH OF W.
 *
 *
 *    @section spec_hwscrt  PROGRAM SPECIFICATIONS 
 *
 *
 *    DIMENSION OF   BDA(N+1),BDB(N+1),BDC(M+1),BDD(M+1),F(IDIMF,N+1),
 *    ARGUMENTS      W(SEE ARGUMENT LIST)
 *
 *    LATEST         JUNE 1, 1976
 *    REVISION
 *
 *    SUBPROGRAMS    HWSCRT,GENBUN,POISD2,POISN2,POISP2,COSGEN,MERGE,
 *    REQUIRED       TRIX,TRI3,PIMACH
 *
 *    SPECIAL        NONE
 *    CONDITIONS
 *
 *    COMMON         NONE
 *    BLOCKS
 *
 *    I/O            NONE
 *
 *    PRECISION      SINGLE
 *
 *    SPECIALIST     ROLAND SWEET
 *
 *    LANGUAGE       FORTRAN
 *
 *    HISTORY        STANDARDIZED SEPTEMBER 1, 1973
 *                   REVISED APRIL 1, 1976
 *
 *    ALGORITHM      THE ROUTINE DEFINES THE FINITE DIFFERENCE
 *                   EQUATIONS, INCORPORATES BOUNDARY DATA, AND ADJUSTS
 *                   THE RIGHT SIDE OF SINGULAR SYSTEMS AND THEN CALLS
 *                   GENBUN TO SOLVE THE SYSTEM.
 *
 *    SPACE          13110(OCTAL) = 5704(DECIMAL) LOCATIONS ON THE NCAR
 *    REQUIRED       CONTROL DATA 7600
 *
 *    TIMING AND        THE EXECUTION TIME T ON THE NCAR CONTROL DATA
 *    ACCURACY       7600 FOR SUBROUTINE HWSCRT IS ROUGHLY PROPORTIONAL
 *                   TO M*N*LOG2(N), BUT ALSO DEPENDS ON THE INPUT
 *                   PARAMETERS NBDCND AND MBDCND.  SOME TYPICAL VALUES
 *                   ARE LISTED IN THE TABLE BELOW.
 *                      THE SOLUTION PROCESS EMPLOYED RESULTS IN A LOSS
 *                   OF NO MORE THAN THREE SIGNIFICANT DIGITS FOR N AND
 *                   M AS LARGE AS 64.  MORE DETAILED INFORMATION ABOUT
 *                   ACCURACY CAN BE FOUND IN THE DOCUMENTATION FOR
 *                   SUBROUTINE GENBUN WHICH IS THE ROUTINE THAT
 *                   SOLVES THE FINITE DIFFERENCE EQUATIONS.
 *
 *
 *                      M(=N)    MBDCND    NBDCND    T(MSECS)
 *                      -----    ------    ------    --------
 *
 *                       32        0         0          31
 *                       32        1         1          23
 *                       32        3         3          36
 *                       64        0         0         128
 *                       64        1         1          96
 *                       64        3         3         142
 *
 *    PORTABILITY    AMERICAN NATIONAL STANDARDS INSTITUTE FORTRAN.
 *                   ALL MACHINE DEPENDENT CONSTANTS ARE LOCATED IN THE
 *                   FUNCTION PIMACH.
 *
 *    REFERENCE      SWARZTRAUBER,P. AND R. SWEET, 'EFFICIENT FORTRAN
 *                   SUBPROGRAMS FOR THE SOLUTION OF ELLIPTIC EQUATIONS'
 *                   NCAR TN/IA-109, JULY, 1975, 138 PP.
 *
 */
  
