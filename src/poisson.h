// -*- C++ -*-
// $Id$

#ifndef POISSON_H
#define POISSON_H

#include <matrix.H>
#include "index.h"
#include "vlasov.h"

/** @brief Wraps the Fortran hwscrt routine in a handy class that
 *  solves poissons equation.
 *
 *  At the moment only periodic boundary conditions are implemented.
 *
 *  See the page @ref hwscrt
 */
class Poisson {
    private:
        /// The physical coordinates of the minimum and maximum points
        PositionD pMiN, pMaX;
        /// The type of boundary condition for the two directions
        PositionI boundary;
        /// A temporary working array needed by hwscrt
        double *WorkArray;
    public:
        enum bcond { periodic=0, fixfix=1, fixder=2, derfix=3, derder=4 };
        Poisson();
        Poisson(const PositionD &pMiN_, const PositionD &pMaX_, 
                const PositionI &steps_,
                bcond boundx, bcond boundy);
        ~Poisson();
        void resize(const PositionD &pMiN_, const PositionD &pMaX_, 
                    const PositionI &steps_,
                    bcond boundx, bcond boundy);
        void resize(const PositionD &pMiN_, const PositionD &pMaX_, 
                    const PositionI &steps_);
        void BCond(bcond boundx, bcond boundy);
        void solve(ScalarField &In, ScalarField &Out);
    private:
        void solve_step(ScalarField &In, ScalarField &Out);
        void do_resize(const PositionD &pMiN_, const PositionD &pMaX_, 
                       const PositionI &steps_);       
};

#endif
