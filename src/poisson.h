#ifndef POISSON_H
#define POISSON_H

#include <matrix.H>
#include "index.h"
#include "vlasov.h"


class Poisson {
    private:
        PositionD pMiN, pMaX;
        PositionI steps;
        PositionI boundary;
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
        void do_resize(const PositionD &pMiN_, const PositionD &pMaX_, 
                       const PositionI &steps_);       
};

#endif
