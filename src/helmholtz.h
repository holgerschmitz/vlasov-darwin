// -*- C++ -*-
// $Id$

#include <matrix.H>

/** @brief Solves Helmholtz equation using a multigrid method.
 *
 *  Assumes that arrays range from 0 to N+1, where 0 and N+1 are
 *  periodic boundaries.
 */
class Helmholtz {
    protected:
        /// physical dimensions of a grid cell. Read from global variables
        double dx[2];

        /** @brief Do nu1 Gauss steps then repeat multigrid refinement
         *  gama times then do nu2 Gauss steps.
         */
        int gama,nu1,nu2;
        
        /// The error to reach.
        double epsilon;
    public:
        Helmholtz();
        void solve(   NumMatrix<double,2> &u,
	                  NumMatrix<double,2> &f,
                      NumMatrix<double,2> &lambda);

        double Dx() { 
            return dx[0];
        }
        double Dy() { 
            return dx[1];
        }
    private:
        void gauss( NumMatrix<double,2> &u,     
	                NumMatrix<double,2> &f,
                    NumMatrix<double,2> &lambda);   
                     
        void defect(NumMatrix<double,2> &u,
                    NumMatrix<double,2> &f,
                    NumMatrix<double,2> &lambda,
                    NumMatrix<double,2> &fn,
                    NumMatrix<double,2> &lambdan);
                                      
        void prolongate(NumMatrix<double,2> &u,
	                    NumMatrix<double,2> &un);
        
        void mgi(   NumMatrix<double,2> &u,     
	                NumMatrix<double,2> &f,
                    NumMatrix<double,2> &lambda);    
                    
        void normalize( NumMatrix<double,2> &f);
                
        double distance(NumMatrix<double,2> &u,
	                    NumMatrix<double,2> &f);
                        
        void boundary(NumMatrix<double,2> &u);
};
