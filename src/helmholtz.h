// -*- C++ -*-
// $Id$

#include "matrix.h"

/** @brief Solves Helmholtz equation using a multigrid method.
 *
 *  The Helmholtz equation is given by
 *  \f$\Delta u(x) -\lambda(x)u(x) = f(x)\f$
 *  where \f$u(x)\f$ is the unknown scalar field.
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
        /** @brief The default constructor reads the attributes from
         *  global variables
         */
        Helmholtz();
    
        /** @brief Solves the Helmholtz equation
         *  \f$\Delta u(x) -\lambda(x)u(x) = f(x)\f$.
         *  The solution is returned in the first parameter.
         */
        void solve(   NumMatrix<double,2> &u,
	                  NumMatrix<double,2> &f,
                      NumMatrix<double,2> &lambda);

        /// The grid spacing in x-direction
        double Dx() { 
            return dx[0];
        }
        
        /// The grid spacing in y-direction
        double Dy() { 
            return dx[1];
        }
    private:
        /// Perform one Gauss-Seidel red-black iteration
        void gauss( NumMatrix<double,2> &u,     
	                NumMatrix<double,2> &f,
                    NumMatrix<double,2> &lambda);   
                     
        /** @brief Calculates the defect and stores the new
         *  fields of \f$f(x)\f$ and \f$\lambda(x)\f$ on the
         *  coarser grid
         */
        void defect(NumMatrix<double,2> &u,
                    NumMatrix<double,2> &f,
                    NumMatrix<double,2> &lambda,
                    NumMatrix<double,2> &fn,
                    NumMatrix<double,2> &lambdan);
                        
        
        /** @brief Prolongates the solution \f$u(x)\f$ onto 
         *  the finer grid resolution.
         */                              
        void prolongate(NumMatrix<double,2> &u,
	                    NumMatrix<double,2> &un);
        
        /** @brief One single multigrid step.
         *  This method calls itself recursively on coarser
         *  and coarser grids
         */
        void mgi(   NumMatrix<double,2> &u,     
	                NumMatrix<double,2> &f,
                    NumMatrix<double,2> &lambda);    
        
        /** @brief Normalizes a scalar field.
         *  This method is not needed for the Helmhotz solver.
         */            
        void normalize( NumMatrix<double,2> &f);
        
        /// Calculates the maximum norm between two scalar fields
        double distance(NumMatrix<double,2> &u,
	                    NumMatrix<double,2> &f);
        
        /** @brief Performs periodic wrapping of the boundaries
         *  @todo This is geometry dependent and should be supplied
         *  by a special boundary object
         */
        void boundary(NumMatrix<double,2> &u);
};
