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
      /** @brief An enumeration that contains the type of boundary conditions.
       */
      enum bcond { periodic=0, fixfix=1, fixder=2, derfix=3, derder=4 };
  
      /// Default constructor
      Poisson();
      
      /** @brief Construct using the physical extent, the number of
       *  steps in each direction and the boundary conditions
       */
      Poisson(const PositionD &pMiN_, const PositionD &pMaX_, 
              const PositionI &steps_,
              bcond boundx, bcond boundy);
      
      /// Destruct and free the WorkArray buffer
      ~Poisson();     
      
      /** @brief Resize the physical extent, the number of
       *  steps in each direction and the boundary conditions
       */
      void resize(const PositionD &pMiN_, const PositionD &pMaX_, 
                  const PositionI &steps_,
                  bcond boundx, bcond boundy);
      
      /** @brief Resize the physical extent, the number of
       *  steps in each direction. Don't change the boundary
       *  conditions
       */
      void resize(const PositionD &pMiN_, const PositionD &pMaX_, 
                  const PositionI &steps_);
      
      /// Set new boundary conditions
      void BCond(bcond boundx, bcond boundy);
      
      /** @brief Solve the Poisson equation.
       *  The result is returned in the second parameter. Internally the
       *  solve_step method is called twice to promote the second
       *  order hwscrt to a fourth order solver
       */
      void solve(ScalarField &In, ScalarField &Out);
  private:
      
      /// Perform a single call to hwscrt and perform wrapping
      void solve_step(ScalarField &In, ScalarField &Out);
  
      /** @brief Perform resizing. Set the parameters and allocate
       *  memory for the WorkArray
       */
      void do_resize(const PositionD &pMiN_, const PositionD &pMaX_, 
                     const PositionI &steps_);       
};

#endif
