// -*- C++ -*-
// $Id$

#ifndef BOUNDARY_H
#define BOUNDARY_H

#include "numeric.h"
#include "vlasov.h"
#include "scalarfield.h"
#include "parameter.h"
#include "task.h"
#include "numboundary.h"

#ifndef SINGLE_PROCESSOR
#include <mpi.h>
#endif

/** @brief Interface for wrapping and exchanging boundaries.
 *
 *  This interface is used to exchange the boundaries of distribution
 *  functions and scalar fields. It can be implemented to define
 *  periodic boundaries or exchange data with other processes.
 *
 *  See the page @ref indices for a discussion on the numerical
 *  ranges of the fields
 */
class Boundary : public Rebuildable {
  public:
      Boundary() {}
      /** @brief We need a virtual destructor because the class has 
       *  virtual methods
       */
      virtual ~Boundary() {}

      /** @brief Exchange the boundaries of the distribution function
       *  in the x-direction
       */
      virtual void exchangeX(VlasovDist &field) = 0;

      /** @brief Exchange the boundaries of the distribution function
       *  in the y-direction
       */
      virtual void exchangeY(VlasovDist &field) = 0;

      /** @brief Reduce the scalar field. 
       *
       *  The resulting field should be the sum of the fields of all 
       *  the processes. Additional wrapping may be done.
       */
      virtual void ScalarFieldReduce(ScalarField &field) const = 0;
      
      virtual const NumBoundary& getNumBoundary(ScalarField &field) const = 0;

      /// Return the average of a single value over all the processes
      virtual double AvgReduce(double) = 0;    

      /// Return the lower bound of the distribution function
      virtual const PhasePositionI &DistLow() const = 0;

      /// Return the upper bound of the distribution function
      virtual const PhasePositionI &DistHigh() const = 0;

      /// Return the lower bound of scalar fields
      virtual const PositionI &scalarLow() const = 0;

      /// Return the upper bound of scalar fields
      virtual const PositionI &scalarHigh() const = 0;

      /// Return true if this is the master process and false otherwise
      virtual bool master() = 0;

      /// Return the process number
      virtual int procnum() const = 0;
      
      virtual int getUniqueId() const = 0;
};


/** @brief Implements Boudary to supply a periodic system running on
 *  a single processor
 */
class SinglePeriodicBoundary : public Boundary {
  private:
      PPBoundary numBound;
  public:   
      /// Wraps the boundaries in x-direction
      void exchangeX(VlasovDist &field);
      /// Wraps the boundaries in y-direction
      void exchangeY(VlasovDist &field);
      /// Only wraps the boundaries of the scalar field
      void ScalarFieldReduce(ScalarField &field) const;
      /// Returns periodic boundary conditions
      const NumBoundary& getNumBoundary(ScalarField &field) const;
      /// There is no average to be calculated
      double AvgReduce(double val) { return val; }
      /// Returns the global lower bound of the distribution function
      const PhasePositionI &DistLow() const;
      /// Returns the global upper bound of the distribution function
      const PhasePositionI &DistHigh() const;
      /// Return the lower bound of scalar fields
      const PositionI &scalarLow()const;
      /// Return the upper bound of scalar fields
      const PositionI &scalarHigh()const;
      /// There is only one process, so master always returns true
      bool master() { return true; }
      /// The process number is always zero
      int procnum() const { return 0; }
      /// The unique id number is always zero
      int getUniqueId() const { return 0; }
};

#ifndef SINGLE_PROCESSOR

/** @brief Implements Boudary to supply a periodic system running on
 *  a single processor
 */
class MPIPeriodicSplitXBoundary : public Boundary {
  protected:
      PPBoundary numBound;

      /// The number of processes
      int ComSize;
  
      /// The rank of the current process
      int ComRank;
      
      /// The Comm object referring to the cartesian process grid
      MPI_Comm comm;
      
      /// The coordinates of this process
      int mycoord;
      
      /// The ranks of the left and the right neighbour processes
      int leftcoord, rightcoord;
      
      /** @brief The size of the array that needs to be exchanged, 
       *  when the exchangeX method is called
       */
      int exchSize;
      
      /** @brief Buffers of the size exchSize that hold the data
       *  to send and receive
       */
      double *sendarr, *recvarr;

      /// The size of the scalar fields when reducing
      int scalarSize;
      
      /** @brief Buffers of the size scalarSize that hold the data
       *  to reduce the scalar fields
       */
      double *scalarsend, *scalarrecv;

      /** @brief The positions of the lower and upper corner of the
       *  local piece of the distribution function in phase space
       */
      PhasePositionI Low, High;
      void init(int argc, char **argv);
  public:
      MPIPeriodicSplitXBoundary();
      /** @brief Constructor using the command line arguments.
       *  Sets up all the local variables and allocates memory
       *  for the buffers
       */
      MPIPeriodicSplitXBoundary(int argc, char **argv);
  
      /// Virtual destructor deleting all the allocated arrays
      ~MPIPeriodicSplitXBoundary();
      
      /** @brief Exchanges the boundaries in x-direction.
       *  The two outmost simulated cells are sent and the surrounding 
       *  two ghost cells are filled with values
       */
      void exchangeX(VlasovDist &field);
      
      /// Wraps the boundaries in y-direction
      void exchangeY(VlasovDist &field);
      
      /// Adds the scalar fields and wraps them
      void ScalarFieldReduce(ScalarField &field) const;
      
      /// Returns periodic boundary conditions
      const NumBoundary& getNumBoundary(ScalarField &field) const;

      /** @brief Use MPIALLReduce to calculate the sum and then divide
       *  by the number of processes.
       */
      double AvgReduce(double val);
      
      /// Returns the global lower bound of the distribution function
      const PhasePositionI &DistLow()const;
      
      /// Returns the global upper bound of the distribution function
      const PhasePositionI &DistHigh()const;
      
      /// Return the lower bound of scalar fields
      const PositionI &scalarLow()const;
      
      /// Return the upper bound of scalar fields
      const PositionI &scalarHigh()const;

      /// The process with the rank zero is designated master process
      bool master() { return ComRank==0; }
      
      /// Returns the comm rank as given by mpi
      int procnum() const { return ComRank; }
      
      int getUniqueId() const { return mycoord; }
};

/** @brief Implements Boundary to supply a periodic system running on
 *  a single processor
 */
class MPIPeriodicSplitXYBoundary : public Boundary {
  private:
      PPBoundary numBound;

      /// The number of processes
      int ComSize;
  
      /// The rank of the current process
      int ComRank;
      
      /// The Comm object referring to the cartesian process grid
      MPI_Comm comm;

      /// The ranks of the left and the right neighbour processes
      int leftcoord, rightcoord;
      /// The ranks of the top and the bottom neighbour processes
      int topcoord, bottomcoord;
      
      int dims[2];
      /// The cartesian coordinates of this process
      int mycoord[2];

      /** @brief The size of the array that needs to be exchanged, 
       *  when the exchangeX or the exchangeY method is called
       */
      int exchSize[2];
      
      /** @brief Buffers of the size exchSize[0] that hold the data
       *  to send and receive in the x--direction
       */
      double *sendarrx, *recvarrx;
      
      /** @brief Buffers of the size exchSize[1] that hold the data
       *  to send and receive in the y--direction
       */
      double *sendarry, *recvarry;

      /// The size of the scalar fields when reducing
      int scalarSize;
      
      /** @brief Buffers of the size scalarSize that hold the data
       *  to reduce the scalar fields
       */
      double *scalarsend, *scalarrecv;

      /** @brief The positions of the lower and upper corner of the
       *  local piece of the distribution function in phase space
       */
      PhasePositionI Low, High;
      void init(int argc, char **argv);
  public:
      MPIPeriodicSplitXYBoundary();
      /** @brief Constructor using the command line arguments.
       *  Sets up all the local variables and allocates memory
       *  for the buffers
       */
      MPIPeriodicSplitXYBoundary(int argc, char **argv);
      
      
      /// Virtual destructor deleting all the allocated arrays
      ~MPIPeriodicSplitXYBoundary();
      
      /** @brief Exchanges the boundaries in x-direction.
       *  The two outmost simulated cells are sent and the surrounding 
       *  two ghost cells are filled with values
       */
      void exchangeX(VlasovDist &field);
      
      /** @brief Exchanges the boundaries in y-direction.
       *  The two outmost simulated cells are sent and the surrounding 
       *  two ghost cells are filled with values
       */
      void exchangeY(VlasovDist &field);
      
      /// Adds the scalar fields and wraps them
      void ScalarFieldReduce(ScalarField &field) const;
      
      /// Returns periodic boundary conditions
      const NumBoundary& getNumBoundary(ScalarField &field) const;

      /** @brief Use MPIALLReduce to calculate the sum and then divide
       *  by the number of processes.
       */
      double AvgReduce(double val);
      
      /// Returns the global lower bound of the distribution function
      const PhasePositionI &DistLow()const;
      
      /// Returns the global upper bound of the distribution function
      const PhasePositionI &DistHigh()const;
      
      /// Return the lower bound of scalar fields
      const PositionI &scalarLow()const;
      
      /// Return the upper bound of scalar fields
      const PositionI &scalarHigh()const;
      
      /// The process with the rank zero is designated master process
      bool master() { return ComRank==0; }
      
      /// Returns the comm rank as given by mpi
      int procnum() const { return ComRank; }

      int getUniqueId() const { return dims[1]*mycoord[0] + mycoord[1]; }
};
#endif // single processor


#endif
