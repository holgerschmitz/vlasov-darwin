// -*- C++ -*-
// $Id$
/** @file boundary.h
    * @brief Interface for wrapping and exchanging boundaries
    *
    *  This interface is used to exchange the boundaries of distribution
 	*  functions and scalar fields. It can be implemented to define
 	*  periodic boundaries or exchange data with other processes.
 	*
 	*  See the page @ref indices for a discussion on the numerical
 	*  ranges of the fields
	*/
//-----------------------------------------------------------------------------
#ifndef BOUNDARY_H
#define BOUNDARY_H

#include "numeric.h"
#include "vlasov.h"
#include "scalarfield.h"
#include "parameter.h"
#include "task.h"
#include "numboundary.h"
//-----------------------------------------------------------------------------
//if this is a mutli processor run, use MPI for interprocessor communicatio
#ifndef SINGLE_PROCESSOR
#include <mpi.h>
#endif
//-----------------------------------------------------------------------------
//Boundary

/** @brief Interface for wrapping and exchanging boundaries .
 *
 *  This interface is used to exchange the boundaries of distribution
 *  functions and scalar fields. It can be implemented to define
 *  periodic boundaries or exchange data with other processes.
 *  This is the (abstract) base class of all boundary classes.
 *
 *  See the page @ref indices for a discussion on the numerical
 *  ranges of the fields
 */
class Boundary : public Rebuildable {
  public:
	/// Default constructor, 
      Boundary() {}
      /** @brief Virtual destructor
	 *  
	 *  We need a virtual destructor because the class has 
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

      /** @brief Exchange the boundaries of the distribution function
       *  in the x-direction
       */
      virtual void exchangeXLine(NumMatrix<double, 1> &field) = 0;

      /** @brief Exchange the boundaries of the distribution function
       *  in the y-direction
       */
      virtual void exchangeYLine(NumMatrix<double, 1> &field) = 0;

      /** @brief Combine the scalar field. 
       *
       *  The resulting field should be the sum of the fields of all 
       *  the processes. Additional wrapping may be done.
       */
      virtual void ScalarFieldCombine(ScalarField &field) const = 0;

      /** @brief Reduce the scalar field. 
       *
       *  Assume the scalar field is present on all processes.
       *  Only perform wrapping.
       */
      virtual void ScalarFieldReduce(ScalarField &field) const = 0;
      
      /** Returns the numeric boundary for a scalar field. */
      virtual const NumBoundary& getNumBoundary(ScalarField &field) const = 0;

      /// Return the average of a single value over all the processes
      virtual double AvgReduce(double) = 0;    

      /// Return the maximum of a single value over all the processes
      virtual double MaxReduce(double) = 0;    

      /// Return the lower bound of the distribution function
      virtual const PhasePositionI &DistLow() const = 0;

      /// Return the upper bound of the distribution function
      virtual const PhasePositionI &DistHigh() const = 0;

      /// Return the lower bound of scalar fields
      virtual const PositionI &scalarLow() const = 0;

      /// Return the upper bound of scalar fields
      virtual const PositionI &scalarHigh() const = 0;

      /// Return true if this is the master process and false otherwise
      virtual bool master() const = 0;

      /// Return the process number
      virtual int procnum() const = 0;
      
	  ///get a unique Id
      virtual int getUniqueId() const = 0;
};
//Boundary
//-----------------------------------------------------------------------------
//SinglePeriodicBoundary

/** @brief Implements Boudary to supply a periodic system running on
 *  a single processor
 */
class SinglePeriodicBoundary : public Boundary {
	private:
	///a periodic and symmetric boundary
      PPBoundary numBound;
  public:   
      /// Wraps the boundaries in x-direction
      void exchangeX(VlasovDist &field);
      /// Wraps the boundaries in y-direction
      void exchangeY(VlasovDist &field);
      /// Wraps the boundaries in x-direction
      virtual void exchangeXLine(NumMatrix<double, 1> &field);
      /// Wraps the boundaries in y-direction
      virtual void exchangeYLine(NumMatrix<double, 1> &field);
      /// Only wraps the boundaries of the scalar field
      void ScalarFieldCombine(ScalarField &field) const;
      /// Only wraps the boundaries of the scalar field
      void ScalarFieldReduce(ScalarField &field) const;
      /// Returns periodic boundary conditions
      const NumBoundary& getNumBoundary(ScalarField &field) const;
      /// There is no average to be calculated
      double AvgReduce(double val) { return val; }
      /// There is no maximum to be calculated
      double MaxReduce(double val) { return val; }
      /// Returns the global lower bound of the distribution function
      const PhasePositionI &DistLow() const;
      /// Returns the global upper bound of the distribution function
      const PhasePositionI &DistHigh() const;
      /// Return the lower bound of scalar fields
      const PositionI &scalarLow()const;
      /// Return the upper bound of scalar fields
      const PositionI &scalarHigh()const;
      /// There is only one process, so master always returns true
      bool master() const { return true; }
      /// The process number is always zero
      int procnum() const { return 0; }
      /// The unique id number is always zero
      int getUniqueId() const { return 0; }
};
//SinglePeriodicBoundary
//-----------------------------------------------------------------------------
//MPIPeriodicSplitBoundary
#ifndef SINGLE_PROCESSOR
// this will be used for multi processor runs, to split the simulation grid 
// into vertical rectangles

/** @brief a boundary class for mutliple processor runs
  *
  * Is design to be exchanged via the MPI protocol
  * This implementation splits only the x-axis into rectangles.
  */
class MPIPeriodicSplitXBoundary : public Boundary {
  protected:
	///a periodic and symmetric boundary
      PPBoundary numBound;

      /// The number of processes
      int ComSize;
  
      /// The rank of the current process
      int ComRank;
      
      /// The Comm object referring to the cartesian process grid
      MPI_Comm comm;
      
      /// The coordinates of this process
      int mycoord;
      
      
	int leftcoord,  ///< The rank of the left neighbour process
	    rightcoord; ///< The rank of the right  neighbour process
      
      /** @brief The size of the array that needs to be exchanged, 
       *  when the exchangeX method is called
       */
      int exchSize;
      
	double *sendarr, ///< buffer holding the data to be send (size: exchSize)
		 *recvarr; ///< buffer holding the received data (size: exchSize)

      /// The size of the scalar fields when reducing
      int scalarSize;
      
	double *scalarsend, ///< send Buffer for the data used to reduce the scalarfields (size: scalarSize)
		 *scalarrecv; ///< receive Buffer for the data used to reduce the scalarfields (size: scalarSize)

PhasePositionI Low,///<The position of the lower corner of the local piece of the distribution function in phase space
	        High;///<The position of the lower corner of the local piece of the distribution function in phase space
	
	///initialize
      void init(int argc, char **argv);
  public:
	///default constructor
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
      
      /// Wraps the boundaries in x-direction
      virtual void exchangeXLine(NumMatrix<double, 1> &field);
      /// Wraps the boundaries in y-direction
      virtual void exchangeYLine(NumMatrix<double, 1> &field);

      /// Adds the scalar fields and wraps them
      void ScalarFieldCombine(ScalarField &field) const;
      /// Wraps the scalar fields
      void ScalarFieldReduce(ScalarField &field) const;
      
      /// Returns periodic boundary conditions
      const NumBoundary& getNumBoundary(ScalarField &field) const;

      /** @brief Use MPIALLReduce to calculate the sum and then divide
       *  by the number of processes.
       */
      double AvgReduce(double val);
      
      /** @brief Use MPIALLReduce to calculate the maximum
       */
      double MaxReduce(double val);
      
      /// Returns the global lower bound of the distribution function
      const PhasePositionI &DistLow()const;
      
      /// Returns the global upper bound of the distribution function
      const PhasePositionI &DistHigh()const;
      
      /// Return the lower bound of scalar fields
      const PositionI &scalarLow()const;
      
      /// Return the upper bound of scalar fields
      const PositionI &scalarHigh()const;

      /// The process with the rank zero is designated master process
      bool master() const { return ComRank==0; }
      
      /// Returns the comm rank as given by mpi
      int procnum() const { return ComRank; }
	///returns an ID, which is identical with the coordinates
      int getUniqueId() const { return mycoord; }
};
//MPIPeriodicSplitXBoundary
//-----------------------------------------------------------------------------
//MPIPeriodicSplitXYBoundary

/** @brief a boundary class for mutliple processor runs
 *
 * Is designed to be exchanged via the MPI protocol.
 * Here splitting is performed in both spatial directions.
 */
class MPIPeriodicSplitXYBoundary : public Boundary {
  protected:
	///a periodic and symmetric boundary
      PPBoundary numBound;

      /// The number of processes
      int ComSize;
  
      /// The rank of the current process
      int ComRank;
      
      /// The Comm object referring to the cartesian process grid
      MPI_Comm comm;

      
	int leftcoord,  ///< The rank of the left neighbour process
	    rightcoord; ///< The rank of the right neighbour process
      
	int topcoord, ///< The rank of the top neighbour process
	    bottomcoord;  ///< The rank of the bottom neighbour processes
	///dimensions
      int dims[2];
      /// The cartesian coordinates of this process
      int mycoord[2];

      /** @brief The size of the array that needs to be exchanged, 
       *  when the exchangeX or the exchangeY method is called
       */
      int exchSize[2];
      
	double *sendarrx, ///< send Buffer for exchanging data in x-direction (size: exchSize[0])
	       *recvarrx; ///< receive Buffer for exchanging data in x-direction (size: exchSize[0])
      
	double *sendarry, ///< send Buffer for exchanging data in y-direction (size: exchSize[1])
	       *recvarry; ///< receive Buffer for exchanging data in y-direction (size: exchSize[1])

      /// The size of the scalar fields when reducing
      int scalarSize;
      
      
	double *scalarsend, ///< send Buffer for exchanging data used for reducing scalar fields (size: scalarSize)
	       *scalarrecv; ///< receive Buffer for exchanging data used for reducing scalar fields (size: scalarSize)

      PhasePositionI Low, 
	///<The positions of the lower corner of the local piece of the distribution function in phase space
			   High;
	///<The positions of the upper corner of the local piece of the distribution function in phase space
			   
	///initialize
      void init(int argc, char **argv);
  public:
	///default constructor
      MPIPeriodicSplitXYBoundary();
	
     /** @brief Constructor using the command line arguments.
	 *
       *  Sets up all the local variables and allocates memory
       *  for the buffers
       */
      MPIPeriodicSplitXYBoundary(int argc, char **argv);
      
      
      /// Virtual destructor deleting all the allocated arrays
      ~MPIPeriodicSplitXYBoundary();
      
     /** @brief Exchanges the boundaries in x-direction.
       *
       *  The two outmost simulated cells are sent and the surrounding 
       *  two ghost cells are filled with values
       */
      void exchangeX(VlasovDist &field);
      
     /** @brief Exchanges the boundaries in y-direction.
	 *
       *  The two outmost simulated cells are sent and the surrounding 
       *  two ghost cells are filled with values
       */
      void exchangeY(VlasovDist &field);
      
      /// Wraps the boundaries in x-direction
      virtual void exchangeXLine(NumMatrix<double, 1> &field);
      /// Wraps the boundaries in y-direction
      virtual void exchangeYLine(NumMatrix<double, 1> &field);

      /// Adds the scalar fields and wraps them
      void ScalarFieldCombine(ScalarField &field) const;
      /// Wraps the scalar fields
      void ScalarFieldReduce(ScalarField &field) const;
      
      /// Returns periodic boundary conditions
      const NumBoundary& getNumBoundary(ScalarField &field) const;

      /// Use MPIALLReduce to calculate the sum and then divide by the number of processes. 
      double AvgReduce(double val);
      
      /// Use MPIALLReduce to calculate the maximum
      double MaxReduce(double val);
      
      /// Returns the global lower bound of the distribution function
      const PhasePositionI &DistLow()const;
      
      /// Returns the global upper bound of the distribution function
      const PhasePositionI &DistHigh()const;
      
      /// Return the lower bound of scalar fields
      const PositionI &scalarLow()const;
      
      /// Return the upper bound of scalar fields
      const PositionI &scalarHigh()const;
      
      /// The process with the rank zero is designated master process
      bool master() const { return ComRank==0; }
      
      /// Returns the comm rank as given by mpi
      int procnum() const { return ComRank; }
	///returns an ID, which consists of the Dimensions and coordinates
      int getUniqueId() const { return dims[1]*mycoord[0] + mycoord[1]; }
	/// returns "true" since this is periodic in the x-direction
      virtual bool periodicX() { return true; } 
	/// returns "true" since this is periodic in the y-direction
      virtual bool periodicY() { return true; }
};
//MPIPeriodicSplitXYBoundary
//-----------------------------------------------------------------------------
#endif // multiple processor

#endif //BOUNDARY_H
