#ifndef BOUNDARY_H
#define BOUNDARY_H

#include "config.h"
#include "numeric.h"
#include <mpi.h>

/** @brief Interface for wrapping and exchanging boundaries.
 *
 *  This interface is used to exchange the boundaries of distribution
 *  functions and scalar fields. It can be implemented to define
 *  periodic boundaries or exchange data with other processes.
 *
 *  See the page @ref indices for a discussion on the numerical
 *  ranges of the fields
 */
class Boundary {
    public:
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
    virtual void ScalarFieldReduce(ScalarField &field) = 0;
    
    /// Return the sum of all the processes
    virtual double AvgReduce(double) = 0;    
    
    /// Return the lower bound of the distribution function
    virtual PhasePositionI &DistLow() = 0;
    
    /// Return the upper bound of the distribution function
    virtual PhasePositionI &DistHigh() = 0;
    
    virtual bool master() = 0;
    virtual int procnum() = 0;
};


/** @brief Implements Boudary to supply a periodic system running on
 *  a single processor
 */
class SinglePeriodicBoundary : public Boundary {
    public:
    /// Wraps the boundaries in x-direction
    void exchangeX(VlasovDist &field);
    /// Wraps the boundaries in y-direction
    void exchangeY(VlasovDist &field);
    /// Only wraps the boundaries of the scalar field
    void ScalarFieldReduce(ScalarField &field);
    /// There is no sum to be calculated
    double AvgReduce(double val) { return val; }
    /// Returns the global lower bound of the distribution function
    PhasePositionI &DistLow();
    /// Returns the global upper bound of the distribution function
    PhasePositionI &DistHigh();
    bool master() { return true; }
    int procnum() { return 0; }
};

/** @brief Implements Boudary to supply a periodic system running on
 *  a single processor
 */
class MPIPeriodicSplitXBoundary : public Boundary {
    private:
    int ComSize;
    int ComRank;
    MPI_Comm comm;
    int leftcoord, mycoord, rightcoord;

    int exchSize;
    double *sendarr, *recvarr;

    int scalarSize;
    double *scalarsend, *scalarrecv;

    PhasePositionI Low, High;
    public:
    MPIPeriodicSplitXBoundary(int argc, char **argv);
    ~MPIPeriodicSplitXBoundary();
    /// Exchanges the boundaries in x-direction
    void exchangeX(VlasovDist &field);
    /// Wraps the boundaries in y-direction
    void exchangeY(VlasovDist &field);
    /// Adds the scalar fields and wraps them
    void ScalarFieldReduce(ScalarField &field);
    /// Use MPIALLReduce to calculate the sum
    double AvgReduce(double val);
    /// Returns the global lower bound of the distribution function
    virtual PhasePositionI &DistLow();
    /// Returns the global upper bound of the distribution function
    virtual PhasePositionI &DistHigh();
    bool master() { return ComRank==0; }
    int procnum() { return ComRank; }
};

/** @brief Implements Boundary to supply a periodic system running on
 *  a single processor
 */
class MPIPeriodicSplitXYBoundary : public Boundary {
    private:
    int ComSize;
    int ComRank;
    MPI_Comm comm;
    int leftcoord, rightcoord;
    int topcoord, bottomcoord;
    int mycoord[2];

    int exchSize[2];
    double *sendarrx, *recvarrx;
    double *sendarry, *recvarry;

    int scalarSize;
    double *scalarsend, *scalarrecv;

    PhasePositionI Low, High;
    public:
    MPIPeriodicSplitXYBoundary(int argc, char **argv);
    ~MPIPeriodicSplitXYBoundary();
    /// Exchanges the boundaries in x-direction
    void exchangeX(VlasovDist &field);
    /// Wraps the boundaries in y-direction
    void exchangeY(VlasovDist &field);
    /// Adds the scalar fields and wraps them
    void ScalarFieldReduce(ScalarField &field);
    /// Use MPIALLReduce to calculate the sum
    double AvgReduce(double val);
    /// Returns the global lower bound of the distribution function
    virtual PhasePositionI &DistLow();
    /// Returns the global upper bound of the distribution function
    virtual PhasePositionI &DistHigh();
    bool master() { return ComRank==0; }
    int procnum() { return ComRank; }
};

#endif
