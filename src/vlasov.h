// -*- C++ -*-
// $Id$
//-----------------------------------------------------------------------------
/** @file vlasov.h
 * @brief main doc file, typedefs
 *
 * Here the main documentation is found, along with global typedefs for classes 
 * helping in solving the equation.
 */
 //-----------------------------------------------------------------------------
#include "index.h"
#include "matrix.h"
//-----------------------------------------------------------------------------
/**
 *  @todo Test gyro-frequency (how? waves?)
 *  @todo Test magnetosonic modes
 *  @todo Update Documentation
 */
//-----------------------------------------------------------------------------
#ifndef VLASOV_H
#define VLASOV_H

///symbolic constant Pi
#define PI 3.1415926535897932384626433832795029

///dimension of the spatial grid
static const int DIMENSION=2;

///grid position
typedef FixedArray<int,DIMENSION> PositionI;
///position in real numbers
typedef FixedArray<double,DIMENSION> PositionD;

///position in phasespace, measured in grid cells 
typedef FixedArray<int,DIMENSION+3> PhasePositionI;
///position in phasespace, measured in real numbers
typedef FixedArray<double,DIMENSION+3> PhasePositionD;

///velocity, measured in grid cell per time step
typedef FixedArray<int,3> VelocityI;
///velocity, measured in real numbers
typedef FixedArray<double,3> VelocityD;

///distribution in phasespace
typedef NumMatrix<double, DIMENSION+3> VlasovDist;

//seems not to be used anymore ==> commenting out
//extern bool VlasovDiagnostic;
//-----------------------------------------------------------------------------
#endif //VLASOV_H
//-----------------------------------------------------------------------------
//=============================================================================
//Main Page documentation & additional pages 
//=============================================================================
//mainpage
/** @mainpage 2d3v Vlasov Darwin Code
 *
 *  @section int Introduction
 *
 * This code solves the Vlasov equation directly in 2 spatial and 3 velocity dimensions for collisionless plasmas in 
 * the nonrelativistic limit.
 * Advantages are the resolution of the high energy tails of the distribution function and low to none numerical 
 * noise, in contrast to a higher computational effort compared to PIC codes.
 *
 * In order to suppress purely electromagnetic modes in Maxwells equations, the Darwin approximation is implemented, 
 * elimating the necessity to 
 * fullfill the CFL condition \f$ \Delta t \leq \frac{\Delta x}{c} \f$.
 *
 * Next to some other schemes, a third order flux conservative scheme is implemented for spatial integration. This 
 * guarantees the conservation of
 * extremal points in the distribution function, which is a characteristic of the Vlasov equation.
 *
 * Various types of fields can be used for computation, from electrostatic to full electromagnetic fields. 
 * Both Darwin and Boris scheme are available as approximation for the electromagnetic force.
 *
 * Regarding the large amounts of data to be processed, the code is implemented to run on multiple processors systems.
 * To communicate results between processes the mpi-protocols are used. 
 * For simple problems it is also possible to run the programm on a single processor system, though this will limit the 
 * choice of global boundary conditions.
 *
 * Each process process computes an individual result on a subdomain of the problem. The exchange of these results is
 * implemented as an exchange of boundary conditions to neighbouring processes.
 * On each node basically runs an identical copy of the programm with individual boundary conditions and distribution 
 * functions.
 * 
 * A single process is promoted the master process, which will perform global diagnostics.
 *  
 * @section install Installing & Running the code
 * @subsection ins Installation
 * 
 * To install, use "make" in the source directory, which will build the appropriate binary files. 
 * The architectures are specified in the "makeinclude" file, which should be found in the directory containing the
 * source.
 *
 * "make" alone will build only the "vlasov"-executable of the available versions (like vlasov, vlasov_esb, ...) as a default.
 * Several options can be specified to build specific versions of the code:
 *  - vlasov (default value)
 *  - vlasov_esb   equals electrostatic field with constant magnetic field
 *  - vlasov_ms    equals magnetostatic field
 *  - vlasov_const constant fields
 *
 *
 * The HDF-libraries are mandatory to install (currently used version: HDF 4.1). In future HDF may be replaced by a  
 * different data exchange protocol.
 * This file should be located in the "lib" directory and is called "libfhdf.a". 
 * 
 * Furthermore, for inter process communication, the MPI header is used and therefor needed at compiling the code.
 * This header should be found in the "lib" directory as well ("mpi.h").
 * For information on MPI see the link in the question and answers section below (\ref as).
 *
 * @subsection run Running the code 
 * 
 * To run this code build the binary and provide a setup file in the working directory (see here: \ref I).
 * If executing the code on a single processor machine simply call the programm via "./vlasov" where "vlasov" stands
 * for any of the programm versions (assuming linux as an operating system).
 * 
 * On a multiprocessor machine the mpi runtime envrioment is needed to execute the binaries. An instance of this is 
 * invoked by the "lamboot" command. Prior to this the execution command will not work.
 * 
 * With the mpi enviroment running, call "mpirun -np X vlasov" where "X" stands for the number of processors the 
 * program is supposed to run on, while "vlasov" again stands for the desired version of the code. 
 * (The option "-n" sees to running it now, while "-p" specifies the number of processors to run on via the number "X")
 * 
 * 
 * @section as Advancers & schemes
 *
 * Schemes take care of spatial integration, while Advancers carry out the timestep. Note that the spatial integration
 * scheme is a parameter which is pluged into the Advancer, which performs the timesplitting and calls the 
 * spatial integrator at the appropiate time.
 *
 * In order to implement different schemes of integration and advancers, the file "simtype.h" has to be edited
 * manually. 
 *
 * Selection is performed by changing the preprocessor variable in "simtype.h" to the appropiate value, 
 * a default is: POS_FLUX_3_BACKSUBS. Possible values are:
 * - POS_FLUX_3_BACKSUBS
 * - POS_FLUX_3_BAKCSUBS_RK
 * - POS_FLUX_3_TESTING
 * - POS_FLUX_3_TIMESPLIT
 * - POS_FLUX_3_TIMESPLIT_FIRST_ORDER
 * - FINITE_DIFF
 * - FINITE_DIFF_RK3
 * - FINITE_DIFF4_RK3
 *
 * There are several advancers and integration schemes to choose  from:
 * - Advancer can be one of the following:
 *    -# Eulerian: SimpleAdvance (works only with FiniteDiffScheme)
 *    -# Leap Frog: LeapFrogAdvance (deprecated)
 *    -# Simplified Leap Frog: SimpleLeapFrogAdvance
 *    -# Runge Kutta: RungeKuttaAdvance
 *
 * - Scheme can be one of the following:
 *    -# Finite Differences: FiniteDiffScheme (works only with  SimpleAdvance)
 *    -# 3rd order positive Flux Conservative: PosFluxCons3rdOrder
 *    -# Piecewise Parabolic: PParabolicScheme (currently not operational)
 *
 * Note further that the finite difference advancers have a build in Runge Kutta integrators of variing order, since 
 * they do not work with the Scheme design. 
 * 
 * The piecewise parabolic scheme is currently not operational, but will be implemented.
 *
 * @section IO Input & Output
 * 
 * @subsection I Input Data & Setup
 * The data for setup is stored in the file "setup.dat", which should be locate in the directory the executable is in.
 * This file is read by the Rebuildable class and its derivatives.
 * 
 * The syntax understood by the parser included in Rebuildable allows comments in C-style (//...), which are ignored.
 * Each parameter is supplied in plain text with the according value separate by a whitespace:
 * \verbatim
  dt  ...
  \endverbatim
 * and by a newline from the next. 
 * 
 * Blocks of related parameters are listed as shown: 
 * \verbatim
       blockname { 
                    parameter1  ...
                    parameter2  ...
                    ...         ...
                }
  \endverbatim
 * Nesting of parameter blocks is possible.
 * 
 * For an example see the \ref faq section
 * 
 * @subsection O Output & Result
 * 
 * The data returned from the simulation is stored in ASCII-files specified in "setup.dat" by the Diagnostic objects.
 * In "setup.dat" the entry :
 * \verbatim 
   diagnosticsname {
                    file       filename.out
                    interval   x
                    append     y
                   }
  \endverbatim
 * sets the mode and destination of output. 
 * Interval is the desired number of timesteps between writings.
 * If "append y" is specified, each output is appended to the file separate by a blank line.
 * Note, that "y" is the default for append and can  be omitted.
 * If "append n" the file will be overwritten, unless the filename is supplied as follows (or equivalent):
  \verbatim 
   diagnosticsname {
                    file       filename#p-#t.out
                    ...
                   }
  \endverbatim
 * In this setup, the token "#p" is replaced with the process number and "#t" with the timestep.
 * 
 * @section multi Multi-Processor Runs
 * To reduce computation time the spatial grid  is split in equal portions, which can be 
 * squares or rectangles of full height.
 *
 * On each subgrid a separate process advances the distribution function in time and (phase-) space. Thus the amount of 
 * data per process is reduced to manageable size.
 *
 * In order to conserve the consistency of the simulation, each process sis needed to communicate its results
 * to the whole network. This is done through ghostcells adjacent to the processing region, which function as boundary 
 * or starting values to the process itself and are send via MPI routines to the neighbouring processes.
 * There these supply the new boundary conditions for this process. But these are not boundary conditions in the sense, 
 * that they remain unchanged by the process, computation works as normal for the ghostcells.
 * 
 * Exchange is done one direction in phasespace after another.
 *
 * For Fields and Densities one set of ghostcells is needed in each direction, while the distribution function has two 
 * sets, due to higher order interpolation.
 *
 * @section furth Further Information
 *
 * See 
 * 
 * \b \em  "Darwin-Vlasov simulations of magnetised plasmas" 
 * 
 *(H.Schmitz, R.Grauer, RUB) 
 *
 * for further details.
 *
 *    http://www.tp1.ruhr-uni-bochum.de/~hs/papers/JCompPhys2006.pdf
 *
 * Published in: Journal of Computational Physics, 14 october 2005
 * 
 * @section faq Questions and Answers
 *
 *  - For HDF library documentation see: http://hdf.ncsa.uiuc.edu
 *
 *  - The  HDF4.1 librarys can be acquired via the ncsa: http://hdf.ncsa.uiuc.edu/hdf4.html
 * 
 *  - For MPI information see: http://www.mpi-forum.org/index.html
 *
 *  - A default makefile is found here: @ref make
 *
 *  - An example for "setup.dat" is found here: \ref setup
 *
 *  - All class template member functions are defined in the files ending *.t or xxx_temp.cpp, to separate these from 
 *    normal implementation files
 * 
 *  - For the use of indices in the arrays see the page @ref indices
 *
 *  - The makeinclude files should include all architecture this code is supposed to run on: 
 *		- HP-UX-9000
 *		- SunOS/SunOS_gnu
 *		- Linux-386_64 (tested)
 *		- Linux-alpha/Linux-alpha.g++
 *		- Linux-i386 (tested)
 *		- IRIX64 
 *		- Linux-alpha.cxx
 *		- Linux-ia64
 */
//mainpage
//-----------------------------------------------------------------------------
//Makefile 

/** @page make A Default Makefile
  * This makefile is supplied as a default and will build all versions of the code if called with "make":
\verbatim
include $(HOME)/src/makeinclude/include.mak

#OFLAGS  = $(X_OFLAGS) -DNDEBUG
OFLAGS  = $(X_OFLAGS) 

BASE_CXXFLAGS = $(OFLAGS) -I. $(X_INC) $(X_NO_AUTO_INCLUDE)

#the implementation files
VLASOV_SOURCES = darwin.cpp helmholtz.cpp numeric.cpp potential.cpp \
		vlasovinit.cpp wrapvlasov.cpp diaghelper.cpp main.cpp \
		poisson.cpp vlasovinit.cpp  vlasovbase.cpp process.cpp \
		derivedfields.cpp forcefield.cpp hdfstream.cpp task.cpp \
		globals.cpp diagnostic.cpp fielddiag.cpp reconnection.cpp \
		scalarfield.cpp magnetostatic.cpp bernstein.cpp

#the headers
VLASOV_HEADERS = darwin.h diaghelper.h index.h numeric.h opar.h \
		parameter.h poisson.h potential.h stlpwrapper.h \
		vlasov.h vlasovinit.h vlasovbase.h wrapvlasov.h scheme.h \
		derivedfields.h forcefield.h hdfstream.h task.h parameter.h \
		process.h globals.h diagnostic.h fielddiag.h reconnection.h \
		numboundary.h magnetostatic.h bernstein.h

#object files, for linking
OBJECTS = boundary.o darwin.o diaghelper.o helmholtz.o numeric.o \
		potential.o poisson.o wrapvlasov.o main.o \
		vlasovinit.o vlasovbase.o forcefield.o process.o\
		derivedfields.o hdfstream.o task.o globals.o diagnostic.o \
		fielddiag.o reconnection.o scalarfield.o magnetostatic.o \
		bernstein.o

HEADERS = $(VLASOV_HEADERS)

LDFLAGS = $(X_LIB)

LOADLIBS = $(X_HDF_LIBS) -lm

BINARY_PREFIX = $(X_ROOT_DIR)/bin/$(X_OSTYPE)

all: vlasov

vlasov:	CXXFLAGS = $(BASE_CXXFLAGS) -DSIM_DARWIN
vlasov: $(OBJECTS) 
	$(CXX)  -g $^ -o $@ $(LDFLAGS) $(LOADLIBS)
	cp $@ $(HOME)/bin/$(X_OSTYPE)/

vlasov_es:	CXXFLAGS = $(BASE_CXXFLAGS) -DSIM_ELECTROSTATIC
vlasov_es: $(OBJECTS) 
	$(CXX)  -g $^ -o $@ $(LDFLAGS) $(LOADLIBS)
	cp $@ $(HOME)/bin/$(X_OSTYPE)/

vlasov_esb:	CXXFLAGS = $(BASE_CXXFLAGS) -DSIM_ELECTROSTATIC_B
vlasov_esb: $(OBJECTS) 
	$(CXX)  -g $^ -o $@ $(LDFLAGS) $(LOADLIBS)
	cp $@ $(HOME)/bin/$(X_OSTYPE)/

vlasov_ms:	CXXFLAGS = $(BASE_CXXFLAGS) -DSIM_MAGNETOSTATIC
vlasov_ms: $(OBJECTS) 
	$(CXX)  -g $^ -o $@ $(LDFLAGS) $(LOADLIBS)
	cp $@ $(HOME)/bin/$(X_OSTYPE)/

vlasov_const:	CXXFLAGS = $(BASE_CXXFLAGS) -DSIM_CONSTFIELD
vlasov_const: $(OBJECTS) 
	$(CXX)  -g $^ -o $@ $(LDFLAGS) $(LOADLIBS)
	cp $@ $(HOME)/bin/$(X_OSTYPE)/

test_helmholtz:  helmholtz.o Vlas2d.o test_helmholtz.o
	$(CXX)  $^ -o $@ $(LDFLAGS) $(LOADLIBS)

 *%.o: %.cpp
 	$(CXX) -c $(CXXFLAGS) $(INCLUDE) $<
 
 clean:
 	-rm -f *.o core vlasov vlasov_es vlasov_const
\endverbatim
*/
//mainpage
//-----------------------------------------------------------------------------
//Array Indices

/** @page indices Array indices
 *
 *  @section general_ind General
 *
 *  We divide the whole simulation area into a number of grid cells. The
 *  interior of the simulation is indexed from \f$1\ldots N\f$. If we have 
 *  periodic boundary conditions, this means that \f$A(0) = A(N)\f$, of in 
 *  general \f$A(i) = A(N+i)\f$.
 *
 *  The physical extent of the simulation box is \f$N\Delta x\f$.
 *
 *  @section disttrib_ind Distribution Function
 *
 *  The distribution function may be distributed among a number
 *  of processors. For the third order advancer we need at least two
 *  cells to every side of the. The whole simulation region in memory thus 
 *  has the extent \f$ -1 \ldots N+2 \f$.
 *  
 *  @section fields_ind Fields and Potentials solver
 *
 *  Internally the Poisson solver uses the Fortran hwscrt routine. 
 *  In the case of periodic boundaries this routine needs one cell to be 
 *  duplicated. This means the region spans \f$1\ldots N+1\f$, where
 *  \f$A(1) = A(N+1)\f$. If the boundaries are fixed, however we still want the 
 *  interior to span \f$1\ldots N\f$, so we need an array \f$0\ldots N+1\f$.
 *  In general potentials will thus span \f$0\ldots N+1\f$ and in the case of
 *  periodic boundary conditions only the sub-array \f$1\ldots N+1\f$ will
 *  be passed to hwscrt. 
 *
 *  The boundary points can be exploited when spatial derivatives of the 
 *  potentials have to be taken.
 *
 */
//Array Indices
//-----------------------------------------------------------------------------
//setup.dat
/** @page setup Example setup.dat
\verbatim
//this will ignored by the parser
// \lambda_D = 0.1*sqrt(2)
// \beta = 1
// \alpha = 0.2

// \lambda = 0.5
// \gamma = T_e / T_i = 0.2

// d_i := 1
// d_e = 0.2
// m_i/m_e = 25

//general parameters
general {
    //physical data of the simulation grid
    x-size  128
    y-size  2
    v-size-x 20
    v-size-y 20
    v-size-z 20
    
    //extension of a single cell
    x-space 0.05
    y-space 0.05
    
    //timestep width
    dt  0.01       
    
    //velocity over speed of light
    vel-ratio   0.1

    // Total run time
    T   1000       
    
    //total density
    bg-density  1
    
    //external magnetic field
    Bx  1
    
    //diagnostics for the magnetic field, each component is appended
    //to a separate file
    field-diag {
        bx {
            file        bx.out
            interval    1
            append      y
            }
        by {
            file        by.out
            interval    1
            append      y
            }
        bz {
            file        bz.out
            interval    1
            append      y
            }
        ex {
            file        ex.out
            interval    1
            append      y
            }
        ey {
            file        ey.out
            interval    1
            append      y
            }
        ez {
            file        ez.out
            interval    1
            append      y
            }
    }
}
//end of general

//type of boundary 
boundary {
    single-periodic
}


//parameters for a species
// ===============
//    Electrons
// ===============
species {
    charge -1
    mass    1
    vx-range    4
    vy-range    4
    vz-range    4

    //initialize with maxwell distribution
    init {
        maxwell {            
            Therm_vx    1.0
            Therm_vy    1.0
            Therm_vz    1.0
            N0          1.0
            N1          0.05
            kx          1
            ky          0
        }
    }

    //diagnostics, each output is written to a new file
    phase-space {
        file        electron#t-#p.out
        interval    100
                }

    fielddiag {
        class       rho
        field       rho
        file        el_rho.out
        interval    1
        append      y
    }

    fielddiag {
        class       vels
        field       Jx
        file        el_vx.out
        interval    1
        append      y
    }

    fielddiag {
        class       vels
        field       Jy
        file        el_vy.out
        interval    1
        append      y
    }

    fielddiag {
        class       vels
        field       Jz
        file        el_vz.out
        interval    1
        append      y
    }
}
\endverbatim
*/
//setup.dat
//-----------------------------------------------------------------------------

