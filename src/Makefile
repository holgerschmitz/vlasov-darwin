include $(HOME)/src/makeinclude/include.mak

OFLAGS  = $(X_OFLAGS) -DNDEBUG

CXXFLAGS = $(OFLAGS) -I. $(X_INC) $(X_NO_AUTO_INCLUDE)

VLASOV_SOURCES = darwin.cpp helmholtz.cpp numeric.cpp potential.cpp \
		vlasovinit.cpp wrapvlasov.cpp diaghelper.cpp main.cpp \
		poisson.cpp Vlas2d.cpp vlasovinit.cpp  vlasovbase.cpp \
		derivedfields.cpp forcefield.cpp

VLASOV_HEADERS = darwin.h diaghelper.h index.h numeric.h opar.h \
		parameter.h poisson.h potential.h stlpwrapper.h Vlas2d.h \
		vlasov.h vlasovinit.h vlasovbase.h wrapvlasov.h scheme.h \
		derivedfields.h forcefield.h

OBJECTS = boundary.o darwin.o diaghelper.o helmholtz.o numeric.o \
		potential.o poisson.o wrapvlasov.o main.o \
		Vlas2d.o vlasovinit.o vlasovbase.o forcefield.o\
		derivedfields.o

HEADERS = $(VLASOV_HEADERS)

LDFLAGS = $(X_LIB)
LOADLIBES = -lutil_mt -lmatrix_mt -lcyclic_mt -lm

BINARY_PREFIX = $(X_ROOT_DIR)/bin/$(X_OSTYPE)

#CXX = /scratch/raid/hs/bin/i686-pc-linux-gnu-g++ -DSINGLE_PROCESSOR
#CXX = /usr/local/mpich/bin/mpiCC -D__USE_STD_IOSTREAM -use_ld_input=templateinst
#CXX = /usr/lib/mpich/bin/mpiCC
CXX = /usr/local/mpich/bin/mpiCC -w1 -fast -arch ev6 -O5 -std gnu -D__USE_STD_IOSTREAM  -noimplicit_include -accept restrict_keyword

all: vlasov

vlasov:  $(OBJECTS)
	$(CXX)  $^ -o $@ $(LDFLAGS) $(LOADLIBES)

test_helmholtz:  helmholtz.o Vlas2d.o test_helmholtz.o
	$(CXX)  $^ -o $@ $(LDFLAGS) $(LOADLIBES)

.s.o:
	$(COMPILE) -c $<

.S.o:
	$(COMPILE) -c $<

clean:
	-rm -f *.o core vlasov


