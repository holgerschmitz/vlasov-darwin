include $(HOME)/src/makeinclude/include.mak

OFLAGS  = $(X_OFLAGS) -DNDEBUG

CXXFLAGS = $(OFLAGS) -I. $(X_INC) $(X_NO_AUTO_INCLUDE)

VLASOV_SOURCES = darwin.cpp helmholtz.cpp numeric.cpp potential.cpp \
		vlasovinit.cpp wrapvlasov.cpp diaghelper.cpp main.cpp \
		poisson.cpp Vlas2d.cpp 

VLASOV_HEADERS = darwin.h diaghelper.h index.h numeric.h opar.h \
		parameter.h poisson.h potential.h stlpwrapper.h Vlas2d.h \
		vlasov.h vlasovinit.h wrapvlasov.h scheme.h

OBJECTS = boundary.o darwin.o diaghelper.o helmholtz.o numeric.o \
		potential.o poisson.o wrapvlasov.o scheme.o main.o \
		Vlas2d.o

HEADERS = $(VLASOV_HEADERS)

LDFLAGS = $(X_LIB)
LOADLIBES = -lutil_mt -lmatrix_mt -lcyclic_mt -lm

BINARY_PREFIX = $(X_ROOT_DIR)/bin/$(X_OSTYPE)

CXX = g++
#CXX = /usr/local/mpich/bin/mpiCC
#CXX = /usr/lib/mpich/bin/mpiCC

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


