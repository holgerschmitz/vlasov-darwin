include $(HOME)/src/makeinclude/include.mak

OFLAGS  = $(X_OFLAGS) -DNDEBUG

CXXFLAGS = $(OFLAGS) -I. $(X_INC) $(X_NO_AUTO_INCLUDE)

VLASOV_SOURCES = darwin.cpp helmholtz.cpp numeric.cpp potential.cpp \
		vlasovinit.cpp wrapvlasov.cpp diaghelper.cpp main.cpp \
		poisson.cpp vlasovinit.cpp  vlasovbase.cpp process.cpp \
		derivedfields.cpp forcefield.cpp hdfstream.cpp task.cpp \
		globals.cpp

VLASOV_HEADERS = darwin.h diaghelper.h index.h numeric.h opar.h \
		parameter.h poisson.h potential.h stlpwrapper.h \
		vlasov.h vlasovinit.h vlasovbase.h wrapvlasov.h scheme.h \
		derivedfields.h forcefield.h hdfstream.h task.h parameter.h \
		process.h globals.h

OBJECTS = boundary.o darwin.o diaghelper.o helmholtz.o numeric.o \
		potential.o poisson.o wrapvlasov.o main.o \
		vlasovinit.o vlasovbase.o forcefield.o process.o\
		derivedfields.o hdfstream.o task.o globals.o

HEADERS = $(VLASOV_HEADERS)

LDFLAGS = $(X_LIB)

HDFLIBS = -lmfhdf -ldf -ljpeg -lz
LOADLIBS = $(HDFLIBS) -lutil_mt -lmatrix_mt -lcyclic_mt -lm

BINARY_PREFIX = $(X_ROOT_DIR)/bin/$(X_OSTYPE)

all: vlasov

vlasov:  $(OBJECTS)
	$(CXX)  $^ -o $@ $(LDFLAGS) $(LOADLIBS)

test_helmholtz:  helmholtz.o Vlas2d.o test_helmholtz.o
	$(CXX)  $^ -o $@ $(LDFLAGS) $(LOADLIBS)

.s.o:
	$(COMPILE) -c $<

.S.o:
	$(COMPILE) -c $<

clean:
	-rm -f *.o core vlasov


