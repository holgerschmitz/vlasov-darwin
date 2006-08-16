include $(HOME)/src/makeinclude/include.mak

#OFLAGS  = $(X_OFLAGS) -DNDEBUG
OFLAGS  = $(X_OFLAGS) 

BASE_CXXFLAGS = $(OFLAGS) -I. $(X_INC) $(X_NO_AUTO_INCLUDE)

VLASOV_SOURCES = darwin.cpp helmholtz.cpp numeric.cpp potential.cpp \
		vlasovinit.cpp wrapvlasov.cpp diaghelper.cpp main.cpp \
		poisson.cpp vlasovinit.cpp  vlasovbase.cpp process.cpp \
		derivedfields.cpp forcefield.cpp hdfstream.cpp task.cpp \
		globals.cpp diagnostic.cpp fielddiag.cpp reconnection.cpp \
		scalarfield.cpp magnetostatic.cpp bernstein.cpp shock.cpp \
		openbound.cpp

VLASOV_HEADERS = darwin.h diaghelper.h index.h numeric.h opar.h \
		parameter.h poisson.h potential.h stlpwrapper.h \
		vlasov.h vlasovinit.h vlasovbase.h wrapvlasov.h scheme.h \
		derivedfields.h forcefield.h hdfstream.h task.h parameter.h \
		process.h globals.h diagnostic.h fielddiag.h reconnection.h \
		numboundary.h magnetostatic.h bernstein.h shock.h \
		openbound.h

OBJECTS = boundary.o darwin.o diaghelper.o helmholtz.o numeric.o \
		potential.o poisson.o wrapvlasov.o main.o \
		vlasovinit.o vlasovbase.o forcefield.o process.o\
		derivedfields.o hdfstream.o task.o globals.o diagnostic.o \
		fielddiag.o reconnection.o scalarfield.o magnetostatic.o \
		bernstein.o shock.o openbound.o

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

vlasov_esb_bs_ex:	CXXFLAGS = $(BASE_CXXFLAGS) -DSIM_ELECTROSTATIC_B
vlasov_esb_bs_ex: $(OBJECTS) 
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

%.o: %.cpp
	$(CXX) -c $(CXXFLAGS) $(INCLUDE) $<


clean:
	-rm -f *.o core vlasov vlasov_es vlasov_const


