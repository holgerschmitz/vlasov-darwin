include $(HOME)/src/makeinclude/include.mak

OFLAGS  = $(X_OFLAGS) -DNDEBUG

BASE_CXXFLAGS = $(OFLAGS) -I. $(X_INC) $(X_NO_AUTO_INCLUDE)

VLASOV_SOURCES = darwin.cpp helmholtz.cpp numeric.cpp potential.cpp \
		vlasovinit.cpp wrapvlasov.cpp diaghelper.cpp main.cpp \
		poisson.cpp vlasovinit.cpp  vlasovbase.cpp process.cpp \
		derivedfields.cpp forcefield.cpp hdfstream.cpp task.cpp \
		globals.cpp diagnostic.cpp fielddiag.cpp reconnection.cpp \
		scalarfield.cpp

VLASOV_HEADERS = darwin.h diaghelper.h index.h numeric.h opar.h \
		parameter.h poisson.h potential.h stlpwrapper.h \
		vlasov.h vlasovinit.h vlasovbase.h wrapvlasov.h scheme.h \
		derivedfields.h forcefield.h hdfstream.h task.h parameter.h \
		process.h globals.h diagnostic.h fielddiag.h reconnection.h \
		numboundary.h

OBJECTS = boundary.o darwin.o diaghelper.o helmholtz.o numeric.o \
		potential.o poisson.o wrapvlasov.o main.o \
		vlasovinit.o vlasovbase.o forcefield.o process.o\
		derivedfields.o hdfstream.o task.o globals.o diagnostic.o \
		fielddiag.o reconnection.o scalarfield.o

HEADERS = $(VLASOV_HEADERS)

LDFLAGS = $(X_LIB)

LOADLIBS = $(X_HDF_LIBS) -lm

BINARY_PREFIX = $(X_ROOT_DIR)/bin/$(X_OSTYPE)

all: vlasov

vlasov:	CXXFLAGS = $(BASE_CXXFLAGS) -DSIM_DARWIN
vlasov: $(OBJECTS) 
	$(CXX)  $^ -o $@ $(LDFLAGS) $(LOADLIBS)
	cp $@ ~/bin/$(X_OSTYPE)

vlasov_es:	CXXFLAGS = $(BASE_CXXFLAGS) -DSIM_ELECTROSTATIC
vlasov_es: $(OBJECTS) 
	$(CXX)  $^ -o $@ $(LDFLAGS) $(LOADLIBS)
	cp $@ ~/bin/$(X_OSTYPE)

vlasov_const:	CXXFLAGS = $(BASE_CXXFLAGS) -DSIM_CONSTFIELD
vlasov_const: $(OBJECTS) 
	$(CXX)  $^ -o $@ $(LDFLAGS) $(LOADLIBS)
	cp $@ ~/bin/$(X_OSTYPE)

test_helmholtz:  helmholtz.o Vlas2d.o test_helmholtz.o
	$(CXX)  $^ -o $@ $(LDFLAGS) $(LOADLIBS)

.s.o:
	$(COMPILE) $(INCLUDE) -c $<

.S.o:
	$(COMPILE) $(INCLUDE) -c $<

clean:
	-rm -f *.o core vlasov vlasov_es vlasov_const


