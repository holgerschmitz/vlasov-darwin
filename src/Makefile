include $(HOME)/src/makeinclude/include.mak

OFLAGS  = $(X_OFLAGS) -DNDEBUG

CXXFLAGS = $(OFLAGS) -I. $(X_INC) $(X_NO_AUTO_INCLUDE)

VLASOV_SOURCES = darwin.cpp helmholtz.cpp numeric.cpp potential.cpp \
		vlasovinit.cpp wrapvlasov.cpp diaghelper.cpp main.cpp \
		poisson.cpp Vlas2d.cpp 
VLASOV_HEADERS = darwin.h diaghelper.h index.h numeric.h opar.h \
		parameter.h poisson.h potential.h stlpwrapper.h Vlas2d.h \
		vlasov.h vlasovinit.h wrapvlasov.h

OBJECTS = darwin.o helmholtz.o numeric.o potential.o \
		vlasovinit.o wrapvlasov.o diaghelper.o main.o \
		poisson.o Vlas2d.o

HEADERS = $(OPAR_HEADERS)

LDFLAGS = $(X_LIB)
LOADLIBES = -lutil_mt -lmatrix_mt -lcyclic_mt -lm

BINARY_PREFIX = $(X_ROOT_DIR)/bin/$(X_OSTYPE)

all: vlasov

vlasov:  $(OBJECTS)
	@rm -f vlasov
	$(CXX) $^ $(LDFLAGS) $(LOADLIBES) -o $@
	cp vlasov $(BINARY_PREFIX)/vlasov


.s.o:
	$(COMPILE) -c $<

.S.o:
	$(COMPILE) -c $<

clean:
	-rm -f *.o core vlasov


