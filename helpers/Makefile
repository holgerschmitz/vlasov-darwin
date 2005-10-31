include $(HOME)/src/makeinclude/include.mak

AR       = ar cqs

#OFLAGS  = $(X_OFLAGS) -DNDEBUG
OFLAGS  = $(X_OFLAGS) 

CXXFLAGS = $(OFLAGS) -I. $(X_INC) $(X_NO_AUTO_INCLUDE)

HELPERLIB_SOURCES = filename.cpp hdfstream.cpp multiproc.cpp

HELPERLIB_HEADERS = filename.h hdfstream.h multiproc.h

HELPERLIB_OBJECTS = filename.o hdfstream.o multiproc.o

PHASEBROP_OBJ = phaseprobe.o
REDIST_OBJ = redist.o

HEADERS = $(VLASOV_HEADERS)

LDFLAGS = $(X_LIB)

LOADLIBS = $(X_HDF_LIBS) -lm

BINARY_PREFIX = $(X_ROOT_DIR)/bin/$(X_OSTYPE)

all: libvlasov.a

libvlasov.a: $(HELPERLIB_OBJECTS) 
	$(AR) $@ $^
	cp $@ $(HOME)/lib/$(X_OSTYPE)/

phaseprobe: libvlasov.a $(PHASEBROP_OBJ)
	$(CXX)  -g $(PHASEBROP_OBJ) -o $@ $(LDFLAGS) $(LOADLIBS) -lvlasov
	cp $@ $(HOME)/bin/$(X_OSTYPE)/
	
redist: libvlasov.a $(REDIST_OBJ)
	$(CXX)  -g $(REDIST_OBJ) -o $@ $(LDFLAGS) $(LOADLIBS) -lvlasov
	cp $@ $(HOME)/bin/$(X_OSTYPE)/
	
%.o: %.cpp
	$(CXX) -c $(CXXFLAGS) $(INCLUDE) $<


clean:
	-rm -f *.o core libvlasov.a


