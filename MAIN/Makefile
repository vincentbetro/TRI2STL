
include ../Makefile.options

LIBRARY_PATH += -L../GEOMETRY -L../IO -L../UTIL -L../TETGEN

LIBS += -lgeometry -lio -lutility -loctree -lSpaceField -lcgns -lhdf5 -ltet
#LIBS += -lgeometry -lio -lutility -loctree -lcgns /usr/local/hdf5-1.6.4-icc/lib/libhdf5.a -ltet

LDFLAGS = $(LIBRARY_PATH) $(LIBS)

OPTC_SRCS = TRI2STL_WRAP.cpp \
			TRI2STL.cpp \

OPTSRCS = $(OPTC_SRCS)
OPTOBJECTS = $(OPTC_SRCS:.cpp=.o)

all: TRI2STL

TRI2STL: geometry io utility $(OPTOBJECTS)
	$(LINK) -o $@ $(OPTOBJECTS) $(LDFLAGS)
	cp $@ ../$@.$(EXE_SUFFIX)

geometry:
	cd ../GEOMETRY; $(MAKE)

io:
	cd ../IO; $(MAKE)

utility:
	cd ../UTIL; $(MAKE)

clean:
	/bin/rm -f *.o
	/bin/rm -f TRI2STL

include ../suffix.rules
