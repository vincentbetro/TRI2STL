
include ../Makefile.options

SRCS=\
	geometry.cpp \

OBJECTS = $(SRCS:.cpp=.o)

geometry.a:	$(OBJECTS)
		$(LINKSTATIC) lib$@ $(OBJECTS)
		$(RANLIB) lib$@

geometry.so:	$(OBJECTS)
		$(LINKSHARED) $@ $(OBJECTS)

clean:
	/bin/rm -f *.o
	/bin/rm -f libgeometry.a
	/bin/rm -f libgeometry.so

include ../suffix.rules
