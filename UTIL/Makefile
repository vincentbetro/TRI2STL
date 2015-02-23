
include ../Makefile.options

SRCS=\
	FV_routines.cpp \
	List.cpp \
	Util.cpp \
	sort.cpp \
	svdcmp.cpp \

OBJECTS = $(SRCS:.cpp=.o)

utility.a:		$(OBJECTS)
		$(LINKSTATIC) lib$@ $(OBJECTS)
		$(RANLIB) lib$@

utility.so:	$(OBJECTS)
		$(LINKSHARED) $@ $(OBJECTS)

clean:
	/bin/rm -f *.o
	/bin/rm -f libutility.a
	/bin/rm -f libutility.so

include ../suffix.rules

