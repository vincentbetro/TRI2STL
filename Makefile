
SRCDIRS = GEOMETRY IO UTIL MAIN

all:
	@for dir in $(SRCDIRS); do \
	  cd $$dir ; \
	  pwd ; \
	  $(MAKE) ; \
	  cd .. ; \
	done

clean:
	@for dir in $(SRCDIRS); do \
	  cd $$dir ; \
	  pwd ; \
	  $(MAKE) clean ; \
	  cd .. ; \
	done
	/bin/rm -f TRI2STL.MACOSX

cleanest:
	@for dir in $(SRCDIRS); do \
	  cd $$dir ; \
	  pwd ; \
	  $(MAKE) clean ; \
	  cd .. ; \
	done
	/bin/rm -f TRI2STL.MACOSX
	/bin/rm -f *.jou
	/bin/rm -f *.params
	/bin/rm -f debug.*
	/bin/rm -f *.uns
