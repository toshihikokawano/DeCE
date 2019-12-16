LDFLAGS		= -lm
CPPFLAGS	= -g -Wall -O3 -Wextra
CPP		= g++
RM      	= rm
LIBOBJS		= endflib.o endfio.o
PRGOBJS		= example1.o example2.o example3.o
PROG		= example1 example2 example3

all:		$(PROG)

example1:	$(LIBOBJS) example1.o
		$(CPP) example1.o $(LIBOBJS) $(LDFLAGS) -o example1

example2:	$(LIBOBJS) example2.o
		$(CPP) example2.o $(LIBOBJS) $(LDFLAGS) -o example2

example3:	$(LIBOBJS) example3.o
		$(CPP) example3.o $(LIBOBJS) $(LDFLAGS) -o example3

clean:
		$(RM) -f $(LIBOBJS) $(PROG) $(PRGOBJS)

example1.o:	../source/endflib.h
example2.o:	../source/endflib.h
example3.o:	../source/endflib.h
endflib.o:	../source/endflib.cpp ../source/endflib.h
		$(CPP) $(CPPFLAGS) ../source/endflib.cpp -c
endfio.o:	../source/endfio.cpp ../source/endflib.h
		$(CPP) $(CPPFLAGS) ../source/endfio.cpp -c
