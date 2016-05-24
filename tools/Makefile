LDFLAGS		= -lm
CPPFLAGS	= -g -Wall -O4 -Wextra
CPP		= g++
RM      	= rm
LIBOBJS		= endflib.o endfio.o
PRGOBJS		= decemf4.o decemf5.o decemf6.o decemf12.o decemacs.o
PROG		= decemf4 decemf5 decemf6 decemf12 decemacs

all:		$(PROG)

decemf4:	$(LIBOBJS) decemf4.o
		$(CPP) decemf4.o $(LIBOBJS) $(LDFLAGS) -o decemf4

decemf5:	$(LIBOBJS) decemf5.o
		$(CPP) decemf5.o $(LIBOBJS) $(LDFLAGS) -o decemf5

decemf6:	$(LIBOBJS) decemf6.o
		$(CPP) decemf6.o $(LIBOBJS) $(LDFLAGS) -o decemf6

decemf12:	$(LIBOBJS) decemf12.o
		$(CPP) decemf12.o $(LIBOBJS) $(LDFLAGS) -o decemf12

decemacs:	$(LIBOBJS) decemacs.o
		$(CPP) decemacs.o $(LIBOBJS) $(LDFLAGS) -o decemacs

clean:
		$(RM) -f $(LIBOBJS) $(PRGOBJS) $(PROG)

decemf4.o:	../source/endflib.h
decemf5.o:	../source/endflib.h
decemf6.o:	../source/endflib.h
decemf12.o:	../source/endflib.h
decemacs.o:	../source/endflib.h
endflib.o:	../source/endflib.cpp ../source/endflib.h
		$(CPP) $(CPPFLAGS) ../source/endflib.cpp -c
endfio.o:	../source/endfio.cpp ../source/endflib.h
		$(CPP) $(CPPFLAGS) ../source/endfio.cpp -c
