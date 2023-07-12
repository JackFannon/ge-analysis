FC	= gfortran
CFLAGS	= -O -g
LD	= gfortran

#
CERN = /usr/local/sklib_gcc8/cern/2005
#CERN = /usr/local/sklib_g77/cern/2005
LIBDIR	= -L$(CERN)/lib
CLIB	= -lgeant321 -ljetset74 -lmathlib \
          -lgraflib -lgrafX11 -lpacklib -lkernlib
GLIB	= -lX11 -lnsl
LIBS	= $(CLIB) $(GLIB)
FFLAGS = -O0 -I$(CERN)/include/geant321

PROGRAM   = gedeo

OBJS      = gedeo.o\
            uginit.o\
	    ugeom120830.o\
	    gukine0117.o\
	    gustep120810.o\
	    guhadr.o\
	    guout120810.o\
	    guphad.o\
	    gutrev120810.o\
	    uffkey.o\
	    uhinit21.o\
	    viewyz.o\
            gsum.o\
            gsum2.o\
	    uglast.o\
            rngaus.o\

.f.o	: ugeant.h gso.h
	  $(FC) $(FFLAGS) -c $<

all	: $(PROGRAM)

$(PROGRAM)	: $(OBJS)
	$(LD) $(CFLAGS) -o $(PROGRAM) $(OBJS) $(LIBDIR) $(LIBS)

clean	:
	rm -f *.o $(PROGRAM) *~ */*~



