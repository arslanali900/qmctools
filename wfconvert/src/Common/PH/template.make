SOURCES = PH.cc Chebyshev.cc Potential.cc Atom.cc Cost.cc CoreTransform.cc    \
	  ScreenedPot.cc QuinticPH.cc CoulombPot.cc SplinePot.cc HeAzizPot.cc \
	  kSpacePH.cc TestkSpacePH.cc

Objs =	PH.o Chebyshev.o Potential.o Atom.o Cost.o CoreTransform.o            \
	ScreenedPot.o QuinticPH.o CoulombPot.o SplinePot.o HeAzizPot.o        \
        kSpacePH.o 

TestObjs  = Potential.o ScreenedPot.o QuinticPH.o CoulombPot.o SplinePot.o    \
            HeAzizPot.o kSpacePH.o TestkSpacePH.o

IOObjs = ../IO/InputOutput.o ../IO/InputOutputASCII.o ../IO/InputOutputHDF5.o \
 ../IO/InputOutputXML.o

SplineObjs = ../Splines/BicubicSpline.o ../Splines/CubicSpline.o              \
             ../Splines/Grid.o ../Splines/MyTricubicSpline.o                  \
             ../Splines/QuinticSpline.o ../Splines/QuinticSplines.o 

MiscObjs =   ../Integration/GKIntegration.o ../Fitting/Fitting.o              \
             ../MatrixOps/MatrixOps.o

all:    TestkSpacePH $(Objs)

TestkSpacePH:  $(TestObjs) $(IOObjs) $(SplineObjs) $(MiscObjs)
	$(LD) -o $@ $(TestObjs) $(IOObjs) $(SplineObjs) $(MiscObjs) $(LIBS)

clean:
	rm -f *.o

.cc.o: $(HEADERS)
	$(CC) $(CCFLAGS) $(DEFS) $(INCL) -o $*.o $< 
.f.o:
	g77 -c -o $*.o $<

newmake:
	+$(MAKE) -f template.make Makefile FRC=force_rebuild

Makefile:	$(FRC)
	rm -f $@
	cp template.make $@
	echo 'Automatically generated dependency list:' >> $@
	$(MAKECC) $(CCFLAGS) $(INCL) -MM $(SOURCES) >> $@
	chmod -w $@


force_rebuild:
