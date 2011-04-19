SOURCES = Atom.cc RadialWF.cc TestRadialWF.cc DFTAtom.cc TestDFTAtom.cc 

all:	Atom.o RadialWF.o TestRadialWF DFTAtom.o TestDFTAtom 

TestWFobjs =                     \
	TestRadialWF.o           \
	RadialWF.o               \
        ../PH/CoulombPot.o       \
	../Splines/CubicSpline.o \
	../IO/InputOutput.o      \
 	../IO/InputOutputHDF5.o  \
	../IO/InputOutputASCII.o \
	../IO/InputOutputXML.o

TestDFTAtomobjs =                   \
	RadialWF.o                  \
	TestDFTAtom.o               \
	DFTAtom.o                   \
	../PH/CoulombPot.o          \
	../PH/QuinticPH.o           \
	../PH/ScreenedPot.o         \
	../PH/SplinePot.o           \
	../PH/HeAzizPot.o           \
	../PH/Potential.o           \
	../Splines/CubicSpline.o    \
	../Splines/QuinticSpline.o  \
	../Splines/QuinticSplines.o \
	../DFT/Functionals.o        \
	../DFT/ExCorr.o             \
	../IO/InputOutput.o         \
	../IO/InputOutputHDF5.o     \
	../IO/InputOutputASCII.o    \
	../IO/InputOutputXML.o


TestRadialWF:  $(TestWFobjs)
	$(LD) -o TestRadialWF $(TestWFobjs) $(LIBS)

TestDFTAtom:  $(TestDFTAtomobjs)
	$(LD) -o TestDFTAtom $(TestDFTAtomobjs) $(LIBS)

clean:
	rm -f *.o

.cc.o: 
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
