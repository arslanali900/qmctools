SOURCES = Hamiltonian.cc ConjGrad.cc TestPW.cc GVecs.cc                       \
          FFTBox.cc Hamiltonian2.cc TestPW2.cc PlaneWaves.cc ConjGrad2.cc     \
          BandStructure.cc


Objs = Hamiltonian.o ConjGrad.o GVecs.o FFTBox.o 

Objs2 = Hamiltonian2.o ConjGrad2.o GVecs.o FFTBox.o PlaneWaves.o

BSObjs = BandStructure.o

TestObjs  = TestPW.o

TestObjs2  = TestPW2.o

PHObjs = ../PH/Potential.o ../PH/ScreenedPot.o ../PH/QuinticPH.o              \
         ../PH/CoulombPot.o ../PH/SplinePot.o ../PH/HeAzizPot.o               \
         ../PH/kSpacePH.o

IOObjs = ../IO/InputOutput.o ../IO/InputOutputASCII.o ../IO/InputOutputHDF5.o \
 ../IO/InputOutputXML.o

SplineObjs = ../Splines/BicubicSpline.o ../Splines/CubicSpline.o              \
             ../Splines/Grid.o ../Splines/MyTricubicSpline.o                  \
             ../Splines/QuinticSpline.o ../Splines/QuinticSplines.o 

MiscObjs =   ../Integration/GKIntegration.o ../Fitting/Fitting.o              \
             ../MatrixOps/MatrixOps.o ../FFT/FFT.o ../MPI/Communication.o

all:    $(Objs) $(Objs2) TestPW TestPW2 BandStructure Hamiltonian.s

TestPW:  $(TestObjs) $(Objs) $(PHObjs) $(IOObjs) $(SplineObjs) $(MiscObjs)
	$(LD) -o $@ $(TestObjs) $(Objs) $(PHObjs) $(IOObjs) $(SplineObjs) \
        $(MiscObjs) $(LIBS)

TestPW2:  $(TestObjs2) $(Objs2) $(PHObjs) $(IOObjs) $(SplineObjs) $(MiscObjs)
	$(LD) -o $@ $(TestObjs2) $(Objs2) $(PHObjs) $(IOObjs) $(SplineObjs) \
        $(MiscObjs) $(LIBS)

TestPW2:  $(TestObjs2) $(Objs2) $(PHObjs) $(IOObjs) $(SplineObjs) $(MiscObjs)
	$(LD) -o $@ $(TestObjs2) $(Objs2) $(PHObjs) $(IOObjs) $(SplineObjs) \
        $(MiscObjs) $(LIBS)

BandStructure:  $(BSObjs) $(Objs2) $(PHObjs) $(IOObjs) $(SplineObjs) \
                $(MiscObjs)
	$(LD) -o $@ $(BSObjs) $(Objs2) $(PHObjs) $(IOObjs) $(SplineObjs) \
        $(MiscObjs) $(LIBS)



Hamiltonian.s:  Hamiltonian.cc
	$(CC) -S -O3 -ffast-math -g  $(DEFS) $(INCL) -o $*.s $<

clean:
	rm -f *.o

.cc.o: $(HEADERS)
	$(CC) $(CCFLAGS) $(DEFS) $(INCL) -o $*.o $< 
.f.o:
	g77 -c -o $*.o $<

PASS_DEFS = "CC=${CC}" "LD=${LD}" "CCFLAGS=${CCFLAGS}" "DEFS=${DEFS}" "INCL=${INCL}" "LIBS=${LIBS}" "F77=${F77}"\
  "MAKECC=${MAKECC}"

newmake:
	+$(MAKE) -f template.make Makefile FRC=force_rebuild $(PASS_DEFS)

Makefile:	$(FRC)
	rm -f $@
	cp template.make $@
	echo 'Automatically generated dependency list:' >> $@
	$(MAKECC) $(CCFLAGS) $(INCL) $(DEFS) -MM $(SOURCES) >> $@
	chmod -w $@


force_rebuild:
