SOURCES = CubicSpline.cc Grid.cc BicubicSpline.cc TestBicubic.cc TestGrid.cc TestTricubic.cc MyTricubicSpline.cc TestMyTricubic.cc QuinticSpline.cc TestQuintic.cc DyutimanTest.cc MultiTricubicSpline.cc TestMultiTricubicSpline.cc TestMultiTricubicSpline2.cc MultiTricubicSpline3.cc TestMultiTricubicSpline3.cc PeriodicSpline.cc TestPeriodic.cc ComplexMultiTricubicSpline.cc TestComplexMultiTricubicSpline.cc ComplexMultiTricubicSpline2.cc TestComplexMultiTricubicSpline2.cc 

IOobjs = ../IO/InputOutput.o ../IO/InputOutputHDF5.o ../IO/InputOutputASCII.o  ../IO/InputOutputXML.o

F77Objs = fortran/evtricub.o  fortran/herm3ev.o  fortran/mktricubw.o  fortran/tcspline.o fortran/ibc_ck.o fortran/splinck.o fortran/zonfind.o fortran/tcspeval.o fortran/v_spline.o fortran/bcspline.o fortran/bcspeval.o

Objs:	MyTricubicSpline.o CubicSpline.o QuinticSpline.o DyutimanTest.o MultiTricubicSpline.o MultiTricubicSpline3.o PeriodicSpline.o ComplexMultiTricubicSpline.o ComplexMultiTricubicSpline2.o 

all:	FortranObjs TestBicubic TestGrid  TestMyTricubic Objs TestQuintic DyutimanTest TestMultiTricubicSpline TestMultiTricubicSpline2 TestMultiTricubicSpline3 TestPeriodic TestComplexMultiTricubicSpline TestComplexMultiTricubicSpline2 

TestQuintic:	QuinticSpline.o QuinticSplines.o TestQuintic.o
	$(LD) -o TestQuintic QuinticSpline.o QuinticSplines.o TestQuintic.o $(LIBS)

TestBicubic:	CubicSpline.o Grid.o BicubicSpline.o TestBicubic.o 
	$(LD) -o TestBicubic CubicSpline.o Grid.o BicubicSpline.o TestBicubic.o $(IOobjs) $(LIBS)
TestPeriodic:  PeriodicSpline.o TestPeriodic.o Grid.o
	$(LD) -o TestPeriodic  PeriodicSpline.o TestPeriodic.o Grid.o $(IOobjs) $(LIBS)


TestMyTricubic:	CubicSpline.o Grid.o  TestMyTricubic.o MyTricubicSpline.o
	$(LD) -o TestMyTricubic CubicSpline.o Grid.o TestMyTricubic.o MyTricubicSpline.o $(IOobjs) -Lfortran -lpspline $(LIBS)

TestMultiTricubicSpline:	Grid.o  TestMultiTricubicSpline.o MyTricubicSpline.o MultiTricubicSpline.o
	$(LD) -o TestMultiTricubicSpline Grid.o MultiTricubicSpline.o MyTricubicSpline.o TestMultiTricubicSpline.o $(IOobjs) $(LIBS)


TestMultiTricubicSpline2:	Grid.o  TestMultiTricubicSpline2.o MyTricubicSpline.o
	$(LD) -o TestMultiTricubicSpline2 Grid.o MyTricubicSpline.o TestMultiTricubicSpline2.o $(IOobjs) $(LIBS)

TestMultiTricubicSpline3:	Grid.o  TestMultiTricubicSpline3.o MyTricubicSpline.o MultiTricubicSpline3.o r3spline.o
	$(LD) -o TestMultiTricubicSpline3 Grid.o MyTricubicSpline.o MultiTricubicSpline3.o TestMultiTricubicSpline3.o r3spline.o $(IOobjs) $(LIBS)

TestComplexMultiTricubicSpline:	Grid.o  TestComplexMultiTricubicSpline.o MyTricubicSpline.o ComplexMultiTricubicSpline.o
	$(LD) -o TestComplexMultiTricubicSpline Grid.o MyTricubicSpline.o ComplexMultiTricubicSpline.o TestComplexMultiTricubicSpline.o $(IOobjs) $(LIBS)

TestComplexMultiTricubicSpline2:	Grid.o  TestComplexMultiTricubicSpline2.o MyTricubicSpline.o ComplexMultiTricubicSpline2.o z3spline.o
	$(LD) -o TestComplexMultiTricubicSpline2 Grid.o MyTricubicSpline.o ComplexMultiTricubicSpline2.o TestComplexMultiTricubicSpline2.o z3spline.o $(IOobjs) $(LIBS)


DyutimanTest:	CubicSpline.o Grid.o  TestMyTricubic.o MyTricubicSpline.o
	$(LD) -o DyutimanTest CubicSpline.o Grid.o DyutimanTest.o MyTricubicSpline.o $(IOobjs) -Lfortran -lpspline $(LIBS)


TestTricubic:	Grid.o TestTricubic.o FortranObjs TestMyTricubic
	$(LD) -o TestTricubic Grid.o TestTricubic.o $(F77Objs) $(IOobjs) $(LIBS)

TestGrid:	 Grid.o  TestGrid.o 
	$(LD) -o TestGrid Grid.o TestGrid.o $(IOobjs) $(LIBS)

FortranObjs:
	cd fortran; $(MAKE)


clean:
	rm -f *.o

.cc.o: 
	$(CC) $(CCFLAGS) $(DEFS) $(INCL) -o $*.o $< 
.f.o:
	$(F77) $(F77FLAGS) -o $*.o $<

newmake:
	+$(MAKE) -f template.make Makefile FRC=force_rebuild

Makefile:	$(FRC)
	rm -f $@
	cp template.make $@
	echo 'Automatically generated dependency list:' >> $@
	$(MAKECC) $(CCFLAGS) $(INCL) $(DEFS) -MM $(SOURCES) >> $@
	chmod -w $@


force_rebuild:
