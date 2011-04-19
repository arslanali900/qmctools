SOURCES = SpecialFunctions.cc HermitePoly.cc TestLegendre.cc TestPoly.cc PolynomialSet.cc TestPolySet.cc

all:	SpecialFunctions.o HermitePoly.o TestLegendre.o TestPoly.o PolynomialSet.o TestPolySet.o TestHermite TestLegendre TestPoly TestPolySet

TestHermite: HermitePoly.o
	$(LD) -o TestHermite HermitePoly.o $(LIBS)

TestLegendre: TestLegendre.o
	$(LD) -o TestLegendre TestLegendre.o SpecialFunctions.o $(LIBS)

TestPoly: TestPoly.o
	$(LD) -o TestPoly TestPoly.o $(LIBS)

TestPolySet: TestPolySet.o PolynomialSet.o
	$(LD) -o TestPolySet TestPolySet.o PolynomialSet.o ../Integration/GKIntegration.o $(LIBS)

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
	$(MAKECC) $(CCFLAGS) $(INCL) -M $(SOURCES) >> $@
	chmod -w $@


force_rebuild:
