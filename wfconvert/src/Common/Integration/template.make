SOURCES = Integrate.cc GKIntegration.cc HermiteQuad.cc TestHermite.cc TestIntegrate.cc

all:	Integrate.o GKIntegration.o HermiteQuad.o TestHermite.o TestIntegrate.o TestHermite TestIntegrate

TestIntegrate: TestIntegrate.o Integrate.o
	$(LD) -o TestIntegrate TestIntegrate.o Integrate.o $(IOobjs) $(LIBS)

TestHermite: TestHermite.o HermiteQuad.o
	$(LD) -o TestHermite  HermiteQuad.o TestHermite.o $(IOobjs) $(LIBS)

clean:
	rm -f *.o


.cc.o: $(HEADERS)
	$(CC) $(CCFLAGS) $(DEFS) $(INCL) -o $*.o $< 
.f.o:
	g77 -c -o $*.o $<

newmake:
	$(MAKE) -f template.make Makefile FRC=force_rebuild

Makefile:	$(FRC)
	rm -f $@
	cp template.make $@
	echo 'Automatically generated dependency list:' >> $@
	$(MAKECC) $(CCFLAGS) $(INCL) -MM $(SOURCES) >> $@
	chmod -w $@


force_rebuild:
