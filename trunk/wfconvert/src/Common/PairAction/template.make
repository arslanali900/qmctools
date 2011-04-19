SOURCES = PAszFit.cc PAcoulombFit.cc PAcoulombBCFit.cc PAclassicalFit.cc PAsFit.cc PAtricubicFit.cc PAzeroFit.cc DavidPAClass.cc PAMonopoleFit.cc PADipoleFit.cc PATripoleFit.cc PAtricubicFit2.cc

#all:	PAszFit.o PAcoulombFit.o PAcoulombBCFit.o PAclassicalFit.o PAsFit.o PAtricubicFit.o PAzeroFit.o DavidPAClass.o
all:	 PAcoulombBCFit.o PAtricubicFit.o PAclassicalFit.o PAzeroFit.o DavidPAClass.o PAMonopoleFit.o PADipoleFit.o PATripoleFit.o PAtricubicFit2.o

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
