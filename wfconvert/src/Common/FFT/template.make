SOURCES = FFT.cc TestFFT.cc

Objs = FFT.o

TestObjs  = TestFFT.o

all:    $(Objs) TestFFT

TestFFT:  $(TestObjs) $(Objs)
	$(LD) -o $@ $(TestObjs) $(Objs) $(LIBS)

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
