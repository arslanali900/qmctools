SOURCES = InputOutput.cc InputOutputHDF5.cc InputFile.cc  InputOutputASCII.cc TestXML.cc TestHDF5.cc

objs = InputOutput.o  InputOutputHDF5.o  InputOutputASCII.o InputOutputXML.o #InputFile.o  
all:	$(objs) TestXML TestHDF5



TestXML:	$(objs) TestXML.o
	$(LD) -o $@ $(objs) TestXML.o $(LIBS)

TestHDF5:	$(objs) TestHDF5.o
	$(LD) -o $@ $(objs) TestHDF5.o $(LIBS)

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
