SOURCES = DistributedMat.cc TestDistributedMat.cc

all:	DistributedMat.o TestDistributedMat.o ../MPI/Communication.o Test

Test: DistributedMat.o TestDistributedMat.o ../MPI/Communication.o
	$(CC) -o Test TestDistributedMat.o DistributedMat.o ../MPI/Communication.o $(LIBS)

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
	$(MAKECC) $(CCFLAGS) $(INCL) $(DEFS) -MM $(SOURCES) >> $@
	chmod -w $@


force_rebuild:
