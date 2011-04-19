ppconvert: CubicSpline.o NLPPClass.o ParseCommand.o ParserClass.o XMLWriterClass2.o
	g++ -g -O2 -o ppconvert  NLPPClass.o CubicSpline.o ParseCommand.o \
                      ParserClass.o XMLWriterClass2.o


CubicSpline.o: CubicSpline.h CubicSpline.cc

NLPPClass.o:   CubicSpline.h ParseCommand.h ParserClass.h XMLWriterClass2.h\
               NLPPClass.h 

clean:
	rm -f *.o

.cc.o:
	g++ -g -O2 -c $<
	