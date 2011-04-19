/////////////////////////////////////////////////////////////
// Copyright (C) 2003-2006 Bryan Clark and Kenneth Esler   //
//                                                         //
// This program is free software; you can redistribute it  //
// and/or modify it under the terms of the GNU General     //
// Public License as published by the Free Software        //
// Foundation; either version 2 of the License, or         //
// (at your option) any later version.  This program is    //
// distributed in the hope that it will be useful, but     //
// WITHOUT ANY WARRANTY; without even the implied warranty //
// of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. //  
// See the GNU General Public License for more details.    //
// For more information, please see the PIMC++ Home Page:  //
//           http://pathintegrals.info                     //
/////////////////////////////////////////////////////////////

#include "IO.h"

using namespace IO;


void TestDoubleWrite(IOSectionClass &io)
{
  io.NewSection("Double");

  double x0 = 2.0;
  io.WriteVar("x0", x0);

  Array<double,1> x1write(5), x1read(5);
  x1write =  1.0, 2.0, 3.0, 4.0, 5.0 ;
  io.WriteVar("x1", x1write);
  io.ReadVar ("x1", x1read);
  for (int i=0; i<5; i++)
    assert (x1read(i) == x1write(i));

  Array<double,2> x2write(5,3), x2read(5,3);
  x2write =  
    1.0,   2.0,  3.0, 
    4.0,   5.0,  5.0,
    7.0,   8.0,  9.0,
    10.0, 11.0, 12.0,
    23.0, 14.0, 15.0;
  io.WriteVar("x2", x2write);
  io.ReadVar ("x2", x2read);
  for (int i=0; i<5; i++)
    for (int j=0; j<3; j++)
      assert(x2read(i,j) == x2write(i,j));
  
  IOVarBase* ioVar = io.GetVarPtr("x2");
  ioVar->Write(x1write,Range::all(), 0);
  x2write(Range::all(),0) = x1write;
  io.ReadVar("x2", x2read);
  for (int i=0; i<5; i++)
    for (int j=0; j<3; j++)
      assert(x2read(i,j) == x2write(i,j));
  
  ioVar->Resize(6);
  Array<double,1> x3write(3);
  x3write = 1.0, 2.0, 3.0;
  ioVar->Write(x3write,5,Range::all());
  io.ReadVar("x2", x2read);
  cerr << "extended x2 = " << x2read << endl;

  io.CloseSection();

  io.FlushFile();
}


void TestDoubleRead(IOSectionClass &io)
{
  io.OpenSection("Double");

  double x0;
  assert(io.ReadVar("x0", x0));
  assert (fabs(x0-2.0) < 1.0e-15);

  Array<double,1> x1;
  assert(io.ReadVar("x1", x1));
  for (int i=0; i<5; i++)
    assert (fabs(x1(i) - (i+1.0)) < 1.0e-15);

  Array<double,2> x2, x2test(5,3);
  x2test =  
    1.0,  2.0,  3.0, 
    2.0,  5.0,  5.0,
    3.0,  8.0,  9.0,
    4.0, 11.0, 12.0,
    5.0, 14.0, 15.0;
  
  assert(io.ReadVar("x2", x2));
  cerr << " x2 = " << x2 << endl;
  cerr << " x2test = " << x2test << endl;
  for (int i=0; i<5; i++)
    for (int j=0; j<3; j++) {
      assert (fabs(x2(i,j)-x2test(i,j)) < 1.0e-15);
    }

  Array<double,1> x2Range;
  IOVarBase* ioVar = io.GetVarPtr("x2");
  ioVar->Read(x2Range,Range::all(), 1);
  
  io.CloseSection();
}

void TestComplexWrite (IOSectionClass &out)
{
  Array<complex<double>,1> vals(100);
  for (int i=0; i<100; i++) {
    double phi = (double)i/99.0 * 2.0*M_PI;
    double re = cos(phi);
    double im = sin(phi);
    vals(i) = complex<double>(re,im);
  }
  out.WriteVar("CompTest", vals);
}


void TestComplexRead (IOSectionClass &out)
{
  Array<complex<double>,1> vals;
  out.ReadVar("CompTest", vals);
}


void
TestAppend (IOSectionClass &out)
{
  out.NewSection("Append");
  Array<double,1> E(1);
  E(0) = 5.0;
  out.WriteVar ("E", E);
  IOVarBase *Evar = out.GetVarPtr("E");
  for (double e=6.0; e<20.0; e++)
    Evar->Append(e);
  out.CloseSection();
}



void TestHDF5()
{
  IOSectionClass out;
  
  out.NewFile("IOtest.h5");
  Array<double,2> myVar(2,3), myVar2;
  myVar(0,0)=1.0; myVar(0,1)=2.0; myVar(0,2)=3.0;
  myVar(1,0)=4.0; myVar(1,1)=5.0; myVar(1,2)=6.0;

  string TestString = "Hello World";
  out.WriteVar("TestString", TestString);

  out.WriteVar("Test", myVar);
  double a = 3.0;
  out.WriteVar ("a", a);

  out.CloseFile();

  Array<double,1> x(3);
  x(0) = 7.0; x(1) = 8.0; x(2)=9.0;

  IOSectionClass in;
  assert(in.OpenFile ("IOtest.h5"));
  IOVarBase *test = in.GetVarPtr("Test");
  nilArraySection n;
  test->Write(x,0,Range::all());
  in.PrintTree();
  in.ReadVar("Test", myVar2);
  cerr << "myVar2 = " << myVar2 << endl;
  Array<double,2> RangeVar;
  test->Read(RangeVar, Range(0,1), Range(0,2,2));
  cerr << "RangeVar = " << RangeVar << endl;
  in.CloseFile();

}


void TestASCII()
{
  cerr << "\n\nASCII Test:\n\n";
  Array<double,2> myVar;
  Array<double,1> y(2);
  Array<double,1> x(3);
  x(0) = 7.0; x(1) = 8.0; x(2)=9.0;

  IOSectionClass in;
  assert(in.OpenFile ("IOtest.txt"));
  in.PrintTree();
  IOVarBase *test = in.GetVarPtr("Test");
  nilArraySection n;
  test->Write(x,0,Range::all());
  test->Read(y,1,Range(1,2));
  cerr << "y = " << y << endl;
  in.PrintTree();
  in.ReadVar("Test", myVar);
  test->Read(myVar);
  cerr << "myVar = " << myVar << endl;
  //  in.CloseFile();

}


main()
{
  TestHDF5();
  TestASCII();
  IOSectionClass ioHDF5;
  ioHDF5.NewFile ("TestHDF5.h5");
  TestDoubleWrite(ioHDF5);
  TestDoubleRead(ioHDF5);
  TestAppend(ioHDF5);
  TestComplexWrite (ioHDF5);
  TestComplexRead  (ioHDF5);


  IOSectionClass ioASCII;
  ioASCII.NewFile ("TestASCII.txt");
  TestDoubleWrite(ioASCII);
  TestDoubleRead(ioASCII);
  TestComplexWrite (ioASCII);
  TestComplexRead  (ioASCII);


//   IOSectionClass hack;
//   hack.OpenFile("/home/esler/PairActions/Na_HF_NLPP.square.beta16.0.h5");
//   hack.OpenSection("Potential");
//   hack.WriteVar("Type", "General");
//   hack.CloseFile();

}
