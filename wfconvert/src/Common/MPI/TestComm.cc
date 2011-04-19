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

#include "Communication.h"

bool TestAllGatherRows()
{
  CommunicatorClass Comm;
  Comm.SetWorld();

  int cols = 3;
  int rows = 2*Comm.NumProcs()+1;
  Array<complex<double>,2> mat(rows,cols);
  mat = 0.0;
  int currRow = 0;
  for (int proc=0; proc<Comm.NumProcs(); proc++) {
    int procRows = rows/Comm.NumProcs() + ((rows%Comm.NumProcs())>proc);
    if (proc == Comm.MyProc()) {
      for (int row=currRow; row<(currRow+procRows); row++)
	for (int col=0; col<cols; col++)
	  mat(row,col) = (double)row + (double)col/10.0;
    }
    currRow += procRows;
  }
//   if (Comm.MyProc() == 0) 
//     cerr << "mat = " << mat << endl;

  Comm.AllGatherRows(mat);
  bool passed = true;
  for (int row=0; row<rows; row++)
    for (int col=0; col<cols; col++)
      passed = passed && 
	(fabs(mat(row,col).real()-((double)row + (double)col/10.0)) < 1.0e-12);
//   if (Comm.MyProc() == 0) 
//     cerr << "mat = " << mat << endl;
  return passed;
}


bool TestAllGatherVec()
{
  CommunicatorClass Comm;
  Comm.SetWorld();

  int elems = 2*Comm.NumProcs()+1;
  Array<double,1> vec(elems);
  vec = 0.0;
  int currElem = 0;
  for (int proc=0; proc<Comm.NumProcs(); proc++) {
    int procElems = elems/Comm.NumProcs() + ((elems%Comm.NumProcs())>proc);
    if (proc == Comm.MyProc()) {
      for (int elem=currElem; elem<(currElem+procElems); elem++)
	  vec(elem) = (double)elem;
    }
    currElem += procElems;
  }

  Comm.AllGatherVec(vec);
  bool passed = true;
  for (int elem=0; elem<elems; elem++)
    passed = passed && 
      (fabs(vec(elem)-(double)elem) < 1.0e-12);

  return passed;
}


#include <time.h>
#include <../Random/Random.h>
void 
TestSumSpeed()
{
  CommunicatorClass Comm;
  Comm.SetWorld();
  RandomClass random(Comm);
  random.Init();

  clock_t start, end;
  start = clock();
  for (int i=0; i<1000000; i++) {
    double x = random.Local();
    double y = Comm.AllSum(x);
  }
  end = clock();
  if (Comm.MyProc() == 0) {
    double t = (double)(end-start)/(double)CLOCKS_PER_SEC;
    cerr << "Sums per second = " << (1.0e6/t) << endl;
  }
}


main(int argc, char *argv[])
{
  COMM::Init(argc, argv);
  CommunicatorClass comm;
  comm.SetWorld();

  if (comm.MyProc()==0)
    cerr << "Number of procs = " << comm.NumProcs() << endl;
  bool passed = TestAllGatherRows();  
  if (comm.MyProc()==0)
    cerr << "AllGatherRows() check:  " 
	 << (passed ? "passed." : "failed.") << endl;

  passed = TestAllGatherVec();  
  if (comm.MyProc()==0)
    cerr << "AllGatherVec() check:   " 
	 << (passed ? "passed." : "failed.") << endl;
  TestSumSpeed();
}
