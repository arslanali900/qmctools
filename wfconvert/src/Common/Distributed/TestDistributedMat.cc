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

#include "DistributedMat.h"

bool TestDistributedSymmMat()
{
  const int N=6;
  CommunicatorClass comm;
  comm.SetWorld();
  
  DistributedSymmMat  dmat(N,comm);
  int row, col;

  for (int i=0; i<dmat.MyNumElements(); i++)
    {
      dmat.MyElement(i,row,col);
      dmat(row,col) = 10.0*row+col;
    }
  dmat.AllGather();
  bool passed=true;
  for (int row=0; row<N; row++)
    for (int col=0; col<=row; col++)
      if (dmat(row,col) != (10.0*row+col))
	passed = false;
  if (!passed)
    cerr << "Error in TestDistributedSymmMat()\n";
  return (passed);
}

bool TestDistributedArray3()
{
  const int L=5;
  const int M=7;
  const int N=3;
  CommunicatorClass comm;
  comm.SetWorld();
  DistributedArray3 d3(L,M,N,comm);
  int row,col;
  for (int i=0; i<d3.MyNumElements(); i++) {
    d3.MyElement(i,row,col);
    for (int j=0; j<N; j++){
      d3(row,col,j) = 100.0*row + 10.0*col + (double)(j+1);
    }
  }
  d3.AllGather();
  bool passed = true;
  for (int i=0; i<L; i++)
    for (int j=0; j<M; j++)
      for (int k=0; k<N; k++) 
	if (d3(i,j,k) != (100.0*i + 10.0*j + (double)(k+1))) {
	  cerr << "Error in TestDistributedArray3 at i="
	       << i << " j=" << j << " k=" << k << endl;
	  cerr << d3(i,j,k) << endl;
	  passed = false;
	}
  if (!passed)
    cerr << "Error in TestDistributedArray3()\n";
  return passed;
}


bool TestDistributedArray3b()
{
  const int L=5;
  const int M=7;
  const int N=3;
  CommunicatorClass comm;
  comm.SetWorld();
  DistributedArray3b d3(L,M,N,comm);
  int row,col;
  for (int elem=0; elem<d3.MyNumElements(); elem++) {
    int i,j,k;
    d3.MyElement(elem,i,j,k);
    d3(i,j,k) = 100.0*i + 10.0*j + (double)(k+1);
  }
  d3.AllGather();
  bool passed = true;
  for (int i=0; i<L; i++)
    for (int j=0; j<M; j++)
      for (int k=0; k<N; k++) 
	if (d3(i,j,k) != (100.0*i + 10.0*j + (double)(k+1))) {
	  cerr << "Error in TestDistributedArray3b at i="
	       << i << " j=" << j << " k=" << k << endl;
	  cerr << d3(i,j,k) << endl;
	  passed = false;
	}
  if (!passed)
    cerr << "Error in TestDistributedArray3b()\n";
  return passed;
}


main(int argc, char **argv)
{
  COMM::Init(argc, argv);
  CommunicatorClass comm;
  comm.SetWorld();
  if (comm.NumProcs() == 1) {
    cerr << "Error:  This test must be run on more than one processor.\n";
    exit(-1);
  }
  bool passed;
  passed = TestDistributedSymmMat();
  passed = passed && TestDistributedArray3();
  passed = passed && TestDistributedArray3b();
  int MyProc = comm.MyProc();
  COMM::Finalize();
  if (MyProc == 0) {
    cout << "Distributed Matrices Test:\n";
    if (passed)
      cout << "  Passed.\n";
    else
      cout << "  Failed.\n";
  }      
  if (passed) 
    return (0);
  else 
    return (-1);
}
