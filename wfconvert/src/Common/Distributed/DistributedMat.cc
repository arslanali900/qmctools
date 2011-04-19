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

void DistributedMat::AllGather()
{
  int maxElements = NumElements(0);
  int numProcs = MyComm.NumProcs();
  int row, col;

  // Allocate send an receive buffers
  Array<double,1> SendBuffer(maxElements),
    RecvBuffer(maxElements*numProcs);
  
  // Now fill send buffer with my elements
  for (int i=0; i<MyNumElements(); i++)
    {
      MyElement(i, row, col);
      SendBuffer(i) = Mat(row,col);
    }

  // Do collective communication
  MyComm.AllGather(SendBuffer, RecvBuffer);

  // Now, put elements where they should be
  for (int proc=0; proc<numProcs; proc++)
    {
      for (int i=0; i<NumElements(proc); i++)
	{
	  Element(proc, i, row, col);
	  Mat(row,col) = RecvBuffer(i+maxElements*proc);
	}
    }
}



void DistributedSymmMat::AllGather()
{
  int maxElements = NumElements(0);
  int numProcs = MyComm.NumProcs();
  int row, col;

  // Allocate send an receive buffers
  Array<double,1> SendBuffer(maxElements),
    RecvBuffer(maxElements*numProcs);
  
  // Now fill send buffer with my elements
  for (int i=0; i<MyNumElements(); i++)
    {
      MyElement(i, row, col);
      SendBuffer(i) = (*this)(row,col);
    }

  // Do collective communication
  MyComm.AllGather(SendBuffer, RecvBuffer);

  // Now, put elements where they should be
  for (int proc=0; proc<numProcs; proc++)
    {
      for (int i=0; i<NumElements(proc); i++)
	{
	  Element(proc, i, row, col);
	  (*this)(row,col) = RecvBuffer(i+maxElements*proc);
	}
    }
}



void DistributedArray3::AllGather()
{
  int maxElements = NumElements(0);
  int numProcs = MyComm.NumProcs();
  int row, col;

  int depth = Mat.extent(2);
  // Allocate send an receive buffers
  Array<double,1> SendBuffer(maxElements*depth),
    RecvBuffer(maxElements*numProcs*depth);
  
  // Now fill send buffer with my elements
  for (int i=0; i<MyNumElements(); i++) {
    MyElement(i, row, col);
    for (int k=0; k<depth; k++)
      SendBuffer(i*depth+k) = Mat(row,col,k);
  }

  // Do collective communication
  MyComm.AllGather(SendBuffer, RecvBuffer);

  // Now, put elements where they should be
  for (int proc=0; proc<numProcs; proc++)
    for (int i=0; i<NumElements(proc); i++) {
      Element(proc, i, row, col);
      for (int k=0; k<depth; k++)
	Mat(row,col,k) = RecvBuffer((i+maxElements*proc)*depth+k);
    }
}



void DistributedArray3b::AllGather()
{
  int maxElements = NumElements(0);
  int numProcs = MyComm.NumProcs();

  // Allocate send an receive buffers
  Array<double,1> SendBuffer(maxElements),
    RecvBuffer(maxElements*numProcs);
  
  // Now fill send buffer with my elements
  for (int elem=0; elem<MyNumElements(); elem++) {
    int i, j, k;
    MyElement(elem, i, j, k);
    SendBuffer(elem) = Mat(i,j,k);
  }

  // Do collective communication
  MyComm.AllGather(SendBuffer, RecvBuffer);

  // Now, put elements where they should be
  for (int proc=0; proc<numProcs; proc++)
    for (int elem=0; elem<NumElements(proc); elem++) {
      int i,j,k;
      Element(proc, elem, i, j, k);
      Mat(i,j,k) = RecvBuffer(elem+maxElements*proc);
    }
}




/* main(int argc, char **argv)
{
  COMM::Init(argc, argv);
  TestDistributedMat();
  COMM::Finalize();
}
*/
