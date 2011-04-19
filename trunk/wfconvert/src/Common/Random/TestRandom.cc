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

#include "Random.h"
#include <time.h>

void SpeedTest()
{
  CommunicatorClass comm;
  comm.SetWorld();
  RandomClass rand(comm);
  rand.Init(314159);
  
  clock_t start1, end1, start2, end2;
  start1 = clock();
  for (int i=0; i<10000000; i++) 
    double d = rand.LocalGaussian(1.0);
  end1 = clock();

  start2 = clock();
  for (int i=0; i<10000000; i++) 
    double d = rand.LocalGaussian2(1.0);
  end2 = clock();

  cerr << "LocalGaussian time:   " << 
    (double)(end1-start1)/(double)CLOCKS_PER_SEC << endl;
  cerr << "LocalGaussian2 time:  " << 
    (double)(end2-start2)/(double)CLOCKS_PER_SEC << endl;
}



int main(int argc, char** argv)
{
#ifdef USE_MPI
  MPI_Init (&argc, &argv);
#endif
  CommunicatorClass comm;

  comm.SetWorld();
  RandomClass rand(comm);
  rand.Init(314159);
  int numLocal;
  if (comm.MyProc()==0) 
    numLocal = 4;
  else
    numLocal = 6;
  for (int i=0; i<numLocal; i++)
      cerr << "Proc " << comm.MyProc() << "  local random num = " 
	   << rand.Local() << endl;

  cerr << "Proc " << comm.MyProc() << " common random num = " 
       << rand.Common() << endl;

  SpeedTest();
#ifdef USE_MPI
  MPI_Finalize();
#endif
}

  

  
