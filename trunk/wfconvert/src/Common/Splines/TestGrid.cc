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

#include "Grid.h"
#include "../IO/IO.h"

void TestOptimal()
{
  OptimalGrid rgrid(600, 15.0, 0.0039);
  double delta = rgrid(1) - rgrid(0);

  for (int i=0; i<rgrid.NumPoints; i++) {
    double r = rgrid(i) + 0.1*delta;
    int j = rgrid.ReverseMap(r);
    if (i != j)
      cerr << "Error in Optimal ReverseMap!\n"
	   << " i = " << i << "j = " << j << endl;
  }
}


void TestOptimal2()
{
  OptimalGrid2 rgrid(1.0e-6, 15, 10.0, 300);
  double delta = rgrid(1) - rgrid(0);

  for (int i=0; i<rgrid.NumPoints; i++) {
    double r = rgrid(i) + 0.1*delta;
    int j = rgrid.ReverseMap(r);
    if (i != j)
      cerr << "Error in Optimal2 ReverseMap!\n"
	   << " i = " << i << "j = " << j << endl;
  }
}


void TestCluster()
{
  ClusterGrid rgrid(1.0e-6, 15, 1.0/10.0, 300);
  double delta = rgrid(1) - rgrid(0);

  for (int i=0; i<rgrid.NumPoints; i++) {
    double r = rgrid(i) + 0.1*delta;
    int j = rgrid.ReverseMap(r);
    if (i != j)
      cerr << "Error in Cluster ReverseMap!\n"
	   << " i = " << i << "j = " << j << endl;
  }
  for (int i=0; i<rgrid.NumPoints; i++)
    fprintf (stdout, "%1.16e\n", rgrid(i));
}


void TestGeneral()
{
  IOSectionClass in;
  in.OpenFile("Na_HF_NLP.h5");
  in.OpenSection("Grid");
  Grid *grid = ReadGrid(in);
  int lo = grid->ReverseMap(3.2);
  cerr << "Grid(lo) = " << (*grid)(lo) << endl;
  cerr << "Grid(hi) = " << (*grid)(lo+1) << endl;
  in.CloseFile();
}

#include <time.h>

void TestSpeed()
{
  clock_t cstart, cend, gstart, gend, rstart, rend;
  CenterGrid center;
  GeneralGrid general;
  const int N = 100000000;
  int csum=0, gsum=0;

  center.Init  (-10.0, 10.0, 15.0, 200);
  general.Init (center.Points());

  rstart = clock();
  for (int i=0; i<N; i++) {
    double r = -10.0 + 20.0*drand48();
  }
  rend = clock();
  
  cstart = clock();
  for (int i=0; i<N; i++) {
    double r = -10.0 + 20.0*drand48();
    csum += center.ReverseMap (r);
  }
  cend = clock();
  
  gstart = clock();
  for (int i=0; i<N; i++) {
    double r = -10.0 + 20.0*drand48();
    gsum += general.ReverseMap (r);
  }
  gend = clock();
  cerr << "csum = " << csum << "  gsum = " << gsum << endl;

  double centerTime = 
    (double)(cend - cstart + rstart - rend)/(double)CLOCKS_PER_SEC;
  double generalTime = 
    (double)(gend - gstart + rstart - rend)/(double)CLOCKS_PER_SEC;

  cerr << "CenterGrid time  = " << centerTime  << " seconds.\n";
  cerr << "GeneralGrid time = " << generalTime << " seconds.\n";
}


main()
{
  TestSpeed();
  //  TestGeneral();
//   TestOptimal();
//   TestOptimal2();
//   TestCluster();
}
