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

#include "FastMultiCubicSpline.h"
#include "CubicSpline.h"
#include "MyTricubicSpline.h"

#include<time.h>


void
TestFastMultiCubic()
{
  Array<double,2> y(500,5);
  for (int i=0; i<y.size(); i++) {
    double x = 2.0*M_PI*(double)i/(double)(y.size()-1);
    for (int j=0; j<5; j++)
      y(i,j) = sin(x)*cos((double)j);
  }
  y(19) = y(0);
  FastMultiCubicSpline spline;
  spline.Init (0.0, 2.0*M_PI, y, true);
  FILE *fout = fopen ("FastCubic.dat", "w");
  Array<double,1> val(5);
  for (double x=0.0; x<= 2.0*M_PI; x+=0.001) {
    spline(x, val);
    fprintf (fout, "%1.18e ", x);
    for (int j=0; j<val.size(); j++) 
      fprintf (fout, "%1.18e %1.18e", val(j), sin(x)*cos((double)j));
    fprintf (fout, "\n");
  }
  fclose(fout);
}

void
FastMultiCubicSpeed()
{
  Array<double,2> y(100,5);
  for (int i=0; i<y.extent(0); i++) {
    for (int j=0; j<y.extent(1); j++) {
      double x = 2.0*M_PI*(double)i/(double)(y.size()-1);
      y(i) = sin(x)*cos((double)j);
    }
  }
  FastMultiCubicSpline spline;
  spline.Init (0.0, 2.0*M_PI, y);

  clock_t start, end, rstart, rend;

  const int numSamples = 100000000;  
  rstart = clock();
  for (int i=0; i<numSamples; i++) {
    double x = 2.0*M_PI*drand48();
  }
  rend = clock();

  Array<double,1> val(y.extent(1));
  start = clock();
  for (int i=0; i<numSamples; i++) {
    double x = 2.0*M_PI*drand48();
    spline(x, val);
  }
  end = clock();

  double time = (double)((end-start)-(rend-rstart))/(double)CLOCKS_PER_SEC;

  fprintf (stderr, "Time for 1 multi spline evaluation = %1.5e\n",
	   time/(double)numSamples);
}


void
CubicSpeed()
{
  Array<double,1> y(500);
  LinearGrid xGrid (0.0, 2.0*M_PI, 500);
  for (int i=0; i<y.size(); i++) {
    double x = 2.0*M_PI*(double)i/(double)(y.size()-1);
    y(i) = sin(x);
  }
  CubicSpline spline (&xGrid, y);

  clock_t start, end, rstart, rend;

  const int numSamples = 100000000;  
  rstart = clock();
  for (int i=0; i<numSamples; i++) {
    double x = 2.0*M_PI*drand48();
  }
  rend = clock();

  start = clock();
  for (int i=0; i<numSamples; i++) {
    double x = 2.0*M_PI*drand48();
    double p = spline(x);
  }
  end = clock();

  double time = (double)((end-start)-(rend-rstart))/(double)CLOCKS_PER_SEC;


  fprintf (stderr, "Time for 1 spline evaluation = %1.5e\n",
	   time/(double)numSamples);
}

void
TricubicSpeed()
{
  Array<double,1> y(500);
  LinearGrid xGrid (0.0, 2.0*M_PI, 500);
  LinearGrid yGrid (0.0, 1.0, 50);
  LinearGrid zGrid (0.0, 1.0, 50);
  Array<double,3> F(500,50,50);
  for (int i=0; i<F.extent(0); i++) {
    double x = xGrid(i);
    for (int j=0; j<F.extent(1); j++) {
      double y = yGrid(j);
      for (int k=0; k<F.extent(2); k++) {
	double z = zGrid(k);
	F(i,j,k) = sin(x)*(1.0+y*y*z*z*z);
      }
    }
  }
  MyTricubicSpline spline (&xGrid, &yGrid, &zGrid, F);

  clock_t start, end, rstart, rend;

  const int numSamples = 10000000;  
  rstart = clock();
  for (int i=0; i<numSamples; i++) {
    double x = 2.0*M_PI*drand48();
    double y = drand48();
    double z = drand48();
  }
  rend = clock();

  start = clock();
  for (int i=0; i<numSamples; i++) {
    double x = 2.0*M_PI*drand48();
    double y = drand48();
    double z = drand48();
    double p = spline(x,y,z);
  }
  end = clock();

  double time = (double)((end-start)-(rend-rstart))/(double)CLOCKS_PER_SEC;


  fprintf (stderr, "Time for 1 spline evaluation = %1.5e\n",
	   time/(double)numSamples);
}


main()
{
  //  TestFastMultiCubic();
  FastMultiCubicSpeed();
  //  CubicSpeed();
  //  TricubicSpeed();
} 
