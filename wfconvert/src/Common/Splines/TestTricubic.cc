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

#include "TricubicSpline.h"

main()
{
  LinearGrid xgrid(0.0, 5.0, 40);
  LinearGrid ygrid(0.0, 5.0, 40);
  LinearGrid zgrid(0.0, 5.0, 40);
  Array<double,3> f(40,40,40);
  
  f = 0.0;

  for (int ix=0; ix<xgrid.NumPoints; ix++)
    for (int iy=0; iy<ygrid.NumPoints; iy++)
      for (int iz=0; iz<zgrid.NumPoints; iz++) {
	double x = xgrid(ix);	
	double y = ygrid(iy);  
	double z = zgrid(iz);
	f(ix,iy,iz) = sin(2.0*x)*sin(2.0*y)*sin(2.0*z);
      }

  TricubicSpline spline(&xgrid, &ygrid, &zgrid, f);  
  LinearGrid xgrid2(0.0, 4.9, 51);
  LinearGrid ygrid2(0.0, 4.9, 51);
  LinearGrid zgrid2(0.0, 4.5, 51);

  for (int ix=0; ix<xgrid2.NumPoints; ix++)
    for (int iy=0; iy<ygrid2.NumPoints; iy++)
      for (int iz=0; iz<zgrid2.NumPoints; iz++) {
	double x = xgrid2(ix);	
	double y = ygrid2(iy);  
	double z = zgrid2(iz);
	double ex = sin(2.0*x)*sin(2.0*y)*sin(2.0*z);
	double sp = spline (x,y,z);  
	//fprintf (stderr, "%1.12e %1.12e\n", ex, sp);
      }
  
  // Speed test
  clock_t start, end;
  start = clock();
  for (int i=0; i<10000000; i++)
    double a = spline(1.0, 1.2, 1.0);
  end = clock();

  double rate = 1.0e7/((double)(end-start)/1.0e6);
  cerr << "Rate = " << rate << " per second.\n";

}
