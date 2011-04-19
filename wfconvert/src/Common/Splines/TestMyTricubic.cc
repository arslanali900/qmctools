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

#include "MyTricubicSpline.h"
#include "TricubicSpline.h"

main()
{
  int Nx = 10;
  int Ny = 13;
  int Nz = 17;
  Array<double,3> f(Nx,Ny,Nz);
  LinearGrid xgrid(0.0, 1.0, Nx);
  LinearGrid ygrid(0.0, 1.0, Ny);
  LinearGrid zgrid(0.0, 1.0, Nz);

  for (int ix=0; ix<Nx; ix++)
    for (int iy=0; iy<Ny; iy++)
      for (int iz=0; iz<Nz; iz++) {
	double x = xgrid(ix);
	double y = ygrid(iy);
	double z = zgrid(iz);
	
	f(ix,iy,iz) = sin(4.32*M_PI*x)*cos(3.32*M_PI*y)*sin(5.32*M_PI*z);
      }

  TricubicSpline TC(&xgrid, &ygrid, &zgrid, f);
  double g = TC(0.5, 0.5, 0.5);
  MyTricubicSpline MyTC(&xgrid, &ygrid, &zgrid, f);
  MyTC.Update();

//   for (int ix=0; ix<Nx; ix++)
//     for (int iy=0; iy<Ny; iy++)
//       for (int iz=0; iz<Nz; iz++) {
// 	double x = xgrid(ix);
// 	double y = ygrid(iy);
// 	double z = zgrid(iz);
	
// 	double dfdx = cos(4.32*M_PI*x)*cos(3.32*M_PI*y)*sin(5.32*M_PI*z);
// 	double mine = MyTC.F(ix,iy,iz)[1];
// 	fprintf (stderr, "%1.16e %1.16e\n", dfdx, mine);
//       }


  clock_t myStart, myEnd, fortStart, fortEnd;
  
  myStart = clock();
  for (int i=0; i<1000000; i++)
    double ex = MyTC(0.2, 0.3, 0.4);
  myEnd = clock();
  
  fortStart = clock();
  for (int i=0; i<1000000; i++)
    double ex = TC(0.2, 0.3, 0.4);
  fortEnd = clock();

  double fortRate = 1e12/(double)(fortEnd-fortStart);
  double myRate = 1e12/(double)(myEnd-myStart);

  cerr << "Fortran rate: " << fortRate << " per second\n";
  cerr << "My rate:      " << myRate << " per second\n";
 


  for (double x=0.0; x < 1.0; x+=0.001) {
    double y = 0.3221;
    double z = 0.481;
    double fort = TC(x, y, z);
    double mine = MyTC(x, y, z);
    double ex = sin(4.32*M_PI*x)*cos(3.32*M_PI*y)*sin(5.32*M_PI*z);
    fprintf (stderr, "%1.16e %1.16e %1.16e %1.16e\n", 
	     x, ex, fort, mine);
  }

}
