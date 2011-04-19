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

#include "PeriodicSpline.h"

void TestPeriodic()
{
  LinearGrid xGrid(0.0, 4.0, 5);
  Array<double,1> vals(5);
  vals(0) = 1.0;
  vals(1) = 0.0;
  vals(2) = -1.0;
  vals(3) = 0.0;
  vals(4) = 1.0;

  PeriodicSpline spline;
  spline.Init (&xGrid, vals);

  for (double x=0.0; x < 4.0; x+= 0.001)
    fprintf (stdout, "%1.16e %1.16e %1.16e\n", x, spline(x), spline.Deriv(x));

}


main()
{
  TestPeriodic();
}
