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

#include "LegendrePoly.h"
#include "SpecialFunctions.h"


void LegendreTest()
{
  Array<double,1> Pn(50);
  for (int n=0; n<50; n++) {
    for (double x=-1.0; x<1.0; x+=0.001) {
      double l1 = Legendre(n,x);
      LegendrePoly(x, Pn);
      double l2 = Pn(n);
      fprintf (stderr, "%7.4f %21.14e %21.14e\n", x, l1, l2);
    }
  }
}

main()
{
  LegendreTest();

}


