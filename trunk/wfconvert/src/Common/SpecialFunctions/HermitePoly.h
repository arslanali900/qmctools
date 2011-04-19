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

#ifndef HERMITE_H
#define HERMITE_H

#include "../Blitz.h"

inline void HermitePoly(double t, Array<double,1> &H)
{
  H(0) = 1.0;
  if (H.rows()>1)
    H(1) = 2.0*t;
  for (int n=2; n < H.size(); n++)
    H(n) = 2.0*t*H(n-1) - 2.0*(n-1)*H(n-2);
}

inline double HermitePoly(int m, double t)
{
  if (m == 0)
    return (1.0);
  else if (m == 1)
    return (2.0*t);
  else {
    double Hnm2 = 1.0;
    double Hnm1 = 2.0*t;
    double Hn;
    for (int n=2; n<=m; n++) {
      Hn = 2.0*t*Hnm1 - 2.0*(n-1)*Hnm2;
      Hnm2 = Hnm1;
      Hnm1 = Hn;
    }
    return (Hn);
  }
}


inline void TestHermitePoly()
{
  Array<double,1> Hm(5);
//   for (double x=-2.0; x<=2.0; x+= 0.01) {
//     HermitePoly(x, Hm);
//     fprintf (stdout, "%1.5e ", x);
//     for (int m=0; m<Hm.size(); m++)
//       fprintf (stdout, "%1.5e ", Hm(m));
//     fprintf (stdout, "\n");
//   }

  for (double x=-2.0; x<=2.0; x+= 0.01) {
    fprintf (stdout, "%1.5e ", x);
    for (int m=0; m<Hm.size(); m++)
      fprintf (stdout, "%1.5e ", HermitePoly(m,x));
    fprintf (stdout, "\n");
  }
}
    
#endif
