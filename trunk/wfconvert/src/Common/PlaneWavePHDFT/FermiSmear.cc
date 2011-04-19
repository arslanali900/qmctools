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

#include "FermiSmear.h"
#include "../SpecialFunctions/HermitePoly.h"

void
MethfesselPaxton::SetOrder(int order)
{
  Order = order;
  A.resize(order+1);
  H.resize(2*(order+1));
  double nfact = 1.0;
  double sign = 1.0;
  double sqrtPi = sqrt(M_PI);
  double four2n = 1.0;
  for (int n=0; n<=Order; n++) {
    A(n) = sign/(nfact*four2n*sqrtPi);
    nfact *= (double)(n+1);
    sign = -sign;
    four2n *= 4.0;
  }
}
  
      

double
MethfesselPaxton::D(double E, double mu)
{
  double x = (E-mu)/Width;
  HermitePoly (x, H);
  double Dn = 0.0;
  for (int n=0; n<=Order; n++)
    Dn += A(n)*H(2*n);
  Dn *= exp(-x*x);
  return Dn;
}

double
MethfesselPaxton::S(double E, double mu)
{
  double x = (E-mu)/Width;
  HermitePoly (x, H);
  double Sn = 0.0;
  for (int n=1; n<=Order; n++)
    Sn += A(n)*H(2*n-1);
  Sn *= exp(-x*x);
  Sn += 0.5*erfc(x);
  
  return Sn;
}
  
MethfesselPaxton& 
MethfesselPaxton::operator=(const MethfesselPaxton& mp)
{
  Order = mp.Order;
  Width = mp.Width;
  A.resize(mp.A.shape());  A = mp.A;
  H.resize(mp.H.shape());  H = mp.H;
  return *this;
}
