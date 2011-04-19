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

#include "../Blitz.h"
#include "HeAzizPot.h"

const double HeAzizPot::bb[7] = 
{ 
   0,
   1.263030e6,
   1.399649e6,
  -8.389601e5,
   7.426020e6,
  -2.006420e6,
   8.570426e5
};

double HeAzizPot::V (double rv)
{
  double d,eps,a,alpha,c6,c8,c10;
  double rm,beta;
  double a2i,a6i;
  double f;
  double pott;
  double beta1;
  d=1.241314;
  eps=10.8;
  a=544850.4;
  alpha=13.353384;
  c6=1.3732412;
  c8=0.4253785;
  c10=0.178100;
  rm=2.9673;
  beta=0.0;
  beta1=3.32316;
  if (rv>1.828){
    double x=rv/rm;
    a2i=1.0/(x*x);
    a6i=a2i*a2i*a2i;
    if (x<d){
      f=exp(-pow((d/x-1.0),2));
    }
    else {
      f=1.0;
    }
    pott=eps*(a*exp(-alpha*x+beta*x*x)-f*a6i*(c6+a2i*(c8+a2i*c10)));
  }
  else {
    double x=rv*1.8897266;
    double sum=0;
    for (int i=1;i<=6;i++){
      sum=bb[7-i]+x*sum;
    }
    pott=sum*exp(-beta1*x)/x;
  }
  return pott;
}

double HeAzizPot::dVdr (double r)
{
  return 0;
}

double HeAzizPot::d2Vdr2 (double r)
{
  return 0;
}

void HeAzizPot::Write(IOSectionClass &out)
{
  out.WriteVar("Type", "HeAziz");
}

void HeAzizPot::Read(IOSectionClass &in)
{
  // do nothing -- no free parameters.
}


