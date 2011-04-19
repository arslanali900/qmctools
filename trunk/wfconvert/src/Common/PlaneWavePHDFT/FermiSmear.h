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

#ifndef FERMI_SMEAR_H
#define FERMI_SMEAR_H

#include "../Blitz.h"

class FermiSmearClass
{
public:
  virtual double D(double E, double mu) = 0;
  virtual double S(double E, double mu) = 0;
};

class MethfesselPaxton : public FermiSmearClass
{
protected:
  int Order;
  double Width;
  Array<double,1> A;
  Array<double,1> H;
  double S0(double x);
public:
  void   SetOrder (int order);
  inline int    GetOrder ()             { return Order; }
  inline void   SetWidth (double width) { Width = width; }
  inline double GetWidth ()             { return Width; }
  double D(double E, double mu);
  double S(double E, double mu);
  MethfesselPaxton& operator=(const MethfesselPaxton& mp);
};

#endif
  
