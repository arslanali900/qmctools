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

#ifndef POLYNOMIAL_SET_CLASS_H
#define POLYNOMIAL_SET_CLASS_H
#include "Polynomial.h"

class WeightFuncClass 
{
public:
  virtual double operator()(double x) = 0;
};


class PolynomialSetClass
{
  Array<PolynomialClass,1> P;
public:
  inline double operator()(int k, double x)
  { return P(k)(x);}
  inline PolynomialClass& operator()(int k)
  { return P(k);}

  void MakeOrthoSet(int order, double xmin, double xmax, 
		    WeightFuncClass &wf);
};


#endif
