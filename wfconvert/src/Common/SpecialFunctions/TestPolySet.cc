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

#include "PolynomialSet.h"
#include "LegendrePoly.h"
#include "HermitePoly.h"
#include "../Integration/GKIntegration.h"

class ConstClass : public WeightFuncClass
{
public:
  double operator()(double x)
  { return 1.0; }
};

class GaussianClass : public WeightFuncClass
{
public:
  double operator()(double x)
  { return exp(-x*x); }
};



void TestPolySet()
{
  const int order=10;
  PolynomialSetClass Pset;
  ConstClass wf;
  Pset.MakeOrthoSet(order, -1.0, 1.0, wf);
  for (double x=-1.0; x<=1.00001; x+=0.001) {
    fprintf (stderr, "%1.4e ", x);
    for (int n=0; n<=order; n++)
      fprintf (stderr, "%1.14e %1.14e ", Pset(n, x)/Pset(n,1.0),
	       LegendrePoly(n, x));
    fprintf (stderr, "\n");
  }
}


void TestPolySet2()
{
  const int order=10;
  PolynomialSetClass Pset;
  GaussianClass wf;
  Pset.MakeOrthoSet(order, -15.0, 15.0, wf);
  for (double x=-2.0; x<=2.00001; x+=0.001) {
    fprintf (stderr, "%1.4e ", x);
    for (int n=0; n<=order; n++)
      fprintf (stderr, "%1.14e %1.14e ", 
	       Pset(n, x)/Pset(n,1.0) * HermitePoly(n,1.0),
	       HermitePoly(n, x));
    fprintf (stderr, "\n");
  }
}


class OrthoCheckClass
{
  int n, m;
  PolynomialSetClass &Pset;
  WeightFuncClass &WF;
public:
  double operator()(double x)
  {
    return (WF(x)*Pset(m)(x)*Pset(n)(x));
  }

  OrthoCheckClass(int M, int N, 
		  WeightFuncClass &wf, 
		  PolynomialSetClass &pset) 
    : Pset(pset), WF(wf), m(M), n(N)
  {}
};

/// Test the orthonormality of a PolynomialSet by constructing an
/// overlap matrix.  The result should be the identity matrix.
void TestPolySet3()
{
  const int order=10;
  PolynomialSetClass Pset;
  GaussianClass wf;
  Pset.MakeOrthoSet(order, -15.0, 15.0, wf);

  for(int m=0; m<=order; m++) {
    for (int n=0; n<=order; n++) {
      OrthoCheckClass orthoCheck(m, n, wf, Pset);
      GKIntegration<OrthoCheckClass, GK15> Integrator(orthoCheck);
      double val = Integrator.Integrate(-15.0, 15.0, 1.0e-12);
      fprintf (stderr, "%1.16e ", val);
    }
    fprintf (stderr, "\n");
  }
}

main()
{
  TestPolySet2();
}
