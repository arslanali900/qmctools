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

#include "Polynomial.h"
#include <vector>

void TestPoly()
{
  PolynomialClass a(3), b(2), c(4);
  a[0]= 2.0; a[1]= -1.0; a[2] = 0.2; a[3]=1.5;
  b[0]= 0.5; b[1]= 2.1;  b[2]=-4.9;
  c[0]= 0.3; c[1]= 1.5;  c[2]=1.5; c[3]=-2.0; c[4]=0.2;
 
  double x = 3.14159;
  double aofx = a(x);
  double bofx = b(x);
  double cofx = c(x);
  
  PolynomialClass ab(a*b);
  PolynomialClass aplusb;
  aplusb = a+b;
 
  fprintf (stderr, "a(x)*b(x) = %1.16e\n", a(x)*b(x));
  fprintf (stderr, "(a*b)(x)  = %1.16e\n", ab(x));
  fprintf (stderr, "a(x)+b(x) = %1.16e\n", a(x)+b(x));
  fprintf (stderr, "(a+b)(x)  = %1.16e\n", aplusb(x));

}

void MakeTable()
{
  vector<PolynomialClass> u(4), d2u(4);
  PolynomialClass p1(3), p2(3), q1(3), q2(3);
  p1[0]=1.0; p1[1]=0.0; p1[2]=-3.0; p1[3]=2.0;
  p2[0]=0.0; p2[1]=0.0; p2[2]= 3.0; p2[3]=-2.0;
  q1[0]=0.0; q1[1]=1.0; q1[2]=-2.0; q1[3]=1.0;
  q2[0]=0.0; q2[1]=0.0; q2[2]=-1.0; q2[3]=1.0;
  
  u[0]=p1; u[1]=p2; u[2]=q1; u[3]=q2;

  PolynomialClass d2p1, d2p2, d2q1, d2q2;
  d2p1 = Deriv(Deriv(p1));
  d2p2 = Deriv(Deriv(p2));
  d2q1 = Deriv(Deriv(q1));
  d2q2 = Deriv(Deriv(q2));

  d2u[0]=d2p1; d2u[1]=d2p2; d2u[2]=d2q1; d2u[3]=d2q2;

  PolynomialClass uij, intuij;
  TinyMatrix<double,4,4> Amat;
  cerr << "Amat = " << endl;
  for (int i=0; i<4; i++) {
    for (int j=0; j<4; j++) {
      uij = u[i]*u[j];
      intuij = Integral(uij);
      Amat(i,j) = intuij(1.0)-intuij(0.0);
      fprintf (stderr, "%24.20f ", Amat(i,j));
    }
    fprintf (stderr, "\n");
  }
  TinyMatrix<double,4,4> d2Amat;
  cerr << "d2Amat = " << endl;
  for (int i=0; i<4; i++) {
    for (int j=0; j<4; j++) {
      uij = u[i]*d2u[j];
      intuij = Integral(uij);
      d2Amat(i,j) = intuij(1.0)-intuij(0.0);
      fprintf (stderr, "%24.20f ", d2Amat(i,j));
    }
    fprintf (stderr, "\n");
  }
  
}

main()
{
  TestPoly();
  MakeTable();
}
