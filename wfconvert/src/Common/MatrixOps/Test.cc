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

#include "MatrixOps.h"

void TestDet()
{
  Array<double,2> C(10,10);
  fprintf (stderr, "C = \n");
  fprintf (stdout, "C = [ ");
  for (int i=0; i<10; i++) {
    for (int j=0; j<10; j++){
      C(i,j) = 2.0*drand48()-1.0;
      fprintf (stdout, "%1.16e ", C(i,j));
    }
    if (i < 9)
      fprintf (stdout, "; ");
    else
      fprintf (stdout, "]\n ");
  }
  fprintf (stderr, "det(C) = %1.16e\n", Determinant(C));
}

#include <time.h>

void TimeDet()
{
  Array<double,2> C(16,16);
  fprintf (stderr, "C = \n");
  fprintf (stdout, "C = [ ");
  for (int i=0; i<16; i++) 
    for (int j=0; j<16; j++)
      C(i,j) = 2.0*drand48()-1.0;
    
  clock_t start, end;

  start  = clock();
  for (int i=0; i<100000; i++) {
    for (int j=0; j<16; j++)
      for (int k=0; k<16; k++)
	C(j,k) = exp(-(double)(k-j)*(k-j));
    Determinant(C);
  }
  end = clock();

  fprintf (stderr, "Time = %5.2f sec.\n", (double)(end-start)/CLOCKS_PER_SEC);
  fprintf (stderr, "det(C) = %1.16e\n", Determinant(C));
}


bool TestInv()
{
  const int N = 5;
  Array<double,2> C(N,N), Cinv(N,N), eye(N,N);
  
  // Random matrix
  for (int i=0; i<N; i++) 
    for (int j=0; j<N; j++)
      C(i,j) = 2.0*drand48()-1.0;
  
  Cinv = C;
  double det = Determinant (C);
  double gjdet = GJInverse(Cinv);

  fprintf (stderr, "det   = %20.16e\n", det);
  fprintf (stderr, "gjdet = %20.16e\n", gjdet);

  eye = 0.0;
  for (int i=0; i<N; i++)
    for (int j=0; j<N; j++)
      for (int k=0; k<N; k++)
	eye(i,j) += Cinv(i,k)*C(k,j);

  eye = Cinv * C;
  bool failed = false;
  for (int i=0; i<N; i++)
    for (int j=0; j<N; j++)
      if (i == j)
	failed = (fabs(eye(i,j) - 1.0) > 1.0e-14);
      else
	failed = (fabs(eye(i,j)) > 1.0e-14);
  cerr << "eye = " << eye << endl;
  if (failed)
    cerr << "Failed inverse test.\n";
  else
    cerr << "Passed inverse test.\n";

  return !failed;
}

bool TestInvPartial()
{
  const int N = 5;
  Array<double,2> C(N,N), Cinv(N,N), eye(N,N);
  
  // Random matrix
  for (int i=0; i<N; i++) 
    for (int j=0; j<N; j++)
      C(i,j) = 2.0*drand48()-1.0;
  
  Cinv = C;
  double det = Determinant (C);
  double gjdet = GJInversePartial(Cinv);

  fprintf (stderr, "det   = %20.16e\n", det);
  fprintf (stderr, "gjdet = %20.16e\n", gjdet);

  eye = 0.0;
  for (int i=0; i<N; i++)
    for (int j=0; j<N; j++)
      for (int k=0; k<N; k++)
	eye(i,j) += Cinv(i,k)*C(k,j);

  eye = Cinv * C;
  bool failed = false;
  for (int i=0; i<N; i++)
    for (int j=0; j<N; j++)
      if (i == j)
	failed = (fabs(eye(i,j) - 1.0) > 1.0e-14);
      else
	failed = (fabs(eye(i,j)) > 1.0e-14);
  cerr << "eye = " << eye << endl;
  if (failed)
    cerr << "Failed inverse test.\n";
  else
    cerr << "Passed inverse test.\n";

  return !failed;
}
  

double MulTest(int N)
{
  Array<double,2> A(N,N), B(N,N), C(N,N);
  for (int i=0; i<N; i++) 
    for (int j=0; j<N; j++) {
      A(i,j) = 2.0*drand48()-1.0;
      B(i,j) = 2.0*drand48()-1.0;
    }

  clock_t tstart, tend;
  tstart = clock();
  for (int i=0; i<10000; i++)
    //    C = A*B;
    MatMult (A, B, C);
  tend = clock();

  double ops = 2*N*N*N + 2*N*N;
  double time = (double)(tend-tstart)/1.0e6;
  double rate = 10000.0*ops/time;
  return rate;
}
  

void TestDetCofactors()
{
  const int N = 5;
  int worksize = DetCofactorsWorksize(N);
  Array<double,1> workspace(worksize);

  Array<double,2> A(N,N), B(N,N);
  
  for (int i=0; i<N; i++)
    for (int j=0; j<N; j++)
      A(i,j) = drand48()-0.5;
  
  B = A;
  double detA = Determinant(A);
  GJInverse(A);
  Transpose(A);
  A = detA * A;
  double detB = DetCofactors(B, workspace);

  cerr << "det A = " << detA << endl;
  cerr << "det B = " << detB << endl;

  cerr << "A cofactors = " << A << endl;
  cerr << "B cofactors = " << B << endl;

}


void 
TimeDetCofactors()
{
  const int N = 16;
  int worksize = DetCofactorsWorksize(N);
  Array<double,1> workspace(worksize);

  Array<double,2> A(N,N), B(N,N);
  
  for (int i=0; i<N; i++)
    for (int j=0; j<N; j++)
      A(i,j) = drand48()-0.5;
  
  clock_t start, end;

  start = clock();
  const int num = 100000;

  for (int i=0; i < num; i++) {
    B = A;
    double detA = Determinant(B);
    GJInverse(B);
    Transpose(B);
    B = detA * B;
  }
  end = clock();
  fprintf (stderr, "Old way time = %1.5e\n", 
	   (double)(end-start)/((double)CLOCKS_PER_SEC * (double)num));
  
  start = clock();
  for (int i=0; i<num; i++) {
    B = A;
    double detB = DetCofactors(B, workspace);
  }
  end = clock();
  fprintf (stderr, "New way time = %1.5e\n", 
	   (double)(end-start)/((double)CLOCKS_PER_SEC * (double)num));
  
}




main()
{
  TestInvPartial();
  //TestDetCofactors();
  //TimeDetCofactors();
  //TimeDet();
//   for (int N=10; N<300; N+=5)
//     cerr << "N = " << N << " Rate = " << MulTest(N) << endl;

//   for (int N=10; N<300; N+=5)
//     cerr << "N = " << N << " Rate = " << MulTest(N) << endl;
//   TestDet();
//   TestInv();

//   Array<double,2> A(3,3), LU(3,3), Ainv(3,3);
//   A = 1.2,  -5.6,  7.3,
//      -2.1,  -9.5,  7.1,
//       4.9,   1.6, -6.6;
//   Array<double,1> b(3), x(3);
  
//   b = 1.0, 2.0, 3.0;

//   LU = A;
//   Array<int,1> perm;
//   double sign;
//   //cerr << "A = " << A << endl;
//   LUdecomp (LU, perm, sign);
//   x = b;
//   LUsolve (LU, perm, x);
//   //cerr << "LU = " << LU << endl;
//   //cerr << "b = " << b << endl;
//   //cerr << "x = " << x << endl;
//   Ainv = Inverse(A);
//   //cerr << "Ainv = " << Ainv << endl;
  
//   Array<double,2> AinvA(3,3);
//   AinvA = 0.0;
//   for (int i=0; i<3; i++)
//     for (int j=0; j<3; j++)
//       for (int k=0; k<3; k++)
// 	AinvA(i,j) += Ainv(i,k)*A(k,j);
//   bool success = true;
//   for (int i=0; i<3; i++)
//     for (int j=0; j<3; j++)
//       if ((i==j)) {
// 	if (fabs(AinvA(i,j) - 1.0) > 1e-14)
// 	  success = false;
//       }
//       else {
// 	if (fabs(AinvA(i,j)) > 1e-14)
// 	  success = false;
//       }
//   cerr << (success ? "Passed Matrix Inverse Test\n" : 
//     "Failed Matrix Inverse Test\n");

//   // SVD test
//   Array<double,2> U(3,3), V(3,3);
//   Array<double,1> S(3);
//   Array<double,2> B(4,3);

//   B = 1.2,  -5.6,  7.3,
//      -2.1,  -9.5,  7.1,
//       4.9,   1.6, -6.6,
//       1.3,   1.8,  9.0;

//   SVdecomp (B, U, S, V);
//   cerr << "B = " << B << endl;
//   cerr << "U = " << U << endl;
//   cerr << "V = " << V << endl;
//   cerr << "S = " << S << endl;


//   // Determinant test
//   double det = Determinant(A);
//   fprintf (stderr, "det(A) = %1.16e\n", det);

}
