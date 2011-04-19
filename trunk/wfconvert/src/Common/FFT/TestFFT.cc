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

#include "FFT.h"

void TestFFT3D()
{
  FFT3D fft;
  fft.resize(2,2,2);
  fft.kBox(0,0,0) = 1.0;  fft.kBox(0,0,1) = 2.0;
  fft.kBox(0,1,0) = 3.0;  fft.kBox(0,1,1) = 4.0;
  fft.kBox(1,0,0) = 5.0;  fft.kBox(1,0,1) = 6.0;
  fft.kBox(1,1,0) = 7.0;  fft.kBox(1,1,1) = 8.0;

  cerr << "Orignal:\n";
  for (int i=0; i<2; i++)
    for (int j=0; j<2; j++)
      for (int k=0; k<2; k++)
	fprintf (stderr, "%1.5f\n", fft.kBox(i,j,k).real());
  fft.k2r();
  //fft.kBox = 0.0;
  fft.r2k();
    cerr << "FFT^-1(FFT(Original)):\n";
  for (int i=0; i<2; i++)
    for (int j=0; j<2; j++)
      for (int k=0; k<2; k++)
	fprintf (stderr, "%1.5f\n", fft.kBox(i,j,k).real());
}

void TestFFTVec3D()
{
  FFTVec3D fft;
  fft.resize(2,2,2);
  fft.kBox(0,0,0)[0]=1.0;  fft.kBox(0,0,0)[1]=1.1;  fft.kBox(0,0,0)[2]=1.2;  
  fft.kBox(0,0,1)[0]=2.0;  fft.kBox(0,0,1)[1]=2.1;  fft.kBox(0,0,1)[2]=2.2;  
  fft.kBox(0,1,0)[0]=3.0;  fft.kBox(0,1,0)[1]=3.1;  fft.kBox(0,1,0)[2]=3.2;  
  fft.kBox(0,1,1)[0]=4.0;  fft.kBox(0,1,1)[1]=4.1;  fft.kBox(0,1,1)[2]=4.2;  
  fft.kBox(1,0,0)[0]=5.0;  fft.kBox(1,0,0)[1]=5.1;  fft.kBox(1,0,0)[2]=5.2;  
  fft.kBox(1,0,1)[0]=6.0;  fft.kBox(1,0,1)[1]=6.1;  fft.kBox(1,0,1)[2]=6.2;  
  fft.kBox(1,1,0)[0]=7.0;  fft.kBox(1,1,0)[1]=7.1;  fft.kBox(1,1,0)[2]=7.2;  
  fft.kBox(1,1,1)[0]=8.0;  fft.kBox(1,1,1)[1]=8.1;  fft.kBox(1,1,1)[2]=8.2;  

  cerr << "Orignal:\n";
  for (int i=0; i<2; i++)
    for (int j=0; j<2; j++)
      for (int k=0; k<2; k++){
	for (int l=0; l<3; l++)
	  fprintf (stderr, "%1.1f ", fft.kBox(i,j,k)[l].real());
	fprintf (stderr, "\n");
      }
  fft.k2r();
  //fft.kBox = complex<double>(0.0,0.0);
  fft.r2k();
  cerr << "FFT^-1(FFT(Original)):\n";
  for (int i=0; i<2; i++)
    for (int j=0; j<2; j++)
      for (int k=0; k<2; k++) {
	for (int l=0; l<3; l++)
	  fprintf (stderr, "%1.1f ", fft.kBox(i,j,k)[l].real());
	fprintf (stderr, "\n");
      }
}

void TestFFTMat3D()
{
  FFTMat3D fft;
  fft.resize(2,2,2);
  fft.kBox(0,0,0)(0,0)=1.0; fft.kBox(0,0,0)(0,1)=1.1; fft.kBox(0,0,0)(0,2)=1.2;
  fft.kBox(0,0,0)(1,0)=1.3; fft.kBox(0,0,0)(1,1)=1.4; fft.kBox(0,0,0)(1,2)=1.5;
  fft.kBox(0,0,0)(2,0)=1.6; fft.kBox(0,0,0)(2,1)=1.7; fft.kBox(0,0,0)(2,2)=1.8;

  fft.kBox(0,0,1)(0,0)=2.0; fft.kBox(0,0,1)(0,1)=2.1; fft.kBox(0,0,1)(0,2)=2.2;
  fft.kBox(0,0,1)(1,0)=2.3; fft.kBox(0,0,1)(1,1)=2.4; fft.kBox(0,0,1)(1,2)=2.5;
  fft.kBox(0,0,1)(2,0)=2.6; fft.kBox(0,0,1)(2,1)=2.7; fft.kBox(0,0,1)(2,2)=2.8;


  fft.kBox(0,1,0)(0,0)=3.0; fft.kBox(0,1,0)(0,1)=3.1; fft.kBox(0,1,0)(0,2)=3.2;
  fft.kBox(0,1,0)(1,0)=3.3; fft.kBox(0,1,0)(1,1)=3.4; fft.kBox(0,1,0)(1,2)=3.5;
  fft.kBox(0,1,0)(2,0)=3.6; fft.kBox(0,1,0)(2,1)=3.7; fft.kBox(0,1,0)(2,2)=3.8;

  fft.kBox(0,1,1)(0,0)=4.0; fft.kBox(0,1,1)(0,1)=4.1; fft.kBox(0,1,1)(0,2)=4.2;
  fft.kBox(0,1,1)(1,0)=4.3; fft.kBox(0,1,1)(1,1)=4.4; fft.kBox(0,1,1)(1,2)=4.5;
  fft.kBox(0,1,1)(2,0)=4.6; fft.kBox(0,1,1)(2,1)=4.7; fft.kBox(0,1,1)(2,2)=4.8;



  fft.kBox(1,0,0)(0,0)=5.0; fft.kBox(1,0,0)(0,1)=5.1; fft.kBox(1,0,0)(0,2)=5.2;
  fft.kBox(1,0,0)(1,0)=5.3; fft.kBox(1,0,0)(1,1)=5.4; fft.kBox(1,0,0)(1,2)=5.5;
  fft.kBox(1,0,0)(2,0)=5.6; fft.kBox(1,0,0)(2,1)=5.7; fft.kBox(1,0,0)(2,2)=5.8;

  fft.kBox(1,0,1)(0,0)=6.0; fft.kBox(1,0,1)(0,1)=6.1; fft.kBox(1,0,1)(0,2)=6.2;
  fft.kBox(1,0,1)(1,0)=6.3; fft.kBox(1,0,1)(1,1)=6.4; fft.kBox(1,0,1)(1,2)=6.5;
  fft.kBox(1,0,1)(2,0)=6.6; fft.kBox(1,0,1)(2,1)=6.7; fft.kBox(1,0,1)(2,2)=6.8;


  fft.kBox(1,1,0)(0,0)=7.0; fft.kBox(1,1,0)(0,1)=7.1; fft.kBox(1,1,0)(0,2)=7.2;
  fft.kBox(1,1,0)(1,0)=7.3; fft.kBox(1,1,0)(1,1)=7.4; fft.kBox(1,1,0)(1,2)=7.5;
  fft.kBox(1,1,0)(2,0)=7.6; fft.kBox(1,1,0)(2,1)=7.7; fft.kBox(1,1,0)(2,2)=7.8;

  fft.kBox(1,1,1)(0,0)=8.0; fft.kBox(1,1,1)(0,1)=8.1; fft.kBox(1,1,1)(0,2)=8.2;
  fft.kBox(1,1,1)(1,0)=8.3; fft.kBox(1,1,1)(1,1)=8.4; fft.kBox(1,1,1)(1,2)=8.5;
  fft.kBox(1,1,1)(2,0)=8.6; fft.kBox(1,1,1)(2,1)=8.7; fft.kBox(1,1,1)(2,2)=8.8;

  cerr << "Orignal:\n";
  for (int i=0; i<2; i++)
    for (int j=0; j<2; j++)
      for (int k=0; k<2; k++){
	for (int l=0; l<3; l++)
	  for (int m=0; m<3; m++)
	    fprintf (stderr, "%1.1f ", fft.kBox(i,j,k)(l,m).real());
	fprintf (stderr, "\n");
      }
  fft.k2r();
//   complex<double> z(0.0, 0.0);
//   cMat3 a;
//   a(0,0)=z; a(0,1)=z; a(0,2)=z;
//   a(1,0)=z; a(1,1)=z; a(1,2)=z;
//   a(2,0)=z; a(2,1)=z; a(2,2)=z;
//    fft.kBox = a;

  fft.r2k();
  cerr << "FFT^-1(FFT(Original)):\n";
  for (int i=0; i<2; i++)
    for (int j=0; j<2; j++)
      for (int k=0; k<2; k++) {
	for (int l=0; l<3; l++)
	  for (int m=0; m<3; m++)
	    fprintf (stderr, "%1.1f ", fft.kBox(i,j,k)(l,m).real());
	fprintf (stderr, "\n");
      }
}


#include <time.h>
#include <stdlib.h>

int N = 32;

void TimeTest()
{
  FFT3D fft;
  fft.resize(N,N,N);
  for (int i=0; i<N; i++)
    for (int j=0; j<N; j++)
      for (int k=0; k<N; k++)
	fft.rBox(i,j,k) = drand48();

  clock_t start,end;
  start = clock();
  for (int n=0; n<10000; n++)
    fft.r2k();
  end = clock();

  fprintf (stderr, "Time for (%d,%d,%d) complex FFT = %1.5f\n",
	   N,N,N, (double)(end-start)/((double)10000.0*CLOCKS_PER_SEC));

}

void TimeTestVec()
{
  FFTVec3D fft;
  fft.resize(N,N,N);
  for (int i=0; i<N; i++)
    for (int j=0; j<N; j++)
      for (int k=0; k<N; k++)
	for (int l=0; l<3; l++)
	  fft.rBox(i,j,k)[l] = drand48();

  clock_t start,end;
  start = clock();
  for (int n=0; n<10000; n++)
    fft.r2k();
  end = clock();

  fprintf (stderr, "Time for (%d,%d,%d)[3] Vec3 complex FFT = %1.5f\n",
	   N,N,N, (double)(end-start)/((double)10000.0*CLOCKS_PER_SEC));

}

void TimeTestMat()
{
  FFTMat3D fft;
  fft.resize(N,N,N);
  for (int i=0; i<N; i++)
    for (int j=0; j<N; j++)
      for (int k=0; k<N; k++)
	for (int l=0; l<3; l++)
	  for (int m=0; m<3; m++)
	    fft.rBox(i,j,k)(l,m) = drand48();

  clock_t start,end;
  start = clock();
  for (int n=0; n<1000; n++)
    fft.r2k();
  end = clock();

  fprintf (stderr, "Time for (%d,%d,%d)[3] Mat3 complex FFT = %1.5f\n",
	   N,N,N, (double)(end-start)/((double)1000.0*CLOCKS_PER_SEC));

}



main()
{
  TestFFT3D();
  TestFFTVec3D();
  TestFFTMat3D();
  TimeTest();
  TimeTestVec();
  TimeTestMat();
}
