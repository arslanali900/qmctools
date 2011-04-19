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

void 
FFT1D::resize(int n)
{
  if (Allocated) {
    FFT_FREE(rData);
    if (!InPlace)
      FFT_FREE(kData);
    FFTNAME(destroy_plan)(r2kPlan);
    FFTNAME(destroy_plan)(k2rPlan);
  }
  rData = (complex<FFT_FLOAT>*) FFT_MALLOC(sizeof(FFTNAME(complex))*n+FFT_EXTRA_MEM);
  if (!InPlace)
    kData = (complex<FFT_FLOAT>*) FFT_MALLOC(sizeof(FFTNAME(complex))*n+FFT_EXTRA_MEM);
  else
    kData = rData;

  complex<FFT_FLOAT>* rAligned = FFTAlign(rData);
  complex<FFT_FLOAT>* kAligned = FFTAlign(kData);

  Array<complex<FFT_FLOAT>,1> *temp;
  temp = new Array<complex<FFT_FLOAT>,1>(rAligned, shape(n), neverDeleteData);
  rBox.reference (*temp);
  delete temp;
  
  temp = new Array<complex<FFT_FLOAT>,1>(kAligned, shape(n), neverDeleteData);
  kBox.reference(*temp);
  delete temp;

  r2kPlan = FFTNAME(plan_dft_1d)
    (n, reinterpret_cast<FFTNAME(complex)*>(rAligned), 
     reinterpret_cast<FFTNAME(complex)*>(kAligned), 1, FFTW_MEASURE);
  k2rPlan = FFTNAME(plan_dft_1d) 
    (n, reinterpret_cast<FFTNAME(complex)*>(kAligned), 
     reinterpret_cast<FFTNAME(complex)*>(rAligned), -1,FFTW_MEASURE);

  Ninv = 1.0/(FFT_FLOAT)n;
  Allocated = true;
}

void 
FFT1D::r2k()
{
  FFTNAME(execute)(r2kPlan);
  rBox *= Ninv;
}


void 
FFT1D::k2r()
{
  FFTNAME(execute)(k2rPlan);
}




void 
FFT3D::resize(int nx, int ny, int nz)
{
  if (Allocated) {
    FFT_FREE(rData);
    if (!InPlace)
      FFT_FREE(kData);
    FFTNAME(destroy_plan)(r2kPlan);
    FFTNAME(destroy_plan)(k2rPlan);
  }
  rData = (complex<FFT_FLOAT>*) FFT_MALLOC(sizeof(FFTNAME(complex))*nx*ny*nz+FFT_EXTRA_MEM);
  if (!InPlace)
    kData = (complex<FFT_FLOAT>*) FFT_MALLOC(sizeof(FFTNAME(complex))*nx*ny*nz+FFT_EXTRA_MEM);
  else
    kData = rData;
  
  complex<FFT_FLOAT>* rAligned = FFTAlign(rData);
  complex<FFT_FLOAT>* kAligned = FFTAlign(kData);
  


  Array<complex<FFT_FLOAT>,3> *temp;
  temp = new Array<complex<FFT_FLOAT>,3>(rAligned, shape(nx,ny,nz), neverDeleteData);
  rBox.reference (*temp);
  delete temp;

  temp = new Array<complex<FFT_FLOAT>,3>(kAligned, shape(nx,ny,nz), neverDeleteData);
  kBox.reference(*temp);
  delete temp;

  r2kPlan = FFTNAME(plan_dft_3d)
    (nx, ny, nz, reinterpret_cast<FFTNAME(complex)*>(rAligned), 
     reinterpret_cast<FFTNAME(complex)*>(kAligned), 1, FFTW_MEASURE);
  k2rPlan =  FFTNAME(plan_dft_3d)
    (nx, ny, nz, reinterpret_cast<FFTNAME(complex)*>(kAligned), 
     reinterpret_cast<FFTNAME(complex)*>(rAligned), -1,FFTW_MEASURE);
  
  Ninv = 1.0/(FFT_FLOAT)(nx*ny*nz);
  Allocated = true;
}


void 
FFT3D::r2k()
{
  FFTNAME(execute)(r2kPlan);
  kBox *= Ninv;
//   for (int i=0; i<kBox.size(); i++)
//     *(kBox.data()+i) *= sqrtNinv;
}


void 
FFT3D::k2r()
{
  FFTNAME(execute)(k2rPlan);
//   for (int i=0; i<rBox.size(); i++)
//     *(rBox.data()+i) *= sqrtNinv;
}


void 
FFT3D_r2c::resize(int nx, int ny, int nz)
{
  if (Allocated) {
    FFT_FREE(rData);
    FFT_FREE(kData);
    FFTNAME(destroy_plan)(r2kPlan);
    FFTNAME(destroy_plan)(k2rPlan);
  }
  int nz_k = nz/2+1;

  rData = (FFT_FLOAT*) FFT_MALLOC(sizeof(FFT_FLOAT)*nx*ny*nz+FFT_EXTRA_MEM);
  kData = (complex<FFT_FLOAT>*) FFT_MALLOC(sizeof(FFTNAME(complex))*nx*ny*nz_k+FFT_EXTRA_MEM);
  
  FFT_FLOAT*          rAligned = FFTAlign(rData);
  complex<FFT_FLOAT>* kAligned = FFTAlign(kData);
  


  Array<FFT_FLOAT,3> *temp;
  temp = new Array<FFT_FLOAT,3>(rAligned, shape(nx,ny,nz), neverDeleteData);
  rBox.reference (*temp);
  delete temp;

  Array<complex<FFT_FLOAT>,3> *ktemp;
  ktemp = new Array<complex<FFT_FLOAT>,3>(kAligned, shape(nx,ny,nz_k), neverDeleteData);
  kBox.reference(*ktemp);
  delete ktemp;

  r2kPlan = FFTNAME(plan_dft_r2c_3d)
    (nx, ny, nz, (FFT_FLOAT*)(rAligned), 
     reinterpret_cast<FFTNAME(complex)*>(kAligned), FFTW_MEASURE);
  k2rPlan =  FFTNAME(plan_dft_c2r_3d)
    (nx, ny, nz, reinterpret_cast<FFTNAME(complex)*>(kAligned), 
     (FFT_FLOAT*)(rAligned),FFTW_MEASURE);
  
  Ninv = 1.0/(FFT_FLOAT)(nx*ny*nz);
  Allocated = true;
}


void 
FFT3D_r2c::r2k()
{
  FFTNAME(execute)(r2kPlan);
  kBox *= Ninv;
}


void 
FFT3D_r2c::k2r()
{
  FFTNAME(execute)(k2rPlan);
}


FFT3D&
FFT3D::operator= (const FFT3D& fft)
{
  InPlace = fft.InPlace;
  if (rBox.shape() != fft.rBox.shape()) 
    resize (fft.rBox.extent(0), fft.rBox.extent(1), fft.rBox.extent(2));
  rBox = fft.rBox;
  kBox = fft.kBox;
  return *this;
}


void 
FFTVec3D::resize(int nx, int ny, int nz)
{
  if (Allocated) {
    FFT_FREE(rData);
    if (!InPlace)
      FFT_FREE(kData);
    FFTNAME(destroy_plan)(r2kPlan);
    FFTNAME(destroy_plan)(k2rPlan);
  }
  rData = ( TinyVector<complex<FFT_FLOAT>,3>*) FFT_MALLOC(3*sizeof(FFTNAME(complex))*nx*ny*nz+FFT_EXTRA_MEM);
  if (!InPlace)
    kData = ( TinyVector<complex<FFT_FLOAT>,3>*) FFT_MALLOC(3*sizeof(FFTNAME(complex))*nx*ny*nz+FFT_EXTRA_MEM);
  else
    kData = rData;

  TinyVector<complex<FFT_FLOAT>,3>* rAligned = 
    (TinyVector<complex<FFT_FLOAT>,3>*) FFTAlign((complex<FFT_FLOAT>*) rData);
  TinyVector<complex<FFT_FLOAT>,3>* kAligned = 
    (TinyVector<complex<FFT_FLOAT>,3>*) FFTAlign((complex<FFT_FLOAT>*) kData);

  Array< TinyVector<complex<FFT_FLOAT>,3>,3> *temp;
  temp = new Array< TinyVector<complex<FFT_FLOAT>,3>,3>(rAligned, shape(nx,ny,nz), neverDeleteData);
  rBox.reference (*temp);
  delete temp;

  temp = new Array< TinyVector<complex<FFT_FLOAT>,3>,3>(kAligned, shape(nx,ny,nz), neverDeleteData);
  kBox.reference(*temp);
  delete temp;

  int n[3] = {nx, ny, nz};
  r2kPlan = FFTNAME(plan_many_dft)(3, n, 3, 
				   reinterpret_cast<FFTNAME(complex)*>(rAligned),
				   n, 3, 1, 
				   reinterpret_cast<FFTNAME(complex)*>(kAligned), 
				   n, 3, 1, 1, FFTW_MEASURE);
  assert (r2kPlan != NULL);
  k2rPlan = 
    FFTNAME(plan_many_dft)(3, n, 3, reinterpret_cast<FFTNAME(complex)*>(kAligned), 
			   n, 3, 1, 
			   reinterpret_cast<FFTNAME(complex)*>(rAligned), n,
			   3, 1, -1, FFTW_MEASURE);
  assert (k2rPlan != NULL);
  Ninv = 1.0/(FFT_FLOAT)(nx*ny*nz);
  Allocated = true;
}


void 
FFTVec3D::r2k()
{
  FFTNAME(execute)(r2kPlan);
  kBox *= Ninv;
}


void 
FFTVec3D::k2r()
{
  FFTNAME(execute)(k2rPlan);
}




void 
FFTMat3D::resize(int nx, int ny, int nz)
{
  if (Allocated) {
    FFT_FREE(rData);
    if (!InPlace)
      FFT_FREE(kData);
    FFTNAME(destroy_plan)(r2kPlan);
    FFTNAME(destroy_plan)(k2rPlan);
  }
  rData = (TinyMatrix<complex<FFT_FLOAT>,3,3>*) 
    FFT_MALLOC(9*sizeof(FFTNAME(complex))*nx*ny*nz+FFT_EXTRA_MEM);
  if (!InPlace)
    kData = (TinyMatrix<complex<FFT_FLOAT>,3,3>*) 
      FFT_MALLOC(9*sizeof(FFTNAME(complex))*nx*ny*nz+FFT_EXTRA_MEM);
  else
    kData = rData;
  
  TinyMatrix<complex<FFT_FLOAT>,3,3>* rAligned = 
    (TinyMatrix<complex<FFT_FLOAT>,3,3>*) FFTAlign((complex<FFT_FLOAT>*)rData);
  TinyMatrix<complex<FFT_FLOAT>,3,3>* kAligned = 
    (TinyMatrix<complex<FFT_FLOAT>,3,3>*) FFTAlign((complex<FFT_FLOAT>*)kData);

  Array<TinyMatrix<complex<FFT_FLOAT>,3,3>,3> *temp;
  temp = new Array<TinyMatrix<complex<FFT_FLOAT>,3,3>,3>
    (rAligned, shape(nx,ny,nz), neverDeleteData);
  rBox.reference (*temp);
  delete temp;
  
  temp = new Array<TinyMatrix<complex<FFT_FLOAT>,3,3>,3>
    (kAligned, shape(nx,ny,nz), neverDeleteData);
  kBox.reference(*temp);
  delete temp;

  int n[3] = {nx, ny, nz};
  r2kPlan = 
    FFTNAME(plan_many_dft)(3, n, 9, 
			   reinterpret_cast<FFTNAME(complex)*>(rAligned),
			   n, 9, 1, 
			   reinterpret_cast<FFTNAME(complex)*>(kAligned), 
			   n, 9, 1, 1, FFTW_MEASURE);
  assert (r2kPlan != NULL);
  k2rPlan = 
    FFTNAME(plan_many_dft)(3, n, 9, reinterpret_cast<FFTNAME(complex)*>(kAligned),
			   n, 9, 1, 
			   reinterpret_cast<FFTNAME(complex)*>(rAligned), 
			n, 9, 1, -1, FFTW_MEASURE);
  assert (k2rPlan != NULL);

  Ninv = 1.0/(FFT_FLOAT)(nx*ny*nz);
  Allocated = true;
}


void 
FFTMat3D::r2k()
{
  FFTNAME(execute)(r2kPlan);
  for (int ix=0; ix<kBox.extent(0); ix++)
    for (int iy=0; iy<kBox.extent(1); iy++)
      for (int iz=0; iz<kBox.extent(2); iz++) {
	kBox(ix,iy,iz)(0,0) *= Ninv;
	kBox(ix,iy,iz)(0,1) *= Ninv;
	kBox(ix,iy,iz)(0,2) *= Ninv;
	kBox(ix,iy,iz)(1,0) *= Ninv;
	kBox(ix,iy,iz)(1,1) *= Ninv;
	kBox(ix,iy,iz)(1,2) *= Ninv;
	kBox(ix,iy,iz)(2,0) *= Ninv;
	kBox(ix,iy,iz)(2,1) *= Ninv;
	kBox(ix,iy,iz)(2,2) *= Ninv;
      }
}


void 
FFTMat3D::k2r()
{
  FFTNAME(execute)(k2rPlan);
}

