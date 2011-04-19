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

#include "LocalPotFFT.h"

void
LocalPotFFTClass::Vmatrix (Array<complex<double>,2> &vmat)
{
  if (!IsSetup)
    Setup();
  double volInv = 1.0/GVecs.GetBoxVol();
  for (int i=0; i<vmat.rows(); i++) 
    for (int j=0; j<=i; j++) {
      Vec3 diff = GVecs(i) - GVecs(j);
      complex<double> s(0.0,0.0);
      for (int zi=0; zi<Rions.size(); zi++) {
	double cosVal, sinVal, phase;
	phase = dot (diff, Rions(zi));
	sincos(phase, &sinVal, &cosVal);
	s += complex<double> (cosVal,sinVal);
      }
      vmat (i,j) = s*kPH.V(kPoint, GVecs(i), GVecs(j))*volInv;
      vmat (j,i) = conj(vmat(i,j));
    }
//   for (int i=0; i<vmat.rows(); i++)
//     cerr << "vmat(" << i << "," << i << ") = " << vmat(i,i) <<
//     endl;
//  cerr << "vmat = " << vmat << endl;
}


void 
LocalPotFFTClass::SetupkPotentials()
{
  kPH.CalcTailCoefs (30.0, 60.0);

  double volInv = 1.0/GVecs.GetBoxVol();
  // Setup V and F tensors in k-space
  //  double gMag, lastMag2;
  double gMag, lastMag2;
  lastMag2 = -1.0;
  double a, bPerp, bPar, V;
  int numCalls = 0;
  for (int i=0; i<GVecs.DeltaSize(); i++) {
    double gMag2 = dot(GVecs.DeltaG(i), GVecs.DeltaG(i));
    if (fabs(gMag2-lastMag2) > 1.0e-12) {
      lastMag2 = gMag2;
      gMag  = sqrt(lastMag2);
      kPH.GetVals(gMag, a, bPerp, bPar, V);
      numCalls++;
    }
    if (gMag2 > 1.0e-10)
      VG(i) = volInv * V;
    else {
      VG(i) = 0.0;
      VG0 = volInv * V;
    }
  }
}

void 
LocalPotFFTClass::SetuprPotentials()
{
  cFFT.kBox = complex<FFT_FLOAT>(0.0, 0.0);
  for (int i=0; i<GVecs.DeltaSize(); i++) {
    Int3 I = GVecs.DeltaI(i);
    cFFT.kBox(I)   = StructureFactor(i)*VG(i);
  }
  cFFT.k2r();
  Vr = cFFT.rBox;
}

void
LocalPotFFTClass::SetIons(const Array<Vec3,1> &rions)
{
  // Calculate the structure factor
  VionBase::SetIons(rions);
  if (IsSetup)
    SetuprPotentials();
}


void 
LocalPotFFTClass::Setup()
{
  kPH.CalcTailCoefs (15.0, 60.0);
  int nx, ny, nz;
  cFFT.GetDims(nx,ny,nz);
  Vr.resize(nx,ny,nz);
  Vc.resize(GVecs.size());
  VG.resize(GVecs.DeltaSize());
//   StructureFactor.resize(GVecs.DeltaSize());
//   StructureFactor = complex<double>(1.0, 0.0);

  SetupkPotentials();
  SetuprPotentials();

  IsSetup = true;
}
  
void
LocalPotFFTClass::Setk (Vec3 k)
{
  kPoint = k;
  
  int nx, ny, nz;
  cFFT.GetDims(nx,ny,nz);
  Vr.resize(nx,ny,nz);
  Vc.resize(GVecs.size());
  VG.resize(GVecs.DeltaSize());

  SetupkPotentials();
  SetIons(Rions);
}


void 
LocalPotFFTClass::Apply (const zVec &c, zVec &Hc)
{
  if (!IsSetup)
    Setup();
  int nx, ny, nz;
  cFFT.GetDims(nx, ny, nz);

  ////////////////////
  // Potential part //
  ////////////////////
  // Transform c into real space
  cFFT.PutkVec (c);
  cFFT.k2r();
  // Multiply by V
  cFFT.rBox *= Vr;
  // Transform back
  cFFT.r2k();

  // Get vector
  cFFT.GetkVec (Vc);
  for (int i=0; i<GVecs.size(); i++)
    Hc(i) += Vc(i);
}

void 
LocalPotFFTClass::Apply (const zVec &c, zVec &Hc,
			 Array<double,3> &VHXC)
{
  if (!IsSetup)
    Setup();
  int nx, ny, nz;
  cFFT.GetDims(nx, ny, nz);

  ////////////////////
  // Potential part //
  ////////////////////
  // Transform c into real space
  cFFT.PutkVec (c);
  cFFT.k2r();
  // Multiply by V
  //  cFFT.rBox *= (Vr+VHXC);
  for (int ix=0; ix<nx; ix++)
    for (int iy=0; iy<ny; iy++)
      for (int iz=0; iz<nz; iz++)
	cFFT.rBox(ix,iy,iz) *= (Vr(ix,iy,iz)+VHXC(ix,iy,iz));

  // Transform back
  cFFT.r2k();

  // Get vector
  cFFT.GetkVec (Vc);
  for (int i=0; i<GVecs.size(); i++)
    Hc(i) += Vc(i);
}

