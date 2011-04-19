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

#include "Hamiltonian.h"
#include "../PH/kSpacePH.h"

void KineticClass::Setup()
{
  halfG2.resize(GVecs.size());
  for (int i=0; i<GVecs.size(); i++)
    halfG2(i) = 0.5*dot(GVecs(i),GVecs(i));
  IsSetup = true;
}

void KineticClass::Apply(const zVec &c, zVec &Kc)
{
  if (!IsSetup)
    Setup();
  for (int i=0; i<c.size(); i++)
    Kc(i) += halfG2(i)*c(i);
}


void CoulombClass::Apply(const zVec &c, zVec &Vc)
{
  double volInv = 1.0/GVecs.GetBoxVol();
  for (int i=0; i<c.size(); i++) {
    for (int j=0; j<i; j++) {
      Vec3 diff = GVecs(i) - GVecs(j);
      Vc(i) += -4.0*volInv*M_PI*Z/dot(diff,diff)*c(j);
    }
    for (int j=i+1; j<c.size(); j++) {
      Vec3 diff = GVecs(i) - GVecs(j);
      Vc(i) += -4.0*volInv*M_PI*Z/dot(diff,diff)*c(j);
    }
  }
}


void CoulombFFTClass::Setup()
{
  FFT.Setup();
  Vec3 dG;
  Int3 dI;
  //double prefact = -4.0*M_PI*Z;// /GVecs.GetBoxVol();
  double prefact = -4.0*M_PI*Z/GVecs.GetBoxVol();
  FFT.kBox = 0.0;
//   for (int m=0; m<GVecs.size(); m++)
//     for (int n=0; n<GVecs.size(); n++) {
//       dG = GVecs(m)-GVecs(n);
//       dI = GVecs.Index(n) - GVecs.Index(m);
//       double val = (m == n) ? 0.0 : prefact/dot(dG,dG);
//       int i = (dI[0] + FFT.Nx) % FFT.Nx;
//       int j = (dI[1] + FFT.Ny) % FFT.Ny;
//       int k = (dI[2] + FFT.Nz) % FFT.Nz;
//       FFT.kBox(i,j,k) = val;
//     }
  Vec3 Gdiff;
  int Nx, Ny, Nz;
  FFT.GetDims(Nx, Ny, Nz);
  cerr << "Nx=" << Nx << " Ny=" << Ny << " Nz=" << Nz << endl;
  int mx = (Nx-1)/2;
  int my = (Ny-1)/2;
  int mz = (Nz-1)/2;
  double cut = 4.0*GVecs.GetkCut() * GVecs.GetkCut();
  for (int ix=-mx; ix<=mx; ix++) {
    int nx = (ix+Nx)%Nx;
    Gdiff[0] = GVecs.GetkBox()[0]*ix;
    for (int iy=-my; iy<=my; iy++) {
      int ny = (iy+Ny)%Ny;
      Gdiff[1] = GVecs.GetkBox()[1]*iy;
      for (int iz=-mz; iz<=mz; iz++) {
	int nz = (iz+Nz)%Nz;
	Gdiff[2] = GVecs.GetkBox()[2]*iz;
	double prod = dot(Gdiff,Gdiff);
	if (prod < cut) {
	  //cerr << "oldval = " << FFT.kBox(nx,ny,nz) <<endl;
	  FFT.kBox(nx,ny,nz) = prefact/prod;
	  //cerr << "newval = " << FFT.kBox(nx,ny,nz) <<endl;
	}
      }
    }
  }
  FFT.kBox(0,0,0) = 0.0;
    
  FFT.k2r();
  FILE *fout = fopen ("Vr.dat", "w");
  for (int i=0; i<Nx; i++) {
    for (int j=0; j<Ny; j++)
      fprintf (fout, "%1.15e ", FFT.rBox(i,j,0).real());
    fprintf (fout, "\n");
  }
  fclose (fout);
  Vr.resize(Nx, Ny, Nz);
  Vr = FFT.rBox;
  IsSetup = true;
  cerr << "Finished CoulombFFT setup.\n";
}

void 
CoulombFFTClass::Apply(const zVec &c, zVec &Vc)
{
  if (!IsSetup)
    Setup();

  int Nx, Ny, Nz;
  FFT.GetDims(Nx, Ny, Nz);
  // First, put c into FFTbox
  FFT.PutkVec (c);
  // Now, transform to real space
  FFT.k2r();
  // Now, multiply by V
  FFT.rBox *= Vr;
//   for (int i=0; i<FFT.Nx; i++)
//     for (int j=0; j<FFT.Ny; j++)
//       for (int k=0; k<FFT.Nz; k++)
// 	FFT.rBox(i,j,k) *= Vr(i,j,k); 

  // Transform back
  FFT.r2k();
  // And put into Vc
  FFT.AddToVec (Vc);
}


void PHPotClass::Setup()
{
  kPH.CalcTailCoefs (30.0, 80.0);
  VGGp.resize (GVecs.size(), GVecs.size());
  Vec3 k (0.0, 0.0, 0.0);
  double volInv = 1.0/GVecs.GetBoxVol();
  for (int i=0; i<GVecs.size(); i++)
    for (int j=i; j<GVecs.size(); j++) {
      VGGp (i,j) = kPH.V(k, GVecs(i), GVecs(j))*volInv;
      VGGp (j,i) = conj (VGGp(i,j));
    }

  IsSetup = true;
}

void PHPotClass::Apply (const zVec &c, zVec &Hc)
{
  if (!IsSetup)
    Setup();
      
  for (int i=0; i<c.size(); i++)
    for (int j=0; j<c.size(); j++)
      Hc(i) += VGGp (i,j) * c(j);
}


void Hamiltonian::Apply(const zVec &c, zVec &Hc)
{
  Hc = 0.0;
  Kinetic.Apply (c, Hc);
  //Coulomb.Apply (c, Hc);
  CoulombFFT.Apply (c, Hc);
  //PH.Apply (c, Hc);
  //PHFFT.Apply (c, Hc);
}


void PHPotFFTClass::Setup()
{
  cerr << "In PHPotFFTClass::Setup().\n";
  kPH.CalcTailCoefs (15.0, 60.0);
  int nx, ny, nz;
  cerr << "Before cFFT Setup.\n";
  cFFT.Setup();
  cerr << "Before VecFFT Setup.\n";
  VecFFT.Setup();
  cerr << "Before MatFFT Setup.\n";
  FFTMatBox MatFFT(GVecs);
  MatFFT.Setup();
  MatFFT.GetDims(nx,ny,nz);
  Fr.resize(nx,ny,nz);
  Vr.resize(nx,ny,nz);
  Gc.resize(GVecs.size());
  Vc.resize(GVecs.size());
  Vk.resize(GVecs.DeltaSize());
  Fk.resize(GVecs.DeltaSize());
  StructureFactor.resize(GVecs.DeltaSize());
  StructureFactor = complex<double>(1.0, 0.0);

  double volInv = 1.0/GVecs.GetBoxVol();

  // First, setup V and F tensors in k-space
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
    Vk(i) = volInv*StructureFactor(i)*V;
    Fk(i) = volInv *  StructureFactor(i) * 
      kPH.Ftensor (GVecs.DeltaG(i), a, bPerp, bPar);
  }
  cerr << "Number of GetVals calls = " << numCalls << endl;

//   /// FFT to calculate V in real space
//   for (int i=0; i<GVecs.DeltaSize(); i++) {
//     Vec3 dG = GVecs.DeltaG(i);
//     double dGmag = sqrt (dot(dG,dG));
//     Vk(i) = StructureFactor(i)*kPH.V(dGmag)*volInv;
//   }
  cFFT.PutkVec (Vk);
  cFFT.k2r();
  for (int ix=0; ix<nx; ix++)
    for (int iy=0; iy<ny; iy++)
      for (int iz=0; iz<nz; iz++)
	Vr(ix,iy,iz) = cFFT.rBox(ix,iy,iz);

  /// FFT to calculate F tensor in real space
//   for (int i=0; i<GVecs.DeltaSize(); i++) 
//     Fk(i) = volInv*StructureFactor(i)*kPH.Ftensor(GVecs.DeltaG(i));
  MatFFT.PutkVec(Fk);
  MatFFT.k2r();
  for (int ix=0; ix<nx; ix++)
    for (int iy=0; iy<ny; iy++)
      for (int iz=0; iz<nz; iz++)
	Fr(ix,iy,iz) = MatFFT.rBox(ix,iy,iz);
  IsSetup = true;
}
  
void PHPotFFTClass::Apply (const zVec &c, zVec &Hc)
{
  if (!IsSetup)
    Setup();
  int nx, ny, nz;
  cFFT.GetDims(nx, ny, nz);

  /////////////////////////
  // Pseudo-kinetic part //
  /////////////////////////
  // First, build G*c vector in VecFFT  
  for (int i=0; i<GVecs.size(); i++) 
    Gc(i) = c(i)*(k+GVecs(i));
  VecFFT.PutkVec (Gc);
  VecFFT.k2r();
  
  // Now, muliply by F(r)
  for (int ix=0; ix<nx; ix++)
    for (int iy=0; iy<ny; iy++)
      for (int iz=0; iz<nz; iz++)
	VecFFT.rBox(ix,iy,iz) = Fr(ix,iy,iz)*VecFFT.rBox(ix,iy,iz);
  // Transform back;
  VecFFT.r2k();
  // Rename for clarity
  zVecVec& FGc = Gc;
  VecFFT.GetkVec (FGc);
  for (int i=0; i<GVecs.size(); i++)
    Hc(i) += 0.5*dot(GVecs(i), FGc(i));

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
