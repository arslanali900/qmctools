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

#include "Hamiltonian2.h"
#include "../PH/kSpacePH.h"


void 
PWKineticClass::Setup()
{
  halfG2.resize(GVecs.size());
  for (int i=0; i<GVecs.size(); i++) {
    Vec3 Gpk = GVecs(i)+kPoint;
    halfG2(i) = 0.5*dot(Gpk, Gpk);
  }
  IsSetup = true;
}

void 
PWKineticClass::Setk (Vec3 k)
{
  kPoint = k;
  Setup();
}

void 
PWKineticClass::Apply(const zVec &c, zVec &Kc)
{
  if (!IsSetup)
    Setup();
  for (int i=0; i<c.size(); i++)
    Kc(i) += halfG2(i)*c(i);
}


void
VionBase::SetIons(const Array<Vec3,1> &rions)
{
  Rions.resize(rions.size());
  Rions = rions;
  if (StructureFactor.size() != GVecs.DeltaSize())
    StructureFactor.resize(GVecs.DeltaSize());
  StructureFactor = 0.0;
  for (int gi=0; gi<GVecs.DeltaSize(); gi++) {
    double c, s, phase;
    for (int i=0; i<rions.size(); i++) {
      phase = dot (GVecs.DeltaG(gi), rions(i));
      sincos (phase, &s, &c);
      StructureFactor(gi) += complex<double>(c, s);
    }
  }
}



void
CoulombClass::Setup()
{

}


void 
CoulombClass::Vmatrix (Array<complex<double>,2> &vmat)
{
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
      complex<double> val = -4.0*volInv*s*M_PI*Z/dot(diff,diff);
      vmat(i,j) = val;
      vmat(j,i) = conj(val);
    }
  for (int i=0; i<vmat.rows(); i++)
    vmat(i,i) = 0.0;
}

void 
CoulombClass::Apply(const zVec &c, zVec &Vc)
{
  double volInv = 1.0/GVecs.GetBoxVol();
  for (int i=0; i<c.size(); i++) {
    for (int j=0; j<i; j++) {
      Vec3 diff = GVecs(i) - GVecs(j);
      complex<double> s(0.0,0.0);
      for (int zi=0; zi<Rions.size(); zi++) {
	double cosVal, sinVal, phase;
	phase = dot (diff, Rions(zi));
	sincos(phase, &sinVal, &cosVal);
	s += complex<double> (cosVal,sinVal);
      }
      complex<double> val = -4.0*volInv*s*M_PI*Z/dot(diff,diff);
      Vc(i) += val*c(j);
      Vc(j) += conj(val)*c(i);
    }
  }
}


void
CoulombClass::SetIons (const Array<Vec3,1> &rions)
{
  if (Rions.size() != rions.size())
    Rions.resize(rions.size());
  Rions = rions;
}

void 
CoulombFFTClass::SetVr()
{
  FFT.kBox = 0.0;
  Vec3 Gdiff;
  int nx, ny, nz;
  FFT.GetDims(nx, ny, nz);
  double prefact = -4.0*M_PI*Z/GVecs.GetBoxVol();
  for (int i=0; i<GVecs.DeltaSize(); i++) {
    Vec3 dG = GVecs.DeltaG(i);
    FFT.kBox(GVecs.DeltaI(i)) = StructureFactor(i)* prefact/dot(dG,dG);
  }
  FFT.kBox(0,0,0) = 0.0;
  FFT.k2r();
  Vr.resize(nx, ny, nz);
  Vr = FFT.rBox;
}

void 
CoulombFFTClass::Setup()
{
  SetVr();
  IsSetup = true;
}


void 
CoulombFFTClass::Vmatrix (Array<complex<double>,2> &vmat)
{
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
      complex<double> val = -4.0*volInv*s*M_PI*Z/dot(diff,diff);
      vmat(i,j) = val;
      vmat(j,i) = conj(val);
    }
  for (int i=0; i<vmat.rows(); i++)
    vmat(i,i) = 0.0;
}


void 
CoulombFFTClass::SetIons(const Array<Vec3,1> &rions)
{
  VionBase::SetIons(rions);
  if (!IsSetup)
    Setup();
  else
    SetVr();
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
  // Transform back
  FFT.r2k();
  // And put into Vc
  FFT.AddToVec (Vc);
}


void 
PHPotClass::CalcStructFact()
{
  StructFact.resize(GVecs.size(), GVecs.size());
  StructFact = 0.0;
  for (int i=0; i<GVecs.size(); i++)
    for (int j=0; j<GVecs.size(); j++) {
      Vec3 Gdiff = GVecs(i) - GVecs(j);
      for (int k=0; k<Rions.size(); k++) {
	double phase, c, s;
	phase = dot (Gdiff, Rions(k));
	sincos (phase, &s, &c);
	StructFact(i,j) += complex<double>(c, s);
      }
    }
  SFIsSet = true;
}


void 
PHPotClass::SetIons (const Array<Vec3,1> &rions)
{
  Rions.resize(rions.size());
  Rions = rions;
  CalcStructFact();
}

void
PHPotClass::Vmatrix (Array<complex<double>,2> &vmat)
{
  if (!IsSetup) {
    cerr << "Calling PHPotClass::Setup.\n";
    Setup();
  }
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
}

void 
PHPotClass::SetVmat()
{
  VGGp.resize (GVecs.size(), GVecs.size());
  double volInv = 1.0/GVecs.GetBoxVol();
  for (int i=0; i<GVecs.size(); i++)
    for (int j=i; j<GVecs.size(); j++) {
      VGGp (i,j) = kPH.V(kPoint, GVecs(i), GVecs(j))*volInv;
      VGGp (j,i) = conj (VGGp(i,j));
    } 
  VmatIsSet = true;
}

void 
PHPotClass::Setk(Vec3 k)
{
  kPoint = k;
  kPH.CalcTailCoefs (30.0, 80.0);
  SetVmat();
}


void 
PHPotClass::Setup()
{
  kPH.CalcTailCoefs (30.0, 80.0);
  if (!VmatIsSet)
    SetVmat();
  if (!SFIsSet)
    CalcStructFact();
  IsSetup = true;
}

void 
PHPotClass::Apply (const zVec &c, zVec &Hc)
{
  if (!IsSetup)
    Setup();
      
  for (int i=0; i<c.size(); i++)
    for (int j=0; j<c.size(); j++)
      Hc(i) += StructFact(i,j)*VGGp (i,j) * c(j);
}



void
PHPotFFTClass::Vmatrix (Array<complex<double>,2> &vmat)
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
PHPotFFTClass::SetupkPotentials()
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
    Vk(i) = volInv * V;
    Fk(i) = volInv * kPH.Ftensor (GVecs.DeltaG(i), a, bPerp, bPar);
  }
}

void 
PHPotFFTClass::SetuprPotentials()
{
  cFFT.kBox = complex<double>(0.0, 0.0);
  complex<double> z(0.0, 0.0);
  cMat3 Mat0;
  Mat0(0,0)=z; Mat0(0,1)=z; Mat0(0,2)=z; 
  Mat0(1,0)=z; Mat0(1,1)=z; Mat0(1,2)=z;
  Mat0(2,0)=z; Mat0(2,1)=z; Mat0(2,2)=z;
  MatFFT.kBox = Mat0;
  for (int i=0; i<GVecs.DeltaSize(); i++) {
    Int3 I = GVecs.DeltaI(i);
    cFFT.kBox(I)   = StructureFactor(i)*Vk(i);
    MatFFT.kBox(I) = StructureFactor(i)*Fk(i);
  }
  cFFT.k2r();
  Vr = cFFT.rBox;
  MatFFT.k2r();
  Fr = MatFFT.rBox;
}

void
PHPotFFTClass::SetIons(const Array<Vec3,1> &rions)
{
  // Calculate the structure factor
  VionBase::SetIons(rions);
  if (IsSetup)
    SetuprPotentials();
}


void 
PHPotFFTClass::Setup()
{
  cerr << "Setting up FFT boxes:\n";
  VecFFT.Setup();
  MatFFT.Setup();
  cerr << "done.\n";

  kPH.CalcTailCoefs (15.0, 60.0);
  int nx, ny, nz;
  MatFFT.GetDims(nx,ny,nz);
  Fr.resize(nx,ny,nz);
  Vr.resize(nx,ny,nz);
  Gc.resize(GVecs.size());
  Vc.resize(GVecs.size());
  Vk.resize(GVecs.DeltaSize());
  Fk.resize(GVecs.DeltaSize());
//   StructureFactor.resize(GVecs.DeltaSize());
//   StructureFactor = complex<double>(1.0, 0.0);

  SetupkPotentials();
  SetuprPotentials();

  IsSetup = true;
}
  
void
PHPotFFTClass::Setk (Vec3 k)
{
  kPoint = k;
  VecFFT.Setup();
  MatFFT.Setup();
  
  int nx, ny, nz;
  MatFFT.GetDims(nx,ny,nz);
  Fr.resize(nx,ny,nz);
  Vr.resize(nx,ny,nz);
  Gc.resize(GVecs.size());
  Vc.resize(GVecs.size());
  Vk.resize(GVecs.DeltaSize());
  Fk.resize(GVecs.DeltaSize());

  SetupkPotentials();
  SetIons(Rions);
}


void 
PHPotFFTClass::Apply (const zVec &c, zVec &Hc)
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
    Gc(i) = c(i)*(kPoint+GVecs(i));
  VecFFT.PutkVec (Gc);
  VecFFT.k2r();
  
  // Now, muliply by F(r)
  for (int ix=0; ix<nx; ix++)
    for (int iy=0; iy<ny; iy++)
      for (int iz=0; iz<nz; iz++) {
	register complex<double> vx = VecFFT.rBox(ix,iy,iz)[0];
	register complex<double> vy = VecFFT.rBox(ix,iy,iz)[1];
	register complex<double> vz = VecFFT.rBox(ix,iy,iz)[2];
	register complex<double> Fvx, Fvy, Fvz;
	Fvx  = Fr(ix,iy,iz)(0,0) * vx;
	Fvy  = Fr(ix,iy,iz)(1,0) * vx;
	Fvz  = Fr(ix,iy,iz)(2,0) * vx;
	Fvx += Fr(ix,iy,iz)(0,1) * vy;
	Fvy += Fr(ix,iy,iz)(1,1) * vy;
	Fvz += Fr(ix,iy,iz)(2,1) * vy;
	Fvx += Fr(ix,iy,iz)(0,2) * vz;
	Fvy += Fr(ix,iy,iz)(1,2) * vz;
	Fvz += Fr(ix,iy,iz)(2,2) * vz;
	VecFFT.rBox(ix,iy,iz) =
	TinyVector<complex<double>,3>(Fvx,Fvy,Fvz);

//	VecFFT.rBox(ix,iy,iz) = Fr(ix,iy,iz)*VecFFT.rBox(ix,iy,iz);

// 	VecFFT.rBox(ix,iy,iz)[0] = 
// 	  Fr(ix,iy,iz)(0,0)*VecFFT.rBox(ix,iy,iz)[0] + 
// 	  Fr(ix,iy,iz)(0,1)*VecFFT.rBox(ix,iy,iz)[1] + 
// 	  Fr(ix,iy,iz)(0,2)*VecFFT.rBox(ix,iy,iz)[2];
// 	VecFFT.rBox(ix,iy,iz)[1] = 
// 	  Fr(ix,iy,iz)(1,0)*VecFFT.rBox(ix,iy,iz)[0] + 
// 	  Fr(ix,iy,iz)(1,1)*VecFFT.rBox(ix,iy,iz)[1] + 
// 	  Fr(ix,iy,iz)(1,2)*VecFFT.rBox(ix,iy,iz)[2];
// 	VecFFT.rBox(ix,iy,iz)[2] = 
// 	  Fr(ix,iy,iz)(2,0)*VecFFT.rBox(ix,iy,iz)[0] + 
// 	  Fr(ix,iy,iz)(2,1)*VecFFT.rBox(ix,iy,iz)[1] + 
// 	  Fr(ix,iy,iz)(2,2)*VecFFT.rBox(ix,iy,iz)[2];
      }
  // Transform back;
  VecFFT.r2k();
  // Rename for clarity
  zVecVec& FGc = Gc;
  VecFFT.GetkVec (FGc);
  for (int i=0; i<GVecs.size(); i++)
    Hc(i) += 0.5*dot(GVecs(i)+kPoint, FGc(i));

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
HamiltonianClass::SetIonPot (double z, bool useFFT)
{
  if (useFFT)
    Vion = new CoulombFFTClass(z, GVecs, FFT);
  else
    Vion = new CoulombClass (z, GVecs);
}


void 
HamiltonianClass::SetIonPot (Potential &ph, bool useFFT)
{
  if (useFFT)
    Vion = new PHPotFFTClass (ph, GVecs, FFT);
  else
    Vion = new PHPotClass (ph, GVecs);
}


void 
HamiltonianClass::SetIons (const Array<Vec3,1>& rions)
{
  Vion->SetIons(rions);
}

void
HamiltonianClass::Setk (Vec3 k)
{
  kPoint = k;
  Kinetic.Setk(k);
  Vion->Setk(k);
}

void
HamiltonianClass::Apply (const zVec &c, zVec &Hc)
{
  Hc = 0.0;
  Kinetic.Apply(c, Hc);
  Vion->Apply (c, Hc);
}
