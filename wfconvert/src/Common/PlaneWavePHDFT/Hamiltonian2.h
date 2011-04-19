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

#ifndef HAMILTONIAN_H
#define HAMILTONIAN_H

#include "VectorOps.h"
#include "../PH/kSpacePH.h"
#include "FFTBox.h"

#ifdef MAC
#define sincos(p,s,c) *(s)=sin(p); *(c)=cos(p);
#endif

class PWKineticClass
{
private:
  GVecsClass &GVecs;
  Array<double,1> halfG2;
  void Setup();
  Vec3 kPoint;
  bool IsSetup;
public:
  void Apply (const zVec &c, zVec &Kc);
  void Setk  (Vec3 k);

  PWKineticClass (GVecsClass &gvecs) : 
    GVecs (gvecs), IsSetup(false), kPoint(0.0, 0.0, 0.0)
  {
  }
};



class VionBase
{
protected:
  GVecsClass &GVecs;
  bool IsSetup;
  Vec3 kPoint;
  zVec StructureFactor;
  Array<Vec3,1> Rions;
public:
  // Adds H*c to the Hc vector.  Does not zero Hc before accumulating
  virtual void Apply   (const zVec &c, zVec &Hc) = 0;
  virtual void SetIons (const Array<Vec3,1> &rions); 
  virtual void Setup() = 0;
  virtual void Vmatrix (Array<complex<double>,2> &vmat) = 0;
  
  virtual void Setk    (Vec3 k) 
  { kPoint = k; }
  
  VionBase (GVecsClass &gvecs) 
    : GVecs(gvecs), IsSetup(false), kPoint (0.0, 0.0, 0.0)
  {
    // do nothing for now
  }
};


class CoulombClass : public VionBase
{
private:
  double Z;
public:
  void Setup();
  void Apply (const zVec &c, zVec &Hc);
  void SetIons (const Array<Vec3,1> &rions);
  void Vmatrix (Array<complex<double>,2> &vmat);

  CoulombClass (double z, GVecsClass &gvecs) : 
    VionBase (gvecs), Z(z)
  {
    // nothing for now
  }
};

class CoulombFFTClass : public VionBase
{
private:
  double Z;
  Array<complex<double>,3> Vr;

  void Setup();
  void SetVr();
  FFTBox &FFT;
public:
  void Apply (const zVec &c, zVec &Hc);
  void SetIons (const Array<Vec3,1> &rions);
  void Vmatrix (Array<complex<double>,2> &vmat);

  CoulombFFTClass (double z, GVecsClass &gvecs, FFTBox &fft) : 
    VionBase (gvecs), Z(z), FFT(fft)
  {
    // nothing for now
  }
};

class PHPotClass : public VionBase
{
private:
  kSpacePH kPH;
  void Setup();
  Array<complex<double>,2> VGGp;
  Array<complex<double>,2> StructFact;
  void CalcStructFact();
  bool VmatIsSet, SFIsSet;
  void SetVmat();
public:
  void Apply (const zVec &c, zVec &Hc);
  void SetIons (const Array<Vec3,1> &rions);
  void Vmatrix (Array<complex<double>,2> &vmat);
  void Setk(Vec3 k);
  
  PHPotClass (Potential &ph, GVecsClass &gvecs) :
    VionBase (gvecs), kPH(ph), VmatIsSet(false),
    SFIsSet(false)
  {

  }
};


class PHPotFFTClass : public VionBase
{
private:
  kSpacePH kPH;
  Array<Mat3,1> Fk;
  Array<cMat3,3> Fr;
  Array<double,1> Vk;
  Array<complex<double>,3> Vr;
  zVec Vc;
  zVecVec Gc;
  FFTBox      &cFFT;
  FFTVecBox VecFFT;
  FFTMatBox MatFFT;
  void SetupkPotentials();
  void SetuprPotentials();
public:
  void Setup();
  void SetIons (const Array<Vec3, 1>& rions);
  void Vmatrix (Array<complex<double>,2> &vmat);
  void Setk(Vec3 k);

  void Apply (const zVec &c, zVec &Hc);
  PHPotFFTClass (Potential &ph, GVecsClass &gvecs, FFTBox &fft) :
    VionBase (gvecs), kPH(ph), cFFT(fft), VecFFT(gvecs), MatFFT(gvecs)
  {

  }
};



class HamiltonianClass
{
private:
public:
  FFTBox &FFT;
  GVecsClass &GVecs;
  PWKineticClass Kinetic;
  VionBase *Vion;
  Vec3 kPoint;

  void Apply (const zVec &c, zVec &Hc);
  
  void SetIonPot (double z, bool useFFT);
  void SetIonPot (Potential &ph, bool useFFT);
  void SetIons (const Array<Vec3,1>& rions);
  void Setk (Vec3 k);

  HamiltonianClass (GVecsClass &gvecs, FFTBox &fft) 
    : GVecs(gvecs), Kinetic(gvecs), Vion(NULL), FFT(fft)
  {

  }
};


#endif
