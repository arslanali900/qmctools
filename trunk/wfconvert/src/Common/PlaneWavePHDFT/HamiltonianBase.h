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

#ifndef HAMILTONIAN_BASE_H
#define HAMILTONIAN_BASE_H

#include "VectorOps.h"
#include "Common/PH/kSpacePH.h"
#include "FFTBox.h"
#include "config.h"

#if !HAVE_SINCOS
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

  PWKineticClass& operator=(const PWKineticClass &kinetic);
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
  /// This is the potential for single atom evaluated at the
  /// reciporical lattice vectors
  Array<double,1> VG;
  Array<complex<FFT_FLOAT>,3> Vr;
  double VG0;
public:
  // Adds H*c to the Hc vector.  Does not zero Hc before accumulating
  virtual void Apply   (const zVec &c, zVec &Hc) = 0;
  /// This version includes the hartree and exchange-correlation pot
  virtual void Apply   (const zVec &c, zVec &Hc, 
			Array<double,3> &VHXC);
  virtual double NonlocalEnergy (const zVec &c);
  virtual void SetIons (const Array<Vec3,1> &rions); 
  virtual void Setup() = 0;
  virtual void Vmatrix (Array<complex<double>,2> &vmat) = 0;
  inline const Array<Vec3,1> GetRions() {
    return Rions;
  }
  virtual void SetProjectors (bool smooth);

  inline const Array<double,1>& GetVG()
  { return VG; }

  inline double GetVG0() { return VG0; }

  inline const Array<complex<FFT_FLOAT>,3> & GetVr()
  { return Vr; }
  
  virtual void Setk    (Vec3 k) 
  { kPoint = k; }
  
  VionBase (GVecsClass &gvecs) 
    : GVecs(gvecs), IsSetup(false), kPoint (0.0, 0.0, 0.0)
  {
    // do nothing for now
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
  void Apply (const zVec &c, zVec &Hc,
	      Array<double,3> &VHXC);
  
  void SetIonPot (double z, bool useFFT);
  void SetIonPot (Potential &ph, bool useFFT);
  void SetIons (const Array<Vec3,1>& rions);
  void Setk (Vec3 k);
  inline const Array<Vec3,1> GetRions()
  {  return Vion->GetRions(); }
  inline const Array<double,1> GetVG()
  {  return Vion->GetVG();    }
  inline double GetVG0()
  {  return Vion->GetVG0();   }
  inline const Array<complex<FFT_FLOAT>,3> GetVr()
  {  return Vion->GetVr();    }

  HamiltonianClass& operator=(const HamiltonianClass &h);
  HamiltonianClass (GVecsClass &gvecs, FFTBox &fft) 
    : GVecs(gvecs), Kinetic(gvecs), Vion(NULL), FFT(fft)
  {

  }
};



#endif
