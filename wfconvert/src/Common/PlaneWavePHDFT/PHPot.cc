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

#include "PHPot.h"

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

