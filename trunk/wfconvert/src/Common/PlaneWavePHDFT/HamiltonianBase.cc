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

#include "HamiltonianBase.h"

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


double
VionBase::NonlocalEnergy(const zVec &c)
{
  return 0.0;
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
      // HACK HACK HACK -- trying minus sign
      phase = dot (GVecs.DeltaG(gi), rions(i));
      sincos (phase, &s, &c);
      StructureFactor(gi) += complex<double>(c, s);
    }
  }
}

void
VionBase::SetProjectors (bool smooth)
{
  // do nothing
}

void
VionBase::Apply   (const zVec &c, zVec &Hc, 
		   Array<double,3> &VHXC)
{
  cerr << "Applying VHXC potential not supported by non-FFT versions"
       << " of the potentials.\n"; 
  abort();
}


HamiltonianClass&
HamiltonianClass::operator=(const HamiltonianClass& h)
{
  Kinetic = h.Kinetic;
  Vion    = h.Vion;
  kPoint  = h.kPoint;
  return *this;
}

PWKineticClass& 
PWKineticClass::operator=(const PWKineticClass &kinetic)
{
  halfG2.resize (kinetic.halfG2.shape());
  halfG2 = kinetic.halfG2;
  kPoint = kinetic.kPoint;
  IsSetup = kinetic.IsSetup;
  return *this;
}
