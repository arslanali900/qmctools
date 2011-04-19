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

#include "Hamiltonians.h"

void
HamiltonianClass::SetIonPot (double z, bool useFFT)
{
  if (useFFT)
    Vion = new CoulombFFTClass(z, GVecs, FFT);
  else
    Vion = new CoulombClass (z, GVecs);
}


void 
HamiltonianClass::SetIonPot (Potential &pot, bool useFFT)
{
  if (useFFT) {
    if (pot.IsNonlocal()) {
      NLPPClass &nlpp = dynamic_cast<NLPPClass&> (pot);
      assert (&nlpp != NULL);
      Vion = new NLPP_FFTClass (nlpp, GVecs, FFT);
    }
    else if (pot.IsPH())
      Vion = new PHPotFFTClass (pot, GVecs, FFT);
    else
      Vion = new LocalPotFFTClass (pot, GVecs, FFT);
  }
  else
    Vion = new PHPotClass (pot, GVecs);
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

void
HamiltonianClass::Apply (const zVec &c, zVec &Hc,
			 Array<double,3> &VHXC)
{
  Hc = 0.0;
  Kinetic.Apply(c, Hc);
  Vion->Apply (c, Hc, VHXC);
}
