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

#ifndef COULOMB_FFT_H
#define COULOMB_FFT_H

#include "HamiltonianBase.h"

class CoulombFFTClass : public VionBase
{
private:
  double Z;

  void Setup();
  void SetVr();
  FFTBox &FFT;
public:
  void Apply (const zVec &c, zVec &Hc);
  void Apply (const zVec &c, zVec &Hc,
	      const Array<double,3> &VHXC);
  void SetIons (const Array<Vec3,1> &rions);
  void Vmatrix (Array<complex<double>,2> &vmat);

  CoulombFFTClass (double z, GVecsClass &gvecs, FFTBox &fft) : 
    VionBase (gvecs), Z(z), FFT(fft)
  {
    // nothing for now
  }
};

#endif
