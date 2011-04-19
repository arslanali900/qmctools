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

#ifndef PH_POT_FFT_H
#define PH_POT_FFT_H

#include "HamiltonianBase.h"


class PHPotFFTClass : public VionBase
{
private:
  kSpacePH kPH;
  Array<Mat3,1> Fk;
  Array<TinyMatrix<complex<FFT_FLOAT>,3,3>,3> Fr;
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
  void Apply (const zVec &c, zVec &Hc, Array<double,3> &VHXC);
  PHPotFFTClass (Potential &ph, GVecsClass &gvecs, FFTBox &fft) :
    VionBase (gvecs), kPH(ph), cFFT(fft), VecFFT(gvecs), MatFFT(gvecs)
  {

  }
};

#endif
