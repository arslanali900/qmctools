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

#ifndef PLANE_WAVES_H
#define PLANE_WAVES_H

#include "ConjGrad2.h"

class SystemClass
{
protected:
  Array<Vec3,1> Rions;
  FFTBox FFT;
  HamiltonianClass H;
  Array<complex<double>,2> Bands;
  ConjGrad CG;
  Vec3 Box;
  double kCut;
  int NumBands;
  Potential *PH;
  bool UseFFT;
public:
  GVecsClass GVecs;
  void Setup(Vec3 box, Vec3 k, double kcut, Potential &ph, bool useFFT=true);
  void Setup(Vec3 box, Vec3 k, double kcut, double z, bool useFFT=true);
  void SetIons (const Array<Vec3,1> &rions);
  inline Vec3 GetIonPos(int i) { return Rions(i); }
  void Setk (Vec3 k);
  void DiagonalizeH();
  inline double GetEnergy(int band) { return CG.Energies(band); }
  inline int GetNumBands() { return NumBands; }

  /// Gets the FFT box dimensions.
  inline void GetBoxDims(int &nx, int &ny, int &nz)
  { FFT.GetDims(nx, ny, nz); }
  /// This FFT's the desired band into real space.
  void SetRealSpaceBandNum(int num);
  inline complex<double> RealSpaceBand (int ix, int iy, int iz)
  { return FFT.rBox(ix, iy, iz); }
  

  void CalcChargeDensity(Array<double,3> &rho);
  void WriteXSFFile(string filename);

  SystemClass(int numBands) 
    : CG(H, Bands), FFT(GVecs), H(GVecs, FFT), NumBands(numBands)
  {

  }
};

#endif
