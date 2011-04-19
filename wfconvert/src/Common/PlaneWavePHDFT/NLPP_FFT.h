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

#ifndef NLPP_FFT_H
#define NLPP_FFT_H

#include "HamiltonianBase.h"
#include "../PH/NLPP.h"
#include "../PH/CoulombPot.h"
#include "../PH/SplinePot.h"

// Implements the King-Smith method for applying nonlocal
// pseudopotentials in real space.  See PRB 44, 13063.

class KingSmithProjector
{
public:
  int l;
  double R0;
  CubicSpline chi_r;
  double El;

  void Setup (int l, CubicSpline &chi, double R0,
	      double El);
};


// This class stores what is necessary to do Kleinman-Bylander
// projection in real space on a single ion for one l channel.
class Ion_l_Projector
{
private:
  FFTBox *FFT;
  complex<double> Ylm (int l, int m, Vec3 r);
public:
  complex<double> Ylm2(int l, int m, Vec3 r);
  int l;
  Array<Int3, 1> FFTIndices;
  Array<int,1> FFTOffsets;
  // The first index specifies the point.  The second index specifies
  // the value of m.  Since chi_r is real and 
  // Y_l{-m} = (-1)^m conj(Y_lm), we only need to store for values of
  // m from 0 to l.  We will store for all -l to l, though.
  Array<complex<double>,2> ChiYlm;
  // The volume of a mesh element
  double MeshVol;
  void Setup(NLPPClass &nlpp, int l_,
	     Vec3 rion, Vec3 k, FFTBox &fft, bool smooth);
  // This returns the application of the projector to the contents of
  // the FFTBox in chi_psi.  Thus, for l=0, chi_psi has one element,
  // for l=1: 3 elements, l=2:  5 elements...
  void Project (Array<complex<double>,1> &chi_psi);

  // This adds the contribution of the projector to VnlPsi
  void AddToVnl (Array<complex<double>,1> &chi_psi,
		 Array<complex<double>,3> &VnlPsi);
};



class NLPP_FFTClass : public VionBase
{
private:
  // This computes the reciprocal-space version of the local potential 
  kSpacePH kPH;
  zVec Vc;
  FFTBox      &cFFT;
  void SetupkProjectors();
  void SetupkPotentials();
  void SetuprPotentials();
  // Apply the nonlocal part in reciprocal space
  void ApplyNonlocal (const zVec &c, zVec &Hc);
  NLPPClass &NLPP;
  SplinePot Vlocal;
  CoulombPot Vouter;
  // The first index is the ion number, the second is the projector number. 
  Array<Ion_l_Projector,2> Ion_l_Projectors;
  // This is where the results of applying the nonlocal parts get accumulated. 
  Array<complex<double>,3> VnlPsi;
  // The first index is ion number, the second is the projector
  // number, and the third is the plane-wave index
  // wave coefficient number
  Array<complex<double>,3> lambda_lm;
  // Apply the nonlocal parts
  void CalcVnlPsi();
  bool SmoothProjectors;
public:
  void Setup();
  void SetProjectors (bool smooth);
  void SetIons (const Array<Vec3, 1>& rions);
  void Vmatrix (Array<complex<double>,2> &vmat);
  void Setk(Vec3 k);

  void Apply (const zVec &c, zVec &Hc);
  void Apply (const zVec &c, zVec &Hc,
	      Array<double,3> &VHXC);
  double NonlocalEnergy (const zVec &c);

  NLPP_FFTClass (NLPPClass &nlpp, GVecsClass &gvecs, FFTBox &fft,
		 bool smoothProjectors=false) :
    VionBase (gvecs), kPH(Vlocal), cFFT(fft), NLPP(nlpp),
    SmoothProjectors(smoothProjectors)
  {

  }
};

#endif
