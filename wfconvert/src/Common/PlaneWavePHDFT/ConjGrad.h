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

#ifndef CONJ_GRAD_H
#define CONJ_GRAD_H

#include "Hamiltonian.h"

class ConjGrad
{
protected:
  Hamiltonian &H;
  zVec c;
  zVec cnext, Hc, Phi;
  zVec Phip, Phipp, Xi, Eta;
  complex<double> EtaXiLast;
  double E0;
  void CalcPhiSD();
  void CalcPhiCG();
  void Setup();
  bool IsSetup;
  void Precondition();
  double T;
  int iter, NumBands;
  int CurrentBand, LastBand;
public:
  Array<complex<double>,2> Bands;
  Array<double,1> Energies;
  void Iterate(int band);
  void PrintOverlaps();

  ConjGrad (Hamiltonian &h, int numBands) : 
    H(h), IsSetup(false), EtaXiLast(0.0, 0.0), iter(0), NumBands(numBands),
    LastBand(-1)
  {
    // Do nothing for now
  }
};

#endif
