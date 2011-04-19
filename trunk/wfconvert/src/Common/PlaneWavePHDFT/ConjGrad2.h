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

#ifndef CONJ_GRAD2_H
#define CONJ_GRAD2_H

#include "Hamiltonians.h"

class ConjGrad
{
protected:
  HamiltonianClass &H;
  zVec c;
  zVec cnext, Hc, Phi;
  zVec Phip, Phipp, Xi, Eta;
  complex<double> EtaXiLast;
  double E0;
  double CalcPhiSD();
  double CalcPhiCG();
  bool IsSetup;
  void Precondition();
  double T;
  int iter;
  int CurrentBand, LastBand;
  Array<complex<double>,2> &Bands;
  double Tolerance;

public:
  void Setup();
  Array<double,1> Energies;
  void InitBands();
  void Solve(int band);
  inline void SetTolerance(double tol) { Tolerance = tol;}
  void PrintOverlaps();

  ConjGrad (HamiltonianClass &h, Array<complex<double>,2> &bands) : 
    H(h), IsSetup(false), EtaXiLast(0.0, 0.0), iter(0),
    LastBand(-1), Bands(bands), Tolerance (1.0e-6)
  {
    // Do nothing for now
  }
};

#endif
