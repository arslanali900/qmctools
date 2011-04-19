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

#ifndef CONJ_GRAD_MPI_H
#define CONJ_GRAD_MPI_H

#include "Hamiltonians.h"
#include "ChargeMixer.h"

class CommunicatorClass;

typedef enum {ORTHO_ALL, ORTHO_LOWER} OrthoType;

class ConjGradMPI
{
protected:
  bool Verbose;
  HamiltonianClass &H;
  zVec c;
  zVec cnext, Hc, Phi;
  zVec Phip, Phipp, Xi, Eta;
  Array<complex<double>,1> EtaXiLast;
  double E0;
  double CalcPhiSD();
  double CalcPhiCG();
  bool IsSetup;
  void Precondition();
  int iter;
  int CurrentBand, LastBand;
  Array<complex<double>,2> &Bands;
  Array<complex<double>,2> &HBands;
  int NumOccupied, MaxIter;
  Array<complex<double>,2> lastPhis;
  double Tolerance;
  CommunicatorClass &BandComm, &kComm;
  /// Returns the maximum residual of all the bands
  double Iterate();
  int MyFirstBand, MyLastBand;
  Array<double,1> Residuals;
  void CheckOverlaps();
  /// For LDA
  Array<double,3> &VHXC;

  /// Orthogonalize CG vector to all bands
  OrthoType Ortho;
  void CollectBands();

public:
  // The kinetic energies of the bands
  Array<double,1> T;

  /// Applies the current hamiltonian to all bands, storing in Hbands
  void ApplyH();

  void Setup();
  inline void SetNumOccupied (int num) { NumOccupied = num; }
  inline void SetMaxIter     (int num) { MaxIter     = num; }

  /// Stores the band eigenenergies
  Array<double,1> Energies;
  void InitBands();
  void Solve();
  inline void SetTolerance(double tol) { Tolerance = tol;}
  inline void SetOrthoType (OrthoType otype) { Ortho = otype; }
  void PrintOverlaps();

  inline int GetFirstBand() { return MyFirstBand; }
  inline int GetLastBand () { return MyLastBand; }
  inline void SetVerbose(bool verb) { Verbose = verb; }
  ConjGradMPI& operator=(const ConjGradMPI& cg);

  ConjGradMPI (HamiltonianClass &h,
	       Array<complex<double>,2> &bands,
	       Array<complex<double>,2> &hbands,
	       int numOccupied,
	       CommunicatorClass &bandComm, 
	       CommunicatorClass &kcomm,
	       Array<double,3> &vhxc) : 
    H(h), IsSetup(false), iter(0),
    LastBand(-1), Bands(bands), HBands(hbands), 
    NumOccupied(numOccupied), Tolerance (1.0e-6),
    BandComm(bandComm), kComm(kcomm), VHXC(vhxc),
    Ortho (ORTHO_ALL), Verbose(false), MaxIter(10)
  {
    // Do nothing for now
  }
};

#endif
