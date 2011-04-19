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

#ifndef GVECS_H
#define GVECS_H

#include "VectorOps.h"

class CommunicatorClass;

class GVecsClass
{
protected:
  // This stores the actual G-Vectors
  Array<Vec3,1> GVecs;
  // This stores the index of this G-vector in the FFT box
  Array<Int3,1> Indices;
  // This stores the integer multipliers of the reciprocal lattice
  // vector for each G-vector
  Array<Int3,1> Multipliers;
  // This stores the differences between g-vectors:
  Array<Vec3,1> GDiff;
  /// Lattice[0,:] gives the first primitive lattice vector, etc.
  Mat3 Lattice, LatticeInv, RecipLattice;
  // The inverse of the norm of the above:  used to compute the
  // Hartree potential in DFT
  Array<double,1> GDiffInv2;
  Array<Int3,1> IDiff;
  Vec3 Box, kBox;
  int Nx, Ny, Nz;
  double kCut;
  Vec3 k;

public:
  Int3 GetFFTBoxSize (Vec3 box, Vec3 kvec, double kcut,
		      double fftFactor=1.0);
  Int3 GetFFTBoxSize (Mat3 lattice, const Array<Vec3,1> &gvecs, 
		      double fftFactor=1.0);

  void Set (Vec3 box, Vec3 kvec, double kcut, double fftFactor=1.0);
  void Set (Vec3 box, Vec3 kvec, double kcut, Int3 boxSize);
  void Set (Mat3 lattice, Vec3 kvec, double kcut, Int3 boxSize);
  void Set (Mat3 lattice, const Array<Vec3,1> &gvecs, 
	    double fftFactor=1.0);
  void Set (Mat3 lattice, const Array<Vec3,1> &gvecs, 
	    Int3 boxSize);
  void Broadcast (CommunicatorClass &comm, int root);

  inline Vec3 operator()(int i) const
  { return GVecs(i); }

  inline const Vec3& operator()(int i) 
  { return GVecs(i); }

  inline Int3 Index(int i) const
  { return Indices(i); }

  inline Int3 Multiplier (int i) const
  { return Multipliers(i); }

  inline const Array<Int3,1>& GetMultipliers() const
  { return Multipliers; }

  inline Vec3 DeltaG (int i) const
  { return GDiff(i); }

  inline const Vec3& DeltaG (int i)
  { return GDiff(i); }

  inline double DeltaGInv2 (int i)
  { return GDiffInv2(i); }

  inline Int3 DeltaI (int i) const
  { return IDiff(i); }

  inline int size()
  { return GVecs.size(); }

  inline int DeltaSize() 
  { return GDiff.size(); }

  inline double GetBoxVol() 
  { return fabs(det(Lattice)); }

  inline double GetkCut() 
  { return kCut; }

  inline Vec3 GetBox()
  { return Box; }

  inline Mat3 GetLattice() const
  { return Lattice; }

  inline Mat3 GetLatticeInv() const
  { return LatticeInv; }

  inline Mat3 GetRecipLattice() const
  { return RecipLattice; }

  inline Vec3 GetkBox()
  { return kBox; }

  inline Vec3 Getk()
  { return k; }

  void GetFFTBoxSize (int &nx, int &ny, int &nz);

  GVecsClass& operator=(const GVecsClass &gvecs);
};


#endif

