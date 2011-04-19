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

#ifndef K_SPACE_PH_H
#define K_SPACE_PH_H

#include "../Blitz.h"

#include "PotentialBase.h"
#include <vector>
#include <map>

class kCachePoint
{
public:
  double a, bPerp, bPar, V, k;
};

/// This just stores the fourier transforms of the pseudeohamiltonian
/// functions for a single k point.  Used for caching these values in
/// an STL map.
struct kMapPoint
{
  double a, bPerp, bPar, V;
  kMapPoint(double aVal, double bPerpVal, double bParVal, double VVal)
  {
    a = aVal; bPerp=bPerpVal; bPar=bParVal; V=VVal;
  }
  kMapPoint() {}
};

/// This is a simple replacement for the "less than" operator to be
/// used in the map cache of fourier transforms.  It essentially says
/// that two doubles are the same if they are within 1e-12 of each
/// other. 
struct FuzzyLess
{
  bool operator()(const double v1, const double v2) const
  {
    return ((v1+1.0e-12)<v2);
  }
};



class kSpacePH;


/// Deprecated.  We use the STL map instead.
class kCache
{
protected:
  vector<kCachePoint> Cache;
  kSpacePH &kPH;
public:
  void GetVals (double k, double &a, double &bPerp, double &bPar,
		double &V);
  void Clear();
  kCache (kSpacePH &kph) : kPH(kph)
  {
    // do nothing for now
  }
};


/// kSpacePH class takes a PH as its constructor argument.  Its only
/// public member function, V, returns the k-space component.  The
/// formulae used in this class come from 
/// Foulkes and Schluter Phys. Rev. B 42, 11505 (1990)
class kSpacePH
{
protected:
  friend class kCache;
  map<double,kMapPoint,FuzzyLess> kMap;
  Potential &PH;
  /// The tail is fit to the potential to the form: 
  /// V = Ctail/r + Ctail2/r^2 + Ctail3/r^3
  double Ctail1, Ctail2, Ctail3;

  //////////////////////
  // Member functions //
  //////////////////////

  /// Returns fourier component of potential part
  double Vk (double k);
  /// Fourier transform of a(r)
  double a  (double k);
  /// The perpendicular component of the b tensor
  double bPerp (double k);
  /// The perpendicular component of the b tensor
  double bPar  (double k);
  bool HaveTailCoefs;
  double R1, R2;
  kCache Cache;
  bool UseCache;
public:
  /// Calculates the values of Ctail1-3.  r1 and r2 specify the start
  /// and end of the region of the potential fit to the form given
  /// above. 
  void CalcTailCoefs(double r1, double r2);

  /// Returns the F tensor for k
  TinyMatrix<double,3,3> Ftensor (Vec3 deltaG, 
				  double aVal, double bPerpVal,
				  double bPparVal);
  TinyMatrix<double,3,3> Ftensor (Vec3 deltaG);

  void GetVals (double dGmag, 
		double &aval, double &bPerpVal, double &bParVal, double &VVal);

  /// Returns to fourier component of V, using cached values
  double V(double deltaGmag);

  /// This returns the k-space representation of the
  /// pseudohamiltonian.  Note that it does not contain the 1/Vcell
  /// volume factor, and that you MUST call CalcTailCoefs before
  /// calling this fuction.
  double V (Vec3 k, Vec3 G, Vec3 Gp);
  

  kSpacePH (Potential &ph) : PH(ph), HaveTailCoefs(false),
			     Cache (*this), UseCache(true)
  {
    // do nothing else for now
  }
};

#endif
