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

#ifndef CORE_TRANSFORM_H
#define CORE_TRANSFORM_H

#include "Potential.h"

/// The core transform class describes a coordinate transformation
/// used in calculating the pair density matrix through matrix
/// squaring for pseudohamiltonians .  It maps between real radial
/// distances, r, and transformed distances, x.  This transformation
/// allows the radial equation for the pseudohamiltonians to have the
/// same canonical form as a local radial potential, thereby allowing
/// us to use the traditional matrix squaring formalism.
class CoreTransform
{
  double rMax, xMax;
  LinearGrid xgrid, rgrid;
  CubicSpline x_of_r;
  CubicSpline r_of_x;

public:
  /// Convert from x coordinate to r coordinate
  inline double x2r(double x)
  {
    if (x >= xMax)
      return ((x-xMax)+rMax);
    else
      return (r_of_x(x));
  }

  /// Convert from r coordinate to x coordinate
  inline double r2x(double r)
  {
    if (r >= rMax)
      return ((r-rMax)+xMax);
    else
      return (x_of_r(r));
  }
  
  /// Calculate the transform from the pseudohamiltonian.
  void Initialize(Potential *pot, int NumPoints);
};

#endif
