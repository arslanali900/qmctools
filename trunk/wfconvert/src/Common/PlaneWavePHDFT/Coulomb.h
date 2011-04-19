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

#ifndef COULOMB_H
#define COULOMB_H

#include "HamiltonianBase.h"

class CoulombClass : public VionBase
{
private:
  double Z;
public:
  void Setup();
  void Apply (const zVec &c, zVec &Hc);
  void SetIons (const Array<Vec3,1> &rions);
  void Vmatrix (Array<complex<double>,2> &vmat);

  CoulombClass (double z, GVecsClass &gvecs) : 
    VionBase (gvecs), Z(z)
  {
    // nothing for now
  }
};

#endif
