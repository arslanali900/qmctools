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

#ifndef PH_POT_H
#define PH_POT_H

#include "HamiltonianBase.h"

class PHPotClass : public VionBase
{
private:
  kSpacePH kPH;
  void Setup();
  Array<complex<double>,2> VGGp;
  Array<complex<double>,2> StructFact;
  void CalcStructFact();
  bool VmatIsSet, SFIsSet;
  void SetVmat();
public:
  void Apply (const zVec &c, zVec &Hc);
  void SetIons (const Array<Vec3,1> &rions);
  void Vmatrix (Array<complex<double>,2> &vmat);
  void Setk(Vec3 k);
  
  PHPotClass (Potential &ph, GVecsClass &gvecs) :
    VionBase (gvecs), kPH(ph), VmatIsSet(false),
    SFIsSet(false)
  {

  }
};


#endif
