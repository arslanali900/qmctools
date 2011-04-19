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

#ifndef BAND_STRUCTURE_H
#define BAND_STRUCTURE_H

#include "PlaneWaves.h"

class BandStructureClass
{
private:
  SystemClass *System;
  Array<Vec3,1> Rions;
  Potential *PH;
  double kCut;
  Array<Vec3,1> kPoints;
  int NumBands;
  Vec3 Box;
  string OutFilename;
  int InterpPoints;
public:
  void Read(IOSectionClass &in);
  void CalcBands();

};


#endif
