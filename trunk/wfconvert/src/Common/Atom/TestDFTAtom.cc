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

#include "NewAtom.h"

void TestDFTAtom()
{
  CoulombPot barePot;
  barePot.Z1Z2 = -1.0;
  OptimalGrid grid(1.0, 50.0);
  DFTAtom atom;
  atom.RadialWFs.resize(1);
  atom.RadialWFs(0).n = 1;
  atom.RadialWFs(0).l = 0;
  atom.RadialWFs(0).Occupancy = 1.0;
  atom.RadialWFs(0).Energy = -0.3;
  atom.SetGrid (&grid);
  atom.SetBarePot (&barePot);

  atom.NewMix = 0.9;
  atom.Solve();
}



main()
{
  TestDFTAtom();
}
  
  
