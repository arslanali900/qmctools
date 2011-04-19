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

void TestSiAtom()
{
  Potential *barePot;
  IOSectionClass in;
  assert (in.OpenFile ("Si_CASINO.h5"));
  barePot = ReadPotential (in);
  in.CloseFile();

  OptimalGrid grid(4.0, 90.0);
  DFTAtom atom;
  atom.RadialWFs.resize(2);
  atom.RadialWFs(0).n = 1;
  atom.RadialWFs(0).l = 0;
  atom.RadialWFs(0).Occupancy = 2.0;
  atom.RadialWFs(0).Energy = -0.3;
  atom.RadialWFs(1).n = 2;
  atom.RadialWFs(1).l = 1;
  atom.RadialWFs(1).Occupancy = 2.0;
  atom.RadialWFs(1).Energy = -0.3;
  atom.SetGrid (&grid);
  atom.SetBarePot (barePot);

  atom.NewMix = 0.9;
  atom.Solve();
  atom.RadialWFs(0).Normalize();
  atom.RadialWFs(1).Normalize();
  FILE *fout = fopen ("SiWF.dat", "w");
  for (int i=0; i<grid.NumPoints; i++) {
    double r = grid(i);
    fprintf (fout, "%1.16e ", r);
    for (int j=0; j<atom.RadialWFs.size(); j++)
      fprintf (fout, "%1.16e ", atom.RadialWFs(j).u(i));
    fprintf (fout, "\n");
  }
}

void TestNaAtom()
{
  Potential *barePot;
  IOSectionClass in;
  assert (in.OpenFile ("Na_CASINO.h5"));
  barePot = ReadPotential (in);
  in.CloseFile();

  OptimalGrid grid(1.0, 90.0);
  DFTAtom atom;
  atom.RadialWFs.resize(2);
  atom.RadialWFs(0).n = 1;
  atom.RadialWFs(0).l = 0;
  atom.RadialWFs(0).Occupancy = 1.0;
  atom.RadialWFs(0).Energy = -0.3;
  atom.RadialWFs(1).n = 2;
  atom.RadialWFs(1).l = 1;
  atom.RadialWFs(1).Occupancy = 0.0;
  atom.RadialWFs(1).Energy = -0.3;
  atom.SetGrid (&grid);
  atom.SetBarePot (barePot);
  atom.RadialWFs(0).Normalize();
  atom.RadialWFs(1).Normalize();

  atom.NewMix = 0.9;
  atom.Solve();
  FILE *fout = fopen ("NaWF.dat", "w");
  for (int i=0; i<grid.NumPoints; i++) {
    double r = grid(i);
    fprintf (fout, "%1.16e ", r);
    for (int j=0; j<atom.RadialWFs.size(); j++)
      fprintf (fout, "%1.16e ", atom.RadialWFs(j).u(i));
    fprintf (fout, "\n");
  }
}



main()
{
  TestSiAtom();
  TestNaAtom();
}
  
  
