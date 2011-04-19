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

#include "RadialWF.h"

void TestRadialWF1()
{
  CoulombPot pot;
  pot.Z1Z2 = -1.0;

  OptimalGrid grid(1.0, 80.0);
  RadialWF wf;
  wf.l = 0;
  wf.n = 1;
  wf.Energy = -0.25;
  wf.Occupancy = 1.0;
  wf.SetPotential (&pot);
  wf.SetGrid (&grid);
  
  wf.Solve (1.0e-14);
  fprintf (stderr, "Energy = %1.12f\n", wf.Energy);
}

void TestRadialWF_TH()
{
  ToppHopfieldPot pot;
  pot.Z = 1.0;
  pot.a = 1.22441639097900;
  pot.b = -0.17904097875311;
  pot.rc = 3.0;
  pot.V0 = 0.179;

  OptimalGrid grid(1.0, 100.0);
  RadialWF wf;
  wf.l = 0;
  wf.n = 1;
  wf.Energy = -0.25;
  wf.Occupancy = 1.0;
  wf.SetPotential (&pot);
  wf.SetGrid (&grid);
  
  wf.Solve (1.0e-14);
  fprintf (stderr, "Topp-Hopfield energy = %1.12f\n", wf.Energy);
  IOSectionClass out;
  out.NewFile ("ToppHopfield_WF.h5");
  out.WriteVar ("u", wf.u.Data());
  out.WriteVar ("r", wf.grid->Points());
  out.CloseFile();

}

void TestRadialWF_Needs()
{
  Potential *pot;
  IOSectionClass in;
  assert(in.OpenFile("Na_HF_NLPP.h5"));
  pot = ReadPotential(in);

  in.CloseFile(); 
  OptimalGrid grid(1.0, 100.0);
  RadialWF wf;
  wf.l = 0;
  wf.n = 1;
  wf.Energy = -0.25;
  wf.Occupancy = 1.0;
  wf.SetPotential (pot);
  wf.SetGrid (&grid);
  
  wf.Solve (1.0e-14);
  fprintf (stderr, "Needs HF NLPP energy = %1.12f\n", wf.Energy);
  IOSectionClass out;
  out.NewFile ("Na_HF_WF.h5");
  out.WriteVar ("u", wf.u.Data());
  out.WriteVar ("r", wf.grid->Points());
  out.CloseFile();
}

void TestRadialWF_LocalPH()
{
  Potential *pot;
  IOSectionClass in;
  assert(in.OpenFile("NaLocalPH.h5"));
  pot = ReadPotential(in);

  in.CloseFile(); 
  OptimalGrid grid(1.0, 100.0);
  RadialWF wf;
  wf.l = 0;
  wf.n = 1;
  wf.Energy = -0.25;
  wf.Occupancy = 1.0;
  wf.SetPotential (pot);
  wf.SetGrid (&grid);
  
  wf.Solve (1.0e-14);
  fprintf (stderr, "NaLocalPH energy     = %1.12f\n", wf.Energy);
  IOSectionClass out;
  out.NewFile ("NaLocalPH_WF.h5");
  out.WriteVar ("u", wf.u.Data());
  out.WriteVar ("r", wf.grid->Points());
  out.CloseFile();
}



main()
{
  TestRadialWF1();
  TestRadialWF_TH();
  TestRadialWF_Needs();
  TestRadialWF_LocalPH();
}
