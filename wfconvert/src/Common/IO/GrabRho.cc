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

#include "IO.h"

using namespace blitz;
using namespace IO;

main()
{
  IOSectionClass in;
  assert (in.OpenFile ("/home/esler/Runs/NaRuns/Langevin/turing/Na_16_rho220_T2503_32kvecs_BandFix.0.h5"));
  assert (in.OpenSection("Moves"));
  assert (in.OpenSection("Langevin"));
  IOVarBase *var = in.GetVarPtr("Rho");
  assert(var);
  Array<double,3> rho;
  var->Read(rho, 200, Range::all(), Range::all(), Range::all());
  
  IOSectionClass out;
  out.NewFile("Rho200.h5");
  out.WriteVar("Rho", rho);
  out.CloseFile();

  out.NewFile("configs.h5");
  Array<double,3> Rion;
  in.ReadVar("R", Rion);
  out.WriteVar("R", Rion);
  out.CloseFile();


}
