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

#include "PH.h"

//////////////////////////////////////////////////////////////////
//                FullCorePotential Routines                    //
//////////////////////////////////////////////////////////////////

void
FullCorePotential::Read(char *FName)
{
  strncpy (FileName, FName, 1000);
  FILE *fin;
  if ((fin = fopen(FileName, "r")) == NULL)
    {
      cerr << "Cannot open " << FileName << " for reading. Exitting.";
      exit(1);
    }

  fscanf (fin, " %lf ", &Z);
  Grid *TempGrid;
  TempGrid = ReadGrid(fin);
  if (GridInitialized)
    delete(V.grid);
  GridInitialized = 1;
  scalar dummy;
  Array<scalar, 1> PotVals(TempGrid->NumPoints);
  for (int i=0; i<TempGrid->NumPoints; i++)
    fscanf(fin, " %lf %lf ", &PotVals(i), &dummy);
  V.Init(TempGrid, PotVals, 5.0e30, 5.0e30);
  fclose(fin);
}
