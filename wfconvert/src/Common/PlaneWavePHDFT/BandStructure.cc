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

#include "BandStructure.h"

void
BandStructureClass::Read(IOSectionClass &in)
{
  // Read the pseudoHamiltonian
  assert (in.OpenSection ("Potential"));
  PH = ReadPotential (in);
  in.CloseSection();

  // Read the box
  Array<double,1> box;
  assert (in.ReadVar ("Box", box));
  Box[0]=box(0);   Box[1]=box(1);    Box[2]=box(2);
  
  // Read ion positions
  Array<double,2> tmpIons;
  assert (in.ReadVar ("Positions", tmpIons));
  Rions.resize(tmpIons.rows());
  for (int i=0; i<Rions.size(); i++)
    for (int j=0; j<3; j++)
      Rions(i)[j] = tmpIons(i,j);

  // Read the number of bands
  assert (in.ReadVar ("NumBands", NumBands));
	  
  // Read the k-point list
  Array<double,2> tmpkPoints;
  assert (in.ReadVar("kPoints", tmpkPoints));
  kPoints.resize(tmpkPoints.rows());
  for (int i=0; i<kPoints.size(); i++)
    for (int j=0; j<3; j++)
      kPoints(i)[j] = tmpkPoints(i,j);

  // Read the number of interpolation points per kPoint
  assert (in.ReadVar ("InterpPoints", InterpPoints));

  // Read the k-cutoff
  assert(in.ReadVar("kCut", kCut));
  System = new SystemClass(NumBands);

  // Read the output file name
  assert (in.ReadVar ("OutFilename", OutFilename));

  // Setup the plane-wave system
  System->Setup (Box, kPoints(0), kCut, *PH, true);
  System->SetIons(Rions);


};


void
BandStructureClass::CalcBands()
{
  FILE *fout = fopen (OutFilename.c_str(), "w");
  assert (fout != NULL);
 
  for (int ki=0; ki<kPoints.size()-1; ki++) 
    for (int i=0; i<InterpPoints; i++) {
      double alpha = (double)i/(double)InterpPoints;
      Vec3 k = (1.0-alpha)*kPoints(ki) + alpha*kPoints(ki+1);
      fprintf (fout, "%1.6e %1.6e %1.6e ", k[0], k[1], k[2]);
      System->Setk(k);
      System->DiagonalizeH();
      for (int band=0; band<NumBands; band++)
	fprintf (fout, "%1.16e ", System->GetEnergy(band));
      fprintf (fout, "\n");
      fflush (fout);
    }
  fclose(fout);
}

#include <time.h>
main(int argc, char **argv)
{
  if (argc < 2) {
    cout << "Usage:\n";
    cout << "BandStructure myfile.in\n"; 
  }
  else {
    clock_t start, end;

    IOSectionClass in;
    assert (in.OpenFile(argv[1]));
    BandStructureClass BandStructure;
    BandStructure.Read(in);
    in.CloseFile();
    start = clock();
    BandStructure.CalcBands();
    end = clock();
    cerr << "Time = " << ((double)(end-start)/(double)CLOCKS_PER_SEC) << endl;
  }
}
