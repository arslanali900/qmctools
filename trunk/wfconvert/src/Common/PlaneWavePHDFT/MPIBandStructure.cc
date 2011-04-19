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

#include "MPIBandStructure.h"

void
MPIBandStructureClass::Read(IOSectionClass &in)
{
  WorldComm.SetWorld();
  // Read the pseudoHamiltonian
  assert (in.OpenSection ("Potential"));
  PH = ReadPotential (in);
  double Z1Z2 = round (50.0*PH->V(50.0));
  cerr << "Z = " << -Z1Z2 << endl;
  IonPot.Z1Z2 = Z1Z2*Z1Z2;
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
  assert (in.ReadVar ("NumElecs", NumElecs));
  assert (in.ReadVar ("UseLDA", UseLDA));
	  
  // Read the k-point list
  Array<double,2> tmpkPoints;
  assert (in.ReadVar("kPoints", tmpkPoints));
  kPoints.resize(tmpkPoints.rows());
  for (int i=0; i<kPoints.size(); i++)
    for (int j=0; j<3; j++)
      kPoints(i)[j] = tmpkPoints(i,j);

  if (UseLDA) {
    assert (in.ReadVar("kMesh", tmpkPoints));
    kMesh.resize(tmpkPoints.rows());
    for (int i=0; i<kMesh.size(); i++) {
      kMesh(i)[0] = tmpkPoints(i,0);
      kMesh(i)[1] = tmpkPoints(i,1);
      kMesh(i)[2] = tmpkPoints(i,2);
    }
    // Setup bandComm
    int procsPerk = WorldComm.NumProcs()/kMesh.size();
    if ((BandComm.NumProcs()%kMesh.size()) != 0) {
      cerr << "The number of processors is not a multiple"
	   << " of the number of k-points.\n";
      abort();
    }
    int color = WorldComm.MyProc()/procsPerk;
    WorldComm.Split (color, BandComm);
    Array<int,1> kCommRanks(kMesh.size());
    for (int ki=0; ki<kMesh.size(); ki++)
      kCommRanks(ki) = ki*procsPerk;
    WorldComm.Subset (kCommRanks, kComm);
    MykPoint = WorldComm.MyProc()/procsPerk;
  }    
  else {
    BandComm.SetWorld();
    kComm.SetWorld();
  }
  

  // Read the number of interpolation points per kPoint
  assert (in.ReadVar ("InterpPoints", InterpPoints));

  // Read the k-cutoff
  assert(in.ReadVar("kCut", kCut));
  System = new MPISystemClass(NumBands, NumElecs,
			      BandComm, kComm);



  // Read the output file name
  assert (in.ReadVar ("OutFilename", OutFilename));

  // Setup the plane-wave system
  if (UseLDA)
    System->Setup (Box, kMesh(MykPoint), kCut, *PH, IonPot, UseLDA);
  else
    System->Setup (Box, kPoints(0), kCut, *PH, IonPot, UseLDA);
  System->SetIons(Rions);

  bool smoothProjectors;
  in.ReadVar("SmoothProjectors", smoothProjectors, false);
  if (smoothProjectors) 
    cerr << "Using King-Smith projector smoothing.\n";
  else
    cerr << "Not using King-Smith projector smoothing.\n";
  System->SetProjectors (smoothProjectors);
  

}

void
MPIBandStructureClass::CalcLDABands()
{
  /// First, do the full, self-consistent LDA calculation with the
  /// k-Mesh points:
  System->Setk(kMesh(MykPoint));
  System->SolveLDA();

  /// Now VHXC is set.  Do a non-self-consistent band-structure sweep,
  /// computing eigenenergies:
  FILE *fout;
  if (MykPoint == 0) {
    if (WorldComm.MyProc() == 0) {
      fout = fopen (OutFilename.c_str(), "w");
      assert(fout != NULL);
    }
    for (int ki=0; ki<kPoints.size()-1; ki++) 
      for (int i=0; i<InterpPoints; i++) {
	double alpha = (double)i/(double)InterpPoints;
	Vec3 k = (1.0-alpha)*kPoints(ki) + alpha*kPoints(ki+1);
	System->Setk_NonSC(k);
	System->DiagonalizeH();
	if (WorldComm.MyProc() == 0) {
	  fprintf (fout, "%1.6e %1.6e %1.6e ", k[0], k[1], k[2]);
	  for (int band=0; band<NumBands; band++)
	    fprintf (fout, "%1.16e ", System->GetEnergy(band));
	  fprintf (fout, "\n");
	  fflush(fout);
	}
      }
    if (WorldComm.MyProc() == 0)
      fclose(fout);
  }
}


void
MPIBandStructureClass::CalcBareBands()
{
  FILE *fout;
  if (BandComm.MyProc()==0) {
    fout = fopen (OutFilename.c_str(), "w");
    assert (fout != NULL);
  }
 
  for (int ki=0; ki<kPoints.size()-1; ki++) 
    for (int i=0; i<InterpPoints; i++) {
      double alpha = (double)i/(double)InterpPoints;
      Vec3 k = (1.0-alpha)*kPoints(ki) + alpha*kPoints(ki+1);
      if (BandComm.MyProc()==0) 
	fprintf (fout, "%1.6e %1.6e %1.6e ", k[0], k[1], k[2]);
      System->Setk(k);
      System->DiagonalizeH();
      if (BandComm.MyProc()==0) {
	for (int band=0; band<NumBands; band++)
	  fprintf (fout, "%1.16e ", System->GetEnergy(band));
	fprintf (fout, "\n");
	fflush (fout);
      }
    }
  if (BandComm.MyProc()==0)
    fclose(fout);
}


void
MPIBandStructureClass::CalcBands()
{
  if (UseLDA)
    CalcLDABands();
  else
    CalcBareBands();
}

#include <time.h>

main(int argc, char **argv)
{
  COMM::Init(argc, argv);
  if (argc < 2) {
    cout << "Usage:\n";
    cout << "BandStructure myfile.in\n"; 
  }
  else {
  clock_t start, end;
    IOSectionClass in;
    assert (in.OpenFile(argv[1]));
    MPIBandStructureClass BandStructure;
    BandStructure.Read(in);
    in.CloseFile();
    start = clock();
    BandStructure.CalcBands();
    end = clock();
    cerr << "Time = " << ((double)(end-start)/(double)CLOCKS_PER_SEC) << endl;
  }
}
