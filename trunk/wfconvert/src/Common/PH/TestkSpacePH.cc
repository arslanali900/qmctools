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

#include "kSpacePH.h"
#include "../MatrixOps/MatrixOps.h"

class kSpaceTest
{
protected:
  Vec3 Box;
  Vec3 GBox;
  Array<Vec3,1> GVecs;
  double kCut;
  Potential &PH;
  kSpacePH kPH;
public:
  Array<double,2> H;

  void SetBox(double x, double y, double z);
  void SetupGVecs(double kcut);
  void CalcHamiltonian (Vec3 k);
  void Diagonalize();
  void PlotRho (Array<double,1> phik);
  kSpaceTest (Potential &ph) : PH(ph), kPH(ph)
  {
    kPH.CalcTailCoefs (10.0, 80.0);
  }
};


void kSpaceTest::SetBox (double x, double y, double z)
{
  Box[0]  = x;      GBox[0] = 2.0*M_PI/x;
  Box[1]  = y;      GBox[1] = 2.0*M_PI/y;
  Box[2]  = z;      GBox[2] = 2.0*M_PI/z;
}

void kSpaceTest::SetupGVecs (double kcut)
{
  kCut = kcut;

  int xMax = (int)ceil(kCut/GBox[0]);  
  int yMax = (int)ceil(kCut/GBox[1]);
  int zMax = (int)ceil(kCut/GBox[2]);

  cerr << "xMax = " << xMax << endl;
  cerr << "yMax = " << yMax << endl;
  cerr << "zMax = " << zMax << endl;
  
  Vec3 G;
  int numVecs = 0;
  for (int ix=-xMax; ix<=xMax; ix++) {
    G[0] = GBox[0] * ix;
    for (int iy=-yMax; iy<=yMax; iy++) {
      G[1] = GBox[1] * iy;
      for (int iz=-zMax; iz<=zMax; iz++) {
	G[2] = GBox[2] * iz;
	if (dot(G,G) < (kCut*kCut))
	  numVecs++;
      }
    }
  }
  GVecs.resize(numVecs);
  numVecs = 0;
  for (int ix=-xMax; ix<=xMax; ix++) {
    G[0] = GBox[0] * ix;
    for (int iy=-yMax; iy<=yMax; iy++) {
      G[1] = GBox[1] * iy;
      for (int iz=-zMax; iz<=zMax; iz++) {
	G[2] = GBox[2] * iz;
	if (dot(G,G) < (kCut*kCut)) {
	  GVecs(numVecs) = G;
	  numVecs++;
	}
      }
    }
  }
}

void kSpaceTest::CalcHamiltonian(Vec3 k)
{
  int numG = GVecs.size();
  cerr << "H will be " << numG << "x" << numG << " in size.\n";
  H.resize(numG, numG);

  double volInv = 1.0/(Box[0] * Box[1] * Box[2]);

  for (int row=0; row<numG; row++) {
    Vec3 &G = GVecs(row);
    for (int col=0; col<numG; col++) {
      H(row,col) = 0.0;
      Vec3 &Gp = GVecs(col);
      if (row == col)
	H(row,col) = 0.5*dot(G+k,G+k);
      H(row,col) += (kPH.V(k, G, Gp)*volInv);
    }
  }
}


void kSpaceTest::PlotRho(Array<double,1> phik)
{
  assert (phik.size() == GVecs.size());
  double vol = Box[0]*Box[1]*Box[2];
  // Now, calculate charge density
  FILE *fout = fopen ("rho.dat", "w");
  for (double x=-0.5*Box[0]; x<=0.5*Box[0]; x+=0.02) {
    Vec3 r(x, 0.0, 0.0);
    complex<double> psi(0.0, 0.0);
    for (int i=0; i<GVecs.size(); i++) {
      double phase = dot (r, GVecs(i));
      complex<double> ex (cos(phase), sin(phase));
      psi += phik(i) * ex;
    }
    double rho = (psi.real()*psi.real() + psi.imag()*psi.imag())/vol;
    fprintf (fout, "%1.12e %1.12e\n", x, rho);
  }
  fclose (fout);
}



void Test (IOSectionClass &in)
{
  string PHname;
  assert (in.ReadVar ("PHfile", PHname));
  IOSectionClass PHin;
  assert (PHin.OpenFile(PHname));
  Potential *PH = ReadPotential(PHin);
  PHin.CloseFile();

  Array<double,1> box;
  assert(in.ReadVar ("Box", box));
  assert(box.size() == 3);

  double kCut;
  assert (in.ReadVar ("kCut", kCut));
  
  kSpaceTest kTest (*PH);
  kTest.SetBox (box(0), box(1), box(2));
  kTest.SetupGVecs (kCut);
  Vec3 k(0.0, 0.0, 0.0);
  kTest.CalcHamiltonian(k);

  int N = kTest.H.rows();
//   cerr << "H = \n";
//   for (int i=0; i<N; i++) {
//     for (int j=0; j<N; j++)
//       fprintf (stderr, "%11.4e ", kTest.H(i,j));
//     fprintf (stderr, "\n");
//   }
  Array<double,1> vals(N);
  Array<double,2> vecs(1,N);
  SymmEigenPairs (kTest.H, 1, vals, vecs);
  cerr << "Energy = " << vals(0) << endl;
  kTest.PlotRho (vecs(0,Range::all()));
}



main(int argc, char **argv)
{
  if (argc < 2) {
    cerr << "Usage:\n"
	 << "TestkSpacePH file.in\n";
  }
  else {
    IOSectionClass in;
    assert (in.OpenFile (argv[1]));
    Test (in);
  }
}
