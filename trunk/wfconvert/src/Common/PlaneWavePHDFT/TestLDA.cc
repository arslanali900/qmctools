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

#include "../MPI/Communication.h"
#include "PlaneWavesMPI.h"
#include "../IO/IO.h"

void TestInitCharge()
{
  CommunicatorClass bandComm, kComm;
  bandComm.SetWorld();
  Array<int,1> root(1);
  root(0) = 0;
  bandComm.Subset (root, kComm);

  IOSectionClass in;
  in.OpenFile("NaLocalPH.h5");
  //in.OpenFile("Na_HF_NLPP.h5");
  Potential *V_elec_ion = ReadPotential(in);
  CoulombPot V_ion_ion;
  V_ion_ion.Z1Z2 = 1.0;
  in.CloseFile();

  int numBands  = 10;
  int numElecs  = 16;
  Vec3 box (26.56, 26.56, 26.56);
  Array<Vec3,1> rions(16);
  rions(0)  = Vec3(  6.7701,   -0.5241,   -0.7601);
  rions(1)  = Vec3( -0.2356,  -11.1819,    0.1041);
  rions(2)  = Vec3(  0.7953,   -3.1228,   -3.3658);
  rions(3)  = Vec3( -9.6992,    2.5853,   -8.7713);
  rions(4)  = Vec3( -1.0770,  -12.9317,  -10.9011);
  rions(5)  = Vec3( -8.4076,   -4.4247,   -6.0560);
  rions(6)  = Vec3(-11.7379,    0.7186,    6.3867);
  rions(7)  = Vec3(  0.7048,    6.5211,   -7.6078);
  rions(8)  = Vec3( -0.9636,   -1.0799,    3.2954);
  rions(9)  = Vec3( 13.1805,   -9.5798,   -8.3901);
  rions(10) = Vec3( -1.8638,   -4.0041,  -11.5047);
  rions(11) = Vec3(  0.4408,    7.2936,    3.7794);
  rions(12) = Vec3(  8.0325,   -6.3283,    4.1347);
  rions(13) = Vec3(  5.4652,    2.5650,    5.4588);
  rions(14) = Vec3( -7.7732,   -6.0467,    1.7018);
  rions(15) = Vec3(  5.5497,  -11.6650,    9.0145);
  Vec3 Gprim = 2.0*M_PI*Vec3(1.0/box[0], 1.0/box[1], 1.0/box[2]);
  Vec3 k = 0.25 * Gprim;
  MPISystemClass system (numBands, numElecs, bandComm, kComm, true, false);
  
  system.Setup (box, k, 4.0, *V_elec_ion, V_ion_ion, true, true);
  system.SetIons(rions);
  IOSectionClass out;
  out.NewFile ("InitDensity.h5");
  out.WriteVar("InitRho", system.GetDensity());
  out.CloseFile();
}


void TestSolveLDA()
{
  CommunicatorClass bandComm, kComm;
  bandComm.SetWorld();
  Array<int,1> root(1);
  root(0) = 0;
  bandComm.Subset (root, kComm);

  IOSectionClass in;
  //in.OpenFile("NaLocalPH.h5");
  in.OpenFile("Na_HF_NLPP.h5");
  Potential *V_elec_ion = ReadPotential(in);
  CoulombPot V_ion_ion;
  V_ion_ion.Z1Z2 = 1.0;
  in.CloseFile();

  int numBands  = 16;
  int numElecs  = 16;
  Vec3 box (26.56, 26.56, 26.56);
  Array<Vec3,1> rions(16);
  rions(0)  = Vec3(  6.7701,   -0.5241,   -0.7601);
  rions(1)  = Vec3( -0.2356,  -11.1819,    0.1041);
  rions(2)  = Vec3(  0.7953,   -3.1228,   -3.3658);
  rions(3)  = Vec3( -9.6992,    2.5853,   -8.7713);
  rions(4)  = Vec3( -1.0770,  -12.9317,  -10.9011);
  rions(5)  = Vec3( -8.4076,   -4.4247,   -6.0560);
  rions(6)  = Vec3(-11.7379,    0.7186,    6.3867);
  rions(7)  = Vec3(  0.7048,    6.5211,   -7.6078);
  rions(8)  = Vec3( -0.9636,   -1.0799,    3.2954);
  rions(9)  = Vec3( 13.1805,   -9.5798,   -8.3901);
  rions(10) = Vec3( -1.8638,   -4.0041,  -11.5047);
  rions(11) = Vec3(  0.4408,    7.2936,    3.7794);
  rions(12) = Vec3(  8.0325,   -6.3283,    4.1347);
  rions(13) = Vec3(  5.4652,    2.5650,    5.4588);
  rions(14) = Vec3( -7.7732,   -6.0467,    1.7018);
  rions(15) = Vec3(  5.5497,  -11.6650,    9.0145);
  Vec3 Gprim = 2.0*M_PI*Vec3(1.0/box[0], 1.0/box[1], 1.0/box[2]);
  Vec3 k = 0.25 * Gprim;
  cerr << "k = " << k << endl;
  MPISystemClass system (numBands, numElecs, bandComm, kComm, true, false);
  
  system.Setup (box, k, 3.0, *V_elec_ion, V_ion_ion, true, true);
  system.SetIons(rions);
  system.SolveLDA();
}


void TestMultiLDA()
{
  CommunicatorClass bandComm, kComm;
  bandComm.SetWorld();
  Array<int,1> root(1);
  root(0) = 0;
  bandComm.Subset (root, kComm);

  IOSectionClass in;
  //in.OpenFile("NaLocalPH.h5");
  //  in.OpenFile("Na_HF_NLPP.h5");
  in.OpenFile("OpiumNaLocal.h5");
  Potential *V_elec_ion = ReadPotential(in);
  in.CloseFile();
  CoulombPot V_ion_ion;
  V_ion_ion.Z1Z2 = 1.0;



  int numBands  = 16;
  int numElecs  = 16;
  Vec3 box (26.56, 26.56, 26.56);
  
  Vec3 Gprim = 2.0*M_PI*Vec3(1.0/box[0], 1.0/box[1], 1.0/box[2]);
  Vec3 k = 0.25 * Gprim;
  cerr << "k = " << k << endl;
  MPISystemClass system (numBands, numElecs, bandComm, kComm, true, false);
  
  system.Setup (box, k, 4.0, *V_elec_ion, V_ion_ion, true, true);
  
  Array<Vec3,1> rions(16);
  Array<double,3> R;
  IOSectionClass configsIn;
  configsIn.OpenFile("configs.h5");
  configsIn.ReadVar("R", R);
  FILE *fout;
  if (bandComm.MyProc() == 0)
    fout = fopen ("Energies1.dat", "w");
  for (int conf=0; conf<1/*R.extent(0)*/; conf++) {
    for (int ri=0; ri<R.extent(1); ri++)
      for (int dim=0; dim<3; dim++)
	rions(ri)[dim] = R(conf,ri,dim);
    system.SetIons(rions);
    system.DoMDExtrap();
    system.SolveLDA();
    if (bandComm.MyProc() == 0) {
      for (int bi=0; bi<numBands; bi++)
	fprintf (fout, "%1.12e ", system.GetEnergy(bi));
      fprintf (fout, "\n");
      fflush (fout);
    }
  }
  




}


void TestLDAForces()
{
  CommunicatorClass bandComm, kComm;
  bandComm.SetWorld();
  Array<int,1> root(1);
  root(0) = 0;
  bandComm.Subset (root, kComm);

  IOSectionClass in;
  //in.OpenFile("NaLocalPH.h5");
  //  in.OpenFile("Na_HF_NLPP.h5");
  in.OpenFile("OpiumNaLocal.h5");
  Potential *V_elec_ion = ReadPotential(in);
  in.CloseFile();
  CoulombPot V_ion_ion;
  V_ion_ion.Z1Z2 = 1.0;

  int numBands  = 16;
  int numElecs  = 16;
  Vec3 box (26.56, 26.56, 26.56);
  
  Vec3 Gprim = 2.0*M_PI*Vec3(1.0/box[0], 1.0/box[1], 1.0/box[2]);
  Vec3 k = 0.25 * Gprim;
  //Vec3 k = 0.0;
  cerr << "k = " << k << endl;
  MPISystemClass system (numBands, numElecs, bandComm, kComm, true, false);
  
  system.Setup (box, k, 4.0, *V_elec_ion, V_ion_ion, true, true);
  
  Array<Vec3,1> rions(16);
  Array<double,3> R;
  IOSectionClass configsIn;
  configsIn.OpenFile("configs.h5");
  configsIn.ReadVar("R", R);
  for (int ri=0; ri<R.extent(1); ri++)
    for (int dim=0; dim<3; dim++)
      rions(ri)[dim] = R(0,ri,dim);
  
  system.SetIons(rions);
  system.DoMDExtrap();
  system.SolveLDA();
  
  Array<Vec3,1> forcesFD(16);
  if (bandComm.MyProc() == 0)
    fprintf (stderr, "Finite difference forces:\n");
  for (int pi=0; pi<16; pi++) {
    for (int dim=0; dim<3; dim++) {
      rions(pi)[dim] += 1.0e-4;
      system.SetIons(rions);
      double Eplus = system.CalcLocalPPEnergy()+system.EwaldEnergy();
      rions(pi)[dim] -= 2.0e-4;
      system.SetIons(rions);
      double Eminus = system.CalcLocalPPEnergy()+system.EwaldEnergy();
      rions(pi)[dim] += 1.0e-4;
      forcesFD(pi)[dim] = -(Eplus-Eminus)/2.0e-4;
      if (bandComm.MyProc() == 0)
	fprintf (stderr, "%18.12f ", forcesFD(pi)[dim]);
    }
    if (bandComm.MyProc() == 0)
      fprintf (stderr, "\n");
  }
  if (bandComm.MyProc() == 0)
    fprintf (stderr, "\nAnalytic forces:\n");
  Array<Vec3,1> forces(16);
  system.CalcIonForces(forces);
  if (bandComm.MyProc() == 0)
    for (int pi=0; pi<16; pi++) 
      fprintf(stderr, "%18.12f %18.12f %18.12f\n",
	      forces(pi)[0], forces(pi)[1], forces(pi)[2]);
}

void TestSolidLDA()
{
  CommunicatorClass bandComm, kComm;
  bandComm.SetWorld();
  Array<int,1> root(1);
  root(0) = 0;
  bandComm.Subset (root, kComm);

  IOSectionClass in;
  //in.OpenFile("NaLocalPH.h5");
  in.OpenFile("Na_HF_NLPP.h5");
  Potential *V_elec_ion = ReadPotential(in);
  in.CloseFile();
  CoulombPot V_ion_ion;
  V_ion_ion.Z1Z2 = 1.0;


  int numBands  = 5;
  int numElecs  = 2;
  Vec3 box (8.11, 8.11, 8.11);
  Array<Vec3,1> rions(2);
  rions(0)  = Vec3( 0.00, 0.00, 0.00); 
  rions(1)  = Vec3( 4.05, 4.05, 4.05);
  Vec3 Gprim = 2.0*M_PI*Vec3(1.0/box[0], 1.0/box[1], 1.0/box[2]);
  Vec3 k = 0.25 * Gprim;
  cerr << "k = " << k << endl;
  MPISystemClass system (numBands, numElecs, bandComm, kComm, true, false);
  
  system.Setup (box, k, 8.0, *V_elec_ion, V_ion_ion, true, true);
  system.SetIons(rions);
  system.SolveLDA();
}


#include "FermiSmear.h"
void 
TestSmear()
{
  FILE *fout = fopen ("smear.dat", "w");
  double mu = 0.0;
  MethfesselPaxton smearer;
  smearer.SetOrder(2);
  smearer.SetWidth(0.5);
  for (double E=-3.0; E<=3.0; E+=0.001) {
    fprintf (fout, "%1.12e ", E);
    for (int order=1; order<6; order++) {
      smearer.SetOrder(order);
      fprintf (fout, "%1.12e ", smearer.S(E,mu));
    }
    fprintf (fout, "\n");
  }
  fclose(fout);
}	     


main(int argc, char **argv)
{
  COMM::Init(argc, argv);
  // TestInitCharge();
  //TestSmear();
  //TestSolveLDA();
  //TestMultiLDA();
  // TestSolidLDA();
  TestLDAForces();
  COMM::Finalize();
}
