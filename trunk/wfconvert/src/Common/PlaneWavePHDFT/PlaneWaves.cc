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

#include "PlaneWaves.h"

void
SystemClass::Setup (Vec3 box, Vec3 k, double kcut, Potential &ph,
		    bool useFFT)
{
  Box = box;
  kCut = kcut;
  GVecs.Set (box, k, kcut);
  FFT.Setup();
  H.SetIonPot (ph, useFFT);
  H.Setk(k);
  Bands.resize (NumBands, GVecs.size());
  PH = &ph;
  UseFFT = useFFT;
}


void 
SystemClass::Setup (Vec3 box, Vec3 k, double kcut, double z,
		    bool useFFT)
{
  Box = box;
  GVecs.Set (box, k, kcut);
  FFT.Setup();
  H.SetIonPot (z, useFFT);
  H.Setk(k);
  Bands.resize (NumBands, GVecs.size());
  UseFFT = useFFT;
}

void
SystemClass::SetIons (const Array<Vec3,1> &rions)
{
  if (Rions.size() != rions.size())
    Rions.resize(rions.size());
  /// This half box addition compensates for the fourier transform
  /// aliasing. 
  Rions = rions;
  H.SetIons(rions+0.5*Box);
}


void 
SystemClass::DiagonalizeH ()
{
  CG.InitBands();
  for (int i=0; i<Bands.rows(); i++) {
    CG.Solve(i);
    fprintf (stderr, "Energy(%d) = %15.12f\n", i, CG.Energies(i)* 27.211383);
  }
}


void
SystemClass::CalcChargeDensity(Array<double,3> &rho)
{
  int nx, ny, nz;
  GVecs.GetFFTBoxSize(nx, ny, nz);
  rho.resize(nx, ny, nz);
  rho = 0.0;
  zVec c;
  for (int band=0; band<Bands.rows(); band++) {
    c.reference (Bands(band,Range::all()));
    FFT.PutkVec (c);
    FFT.k2r();
    for (int ix=0; ix<nx; ix++)
      for (int iy=0; iy<ny; iy++)
	for (int iz=0; iz<nz; iz++) {
	  complex<double> z = FFT.rBox(ix,iy,iz);
	  rho(ix,iy,iz) += z.real()*z.real()+z.imag()*z.imag();
	}
  }    

}

void 
SystemClass::WriteXSFFile (string filename)
{
  const double bohr2angst = 0.52917721;
  FILE *fout = fopen(filename.c_str(), "w");
  assert (fout !=NULL );

  // First, write out box geometry
  fprintf (fout, "CRYSTAL\n");
  fprintf (fout, "PRIMVEC\n");
  for (int i=0; i<3; i++) {
    for (int j=0; j<3; j++) {
      double val = (i==j) ? Box[i] : 0.0;
      fprintf (fout, "%15.10f ", val*bohr2angst);
    }
    fprintf (fout, "\n");
  }
  fprintf (fout, "CONVVEC\n");
  for (int i=0; i<3; i++) {
    for (int j=0; j<3; j++) {
      double val = (i==j) ? Box[i] : 0.0;
      fprintf (fout, "%15.10f ", val*bohr2angst);
    }
    fprintf (fout, "\n");
  }
 
  // Now write atom positions
  fprintf (fout, "PRIMCOORD\n");
  fprintf (fout, "%d 1\n", Rions.size());
  for (int i=0; i<Rions.size(); i++) {
    Vec3 pos = Rions(i);
    for (int j=0; j<3; j++)
      while (pos[j] < 0.0)
	pos[j] += Box[j];
    fprintf (fout, "%d %15.12f %15.12f %15.12f\n", 11,
	     pos[0]*bohr2angst, pos[1]*bohr2angst, pos[2]*bohr2angst);
  }
  
  // Now write out the charge density
  Array<double,3> rho;
  CalcChargeDensity(rho);
  fprintf (fout, "BEGIN_BLOCK_DATAGRID_3D\n");
  fprintf (fout, "ChargeDensity\n");
  fprintf (fout, "BEGIN_DATAGRID_3D_RHO\n");
  int nx = rho.extent(0);
  int ny = rho.extent(1);
  int nz = rho.extent(2);

  fprintf (fout, "%d %d %d\n", nx+1, ny+1, nz+1);
  /// Write origin
  fprintf (fout, "0.0 0.0 0.0\n");

  /// Write spanning vectors
  for (int i=0; i<3; i++) {
    for (int j=0; j<3; j++) {
      double val = (i==j) ? Box[i] : 0.0;
      fprintf (fout, "%15.10f ", val*bohr2angst);
    }
    fprintf (fout, "\n");
  }

  int n = 0;
  for (int k=0; k<=nz; k++)
    for (int j=0; j<=ny; j++)
      for (int i=0; i<=nx; i++) {
	fprintf (fout, "%10.8f ", rho(i%nx,j%ny,k%nz));
	n++;
	if ((n%4) == 0)
	  fprintf (fout, "\n");
      }
  fprintf (fout, "\n");
  fprintf (fout, "END_DATAGRID_3D\n");
  fprintf (fout, "END_BLOCK_DATAGRID_3D\n");
  fclose(fout);
}

	   
	
void 
SystemClass::SetRealSpaceBandNum(int num)
{
  zVec c;
  c.reference (Bands(num,Range::all()));
  FFT.PutkVec (c);
  FFT.k2r();
}


void 
SystemClass::Setk(Vec3 k)
{
  GVecs.Set (Box, k, kCut);
  FFT.Setup();
  //  H.SetIonPot (*PH, UseFFT);
  H.Setk(k);
  /// The half box is necessary to compensate for fourier transform
  /// aliasing.  
  H.SetIons(Rions+0.5*Box);
  Bands.resize (NumBands, GVecs.size());
  CG.Setup();
}
