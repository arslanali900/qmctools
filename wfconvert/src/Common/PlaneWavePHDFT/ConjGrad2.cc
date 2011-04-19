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

#include "ConjGrad2.h"
#include "../MatrixOps/MatrixOps.h"

void ConjGrad::Setup()
{
  int numBands = Bands.rows();
  int N = H.GVecs.size();
  Bands.resize (numBands, N);
  Energies.resize(numBands);
  cnext.resize (N);
  Hc.resize (N);
  Phi.resize (N);
  Phip.resize(N);
  Phipp.resize(N);
  Xi.resize(N);
  Eta.resize(N);
  Bands = 0.0;
  Vec3 kBox = H.GVecs.GetkBox();
  double maxk = max(kBox[0], max(kBox[1], kBox[2]));
  double maxE = 1.0*maxk*maxk;
  int numNonZero;
  for (int band=0; band<numBands; band++) {
    numNonZero = 0;
    c.reference(Bands(band,Range::all()));
    c = 0.0;
    for (int i=0; i<N; i++) {
      double dist2 = dot(H.GVecs(i), H.GVecs(i));
      c(i) = exp(-2.0*dist2);
//       if (dot (H.GVecs(i), H.GVecs(i)) < maxE) {
// 	c(i) = 1.0;
// 	numNonZero++;
//       }
    }
    Normalize(c);
  }
  cerr << "numNonZero = " << numNonZero << endl;
  InitBands();
//   Orthogonalize (Bands);
//   for (int band=0; band<numBands; band++) {
//     c.reference(Bands(band,Range::all()));
//     Normalize (c);
//   }

//  PrintOverlaps();
  IsSetup = true;
}


void ConjGrad::InitBands()
{
  int numBands = Bands.rows();
  //  int numVecs = 6 * numBands;
  int numVecs = 4 * numBands;
  assert (numVecs <= H.GVecs.size());
  Array<complex<double>,2> Hmat(numVecs, numVecs);
  Array<complex<double>,2> EigVecs(numBands, numVecs);
  Array<double,1> EigVals (numBands);

  H.Vion->Vmatrix(Hmat);
  cerr << "kPoint = " << H.kPoint << endl;
  for (int i=0; i<numVecs; i++) {
    Vec3 Gpk = H.GVecs(i) + H.kPoint;
    Hmat(i,i) += 0.5 * dot (Gpk,Gpk);
  }
  /// DEBUG DEBUG DEBUG
//   CommunicatorClass comm;
//   comm.SetWorld();
//   int myProc = comm.MyProc();
//   char fname[100];
//   snprintf (fname, 100, "Hmat%d.dat", myProc);

//   FILE *fout = fopen (fname, "w");
//   for (int i=0; i<numVecs; i++) {
//     for (int j=0; j<numVecs; j++) 
//       fprintf (fout, "%1.16e %1.16e ", Hmat(i,j).real(), Hmat(i,j).imag());
//     fprintf (fout, "\n");
//   }
//   fclose(fout);
  // Now diagonalize
  SymmEigenPairs (Hmat, numBands, EigVals, EigVecs);

//   snprintf (fname, 100, "InitBand%d.dat", myProc);
//   fout = fopen (fname, "w");
//   for (int i=0; i<numVecs; i++) {
//     for (int band=0; band<numBands; band++) 
//       fprintf (fout, "%1.16e ", EigVecs(band, i));
//     fprintf (fout, "\n");
//   }
//   fclose (fout);
  
  for (int i=0; i<numBands; i++)
    cerr << "Mini energy(" << i << ") = " << 27.211383*EigVals(i) << endl;

  // Now put results in Bands
  Bands = 0.0;
  for (int band=0; band<numBands; band++) 
    for (int i=0; i<numVecs; i++)
      Bands(band, i) = EigVecs(band, i);
}

double ConjGrad::CalcPhiSD()
{
  H.Apply(c, Hc);
  E0 = realconjdot (c, Hc);
  cerr << "E = " << E0 << endl;
  Phip = E0*c - Hc;
  double nrm = norm (Phip);
  Normalize (Phip);  
  return nrm;
}

void ConjGrad::Precondition()
{
  double Tinv = 1.0/T;
  Vec3 k = H.kPoint;
  for (int i=0; i<c.size(); i++) {
    double x = 0.5*dot(H.GVecs(i)+k, H.GVecs(i)+k)*Tinv;
    double num = 27.0 + 18.0*x +12.0*x*x + 8.0*x*x*x;
    double denom = num + 16.0*x*x*x*x;
    Eta(i) = (num/denom)* Xi(i);
  }
}


// Returns the norm of the residual Hc - Ec
double ConjGrad::CalcPhiCG()
{
  Hc = 0.0;
  H.Kinetic.Apply (c, Hc);
  T = realconjdot (c, Hc);
  H.Vion->Apply (c, Hc);
  E0 = realconjdot (c, Hc);
  Xi = E0*c - Hc;
  double residualNorm = norm (Xi);
  /// Orthonalize to other bands here
  zVec &Xip = Xi;
  Orthogonalize2 (Bands, Xip, CurrentBand);
  Precondition();

  // Now, orthogonalize to all bands
  // rename for clarity
  zVec &Etap = Eta;
  Orthogonalize2 (Bands, Etap, CurrentBand);
  // Etap = Eta - conjdot (c, Eta)*c;
  // Compute conjugate direction
  complex<double> etaxi = conjdot(Etap, Xip);
  complex<double> gamma; 
  if (EtaXiLast != complex<double>(0.0, 0.0)) 
    gamma = etaxi/EtaXiLast;
  else
    gamma = 0.0;
  EtaXiLast = etaxi;
  
  Phi = Etap + gamma * Phi;
  
  // Orthogonalize to present band
  complex<double> cPhi = conjdot(c,Phi);
  Phip = Phi - cPhi * c;
  Normalize (Phip);
  Energies(CurrentBand) = E0;
  return residualNorm;
}


void ConjGrad::Solve(int band)
{
  CurrentBand = band;
  if (!IsSetup)
    Setup();
  c.reference (Bands(band,Range::all()));
  EtaXiLast = 0.0;
  Orthogonalize2 (Bands, c, band-1);
  Normalize(c);

  double Elast = 1.0e100;
  double residualNorm = 1.0;
  //  while (fabs (Elast - Energies(band)) > Tolerance) {
  int iter=0;
  while ((residualNorm > 1.0e-6) && (iter < 50)) {
    //    cerr << "Energy = " << 27.211383*Energies(band) << endl;
    Elast = Energies(band);
    // First, calculate conjugate gradient direction
    residualNorm = CalcPhiCG();
    
    // Now, pick optimal theta for 
    double dE_dtheta = 2.0*realconjdot(Phip, Hc);

    H.Apply (Phip, Hc);
    double d2E_dtheta2 = 2.0*(realconjdot(Phip, Hc) - E0);
    double thetaMin = 0.5*atan(-dE_dtheta/(0.5*d2E_dtheta2));

//     double theta1 = M_PI/300.0;
//     cnext = cos(theta1)*c + sin(theta1)*Phip;
//     H.Apply (cnext, Hc);
//     double E1 = realconjdot (cnext, Hc);
//     double A1=(E0 - E1 + 0.5*sin(2.0*theta1)*dE_dtheta)/(1.0-cos(2.0*theta1));
//     double B1 = 0.5*dE_dtheta;
//     double thetaMin = 0.5*atan (B1/A1);
    double costhetaMin, sinthetaMin;
    sincos(thetaMin, &sinthetaMin, &costhetaMin);
    c = costhetaMin*c + sinthetaMin*Phip;
    iter++;
  }
  if (residualNorm > 1.0e-8)
    cerr << "Warning:  conjugate gradient residual norm = " 
	 << residualNorm << endl;
  cerr << "# of iterations = " << iter << endl;
}

void
ConjGrad::PrintOverlaps()
{
  cerr << "Overlaps = \n";
  zVec x, y;
  for (int i=0; i<Bands.rows(); i++) {
    x.reference (Bands(i, Range::all()));
    for (int j=0; j<Bands.rows(); j++) {
      y.reference (Bands(j,Range::all()));
      double s = realconjdot(x,y);
      fprintf (stderr, "%12.6e ", s);
    }
    fprintf (stderr, "\n");
  }
}
