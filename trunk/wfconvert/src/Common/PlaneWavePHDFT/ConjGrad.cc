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

#include "ConjGrad.h"


void ConjGrad::Setup()
{
  int N = H.GVecs.size();
  Bands.resize (NumBands, N);
  Energies.resize(NumBands);
  cnext.resize (N);
  Hc.resize (N);
  Phi.resize (N);
  Phip.resize(N);
  Phipp.resize(N);
  Xi.resize(N);
  Eta.resize(N);
  Bands = 0.0;
  for (int band=0; band<NumBands; band++) {
    c.reference(Bands(band,Range::all()));
    c = 0.0;
    for (int i=0; i<N; i++) {
      if (dot (H.GVecs(i), H.GVecs(i)) < 1.0)
	c(i) = 1.0;
    }
    Normalize(c);
  }
//   Orthogonalize (Bands);
//   for (int band=0; band<NumBands; band++) {
//     c.reference(Bands(band,Range::all()));
//     Normalize (c);
//   }

  PrintOverlaps();
  IsSetup = true;
}

void ConjGrad::CalcPhiSD()
{
  H.Apply(c, Hc);
  E0 = realconjdot (c, Hc);
  cerr << "E = " << E0 << endl;
  Phip = E0*c - Hc;
  Normalize (Phip);  
}

void ConjGrad::Precondition()
{
  double Tinv = 1.0/T;
  for (int i=0; i<c.size(); i++) {
    double x = 0.5*dot(H.GVecs(i), H.GVecs(i))*Tinv;
    double num = 27.0 + 18.0*x +12.0*x*x + 8.0*x*x*x;
    double denom = num + 16.0*x*x*x*x;
    Eta(i) = (num/denom)* Xi(i);
  }
}

void ConjGrad::CalcPhiCG()
{
  if (LastBand != CurrentBand) {
    LastBand = CurrentBand;
    EtaXiLast = 0.0;
    Orthogonalize2(Bands, c, CurrentBand-1);
    Normalize(c);
  }
  Hc = 0.0;
  H.Kinetic.Apply (c, Hc);
  T = realconjdot (c, Hc);
  //H.Coulomb.Apply (c, Hc);
  H.CoulombFFT.Apply (c, Hc);
  //H.PH.Apply (c, Hc);
  //H.PHFFT.Apply (c,Hc);
  E0 = realconjdot (c, Hc);
  Xi = E0*c - Hc;
  /// Orthonalize to other bands here
  //Xip = Xi;
  zVec &Xip = Xi;
  Orthogonalize2 (Bands, Xip, CurrentBand);
  //Orthogonalize (Bands, Xip);
  //  CheckOrthog (Bands, Xip);

  Precondition();
  //Eta = Xip;

  // Now, orthogonalize to psi
  // rename for clarity
  zVec &Etap = Eta;
  Orthogonalize2 (Bands, Etap, CurrentBand);
  //Orthogonalize (Bands, Etap);
  //  CheckOrthog (Bands, Etap);
  // Etap = Eta - conjdot (c, Eta)*c;
  complex<double> etaxi = conjdot(Etap, Xip);
  complex<double> gamma; 
  if (EtaXiLast != complex<double>(0.0, 0.0)) 
    gamma = etaxi/EtaXiLast;
  else
    gamma = 0.0;
  EtaXiLast = etaxi;
  
  Phi = Etap + gamma * Phi;
  
  Phipp = Phi - conjdot(c, Phi) * c;
  Phip = Phipp;
  Normalize (Phip);

  Energies(CurrentBand) = E0;
}


void ConjGrad::Iterate(int band)
{
  //  PrintOverlaps();
  CurrentBand = band;
  if (!IsSetup)
    Setup();
  c.reference (Bands(band,Range::all()));
  // First, calculate steepest descent vector:
  CalcPhiCG();
  //CalcPhiSD();

  // Now, pick optimal theta for 
  // cnext = cos(theta)*c + sin(theta)*Phi;
  double dE_dtheta = 2.0*realconjdot(Phip, Hc);
  double theta1 = M_PI/300.0;
  cnext = cos(theta1)*c + sin(theta1)*Phip;
  H.Apply (cnext, Hc);
  double E1 = realconjdot (cnext, Hc);
  double A1 = (E0 - E1 + 0.5*sin(2.0*theta1)*dE_dtheta)/(1.0-cos(2.0*theta1));
  double B1 = 0.5*dE_dtheta;
  double thetaMin = 0.5*atan (B1/A1);
  //  cerr << "thetaMin = " << thetaMin << endl;

  c = cos(thetaMin)*c + sin(thetaMin)*Phip;
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
