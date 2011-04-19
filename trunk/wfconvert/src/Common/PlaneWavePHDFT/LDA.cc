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
#include "../MatrixOps/MatrixOps.h"
#include "../DFT/Functionals.h"
#include "../Ewald/OptimizedBreakup.h"
#include "../Integration/GKIntegration.h"
#include <algorithm>

void
MPISystemClass::InitLDA()
{
  ChargeMixer = new KerkerMixerClass (FFT);
  int NDelta = GVecs.DeltaSize();
  FFT.GetDims    (Nx, Ny, Nz);
  NewRho.resize  (Nx, Ny, Nz);
  Rho_r1.resize  (Nx, Ny, Nz);
  Rho_r2.resize  (Nx, Ny, Nz);
  TempRho.resize (Nx, Ny, Nz);
  VH.resize      (Nx, Ny, Nz);
  VXC.resize     (Nx, Ny, Nz);
  VHXC.resize    (Nx, Ny, Nz);
  Rho_r.resize   (Nx, Ny, Nz);
  h_G.resize     (NDelta);
  Rho_G.resize   (NDelta);
  Hmat.resize(NumBands, NumBands);
  EigVecs.resize(NumBands, NumBands);
  EigVals.resize(NumBands);
  RotBands.resize(NumBands, GVecs.size());
  Occupancies.resize(NumBands);
  CalcRadialChargeDensity();
  Smearer.SetOrder(1);
  Smearer.SetWidth(0.001);
  DoOptimizedBreakup();
  if (BandComm.MyProc() == 0) {
    Numk = kComm.NumProcs();
    Myk  = kComm.MyProc();
  }
  BandComm.Broadcast (0, Numk);
  BandComm.Broadcast (0, Myk);
}

#include "Hamiltonians.h"


void
MPISystemClass::SubspaceRotate()
{
  //  ApplyH();
  zVec c, Hc;
  for (int bi=0; bi<Bands.rows(); bi++) {
    c.reference (Bands(bi,Range::all()));
    for (int bj=0; bj<Bands.rows(); bj++) {
      Hc.reference (HBands(bj, Range::all()));
      Hmat(bi,bj) = conjdot (c, Hc);
    }
  }
  SymmEigenPairs(Hmat, Hmat.rows(), EigVals, EigVecs);
  // Now, compute the new rotated bands.
  RotBands = 0.0;
  for (int bi=0; bi<NumBands; bi++)
    for (int bj=0; bj<NumBands; bj++)
      RotBands(bi,Range::all()) += EigVecs(bi,bj)*Bands(bj,Range::all());
  
  Bands = RotBands;
  // And compute rotated HBands
  RotBands = 0.0;
  for (int bi=0; bi<NumBands; bi++)
    for (int bj=0; bj<NumBands; bj++)
      RotBands(bi,Range::all()) += EigVecs(bi,bj)*HBands(bj,Range::all());
  HBands = RotBands;

  zVec Tc(Bands.extent(1));
  for (int bi=0; bi<NumBands; bi++) {
    Tc = 0.0;
    c.reference (Bands(bi, Range::all()));
    H.Kinetic.Apply (c, Tc);
    CG.T(bi) = realconjdot (c, Tc);
  }
  // HACK HACK HACK
  // CG.ApplyH();

//   for (int bi=0; bi<NumBands; bi++) {
//     c.reference (Bands(bi,Range::all()));
//     Hc.reference (HBands(bi, Range::all()));
//     double E = realconjdot(c, Hc);
//     perr << "EigVal(bi) = " << EigVals(bi) << "cHc = " << E << endl;
//   }
}


void
MPISystemClass::SolveLDA()
{
  if (UseSubspaceRotation)
    CG.SetOrthoType (ORTHO_ALL);
  else
    CG.SetOrthoType (ORTHO_LOWER);

  // Make sure we have enough bands
  assert (NumBands >= (NumElecs+1)/2);
  SCiter = 0;
  bool SC = false;
  int highestOcc = (NumElecs+1)/2-1;
  double lastE = 500.0;
  ChargeMixer->Reset();
  
  double startLambda = ChargeMixer->GetLambda();

  while (!SC) {
    MixChargeDensity();
    CalcVHXC();
    CG.SetTolerance(1.0e-4);
    if (SCiter == 0 && ConfigNum == 0) {
      CG.InitBands();
      CG.ApplyH();
      SubspaceRotate();
    }
    CG.ApplyH();
    if (UseSubspaceRotation) 
      SubspaceRotate();

    CG.Solve();
    SC = (fabs(CG.Energies(highestOcc)-lastE) < 1.0e-8);
    /// Now make sure we all k-points agree that we're
    /// self-consistent.
    if (BandComm.MyProc() == 0) 
      kComm.AllAnd(SC);
    BandComm.Broadcast(0, SC);

    lastE = CG.Energies(highestOcc);
    CalcOccupancies();
    CalcChargeDensity();
    SCiter++;
    /// Stability control:  if convergence is taking forever, we're
    /// probably  charge sloshing. Then increase the Kerker
    /// parameter
    if ((SCiter % 50) == 0)
      ChargeMixer->SetLambda (1.2*ChargeMixer->GetLambda());
  }
  // Reset initial lambda
  ChargeMixer->SetLambda(startLambda);

//   if (BandComm.MyProc() == 0) {
//     if (kComm.MyProc() == 0) {
//       IOSectionClass out;
//       char fname[100];
//       snprintf (fname, 100, "VHXC%d.h5", ConfigNum);
//       out.NewFile (fname);
//       out.WriteVar("VH", VH);
//       out.WriteVar("VXC", VXC);
//       out.WriteVar("VHXC", VHXC);
//       out.WriteVar("Rho", Rho_r);
//       out.WriteVar ("E", CG.Energies);
//       LocalPotFFTClass &pot = *((LocalPotFFTClass*)H.Vion);
//       Array<double,3> magV(VH.shape());
//       for (int ix=0; ix<magV.extent(0); ix++)
// 	for (int iy=0; iy<magV.extent(1); iy++)
// 	  for (int iz=0; iz<magV.extent(2); iz++)
// 	    magV(ix,iy,iz) = mag(pot.Vr(ix,iy,iz));
//       out.WriteVar("Vion", magV);
//       out.CloseFile();
//     }
//   }
  Bands2 = Bands1;
  Rions2 = Rions1;
  Bands1 = Bands;
  Rions1 = Rions;
  Rho_r2 = Rho_r1;
  Rho_r1 = Rho_r;
  ConfigNum++;
}


double
Residue (Array<Vec3,1> &R0, Array<Vec3,1> &R1, Array<Vec3,1> &R2,
	 double alpha)
{
  double residue = 0;
  for (int i=0; i<R0.size(); i++)
    for (int dim=0; dim < 3; dim++) {
      double diff = R0(i)[dim] -(alpha*R1(i)[dim] + (1.0-alpha)*R2(i)[dim]);
      residue += diff*diff;
    }
  return residue;
}

/// See Dario Alfe, Computer Physics Communications 118, 31-33 (1999)
void
MPISystemClass::ExtrapolateCharge()
{


}

/// See Arias et al. Phys. Rev. B 45, 1538 (1992)
void 
MPISystemClass::DoMDExtrap()
{
  if (ConfigNum > 1) {
    double e0, e1, e2;
    e0  = Residue(Rions, Rions1, Rions2,  0.0);
    e1  = Residue(Rions, Rions1, Rions2,  1.0);
    e2  = Residue(Rions, Rions1, Rions2,  2.0);
    
    double a = 0.5*(e2-2.0*e1+e0);
    double b = 0.5*(4.0*e1 - e2 - 3.0*e0);
    double c = e0;
    double alpha;
    if (a > 0.0)  // we have a minimum somewhere
      alpha = -b/(2.0*a);
    else if (e0>e2)
      alpha = 2.0;
    else
      alpha = 0.0;
    
    alpha = max(0.0, alpha);
    alpha = min(2.0, alpha);

    if (Verbose)
      cerr << "alpha = " << alpha << endl;
    Bands = alpha * Bands1 + (1.0-alpha)*Bands2;
    GramSchmidt(Bands);

    // Now extrapolate density
    Rho_r = alpha * Rho_r1 + (1.0-alpha)*Rho_r2;
    // For the benefit of the charge mixer
    NewRho = Rho_r;
  }
}


void
MPISystemClass::CalcOccupancies()
{
  ///////////////////////////////////////
  // Calculate mu from 0.5*(HOMO+LUMO) //
  ///////////////////////////////////////
  // First, gather all the energies from the k-Point clones.

  Array<double,2> allEnergies(Numk, NumBands), occ(Numk, NumBands);
  if (BandComm.MyProc() == 0) {
    allEnergies(Myk, Range::all()) = CG.Energies;
    kComm.AllGatherRows(allEnergies);
  }
  BandComm.Broadcast (0, allEnergies);
  // Now processor zero has all the energies;
  vector<double> Esort;
  for (int i=0; i<allEnergies.extent(0); i++)
    for (int j=0; j<allEnergies.extent(1); j++)
      Esort.push_back(allEnergies(i,j));
  sort (Esort.begin(), Esort.end());
  int HOMO = (NumElecs*Numk+1)/2-1;
  assert ((HOMO+1) < Esort.size());
  double muLo = Esort[0];
  double muHi = Esort[Esort.size()-1];
  double numkInv = 1.0/(double)Numk;
  double mu, occSum;
  do {
    occSum = 0.0;
    mu = 0.5*(muHi + muLo);
    for (int ki=0; ki<Numk; ki++)
      for (int bi=0; bi<NumBands; bi++) {
	occ(ki, bi) = 2.0*numkInv*Smearer.S(allEnergies(ki, bi), mu);
	occSum += occ(ki, bi);
      }
    if (occSum < (double) NumElecs)
      muLo = mu;
    else if (occSum > (double) NumElecs)
      muHi = mu;
  } while (fabs(occSum-(double)NumElecs) > 1.0e-10);
  

//   int LUMO = (NumElecs*Numk)/2;
//   double mu = 0.5*(Esort[HOMO] + Esort[LUMO]);
//   /// Broadcast mu to all procs in kBand;
//   if (BandComm.MyProc() == 0)
//     kComm.Broadcast(0, mu);
//   /// Now proc 0 of each k-point clone has mu.
//   /// Now broadcast mu from each BandComm Proc 0 to the rests of the
//   /// band holders
//   BandComm.Broadcast (0, mu);
//   /// Now everybody has mu!  Let's do the Methfessel-Paxton occupation 
//   double totalOcc = 0.0;
//   for (int ki=0; ki<Numk; ki++)
//     for (int bi=0; bi<NumBands; bi++) {
//       occ(ki, bi) = Smearer.S(allEnergies(ki, bi), mu);
//       totalOcc += occ(ki, bi);
//     }
  /// Now normalize to to make sure we have exactly the right number
  /// of electrons 
  for (int bi=0; bi<NumBands; bi++) 
    Occupancies(bi) = (double)(NumElecs)/occSum * occ(Myk, bi);
  if (Verbose)
    perr << "Occupancies = " << Occupancies << endl;
}

void
MPISystemClass::CalcChargeDensity()
{
  int nx, ny, nz;
  FFT.GetDims(nx,ny,nz);
  if (TempRho.shape() != TinyVector<int,3>(nx,ny,nz))
    TempRho.resize(nx,ny,nz);
  if (NewRho.shape() != TinyVector<int,3>(nx,ny,nz))
    NewRho.resize(nx,ny,nz);
  TempRho = 0.0;
  
  zVec band(Bands.extent(1));
  /// We assume that each of the bands is already normalized.
  for (int bi=CG.GetFirstBand(); bi<=CG.GetLastBand(); bi++) {
    band = Bands(bi, Range::all());
    FFT.PutkVec (band);
    FFT.k2r();
    for (int ix=0; ix<nx; ix++)
      for (int iy=0; iy<ny; iy++)
	for (int iz=0; iz<nz; iz++)
 	  TempRho(ix,iy,iz) += Occupancies(bi)*norm (FFT.rBox(ix,iy,iz));
  }
  /// First, sum over all the band procs
  BandComm.AllSum (TempRho, NewRho);
  TempRho = NewRho;
  // Now sum over all the kPoints
  if (BandComm.MyProc() == 0) 
    kComm.AllSum(TempRho, NewRho);
  // And broadcast;
  BandComm.Broadcast(0, NewRho);
  /// Multiply by the right prefactor;
  double volInv = 1.0/H.GVecs.GetBoxVol();
  double prefactor = volInv; // /(double)(nx*ny*nz);
  NewRho *= prefactor;
  double vCell = H.GVecs.GetBoxVol()/(double)(nx*ny*nz);
  double totalCharge = 0.0;
  for (int ix=0; ix<nx; ix++)
    for (int iy=0; iy<ny; iy++)
      for (int iz=0; iz<nz; iz++)
	totalCharge += vCell * NewRho(ix,iy,iz);
  if (Verbose)
    perr << "total NewRho charge = " << totalCharge << endl;
}



/// Calculate the hartree, exchange, and correlation potentials.
/// We assume that upon entry, Rho_r contains the charge density 
/// in real space and Rho_G contains it in reciporical space.
/// On exit, VH will contain the Hartree potential and VXC will
/// contain the exchange-correlation potential, both in real space.
void 
MPISystemClass::CalcVHXC()
{
  double totalCharge = 0.0;
  double vol = GVecs.GetBoxVol();
  double Ngrid = Rho_r.size();

  ///////////////////////
  // Hartree potential //
  ///////////////////////
  h_G = Rho_G;
  // Compute V_H in reciporical space
  double hartreeTerm = 0.0;
  //  double prefact = (4.0*M_PI/H.GVecs.GetBoxVol());
  double prefact = 4.0*M_PI;
  EH = 0.0;
  for (int i=0; i<h_G.size(); i++) {
    EH += norm(h_G(i)) * H.GVecs.DeltaGInv2(i);
    h_G(i) *= prefact*H.GVecs.DeltaGInv2(i);
  }
  EH *= 0.5*prefact*vol;
  // FFT back to real space
  FFT.PutkVec(h_G);
  FFT.k2r();
  copy (real(FFT.rBox), VH);

  ////////////////////////////////////
  // Exchange-correlation potential //
  ////////////////////////////////////
  EXC = 0.0;
  double exc;
  for (int ix=0; ix<Rho_r.extent(0); ix++)
    for (int iy=0; iy<Rho_r.extent(1); iy++)
      for (int iz=0; iz<Rho_r.extent(2); iz++) {
	FortranExCorr (Rho_r(ix,iy,iz), exc, VXC(ix,iy,iz));
	EXC += exc*Rho_r(ix,iy,iz);
      }
  EXC *= vol/(double)Ngrid;
  VHXC = VH + VXC;

  ///////////////////////////////////
  // Calculate electron-ion energy //
  ///////////////////////////////////
//   Eelec_ion = 0.0;
//   const Array<double,1> &VG = H.GetVG();
//   for (int i=0; i<H.GVecs.DeltaSize(); i++) {
//     complex<double> S(0.0, 0.0);
//     for (int ion=0; ion<Rions.size(); ion++) {
//       Vec3 r = Rions(ion) + 0.5*Box;
//       const Vec3 &G = H.GVecs.DeltaG(i);
//       double phase, c, s;
//       phase = dot (r, G);
//       sincos(phase, &s, &c);
//       S += complex<double> (c, s);
//     }
//     Eelec_ion += real(conj(Rho_G(i))*S*VG(i));
//   }
//   Eelec_ion *= vol;
//   double Eelec_ion2 = 0.0;
//   for (int ix=0; ix<Rho_r.extent(0); ix++)
//     for (int iy=0; iy<Rho_r.extent(1); iy++)
//       for (int iz=0; iz<Rho_r.extent(2); iz++)
// 	Eelec_ion2 += Rho_r(ix,iy,iz)*H.GetVr()(ix,iy,iz).real();

//   Eelec_ion2 *= vol/(double)Ngrid;
  Elocal    = CalcLocalPPEnergy();
  Enonlocal = CalcNonlocalPPEnergy();
  Ekinetic  = CalcKineticEnergy();
  double Ecore = NumElecs*(double)Rions.size() * H.GetVG0();
  
  if (BandComm.MyProc() == 0)
    if (kComm.MyProc() == 0) {
      fprintf (stderr, 
	       "Iter      Ekinetic EH        EXC       Elocal    Enonlocal Eewald    Ecore\n");
      fprintf (stderr, 
	       "%3d    %9.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f\n", 
	       SCiter, Ekinetic, EH, EXC, Elocal, Enonlocal, EwaldEnergy(), 
	       Ecore);
    }	       
//   perr << "Energies:\n"
//        << "   EH         = " << EH  << endl
//        << "   EXC        = " << EXC << endl 
//        << "   Eelec_ion  = " << Eelec_ion << endl
//        << "   Eewald     = " << EwaldEnergy() << endl
//        << "   Ecore      = " << Ecore << endl;
}


double
MPISystemClass::CalcLocalPPEnergy()
{
  double E = 0.0;
  for (int ix=0; ix<Rho_r.extent(0); ix++)
    for (int iy=0; iy<Rho_r.extent(1); iy++)
      for (int iz=0; iz<Rho_r.extent(2); iz++)
	E += Rho_r(ix,iy,iz)*H.GetVr()(ix,iy,iz).real();
  double vol = GVecs.GetBoxVol();
  double Ngrid = Rho_r.size();
  E *= vol/Ngrid;
  return E;
}

double
MPISystemClass::CalcNonlocalPPEnergy()
{
  zVec band(Bands.extent(1));
  double Enl = 0.0;
  /// We assume that each of the bands is already normalized.
  for (int bi=CG.GetFirstBand(); bi<=CG.GetLastBand(); bi++) {
    band = Bands(bi, Range::all());
    Enl += H.Vion->NonlocalEnergy(band) * Occupancies(bi);
  }
  BandComm.AllSum(Enl);
  return Enl;
}

double
MPISystemClass::CalcKineticEnergy()
{
  zVec band(Bands.extent(1));
  double KE = 0.0;
  /// We assume that each of the bands is already normalized.
  for (int gi=0; gi<GVecs.size(); gi++) {
    Vec3 GplusK = GVecs(gi) + GVecs.Getk();
    double E = 0.5 * dot(GplusK, GplusK);
    for (int bi=CG.GetFirstBand(); bi<=CG.GetLastBand(); bi++) 
      KE += Occupancies(bi)*norm (Bands(bi, gi)) * E;
  }
  BandComm.AllSum(KE);
  return KE;
}


void
MPISystemClass::CalcRadialChargeDensity()
{
  perr << "Calculating atomic charge density:\n";
  double maxBox = max(Box[0], max(Box[1], Box[2]));
  double rmax = max (maxBox, 100.0);
  AtomGrid.Init (1.0, rmax);
  DFTAtom atom;
  atom.NewMix = 0.5;
  atom.RadialWFs.resize(1);
  atom.RadialWFs(0).n = 1;
  atom.RadialWFs(0).l = 0;
  atom.RadialWFs(0).Occupancy = 1.0;
  atom.RadialWFs(0).Energy = -0.15;
  atom.SetGrid (&AtomGrid);
  atom.SetBarePot (V_elec_ion);
  atom.Solve();
  Array<double,1> rho(AtomGrid.NumPoints);
  for (int i=0; i<rho.size(); i++) {
    double r = AtomGrid(i);
    double u = atom.RadialWFs(0).u(i);
    rho(i) = u*u/(r*r);
  }
  RadialChargeDensity.Init (&AtomGrid, rho);
  perr << "Done atomic charge density.\n";
}

/// Initializes the density to a superposition of the atomic charge
/// densities from an Atomic DFT calculation.
void
MPISystemClass::InitNewRho()
{
  NewRho = 0.0;
  double nxInv = 1.0/(double)(Nx-1);
  double nyInv = 1.0/(double)(Ny-1);
  double nzInv = 1.0/(double)(Nz-1);
  Vec3 boxInv = Vec3(1.0/Box[0], 1.0/Box[1], 1.0/Box[2]);
  Vec3 r;
  double chargePerAtom = (double)NumElecs/(double)Rions.size();
  perr << "chargePerAtom = " << chargePerAtom << endl;
  double totalCharge = 0.0;
  double cellVol = Box[0]*Box[1]*Box[2]/(double)(Nx*Ny*Nz);
  for (int ix=0; ix<Nx; ix++) {
    r[0] = (nxInv*(double)ix-0.5)*Box[0];
    for (int iy=0; iy<Ny; iy++) {
      r[1] = (nyInv* (double)iy-0.5)*Box[1];
      for (int iz=0; iz<Nz; iz++) {
	r[2] = (nzInv*(double)iz-0.5)*Box[2];
	for (int ion=0; ion<Rions.size(); ion++) {
	  Vec3 disp = r-Rions(ion);
	  /// Find nearest image
	  disp[0] -= nearbyint(disp[0]*boxInv[0])*Box[0];
	  disp[1] -= nearbyint(disp[1]*boxInv[1])*Box[1];
	  disp[2] -= nearbyint(disp[2]*boxInv[2])*Box[2];
	  double dist = sqrt(dot(disp, disp));
	  double rho = RadialChargeDensity(dist);
	  if (rho < 0.0)
	    perr << "Negative rho at (ix,iy,iz) = (" 
		 << ix << ", " << iy << ", " << iz << ")\n";
	  NewRho(ix,iy,iz) += chargePerAtom*rho;
	}
	totalCharge += NewRho(ix,iy,iz);
      }
    }
  }
  totalCharge *= cellVol;
  NewRho *= (double)NumElecs/totalCharge;
}
	 
void 
MPISystemClass::MixChargeDensity()
{
  ChargeMixer->Mix(NewRho, Rho_r, Rho_G);
}

inline complex<double> e2iphi(double phi)
{
  double s,c;
  sincos(phi, &s, &c);
  return complex<double>(c,s);
}

// Note, this returns only the forces of the electrons on the ions,
// not the ions on the ions.
void
MPISystemClass::CalcIonForces (Array<Vec3,1> &F)
{

  if (F.size() != Rions.size())
    F.resize(Rions.size());
  ///////////////////////
  // Electron-ion part //
  ///////////////////////
  const Array<double,1> &VG = H.GetVG();
  F = 0.0;
  for (int ion=0; ion<Rions.size(); ion++) {
    Vec3 r = Rions(ion) + 0.5*Box;
    for (int i=0; i<H.GVecs.DeltaSize(); i++) {
      const Vec3 &G = H.GVecs.DeltaG(i);
      complex<double> z = e2iphi(dot(G,r));
      F(ion) += G*imag(conj(Rho_G(i))*z)*VG(i);
    }
  }
  F *= H.GVecs.GetBoxVol();

  Vec3 boxInv (1.0/Box[0], 1.0/Box[1], 1.0/Box[2]);

  //////////////////
  // Ion-ion part //
  //////////////////
  // Short-range
  for (int i=0; i<Rions.size(); i++)
    for (int j=0; j<Rions.size(); j++)
      if (i != j) {
	Vec3 disp = Rions(j)-Rions(i);
	disp[0] -= floor(disp[0]*boxInv[0]+0.5)*Box[0];
	disp[1] -= floor(disp[1]*boxInv[1]+0.5)*Box[1];
	disp[2] -= floor(disp[2]*boxInv[2]+0.5)*Box[2];
	double dist = sqrt(dot(disp,disp));
	disp = (1.0/dist)*disp;
	double dV = V_ion_ion->dVdr(dist) - Vlong.Deriv(dist);
	F(i) += dV*disp;
      }
  // Long-range
  // Compute structure factor
  Array<complex<double>,1> S(GVecs.size());
  S = 0.0;
  for (int gi=0; gi<H.GVecs.size(); gi++) {
    Vec3 G = H.GVecs(gi);
    for (int ri=0; ri<Rions.size(); ri++) {
      Vec3 r = Rions(ri);
      S(gi) += e2iphi(dot(G,r));
    }
  }
  // Now compute force
  for (int ri=0; ri<Rions.size(); ri++) {
    Vec3 r = Rions(ri);
    for (int gi=0; gi<H.GVecs.size(); gi++) {
      Vec3 G = H.GVecs(gi);
      F(ri) += imag(conj(S(gi))*e2iphi(dot(G,r)))*Vlong_G(gi)*G;
    }
  }
}


double
MPISystemClass::EwaldEnergy()
{
  double Eewald = 0.0;
  Vec3 boxInv (1.0/Box[0], 1.0/Box[1], 1.0/Box[2]);

  for (int i=0; i<Rions.size(); i++)
    for (int j=0; j<Rions.size(); j++) 
      if (i != j) {
	Vec3 disp = Rions(j)-Rions(i);
	disp[0] -= floor(disp[0]*boxInv[0]+0.5)*Box[0];
	disp[1] -= floor(disp[1]*boxInv[1]+0.5)*Box[1];
	disp[2] -= floor(disp[2]*boxInv[2]+0.5)*Box[2];
	double dist = sqrt(dot(disp,disp));
	Eewald += 0.5*(V_ion_ion->V(dist)-Vlong(dist));
      }

  for (int gi=0; gi < H.GVecs.size(); gi++) {
    complex<double> S(0.0, 0.0);
    Vec3 G = H.GVecs(gi);
    for (int ri=0; ri<Rions.size(); ri++) {
      Vec3 r = Rions(ri);
       S += e2iphi(dot(r,G));
    }
    Eewald += 0.5*norm(S)*Vlong_G(gi);
  }
  double N = Rions.size();
  Eewald -= 0.5 * Vlong(0.0) * N;
  Eewald -= 0.5 * Vshort_G0  * N * N;
  return Eewald;
}

class VshortIntegrand
{
private:
  Potential &V;
  CubicSpline &Vlong;
public:
  inline double operator()(double r)
  {
    double vshort = V.V(r) - Vlong(r);
    return r*r*vshort;
  }

  VshortIntegrand(Potential &v, CubicSpline &vlong) :
    V(v), Vlong(vlong)
  { /* do nothing */ }
};

void
MPISystemClass::DoOptimizedBreakup()
{
  double rcut = 0.5*min(Box[0],min(Box[1],Box[2]));
  double vol = Box[0]*Box[1]*Box[2];
  double volInv = 1.0/vol;
  double kvol = (8.0*M_PI*M_PI*M_PI)/vol;
  double kavg = pow(kvol, 1.0/3.0);
  
  LPQHI_BasisClass basis;
  basis.Set_rc(rcut);
  basis.SetBox(Box);
  basis.SetNumKnots (10);
  double kCont = 50.0 * kavg;
  double delta = basis.GetDelta();
  double kMax = 20.0*M_PI/delta;
  OptimizedBreakupClass breakup(basis);
  breakup.SetkVecs (kCut, kCont, kMax);
  int numk = breakup.kpoints.size();
  int N = basis.NumElements();
  Array<double,1> t(N);
  Array<bool,1>   adjust (N);
  Array<double,1> Xk(numk);

  // Setup the grid for storing the long range potential
  double rmax = 0.500001*sqrt(dot(Box,Box));
  const int numPoints = 1000;
  VlongGrid.Init(0.0, rmax, numPoints);
  for (int ki=0; ki<numk; ki++)
    Xk(ki) = volInv * V_ion_ion->X_k(rcut, breakup.kpoints(ki)[0]);
  
  adjust = true;
  perr << "Doing ion-ion breakup...\n";
  breakup.DoBreakup (Xk, t, adjust);
  perr << "Completed.\n";
  
  ///////////////////////////////
  // Calculate real-space part //
  ///////////////////////////////
  Array<double,1> Vlong_r(numPoints);
  Vlong_r = 0.0;
  for (int i=0; i<numPoints; i++) {
    double r = VlongGrid(i);
    if (r <= rcut) {
      // Sum over basis functions
      for (int n=0; n<N; n++) 
	Vlong_r(i) += t(n) * basis.h(n, r);
    }
    else
      Vlong_r(i) = V_ion_ion->V(r);
  }
  Vlong.Init(&VlongGrid, Vlong_r);

  //////////////////////////////////////
  // Calculate reciporical-space part //
  //////////////////////////////////////
  Vlong_G.resize(H.GVecs.size());
  Vlong_G = 0.0;
  for (int ki=0; ki < H.GVecs.size(); ki++) {
    const Vec3 &kv = H.GVecs(ki);
    double k = sqrt (dot(kv,kv));
    // Sum over basis functions
    for (int n=0; n<N; n++)
      Vlong_G(ki) += t(n) * basis.c(n,k);
    // Now add on part from rc to infinity
    Vlong_G(ki) -= volInv*V_ion_ion->X_k(rcut, k);
    if (k < 1.0e-10)
      Vlong_G(ki) = 0.0;
  }

  // Now, write out a file for checking purposes
  if (BandComm.MyProc()==0)
    if (kComm.MyProc()==0) {
      FILE *fout = fopen ("Vion.dat", "w");
      for (int i=0; i<numPoints; i++)
// 	fprintf (fout, "%1.12e %1.12e %1.12e\n",
// 		 VlongGrid(i), V_ion_ion->V(VlongGrid(i)), Vlong_r(i));
 	fprintf (fout, "%1.12e %1.12e %1.12e\n",
 		 VlongGrid(i), V_ion_ion->V(VlongGrid(i)), Vlong(VlongGrid(i)));
      fclose(fout);
    }
  
  ///////////////////////////////////////////// 
  // Compute the fourier transform of Vshort //
  // for G=0.                                //
  /////////////////////////////////////////////
  VshortIntegrand integrand (*V_ion_ion, Vlong);
  GKIntegration<VshortIntegrand, GK31> integrator(integrand);
  Vshort_G0 = 4.0*M_PI/H.GVecs.GetBoxVol() *
    integrator.Integrate(1.0e-100, rcut, 1.0e-10);
}



