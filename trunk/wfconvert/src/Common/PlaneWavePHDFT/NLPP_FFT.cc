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

#include "NLPP_FFT.h"


////////////////////////////////////////////////////////////
//                  IonProjector stuff                    //
////////////////////////////////////////////////////////////

complex<double> 
Ion_l_Projector::Ylm(int l, int m, Vec3 r)
{
  if (l == 0)
    return complex<double>(1.0/sqrt(4.0*M_PI), 0.0);

  if (m < 0) {
    double sign = 1.0;
    for (int i=0; i<m; i++)
      sign *= -1.0;
    return sign * conj(Ylm(l, -m, r));
  }

  Vec3 omega = 1.0/ sqrt(dot(r, r)) * r;
  double costheta = omega[2];
  double sintheta = sqrt(1.0-costheta*costheta);
  double cosphi, sinphi;
      
  if (sintheta != 0.0) {
    double sinthetaInv = 1.0/sintheta;
    cosphi = sinthetaInv*omega[0];
    sinphi = sinthetaInv*omega[1];
  }
  else { // We're at the poles, phi is undefined
    cosphi = 1.0;
    sinphi = 0.0;
  }
  if (l == 1) {
    if (m==0)
      return -sqrt(3.0/(4.0*M_PI)) * costheta *	complex<double>(1.0, 0.0);
    else if (m==1)
      return -sqrt(3.0/(8.0*M_PI))*sintheta * 
	complex<double> (cosphi, sinphi);
    else {
      cerr << "Invalid l,m\n";
      abort();
    }
  }
  else if (l==2) {
    if (m==0)
      return sqrt(5.0/(4.0*M_PI)) * (1.5*costheta *costheta - 0.5) *
	complex<double>(1.0, 0.0);
    else if (m == 1)
      return -sqrt(15.0/(8.0*M_PI))*sintheta*costheta*
	complex<double>(cosphi, sinphi);
    else if (m==2){
      double cos2phi = cosphi*cosphi - sinphi*sinphi;
      double sin2phi = 2.0*sinphi*cosphi;
      return 0.25*sqrt(15.0/(2.0*M_PI))*sintheta*sintheta *
	complex<double>(cos2phi, sin2phi);
    }
  }
  else {
    cerr << "Ylm not implemented for l > 2.\n";
    abort();
  }
  return 0.0;
}


complex<double> 
Ion_l_Projector::Ylm2(int l, int m, Vec3 r)
{
  double r2 = dot(r,r);
  if (l == 0 || r2==0.0)
    return 0.5*sqrt(1.0/M_PI);

  if (m < 0) {
    double sign = ((-m)&1) ? -1.0 : 1.0;
    return sign * conj(Ylm(l, -m, r));
  }

  double nrm = sqrt (1.0/dot(r, r));
  double x = nrm*r[0];
  double y = nrm*r[1];
  double z = nrm*r[2];
  if (l==1) {
    if (m==0) 
      return 0.5*sqrt(3.0/M_PI)*z;
    else if (m==1)
      return -0.5*sqrt(3.0/(2.0*M_PI))*complex<double>(x,y);
  }
  else if (l==2) {
    if (m==0)
      return 0.25*sqrt(5.0/M_PI)*(3.0*z*z - 1.0);
    else if (m==1)
      return -0.5*sqrt(15.0/(2.0*M_PI))*z*complex<double>(x,y);
    else if (m==2) {
      complex<double> xy(x,y);
      return 0.25*sqrt(15.0/(2.0*M_PI))*xy*xy;
    }
  }
  else {
    cerr << "Ylm not implemented for l > 2.\n";
    abort();
  }
  return 0.0;
}


// void
// Ion_l_Projector::Setup(NLPPClass &nlpp, int l_, 
// 		       Vec3 rion, Vec3 k, FFTBox &fft, bool smooth)
// {
//   FFT = &fft;
//   l = l_;
//   double R0 = smooth ? nlpp.GetR0(l) : nlpp.Getrc(l);
//   int nx, ny, nz;
//   fft.GetDims (nx, ny, nz);
//   Vec3 box = fft.GVecs.GetBox();
//   Vec3 boxInv(1.0/box[0], 1.0/box[1], 1.0/box[2]);
//   double dx = box[0]/(double)nx;
//   double dy = box[1]/(double)ny;
//   double dz = box[2]/(double)nz;
//   Vec3 r;
//   vector<Int3> indices;
//   vector<Vec3> disps;
//   for (int ix=0; ix<nx; ix++) {
//     r[0] = dx * ix;
//     for (int iy=0; iy<ny; iy++) {
//       r[1] = dy * iy;
//       for (int iz=0; iz<nz; iz++) {
// 	r[2] = dz * iz;
// 	Vec3 disp = r - rion;
// 	disp[0] -= round(disp[0]*boxInv[0])*box[0];
// 	disp[1] -= round(disp[1]*boxInv[1])*box[1];
// 	disp[2] -= round(disp[2]*boxInv[2])*box[2];
// 	if (dot(disp, disp) <= R0*R0) {
// 	  // This point is inside the projection core
// 	  indices.push_back(Int3(ix,iy,iz));
// 	  disps.push_back(disp);
// 	}
//       }
//     }
//   }
//   // Now, compute projector inside the core
//   ChiYlm.resize(indices.size(), 2*l+1);
//   FFTIndices.resize(indices.size());
//   // Psi is normalized such that \sum_{rFFT} = Nx*Ny*Nz.
//   // We adjust the normalization of the projector so that it has the
//   // same convention. 
//   double normFactor = sqrt(box[0]*box[1]*box[2]);
//   for (int i=0; i<indices.size(); i++) {
//     FFTIndices(i) = indices[i];
//     double r = sqrt(dot(disps[i], disps[i]));
//     for (int m=-l; m<=l; m++) {
//       double chi_r = smooth ? nlpp.GetChi_r(l,r) : nlpp.GetZeta_r(l, r);
//       ChiYlm(i,l+m) = normFactor*Ylm2(l,m,disps[i]) * chi_r;
//       if (isnan(real(ChiYlm(i,l+m))) || isnan(imag(ChiYlm(i,l+m))))
// 	cerr << "NAN at r = " << r << " m=" << m << "  l=" << l << endl;
//     }
//   }
//   // Renormalize the real-space projectors
// //   for (int m=-l; m<=l; m++) {
// //     double nrm = 0.0;
// //     for (int i=0; i<indices.size(); i++)
// //       nrm += norm (ChiYlm(i,l+m));
// //     ChiYlm(Range::all(),l+m) *= sqrt((double)(nx*ny*nz)/nrm);
// //   }
//   MeshVol = dx*dy*dz;
// }


void
Ion_l_Projector::Setup(NLPPClass &nlpp, int l_, 
		       Vec3 rion, Vec3 k, FFTBox &fft, bool smooth)
{
  FFT = &fft;
  l = l_;
  double R0 = smooth ? nlpp.GetR0(l) : nlpp.Getrc(l);
  int nx, ny, nz;
  fft.GetDims (nx, ny, nz);
  Vec3 box = fft.GVecs.GetBox();
  Vec3 boxInv(1.0/box[0], 1.0/box[1], 1.0/box[2]);
  double dx = box[0]/(double)nx;
  double dy = box[1]/(double)ny;
  double dz = box[2]/(double)nz;
  Vec3 r;
  vector<Int3> indices;
  vector<Vec3> disps;
  for (int ix=0; ix<nx; ix++) {
    r[0] = dx * ix;
    for (int iy=0; iy<ny; iy++) {
      r[1] = dy * iy;
      for (int iz=0; iz<nz; iz++) {
	r[2] = dz * iz;
	Vec3 disp = r - rion;
	disp[0] -= round(disp[0]*boxInv[0])*box[0];
	disp[1] -= round(disp[1]*boxInv[1])*box[1];
	disp[2] -= round(disp[2]*boxInv[2])*box[2];
	if (dot(disp, disp) <= R0*R0) {
	  // This point is inside the projection core
	  indices.push_back(Int3(ix,iy,iz));
	  disps.push_back(disp);
	}
      }
    }
  }
  
  ChiYlm.resize(indices.size(), 2*l+1);
  FFTIndices.resize(indices.size());
  FFTOffsets.resize(indices.size());

  GVecsClass &gvecs = fft.GVecs;
  // Now, compute projectors in reciprocal space
  int nG = gvecs.DeltaSize();
  complex<double> i2l (1.0);
  for (int i=0; i<l; i++)
    i2l *= complex<double>(0.0, 1.0);
  double prefactor = 4.0*M_PI/sqrt(gvecs.GetBoxVol());
 
  zVec kProjector(nG);
  kProjector = complex<double>();
  for (int m=-l; m<=l; m++) {
    for (int ig=0; ig<nG; ig++) {
      Vec3 g = gvecs.DeltaG(ig) + k;
      double gmag = sqrt(dot(g,g));
      if (gmag < 3.0*gvecs.GetkCut()) {
	double chi_q = smooth ? 
	  nlpp.GetChi_q(l, gmag) : nlpp.GetZeta_q(l, gmag);
	complex<double> ylm = Ylm2 (l, m, g);
	double phase = dot (g, rion);
	kProjector(ig) = 
	  prefactor*i2l*ylm*chi_q*complex<double>(cos(phase), sin(phase));
      }
    }
    fft.PutkVec (kProjector);
    fft.k2r();
    for (int i=0; i<indices.size(); i++) {
      FFTIndices(i) = indices[i];
      FFTOffsets(i) = (&fft.rBox(indices[i]) - &fft.rBox(0,0,0));
      ChiYlm(i,l+m) = fft.rBox(indices[i]);
    }
  }
  MeshVol = dx*dy*dz;
}



void
NLPP_FFTClass::SetupkProjectors()
{
  GVecsClass &gvecs = cFFT.GVecs;
  int numproj = 0;
  for (int l=0; l<NLPP.NumChannels(); l++)
    if (l != NLPP.LocalChannel())
      numproj += 2*l+1;

  lambda_lm.resize(gvecs.size(), Rions.size(), numproj);
  Ion_l_Projector proj;
  double prefactor = 4.0*M_PI/sqrt(gvecs.GetBoxVol());

  for (int ig=0; ig<gvecs.size(); ig++) {
    Vec3 g = gvecs(ig) + kPoint;
    double gmag = sqrt(dot(g,g));
    int iproj = 0;
    complex<double> i2l(1.0, 0.0);
    for (int l=0; l<NLPP.NumChannels(); l++) {
      if (l != NLPP.LocalChannel()) {
      double zeta_q = NLPP.GetZeta_q(l, gmag);
	for (int m=-l; m<=l; m++) {
	  complex<double> ylm = proj.Ylm2 (l, m, g);
	  for (int ion=0; ion<Rions.size(); ion++) {
	    double phase = dot (g, Rions(ion));
	    lambda_lm(ig, ion, iproj) = 
	      prefactor*i2l*ylm*zeta_q*complex<double>(cos(phase), sin(phase));
	  }
	  iproj++;
	}
      }
      i2l *= complex<double>(0.0, 1.0);
    }
  }
  for (int ion=0; ion < Rions.size(); ion++) {
    int iproj = 0;
    for (int l=0; l<NLPP.NumChannels(); l++) {
      if (l != NLPP.LocalChannel()) {
	for (int m=-l; m<=l; m++) {
	  double nrm = 0.0;
	  for (int ig=0; ig<gvecs.size(); ig++) 
	    nrm += norm(lambda_lm(ig, ion, iproj));
	  // lambda_lm(Range::all(), ion, iproj) *= 1.0/sqrt(nrm);
	  iproj++;
	}
      }
    }
  }

}






// Note:  this presently only includes the local potential
void
NLPP_FFTClass::Vmatrix (Array<complex<double>,2> &vmat)
{
  if (!IsSetup)
    Setup();
  double volInv = 1.0/GVecs.GetBoxVol();
  for (int i=0; i<vmat.rows(); i++) 
    for (int j=0; j<=i; j++) {
      Vec3 diff = GVecs(i) - GVecs(j);
      complex<double> s(0.0,0.0);
      for (int zi=0; zi<Rions.size(); zi++) {
	double cosVal, sinVal, phase;
	phase = dot (diff, Rions(zi));
	sincos(phase, &sinVal, &cosVal);
	s += complex<double> (cosVal,sinVal);
      }
      vmat (i,j) = s*kPH.V(kPoint, GVecs(i), GVecs(j))*volInv;
      vmat (j,i) = conj(vmat(i,j));
    }
}


void 
NLPP_FFTClass::SetupkPotentials()
{
  // First, create Vlocal from the local potential of the NLPP
  Vlocal.Spline = NLPP.GetLocalSpline();
  Vouter.Z1Z2   = -NLPP.GetValenceCharge();
  Vlocal.Vouter = &Vouter;

  // Now, compute the tail coefficients of the local potential
  kPH.CalcTailCoefs (30.0, 60.0);

  // Compute local potential in reciprocal space
  double volInv = 1.0/GVecs.GetBoxVol();
  // Setup V and F tensors in k-space
  //  double gMag, lastMag2;
  double gMag, lastMag2;
  lastMag2 = -1.0;
  double a, bPerp, bPar, V;
  int numCalls = 0;
  for (int i=0; i<GVecs.DeltaSize(); i++) {
    double gMag2 = dot(GVecs.DeltaG(i), GVecs.DeltaG(i));
    if (fabs(gMag2-lastMag2) > 1.0e-12) {
      lastMag2 = gMag2;
      gMag  = sqrt(lastMag2);
      kPH.GetVals(gMag, a, bPerp, bPar, V);
      numCalls++;
    }
    if (gMag2 > 1.0e-10)
      VG(i) = volInv * V;
    else {
      VG(i) = 0.0;
      VG0 = volInv * V;
    }
  }
}


void
NLPP_FFTClass::SetuprPotentials()
{
  // Setup local part
  cFFT.kBox = complex<FFT_FLOAT>(0.0, 0.0);
  for (int i=0; i<GVecs.DeltaSize(); i++) {
    Int3 I = GVecs.DeltaI(i);
    cFFT.kBox(I)   = StructureFactor(i)*VG(i);
  }
  cFFT.k2r();
  Vr = cFFT.rBox;

  // Setup nonlocal part
  cerr << "Setting up nonlocal projectors.\n";
  int numProj = NLPP.NumChannels()-1;
  Ion_l_Projectors.resize( Rions.size(), numProj);
  for (int ri=0; ri<Rions.size(); ri++) {
    int iProj = 0;
    for (int l=0; l<NLPP.NumChannels(); l++) {
      if (l != NLPP.LocalChannel()) {
	Ion_l_Projectors(ri, iProj).Setup 
	  (NLPP, l, Rions(ri), kPoint, cFFT, SmoothProjectors);
	iProj++;
      }
    }
  }
}

void
NLPP_FFTClass::SetProjectors(bool smooth)
{
  if (smooth != SmoothProjectors) {
    SmoothProjectors = smooth;
    SetupkProjectors();
    SetuprPotentials();
  }
}

void
NLPP_FFTClass::SetIons(const Array<Vec3,1> &rions)
{
  // Calculate the structure factor
  VionBase::SetIons(rions);
  if (IsSetup) {
    SetupkProjectors();
    SetuprPotentials();
  }
}



void
NLPP_FFTClass::Setup()
{
  int nx, ny, nz;
  cFFT.GetDims(nx,ny,nz);
  Vr.resize(nx,ny,nz);
  Vc.resize(GVecs.size());
  VG.resize(GVecs.DeltaSize());
  VnlPsi.resize(nx,ny,nz);

  // Compute the Kleinmain-Bylander projectors:
  double kc = cFFT.GVecs.GetkCut();
  NLPP.SetupProjectors(kc, 4.0*kc);

  SetupkPotentials();
  SetupkProjectors();
  SetuprPotentials();

  IsSetup = true;
}

void
NLPP_FFTClass::Setk (Vec3 k)
{
  kPoint = k;
  
  int nx, ny, nz;
  cFFT.GetDims(nx,ny,nz);
  Vr.resize(nx,ny,nz);
  Vc.resize(GVecs.size());
  VG.resize(GVecs.DeltaSize());

  SetupkPotentials();
  SetIons(Rions);
}

void
Ion_l_Projector::Project (Array<complex<double>,1> &chi_psi)
{
  chi_psi = complex<double>();
  int nx, ny, nz;
  FFT->GetDims(nx,ny,nz);
  double normFactor = 1.0/(double)(nx*ny*nz);

  int numPoints = FFTIndices.size();
  for (int i=0; i<numPoints; i++) {
    //    complex<double> psi = FFT->rBox(FFTIndices(i));
    complex<double> psi = FFT->rBox.data()[FFTOffsets(i)];
    for (int m=-l; m<=l; m++) 
      chi_psi(m+l) += conj(ChiYlm(i,m+l)) * psi;
  }
  chi_psi *= normFactor;
}

void
Ion_l_Projector::AddToVnl (Array<complex<double>,1> &Echi_psi,
			   Array<complex<double>,3> &VnlPsi)
{
  int N = FFTIndices.size();
  for (int i=0; i<N; i++) {
    //    Int3 index = FFTIndices(i);
    int offset = FFTOffsets(i);
    for (int m=-l; m<=l; m++)
      //      VnlPsi(index) += ChiYlm(i,m+l)*Echi_psi(m+l);
      VnlPsi.data()[offset] += ChiYlm(i,m+l)*Echi_psi(m+l);
  }
}


double
NLPP_FFTClass::NonlocalEnergy(const zVec &c)
{
  cFFT.PutkVec (c);
  cFFT.k2r();

  int lmax = NLPP.NumChannels()-1;
  Array<complex<double>,1> chi_psi(2*lmax+1);
  double E_nl = 0.0;
  
  for (int ri=0; ri<Rions.size(); ri++) {
    int iProj = 0;
    for (int l=0; l<NLPP.NumChannels(); l++) 
      if (l != NLPP.LocalChannel()) {
	double E_KB = NLPP.GetE_KB(l);
	Ion_l_Projector &proj = Ion_l_Projectors(ri, iProj);
	proj.Project(chi_psi);
	for (int m=-l; m<=l; m++) 
	  E_nl += norm (chi_psi(m+l))*E_KB;
	iProj++;
      }
  }
  return E_nl;
}

// double
// NLPP_FFTClass::NonlocalEnergy(const zVec &c)
// {
//   int lmax = NLPP.NumChannels()-1;
//   double E_nl = 0.0;
//   GVecsClass &gvecs = cFFT.GVecs;
//   for (int ion=0; ion<Rions.size(); ion++) {
//     int iproj = 0;
//     for (int l=0; l<NLPP.NumChannels(); l++) {
//       if (l != NLPP.NumChannels()) {
// 	for (int m=-l; m<=l; m++) {
// 	  complex<double> Zlm(0.0, 0.0);
// 	  for (int ig=0; ig<gvecs.size(); ig++)
// 	    Zlm += conj (lambda_lm(ig, ion, iproj)) * c(ig);
// 	  E_nl += norm(Zlm) * NLPP.GetE_KB(l);
// 	  iproj++;
// 	}
//       }
//     }
//   }
//   return E_nl;
// }


void
NLPP_FFTClass::CalcVnlPsi()
{
  VnlPsi = complex<double>(0.0, 0.0);
  int lmax = NLPP.NumChannels()-1;
  Array<complex<double>,1> chi_psi(2*lmax+1), Echi_psi(2*lmax+1);

  double E_nl = 0.0;

  for (int ri=0; ri<Rions.size(); ri++) {
    int iProj = 0;
    for (int l=0; l<NLPP.NumChannels(); l++) 
      if (l != NLPP.LocalChannel()) {
	double E_KB = NLPP.GetE_KB(l);
	Ion_l_Projector &proj = Ion_l_Projectors(ri, iProj);
	proj.Project(chi_psi);
	Echi_psi = E_KB * chi_psi;
	proj.AddToVnl (Echi_psi, VnlPsi);
	for (int m=-l; m<=l; m++) {
	  E_nl += norm (chi_psi(m+l))*E_KB;
	}
	iProj++;
      }
  }
  //cerr << "E_nl = " << E_nl << endl;
}


void
NLPP_FFTClass::ApplyNonlocal (const zVec &c, zVec &Hc)
{
  GVecsClass &gvecs = cFFT.GVecs;
  for (int ion=0; ion<Rions.size(); ion++) {
    int iproj = 0;
    for (int l=0; l<NLPP.NumChannels(); l++) {
      if (l != NLPP.NumChannels()) {
	for (int m=-l; m<=l; m++) {
	  complex<double> Zlm(0.0, 0.0);
	  for (int ig=0; ig<gvecs.size(); ig++)
	    Zlm += conj(lambda_lm(ig, ion, iproj)) * c(ig);
	  Zlm *= NLPP.GetE_KB(l);
	  for (int ig=0; ig<gvecs.size(); ig++)
	    Hc(ig) += Zlm * (lambda_lm(ig, ion, iproj));
	  iproj++;
	}
      }
    }
  }
}


void 
NLPP_FFTClass::Apply (const zVec &c, zVec &Hc)
{
  if (!IsSetup)
    Setup();
  int nx, ny, nz;
  cFFT.GetDims(nx, ny, nz);

  //////////////////////////
  // Local potential part //
  //////////////////////////
  // Transform c into real space
  cFFT.PutkVec (c);
  cFFT.k2r();

  ////////////////////
  // Nonlocal parts //
  ////////////////////
  //ApplyNonlocal (c, Hc);
  CalcVnlPsi();

  // Multiply by V
  cFFT.rBox *= Vr;
  // Add nonlocal parts
  cFFT.rBox += VnlPsi;

  // Transform back
  cFFT.r2k();

  // Get vector
  cFFT.GetkVec (Vc);
  for (int i=0; i<GVecs.size(); i++)
    Hc(i) += Vc(i);
}

void 
NLPP_FFTClass::Apply (const zVec &c, zVec &Hc,
		      Array<double,3> &VHXC)
{
  if (!IsSetup)
    Setup();
  int nx, ny, nz;
  cFFT.GetDims(nx, ny, nz);

  ////////////////////
  // Potential part //
  ////////////////////
  // Transform c into real space
  cFFT.PutkVec (c);
  cFFT.k2r();
  ////////////////////
  // Nonlocal parts //
  ////////////////////
  //ApplyNonlocal(c, Hc);
  CalcVnlPsi();

  // Apply local potential and VHXC
  //  cFFT.rBox *= (Vr+VHXC);
  for (int ix=0; ix<nx; ix++)
    for (int iy=0; iy<ny; iy++)
      for (int iz=0; iz<nz; iz++)
	cFFT.rBox(ix,iy,iz) *= (Vr(ix,iy,iz)+VHXC(ix,iy,iz));

  // Add nonlocal parts
  cFFT.rBox += VnlPsi;

  // Transform back
  cFFT.r2k();

  // Get vector
  cFFT.GetkVec (Vc);
  for (int i=0; i<GVecs.size(); i++)
    Hc(i) += Vc(i);
}
