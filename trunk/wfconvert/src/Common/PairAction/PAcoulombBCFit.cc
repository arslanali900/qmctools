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

#include "PAcoulombBCFit.h"
#include "../Splines/BicubicSpline.h"
#include "../Fitting/Fitting.h"

/// The following routines are used only if we are creating fits, not
/// using them.

void PAcoulombBCFitClass::ReadParams(IOSectionClass &inSection)
{
  assert(inSection.OpenSection ("qGrid"));
  qgrid = ReadGrid (inSection);
  inSection.CloseSection();
  assert(inSection.OpenSection ("tGrid"));
  tgrid = ReadGrid (inSection);
  inSection.CloseSection();
  GridsAreMine = true;
}

void PAcoulombBCFitClass::WriteBetaIndependentInfo (IOSectionClass &outSection)
{
  outSection.WriteVar ("Type", "coulombBCfit");
  outSection.NewSection("qGrid");
  qgrid->Write(outSection);
  outSection.CloseSection();
  outSection.NewSection("tGrid");
  tgrid->Write(outSection);
  outSection.CloseSection();
}



void PAcoulombBCFitClass::DoFit (Rho &rho)
{
  NumBetas = 1;
  Usplines.resize(1);
  dUsplines.resize(1);

  SmallestBeta = rho.Beta();

  int numq = qgrid->NumPoints;
  int numt = tgrid->NumPoints;
  Array<double,2> Umat(numq, numt), dUmat(numq, numt);
  for (int qi=0; qi<numq; qi++) {
    double q = (*qgrid)(qi);
    for (int ti=0; ti<numt; ti++) {
      double t = (*tgrid)(ti);
      double s = 2.0*q*t;

      double costheta;
      if (q == 0.0)
	costheta = 1.0;
      else
	costheta = 1.0 - s*s/(2.0*q*q);
      costheta = min(1.0, costheta);
      costheta = max(-1.0, costheta);
      
      double U, dU;
      rho.UdU_Coulomb(q,q,costheta, U, dU);
      Umat(qi,ti) = U;
      dUmat(qi,ti) = dU;
    }
  }
  // Initialize the bicubic splines
  Usplines(0).Init(qgrid,tgrid,Umat);
  dUsplines(0).Init(qgrid,tgrid,dUmat);
}



void PAcoulombBCFitClass::Error(Rho &rho, double &Uerror, double &dUerror)
{
  double U2err = 0.0;
  double dU2err = 0.0;
  double weight = 0.0;
  FILE *Uxdat = fopen ("Ux.dat", "w");
  FILE *Ufdat = fopen ("Uf.dat", "w");
  FILE *dUxdat = fopen ("dUx.dat", "w");
  FILE *dUfdat = fopen ("dUf.dat", "w");
  FILE *sdat = fopen ("s.dat", "w");
  FILE *costhetadat = fopen ("costheta.dat", "w");
  FILE *qdat = fopen ("q.dat", "w");
  LinearGrid qgrid2(qgrid->Start, qgrid->End, 315);
  for (int qi=0; qi<qgrid2.NumPoints; qi++) {
    double q = qgrid2(qi);
    LinearGrid sgrid(0.0, 2.0*q, 350);
    for (int si=0; si<sgrid.NumPoints; si++) {
      double s = sgrid(si);
      double w = exp(-s*s/(4.0*rho.lambda*rho.Beta()));
      double Uex, dUex, Ufit, dUfit;
      double costheta;
      if (q == 0.0)
	costheta = 1.0;
      else
	costheta = 1.0 - s*s/(2.0*q*q); 
      rho.UdU_Coulomb(q, q, costheta, Uex, dUex);
      Ufit = U(q, 0.0, s*s, 0);
      dUfit = dU(q, 0.0, s*s, 0);
      U2err += w*(Uex-Ufit)*(Uex-Ufit);
      dU2err += w*(dUex-dUfit)*(dUex-dUfit);
      weight += w;
      fprintf (Uxdat, "%1.16e ", Uex);
      fprintf (Ufdat, "%1.16e ", Ufit);
      fprintf (dUxdat, "%1.16e ", dUex);
      fprintf (dUfdat, "%1.16e ", dUfit);
      fprintf (sdat, "%1.16e ", s);
      fprintf (costhetadat, "%1.16e ", costheta);
      fprintf (qdat, "%1.16e ", q);
    }
    fprintf (Uxdat, "\n");
    fprintf (Ufdat, "\n");
    fprintf (dUxdat, "\n");
    fprintf (dUfdat, "\n");
    fprintf (sdat, "\n");
    fprintf (qdat, "\n");
    fprintf (costhetadat, "\n");
  }
  fclose (Uxdat); fclose(Ufdat); fclose(sdat); fclose(qdat); 
  fclose(costhetadat);
  Uerror = sqrt(U2err/weight);
  dUerror = sqrt(dU2err/weight);
}


void PAcoulombBCFitClass::WriteFit (IOSectionClass &outSection)
{
  Array<double,2> Umat(qgrid->NumPoints, tgrid->NumPoints); 
  Array<double,2> dUmat(qgrid->NumPoints, tgrid->NumPoints); 
  double beta = SmallestBeta;
  outSection.NewSection("Fit");
  outSection.WriteVar ("beta", beta);
  for (int qi=0; qi<qgrid->NumPoints; qi++)
    for (int ti=0; ti<tgrid->NumPoints; ti++) {
      Umat(qi,ti) = Usplines(0)(qi,ti);
      dUmat(qi,ti) = dUsplines(0)(qi,ti);
    }
  outSection.WriteVar ("Umat", Umat);
  outSection.WriteVar ("dUmat", dUmat);
  outSection.CloseSection();
}



double PAcoulombBCFitClass::U(double q, double z, double s2, int level)
{
  if (q <= (qgrid->End*1.0000001)) {
    if (q == 0)
      return (Usplines(level)(0.0,0.0));
    else {
      double t = sqrt(s2)/(2.0*q);
      return (Usplines(level)(q,t));
    }
  }
  else {
    // Coulomb action is independent of z
    double beta = SmallestBeta;
    for (int i=0; i<level; i++)
      beta *= 2.0;
    return (beta*Pot->V(q));
  }
}

double PAcoulombBCFitClass::dU(double q, double z, double s2, int level)
{
  if (q <= (qgrid->End*1.0000001)) {
    if (q == 0)
      return (dUsplines(level)(0.0,0.0));
    else {
      double t = sqrt(s2)/(2.0*q);
      return (dUsplines(level)(q,t));
    }
  }
  else {
    // Coulomb action is independent of z
    return (Pot->V(q));
  }
}

void 
PAcoulombBCFitClass::Derivs (double q, double z, double s2, int level,
			     double &d_dq, double &d_dz)
{
  d_dz = 0.0;
  if (q <= (qgrid->End*1.0000001)) {
    if (q == 0)
      d_dq = Usplines(level).d_dx(0.0,0.0);
    else {
      double t = sqrt(s2)/(2.0*q);
      double dq = Usplines(level).d_dx(q,t);
      double dt = Usplines(level).d_dy(q,t);
      d_dq = dq - (t/q)*dt;
    }
  }
  else {
    // Coulomb action is independent of z
    d_dq = ldexp(SmallestBeta,level)*Pot->dVdr(q);
  }
}

void
PAcoulombBCFitClass::Derivs (double q, double z, double s2, int level,
			     double &d_dq, double &d_dz, double &d_ds)
{
  d_dz = 0.0;
  if (q <= (qgrid->End*1.0000001)) {
    if (q == 0) {
      d_dq = Usplines(level).d_dx(0.0,0.0);
      d_ds = 0.0;
    }
    else {
      double t = sqrt(s2)/(2.0*q);
      double dq = Usplines(level).d_dx(q,t);
      double dt = Usplines(level).d_dy(q,t);
      d_dq = dq - (t/q)*dt;
      d_ds = dt/(2.0*q);
    }
  }
  else {
    // Coulomb action is independent of z
    d_dq = ldexp(SmallestBeta,level)*Pot->dVdr(q);
    d_ds = 0.0;
  }
}



bool PAcoulombBCFitClass::Read (IOSectionClass &in,
				double smallestBeta, int numBetas)
{
  NumBetas = numBetas;
  SmallestBeta = smallestBeta;
  // Resize
  Usplines.resize(NumBetas);
  dUsplines.resize(NumBetas);
  Array<double,2> temp;

  // Read Particles;
  assert(in.OpenSection("Particle1"));
  Particle1.Read(in);
  in.CloseSection();
  assert(in.OpenSection("Particle2"));
  Particle2.Read(in);
  in.CloseSection();
  lambda = Particle1.lambda + Particle2.lambda;

  // Read Potential;
  assert(in.OpenSection("Potential"));
  Pot = ReadPotential(in);
  in.CloseSection();

  Z1Z2 = ((CoulombPot *)Pot)->Z1Z2;

  // Read the fits
  assert(in.OpenSection("Fits"));
  // Read the qgrid
  assert(in.OpenSection("qGrid"));
  qgrid = ReadGrid (in);
  in.CloseSection();
  assert(in.OpenSection("tGrid"));
  tgrid = ReadGrid (in);
  in.CloseSection();
  GridsAreMine=true;
  // Read Order

  double desiredBeta = smallestBeta;
  for (int betaIndex=0; betaIndex<NumBetas; betaIndex++) {
    int i=0;
    int numTemps = in.CountSections("Fit");
    bool found=false;
    while ((i<numTemps) && !found) {
      double beta;
      assert (in.OpenSection("Fit", i));
      assert (in.ReadVar("beta", beta));
      if ((fabs(beta-desiredBeta)/desiredBeta) < 1.0e-12)
	found = true;
      else {
	in.CloseSection();
	i++;
      }
    }
    if (!found) {
    cerr << "Couldn't find beta = " << desiredBeta 
	 << " in fit file.  Exitting\n";
    exit(1);
    }
    // Now read the fit coefficents
    assert(in.ReadVar("Umat", temp));
    Usplines(betaIndex).Init(qgrid,tgrid,temp);
    assert(in.ReadVar("dUmat", temp));
    dUsplines(betaIndex).Init(qgrid,tgrid,temp);
    in.CloseSection(); // "Fit"
    desiredBeta *= 2.0;
  }
  in.CloseSection(); // "Fits"
  return true;
}


bool PAcoulombBCFitClass::IsLongRange()
{
  return true;
}


double PAcoulombBCFitClass::Vlong(double q, int level)
{
  if (q <= 0.0)
    return 2.0/sqrt(M_PI)*Z1Z2*alpha;
  else 
    return Z1Z2/q*erf(alpha*q);
}

double PAcoulombBCFitClass::Vlong_k(double boxVol, double k, int level)
{
  if (k <= 0.0)
    k = 1.0e-30;
  return 4.0*M_PI/(boxVol*k*k)*exp(-k*k/(4.0*alpha*alpha));
}

double PAcoulombBCFitClass::dVlong(double q, int level)
{
  return 0.0;
}


// void PAcoulombBCFitClass::DoBreakup(const dVec &box, const Array<dVec,1> &kVecs)
// {
//   // Calculate the cutoff parameter
//   double minL, boxVol;
//   boxVol = minL = box[0];
//   for (int i=1; i<NDIM; i++) {
//     minL = min(minL, box[i]);
//     boxVol *= box[i];
//   }
//   alpha = 7.0/minL;

//   // First, subtract off long-ranged part from the bicubic spline to
//   // make them short-ranged.
//   for (int level=0; level<NumBetas; level++) {
//     double tau = SmallestBeta * pow(2.0, level);
//     for (int qi=0; qi<Usplines(level).Nx; qi++) {
//       double q = (*Usplines(level).Xgrid)(qi);
//       double v = Vlong(q, level);
//       for (int ti=0; ti<Usplines(level).Ny; ti++) {
// 	Usplines(level)(qi,ti) -= tau*v;
// 	dUsplines(level)(qi,ti) -= v;
//       }
//     }
//   }


//   // Now, calculate the k-space parts
//   Ulong_k.resize(NumBetas, kVecs.size());
//   dUlong_k.resize(NumBetas, kVecs.size());
//   U_RPA_long_k.resize(NumBetas, kVecs.size());
//   dU_RPA_long_k.resize(NumBetas, kVecs.size());
//   for (int level=0; level<NumBetas; level++) {
//     double tau = SmallestBeta * pow(2.0, level);
//     Ulong_0(level)  = tau*Vlong(0.0, level);
//     dUlong_0(level) = Vlong(0.0, level);
//     for (int ki=0; ki<kVecs.size(); ki++) {
//       double k = sqrt(dot(kVecs(ki), kVecs(ki)));
//       Ulong_k = tau*Vlong_k(boxVol, k, level);
//       dUlong_k = Vlong_k(boxVol, k, level);
//     }
//   }
// }

double PAcoulombBCFitClass::Udiag (double r, int level)
{
  double rend = qgrid->End;

  if (r <= (rend*1.0000001)) 
    return Usplines(level)(r, 0.0);
  else {
    // Coulomb action is independent of z
    double beta = SmallestBeta;
    for (int i=0; i<level; i++)
      beta *= 2.0;
    double V = beta * Z1Z2/r;
    // Now add residual term
    double Uend = Usplines(level)(rend,0.0);
    double Vend = beta * Z1Z2/qgrid->End;
    double delta = Uend - Vend;
    double correction = delta * rend*rend*rend/(r*r*r);
    return (V + correction);
  }
}

double PAcoulombBCFitClass::Udiag_p (double r, int level)
{
  if (r <= (qgrid->End*1.0000001)) 
    return Usplines(level).d_dx(r, 0.0);
  else {
    // Coulomb action is independent of z
    double beta = SmallestBeta;
    for (int i=0; i<level; i++)
      beta *= 2.0;
    return (beta*Pot->dVdr(r));
  }
}

double PAcoulombBCFitClass::Udiag_pp (double r, int level)
{
  if (r <= (qgrid->End*1.0000001)) 
    return Usplines(level).d2_dx2(r, 0.0);
  else {
    // Coulomb action is independent of z
    double beta = SmallestBeta;
    for (int i=0; i<level; i++)
      beta *= 2.0;
    return (beta*Pot->d2Vdr2(r));
  }
}



double PAcoulombBCFitClass::dUdiag (double r, int level)
{
  double rend = qgrid->End;
  if (r <= (rend*1.0000001)) 
    return dUsplines(level)(r, 0.0);
  else {
    double V = Z1Z2/r;
    double dUend = dUsplines(level)(rend, 0.0);
    double Vend  = Z1Z2/rend;
    double delta = dUend-Vend;
    double correction = delta * (rend*rend*rend)/(r*r*r);
    return (V + correction);
  }
}

double PAcoulombBCFitClass::dUdiag_p (double r, int level)
{
  if (r <= (qgrid->End*1.0000001)) 
    return dUsplines(level).d_dx(r, 0.0);
  else // Coulomb action is independent of z
    return (Pot->dVdr(r));
}

double PAcoulombBCFitClass::dUdiag_pp (double r, int level)
{
  if (r <= (qgrid->End*1.0000001)) 
    return dUsplines(level).d2_dx2(r, 0.0);
  else  // Coulomb action is independent of z
    return (Pot->d2Vdr2(r));
}

void PAcoulombBCFitClass::Setrc (double rc)
{
  rcut = rc;
  const int numPoints=1000;

  int numLevels = Usplines.size();
  
  Array<double,1> coefs, error;
  Array<double,1> y(numPoints), sigma(numPoints);
  Array<double,2> F(numPoints,3);
  double deltar = (qgrid->End - rcut)/(numPoints-1);
  
  // Setup basis functions:
  for (int i=0; i<numPoints; i++) {
    double r = rcut + i*deltar;
    F(i,0) = 1.0/r;
    F(i,1) = 1.0/(r*r);
    F(i,2) = 1.0/(r*r*r);
    //    F(i,3) = 1.0/(r*r*r*r);
    sigma(i) = 1.0;
  }

  // Do U fits;
  Ucoefs.resize(3, numLevels);
  for (int level=0; level<numLevels; level++) {
    for (int i=0; i<numPoints; i++) {
      double r = rcut + i*deltar;
      y(i) = Udiag (r, level);
    }
    LinFitSVD (y, sigma, F, coefs, error, 1.0e-15);
    for (int i=0; i<coefs.size(); i++)
      Ucoefs(i, level) = coefs(i);
    //    cerr << "Ucoefs = " << coefs << endl;
  }

  // Do dU fits;
  dUcoefs.resize(3, numLevels);
  for (int level=0; level<numLevels; level++) {
    for (int i=0; i<numPoints; i++) {
      double r = rcut + i*deltar;
      y(i) = dUdiag (r, level);
    }
    LinFitSVD (y, sigma, F, coefs, error, 1.0e-15);
    for (int i=0; i<coefs.size(); i++)
      dUcoefs(i, level) = coefs(i);
  }

  // Do V fits;
  Vcoefs.resize(3);
  for (int i=0; i<numPoints; i++) {
    double r = rcut + i*deltar;
    y(i) = V(r);
  }
  LinFitSVD (y, sigma, F, coefs, error, 1.0e-15);
  for (int i=0; i<coefs.size(); i++)
    Vcoefs(i) = coefs(i);
  Vcoefs = 0.0;
  Vcoefs(0) = Z1Z2; Vcoefs(1) = 0.0; Vcoefs(2) = 0.0;
}

#include <gsl/gsl_sf.h>

double PAcoulombBCFitClass::Xk_U (double k, int level)
{
  double C0 = -4.0*M_PI/(k*k) * cos(k*rcut);
  double C1 =  4.0*M_PI/k * (gsl_sf_Si(k*rcut) - 0.5*M_PI);
  double C2 =  4.0*M_PI/k * (k*gsl_sf_Ci(k*rcut) - sin(k*rcut)/rcut);
  double C3 = -4.0*M_PI/k * (k*cos(k*rcut)/(2.0*rcut) + sin(k*rcut)/(2.0*rcut*rcut) +
			     0.5*k*k*(gsl_sf_Si(k*rcut)-0.5*M_PI));

  return (Ucoefs(0,level)*C0  + Ucoefs(1,level)*C1 + 
	  Ucoefs(2,level)*C2 /*+ Ucoefs(3,level)*C3*/);
}

double PAcoulombBCFitClass::dXk_U_dk  (double k, int level)
{
  double kr = k*rcut;
  double dC0 = 4.0*M_PI/(k*k*k) * (2.0*cos(kr) + kr*sin(kr));
  double dC1 = 4.0*M_PI/(k*k) * (0.5*M_PI + sin(kr) - gsl_sf_Si(kr));
  double dC2 = 4.0*M_PI/(k*k*rcut) * sin(kr);
  double dC3 = 2.0*M_PI*(0.5*M_PI - cos(kr)/kr + sin(kr)/(kr*kr) - gsl_sf_Si(kr));

  return (Ucoefs(0,level)*dC0 + Ucoefs(1,level)*dC1 + 
	  Ucoefs(2,level)*dC2 /*+ Ucoefs(3,level)*dC3*/);
}

double PAcoulombBCFitClass::Xk_dU (double k, int level)
{
  double C0 = -4.0*M_PI/(k*k) * cos(k*rcut);
  double C1 =  4.0*M_PI/k * (gsl_sf_Si(k*rcut) - 0.5*M_PI);
  double C2 =  4.0*M_PI/k * (k*gsl_sf_Ci(k*rcut) - sin(k*rcut)/rcut);
  double C3 = -4.0*M_PI/k * (k*cos(k*rcut)/(2.0*rcut) + sin(k*rcut)/(2.0*rcut*rcut) +
			     0.5*k*k*(gsl_sf_Si(k*rcut)-0.5*M_PI));
  return (dUcoefs(0,level)*C0 + dUcoefs(1,level)*C1 +
	  dUcoefs(2,level)*C2 /*+ dUcoefs(3,level)*C3*/);
}
 
double PAcoulombBCFitClass::Xk_V (double k)
{
  double C0 = -4.0*M_PI/(k*k) * cos(k*rcut);

  return (Z1Z2*C0); 
}


double PAcoulombBCFitClass::Vk (double k)
{
  double C0 = 4.0*M_PI/(k*k);
  return (Z1Z2*C0); 
}
