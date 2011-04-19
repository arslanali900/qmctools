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
#include "../Integration/GKIntegration.h"
#include "../Fitting/Fitting.h"
#include <gsl/gsl_sf_expint.h>


void kSpacePH::CalcTailCoefs (double r1, double r2)
{
  R1 = r1; R2 = r2;
  const int numPoints=1000;
  double nInv = 1.0/(double)(numPoints-1);

  Array<double,1> coefs, error;
  Array<double,1> y(numPoints), sigma(numPoints);
  Array<double,2> F(numPoints,3);
  
  // Setup basis functions:
  for (int i=0; i<numPoints; i++) {
    double r = r1 + (r2-r1)*(double)i*nInv;
    F(i,0) = 1.0/r;
    F(i,1) = 1.0/(r*r);
    F(i,2) = 1.0/(r*r*r);
    sigma(i) = 1.0;
    y(i) = PH.V(r);
    sigma(i) = 1.0;
  }
  LinFitSVD (y, sigma, F, coefs, error, 1.0e-15);
  Ctail1 = coefs(0);
  Ctail2 = coefs(1);
  Ctail3 = coefs(2);
  HaveTailCoefs = true;
}


class aIntegrand
{
protected:
  Potential &PH;
  double k, kinv;
public:
  inline void Setk (double newk) 
  {    k = newk;    kinv = 1.0/newk; }
   
  inline double operator()(double r)
  { 
    if (k < 1.0e-12)
      return 4.0*M_PI*(PH.A(r)-1.0)*r*r;
    else
      return 4.0*M_PI*(PH.A(r)-1.0)*r*sin(k*r)*kinv; 
  }

  aIntegrand (Potential &ph) : PH(ph) 
  { /* do nothing for now */ }
};

double kSpacePH::a(double k)
{
  aIntegrand integrand(PH);
  integrand.Setk(k);
  GKIntegration<aIntegrand,GK31> integrator(integrand);
  //  integrator.SetRelativeErrorMode();
  return (integrator.Integrate(0.0, PH.GetCoreRadius(), 1.0e-10));
}


class bPerpIntegrand
{
protected:
  Potential &PH;
  double k, kinv;
  const double Third, TwoThirds;
public:
  inline void Setk (double newk) 
  {    k = newk;    kinv = 1.0/newk; }
   
  inline double operator()(double r)
  { 
    double kr = k*r;
    double krInv = 1.0/kr;
    double b = PH.B(r) - 1.0;
    double j0, j2;
    if (kr < 1.0e-10) {
      j0 = 1.0;
      j2 = 0.0;
    }
    else {
      j0 = sin(kr)/kr;
      j2 = (3.0*krInv*krInv*krInv-krInv)*sin(kr) - 
      3.0*krInv*krInv*cos(kr);
    }
    return (4.0*M_PI*r*r*b*(TwoThirds*j0-Third*j2));
  }
  
  bPerpIntegrand (Potential &ph) : PH(ph), Third(1.0/3.0), TwoThirds(2.0/3.0)
  { /* do nothing for now */ }
};


double kSpacePH::bPerp(double k)
{
  bPerpIntegrand integrand(PH);
  integrand.Setk(k);
  GKIntegration<bPerpIntegrand,GK31> integrator(integrand);
  //  integrator.SetRelativeErrorMode();
  return (integrator.Integrate(0.0, PH.GetCoreRadius(), 1.0e-10));
}



class bParIntegrand
{
protected:
  Potential &PH;
  double k, kinv;
  const double Third, TwoThirds;
public:
  inline void Setk (double newk) 
  {    k = newk;    kinv = 1.0/newk; }
   
  inline double operator()(double r)
  { 
    double kr = k*r;
    double krInv = 1.0/kr;
    double b = PH.B(r) - 1.0;
    double j0, j2;
    if (kr < 1.0e-10) {
      j0 = 1.0;
      j2 = 0.0;
    }
    else {
      j0 = sin(kr)/kr;
      j2 = (3.0*krInv*krInv*krInv-krInv)*sin(kr) - 
      3.0*krInv*krInv*cos(kr);
    }
    return (4.0*M_PI*r*r*b*(TwoThirds*j0+TwoThirds*j2));
  }
  
  bParIntegrand (Potential &ph) : PH(ph), Third(1.0/3.0), TwoThirds(2.0/3.0)
  { /* do nothing for now */ }
};


double kSpacePH::bPar(double k)
{
  bParIntegrand integrand(PH);
  integrand.Setk(k);
  GKIntegration<bParIntegrand,GK31> integrator(integrand);
  //  integrator.SetRelativeErrorMode();
  return (integrator.Integrate(0.0, PH.GetCoreRadius(), 1.0e-10));
}


class VIntegrand
{
protected:
  Potential &PH;
  double k, kinv;
public:
  inline void Setk (double newk) 
  {    k = newk;    kinv = 1.0/newk; }
   
  inline double operator()(double r)
  { return 4.0*M_PI*PH.V(r)*r*sin(k*r)*kinv; }

  VIntegrand (Potential &ph) : PH(ph) 
  { /* do nothing for now */ }
};

class V0Integrand
{
protected:
  Potential &PH;
  double Z;
public:
  inline double operator()(double r)
  { return 4.0*M_PI*(PH.V(r)*r*r+Z*r); }

  V0Integrand (Potential &ph, double rc) : PH(ph)
  { 
    Z = -round(rc*PH.V(rc));
  }
};


double kSpacePH::Vk (double k)
{
  if (fabs(k) < 1.0e-10) {
    V0Integrand integrand(PH, R1);
    GKIntegration<V0Integrand,GK31> integrator(integrand);
    return integrator.Integrate(0.0, R1, 1.0e-10);
    //    return 0.0;
  }
  // First, do the part of the integral up to R1 numerically
  VIntegrand integrand(PH);
  integrand.Setk(k);
  ///////////
  // DEBUG //
  ///////////
//   if (k < 0.4) {
//     FILE *fout = fopen ("Vint.dat", "w");
//     for (double r=0; r<R1; r+=0.001)
//       fprintf (fout, "%1.12e %1.12e %1.12e\n", k, r, integrand(r));
//     fclose (fout);
//   }

  GKIntegration<VIntegrand,GK31> integrator(integrand);  
  //integrator.SetRelativeErrorMode();
  double result = integrator.Integrate(0.0, R1, 1.0e-10);

  // Now, do the remaining part up to infinity using analytic
  // integratio of our fitted form.
  assert (HaveTailCoefs);
  
  result += 4.0*M_PI*Ctail1/(k*k) * cos(k*R1);
  result -= 4.0*M_PI*Ctail2/k*(gsl_sf_Si(k*R1) - 0.5*M_PI);
  result -= 4.0*M_PI*Ctail3/k*(k*gsl_sf_Ci(k*R1) - sin (k*R1)/R1);
//   fprintf (stderr, "k = %1.5e, Vk = %1.12e\n",
// 	   k, result);
  return result;
}


TinyMatrix<double,3,3> kSpacePH::Ftensor (Vec3 deltaG)
{
  double Gmag = sqrt(dot(deltaG, deltaG));
  Vec3 g; 
  if (fabs(Gmag) < 1.0e-13)
    g = Vec3 (1.0, 0.0, 0.0);
  else
    g = deltaG / Gmag;

  double aval, bPerpval, bParval, Vval;
  if (UseCache)
    Cache.GetVals (Gmag, aval, bPerpval, bParval, Vval);
  else {
    aval = a(Gmag);
    bPerpval = bPerp(Gmag);
    bParval = bPar(Gmag);
  }

  TinyMatrix<double,3,3> F, G;
  F(0,0)=0.0;   F(0,1)=0.0; F(0,2)=0.0;
  F(1,0)=0.0;   F(1,1)=0.0; F(1,2)=0.0;  
  F(2,0)=0.0;   F(2,1)=0.0; F(2,2)=0.0;
  for (int i=0; i<3; i++)
    F(i,i) += aval + bPerpval;
  for (int i=0; i<3; i++)
    for (int j=0; j<3; j++)
      G(i,j) = g[i]*g[j];
  F = F + (bParval - bPerpval)*G;
  return F;
}

TinyMatrix<double,3,3> 
kSpacePH::Ftensor (Vec3 deltaG, double aVal, double bPerpVal,
				double bParVal)
{
  Vec3 g;
  double Gmag = sqrt (dot(deltaG, deltaG));
  if (fabs(Gmag) < 1.0e-13)
    g = Vec3 (1.0, 0.0, 0.0);
  else
    g = deltaG / Gmag;

  TinyMatrix<double,3,3> F, G;
  F(0,0)=aVal+bPerpVal;   F(0,1)=0.0;           F(0,2)=0.0;
  F(1,0)=0.0;             F(1,1)=aVal+bPerpVal; F(1,2)=0.0;  
  F(2,0)=0.0;             F(2,1)=0.0;           F(2,2)=aVal+bPerpVal;
  G(0,0)=g[0]*g[0]; G(0,1)=g[0]*g[1]; G(0,2)=g[0]*g[2];
  G(1,0)=g[1]*g[0]; G(1,1)=g[1]*g[1]; G(1,2)=g[1]*g[2];
  G(2,0)=g[2]*g[0]; G(2,1)=g[2]*g[1]; G(2,2)=g[2]*g[2];
  F = F + (bParVal - bPerpVal)*G;

  return F;
}

void
kSpacePH::GetVals (double dGmag, double &aVal, 
		   double &bPerpVal, double &bParVal, double &VVal)
{
  map<double,kMapPoint,FuzzyLess>::iterator iter=kMap.find(dGmag);
  if (iter == kMap.end()) {
    aVal     = a(dGmag);
    bPerpVal = bPerp(dGmag);
    bParVal  = bPar (dGmag);
    VVal     = Vk (dGmag);
    kMap[dGmag] = kMapPoint(aVal, bPerpVal, bParVal, VVal);
    //    cerr << "Adding new val.  Map size = " << kMap.size() << endl;
  }
  else {
    aVal     = (*iter).second.a;
    bPerpVal = (*iter).second.bPerp;
    bParVal  = (*iter).second.bPar;
    VVal     = (*iter).second.V;
  }
}

double 
kSpacePH::V (double deltaGmag)
{
  double aval, bPerpval, bParval, Vval;
  if (UseCache)
    Cache.GetVals (deltaGmag, aval, bPerpval, bParval, Vval);
  else
    Vval = Vk(deltaGmag);
  return Vval;
}
		   

double kSpacePH::V (Vec3 k, Vec3 G, Vec3 Gp)
{
  Vec3 deltaG = G-Gp;
  double Gmag = sqrt(dot(deltaG, deltaG));

  double aval, bPerpval, bParval, Vval;
  if (UseCache) 
    Cache.GetVals (Gmag, aval, bPerpval, bParval, Vval);
  else
    Vval = Vk(Gmag);

  TinyMatrix<double,3,3> F = Ftensor (deltaG);
//   if (G == Gp)
//     cerr << "F = " << F << endl;
  
  Vec3 Gk = G + k;
  Vec3 Gpk = Gp + k;
  Vec3 FGpk (0.0, 0.0, 0.0);
  FGpk[0] = F(0,0)*Gpk[0] + F(0,1)*Gpk[1] + F(0,2)*Gpk[2];
  FGpk[1] = F(1,0)*Gpk[0] + F(1,1)*Gpk[1] + F(1,2)*Gpk[2];
  FGpk[2] = F(2,0)*Gpk[0] + F(2,1)*Gpk[1] + F(2,2)*Gpk[2];
  
  return 0.5*dot(Gk, FGpk) + Vval;
}


void kCache::GetVals(double k, double &a, double &bPerp, double &bPar,
		     double &V)
{
  int i=0;
  while ((i<Cache.size()) && (fabs(Cache[i].k-k) > 1.0e-12))
    i++;
  if (i<Cache.size()) { // It's already in the cache
    a     = Cache[i].a;
    bPerp = Cache[i].bPerp;
    bPar  = Cache[i].bPar;
    V     = Cache[i].V;
  }
  else {
    a     = kPH.a(k);
    bPerp = kPH.bPerp(k);
    bPar  = kPH.bPar (k);
    V     = kPH.Vk    (k);
    kCachePoint newPoint;
    newPoint.a = a;
    newPoint.bPerp = bPerp;
    newPoint.bPar  = bPar;
    newPoint.V     = V;
    newPoint.k = k;
    Cache.push_back(newPoint);
  }
}
