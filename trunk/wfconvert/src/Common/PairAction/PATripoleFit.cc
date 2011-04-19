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

#include "PATripoleFit.h"

/// The following routines are used only if we are creating fits, not
/// using them.

void PATripoleFitClass::ReadParams(IOSectionClass &inSection)
{
}

void PATripoleFitClass::WriteBetaIndependentInfo (IOSectionClass &outSection)
{ }


void PATripoleFitClass::DoFit (Rho &rho)
{
}


void PATripoleFitClass::Error(Rho &rho, double &Uerror, double &dUerror)
{
  Uerror = 0.0;
  dUerror = 0.0;
}


void PATripoleFitClass::WriteFit (IOSectionClass &outSection)
{
}



double PATripoleFitClass::U(double q, double z, double s2, int level)
{
  double beta = SmallestBeta;
  for (int i=0; i<level; i++)
    beta *= 2.0;
  double r = q + 0.5*z;
  double rp = q - 0.5*z;
  double V = 0.5*Z1Z2*(1.0/(r*r*r) + 1.0/(rp*rp*rp));
  return (beta*V);
}

double PATripoleFitClass::dU(double q, double z, double s2, int level)
{
  double r = q + 0.5*z;
  double rp = q - 0.5*z;
  double V = 0.5*Z1Z2*(1.0/(r*r*rp) + 1.0/(rp*rp*rp));
  return (V);
}


bool PATripoleFitClass::Read (IOSectionClass &in,
			      double smallestBeta, int numBetas)
{
  SmallestBeta = smallestBeta;
  NumBetas=numBetas;
  // Read Particles;
  assert(in.OpenSection("Fits"));
  assert(in.OpenSection("Particle1"));
  Particle1.Read(in);
  in.CloseSection();
  assert(in.OpenSection("Particle2"));
  Particle2.Read(in);
  in.CloseSection();
  lambda = Particle1.lambda + Particle2.lambda;
  assert (lambda == 0.0);

  // Read Potential;
  in.OpenSection("Potential");
  if (!in.ReadVar ("Z1Z2", Z1Z2))
    Z1Z2 = 0.0;
  in.CloseSection();
  in.CloseSection();
  return true;
}



/////////////////////////
/// Long-ranged stuff ///
/////////////////////////
bool PATripoleFitClass::IsLongRange()
{
  return (Z1Z2 != 0.0);
}


/// The diagonal action only -- used for long-range breakup
double PATripoleFitClass::Udiag(double q, int level)
{  
  double beta = SmallestBeta;
  for (int i=0; i<level; i++)
    beta *= 2.0;
  return beta*Z1Z2/(q*q*q);
}

/// The q-derivative of the above
double PATripoleFitClass::Udiag_p(double q, int level) 
{  
  double beta = SmallestBeta;
  for (int i=0; i<level; i++)
    beta *= 2.0;
  return -3.0*beta*Z1Z2/(q*q*q*q);
}

/// The q-derivative of the above
double PATripoleFitClass::Udiag_pp(double q, int level) 
{  
  double beta = SmallestBeta;
  for (int i=0; i<level; i++)
    beta *= 2.0;
  return 12.0*beta*Z1Z2/(q*q*q*q*q);
}

/// The beta-derivative of the diagonal action
double PATripoleFitClass::dUdiag    (double q, int level) 
{
  return Z1Z2/(q*q*q);
}

/// The q-derivative of the above
double PATripoleFitClass::dUdiag_p  (double q, int level) 
{
  return -3.0*Z1Z2/(q*q*q*q);
}

/// The q-derivative of the above
double PATripoleFitClass::dUdiag_pp (double q, int level) 
{
  return 12.0*Z1Z2/(q*q*q*q*q);
}

/// The potential to which this action corresponds.
double PATripoleFitClass::V  (double r) 
{
  return Z1Z2/(r*r*r);
}

/// The q-derivative of the above
double PATripoleFitClass::Vp (double r)
{
  return -3.0*Z1Z2/(r*r*r*r);
}

/// The q-derivative of the above
double PATripoleFitClass::Vpp(double r) 
{
  return 12.0*Z1Z2/(r*r*r*r*r);
}

void PATripoleFitClass::Setrc(double rc)
{
  rcut = rc;
}

#include <gsl/gsl_sf.h>

double PATripoleFitClass::Xk_U(double k, int level)
{
  double beta = SmallestBeta;
  for (int i=0; i<level; i++)
    beta *= 2.0;
  return 4.0*M_PI*beta*Z1Z2/k*(k*gsl_sf_Ci(k*rcut)-sin(k*rcut)/rcut);
}

double PATripoleFitClass::dXk_U_dk(double k, int level)
{
  double beta = ldexp(SmallestBeta, level);
  return 4.0*M_PI*beta*Z1Z2/(k*k*rcut)*sin(k*rcut);
}

double PATripoleFitClass::Xk_dU(double k, int level)
{
  return 4.0*M_PI*Z1Z2/k*(k*gsl_sf_Ci(k*rcut)-sin(k*rcut)/rcut);
}

double PATripoleFitClass::Xk_V(double k)
{
  return 4.0*M_PI*Z1Z2/k*(k*gsl_sf_Ci(k*rcut)-sin(k*rcut)/rcut);
}
