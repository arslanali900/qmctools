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

#ifndef RADIAL_WF_H
#define RADIAL_WF_H

#include "../PH/Potential.h"
#include "../Splines/CubicSpline.h"
#include "../IO/IO.h"

class RadialWF
{
private:
  int TurningIndex();
  double IntegrateInOut(int &tindex);
  void OriginBC(double r0, double &u0, double &du0);
  void InfinityBC(double rend, double &uend, double &duend);
  Array<Vec2,1> uduVec;
  Array<double,1> normVec;
  Potential *pot;  
public:
  Grid *grid;
  CubicSpline u, dudr;
  int n, l, CoreNodes;
  double Energy, Occupancy, Weight;
  string Label;

  inline double NormDeriv(double r, double u);
  inline Vec2 PseudoDerivs (double r, Vec2 &u_and_du);
  inline Vec2 NormalDerivs (double r, Vec2 &u_and_du);
  inline Vec2 ScalarRelDerivs (double r, Vec2 &u_and_du);
  int CountNodes();
  void IntegrateOut();
  double PartialNorm();
  double LogDerivative();
  void Solve (double tolerance=1.0e-8);
  void Normalize();
  void SetGrid(Grid *newgrid);
  void SetPotential (Potential *newPot);
  Potential *GetPotential ();
  void Write (IOSectionClass &out);
  void Read  (IOSectionClass &in);
};



// These are wrappers so we can inline integration routines

/// Derivative used for normalization of the radial function
class NormalizeDeriv
{
private:
  RadialWF &WF;
public:
  inline double operator()(double r, double u)
  { return WF.NormDeriv(r, u); }
  NormalizeDeriv (RadialWF &wf) : WF(wf)
  { /* Do nothing */ }
};

/// Derivatives for the radial equation given a pseudoHamiltonian  
class PHDerivs 
{
private:
  RadialWF &WF;
public:
  inline Vec2 operator() (double r, Vec2& u_and_du)
  { return WF.PseudoDerivs(r, u_and_du); }
  PHDerivs (RadialWF &wf) : WF(wf)
  { /* Do nothing */ }
};

/// Derivatives for the radial equation given a pseudoHamiltonian  
class RegularDerivs 
{
private:
  RadialWF &WF;
public:
  inline Vec2 operator() (double r, Vec2& u_and_du)
  { return WF.NormalDerivs(r, u_and_du); }
  RegularDerivs (RadialWF &wf) : WF(wf)
  { /* Do nothing */ }
};



/// Derivatives for the scalar relativistic radial equation
class NonPHDerivs
{
private:
  RadialWF &WF;
public:
  inline Vec2 operator() (double r, Vec2& u_and_du)
  { return WF.ScalarRelDerivs(r, u_and_du); }
  NonPHDerivs (RadialWF &wf) : WF(wf)
  { /* Do nothing */ }
};


/////////////////////////////////////////////////////////////////
//                     Inline functions                        //
/////////////////////////////////////////////////////////////////

inline double RadialWF::NormDeriv (double r, double cumulative)
{
  double uval = u(r);
  return (uval*uval);
}

inline Vec2 RadialWF::PseudoDerivs (double r, Vec2 &u_and_du)
{
  Vec2 derivs;  
  derivs[0] = u_and_du[1];
  double A = pot->A(r);
  double B = pot->B(r);
  double V = pot->V(l,r);
  double dAdr = pot->dAdr(r);
  derivs[1] =  1.0/A*
    (-dAdr*u_and_du[1] + (dAdr/r + (double)(l*(l+1))*B/(r*r) 
			  + 2.0*(V-Energy))*u_and_du[0]);
  return derivs;
}


inline Vec2 RadialWF::NormalDerivs(double r, Vec2 &u_and_du)
{
  Vec2 derivs;
  double V = pot->V(l,r);
  derivs[0] = u_and_du[1];
  derivs[1] = ((double)(l*(l+1))/(r*r) + 2.0*(V-Energy))*u_and_du[0];
  return derivs;
}


inline Vec2 RadialWF::ScalarRelDerivs (double r, Vec2 &u_and_du)
{
  const double alpha = 1.0/137.036;
  const double kappa = -1.0;

  Vec2 derivs;
  derivs[0] = u_and_du[1];
  double V = pot->V(l,r);
  double dVdr = pot->dVdr(r);
  double M = 1.0 - alpha*alpha*0.5*(V-Energy);
  
  derivs[1] = ((double)(l*(l+1))/(r*r) + 2.0*M*(V-Energy))*u_and_du[0] 
    - 0.5*alpha*alpha/M*dVdr*(u_and_du[1] + u_and_du[0]*kappa/r);
  
  return derivs;
}

#endif
