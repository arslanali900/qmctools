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

#ifndef QUINTIC_PH_H
#define QUINTIC_PH_H

#include "PotentialBase.h"
#include "../Splines/QuinticSpline.h"

class QuinticPH : public Potential
{
private:
  double CoreRadius;
  // The minimum value of A and B
  double ABmin;
  inline void Copy (const QuinticPH &ph);
public:
  QuinticSpline pA, pB, Vcore;
  LinearGrid Agrid, Bgrid, Vgrid;
  Potential *Vouter;
  // This is true if this represents a bare potential.  Otherwise, it
  // represents a total potential:  Total = Bare + XC
  bool IsBare;
  // True if we want to use Vcore spline
  bool UseVcore;

  bool IsPH();

  double GetCoreRadius ();
  inline void   SetCoreRadius (double coreRadius);

  // These are used to get and set the minimum value of A and B,
  // adjusting PA and PB accordingly.
  inline double GetABmin();
  inline void SetABmin(double newmin);

  // These are used to get and set the values of A, B and V
  inline int    GetNumAParams();
  inline void   SetNumAParams(int num);
  inline int    GetNumBParams();
  inline void   SetNumBParams(int num);
  inline int    GetNumVParams();
  inline void   SetNumVParams(int num);
  inline double GetAval (int i);
  inline void   SetAval (int i, double Aval, bool isNegative);
  inline double GetBval (int i);
  inline void   SetBval (int i, double Bval, bool isNegative);
  inline double GetVval (int i);
  inline void   SetVval (int i, double val);

  double A      (double r);
  double dAdr   (double r);
  double d2Adr2 (double r);
  double B      (double r);
  double V      (double r);
  double dVdr   (double r);
  double d2Vdr2 (double r);
  inline double  GetParam (int i) const;
  inline void SetParam (int i, double val);
  inline int NumParams();
  inline void Init (Potential *outer, double coreRadius,
		    double abmin, int numA, int numB, int numV);
  void Read (IOSectionClass &in);
  void WriteWithoutVouter (IOSectionClass &out);
  void Write (IOSectionClass &out); 
  inline QuinticPH& operator=(const QuinticPH &ph);

  QuinticPH() : ABmin(0.0), UseVcore(true)
  { /* do nothing for now */ }
};


inline QuinticPH& QuinticPH::operator=(const QuinticPH &ph)
{
  Copy(ph);
  return (*this);
}

inline void QuinticPH::Copy (const QuinticPH &ph)
{
  CoreRadius = ph.CoreRadius;
  ABmin = ph.ABmin;
  Agrid.Init (0.0, ph.CoreRadius, ph.Agrid.NumPoints);
  Bgrid.Init (0.0, ph.CoreRadius, ph.Bgrid.NumPoints);
  Vgrid.Init (0.0, ph.CoreRadius, ph.Vgrid.NumPoints);
  Vouter = ph.Vouter;
  IsBare = ph.IsBare;
  UseVcore = ph.UseVcore;

  Array<double,1> Aparams(ph.pA.NumPoints());
  Aparams = ph.pA.Data();
  pA.Init (&Agrid, Aparams, 0.0, 0.0, 0.0, 0.0);

  Array<double,1> Bparams(ph.pB.NumPoints());
  Bparams = ph.pB.Data();
  pB.Init (&Bgrid, Bparams, 0.0, 0.0, 0.0, 0.0);

  Array<double,1> Vparams(ph.Vcore.NumPoints());
  Vparams = ph.Vcore.Data();
  double Vend = Vouter->V(CoreRadius);
  double dVend = Vouter->dVdr(CoreRadius);
  double d2Vend = Vouter->d2Vdr2(CoreRadius);
  Vcore.Init (&Vgrid, Vparams, NAN, dVend, NAN, d2Vend);
}

inline void QuinticPH::SetCoreRadius(double coreRadius)
{
  CoreRadius = coreRadius;
  Array<double,1> Atmp, Btmp, Vtmp;
  if (pA.NumPoints() != 0) {
    Atmp.resize(   pA.NumPoints());  Atmp  = pA.Data();
    Btmp.resize(   pB.NumPoints());  Btmp  = pB.Data();
    Vtmp.resize(Vcore.NumPoints());  Vtmp = Vcore.Data();
  }
  else {
    Atmp.resize(4);  Atmp = sqrt(1.0-ABmin);
    Btmp.resize(4);  Btmp = sqrt(1.0-ABmin);
    Vtmp.resize(4);  Vtmp = 0.0;
  }
  Agrid.Init (0.0, CoreRadius, Atmp.size());
  Bgrid.Init (0.0, CoreRadius, Btmp.size());
  Vgrid.Init (0.0, CoreRadius, Vtmp.size());
  pA.Init (&Agrid, Atmp, 0.0, 0.0, 0.0, 0.0);
  pB.Init (&Bgrid, Btmp, 0.0, 0.0, 0.0, 0.0);
  double dVend = Vouter->dVdr(CoreRadius);
  double d2Vend = Vouter->d2Vdr2(CoreRadius);
  Vcore.Init(&Vgrid, Vtmp, NAN, dVend, NAN, d2Vend);
}

inline double QuinticPH::GetParam(int i) const
{
  // A is constrained at core radius
  int numA = pA.NumPoints()-1;
  // B is constrained at origin at core radius
  int numB = pB.NumPoints()-2;
  // V is constrained at core radius
  int numV = Vcore.NumPoints()-1;
  if (i < numA)
    return (pA(i));
  else if (i < (numA+numB))
    return (pB(i-numA+1));
  else 
    return (Vcore(i-numA-numB));
}

inline void QuinticPH::SetParam(int i, double val)
{
  // A is constrained at core radius
  int numA = pA.NumPoints()-1;
  // B is constrained at origin at core radius
  int numB = pB.NumPoints()-2;
  // V is constrained at core radius
  int numV = Vcore.NumPoints()-1;

  // B(0) is locked to be A(0)
  if (i == 0) {
    pA(0) = val; 
    pB(0) = val;
  }
  else if (i < numA)
    pA(i) = val;
  else if (i < (numA+numB))
    pB(i-numA+1) = val;
  else 
    Vcore(i-numA-numB) = val;
}

    
inline double QuinticPH::GetAval (int i)
{
  return (pA(i)*pA(i) + ABmin);
}

inline void QuinticPH::SetAval (int i, double Aval, bool isNegative)
{
  double pAval = sqrt (Aval - ABmin);
  if (isNegative)
    pAval = -pAval;
  if (i == 0)
    pB(0) = pAval;
  pA(i) = pAval;
}


inline double QuinticPH::GetBval (int i)
{
  return (pB(i+1)*pB(i+1) + ABmin);
}

inline void QuinticPH::SetBval (int i, double Bval, bool isNegative)
{
  double pBval = sqrt (Bval - ABmin);
  if (isNegative)
    pBval = -pBval;

  pB(i+1) = pBval;
}

inline double QuinticPH::GetVval (int i)
{
  return (Vcore(i));
}

inline void QuinticPH::SetVval (int i, double val)
{
  Vcore(i) = val;
}

inline int QuinticPH::GetNumAParams()
{
  return (pA.NumPoints()-1);
}

inline void QuinticPH::SetNumAParams(int num)
{
  int oldNum = pA.NumPoints();
  Array<double,1> newpA(num+1);
  for (int i=0; i<min(oldNum,num+1); i++)
    newpA(i) = pA(i);
  for (int i=oldNum; i<num+1; i++)
    newpA(i) = sqrt(1.0-ABmin);
  Agrid.Init(0.0, CoreRadius, num+1);
  pA.Init(&Agrid, newpA, 0.0, 0.0, 0.0, 0.0);
}

inline int QuinticPH::GetNumBParams()
{
  return (pB.NumPoints()-2);
}

inline void QuinticPH::SetNumBParams(int num)
{
  int oldNum = pB.NumPoints();
  Array<double,1> newpB(num+2);
  for (int i=0; i<min(oldNum,num+2); i++)
    newpB(i) = pB(i);
  for (int i=oldNum; i<num+2; i++)
    newpB(i) = sqrt(1.0-ABmin);
  Bgrid.Init(0.0, CoreRadius, num+2);
  pB.Init(&Bgrid, newpB, 0.0, 0.0, 0.0, 0.0);
}

inline int QuinticPH::GetNumVParams()
{
  return (Vcore.NumPoints()-1);
}

inline void QuinticPH::SetNumVParams(int num)
{
  int oldNum = Vcore.NumPoints();
  cerr << "CoreRadius = " << CoreRadius << endl;
  double Vend = Vouter->V(CoreRadius);
  cerr << "Vend = " << Vend << endl;
  double dVend = Vouter->dVdr(CoreRadius);
  double d2Vend = Vouter->d2Vdr2(CoreRadius);
  Array<double,1> newV(num+1);
  for (int i=0; i<min(oldNum, num+1); i++)
    newV(i) = Vcore(i);
  for (int i=oldNum; i<num+1; i++)
    newV(i) = Vend;
  newV(num) = Vend;
  Vgrid.Init(0.0, CoreRadius, num+1);
  Vcore.Init(&Vgrid, newV, NAN, dVend, NAN, d2Vend);
}

inline double QuinticPH::GetABmin()
{
  return ABmin;
}

inline void QuinticPH::SetABmin(double newmin)
{
  for (int i=0; i<pA.NumPoints(); i++) {
    double val = pA(i)*pA(i)+ABmin;
    if (pA(i) < 0.0)
      pA(i) = (val>newmin) ? -sqrt(val-newmin) : 0.0;
    else
      pA(i) = (val>newmin) ?  sqrt(val-newmin) : 0.0;
  }
  for (int i=0; i<pB.NumPoints(); i++) {
    double val = pB(i)*pB(i)+ABmin;
    if (pB(i) < 0.0)
      pB(i) = (val>newmin) ? -sqrt(val-newmin) : 0.0;
    else
      pB(i) = (val>newmin) ?  sqrt(val-newmin) : 0.0;
  }
  ABmin = newmin;
}

void QuinticPH::Init (Potential *outer, double coreRadius,
		      double abmin, int numA, int numB, int numV)
{
  Vouter = outer;
  CoreRadius = coreRadius;
  ABmin = abmin;
  Agrid.Init (0.0, CoreRadius, numA+1);
  Bgrid.Init (0.0, CoreRadius, numB+2);
  Vgrid.Init (0.0, CoreRadius, numV+1);
  Array<double,1> initA(numA+1);
  Array<double,1> initB(numB+2);
  Array<double,1> initV(numV+1);
  initA = sqrt(1.0-ABmin);
  initB = sqrt(1.0-ABmin);
  initV = 0.0;
  double Vend = Vouter->V(CoreRadius);
  cerr << "Vend = " << Vend << endl;
  double dVend = Vouter->dVdr(CoreRadius);
  double d2Vend = Vouter->d2Vdr2(CoreRadius);
  initV(numV) = Vend;
  pA.Init (&Agrid, initA, 0.0, 0.0, 0.0, 0.0);
  pB.Init (&Bgrid, initB, 0.0, 0.0, 0.0, 0.0);
  Vcore.Init (&Vgrid, initV, NAN, dVend, NAN, d2Vend);
}

  


#endif
