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

#include "QuinticPH.h"

bool QuinticPH::IsPH()
{ 
  /// Check to see if we really have position-dependent masses.
  bool isph = false;
  for (double y=0.0; y<=1.0; y+=0.001) {
    double r = CoreRadius * y;
    isph = isph || (fabs(A(r)-1.0)>1.0e-6);
    isph = isph || (fabs(B(r)-1.0)>1.0e-6);
  }
  return isph;
}

double QuinticPH::GetCoreRadius()
{
  return (CoreRadius);
}

double QuinticPH::A(double r)
{
  double val;
  if (r < CoreRadius) {
    double pa = pA(r);
    val = pa*pa + ABmin;
  }
  else
    val = 1.0;
  return val;
}


double QuinticPH::dAdr(double r)
{
  double val;
  if (r < CoreRadius) 
    val = 2.0*pA(r)*pA.Deriv(r);
  else
    val = 0.0;
  return val;
}


double QuinticPH::d2Adr2(double r)
{
  double val;
  if (r < CoreRadius) {
    double dpa = pA.Deriv(r);
    val = 2.0*(dpa*dpa + pA(r)*pA.Deriv2(r));
  }
  else
    val = 0.0;
  return val;
}


double QuinticPH::B(double r)
{
  double val;
  if (r < CoreRadius) {
    double pb = pB(r);
    val = pb*pb + ABmin;
  }
  else
    val = 1.0;
  return val;
}

double QuinticPH::V(double r)
{
  double val;
  if (r < CoreRadius && UseVcore)
    val = Vcore(r);
  else
    val = Vouter->V(r);
  return val;
}


double QuinticPH::dVdr(double r)
{
  double val;
  if (r < CoreRadius && UseVcore)
    val = Vcore.Deriv(r);
  else
    val = Vouter->dVdr(r);
  return val;
}


double QuinticPH::d2Vdr2(double r)
{
  double val;
  if (r < CoreRadius && UseVcore)
    val = Vcore.Deriv2(r);
  else
    val = Vouter->d2Vdr2(r);
  return val;
}


void QuinticPH::WriteWithoutVouter(IOSectionClass &out)
{
  out.WriteVar ("Type", "QuinticPH");
  out.NewSection("Agrid");  Agrid.Write(out);  out.CloseSection();
  out.NewSection("Bgrid");  Bgrid.Write(out);  out.CloseSection();
  out.NewSection("Vgrid");  Vgrid.Write(out);  out.CloseSection();

  out.WriteVar ("ABmin", ABmin);
  out.WriteVar ("pA", pA.Data());
  out.WriteVar ("pB", pB.Data());
  out.WriteVar ("Vcore", Vcore.Data());
}

void QuinticPH::Write(IOSectionClass &out)
{
  WriteWithoutVouter(out);
  out.NewSection ("Vouter");
  Vouter->Write(out);
  out.CloseSection();
  out.WriteVar ("UseVcore", UseVcore);
}

/// Note:  Vouter must be set before calling "Read".
void QuinticPH::Read(IOSectionClass &in)
{
  assert (in.OpenSection ("Agrid")); Agrid.Read(in); in.CloseSection();
  assert (in.OpenSection ("Bgrid")); Bgrid.Read(in); in.CloseSection();
  assert (in.OpenSection ("Vgrid")); Vgrid.Read(in); in.CloseSection();
  CoreRadius = Agrid.End;

  assert(in.ReadVar("ABmin", ABmin));
  Array<double,1> temp;
  assert (in.ReadVar ("pA", temp));
  pA.Init (&Agrid, temp, 0.0, 0.0, 0.0, 0.0);
  assert (in.ReadVar ("pB", temp));
  pB.Init (&Bgrid, temp, 0.0, 0.0, 0.0, 0.0);

  bool success = in.OpenSection ("Vouter");
  double dVend, d2Vend;
  if (success) {
    Vouter = ReadPotential(in);
    in.CloseSection();
    dVend  = Vouter->dVdr(CoreRadius);
    d2Vend = Vouter->d2Vdr2(CoreRadius);
  }
  in.ReadVar ("UseVcore", UseVcore);
//   cerr << "UseVcore = " << ( UseVcore ? "true" : "false") << endl;
//   cerr << "sizeof(bool) = " << sizeof(bool) << endl;

  assert (in.ReadVar ("Vcore", temp));
  Vcore.Init (&Vgrid, temp, NAN, dVend, NAN, d2Vend);

  CoreRadius = Agrid.End;
}

