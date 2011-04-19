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

#include "ToppHopfield.h"

double 
ToppHopfieldPot::V(double r)
{
  if (r < rc)
    return (V0*cos(a*r)+b);
  else
    return (-Z/r);
}

double
ToppHopfieldPot::dVdr(double r)
{
  if (r < rc)
    return (-V0*a*sin(a*r));
  else
    return (Z/(r*r));
}

double
ToppHopfieldPot::d2Vdr2(double r)
{
  if (r < rc)
    return (-V0*a*a*cos(a*r));
  else
    return (-2.0*Z/(r*r*r));
}

void
ToppHopfieldPot::Read(IOSectionClass &in)
{
  assert (in.ReadVar("Z",  Z));
  assert (in.ReadVar("V0", V0));
  assert (in.ReadVar("a",  a));
  assert (in.ReadVar("b",  b));
  assert (in.ReadVar("rc", rc));
}

void
ToppHopfieldPot::Write(IOSectionClass &out)
{
  out.WriteVar("Z",  Z);
  out.WriteVar("Type", "ToppHopfield");
  out.WriteVar("V0", V0);
  out.WriteVar("a",  a);
  out.WriteVar("b",  b);
  out.WriteVar("rc", rc);
}
