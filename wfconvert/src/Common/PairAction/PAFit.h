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

#ifndef PA_FIT_H
#define PA_FIT_H

#include "PAcoulombBCFit.h"
#include "PAMonopoleFit.h"
#include "PADipoleFit.h"
#include "PATripoleFit.h"
#include "PAclassicalFit.h"
#include "PAtricubicFit.h"
#include "PAtricubicFit2.h"
#include "PAzeroFit.h"
#include "DavidPAClass.h"

inline PairActionFitClass *ReadPAFit (IOSectionClass &in, 
				      double smallestBeta, int numLevels)
{
  assert (in.OpenSection("Fits"));
  string type;
  assert (in.ReadVar("Type", type));
  in.CloseSection (); // "Fits"
  PairActionFitClass *fit;
  //  if (type == "szfit")
  //  fit = new PAszFitClass;
  //else if (type == "coulombfit")
  //  fit = new PAcoulombFitClass;
  /*else*/ if (type == "coulombBCfit")
    fit = new PAcoulombBCFitClass;
  else if (type == "classical")
    fit = new PAclassicalFitClass;
  else if (type == "coulomb")
    fit = new PACoulombFitClass;
  else if (type == "dipole")
    fit = new PADipoleFitClass;
  else if (type == "tripole")
    fit = new PATripoleFitClass;
  //else if (type == "sfit")
  //  fit = new PAsFitClass;
  else if (type == "tricubicfit")
    fit = new PAtricubicFitClass;
  else if (type == "tricubicfit2")
    fit = new PAtricubicFit2Class;
  else if (type == "zerofit")
    fit=new PAzeroFitClass;
  else if (type=="DavidFit")
    fit = new DavidPAClass;
  else {
    cerr << "Unrecognize pair action fit type \"" 
	 << type << "\".  Exitting.\n";
    exit(1);
  }
  fit->Read(in, smallestBeta, numLevels);
  return (fit);
}


#endif
