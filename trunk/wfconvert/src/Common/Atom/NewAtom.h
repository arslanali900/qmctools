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

#ifndef NEW_ATOM_H
#define NEW_ATOM_H

#include "DFTAtom.h"

inline Atom* ReadAtom (IOSectionClass &in)
{
  Atom *atom;
  string type;
  assert (in.ReadVar ("Type", type));
  if (type == "DFT")
    atom = new DFTAtom;
  else {
    cerr << "Unrecognized atom type \"" << type 
	 << "\" in ReadAtom.  Exitting.\n";
    exit(1);
  }
  atom->Read(in);
  return (atom);
}

#endif
