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

#ifndef ENERGIES_H
#define ENERGIES_H

#include "Atom.h"


scalar KEderiv(scalar r, scalar E, void *WFptr);
scalar PEderiv(scalar r, scalar E, void *Atomptr);
scalar HEderiv(scalar r, scalar E, void *Atomptr);
scalar XCEderiv(scalar r, scalar E, void *Atomptr);
scalar NUCderiv(scalar r, scalar E, void *Atomptr);

#endif
