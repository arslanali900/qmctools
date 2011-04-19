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

#ifndef INTEGRATE_H
#define INTEGRATE_H

#include "../Blitz.h"
#include "../Splines/Grid.h"

void RungeKutta4 (const Grid &grid, 
		  int StartPoint, int EndPoint,
		  Array<scalar, 2> &Result,
		  Array<scalar,1> (*DerivFunc)(scalar x,
					       Array<scalar,1> y,
					       void *Params),
		  void *Params);



void IntegrateFirstOrder (const Grid &grid,
			   int StartPoint, int EndPoint,
			   Array<scalar,1> &Result,
			   scalar (*DerivFunc)(scalar x, scalar y,
					     void *Params),
			  void *Params);
void IntegrateFirstOrderNS (const Grid &grid,
			    int StartPoint, int EndPoint,
			    Array<scalar,1> &Result,
			    scalar (*DerivFunc)(scalar x, scalar y,
						void *Params),
			    void *Params);
void IntegrateSecondOrder (const Grid &grid,
			   int StartPoint, int EndPoint,
			   Array<Vec2,1> &Result,
			   Vec2 (*DerivFunc)(scalar x, Vec2 y,
					     void *Params),
			   void *Params);


#endif
