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

#ifndef SPECFUNC_H
#define SPECFUNC_H
#include "../Blitz.h"

double ModBesselScaled(int l, double x);
double RegModBesselFracScaled(double nu, double x);
double RegModBesselFrac(double nu, double x);
double Legendre (int l, double CosTheta);
double LogRegModBesselScaled(int l, double x);

double il(int l, double z);
double dil_dz (int l, double z);
double d2il_dz2 (int l, double z);

double il_scaled(int l, double z);
double log_il_scaled(int l, double z);
double dil_dz_scaled (int l, double z);
double d2il_dz2_scaled (int l, double z);

void Test_il_Derivs();
void Test_il_series_large();
void Test_il_series_small();
double jl(int l, double z);
#endif
