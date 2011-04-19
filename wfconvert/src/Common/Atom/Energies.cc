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

#include "Energies.h"
#include "Integrate.h"
#include "Functionals.h"

scalar KEderivRel (scalar r, scalar KE, void *WFptr)
{
  RadialWF &WF = *(RadialWF *)WFptr;
  PseudoHamiltonian &PH = *WF.PH;

  const scalar alpha = 1.0/137.036;
  const scalar kappa = -1.0;

  scalar u = WF(r);
  scalar du = WF.Deriv(r);
  scalar d2u = WF.Deriv2(r);
  scalar sl = (scalar) WF.l;
  scalar a2 = alpha * alpha;
  scalar V = WF.PH->V(r);
  scalar dV = WF.PH->dVdr(r);
  scalar E = WF.Energy;
  scalar M = 1.0-0.5*a2*(V-E);

  scalar Tu = -1.0/(2.0*M)*(d2u  +(-sl*(sl+1.0)/(r*r) /*- 2.0*a2*(V-E)*/)*u)
    -0.25*a2/(M*M)*dV*(du+u*kappa/r)  ;

  //Tu = -1.0/(2.0*M) * (d2u /*+ 2.0/r*du*/ - sl*(sl+1.0)*u/(r*r))
  //  -dV*du*a2/(4.0*M*M) - (kappa/*+1.0*/)/r * dV*u*a2/(4.0*M*M);
  scalar uTu = u*Tu;
  return (4.0 * M_PI * uTu);

}


scalar KEderiv(scalar r, scalar E, void *WFptr)
{
  RadialWF &WF = *(RadialWF *)WFptr;
  PseudoHamiltonian &PH = *WF.PH;

  scalar A, B, V, dAdr;
  scalar l = WF.l;
  PH.ABV(r, A, B, V, dAdr);

  scalar u = WF(r);
  scalar up = WF.Deriv(r);
  scalar eps = 1e-8;

  scalar upp = WF.Deriv2(r);
  
  // HACK
  if ((r< 49.0) && (r>2.1e-8))
    upp = (WF.Deriv(r+eps) - WF.Deriv(r-eps)) / (2.0*eps);
  else
    upp = 0.0;
  
  scalar Tu = -0.5 * ((up - u/r) * dAdr + upp * A -
		      l*(l+1.0) * B * u/ (r*r));

  scalar uTu = u*Tu;
  return (4.0 * M_PI * uTu);
}


scalar
RadialWF::KineticEnergy()
{
  int N = grid->NumPoints;
  Array<scalar,1> Temp(N);
  Temp(0) = 0.0;
  scalar (*DerivFunc)(scalar r, scalar E, void *WFptr);

  if (IsRelativistic)
    DerivFunc = KEderivRel;
  else
    DerivFunc = KEderiv;


  IntegrateFirstOrder (*grid, 0, N-1, Temp, DerivFunc, this);

  cerr << "KE = " << Temp(N-1) << "\n";

  return (Temp(N-1));
}




scalar PEderiv(scalar r, scalar E, void *Atomptr)
{
  Atom &atom = *(Atom *)Atomptr;

  scalar V = //grid.Interp(atom.Hartree, r, 4) +
    //grid.Interp(atom.ExCorr, r, 4) +
    atom.PH->V(r) - 0.5 * atom.Hartree(r);
  
  scalar rho = atom.ChargeDensity(r);

  scalar nup = 0.5 * rho;
  scalar ndown = nup;
  scalar Vex = atom.ExCorr(r);
  scalar Eex = FortranXCE (nup, ndown);
  
  V += Eex - Vex;

  
  return (4.0 * M_PI * r * r * rho * V);
  // HACK
  //return (4.0 * M_PI * r * r * rho * atom.PH->V(r));
}


scalar HEderiv(scalar r, scalar E, void *Atomptr)
{
  Atom &atom = *(Atom *)Atomptr;

  scalar V = atom.Hartree(r);
  scalar rho = atom.ChargeDensity (r);
  return (2.0 * M_PI * r * r * rho * V);
}

scalar XCEderiv(scalar r, scalar E, void *Atomptr)
{
  Atom &atom = *(Atom *)Atomptr;

  scalar rho = atom.ChargeDensity (r);

  scalar nup = 0.5 * rho;
  scalar ndown = nup;
  scalar Vex = atom.ExCorr(r);
  scalar Eex = FortranXCE (nup, ndown);

  return (4.0 * M_PI * r * r * rho * Eex);
}


scalar NUCderiv(scalar r, scalar E, void *Atomptr)
{
  Atom &atom = *(Atom *)Atomptr;
  Grid &grid = *atom.RadialWFs(0).grid;

  scalar V = -atom.PH->FullCoreV->Z / r;
  scalar rho = atom.ChargeDensity(r);
  return (4.0 * M_PI * r * r * rho * V);
}



scalar 
Atom::TotalEnergy()
{
  scalar TotalE;
  for (int i=0; i<NumRadialWFs; i++)
    TotalE += RadialWFs(i).KineticEnergy() * RadialWFs(i).Occupancy;

  fprintf (stderr, "Kinetic Energy = %1.6f\n", TotalE);

  CalcChargeDensity();
  CalcHartreePot();
  CalcExCorrPot();


  Grid &grid = *RadialWFs(0).grid;
  int N = grid.NumPoints;

  Array<scalar,1> Temp(N);
  Temp(0) = 0.0;
  IntegrateFirstOrderNS(grid, 0, N-1, Temp, PEderiv, this);
  scalar PE = Temp(N-1);
  fprintf (stderr, "Potential Energy = %1.6f\n", PE);
  
  TotalE += Temp(N-1);
  
  Temp(0) = 0.0;
  IntegrateFirstOrderNS(grid, 0, N-1, Temp, HEderiv, this);
  scalar HartreeE = Temp(N-1);

  Temp(0) = 0.0;
  IntegrateFirstOrderNS(grid, 0, N-1, Temp, XCEderiv, this);
  scalar XCE = Temp(N-1);

  //  Temp(0) = 0.0;
  // IntegrateFirstOrderNS(grid, 0, N-1, Temp, NUCderiv, this);
  //scalar NUC = Temp(N-1);

  scalar NucE = PE - HartreeE - XCE;
  fprintf (stderr, "HartreeE Energy = %1.6f\n", HartreeE);
  fprintf (stderr, "ExCorrE Energy = %1.6f\n", XCE);
  fprintf (stderr, "Nuclear Energy = %1.6f\n", NucE);
  //cerr << "NUC     Energy = " << NUC << "\n";



  return (TotalE);
}
