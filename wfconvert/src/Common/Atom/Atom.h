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

#ifndef ATOM_H
#define ATOM_H

#include "../PH/PH.h"
#include "../IO/IO.h"
#include <string>

inline double SmoothStep (double x, double CutOff, double Width)
{
  Width /= 10.0;
  return (1.0/(1 + exp(-(x-CutOff)/Width)));
}




class RadialWF
{
private:
  int TurningIndex();
  double IntegrateInwardOutward(int &Tindex);
  int CountNodes();
  double KineticEnergy();
  void OriginBC(double r0, double &u0, double &du0);
  void InfinityBC(double r0, double &u0, double &du0);

public:
  Grid *grid;
  PseudoHamiltonian *PH;

  bool IsRelativistic;
  double Energy, PartialNorm;
  int l;
  int CoreNodeNum;
  int DesiredNodeNum;
  double Weight;
  double CutOffRadius;
  double Occupancy;
  double MaxCurve;      // Stores the limit on the curvature allowed.
  CubicSpline u, dudr;
  string Label;
  //Array<double,1> R;
  //Array<double,1> dRdr;

  void Write (IOSectionClass &out);
  void Read  (IOSectionClass &in);

  inline double operator()(double r)
  {
    return (u(r));
  }
  inline double Deriv (double r)
  {
    return (dudr(r));
  }
  inline double Deriv2 (double r)
  {
    return (dudr.Deriv(r));
  }

  inline double WeightFunc(double r)
  {
    double rend = grid->End;
    double weight;
    double width1 = rend/5.0;
    double width2 = rend/10.0;
    weight = SmoothStep (r, CutOffRadius, width1);
    weight += 10.0 * SmoothStep (r, rend, width2);
    double norm = (rend-CutOffRadius) + 0.3634;
    weight *= Weight/norm;
    return (weight);
  }

  void IntegrateOutward();
  void SolveRadialEquation();
  void Normalize();

  inline void Init(int L, double E, PseudoHamiltonian *PseudoH, Grid *NewGrid) 
  {
    grid = NewGrid;
    l = L;
    Energy = E;
    PH = PseudoH;
    Array<double,1> temp(grid->NumPoints);
    temp = 0.0;
    u.Init(grid, temp, 5.0e30, 5.0e50);
    dudr.Init(grid,temp, 5.0e30, 5.0e30);
    IsRelativistic = false;
  }
  RadialWF()
  {
    IsRelativistic = false;
  }
};



class Atom
{
public:
  PseudoHamiltonian *PH;
  Grid *grid;
  int NumRadialWFs;
  bool DoSelfConsistent;
  Array<RadialWF,1> RadialWFs;
  CubicSpline ChargeDensity;
  CubicSpline Hartree;          // Holds the hartee potential
  CubicSpline ExCorr;           // Holds the exchange correlation potential
  double NewMix;                // Fraction of new charge density to mix
                                // in in self-conistent loop.
  double MaxChargeChange;       // What the maximum fractional charge
                                // density change can be before we're
                                // converged in the SC loop.


  void CalcHartreePot();
  void CalcExCorrPot();
  void UpdateVHXC();
  void CalcPartialNorms();
  void CalcChargeDensity();
  void CalcChargeModes();
  void SelfConsistentSolve();
  double TotalEnergy();
  double TotalCharge();
  void IntegrateRadialEquations();
  void NormalizeRadialWFs();
  void Read (IOSectionClass &IO);
  void ReadInput (IOSectionClass &in);
  void Write(IOSectionClass &out);
  void WritePH (char *FileName, Grid *OutputGrid);
  void WriteDavidsPH(char *FileName);
  void WriteWFs(char *FileName);
  // Write the States section of an input file
  void WriteStates(FILE *fout);
  PH_CubicSplineXC *CalcPH_XC(Grid *Agrid, Grid *Bgrid, Grid *Vgrid);

  void Init(PseudoHamiltonian *ph, Grid *newgrid, 
	    Array<int,1> l, Array<double,1> Energies, 
	    Array<double,1> Weights, Array<double,1> Occupancies, 
	    Array<int,1> CoreNodes, Array<int,1> Nodes, 
	    Array<double,1> CutOffRadii,  Array<double,1> MaxCurve,
	    Array<string,1> Labels, bool IsRel)
  {
    PH = ph;
    NumRadialWFs = l.rows();
    RadialWFs.resize(l.rows());
    grid = newgrid;
    Array<double,1> temp(grid->NumPoints);
    temp = 0.0;
    for (int i=0; i<l.rows(); i++)
      {
	RadialWFs(i).l = l(i);
	RadialWFs(i).Energy = Energies(i);
	RadialWFs(i).CoreNodeNum = CoreNodes(i);
	RadialWFs(i).DesiredNodeNum = Nodes(i);
	RadialWFs(i).Weight = Weights(i);
	RadialWFs(i).Occupancy = Occupancies(i);
	RadialWFs(i).PH = ph;
	RadialWFs(i).grid = grid;
	// Initialize cubic splines to 0
	RadialWFs(i).u.Init(grid, temp, 5.0e30, 5.0e30);
	RadialWFs(i).dudr.Init(grid, temp, 5.0e30, 5.0e30);
	RadialWFs(i).IsRelativistic = IsRel;
	RadialWFs(i).CutOffRadius = CutOffRadii(i);
	RadialWFs(i).MaxCurve = MaxCurve(i);
	RadialWFs(i).Label = Labels(i);
      }
  } 
    
  Atom(PseudoHamiltonian *PseudoH, Array<double,1> Energies, Grid *grid)
  {
    PH = PseudoH;
    NumRadialWFs = Energies.rows();
    RadialWFs.resize(NumRadialWFs);
    for (int l=0; l<NumRadialWFs; l++)
      RadialWFs(l).Init(l, Energies(l), PseudoH, grid);
  }
  Atom()
  {
    // Do nothing
  }
};



Array<double, 2> Inverse (Array<double,2> A);
Array<double,1> Prod(Array<double,2> A, Array<double,1>b);
Array<double,1> PseudoRadialDerivs(double r, Array<double,1> u_and_du,
				   void *WFptr);
Vec2 PseudoRadialDerivs(double r, Vec2 u_and_du,
			void *WFptr);


#endif
