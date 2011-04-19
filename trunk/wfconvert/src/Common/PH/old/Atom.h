#ifndef ATOM_H
#define ATOM_H
#include "PH.h"
#include <string>

inline scalar SmoothStep (scalar x, scalar CutOff, scalar Width)
{
  Width /= 10.0;
  return (1.0/(1 + exp(-(x-CutOff)/Width)));
}




class RadialWF
{
public:
  bool IsRelativistic;
  Grid *grid;
  PseudoHamiltonian *PH;
  scalar Energy, PartialNorm;
  int l;
  int CoreNodeNum;
  int DesiredNodeNum;
  scalar Weight;
  scalar CutOffRadius;
  scalar Occupancy;
  scalar MaxCurve;      // Stores the limit on the curvature allowed.
  CubicSpline u, dudr;
  string Label;
  //Array<scalar,1> R;
  //Array<scalar,1> dRdr;


  inline scalar operator()(scalar r)
  {
    return (u(r));
  }
  inline scalar Deriv (scalar r)
  {
    return (dudr(r));
  }
  inline scalar Deriv2 (scalar r)
  {
    return (dudr.Deriv(r));
  }

  inline scalar WeightFunc(scalar r)
  {
    scalar rend = grid->End;
    scalar weight;
    scalar width1 = rend/5.0;
    scalar width2 = rend/10.0;
    weight = SmoothStep (r, CutOffRadius, width1);
    weight += 10.0 * SmoothStep (r, rend, width2);
    scalar norm = (rend-CutOffRadius) + 0.3634;
    weight *= Weight/norm;
    return (weight);
  }

  int CountNodes();
  int TurningIndex();
  void IntegrateOutward();
  scalar IntegrateInwardOutward(int &Tindex);
  void SolveRadialEquation();
  scalar KineticEnergy();
  void Normalize();
  void OriginBC(scalar r0, scalar &u0, scalar &du0);
  void InfinityBC(scalar r0, scalar &u0, scalar &du0);

  inline void Init(int L, scalar E, PseudoHamiltonian *PseudoH, 
		   Grid *NewGrid) 
  {
    grid = NewGrid;
    l = L;
    Energy = E;
    PH = PseudoH;
    Array<scalar,1> temp(grid->NumPoints);
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
  scalar NewMix;                // Fraction of new charge density to mix
                                // in in self-conistent loop.
  scalar MaxChargeChange;       // What the maximum fractional charge
                                // density change can be before we're
                                // converged in the SC loop.


  void CalcHartreePot();
  void CalcExCorrPot();
  void UpdateVHXC();
  void CalcPartialNorms();
  void CalcChargeDensity();
  void CalcChargeModes();
  void SelfConsistentSolve();
  scalar TotalEnergy();
  scalar TotalCharge();
  void IntegrateRadialEquations();
  void NormalizeRadialWFs();
  void WritePH (char *FileName, Grid *OutputGrid);
  void WriteDavidsPH(char *FileName);
  void WriteWFs(char *FileName);
  // Write the States section of an input file
  void WriteStates(FILE *fout);
  PH_CubicSplineXC *CalcPH_XC(Grid *Agrid, Grid *Bgrid, Grid *Vgrid);

  void Init(PseudoHamiltonian *ph, Grid *newgrid, 
	    Array<int,1> l, Array<scalar,1> Energies, 
	    Array<scalar,1> Weights, Array<scalar,1> Occupancies, 
	    Array<int,1> CoreNodes, Array<int,1> Nodes, 
	    Array<scalar,1> CutOffRadii,  Array<scalar,1> MaxCurve,
	    Array<string,1> Labels, bool IsRel)
  {
    PH = ph;
    NumRadialWFs = l.rows();
    RadialWFs.resize(l.rows());
    grid = newgrid;
    Array<scalar,1> temp(grid->NumPoints);
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
    
  Atom(PseudoHamiltonian *PseudoH, Array<scalar,1> Energies, Grid *grid)
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



Array<scalar, 2> Inverse (Array<scalar,2> A);
Array<scalar,1> Prod(Array<scalar,2> A, Array<scalar,1>b);
Array<scalar,1> PseudoRadialDerivs(scalar r, Array<scalar,1> u_and_du,
				   void *WFptr);
Vec2 PseudoRadialDerivs(scalar r, Vec2 u_and_du,
			void *WFptr);


#endif
