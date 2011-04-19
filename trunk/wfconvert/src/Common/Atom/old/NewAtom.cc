#include "NewAtom.h"

void DFTAtom::CalcVHXC()
{
  UpdateChargeDensity();
  UpdateHartree();
  UpdateExCorr();

  for (int i=0; i < grid->NumPoints; i++)
    xcpot.XC(i) = Hartree(i) + ExCorr(i);
}

/// Radial WFs must be normalized before calling this function
void DFTAtom::UpdateChargeDensity()
{
  for (int i=0; i<grid->NumPoints; i++) {
    temp(i) = 0.0;
    double r = (*grid)(i);
    double rinv2 = 1.0/(r*r);
    for (int j=0; j<RadialWFs.size(); j++) {
      double u = RadialWFs(i).u(j);
      temp(i) += RadialWFs(i).Occupancy * u*u * rinv2;
    }
  }
}


class HartreeDeriv1
{
  DFTAtom &atom;
public:
  inline double operator()(double r, double sum)
  { return atom.Hartree1(r,sum); }

  HartreeDeriv1(DFTAtom &newAtom) : atom(newAtom) 
  { /* Do nothing */ }
};

class HartreeDeriv2
{
  DFTAtom &atom;
public:
  inline double operator()(double r, double sum)
  { return atom.Hartree2(r,sum); }

  HartreeDeriv2(DFTAtom &newAtom) : atom(newAtom) 
  { /* Do nothing */ }
};


void DFTAtom::UpdateHartree()
{


}

void DFTAtom::SetGrid(Grid *newgrid)
{
  grid = newgrid;
  int N = grid->NumPoints;
  temp.resize(N);
  temp = 0.0;
  xcpot.XC.Init(grid, temp);
  ChargeDensity.Init(grid,temp);
  Hartree.Init(grid, temp);
  ExCorr.Init(grid,temp);
  for (int i=0; i<RadialWFs.size(); i++)
    RadialWFs(i).SetGrid(grid);
}

void DFTAtom::SetBarePot(Potential *newPot)
{
  BarePot = newPot;
  xcpot.BarePot = BarePot;
}
