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

#include "Atom.h"
#include "../DFT/Functionals.h"

//////////////////////////////////////////////////////////////////
//                    Generic Atom Routines                     // 
//////////////////////////////////////////////////////////////////
double NormalizationDeriv (double r, double sum,
			   void *WFptr)
{
  double deriv;
  RadialWF &WF = *((RadialWF *)WFptr);
  double u = WF(r);
  deriv = u*u;
  return (deriv);
}


void
RadialWF::Normalize()
{
  Array<double,1> Temp(grid->NumPoints);
  Temp(0) = 0.0;
  IntegrateFirstOrderNS(*grid, 0, grid->NumPoints-1, Temp,
			NormalizationDeriv, this);

  double norm = 
	sqrt(1.0/(4.0*M_PI*Temp(grid->NumPoints-1)));

  for (int i=0; i<grid->NumPoints; i++)
    {
      u.Params(i) *= norm;
      dudr.Params(i) *= norm;
    }
}



void
Atom::NormalizeRadialWFs()
{
  for (int l=0; l<NumRadialWFs; l++)
    RadialWFs(l).Normalize();
}



double NormDeriv (double r, double sum,
		  void *WFptr)
{
  RadialWF &WF = *((RadialWF *)WFptr);
  double u = WF(r);
  return(u*u);
}



void
Atom::CalcPartialNorms()
{
  Array<double,1> Temp(grid->NumPoints);
  for (int i=0; i<NumRadialWFs; i++)
    {
      Temp(0) = 0.0;
      IntegrateFirstOrder(*grid, 0, grid->NumPoints-1, Temp, 
			  NormDeriv, &(RadialWFs(i)));
      double uend = RadialWFs(i).u.Params(grid->NumPoints-1);
      RadialWFs(i).PartialNorm = Temp(grid->NumPoints-1) / (uend*uend);
    }
}



Vec2 ModeDerivs (double r, Vec2 Sum,
		 void *Atomptr)
{

  Atom &atom = *(Atom *)Atomptr;
  double Charge = 0.0;
  for (int i=0; i<atom.NumRadialWFs; i++)
    {
      double u = atom.RadialWFs(i)(r);
      double uend = atom.RadialWFs(i).u.Params(atom.grid->NumPoints-1);
      Charge += atom.RadialWFs(i).Occupancy*u*u/(uend*uend);
    }
  Vec2 derivs;
  derivs[0] = r*Charge;
  derivs[1] = r*derivs[0];
  return (derivs);
}


//////////////////////////////////////////////////////////////////
//                 Input and output routines                    //
//////////////////////////////////////////////////////////////////
void Atom::ReadInput (IOSectionClass &IO)
{
  assert (IO.OpenSection("Potential"));
  PH = ReadPH(IO);
  IO.CloseSection();
  assert (IO.OpenSection("RadialWFs"));
  assert (IO.OpenSection("Grid"));
  grid = ReadGrid(IO);
  IO.CloseSection();
  bool IsPseudo;

  assert(IO.ReadVar("NumRadialWFs", NumRadialWFs));
  assert(IO.ReadVar("IsPseudo", IsPseudo));

  Array <double,1> Energies, Weights, Occupancies, MaxCurves;
  Array <int,1> ls, CoreNodeNums, TotalNodeNums;
  Array <string,1> Labels;
  
  bool HaveEnergies      = IO.ReadVar("Energies", Energies);
  bool HaveWeights       = IO.ReadVar("Weights", Weights);
  bool HaveOccupancies   = IO.ReadVar("Occupancies", Occupancies);
  bool HaveMaxCurves     = IO.ReadVar("MaxCurves", MaxCurves);
  bool Havels            = IO.ReadVar("ls", ls);
  bool HaveCoreNodeNums  = IO.ReadVar("CoreNodeNums", CoreNodeNums);
  bool HaveTotalNodeNums = IO.ReadVar("TotalNodeNums", TotalNodeNums);
  bool HaveLabels        = IO.ReadVar("Labels", Labels);
  IO.CloseSection();
  
  RadialWFs.resize(NumRadialWFs);
  Array<double,1> temp(grid->NumPoints);
  temp = 0.0;
  for (int i=0; i<NumRadialWFs; i++) {
    RadialWFs(i).IsRelativistic = !IsPseudo;
    RadialWFs(i).grid = grid;
    RadialWFs(i).PH = PH;
    RadialWFs(i).u.Init(grid, temp, 5.0e50, 5.0e50);
    RadialWFs(i).dudr.Init(grid, temp, 5.0e50, 5.0e50);
    if (HaveEnergies) RadialWFs(i).Energy = Energies(i);
    if (HaveWeights)  RadialWFs(i).Weight = Weights(i);
    if (HaveOccupancies) RadialWFs(i).Occupancy = Occupancies(i);
    if (HaveMaxCurves)   RadialWFs(i).MaxCurve  = MaxCurves(i);
    if (Havels)       RadialWFs(i).l      = ls(i);
    else { cerr << "Cannot find variable ls\n"; exit (1); }
    if (HaveCoreNodeNums) RadialWFs(i).CoreNodeNum = CoreNodeNums(i);
    if (HaveTotalNodeNums) RadialWFs(i).DesiredNodeNum = TotalNodeNums(i);
    if (HaveLabels) RadialWFs(i).Label = Labels(i);
  }    
}


//////////////////////////////////////////////////////////////////
//              PseudoAtom Integration Routines                 //
//////////////////////////////////////////////////////////////////



Vec2 PseudoRadialDerivs(double r, Vec2 u_and_du,
			void *WFptr)
{
  Vec2 deriv;
  RadialWF &WF = *((RadialWF *)WFptr);
  PseudoHamiltonian &PH = *(WF.PH);

  // d/dr u = u'
  deriv[0] = u_and_du[1];

  double A, B, V, dAdr;

  PH.ABV(r, A, B, V, dAdr);
  
  double u = u_and_du[0];
  double dudr = u_and_du[1];
  double l = (double)WF.l;
  double E = WF.Energy;

  deriv[1] = 1.0/A * 
    (-dAdr*dudr + (dAdr/r + l*(l+1.0)*B/(r*r) + 2.0*(V-E))*u);

  return (deriv);
}



Vec2 RelativisticRadialDerivs(double r, Vec2 u_and_du,
			      void *WFptr)
{
  RadialWF &WF = *((RadialWF *)WFptr);
  PseudoHamiltonian &PH = *(WF.PH);

  const double alpha = 1.0/137.036;
  const double kappa = -1.0;
  Vec2 derivs;
  double u, dudr;
  
  double V, dVdr;

  //V = (*PH.FullCoreV)(r);
  //dVdr = PH.FullCoreV->Deriv(r);
  V = PH.V(r);
  dVdr = PH.dVdr(r);

  u = u_and_du[0];
  dudr = u_and_du[1];

  derivs[0] = dudr;
  double E = WF.Energy;
  double M = 1.0 - alpha*alpha*0.5*(V-E);
  double l = (double) WF.l;

  derivs[1] = (l*(l+1.0)/(r*r) + 2.0 * M *(V-E)) * u
    - 0.5*alpha*alpha/M*dVdr*(dudr + u*kappa/r);
  
  return (derivs);
}



// Setup boundary conditions at the origin
void 
RadialWF::OriginBC(double r0, double &u0, double &du0)
{
  if (IsRelativistic)
    {
      const double alpha = 1.0/137.036;
      double Z = PH->Z;
      double a2Z2 = alpha*alpha*(double)Z*(double)Z;

      double sl = (double) l;
      double gamma = 
	(sl+1.0)*sqrt((sl+1.0)*(sl+1.0)-a2Z2) /
	(2.0*sl + 1.0); 
      if (l!=0)
	gamma += sl*sqrt(sl*sl-a2Z2)/(2.0*sl+1.0);
      
      u0 = pow(r0, gamma);
      du0 = gamma * pow(r0,gamma-1.0);
    }
  else
    {
      u0 = 1.0;
      du0 = (double)(l+1)/r0;
      if (l == 0)
	{
	  double E = Energy;
	  double A0, B0, V0, dAdr0;
	  double A, B, V, dAdr;
	  PH->ABV(r0, A, B, V, dAdr);
	  PH->ABV(0.0, A0, B0, V0, dAdr0);
	  if (V0 > Energy)
	    {
	      double kappa = sqrt(2.0 * (V-E)/A);
	      u0 = sinh(kappa * r0);
	      du0 = kappa * cosh (kappa * r0);
	    }
	  else
	    {
	      double k = sqrt(2.0 * (E-V)/A);
	      u0 = sin(k * r0);
	      du0 = k * cos(k * r0);
	    }
	}      
    }  

  if ((DesiredNodeNum % 2) == 1)
    {
      u0 *= -1.0;
      du0 *= -1.0;
    }

}



void 
RadialWF::InfinityBC(double rend, double &uend, double &duend)
{
  double Aend,Bend,Vend,dAdrend;
  PH->ABV(rend, Aend,Bend,Vend,dAdrend);
  double dVdrend = PH->dVdr(rend);
  double E = Energy;
  double k = sqrt(l*(l+1.0)/(rend*rend) + (Vend-E));
  uend = 1.0;
  duend = -uend * (k + 0.5*rend/k*(dVdrend - 2.0*l*(l+1.0)/(rend*rend*rend)));
}




void
RadialWF::IntegrateOutward()
{
  Vec2 (*DerivFunc)(double r, Vec2 u_and_du, void *WFptr);

  if (IsRelativistic)
    DerivFunc = RelativisticRadialDerivs;
  else
    DerivFunc = PseudoRadialDerivs;

  Array<Vec2,1> Temp(grid->NumPoints);

  // Set up initial conditions:
  double r0 = (*grid)(0); 
  OriginBC(r0, Temp(0)[0], Temp(0)[1]);

  // Now do integration

  IntegrateSecondOrder(*grid, 0, grid->NumPoints-1, Temp, DerivFunc, this);

  // Copy results of integration into RadialWF
  for (int i=0; i < grid->NumPoints; i++)
    {
      u.Params(i) = Temp(i)[0];
      dudr.Params(i) = Temp(i)[1];
    }
  Normalize();
}



void
Atom::IntegrateRadialEquations()
{
  for (int i=0; i<NumRadialWFs; i++)
    {
      RadialWFs(i).IntegrateOutward();
      double uend = RadialWFs(i).u.Params(grid->NumPoints-1);
      if (uend == 0.0)
	cerr << "uend = 0!!!!!!!!!!!!!!!!!!!!!1\n";
      for (int j=0; j<RadialWFs(i).grid->NumPoints; j++)
	{
	  RadialWFs(i).u.Params(j) /= uend;
	  RadialWFs(i).dudr.Params(j) /= uend;
	}
    }
  CalcPartialNorms();  
}





//////////////////////////////////////////////////////////////////
//                 EigenState related routines                  //
//////////////////////////////////////////////////////////////////

int RadialWF::CountNodes()
{
  int nodes = 0;
  double sign = u.Params(0);
  for (int i=1; i<grid->NumPoints; i++)
    if(sign*u.Params(i) < 0.0)
      {
	nodes++;
	sign = u.Params(i);
      }
  return(nodes);
}


int RadialWF::TurningIndex()
{
  // Start in classically forbidden region and search inward toward
  // origin looking for classically allowed region

  int index = grid->NumPoints-1;
  int done = 0;
  while (!done && (index >=0))
    {
      double r = (*grid)(index);
      
      double A,B,V,dAdr;
      PH->ABV(r,A,B,V,dAdr);
      double E = Energy;
      double sl = (double) l;
      
      double Deriv2 = dAdr/r + sl*(sl+1.0)*B/(r*r) + 2.0*(V-E);
      
      if (Deriv2 < 0.0)
	done = 1;
      else
	index--;
    }
  if (index == -1)
    index = grid->ReverseMap(0.5*(*grid)(grid->NumPoints-1));

  return (index);
}




double 
RadialWF::IntegrateInwardOutward(int &TIndex)
{
  Vec2 (*DerivFunc)(double r, Vec2 u_and_du,
		      void *WFptr);
  if (IsRelativistic)
    DerivFunc = RelativisticRadialDerivs;
  else
    DerivFunc = PseudoRadialDerivs;


  // Allocate a temporary array for integration
  Array<Vec2,1> Temp(grid->NumPoints);

  // Find classical turning point
  TIndex = TurningIndex();

  // Compute starting value and derivative at origin
  double r0 = (*grid)(0);
  OriginBC (r0, Temp(0)[0], Temp(0)[1]);
 
  // Do integration from the origin to the turning point
  IntegrateSecondOrder(*grid, 0, TIndex, Temp, DerivFunc, this);
  for (int j=0; j<=TIndex; j++)
    {
      u.Params(j) = Temp(j)[0];
      dudr.Params(j) = Temp(j)[1];
    }

  
  // Now initilize value and derivative at rmax
  int EndPoint = grid->NumPoints-1;
  double rend = (*grid)(EndPoint);

  InfinityBC(rend, Temp(EndPoint)[0], Temp(EndPoint)[1]);
  
  // Do integration from the right to the turning point
  IntegrateSecondOrder(*grid,EndPoint, TIndex, Temp, DerivFunc, this);  
  
  if ((u.Params(TIndex)*Temp(TIndex)[0]) < 0.0)
    for (int i=TIndex; i<=EndPoint; i++)
      {
	Temp(i)[0] *= -1.0;
	Temp(i)[1] *= -1.0;
      }

  double CuspValue = Temp(TIndex)[1]/Temp(TIndex)[0] -
    dudr.Params(TIndex)/u.Params(TIndex);

  // Copy values in, normalizing to make R continuous
  double factor = u.Params(TIndex) / Temp(TIndex)[0];
  for (int i=TIndex+1; i<grid->NumPoints; i++)
    {
      u.Params(i)    = factor*Temp(i)[0];
      dudr.Params(i) = factor*Temp(i)[1];
    }
  return (CuspValue);
}





/*void 
RadialWF::SolveRadialEquation()
{
  int Tindex;
  const double Tolerance = 1e-12;
  double Ehigh, Elow, Etrial;
  Ehigh = 0.0;
  
  double r = (*grid)(0);
  double A,B,V,dAdr;
  PH->ABV(r,A,B,V,dAdr);

  if (IsRelativistic)
    {
      double n = DesiredNodeNum + 1.0;
      Elow = -1.5*PH->Z*PH->Z/(n*n);
    }
  else
    {
      // Eigenenergy can't be lower than lowest potential -- otherwise
      // we'd have no turning point.
      Elow = PH->V((*grid)(0));
      for (int i=1; i<grid->NumPoints; i++)
	{
	  double r = (*grid)(i);
	  if (PH->V(r) < Elow)
	    Elow = PH->V(r);
	}
    }

  // HACK
  //Elow = -300.0;
  
  while ((Ehigh-Elow) > Tolerance)
    {
      double Etrial = 0.5 * (Ehigh + Elow);
      Energy = Etrial;
      double CuspValue = IntegrateInwardOutward(Tindex);
      int NumNodes = CountNodes();
      if (NumNodes > DesiredNodeNum)
	Ehigh = Etrial;
      else if (NumNodes < DesiredNodeNum)
	Elow = Etrial;
      else if  (CuspValue < 0.0)
	Elow = Etrial;
      else if (CuspValue > 0.0)
	Ehigh = Etrial;
      else if (isnan(CuspValue))
	Elow = Etrial;
    }
  Energy = 0.5*(Ehigh + Elow);
}
*/

void 
RadialWF::SolveRadialEquation()
{
  const double Tolerance = 1e-12;
  int Tindex;

  IntegrateInwardOutward(Tindex);
  int NumNodes = CountNodes();


  double Ehigh, Elow, Etrial;
  Ehigh = 0.0;
  
  if (IsRelativistic)
    {
      double n = DesiredNodeNum + 1.0;
      double Z = PH->FullCoreV->Z;
      Elow = -1.5*Z*Z/(n*n);
    }
  else
    {
      // Eigenenergy can't be lower than lowest potential -- otherwise
      // we'd have no turning point.
      Elow = PH->V((*grid)(0));
      for (int i=1; i<grid->NumPoints; i++)
	{
	  double r = (*grid)(i);
	  if (PH->V(r) < Elow)
	    Elow = PH->V(r);
	}
    }

  /*if (NumNodes != DesiredNodeNum)
    {  
      Etrial = 0.5 * (Ehigh + Elow);
      Energy = Etrial;
      IntegrateInwardOutward(Tindex);
      int NumNodes = CountNodes();
      while (NumNodes != DesiredNodeNum)
	{
	  
	  if (NumNodes > DesiredNodeNum)
	    Ehigh = Etrial;
	  else if (NumNodes < DesiredNodeNum)
	    Elow = Etrial;
	  Etrial = 0.5 *(Ehigh + Elow);
	  Energy = Etrial;
	  IntegrateInwardOutward(Tindex);
	  NumNodes = CountNodes();
	}
	}*/

  double Eold;
  Eold = Etrial = Energy;
  int done = 0;
  double A, B, V, dAdr;
  while (!done)
    {
      //cerr << "Etrial = " << Etrial << endl;
      //cerr << "Ehigh = " << Ehigh << " Elow = " << Elow << endl;
      double CuspValue = IntegrateInwardOutward(Tindex);
      //cerr << "Cusp value = " << CuspValue << "\n";
      
      NumNodes = CountNodes();
      //cerr << "NumNodes = " << NumNodes << "\n";

      if (NumNodes > DesiredNodeNum)
	Ehigh = Etrial;
      else if (NumNodes < DesiredNodeNum)
	Elow = Etrial;
      else if (CuspValue < 0.0)
	Elow = Etrial;
      else if (CuspValue > 0.0)
	Ehigh = Etrial;
      else if (isnan(CuspValue))
	Elow = Etrial;
     
      Normalize();
      PH->ABV((*grid)(Tindex), A, B, V, dAdr);
      double C = 0.5 * A * CuspValue;
      double u0 = u.Params(Tindex);
      double Etest = Eold - 4.0*M_PI*C*u0*u0;
      if ((Etest > Elow) && (Etest < Ehigh) && (NumNodes == DesiredNodeNum))
	Etrial = Etest;
      else 
	Etrial = 0.5 * (Ehigh + Elow);
			  
      Energy = Etrial;

      if ((NumNodes == DesiredNodeNum) && (fabs(Etrial - Eold) < Tolerance))
	 done = 1;
       Eold = Etrial;       
    }
  IntegrateInwardOutward(Tindex);
  Normalize();
  NumNodes = CountNodes();
  if (NumNodes != DesiredNodeNum)
    cerr << "Node number error!!!!!!!!!!.\n";

}



//////////////////////////////////////////////////////////////////
//            Exchange/Correlation related routines             //
//////////////////////////////////////////////////////////////////


double ChargeDensityDeriv(double r, double sum, void *atomptr)
{
  Atom &atom = *(Atom *)atomptr;
  return (4.0*M_PI * r * r *atom.ChargeDensity(r));
}

double
Atom::TotalCharge()
{
  Array<double,1> temp(grid->NumPoints);
  temp(0) = 0.0;
  IntegrateFirstOrderNS(*grid, 0, grid->NumPoints-1, temp,
			ChargeDensityDeriv, this);
  return (temp(grid->NumPoints-1));
}


void
Atom::CalcChargeDensity()
{
  Array<double,1> Rho(grid->NumPoints);

  Rho = 0.0;

  NormalizeRadialWFs();
  for (int i=0; i<grid->NumPoints; i++)
    {
      Rho(i) = 0.0;
      double r = (*grid)(i);
      double rinv2 = 1.0/(r*r);
      for (int n=0; n<NumRadialWFs; n++)
	{
	  double u = RadialWFs(n)(r);
	  Rho(i) += RadialWFs(n).Occupancy * u * u * rinv2;
	}
    }
  ChargeDensity.Init(grid, Rho, 5.0e30, 5.0e30);
  //  FILE *fout;
  //fout = fopen ("charge.dat", "w");
  //for (int i=0; i<grid->NumPoints; i++)
  // fprintf (fout, "%1.16e %1.16e\n", (*grid)(i), Rho(i));
  //fclose (fout);
}



double Hartree1 (double r, double N,
		 void *AtomPtr)
{
  Atom &atom = *((Atom *)AtomPtr);

  double rho= atom.ChargeDensity(r);
  return(rho*r*r); // dV/dr = V'
}

double Hartree2 (double r, double N,
		 void *AtomPtr)
{
  Atom &atom = *((Atom *)AtomPtr);
  double rho= atom.ChargeDensity(r);
  return(rho*r); // dV/dr = V'
}




void
Atom::CalcHartreePot ()
{
  int NumPoints = grid->NumPoints;
  int EndPoint = NumPoints-1;
  Array<double,1> term1(NumPoints), term2(NumPoints);
  Array<double,1> HV(NumPoints);

  term1(0) = 0.0;
  IntegrateFirstOrder(*grid, 0, EndPoint, term1, Hartree1, this);
  term2(0) = 0.0;
  IntegrateFirstOrder(*grid, 0, EndPoint, term2, Hartree2, this);
  
  double term2infty = term2(EndPoint);
  for (int i=0; i<=EndPoint; i++)
    {
      double r = (*grid)(i);
      HV(i) = 4.0*M_PI*(term1(i)/r + (term2infty-term2(i)));
    }
  Hartree.Init(grid, HV, 5.0e30, 5.0e30);
}

void
Atom::CalcExCorrPot ()
{
  int NumPoints = grid->NumPoints;
  Array<double,1> VXC(NumPoints);

  for (int i=0; i<NumPoints; i++)
    {
      double Vxcup, Vxcdown;
      double UpCharge = 0.5 * ChargeDensity.Params(i);
      double DownCharge = UpCharge;
      FortranExCorr(UpCharge,DownCharge, Vxcup, Vxcdown);
      VXC(i) = Vxcup;
    }
  ExCorr.Init(grid,VXC, 5.0e30, 5.0e30);
}


void
Atom::UpdateVHXC()
{
  CalcHartreePot();
  CalcExCorrPot();

  Array<double,1> VHXC(grid->NumPoints);
  for (int i=0; i<grid->NumPoints; i++)
    VHXC(i) = Hartree.Params(i) + ExCorr.Params(i);

  double NumElecs = 0.0;
  for (int i=0; i<NumRadialWFs; i++)
    NumElecs += RadialWFs(i).Occupancy;
  PH->VHXC.Init(grid, VHXC, 5.0e30, 5.0e30);
  PH->UseVHXC = 1;
  PH->NumElecs = NumElecs;
  //PH->UpdateVHXC(VHXCSpline, NumElecs);

}



void 
Atom::WriteWFs(char *FileName)
{
  FILE *fout;
  if ((fout = fopen (FileName, "w")) == NULL)
    {
      cerr << "Can't open " << FileName << " for writing.\n";
      exit(1);
    }
  for (int i=0; i < grid->NumPoints; i++)
    {
      fprintf (fout, "%1.16e ", (*grid)(i));
      for (int j=0; j < NumRadialWFs; j++)
	fprintf (fout, "%1.16e ", RadialWFs(j).u.Params(i));
      fprintf (fout, "\n");
    }
  fclose(fout);
}



void
Atom::WriteStates(FILE *fout)
{
  fprintf (fout, "States\n{\n");

  fprintf (fout, "  NumRadialWFs = %d;\n", NumRadialWFs);
  fprintf (fout, "  Energies = [ %1.16e\n", RadialWFs(0).Energy);
  for (int i=1; i<(NumRadialWFs-1); i++)
    fprintf (fout, "               %1.16e\n", RadialWFs(i).Energy);
  fprintf (fout, "               %1.16e ];\n", 
	   RadialWFs(NumRadialWFs-1).Energy);
  fprintf (fout, "  l = [ ");
  for (int i=0; i<NumRadialWFs; i++)
    fprintf (fout, "%d ", RadialWFs(i).l);
  fprintf (fout, "];\n");
  fprintf (fout, "  CoreNodes = [ ");
  for (int i=0; i<NumRadialWFs; i++)
    fprintf (fout, "%d ", RadialWFs(i).CoreNodeNum);
  fprintf (fout, " ];\n");
  fprintf (fout, "  PHNodes = [ ");
  for (int i=0; i<NumRadialWFs; i++)
    fprintf (fout, "%d ", RadialWFs(i).DesiredNodeNum);
  fprintf (fout, " ];\n");
  fprintf (fout, "  Weights = [ ");
  for (int i=0; i<NumRadialWFs; i++)
    fprintf (fout, "%1.5e ", RadialWFs(i).Weight);
  fprintf (fout, " ];\n");
  fprintf (fout, "  Occupancies = [ ");
  for (int i=0; i<NumRadialWFs; i++)
    fprintf(fout, "%1.3f ", RadialWFs(i).Occupancy);
  fprintf (fout, " ];\n");
  fprintf (fout, "  CutOffRadii = [ ");
  for (int i=0; i<NumRadialWFs; i++)
    fprintf (fout, "%1.3f ", RadialWFs(i).CutOffRadius);
  fprintf (fout, " ];\n");
  fprintf (fout, "}\n");
}


void
Atom::WriteDavidsPH(char *FileName)
{
  Grid &AtomGrid = *RadialWFs(0).grid;
  LogGrid DavidsGrid(5000, 26.0, 1.0e-3 * 26.0, 1.0024);

  //double IonCharge=0.0;
  //for (int i=0; i<NumRadialWFs; i++)
  // IonCharge += RadialWFs(i).Occupancy;

  CalcChargeDensity();
  CalcHartreePot();
  CalcExCorrPot();

  // Open the file for writing
  FILE *fout;
  if ((fout = fopen (FileName, "w")) == NULL)
    {
      cerr << "Cannot open " << FileName << " for writing.\n";
      exit(1);
    }
  // First, write out preliminaries.
  fprintf (fout, "%s\n", FileName);
  fprintf (fout, "%d ", DavidsGrid.NumPoints);
  fprintf (fout, "%1.16e %1.16e %1.16e \n", 
	   DavidsGrid.Z, DavidsGrid.r0, DavidsGrid.Spacing);
  // Next 3 parameters aren't used
  //fprintf (fout, "1.0 1.0 1.0\n");
  fprintf (fout, "%d ", NumRadialWFs);
  for (int n=0; n<NumRadialWFs; n++)
    fprintf (fout, "%d ", RadialWFs(n).l);
  fprintf (fout, "\n");
  
  for (int i=0; i<DavidsGrid.NumPoints; i++)
    {
      double r = DavidsGrid(i);
      double a, dadr, b, Vion;
      if (r <= AtomGrid.End)
	{
	  double A,B,Vtot,dAdr;
	  PH->ABV(r,A,B,Vtot,dAdr);
	  
	  a = A - 1.0;
	  dadr = dAdr;
	  b = B - a - 1.0;
	  double VHartree = Hartree (r);
	  double VExCorr = ExCorr(r);
	  //if (DoSelfConsistent)
	  Vion = Vtot - VHartree - VExCorr;
	  //else
	  //Vion = Vtot;
	}
      else
	{
	  a = 0.0; dadr = 0.0; b = 0.0;
	  Vion = -/*IonCharge*/PH->Zion / r;
	  // HACK
	  //Vion = 0.0;
	}

      fprintf (fout, "%1.16e %1.16e %1.16e %1.16e %1.16e ", 
	       r, b, Vion, a, dadr);
      if (r <= AtomGrid.End)
	for (int j=0; j<NumRadialWFs; j++)
	  // HACK
	  if (RadialWFs(j).l == 2)
	    fprintf (fout, "%1.16e ", RadialWFs(j)(r)); // /r);
	  else
	    fprintf (fout, "%1.16e ", RadialWFs(j)(r));
      else
	for (int j=0; j<NumRadialWFs; j++)
	  fprintf (fout, "0.0 ");
      fprintf (fout, "\n");	    
    } 
  fclose (fout);
}



PH_CubicSplineXC *
Atom::CalcPH_XC(Grid *Agrid, Grid *Bgrid, Grid *Vgrid)
{
  if (fabs(Agrid->End - PH->CoreRadius) > 1.0e-10)
    {
      cerr << "Agrid does not match CoreRadius in CalcPH_XC\n";
      exit(1);
    }
  if (fabs(Bgrid->End - PH->CoreRadius) > 1.0e-10)
    {
      cerr << "Bgrid does not match CoreRadius in CalcPH_XC\n";
      exit(1);
    }

  PH_CubicSplineXC &PHXC = *(new PH_CubicSplineXC);
  
  Array<double, 1> Ainit(Agrid->NumPoints-1);
  Array<double, 1> Binit(Bgrid->NumPoints-2);
  Array<double, 1> Vinit(Vgrid->NumPoints-1);

  for (int i=0; i<Agrid->NumPoints-1; i++)
    {
      double r = (*Agrid)(i);
      double A, B, V, dAdr;
      PH->ABV(r, A, B, V, dAdr);
      Ainit(i) = sqrt(A);
    }

  for (int i=0; i<Bgrid->NumPoints-2; i++)
    {
      double r = (*Bgrid)(i+1);
      double A, B, V, dAdr;
      PH->ABV(r, A, B, V, dAdr);
      Binit(i) = sqrt(B);
    }

  CalcChargeDensity();
  CalcHartreePot();
  CalcExCorrPot();

  double Zion;
  for (int i=0; i<NumRadialWFs; i++)
    Zion += RadialWFs(i).Occupancy;

  Grid &AtomGrid = *RadialWFs(0).grid;
  for (int i=0; i<Vgrid->NumPoints-1; i++)
    {
      double r = (*Vgrid)(i);
      double A, B, V, dAdr;

      if (r < AtomGrid.End)
	{
	  PH->ABV(r, A, B, V, dAdr);
	  double VHartree = Hartree (r);
	  double VExCorr = ExCorr(r);
	  double Vion = V - VHartree - VExCorr;
	  Vinit(i) = Vion;
	  // HACK
	  //Vinit(i) = V;
	}
      else
	Vinit(i) = -Zion/r;
    }


  PHXC.Init(Ainit, Binit, Vinit, Agrid, Bgrid, Vgrid, Zion, PH->CoreRadius);

  return (&PHXC);

}


void RadialWF::Write(IOSectionClass &out)
{
  out.WriteVar ("IsRelativistic", IsRelativistic);
  out.WriteVar ("Energy", Energy);
  out.WriteVar ("l", l);
  out.WriteVar ("CoreNodeNum", CoreNodeNum);
  out.WriteVar ("DesiredNodeNum", DesiredNodeNum);
  out.WriteVar ("Weight", Weight);
  out.WriteVar ("CutOffRadius", CutOffRadius);
  out.WriteVar ("Occupancy", Occupancy);
  out.WriteVar ("MaxCurve", MaxCurve);
  out.WriteVar ("Label", Label);
  Array<double,1> uVec(grid->NumPoints), duVec(grid->NumPoints);
  for (int i=0; i<grid->NumPoints; i++) {
    uVec(i) = u(i);
    duVec(i) = dudr(i);
  } 
  out.WriteVar ("u", uVec);
  out.WriteVar ("dudr", duVec);
}


// Nota bene:  grid must be set before calling this function
void RadialWF::Read (IOSectionClass &in)
{
  assert(in.ReadVar("IsRelativistic", IsRelativistic));
  assert(in.ReadVar("Energy", Energy));
  assert(in.ReadVar("l", l));
  assert(in.ReadVar("CoreNodeNum", CoreNodeNum));
  assert(in.ReadVar("DesiredNodeNum", DesiredNodeNum));
  assert(in.ReadVar("Weight", Weight));
  assert(in.ReadVar("CutOffRadius", CutOffRadius));
  assert(in.ReadVar("Occupancy", Occupancy));
  assert(in.ReadVar("MaxCurve", MaxCurve));
  assert(in.ReadVar("Label", Label));
  Array<double,1> uVec, duVec;
  assert(in.ReadVar ("u", uVec));
  assert(in.ReadVar ("dudr", duVec));
  u.Init(grid, uVec);
  dudr.Init(grid, duVec);
}


void Atom::Write(IOSectionClass &out)
{
  out.NewSection("PH");
  PH->Write(out);
  out.CloseSection();

  out.NewSection("Grid");
  grid->Write(out);
  out.CloseSection();

  for (int i=0; i<NumRadialWFs; i++) {
    out.NewSection ("RadialWF");
    RadialWFs(i).Write(out);
    out.CloseSection();
  }
}


void Atom::Read(IOSectionClass &in)
{
  in.OpenSection("PH");
  PH = ReadPH (in);
  in.CloseSection();

  in.OpenSection("Grid");
  grid = ReadGrid (in);
  in.CloseSection();

  NumRadialWFs = in.CountSections("RadialWF");
  RadialWFs.resize(NumRadialWFs);
  for (int i=0; i<NumRadialWFs; i++) {
    RadialWFs(i).PH = PH;
    RadialWFs(i).grid = grid;
    in.OpenSection ("RadialWF", i);
    RadialWFs(i).Read(in);
    in.CloseSection();
  }
}
