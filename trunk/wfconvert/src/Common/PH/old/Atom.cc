#include "Atom.h"
#include "../DFT/Functionals.h"

//////////////////////////////////////////////////////////////////
//                    Generic Atom Routines                     // 
//////////////////////////////////////////////////////////////////
scalar NormalizationDeriv (scalar r, scalar sum,
			   void *WFptr)
{
  scalar deriv;
  RadialWF &WF = *((RadialWF *)WFptr);
  scalar u = WF(r);
  deriv = u*u;
  return (deriv);
}


void
RadialWF::Normalize()
{
  Array<scalar,1> Temp(grid->NumPoints);
  Temp(0) = 0.0;
  IntegrateFirstOrderNS(*grid, 0, grid->NumPoints-1, Temp,
			NormalizationDeriv, this);

  scalar norm = 
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



scalar NormDeriv (scalar r, scalar sum,
		  void *WFptr)
{
  RadialWF &WF = *((RadialWF *)WFptr);
  scalar u = WF(r);
  return(u*u);
}



void
Atom::CalcPartialNorms()
{
  Array<scalar,1> Temp(grid->NumPoints);
  for (int i=0; i<NumRadialWFs; i++)
    {
      Temp(0) = 0.0;
      IntegrateFirstOrder(*grid, 0, grid->NumPoints-1, Temp, 
			  NormDeriv, &(RadialWFs(i)));
      scalar uend = RadialWFs(i).u.Params(grid->NumPoints-1);
      RadialWFs(i).PartialNorm = Temp(grid->NumPoints-1) / (uend*uend);
    }
}



Vec2 ModeDerivs (scalar r, Vec2 Sum,
		 void *Atomptr)
{

  Atom &atom = *(Atom *)Atomptr;
  scalar Charge = 0.0;
  for (int i=0; i<atom.NumRadialWFs; i++)
    {
      scalar u = atom.RadialWFs(i)(r);
      scalar uend = atom.RadialWFs(i).u.Params(atom.grid->NumPoints-1);
      Charge += atom.RadialWFs(i).Occupancy*u*u/(uend*uend);
    }
  Vec2 derivs;
  derivs[0] = r*Charge;
  derivs[1] = r*derivs[0];
  return (derivs);
}


//////////////////////////////////////////////////////////////////
//              PseudoAtom Integration Routines                 //
//////////////////////////////////////////////////////////////////



Vec2 PseudoRadialDerivs(scalar r, Vec2 u_and_du,
			void *WFptr)
{
  Vec2 deriv;
  RadialWF &WF = *((RadialWF *)WFptr);
  PseudoHamiltonian &PH = *(WF.PH);

  // d/dr u = u'
  deriv[0] = u_and_du[1];

  scalar A, B, V, dAdr;

  PH.ABV(r, A, B, V, dAdr);
  
  scalar u = u_and_du[0];
  scalar dudr = u_and_du[1];
  scalar l = (scalar)WF.l;
  scalar E = WF.Energy;

  deriv[1] = 1.0/A * 
    (-dAdr*dudr + (dAdr/r + l*(l+1.0)*B/(r*r) + 2.0*(V-E))*u);

  return (deriv);
}



Vec2 RelativisticRadialDerivs(scalar r, Vec2 u_and_du,
			      void *WFptr)
{
  RadialWF &WF = *((RadialWF *)WFptr);
  PseudoHamiltonian &PH = *(WF.PH);

  const scalar alpha = 1.0/137.036;
  const scalar kappa = -1.0;
  Vec2 derivs;
  scalar u, dudr;
  
  scalar V, dVdr;

  //V = (*PH.FullCoreV)(r);
  //dVdr = PH.FullCoreV->Deriv(r);
  V = PH.V(r);
  dVdr = PH.dVdr(r);

  u = u_and_du[0];
  dudr = u_and_du[1];

  derivs[0] = dudr;
  scalar E = WF.Energy;
  scalar M = 1.0 - alpha*alpha*0.5*(V-E);
  scalar l = (scalar) WF.l;

  derivs[1] = (l*(l+1.0)/(r*r) + 2.0 * M *(V-E)) * u
    - 0.5*alpha*alpha/M*dVdr*(dudr + u*kappa/r);
  
  return (derivs);
}



// Setup boundary conditions at the origin
void 
RadialWF::OriginBC(scalar r0, scalar &u0, scalar &du0)
{
  if (IsRelativistic)
    {
      const scalar alpha = 1.0/137.036;
      scalar Z = PH->Z;
      scalar a2Z2 = alpha*alpha*(scalar)Z*(scalar)Z;

      scalar sl = (scalar) l;
      scalar gamma = 
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
      du0 = (scalar)(l+1)/r0;
      if (l == 0)
	{
	  scalar E = Energy;
	  scalar A0, B0, V0, dAdr0;
	  scalar A, B, V, dAdr;
	  PH->ABV(r0, A, B, V, dAdr);
	  PH->ABV(0.0, A0, B0, V0, dAdr0);
	  if (V0 > Energy)
	    {
	      scalar kappa = sqrt(2.0 * (V-E)/A);
	      u0 = sinh(kappa * r0);
	      du0 = kappa * cosh (kappa * r0);
	    }
	  else
	    {
	      scalar k = sqrt(2.0 * (E-V)/A);
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
RadialWF::InfinityBC(scalar rend, scalar &uend, scalar &duend)
{
  scalar Aend,Bend,Vend,dAdrend;
  PH->ABV(rend, Aend,Bend,Vend,dAdrend);
  scalar dVdrend = PH->dVdr(rend);
  scalar E = Energy;
  scalar k = sqrt(l*(l+1.0)/(rend*rend) + (Vend-E));
  uend = 1.0;
  duend = -uend * (k + 0.5*rend/k*(dVdrend - 2.0*l*(l+1.0)/(rend*rend*rend)));
}




void
RadialWF::IntegrateOutward()
{
  Vec2 (*DerivFunc)(scalar r, Vec2 u_and_du, void *WFptr);

  if (IsRelativistic)
    DerivFunc = RelativisticRadialDerivs;
  else
    DerivFunc = PseudoRadialDerivs;

  Array<Vec2,1> Temp(grid->NumPoints);

  // Set up initial conditions:
  scalar r0 = (*grid)(0); 
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
      scalar uend = RadialWFs(i).u.Params(grid->NumPoints-1);
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
  scalar sign = u.Params(0);
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
      scalar r = (*grid)(index);
      
      scalar A,B,V,dAdr;
      PH->ABV(r,A,B,V,dAdr);
      scalar E = Energy;
      scalar sl = (scalar) l;
      
      scalar Deriv2 = dAdr/r + sl*(sl+1.0)*B/(r*r) + 2.0*(V-E);
      
      if (Deriv2 < 0.0)
	done = 1;
      else
	index--;
    }
  if (index == -1)
    index = grid->ReverseMap(0.5*(*grid)(grid->NumPoints-1));

  return (index);
}




scalar 
RadialWF::IntegrateInwardOutward(int &TIndex)
{
  Vec2 (*DerivFunc)(scalar r, Vec2 u_and_du,
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
  scalar r0 = (*grid)(0);
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
  scalar rend = (*grid)(EndPoint);

  InfinityBC(rend, Temp(EndPoint)[0], Temp(EndPoint)[1]);
  
  // Do integration from the right to the turning point
  IntegrateSecondOrder(*grid,EndPoint, TIndex, Temp, DerivFunc, this);  
  
  if ((u.Params(TIndex)*Temp(TIndex)[0]) < 0.0)
    for (int i=TIndex; i<=EndPoint; i++)
      {
	Temp(i)[0] *= -1.0;
	Temp(i)[1] *= -1.0;
      }

  scalar CuspValue = Temp(TIndex)[1]/Temp(TIndex)[0] -
    dudr.Params(TIndex)/u.Params(TIndex);

  // Copy values in, normalizing to make R continuous
  scalar factor = u.Params(TIndex) / Temp(TIndex)[0];
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
  const scalar Tolerance = 1e-12;
  scalar Ehigh, Elow, Etrial;
  Ehigh = 0.0;
  
  scalar r = (*grid)(0);
  scalar A,B,V,dAdr;
  PH->ABV(r,A,B,V,dAdr);

  if (IsRelativistic)
    {
      scalar n = DesiredNodeNum + 1.0;
      Elow = -1.5*PH->Z*PH->Z/(n*n);
    }
  else
    {
      // Eigenenergy can't be lower than lowest potential -- otherwise
      // we'd have no turning point.
      Elow = PH->V((*grid)(0));
      for (int i=1; i<grid->NumPoints; i++)
	{
	  scalar r = (*grid)(i);
	  if (PH->V(r) < Elow)
	    Elow = PH->V(r);
	}
    }

  // HACK
  //Elow = -300.0;
  
  while ((Ehigh-Elow) > Tolerance)
    {
      scalar Etrial = 0.5 * (Ehigh + Elow);
      Energy = Etrial;
      scalar CuspValue = IntegrateInwardOutward(Tindex);
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
  const scalar Tolerance = 1e-12;
  int Tindex;

  IntegrateInwardOutward(Tindex);
  int NumNodes = CountNodes();


  scalar Ehigh, Elow, Etrial;
  Ehigh = 0.0;
  
  if (IsRelativistic)
    {
      scalar n = DesiredNodeNum + 1.0;
      scalar Z = PH->FullCoreV->Z;
      Elow = -1.5*Z*Z/(n*n);
    }
  else
    {
      // Eigenenergy can't be lower than lowest potential -- otherwise
      // we'd have no turning point.
      Elow = PH->V((*grid)(0));
      for (int i=1; i<grid->NumPoints; i++)
	{
	  scalar r = (*grid)(i);
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

  scalar Eold;
  Eold = Etrial = Energy;
  int done = 0;
  scalar A, B, V, dAdr;
  while (!done)
    {
      scalar CuspValue = IntegrateInwardOutward(Tindex);
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
      scalar C = 0.5 * A * CuspValue;
      scalar u0 = u.Params(Tindex);
      scalar Etest = Eold - 4.0*M_PI*C*u0*u0;
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


scalar ChargeDensityDeriv(scalar r, scalar sum, void *atomptr)
{
  Atom &atom = *(Atom *)atomptr;
  return (4.0*M_PI * r * r *atom.ChargeDensity(r));
}

scalar
Atom::TotalCharge()
{
  Array<scalar,1> temp(grid->NumPoints);
  temp(0) = 0.0;
  IntegrateFirstOrderNS(*grid, 0, grid->NumPoints-1, temp,
			ChargeDensityDeriv, this);
  return (temp(grid->NumPoints-1));
}


void
Atom::CalcChargeDensity()
{
  Array<scalar,1> Rho(grid->NumPoints);

  Rho = 0.0;

  NormalizeRadialWFs();
  for (int i=0; i<grid->NumPoints; i++)
    {
      Rho(i) = 0.0;
      scalar r = (*grid)(i);
      scalar rinv2 = 1.0/(r*r);
      for (int n=0; n<NumRadialWFs; n++)
	{
	  scalar u = RadialWFs(n)(r);
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



scalar Hartree1 (scalar r, scalar N,
		 void *AtomPtr)
{
  Atom &atom = *((Atom *)AtomPtr);

  scalar rho= atom.ChargeDensity(r);
  return(rho*r*r); // dV/dr = V'
}

scalar Hartree2 (scalar r, scalar N,
		 void *AtomPtr)
{
  Atom &atom = *((Atom *)AtomPtr);
  scalar rho= atom.ChargeDensity(r);
  return(rho*r); // dV/dr = V'
}




void
Atom::CalcHartreePot ()
{
  int NumPoints = grid->NumPoints;
  int EndPoint = NumPoints-1;
  Array<scalar,1> term1(NumPoints), term2(NumPoints);
  Array<scalar,1> HV(NumPoints);

  term1(0) = 0.0;
  IntegrateFirstOrder(*grid, 0, EndPoint, term1, Hartree1, this);
  term2(0) = 0.0;
  IntegrateFirstOrder(*grid, 0, EndPoint, term2, Hartree2, this);
  
  scalar term2infty = term2(EndPoint);
  for (int i=0; i<=EndPoint; i++)
    {
      scalar r = (*grid)(i);
      HV(i) = 4.0*M_PI*(term1(i)/r + (term2infty-term2(i)));
    }
  Hartree.Init(grid, HV, 5.0e30, 5.0e30);
}

void
Atom::CalcExCorrPot ()
{
  int NumPoints = grid->NumPoints;
  Array<scalar,1> VXC(NumPoints);

  for (int i=0; i<NumPoints; i++)
    {
      scalar Vxcup, Vxcdown;
      scalar UpCharge = 0.5 * ChargeDensity.Params(i);
      scalar DownCharge = UpCharge;
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

  Array<scalar,1> VHXC(grid->NumPoints);
  for (int i=0; i<grid->NumPoints; i++)
    VHXC(i) = Hartree.Params(i) + ExCorr.Params(i);

  scalar NumElecs = 0.0;
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

  //scalar IonCharge=0.0;
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
      scalar r = DavidsGrid(i);
      scalar a, dadr, b, Vion;
      if (r <= AtomGrid.End)
	{
	  scalar A,B,Vtot,dAdr;
	  PH->ABV(r,A,B,Vtot,dAdr);
	  
	  a = A - 1.0;
	  dadr = dAdr;
	  b = B - a - 1.0;
	  scalar VHartree = Hartree (r);
	  scalar VExCorr = ExCorr(r);
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
  
  Array<scalar, 1> Ainit(Agrid->NumPoints-1);
  Array<scalar, 1> Binit(Bgrid->NumPoints-2);
  Array<scalar, 1> Vinit(Vgrid->NumPoints-1);

  for (int i=0; i<Agrid->NumPoints-1; i++)
    {
      scalar r = (*Agrid)(i);
      scalar A, B, V, dAdr;
      PH->ABV(r, A, B, V, dAdr);
      Ainit(i) = sqrt(A);
    }

  for (int i=0; i<Bgrid->NumPoints-2; i++)
    {
      scalar r = (*Bgrid)(i+1);
      scalar A, B, V, dAdr;
      PH->ABV(r, A, B, V, dAdr);
      Binit(i) = sqrt(B);
    }

  CalcChargeDensity();
  CalcHartreePot();
  CalcExCorrPot();

  scalar Zion;
  for (int i=0; i<NumRadialWFs; i++)
    Zion += RadialWFs(i).Occupancy;

  Grid &AtomGrid = *RadialWFs(0).grid;
  for (int i=0; i<Vgrid->NumPoints-1; i++)
    {
      scalar r = (*Vgrid)(i);
      scalar A, B, V, dAdr;

      if (r < AtomGrid.End)
	{
	  PH->ABV(r, A, B, V, dAdr);
	  scalar VHartree = Hartree (r);
	  scalar VExCorr = ExCorr(r);
	  scalar Vion = V - VHartree - VExCorr;
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
