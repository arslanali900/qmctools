#include "Cost.h"


//////////////////////////////////////////////////////////////////
//                    CostFunction Routines                     //
//////////////////////////////////////////////////////////////////

scalar NegativityDeriv (scalar r, scalar sum, 
				  void *WFptr)
{
  RadialWF &WF = *((RadialWF *)WFptr);
  scalar u = WF(r);
  return(fabs(u)-u);
}



scalar CurvatureDeriv (scalar r, scalar sum, void *WFptr)
{
  RadialWF &WF = *(RadialWF *)WFptr;
  
  scalar uend = WF.u.Params(WF.grid->NumPoints-1);
  scalar d2udr2 = WF.dudr.Deriv(r);
  scalar curvature = fabs((d2udr2*d2udr2) / (uend*uend));
  scalar diff = curvature - WF.MaxCurve;
  scalar over = 0.5 * (diff + fabs(diff));
  //  return ((d2udr2*d2udr2) / (uend*uend));
  return (over * over);
}


scalar Curvature(RadialWF &WF)
{  
  Array<scalar,1> curve(WF.grid->NumPoints);
  curve(0) = 0.0;
  IntegrateFirstOrder(*WF.grid, 0, WF.grid->NumPoints-1, curve,
		      CurvatureDeriv, &WF);
  return (curve(WF.grid->NumPoints-1));
}


typedef struct
{
  RadialWF *WF1, *WF2;
} WFpair;


scalar WFMatchDeriv(scalar r, scalar sum, void *WFpairptr)
{
  WFpair &TwoWFs = *(WFpair *)WFpairptr;
  RadialWF &WF1 = *(RadialWF *) TwoWFs.WF1;
  RadialWF &WF2 = *(RadialWF *) TwoWFs.WF2;

  scalar diff = WF1(r) - WF2(r);
  return ((diff * diff) * WF1.WeightFunc(r));
}


scalar PHCost::WFMatchCost()
{
  WFpair pair;
  scalar cost = 0.0;
  Array<scalar,1> Temp(Patom->grid->NumPoints);

  for (int n=0; n < Patom->NumRadialWFs; n++)
    {
      pair.WF1 = &Patom->RadialWFs(n);
      pair.WF2 = &FCatom->RadialWFs(n);
      Temp(0) = 0.0;
      IntegrateFirstOrderNS (*(Patom->grid), 0, Patom->grid->NumPoints-1, 
			     Temp, WFMatchDeriv, &pair);
      cost += pair.WF1->Weight * Temp(Patom->grid->NumPoints-1);
    }
  return (cost);
}




scalar PHCost::Cost()
{
  Patom->IntegrateRadialEquations();

  scalar cost = 0.0;
  Grid &PseudoGrid = *(Patom->RadialWFs(0).grid);
  Grid &FCGrid = *(FCatom->RadialWFs(0).grid);
  // First, ensure that the various u_l's are positive
  
  Array<scalar,1> PseudoTemp(Patom->RadialWFs(0).grid->NumPoints);
  Array<scalar,1> FCTemp(FCatom->RadialWFs(0).grid->NumPoints);

  scalar NegativityCost = 0.0;
  for (int l=0; l<Patom->NumRadialWFs; l++)
    {
      PseudoTemp(0) = 0.0;
      IntegrateFirstOrder(PseudoGrid, 0, PseudoGrid.NumPoints-1, PseudoTemp,
			  NegativityDeriv, &(Patom->RadialWFs(l)));
      NegativityCost += NegativityCoef * PseudoTemp(PseudoGrid.NumPoints-1);
    }

  scalar CurveCost = 0.0;
  for (int l=0; l<Patom->NumRadialWFs; l++)
    {
      //Patom->RadialWFs(l).MaxCurve = MaxCurve;
      CurveCost += CurveCoef * Patom->RadialWFs(l).Weight *
	Curvature(Patom->RadialWFs(l));
    }

  scalar NodeCost = 0.0;
  for (int l=0; l<Patom->NumRadialWFs; l++)
    {
      if (Patom->RadialWFs(l).CountNodes() !=
	  Patom->RadialWFs(l).CoreNodeNum)
	NodeCost+=NodeCoef*Patom->RadialWFs(l).Weight;
    }

  // Add a term to the cost to keep reasonable values for A, B, and V
  // potentials
  scalar OverMaxCost=0.0;
  for (int i=0; i<100; i++)
    {
      scalar r,A,B,V,dAdr;
      r = (scalar)i/99.0 * Patom->PH->CoreRadius;
      Patom->PH->ABV(r, A, B, V, dAdr);
      if (A > Patom->PH->Amax)
	OverMaxCost += OverMaxCoef * (fabs(A)-Patom->PH->Amax);
      if (fabs(B) > Patom->PH->Bmax)
	OverMaxCost += OverMaxCoef * (fabs(B)-Patom->PH->Bmax);
      if (fabs(V) > Patom->PH->Vmax)
	OverMaxCost += OverMaxCoef * (fabs(V)-Patom->PH->Vmax);
    }

  // Next ensure that their logorithmic derivatives match.
  scalar CoreRadius = Patom->PH->CoreRadius;
  scalar LogDerivativesCost = 0.0;
  for (int l=0; l<Patom->NumRadialWFs; l++)
    {
      scalar PseudoLogDeriv = Patom->RadialWFs(l).Deriv(CoreRadius) /
	(Patom->RadialWFs(l))(CoreRadius);
      scalar FCLogDeriv = FCatom->RadialWFs(l).Deriv(CoreRadius) /
	(FCatom->RadialWFs(l))(CoreRadius);
      scalar diff = PseudoLogDeriv - FCLogDeriv;
      LogDerivativesCost += LogDerivCoef * Patom->RadialWFs(l).Weight 
	* diff * diff;
    }
  
  // Next ensure that their partial norms match.
  scalar PartialNormsCost = 0.0;
  for (int l=0; l<Patom->NumRadialWFs; l++)
    {
      scalar diff = Patom->RadialWFs(l).PartialNorm - 
	FCatom->RadialWFs(l).PartialNorm;
      PartialNormsCost += PartialNormCoef * Patom->RadialWFs(l).Weight 
	* diff * diff;
    }

  scalar MatchCost = MatchCoef * WFMatchCost();

  /*
  cerr << "NegativityCost = " << NegativityCost << ".\n";
  cerr << "LogDerivativesCost = " << LogDerivativesCost << ".\n";
  cerr << "PartialNormCost = " << PartialNormsCost << ".\n";
  cerr << "NodeCost = " << NodeCost << ".\n";
  cerr << "OverMaxCost = " << OverMaxCost << ".\n";
  cerr << "CurveCost = " << CurveCost << ".\n";
  cerr << "MatchCost = " << MatchCost << ".\n";
  */

  cost = NegativityCost + LogDerivativesCost + PartialNormsCost +
    OverMaxCost + NodeCost + CurveCost + MatchCost;

  /*cerr << "log10(cost) = " << log10(cost) << "\n"; */
  if (UseLog)
    return (log10(cost));
  else
    return (cost);
}


