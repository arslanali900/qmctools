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
#include "../MPI/Communication.h" //include not needed
#include "U_l.h"
#include "../Integration/GKIntegration.h"
#include "../Integration/HermiteQuad.h"
#include "../Distributed/DistributedMat.h"

void
U_l::Initialize(int l_, double lambda_, double FinalBeta, int NumSquares, 
		Grid *grid_, CoreTransform *transform_, 
		Potential *Pot_, 
		CommunicatorClass groupComm) 
{
  lambda = lambda_;
  l = l_;
  //beta = FinalBeta * pow(0.5, NumSquares);
  beta = ldexp (FinalBeta, -NumSquares);
  grid = grid_;
  transform = transform_;
  Pot = Pot_;
  GroupComm = groupComm;

  NumPoints = grid->NumPoints;
  xMax = (*grid)(NumPoints-1);
  xMaxSC = xMax+20.0*sqrt(lambda*FinalBeta);

  // Initialize the SemiClassical (SC) splines
  // Set the spacing to be the maximum spacing of the U grid
  double spacing = (*grid)(NumPoints-1) - (*grid)(NumPoints-2);
  int SCPoints = (int)ceil((xMaxSC - xMax)/spacing) + 1;
  SCgrid.Init(xMax, xMaxSC, SCPoints);

  Array<double,2> SCMat(SCPoints, NumPoints);
  SCMat = 0.0;
  SCspline.Init(&SCgrid, grid, SCMat);

  Array<double,2> TempMat(NumPoints, NumPoints);  
  TempMat = 0.0;
  Uspline.Init(grid,grid,TempMat);
  dUspline.Init(grid,grid,TempMat);  
}



void U_l::FillWithSC()
{
  int SCPoints = SCgrid.NumPoints;
  DistributedMat SCMat(SCPoints, NumPoints, GroupComm);

  if (COMM::WorldProc() == 0)
    cerr << "Initializing Semiclassical splines for l = " << l << ":\n"; 
  // Put elemenents I'm responsible for the the matrix
  for (int i = 0; i< SCMat.MyNumElements(); i++) {
    int row, col;
    SCMat.MyElement(i, row, col);
    double x = (*grid)(col);
    double xp = SCgrid(row);
    SCMat(row,col) = USC(x,xp)/beta;
  }
  // Now do collective communication gathering elements
  SCMat.AllGather();

  // Now copy into the splines
  for (int row=0; row<SCPoints; row++)
    for (int col=0; col<NumPoints; col++)
      SCspline(row,col) = SCMat(row,col);

  if (COMM::WorldProc() == 0)
    cerr << "Initializing U and dU splines for l = " << l << ":\n";

  DistributedSymmMat UMat(NumPoints, GroupComm);

  // Put elements I am responsible in the matrix
  for (int i = 0; i< UMat.MyNumElements(); i++) {
    int row, col;
    UMat.MyElement(i, row, col);
    double x = (*grid)(col);
    double xp =(*grid)(row);
    UMat(row,col) = USC(x,xp);//Primitive(x,xp);
  }
  // Now do collective gather
  UMat.AllGather();
  // Now copy into the splines
  double betainv = 1.0/beta;
  for (int row=0; row<NumPoints; row++)
    for (int col=0; col<NumPoints; col++) {
      Uspline(row,col) = UMat(row,col);
      dUspline(row,col) = betainv * UMat(row,col);
    }
  
  // Fill in Yl;
  Yl.resize(NumPoints);
  for (int row=0; row<NumPoints; row++)
    Yl(row) = Y_l((*grid)(row));
}









void U_l::CalcElement(int row, int col, 
		      double &U, double &dU)
{
  USquaringIntegrand UIntegrand;
  dUSquaringIntegrand dUIntegrand;
  UIntegrand.Ul = this;

//   if (beta >= 8.0) {
//     absTol *= 100.0;
//     relTol *= 100.0;
//   }

  double sigma = sqrt(lambda*beta);
  double x = (*grid)(row);
  double xp = (*grid)(col);
  double xbar = 0.5*(x+xp);
  // This expression is valid for large l and small (x+xp)
  // The centripetal term pushes the center of the free-particle
  // integrand outward.
  double center_push = sqrt(2.0*lambda*beta*(double)(l+1));
  // center is the center of the free-particle part of the integrand
  // Make sure our center is far enough out.
  double center = max(xbar, center_push);

  double Start = center - 9.0*sigma;
  if (Start < 1.0e-10)
    Start = 1.0e-10;
  double End = center + 18.0*sigma;
  double deltaStart = Start - xbar;
  double deltaEnd = End - xbar;

  UIntegrand.row = row;
  UIntegrand.col = col;
  UIntegrand.SetScaleExp();  
  //cerr << "Doing U:\n";
  GKIntegration<USquaringIntegrand,GK31> U_gk(UIntegrand);
  U_gk.SetRelativeErrorMode();
//   double Uval = U_gk.Integrate(deltaStart, deltaEnd, Tolerance);
//   U = -(log(Uval)-UIntegrand.ScaleExp);
  if (Use_m1) {
    double expmUm1 = 
      U_gk.Integrate(deltaStart, deltaEnd, AbsTol*beta, RelTol, false);
    U = -log1p(expmUm1);
  }
  else {
    double expmU = U_gk.Integrate(deltaStart, deltaEnd, AbsTol, RelTol, false);
    U = -log(expmU);
  }

  double alpha;
  dUIntegrand.Ul = this;
  dUIntegrand.row = row;
  dUIntegrand.col = col;
  dUIntegrand.SetScaleExp();
  GKIntegration<dUSquaringIntegrand,GK31> dU_gk(dUIntegrand);  
  // HACK  
  //dU_gk.SetRelativeErrorMode();
  //dU_gk.SetAbsoluteErrorMode();
  //cerr << "Doing dU:\n";
  // We must allow both relative and absolute tolerances since we can
  // have dUvals that are very close to zero.  dU/dbeta should be of
  // order unity independent of beta.
  double dUval = 
    dU_gk.Integrate(deltaStart, deltaEnd, 50.0*AbsTol, 50.0*RelTol, false)*
    exp(-dUIntegrand.ScaleExp);

  dU = -exp(U) * dUval;
}


void U_l::CalcElement_Hermite(int row, int col, 
			      double &U, double &dU)
{
  double sigma = sqrt(lambda*beta);
  double x = (*grid)(row);
  double xp = (*grid)(col);

  USquaringIntegrand_Hermite UIntegrand;
  dUSquaringIntegrand_Hermite dUIntegrand;

  UIntegrand.Ul = this;
  UIntegrand.row = row;
  UIntegrand.col = col;
  UIntegrand.SetScaleExp();  

  dUIntegrand.Ul = this;  
  dUIntegrand.row = row;
  dUIntegrand.col = col;  
  dUIntegrand.SetScaleExp();  

  HermiteQuadClass<Hermite30, USquaringIntegrand_Hermite> UIntegrator;
  HermiteQuadClass<Hermite30, dUSquaringIntegrand_Hermite> dUIntegrator;

  UIntegrator.SetSigma(sigma);
  dUIntegrator.SetSigma(sigma);

//   double Uval = UIntegrator.Integrate(UIntegrand);
//   U = -(log(Uval)-UIntegrand.ScaleExp);

  if (Use_m1) {
    double expmUm1 = UIntegrator.Integrate(UIntegrand);
    U = -log1p(expmUm1);
  }
  else {
    double expmU = UIntegrator.Integrate(UIntegrand);
    U = -log(expmU);
  }


  //cerr << "Doing dU:\n";
  double dUval = dUIntegrator.Integrate(dUIntegrand)*
    exp(-dUIntegrand.ScaleExp);

  dU = -exp(U) * dUval;
}







void CompletionBar(double frac)
{
  const int cols = 60;
  static int OldNumEquals = -1;
  static int OldPercent = -1;
  static bool FirstTime = true;
  
  int NumEquals = (int)ceil(frac*(cols-3));
  if (NumEquals < OldNumEquals)
    FirstTime = true;
  int Percent = (int)floor (frac * 100.0 + 0.5);
  if ((NumEquals != OldNumEquals) || (Percent != OldPercent))
    {
      OldNumEquals = NumEquals;      
      OldPercent = Percent;
      if (FirstTime)
	FirstTime = false;
      else
	fprintf (stderr, "%c[2A\n", 27);      
      fprintf (stderr, "[");

      int NumSpaces = (cols-3-NumEquals);
      for (int i=0; i<NumEquals; i++)
	fprintf (stderr, "=");
      fprintf (stderr, ">");
      for (int i=0; i<NumSpaces; i++)
	fprintf (stderr, " ");
      fprintf (stderr, "]");

      fprintf (stderr, " %3d%c\n", Percent,'%');
    }      
}



void
U_l::Square()
{
  // Clamp high actions to prevent underflow
  for (int i=0; i<NumPoints; i++)
    for (int j=0; j<NumPoints; j++)
      if (U(i,j) > 200.0)
	U(i,j) = 200.0;

  CommunicatorClass WorldComm;
  WorldComm.SetWorld();
  // Divide work up between nodes
  DistributedSymmMat UMat(NumPoints, GroupComm), 
    dUMat(NumPoints, GroupComm);
  int lastrow = -1;
  for (int i=0; i<UMat.MyNumElements(); i++) {
    int row, col;
    UMat.MyElement(i,row,col);
    if (row != lastrow && WorldComm.MyProc() == 0) {
      lastrow = row;
      cerr << "   row = " << row << endl;
    }
    // This ensures that we don't get exactly -1 for expm1(-U).
    // Also ensure precision at lower temperatures.
    Use_m1 = (beta < 1.0e-3) && (U(row,col) < 10.0);


    double x = (*grid)(row);
    double xp = (*grid)(col);
    double xbar = 0.5*(x+xp);
    
    //cerr << "row = " << row << " col = " << col << endl;
    //      if (xbar > (7.0*sqrt(beta)))
    double xmin = sqrt(2.0*lambda*beta*(double)(l+1));
    
    if ((xbar > (9.0*sqrt(2.0*lambda*beta))) && (xbar > 3.0*xmin))
      CalcElement_Hermite(row,col, UMat(row,col), dUMat(row,col));
    else
      CalcElement (row, col, UMat(row,col), dUMat(row,col));
    // Output a visual indicator of how far we are
    double FracDone = (double)i/(double)UMat.MyNumElements();
    //      if (COMM::WorldProc()==0)
      //CompletionBar(FracDone);
  }
  // Now gather up all the results with collective communication
  UMat.AllGather();
  dUMat.AllGather();

  // Put new data into the permanent locations
  for (int row=0; row<NumPoints; row++)
    for (int col=0; col<NumPoints; col++) {
      U(row,col) = UMat(row,col);
      dU(row,col) = dUMat(row,col);
    }

  UseSC = false;
  beta *= 2.0;
}



double U_l::rho_l(int row, int col)
{
  double x = (*grid)(row);
  double xp = (*grid)(col);
  
  double rho0 = exp(Log_rho0_l(l, x, xp, lambda, beta));
  return (rho0 * exp(-U(row,col)));
}

double U_l::rho_l(double x, double xp)
{
  double rho0 = exp(Log_rho0_l(l,x,xp,lambda,beta));
  double Uval = U(x,xp);
  return (rho0*exp(-Uval));
}


double U_l::drho_l_dbeta(int row, int col)
{
  double x = (*grid)(row);
  double xp = (*grid)(col);
  
  double rho0 = exp(Log_rho0_l(l, x, xp, lambda, beta));
  double drho0 = drho0_dbeta_l(l, x, xp, lambda, beta);
  double Uval = U(row,col);
  double dUval = dU(row,col);
  double drho = exp(-Uval) * (drho0 - rho0*dUval);
  // HACK
  //return (drho0);
  return (drho);
}

double U_l::drho_l_dbeta(double x, double xp)
{
  double rho0 = exp(Log_rho0_l(l, x, xp, lambda, beta));
  double drho0 = drho0_dbeta_l(l, x, xp, lambda, beta);
  double Uval = U(x, xp);
  double dUval = dU(x,xp);
  double drho = exp(-Uval) * (drho0 - rho0*dUval);
  // HACK
  //  return (drho0);
  return (drho);
}





double U_l::H_rho_l(int row, int col)
{
  double x = (*grid)(row);
  double xp = (*grid)(col);

  // First, calculate the second derivative
  double d2rho_dx2;
  double rho0, drho0_dx, d2rho0_dx2;
  double U, dU_dx, d2U_dx2;
  U       = Uspline(col,row);
  dU_dx   = Uspline.Deriv(col,x);
  d2U_dx2 = Uspline.Deriv2(col,x);

  rho0        = exp(Log_rho0_l(l, x, xp, lambda, beta));
  drho0_dx   = drho0_dx_l(l, x, xp, lambda, beta);
  d2rho0_dx2 = d2rho0_dx2_l(l, x, xp, lambda, beta);

  double rho = rho_l(row,col);

  double W = W_l(x);

  d2rho_dx2 = exp(-U)*(-rho0*(d2U_dx2 - dU_dx*dU_dx) 
		       -2.0*drho0_dx*dU_dx + d2rho0_dx2);
  //fprintf (stderr, "analytic: d2U_dx2 = %1.8e\n", d2U_dx2);
  //cerr << "analytic: d2rho_dx2 = " << d2rho_dx2 << "\n";

  return (-0.5*d2rho_dx2 + W*rho);
}

double U_l::dU_Bloche (int row, int col)
{
  double x = (*grid)(row);
  double xp = (*grid)(col);
  double rho0 = exp(Log_rho0_l(l, x, xp, lambda, beta));
  double U = Uspline(row,col);
  double Hrho = H_rho_l(row,col);
  double drho0 = drho0_dbeta_l(l, x, xp, lambda, beta);

  return (1.0/rho0 * (drho0 + exp(U)*Hrho));
}



	  
double U_l::H_rho_l_FD(int row, int col)
{
  double delta = 2.0e-4;

  double x = (*grid)(row);
  double xp = (*grid)(col);

  double rho = rho_l(x,xp);
  double rhoplus = rho_l(x+delta,xp);
  double rhominus = rho_l(x-delta,xp);
  double d2rho_dx2 = 1.0/(delta*delta) * (rhoplus-2.0*rho+rhominus);
  double Uval = Uspline(row,col);
  double W = W_l(x);

  double Uplus = U(x+delta,xp);
  double Uminus = U(x-delta,xp);
  double dU_dx = (Uplus-Uminus)/(2.0*delta);
  double d2U_dx2 = (Uplus-2.0*Uval+Uminus)/(delta*delta);

  //cerr << "FD:       du_dx = " << du_dx << "\n";
  //  fprintf (stderr, "FD:       du2_dx2 = %1.8e\n", d2u_dx2);
  // cerr << "FD:       d2rho_dx2 = " << d2rho_dx2 << "\n";

  return (-0.5*d2rho_dx2 + W*rho);
}





void U_l::Write(IOSectionClass &outSection)
{
  // Write out U matrix
  Array<double,2> temp(NumPoints,NumPoints);
  for (int i=0; i<NumPoints; i++)
    for (int j=0; j<NumPoints; j++)
      temp(i,j) = U(i,j);
  outSection.WriteVar ("U", temp);

  // Write out dU
  for (int i=0; i<NumPoints; i++)
    for (int j=0; j<NumPoints; j++)
      temp(i,j) = dU(i,j);
  outSection.WriteVar ("dU", temp);

  // Write out Yl
  outSection.WriteVar ("Yl", Yl);

  // Write out SCspline
  outSection.NewSection ("SCgrid");
  SCgrid.Write (outSection);
  outSection.CloseSection(); // "SCgrid"
  temp.resize (SCgrid.NumPoints, NumPoints);
  for (int i=0; i<SCgrid.NumPoints; i++)
    for (int j=0; j<NumPoints; j++)
      temp(i,j) = SCspline(i,j);
  outSection.WriteVar ("SC_splines", temp);
}


bool U_l::Read (int this_l, double this_lambda, 
		double this_beta, Grid *Ugrid,
		CoreTransform *this_transform,
		Potential *this_Pot,
		CommunicatorClass comm,
		IOSectionClass &inSection)
{
  l = this_l;
  beta = this_beta;
  lambda = this_lambda;
  grid = Ugrid;
  NumPoints = grid->NumPoints;
  transform = this_transform;
  Pot = this_Pot;
  UseSC = false;
  GroupComm = comm;
  // Read in U
  Array <double,2> temp;
  assert (inSection.ReadVar ("U", temp));
  Uspline.Init (grid, grid, temp);
  assert (inSection.ReadVar ("dU", temp));
  dUspline.Init (grid, grid, temp);
  xMax = grid->End;
  assert (inSection.OpenSection("SCgrid"));
  SCgrid.Read (inSection);
  inSection.CloseSection();
  xMaxSC = SCgrid.End;
  assert (inSection.ReadVar ("SC_splines", temp));
  SCspline.Init(&SCgrid, grid, temp);
  return true;
}
