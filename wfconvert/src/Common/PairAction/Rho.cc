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
#include "Rho.h"
#include "../SpecialFunctions/SpecialFunctions.h"
#include "../nan.h"
#include <unistd.h>


/// This function allocates the U_l's up to lmax and distributes them
/// across the processes.  It does this in such a way as to be as to
/// make the load as uniform as possible.  It currently restricts the
/// number of processes to be an integer multiple of the number of l
/// channels or the number of l channels to be an integer multiple of
/// the number of processes.
void Rho::Initialize(int lmax, double Lambda, scalar FinalBeta, int NumSquares,
		     Grid *newGrid, Potential *newPot, double tailPower,
		     double AbsTol, double RelTol)
{
  lambda = Lambda;
  grid = newGrid;
  Pot = newPot;
  Transform.Initialize(Pot, 1000);
  //beta = FinalBeta * pow(0.5, NumSquares);
  beta = ldexp(FinalBeta, -NumSquares);
  WorldComm.SetWorld();
  int numProcs = WorldComm.NumProcs();
  NumChannels = lmax + 1;
  int myProc = WorldComm.MyProc();
  if (WriteDebug) {
    char debugName[50];
    snprintf (debugName, (size_t)50, "debug.%d", WorldComm.MyProc());
    DebugFile = fopen (debugName, "w");
  }

  if ((numProcs % NumChannels) == 0) {
    // We have 1 or more processors per channel
    int ProcsPerChannel = numProcs/NumChannels;
    ChannelsPerGroup = 1;
    // Now construct processor groups
    MyGroup = myProc/ProcsPerChannel;
    NumGroups = NumChannels;
    ProcsPerGroup = ProcsPerChannel;
    WorldComm.Split(MyGroup, GroupComm);
  }
  else if ((NumChannels % numProcs) == 0) {
    // We have to split up the channels
    ChannelsPerGroup = NumChannels/numProcs;
    MyGroup = myProc;
    NumGroups = numProcs;
    ProcsPerGroup = 1;
    WorldComm.Split(MyGroup, GroupComm);
  }
  else {
    if (WorldComm.MyProc() == 0)
      cerr << "NumProcs (" <<numProcs<< ") is not an integer "
	   << "multiple of NumChannels (" <<NumChannels
	   << ") or vice-versa.\nThis is an inefficient use "
	   << "of CPU time.  Exitting.\n";
    COMM::Finalize();
    exit(1);
  }
  // Now, let's allocate the U_l's
  if (myProc != 0) {
    U_ls.resize(ChannelsPerGroup);
    for (int lindex=0; lindex<ChannelsPerGroup; lindex++) {
      int l = NumGroups*lindex + MyGroup;
      U_ls(lindex).Initialize(l, lambda, FinalBeta, NumSquares,
			      grid, &Transform, Pot, GroupComm);
      U_ls(lindex).TailPower = tailPower;
      U_ls(lindex).AbsTol = AbsTol;
      U_ls(lindex).RelTol = RelTol;
      U_ls(lindex).FillWithSC();
    }
  }
  else {
    U_ls.resize(NumChannels);
    for (int l=0; l<=lmax; l++) {
      U_ls(l).Initialize(l, lambda, FinalBeta, NumSquares,
			 grid, &Transform, Pot, GroupComm);
      U_ls(l).TailPower = tailPower;
      U_ls(l).AbsTol = AbsTol;
      U_ls(l).RelTol = RelTol;
    }
    // Init values for only the U_ls I will do.
    for (int l=0; l<=lmax; l+=NumGroups)
      U_ls(l).FillWithSC();
  }
}

	     

/// This function allocates the U_l's up to lmax.
void Rho::Initialize(int lmax, double Lambda,
		     Grid *newGrid, Potential *newPot)
{
  lambda = Lambda;
  grid = newGrid;
  Pot = newPot;
  Transform.Initialize(Pot, 1000);
  cerr << "x(r=2.0) = " << Transform.r2x(2.0) << endl;
  cerr << "r(x=2.0) = " << Transform.x2r(2.0) << endl;
  WorldComm.SetWorld();
  U_ls.resize(lmax+1);
  for (int l=0; l<=lmax; l++)
    U_ls(l).Initialize(l, lambda, 1.0, 1,
		       grid, &Transform, Pot, WorldComm);
}




/// Collect the various U_l's from the different processor groups.
void 
Rho::GatherU_ls()
{
  int myWorldProc = WorldComm.MyProc();
  cerr << "MyProc = " << myWorldProc << " and I'm in GatherU_ls().\n";
  int myGroupProc = GroupComm.MyProc();
  /// This will be a receive buffer for proc 0 and a send buffer
  /// for all others.
  int NumPoints = grid->NumPoints;
  Array<double,1> Buff((NumPoints*(NumPoints+1))/2);
  int NumSCPoints = U_ls(0).SCgrid.NumPoints;
  Array<double,1> SCBuff(NumSCPoints*NumPoints);
  if (myWorldProc == 0) {
    for (int group=1; group < NumGroups; group++)
      for (int lindex=0; lindex<ChannelsPerGroup; lindex++) {
	int l = lindex*NumGroups + group;
	int groupLeader = group * ProcsPerGroup;
	// First, receive the U matrix
	WorldComm.Receive (groupLeader, Buff);
	for (int row=0; row<NumPoints; row++)
	  for (int col=0; col<=row; col++) {
	    int index = Index(row,col);
	    U_ls(l).U(row,col) = Buff(index);
	    U_ls(l).U(col,row) = Buff(index);
	  }
	// Now receive the dU matrix
	WorldComm.Receive (groupLeader, Buff);
	for (int row=0; row<NumPoints; row++)
	  for (int col=0; col<=row; col++) {
	    int index = Index(row,col);
	    U_ls(l).dU(row,col) = Buff(index);
	    U_ls(l).dU(col,row) = Buff(index);
	  }
	// Now receive the SC matrix
	WorldComm.Receive (groupLeader, SCBuff);
	for (int row=0; row<NumSCPoints; row++)
	  for (int col=0; col<NumPoints; col++) {
	    int index = row*NumPoints+col;
	    U_ls(l).SC(row,col) = SCBuff(index);
	  }
	// Now receive the Yl array
	U_ls(l).Yl.resize(NumPoints);
	WorldComm.Receive (groupLeader, U_ls(l).Yl);
	// Make sure beta is updated
	U_ls(l).beta = U_ls(0).beta;
      }
  }
  else if (myGroupProc == 0) // I'm the group leader {
    for (int lindex=0; lindex<ChannelsPerGroup; lindex++) {
      // First, send the U matrix
      for (int row=0; row<NumPoints; row++)
	for (int col=0; col<=row; col++) {
	  //fprintf (stderr, "row = %4d col = %4d index = %4d\n", 
	  //	 row, col, Index(row,col));
	  Buff(Index(row,col)) = U_ls(lindex).U(row,col);
	}
      
      WorldComm.Send(0, Buff);
      // Next, send the dU matrix
      for (int row=0; row<NumPoints; row++)
	for (int col=0; col<=row; col++)
	  Buff(Index(row,col)) = U_ls(lindex).dU(row,col);
      WorldComm.Send(0, Buff);
      // Now send the SC matrix
      for (int row=0; row<NumSCPoints; row++)
	for (int col=0; col<NumPoints; col++) {
	  int index = row*NumPoints+col;
	  SCBuff(index) = U_ls(lindex).SC(row,col);
	}
      WorldComm.Send(0, SCBuff);
      // Now send the Yl array
      WorldComm.Send(0, U_ls(lindex).Yl);
    }
  cerr << "MyProc = " << myWorldProc << " and I'm in GatherU_ls().\n";
}

/// Broadcast all the U_l's from processor 0 to all the other processors.
/// Note that the semiclassical information is sent
void Rho::BroadcastU_ls()
{
  int myWorldProc = WorldComm.MyProc();
  int myGroupProc = GroupComm.MyProc();
  /// This will be a send buffer for proc 0 and a receive buffer
  /// for all others.
  int NumPoints = grid->NumPoints;
  Array<double,1> UBuff((NumPoints*(NumPoints+1))/2);
  int NumSCPoints = U_ls(0).SCgrid.NumPoints;
  Array<double,1> SCBuff(NumSCPoints*NumPoints);
  if (myWorldProc == 0) {
    for (int l=0; l<U_ls.size(); l++) {
      // Fill the U buffer
      for (int row=0; row<NumPoints; row++)
	for (int col=0; col<NumPoints; col++)
	  UBuff(Index(row,col)) = U_ls(l).U(row,col);
      // Now Broadcast it
      WorldComm.Broadcast(0,UBuff);

      // Fill the U buffer with dU
      for (int row=0; row<NumPoints; row++)
	for (int col=0; col<NumPoints; col++)
	  UBuff(Index(row,col)) = U_ls(l).dU(row,col);
      // Now Broadcast it
      WorldComm.Broadcast(0,UBuff);

      // Fill the SC buffer
      for (int row=0; row<NumSCPoints; row++)
	for (int col=0; col<NumPoints; col++) {
	  int index = row*NumPoints+col;
	  SCBuff(index) = U_ls(l).SC(row,col);
	}
      // Broadcast it
      WorldComm.Broadcast(0,SCBuff);
    }
  }
  else {
    // I'm receiving
    for (int l=0; l<U_ls.size(); l++) {
      WorldComm.Broadcast(0,UBuff);
      for (int row=0; row<NumPoints; row++)
	for (int col=0; col<NumPoints; col++)
	  U_ls(l).U(row,col) = UBuff(Index(row,col));
      
      WorldComm.Broadcast(0,UBuff);
      for (int row=0; row<NumPoints; row++)
	for (int col=0; col<NumPoints; col++)
	  U_ls(l).dU(row,col) = UBuff(Index(row,col));
      
      WorldComm.Broadcast(0,SCBuff);
      for (int row=0; row<NumSCPoints; row++)
	for (int col=0; col<NumPoints; col++) {
	  int index = row*NumPoints+col;
	  U_ls(l).SC(row,col) = SCBuff(index);
	}
    }
  }
}

/// Square all the l-channel actions, halving the temperature of the
/// density matrix.
void Rho::Square()
{
  if (lambda != 0.0)
    for (int lindex=0; lindex<ChannelsPerGroup; lindex++)
      {
	int l;
	if (COMM::WorldProc() == 0)
	  l = lindex*NumGroups;
	else
	  l = lindex;
	if (WriteDebug) {
	  fprintf (DebugFile, "Doing l %d, beta = %1.5e\n",
		   U_ls(l).lchannel(), beta);
	  fflush (DebugFile);
	}
	if (COMM::WorldProc() == 0)
	  cerr << "Doing l = " << U_ls(l).lchannel() << ":\n";
	U_ls(l).Square();
      }
  beta *= 2.0;
}
    

// double Rho::U(double r, double rp, double theta)
// {
//   double LogTotalRho0r = Log_rho0 (r, rp, theta, lambda, beta);
//   double x = Transform->r2x(r);
//   double xp = Transform->r2x(rp);
//   double LogTotalRho0x = Log_rho0 (x, xp, theta, lambda, beta);
//   double A, Ap, B, V, dA;
//   Pot->ABV(r, A, B, V, dA);
//   Pot->ABV(rp, Ap, B, V, dA);
//   double Prefact = 1.0/(4.0*M_PI*r*rp*pow(A*Ap,0.25));
//   double sum=0.0;
//   double costheta = cos(theta);
//   int lmax = U_ls.size()-1;
//   double Ulmax = exp (-U_ls(lmax).U(x,xp));
//   for (int l=0; l<=lmax; l++)
//     {
//       //double Exp = LogFreelDMValue (l, r, rp, beta) -
//       //	LogTotalRho0;
//       double Exp = LogFreelDMValue (l, x, xp, beta) -
// 	LogTotalRho0x;
//       double Uval = exp(-U_ls(l).U(x,xp)) ;
//       if (Extrapolate_ls)
// 	Uval -= Ulmax;
//       //Exp -= U_ls(l)(x,xp);
//       sum += (2.0*l+1.0) * Legendre(l,costheta) * Uval * exp(Exp);
//     }
//   sum *= Prefact;
//   if (Extrapolate_ls)
//     sum += Ulmax;

//   double U;
//   if (sum < 0.0)
//     U = 0.0;
//   else
//     U = log(sum) + LogTotalRho0x - LogTotalRho0;
//   U *= -1.0;

//   return (U);
// }

// double Rho::U(double r, double rp, double theta,
// 	      Array<double,1> &Ularray)
// {
//   double LogTotalRho0r = Log_rho0 (r, rp, theta, lambda, beta);
//   double x = Transform->r2x(r);
//   double xp = Transform->r2x(rp);
//   double LogTotalRho0x = Log_rho0 (x, xp, theta, lambda, beta);
//   double A, Ap, B, V, dA;
//   Pot->ABV(r, A, B, V, dA);
//   Pot->ABV(rp, Ap, B, V, dA);
//   double Prefact = 1.0/(4.0*M_PI*r*rp*pow(A*Ap,0.25));
//   double sum=0.0;
//   double costheta = cos(theta);
//   int lmax = U_ls.size()-1;
//   double Ulmax = exp (-Ularray(lmax));
//   for (int l=0; l<=lmax; l++)
//     {
//       //double Exp = LogFreelDMValue (l, r, rp, beta) -
//       //	LogTotalRho0;
//       double Exp = Log_Rho0_l (l, x, xp, beta) -
// 	LogTotalRho0x;
//       double Uval = exp(-Ularray(l));
//       if (Extrapolate_ls)
// 	Uval -= Ulmax;
//       //Exp -= U_ls(l)(x,xp);
//       sum += (2.0*l+1.0) * Legendre(l,costheta) * Uval * exp(Exp);
//     }
//   if (Extrapolate_ls)
//     sum += 4.0*M_PI*x*xp*Ulmax;

//   sum *= Prefact;

//   double U;
//   if (sum <= 0.0)
//     U = 0.0;
//   else
//     U = log(sum) + LogTotalRho0x - LogTotalRho0r;
//   U *= -1.0;

//   return (U);
// }



void Rho::U_lArray(double r, double rp, 
		   Array<double,1> &Ulvec,
		   Array<double,1> &dUlvec)
{
  Ulvec.resize(U_ls.size());
  dUlvec.resize(U_ls.size());
  double x = Transform.r2x (r);
  double xp = Transform.r2x (rp);
  for (int i=0; i<U_ls.size(); i++) {
    Ulvec(i) = U_ls(i).U(x,xp);
    dUlvec(i) = U_ls(i).dU(x,xp);
  }
}



double Rho::U(double r, double rp, double costheta,
	      Array<double,1> &Ularray)
{
  double LogTotalRho0r = Log_rho0 (r, rp, costheta, lambda, beta);
  double x = Transform.r2x(r);
  double xp = Transform.r2x(rp);
  double LogTotalRho0x = Log_rho0 (x, xp, costheta, lambda, beta);
  double A, Ap;
  A = Pot->A(r); Ap = Pot->A(rp); 
  double Prefact = 1.0/(4.0*M_PI*r*rp*pow(A*Ap,0.25));
  double sum=0.0;
  int lmax = U_ls.size()-1;
  double Ulmax = exp (-Ularray(lmax));
  for (int l=0; l<=lmax; l++)
    {
      //double Exp = LogFreelDMValue (l, r, rp, beta) -
      //	LogTotalRho0;
      double Exp = Log_rho0_l (l, x, xp, lambda, beta) -
	LogTotalRho0x;
      double Uval = exp(-Ularray(l));
      if (Extrapolate_ls)
	Uval -= Ulmax;
      //Exp -= U_ls(l)(x,xp);
      sum += (2.0*l+1.0) * Legendre(l,costheta) * Uval * exp(Exp);
    }
  if (Extrapolate_ls)
    sum += 4.0*M_PI*x*xp*Ulmax;

  sum *= Prefact;

  double U;
  //  if (sum <= 0.0)
  //  U = 0.0;
  //  else
    U = log(sum) + LogTotalRho0x - LogTotalRho0r;
  U *= -1.0;

  return (U);
}



void Rho::UdU_local(double r, double rp, double costheta,
		    Array<double,1> &Ulvec, Array<double,1> &dUlvec,
		    double &U, double &dU)
{
  double rho0r = rho0 (r, rp, costheta, lambda, beta);
  double drho0r = drho0_dbeta (r, rp, costheta, lambda, beta);
  double x = Transform.r2x(r);
  double xp = Transform.r2x(rp);
  double rho0x = rho0 (x, xp, costheta, lambda, beta);
  double drho0x = drho0_dbeta (x, xp, costheta, lambda, beta);
  double A, Ap;
  A = Pot->A(r);
  Ap = Pot->A(rp);
  double Prefact = 1.0/(4.0*M_PI*r*rp*pow(A*Ap,0.25));
  double rho=0.0;
  double drho=0.0;
  int lmax = U_ls.size()-1;
  double expmUlmax = exp (-Ulvec(lmax));
  double dUlmax = dUlvec(lmax);
  for (int l=0; l<=lmax; l++)
    {
      double rho0l = exp(Log_rho0_l (l, x, xp, lambda, beta));
      double drho0l = drho0_dbeta_l(l, x, xp, lambda, beta);
      double expmUval = exp(-Ulvec(l));
      double dUexp = dUlvec(l)*exp(-Ulvec(l));
      if (Extrapolate_ls /*&& ((r+rp)>5.0)*/) {
	expmUval -= expmUlmax;
	dUexp -= dUlmax * exp(-Ulvec(lmax));
      }
      double LegTerm = (2.0*l+1.0)* Legendre(l, costheta);
      rho += LegTerm * expmUval * rho0l;
      drho += LegTerm * (expmUval * drho0l - rho0l*dUexp);
//       fprintf (stderr, "rho0l = %1.6e\n", rho0l);
//       fprintf (stderr, "expmUval = %1.6e\n", expmUval);
//       fprintf (stderr, "Legendre(%d) = %1.6e\n", l, LegTerm);
//       fprintf (stderr, "U_%d = %1.6e\n", l, Ulvec(l));
    }
  if (Extrapolate_ls /*&& ((r+rp)>5.0)*/) {
    rho += 4.0*M_PI*x*xp*expmUlmax*rho0x;
    drho += 4.0*M_PI*x*xp*(drho0x*expmUlmax - rho0x*expmUlmax*dUlmax);
  }

  rho *= Prefact;
  drho*= Prefact;

//   fprintf (stderr, "rho0r = %1.6e\n", rho0r);
//   fprintf (stderr, "rho   = %1.6e\n", rho);
  

  U = -log(rho/rho0r);
  dU = (drho0r - exp(U)*drho)/rho0r;
}




void Rho::UdU_PH(double r, double rp, double costheta,
		 Array<double,1> &Ulvec, Array<double,1> &dUlvec,
		 double &U, double &dU)
{
  double rho0r = rho0 (r, rp, costheta, lambda, beta);
  double drho0r = drho0_dbeta (r, rp, costheta, lambda, beta);
  double x = Transform.r2x(r);
  double xp = Transform.r2x(rp);
  double rho0x = rho0 (x, xp, costheta, lambda, beta);
  double drho0x = drho0_dbeta (x, xp, costheta, lambda, beta);
  double A, Ap;
  A = Pot->A(r); Ap = Pot->A(rp);
  double Prefact = 1.0/(4.0*M_PI*r*rp*pow(A*Ap,0.25));
  double rho=0.0;
  double drho=0.0;
  int lmax = U_ls.size()-1;
  double Ulmax = Ulvec(lmax);
  double dUlmax = dUlvec(lmax);
  for (int l=0; l<=3*lmax; l++) {
    double rho0l = exp(Log_rho0_l (l, x, xp, lambda, beta));
    double drho0l = drho0_dbeta_l(l, x, xp, lambda, beta);
    double Ul, dUl;
    if ((r+rp) < 3.8) {
      Ul  = (l > lmax) ? (Ulmax*(double)l/(double)lmax) : Ulvec(l);
      dUl = (l > lmax) ? (dUlmax*(double)l/(double)lmax) : dUlvec(l);
    }
    else {
      Ul  = (l > lmax) ?  Ulmax :  Ulvec(l);
      dUl = (l > lmax) ? dUlmax : dUlvec(l);
    }
    double expmUval = exp(-Ul);
    double dUexp = dUl*exp(-Ul);
    double LegTerm = (2.0*l+1.0)* Legendre(l, costheta);
    rho += LegTerm * expmUval * rho0l;
    drho += LegTerm * (expmUval * drho0l - rho0l*dUexp);
  }

  rho *= Prefact;
  drho*= Prefact;
//   if (rho < 0.0) {
//     cerr << "rho = " << rho << endl;
//     cerr << "rho0r = " << rho0r << endl;
//     cerr << "r = " << r << endl;
//     cerr << "rp = " << rp << endl;
//     cerr << "costheta = " << costheta << endl;
//   }

//   fprintf (stderr, "rho0r = %1.6e\n", rho0r);
//   fprintf (stderr, "rho   = %1.6e\n", rho);
  

  U = -log(rho/rho0r);
  dU = (drho0r - exp(U)*drho)/rho0r;
}


void Rho::UdU(double r, double rp, double costheta,
	      Array<double,1> &Ulvec, Array<double,1> &dUlvec,
	      double &U, double &dU)
{
  // HACK HACK HACK
  // Inside the core, it seems that the approximation that 
  // U_l \approx l/lmax *U_{lmax} for l>lmax. 
  //  Outside the core, it seems that U_l \approx U_{lmax} for l>lmax.
  if (Pot->IsPH()&&((r+rp)<3.8)) {
    cerr << "Using PH version of UdU.\n";
    UdU_PH(r, rp, costheta, Ulvec, dUlvec, U, dU);
  }
  else
    UdU_local (r, rp, costheta, Ulvec, dUlvec, U, dU);
}


void Rho::UdU(double r, double rp, double costheta,
	      double &U, double &dU)
{
  Array<double,1> Ulvec, dUlvec;
  U_lArray(r, rp, Ulvec, dUlvec);
  UdU(r, rp, costheta, Ulvec, dUlvec, U, dU);
//   int numls = U_ls.size();
//   Array<double,1> Ulvec, dUlvec;
//   U_lArray(r, rp, Ulvec, dUlvec);

//   double rho0r = rho0 (r, rp, costheta, lambda, beta);
//   double drho0r = drho0_dbeta (r, rp, costheta, lambda, beta);
//   double x = Transform.r2x(r);
//   double xp = Transform.r2x(rp);

//   double rho0x = rho0 (x, xp, costheta, lambda, beta);
//   double drho0x = drho0_dbeta (x, xp, costheta, lambda, beta);
//   double A, Ap;
//   A = Pot->A(r); Ap = Pot->A(rp);
//   double Prefact = 1.0/(4.0*M_PI*r*rp*pow(A*Ap,0.25));
//   double rho=0.0;
//   double drho=0.0;
//   int lmax = U_ls.size()-1;
//   double expmUlmax = exp (-Ulvec(lmax));
//   double dUlmax = dUlvec(lmax);
//   for (int l=0; l<=lmax; l++)
//     {
//       double rho0l = exp(Log_rho0_l (l, x, xp, lambda, beta));
//       double drho0l = drho0_dbeta_l(l, x, xp, lambda, beta);
//       double expmUval = exp(-Ulvec(l));
//       double dUexp = dUlvec(l)*exp(-Ulvec(l));
//       if (Extrapolate_ls && ((r+rp)>5.0)) {
// 	expmUval -= expmUlmax;
// 	dUexp -= dUlmax * exp(-Ulvec(lmax));
//       }
//       double LegTerm = (2.0*l+1.0)* Legendre(l, costheta);
//       rho += LegTerm * expmUval * rho0l;
//       drho += LegTerm * (expmUval * drho0l - rho0l*dUexp);
//     }
//   if (Extrapolate_ls && ((r+rp)>5.0)) {
//     rho += 4.0*M_PI*x*xp*expmUlmax*rho0x;
//     drho += 4.0*M_PI*x*xp*(drho0x*expmUlmax - rho0x*expmUlmax*dUlmax);
//   }

//   rho *= Prefact;
//   drho*= Prefact;

//   U = -log(rho/rho0r);
//   dU = (drho0r - exp(U)*drho)/rho0r;
}




// Uses special symmetry of the coulomb potential to calculate U and
// dU from the l=0 partial wave.
void Rho::UdU_Coulomb (double r, double rp, double costheta,
		       double &U, double &dU)
{
  double q = 0.5*(r+rp);
  double s2 = max(0.0, r*r + rp*rp - 2.0*r*rp*costheta);
  double s = sqrt(s2);
  s = min (2.0*q, s);
  double u = max(1.0e-10,q + 0.5*s);
  double v = max(1.0e-10,q - 0.5*s);
  
  SymmBicubicSpline &U0  = U_ls(0).Uspline;
  SymmBicubicSpline &dU0 = U_ls(0).dUspline;

  U = U0(u,v);
  dU = dU0(u,v);
  double log_rho_0 = Log_rho0 (r,rp,costheta,lambda,beta);
  //double rho_0 = rho0(r,rp,costheta, lambda, beta);
  double drho_0_dbeta_scaled = drho0_dbeta_scaled(r,rp,costheta,lambda,beta);
  //double drho_0_dbeta = drho0_dbeta(r,rp,costheta,lambda,beta);
  double log_rho_00 = Log_rho0_l(0, u, v, lambda, beta);
  //double rho_00 = exp(Log_rho0_l(0, u, v, lambda, beta));

  double log_scale = -(u-v)*(u-v)/(4.0*lambda*beta);
  double drho_00_dbeta_scaled = drho0_dbeta_l_scaled(0, u, v, lambda, beta);
  //double drho_00_dbeta = drho0_dbeta_l(0, u, v, lambda, beta);
  if (s > 1.0e-7) {  // Use off-diagonal formula
    double C1 = 1.0/(4.0*M_PI*s)*(U0.d_dx(u,v) - U0.d_dy(u,v));
    double C2 = 1.0/(4.0*M_PI*s)*(dU0.d_dx(u,v) - dU0.d_dy(u,v));
    //    double C3 = 1.0+rho_00/rho_0*C1;
    double C3 = 1.0+exp(log_rho_00-log_rho_0)*C1;
    U -= log(C3);
    //double C4 = (drho_00_dbeta/rho_0 -
    //rho_00*drho_0_dbeta/(rho_0*rho_0));
    double C4 = drho_00_dbeta_scaled*exp(log_scale-log_rho_0) -
      drho_0_dbeta_scaled*exp(log_rho_00-log_rho_0);
    //double C5 = C4*C1 + rho_00/rho_0*C2;
    double C5 = C4*C1 + exp(log_rho_00-log_rho_0)*C2;
    dU -= C5/C3;
    if (isnan(C4)) {
      fprintf (stderr, "u = %1.16e v = %1.16e\n", u, v);
      fprintf (stderr, 
	       "C1=%1.16e\nC2=%1.16e\nC3=%1.16e\nC4=%1.16e\nC5=%1.16e\n", 
	       C1, C2, C3, C4, C5);
      fprintf (stderr, "log_rho0 = %1.16e\n", log_rho_0);
      fprintf (stderr, "drho00_scaled = %1.16e\n", drho_00_dbeta_scaled);
      fprintf (stderr, "drho0_scaled = %1.16e\n", drho_0_dbeta_scaled);
      fprintf (stderr, "costheta = %1.16e\n", costheta);
    }
  }
  else {             // Use diagonal formula
    double C1 = 1.0/(8.0*M_PI)*(U0.d2_dx2(u,v) + U0.d2_dy2(u,v) - 
				2.0*U0.d2_dxdy(u,v));
    double C2 = 1.0/(8.0*M_PI)*(dU0.d2_dx2(u,v) + dU0.d2_dy2(u,v) -
				2.0*dU0.d2_dxdy(u,v));
    //double C3 = 1.0+rho_00/rho_0*C1;
    double C3 = 1.0+exp(log_rho_00-log_rho_0)*C1;
    U -= log(C3);
    //double C4 = (drho_00_dbeta/rho_0 - rho_00*drho_0_dbeta/(rho_0*rho_0));
    double C4 = (drho_00_dbeta_scaled*exp(log_scale-log_rho_0) - 
		 drho_0_dbeta_scaled*exp(log_rho_00-log_rho_0));
    //double C5 = C4*C1 + rho_00/rho_0*C2;
    double C5 = C4*C1 + exp(log_rho_00-log_rho_0)*C2;
    dU -= C5/C3;
    // U -= log(1.0 + rho_00/rho_0*(8.0*M_PI)*(U0.d2_dx2(u,v) + 
    //     					    U0.d2_dy2(u,v) -
    //     					    2.0*U0.d2_dxdy(u,v)));
  }
}
  



double Rho::U(double r, double rp, double costheta)
{
//   cerr << "r = " << r << " rp = " << rp << " costheta = " 
//        << costheta << endl;

  double LogTotalRho0r = Log_rho0 (r, rp, costheta, lambda, beta);
  double x = Transform.r2x(r);
  double xp = Transform.r2x(rp);
  double LogTotalRho0x = Log_rho0 (x, xp, costheta, lambda, beta);
  double A, Ap;
  A = Pot->A(r);  Ap = Pot->A(rp);
  double Prefact = 1.0/(4.0*M_PI*r*rp*pow(A*Ap,0.25));
  double sum=0.0;
  int lmax = U_ls.size()-1;

  double Ulmax = exp (-U_ls(lmax).U(r,rp));
  for (int l=0; l<=lmax; l++)
    {
      double Exp = Log_rho0_l (l, x, xp, lambda, beta) -
	LogTotalRho0x;
      double Uval = exp(-U_ls(l).U(x,xp));
      if (Extrapolate_ls)
	Uval -= Ulmax;
      sum += (2.0*l+1.0) * Legendre(l,costheta) * Uval * exp(Exp);
    }
  if (Extrapolate_ls)
    sum += 4.0*M_PI*x*xp*Ulmax;

  sum *= Prefact;

  double U;
  //  if (sum <= 0.0)
  //  U = 0.0;
  //  else
  U = log(sum) + LogTotalRho0x - LogTotalRho0r;
  U *= -1.0;

  return (U);
}


			    

void Rho::Write(IOSectionClass &outSection)
{
  outSection.WriteVar ("version", "1.0");
  outSection.WriteVar ("beta", beta);
  outSection.WriteVar ("lambda", lambda);
    
  outSection.NewSection ("Potential");
  Pot->Write (outSection);
  outSection.CloseSection();

  outSection.NewSection ("U_ls");
  outSection.NewSection ("Grid");
  grid->Write(outSection);
  outSection.CloseSection();  // "Grid"
  
  for (int l=0; l<U_ls.size(); l++)
    {
      outSection.NewSection("U_l");
      U_ls(l).Write(outSection);
      outSection.CloseSection (); // "U_l"
    }
    outSection.CloseSection (); // "U_ls"
}




void Rho::WriteU_ls(IOSectionClass &outSection)
{
  outSection.NewSection ("U_ls");
  outSection.WriteVar ("beta", beta);
  for (int l=0; l<U_ls.size(); l++)
    {
      outSection.NewSection("U_l");
      U_ls(l).Write(outSection);
      outSection.CloseSection (); // "U_l"
    }
  outSection.CloseSection (); // "U_ls"
  outSection.FlushFile();
}


void Rho::ReadU_ls(IOSectionClass &inSection)
{
  inSection.ReadVar ("beta", beta);
  for (int l=0; l<U_ls.size(); l++)
    {
      assert(inSection.OpenSection("U_l",l));
      U_ls(l).Read(l, lambda, beta, grid, &Transform,
		   Pot, WorldComm, inSection);
      inSection.CloseSection (); // "U_l"
    }
}


void Rho::Read(IOSectionClass &inSection)
{
  string version;
  assert(inSection.ReadVar("version", version));
  assert(version == "1.0");
  assert(inSection.ReadVar("beta", beta));
  assert(inSection.ReadVar("lambda", lambda));
  assert(inSection.OpenSection("Potential"));
  Pot=ReadPotential(inSection);
  inSection.CloseSection();
  Transform.Initialize(Pot, 1000);
  cerr << "x(r=2.0) = " << Transform.r2x(2.0) << endl;
  cerr << "r(x=2.0) = " << Transform.x2r(2.0) << endl;
  assert(inSection.OpenSection("U_ls"));
  assert(inSection.OpenSection("Grid"));
  grid = ReadGrid (inSection);
  inSection.CloseSection(); // "Grid"
  int NumU_ls = inSection.CountSections ("U_l");
  U_ls.resize(NumU_ls);
  for (int l=0; l<U_ls.size(); l++)
    {
      assert(inSection.OpenSection("U_l", l));
      U_ls(l).Read (l, lambda, beta, grid, &Transform, Pot,
		    GroupComm, inSection);
      inSection.CloseSection(); // "U_l"
    }
  inSection.CloseSection(); // "U_ls"
}



