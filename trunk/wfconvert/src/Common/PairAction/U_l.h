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

#ifndef U_L_H
#define U_L_H

#include "../PH/CoreTransform.h"
#include "../Splines/BicubicSpline.h"
#include "../MPI/Communication.h"
#include "../Integration/GKIntegration.h"
#include "../SpecialFunctions/SpecialFunctions.h"
#include "../IO/IO.h"
#include "FreeParticles.h"

////////////////////////////////////////////////////////////
//                    The U_l class                       //
// The U_l class holds the potential action part of each  //
// partial wave.  l refers to the angular momentum        //
// channel.  It also stores the inverse temperature to    //
// which it corresponds and the grid on which it is       //
// represented.                                           // 
////////////////////////////////////////////////////////////

class U_l
{
private:
  const double MaxRatio;
  /// Stores the semiclassical approximation to the U action outside
  /// the U grid;
  double xMaxSC;
  double xMax;
  CoreTransform *transform;
  Potential *Pot;

  int l, NumPoints;
  CommunicatorClass GroupComm;
  // Use the expm1 and log1p functions to preserve accuracy at high
  // temperature.
public:
  SymmBicubicSpline Uspline;
  SymmBicubicSpline dUspline;
  blitz::Array<double,1> Yl;
  BicubicSpline SCspline;

  double beta, lambda, TailPower, AbsTol, RelTol;
  Grid *grid;
  LinearGrid SCgrid;
  bool UseSC, Use_m1;
      
  inline double W_l(double x) 
  {
    double r = transform->x2r(x);
    double A, B, V, dAdr, d2Adr2;
    A=Pot->A(r); 
    B=Pot->B(r); 
    V=Pot->V(r); 
    dAdr=Pot->dAdr(r); 
    d2Adr2=Pot->d2Adr2(r);
    double sl = (double)l;
    double Wl = sl*(sl+1.0)*B / (2.0*r*r);
    Wl -= dAdr*dAdr / (32.0*A);
    Wl += 0.125 * d2Adr2 + 0.5 *dAdr/r;
    Wl += V;
    // HACK to cut out coulomb singularity
    if (beta * Wl < -200.0)
      Wl = -200.0/beta;
    return (Wl);
  }

  inline double Y_l(double x) {
    double Wl = W_l(x);
    double A0, B0;
    A0 = Pot->A(0.0); 
    B0 = Pot->B(0.0);
    double llplus1 = l*(l+1);
    
    return (Wl - B0/A0*0.5*llplus1/(x*x));
  }


  inline double USC(double x, double xp)
  {
    if (xp<x)
      {
	double temp = x;
	x = xp;
	xp = temp;
      }

    double delta = (xp-x);
    if (delta == 0.0)
      return (beta*Y_l(x));
    else
      {
	GKIntegration<U_l,GK31> Yl(*this);
	//Yl.SetRelativeErrorMode();
	return (beta*Yl.Integrate(x,xp,1.0e-6,1.0e-6,false) / delta);
      }
  }

  inline double dUSC (double x, double xp)
  {
    return (USC(x,xp)/beta);
  }


  inline double U(int row, double x) 
  {
    x = max (x, grid->Start);
    double UVal;

    if (UseSC)
      UVal = USC((*grid)(row), x);
    else if (x<=xMax)
      UVal = (Uspline(row,x));
    else {
      double Umax = Uspline(row, NumPoints-1);
      double U2 = Uspline (row, 0.99*xMax);
      double a = (U2-Umax)/(pow(0.99,TailPower)-1.0);
      double c = Umax - a;
      //UVal = Umax*xMax/x;
      UVal = c + a *pow(x/xMax,TailPower);

//       double U_SCmax = SCspline(0,row);
//       double U_SC;
//       if (x<=xMaxSC)
// 	U_SC = SCspline(x,row);
//       //	  U_SC = SCspline(row,x);
//       else {
// 	//HACK
// 	cerr << "Warning:  Outside SC table in dU.  Very slow.\n";
// 	U_SC = USC((*grid)(row),x);
//       }


//       //double ratio = U_SC/U_SCmax;      
//       // Check for singularities   
// //       if ((ratio > MaxRatio) || (U_SCmax == 0.0))
// // 	ratio = MaxRatio;
// //       if (ratio < (-MaxRatio))
// // 	ratio = -MaxRatio;
// //       UVal = ratio * Umax;
//       UVal = Umax + SCscale(row)*(U_SC-U_SCmax);
      

// //       if ((U_SCmax == 0.0) || fabs(ratio)>50.0)
// // 	UVal = Umax + beta*(U_SC-U_SCmax);
// //       else
// // 	UVal = U_SC/U_SCmax * Umax;
    }
    return (UVal);
  }

  inline double dU(int row, double x) 
  {
    x = max (x, grid->Start);
    double dUVal;

    if (UseSC)
      dUVal = dUSC((*grid)(row), x);
    else if (x<=xMax)
      dUVal = (dUspline(row,x));
    else {
      double dUmax = dUspline(row, NumPoints-1);
      //dUVal = dUmax*xMax/x;
      double dU2 = dUspline (row, 0.99*xMax);
      double a = (dU2-dUmax)/(pow(0.99,TailPower)-1.0);
      double c = dUmax - a;
      dUVal = c + a *pow(x/xMax,TailPower);


//       double U_SCmax = SCspline(0,row);
//       double U_SC;
//       if (x<=xMaxSC)
// 	//	  U_SC = SCspline(row,x);
// 	U_SC = SCspline(x,row);
//       else {
// 	// HACK
// 	cerr << "Warning:  Outside SC table in dU.  Very slow.\n";
// 	U_SC = USC((*grid)(row),x);
//       }

// //       double ratio = U_SC/U_SCmax;      
// //       // Check for singularities   
// //       if ((ratio > MaxRatio) || (U_SCmax == 0.0))
// // 	ratio = MaxRatio;
// //       if (ratio < (-MaxRatio))
// // 	ratio = -MaxRatio;
// //       dUVal = ratio * dUmax;

//       dUVal = dUmax + SCscale(row)/beta*(U_SC - U_SCmax);
      
//       //       if ((U_SCmax == 0.0) || fabs(ratio)>50.0)
//       // 	dUVal = dUmax;
//       //       else
//       // 	dUVal = U_SC/U_SCmax * dUmax;
    }
    return (dUVal);
  }

  inline int lchannel()
  {
    return (l);
  }

  // U_l is symmetric in x, x';
  inline double U(double x, int col) 
  { return (U(col, x)); }
  inline double dU(double x, int col) 
  { return (dU(col, x)); }
  inline double U(int i, int j) const
  { return (Uspline(i,j)); }
  inline double dU(int i, int j) const
  { return (dUspline(i,j)); }
  inline double & U(int i, int j) 
  { return (Uspline(i,j)); }
  inline double & dU(int i, int j) 
  { return (dUspline(i,j)); }
  inline double SC (int i, int j) const
  { return SCspline(i,j); }
  inline double& SC (int i, int j)
  { return SCspline(i,j); }

  inline double U(double x, double xp)
  {
    //      cerr << "x = " << x << " xp = " << xp 
    //   << "l = " << l << endl;
    if ((x > xMax) || (xp > xMax) || UseSC) {
      //cerr << "Using SC!!!\n";
      return (USC(x,xp));
    }
    else
      return (Uspline(x, xp));
  }
  
  inline double dU(double x, double xp)
  {
    if ((x > xMax) || (xp > xMax) || UseSC)
      return (dUSC(x,xp));
    else
      return (dUspline(x, xp));
  }
  
  /// This operator is used in integration of Y_l(x) for
  /// the semiclassical approximation
  inline double operator()(double x)
  {
    return (Y_l(x));
  }

  void Initialize(int l_, double lambda_, double FinalBeta, int NumSquares, 
		  Grid *grid_, CoreTransform *transform_, 
		  Potential *Pot_, CommunicatorClass groupComm);

  void FillWithSC();
  
  void Read (char *FileName);
  // Square correctly reduces the temperature of the partial wave
  // by a factor of 2 (increases beta by a factor or 2).
  void Square();
  void CalcElement(int row, int col, double &U, double &dU);
  void CalcElement_Hermite(int row, int col, double &U, double &dU);

  double rho_l(int row, int col);
  double rho_l(double x, double xp);
  double H_rho_l(int row, int col);
  double H_rho_l_FD(int row, int col);
  double dU_Bloche(int row, int col);
  double drho_l_dbeta(int row, int col);
  double drho_l_dbeta(double x, double xp);
  void Write (IOSectionClass &outSection);
  bool Read (int this_l, double this_lambda, 
	     double this_beta, Grid *Ugrid, CoreTransform *this_transform, 
	     Potential *this_Pot, CommunicatorClass comm,
	     IOSectionClass &inSection);
  U_l() : 
    MaxRatio (50.0), 
    UseSC(false)
  { /* do nothing */ }
};





class USquaringIntegrand
{
public:
  int row, col;
  U_l *Ul;
  double ScaleExp;
  int Ncalls;

  inline double Exp (double delta)
  {
    double beta = Ul->beta;
    double lambda = Ul->lambda;
    int l = Ul->lchannel();
    double x = (*(Ul->grid))(row);
    double xp = (*(Ul->grid))(col);
    double xbar = 0.5*(x+xp);
    double xpp = xbar + delta;
    double Log = LogFreeIntegrand(l, x,xp,delta, lambda, beta);
    double Ul1 = Ul->U(row,xpp);
    double Ul2 = Ul->U(col,xpp);
    Log -= (Ul1 + Ul2);
    return (Log);
  }

  inline void SetScaleExp()
  {
//     const int NumSamples = 200;
//     double x = (*(Ul->grid))(row);
//     double xp = (*(Ul->grid))(col);
//     double center = 0.5*(x+xp);
//     double sigma = sqrt(0.5*Ul->beta);
//     double low =  -9.0*sigma;
//     double high = 9.0*sigma;
//     if ((low+center) < 1.0e-10)
//       low = 1.0e-10 - center;
//     //    double Max = Exp(1.0e-10-center);
//     // if (Exp(low) > Max)
//     double Max = Exp(low);
//     double deltaMax = low;
//     int imax = 0;
//     for (int i=0; i<NumSamples; i++)
//       {
// 	double delta = low + (high-low)*i/(NumSamples-1);
// 	double val = Exp(delta);
// 	if (val > Max)
// 	  {
// 	    Max = val;
// 	    deltaMax = delta;
// 	    imax = i;
// 	  }
//       }
//     // Now, zero in more with binary search
//     double a = low + (high-low)*(imax-1)/(NumSamples-1);
//     if ((a+center) < 1.0e-10)
//       a = 1.0e-10 - center;
//     double b = low + (high-low)*(imax+1)/(NumSamples-1);
//     if ((b+center) < 1.0e-10)
//       b = 1.0e-10 - center;
//     double mid, midb, mida;
//     int done = 0;
//     while (!done)
//       {
// 	mid = 0.5 *(a+b);
// 	mida = 0.5 *(a+mid);
// 	midb = 0.5 *(b+mid);
	
// 	double valmid = Exp(mid);
// 	double vala = Exp(mida);
// 	double valb = Exp(midb);
   
// 	if (vala > valb)
// 	  b = midb;
// 	else
// 	  a = mida;
// 	if ((b-a)/sigma < 1.0e-10)
// 	  done = 1;
// 	if (valmid > Max)
// 	  Max = valmid;
//       }
    

//     ScaleExp = -Max;
    // HACK
    ScaleExp = 0.0;
  }

  inline double operator()(double delta)
  {
    double beta = Ul->beta;
    double lambda = Ul->lambda;
    int l = Ul->lchannel();
    double x = (*(Ul->grid))(row);
    double xp = (*(Ul->grid))(col);
    double xbar = 0.5*(x+xp);
    double xpp = xbar + delta;
    double LogI = LogFreeIntegrand(l, x,xp,delta, lambda, beta);
    double Ul1 = Ul->U(row,xpp);
    double Ul2 = Ul->U(col,xpp);
    if (Ul->Use_m1)
      return (exp(LogI)*expm1(-(Ul1+Ul2)));
    else
      return (exp(LogI-(Ul1+Ul2)));
  }


//   inline double operator()(double delta)
//   {
//     double Log = Exp(delta);
//     Log += ScaleExp;
    
//     Ncalls++;

//     if (isnan(exp(Log)))
//       {
// 	cerr << "We have an NAN! row = " << row << " col = " << col << "\n";
// 	cerr << "Log = " << Log << "\n";
//       }
//     return (exp(Log));
//   }
  USquaringIntegrand()
  {
    Ncalls = 0;
  }
};




class USquaringIntegrand_Hermite
{
public:
  int row, col;
  U_l *Ul;
  double ScaleExp;
  int Ncalls;

  inline double Exp (double delta)
  {
    double beta = Ul->beta;
    double lambda = Ul->lambda;
    int l = Ul->lchannel();
    double x = (*(Ul->grid))(row);
    double xp = (*(Ul->grid))(col);
    double xbar = 0.5*(x+xp);
    double xpp = xbar + delta;
    double Log = LogFreeIntegrand_Hermite(l, x,xp,delta, lambda, beta);
    double Ul1 = Ul->U(row,xpp);
    double Ul2 = Ul->U(col,xpp);
    Log -= (Ul1 + Ul2);
    return (Log);
  }

  inline void SetScaleExp()
  {
    ScaleExp = 0.0;
  }


  inline double operator()(double delta)
  {
    double beta = Ul->beta;
    double lambda = Ul->lambda;
    int l = Ul->lchannel();
    double x = (*(Ul->grid))(row);
    double xp = (*(Ul->grid))(col);
    double xbar = 0.5*(x+xp);
    double xpp = xbar + delta;
    if (xpp < 0.0)
      return (0.0);
    double LogI = LogFreeIntegrand_Hermite(l, x,xp,delta, lambda, beta);
    double Ul1 = Ul->U(row,xpp);
    double Ul2 = Ul->U(col,xpp);
    if (Ul->Use_m1)
      return (exp(LogI)*expm1(-(Ul1+Ul2)));
    else
      return (exp(LogI-(Ul1+Ul2)));
  }


//   inline double operator()(double delta)
//   {
//     double x = (*(Ul->grid))(row);
//     double xp = (*(Ul->grid))(col);
//     double xbar = 0.5*(x+xp);
//     double xpp = xbar + delta;
//     if (xpp < 0.0)
//       return (0.0);

//     double Log = Exp(delta);
//     Log += ScaleExp;
    
//     Ncalls++;

//     if (isnan(exp(Log)))
//       {
// 	cerr << "We have an NAN! row = " << row << " col = " << col << "\n";
// 	cerr << "Log = " << Log << "\n";
//       }
//     return (exp(Log));
//   }
  USquaringIntegrand_Hermite()
  {
    Ncalls = 0;
  }
};





class dUSquaringIntegrand
{
public:
  int row, col;
  U_l *Ul;
  double ScaleExp;
  double sign;

  int Ncalls;
  
  inline double Exp (double delta)
  {
    double beta = Ul->beta;
    double lambda = Ul->lambda;
    int l = Ul->lchannel();
    double x = (*(Ul->grid))(row);
    double xp = (*(Ul->grid))(col);
    double xbar = 0.5*(x+xp);
    double xpp = xbar + delta;
    double Log = LogFreeIntegrand(l, x,xp,delta,lambda,beta);
    double U1  = Ul->U (row,xpp);
    double U2  = Ul->U (col,xpp);
    double dU1 = Ul->dU(row,xpp);
    double dU2 = Ul->dU(col,xpp);
    Log -= (U1 + U2);

    double lbetainv = 1.0/(lambda*beta);
    double y = 0.5*x*xpp*lbetainv;
    double yp = 0.5*xp*xpp*lbetainv;
    double z = 0.25*x*xp*lbetainv;
    double Cpart;
//     if (z < 20.0)
//       Cpart = betainv *(0.25*betainv*(x*x+xp*xp+4.0*xpp*xpp) 
// 			-(1.5 + (double)l)
// 			- y*il_scaled(l+1,y)/il_scaled(l,y)
// 			- yp*il_scaled(l+1,yp)/il_scaled(l,yp)
// 			+ z*il_scaled(l+1,z)/il_scaled(l,z));
    if ((((double)l/z) > 0.005)||(z<20.0))//(z < 5000.0)
      Cpart = 0.5*lbetainv*(0.125*lbetainv*(x*x+xp*xp+4.0*xpp*xpp) 
			-(1.5 + (double)l)
			- y*exp(log_il_scaled(l+1,y)-log_il_scaled(l,y))
			- yp*exp(log_il_scaled(l+1,yp)-log_il_scaled(l,yp))
			+ z*exp(log_il_scaled(l+1,z)-log_il_scaled(l,z)));
    else
      {
	double sl = (double)l;
	double alpha = 0.5 * sl*(sl+1.0);
	double gamma = -0.125*(sl-2.0)*sl*(sl+1.0)*(sl+3.0);
	double epsilon = -0.5*sl*(sl+1.0)*(sl*sl+sl-3.0);
	Cpart = 0.5*lbetainv*(delta*delta*0.5*lbetainv - 0.5);
	Cpart += 2.0*alpha *delta/(x*xp*xpp);
	Cpart += 8.0*alpha*lambda*beta/(x*x*xp*xp*xpp*xpp)*
	  (delta*(xpp+xbar) + 0.5*x*xp);
	//Cpart += betainv*alpha*(-1.0/(y*y) - 1.0/(yp*yp) + 1.0/(z*z));
	Cpart += 0.5*lbetainv*gamma*(-1.0/(y*y*y) - 1.0/(yp*yp*yp) + 1.0/(z*z*z));
	Cpart += 0.5*lbetainv*epsilon*(-1.0/(y*y*y*y) - 1.0/(yp*yp*yp*yp) 
				  + 1.0/(z*z*z*z));
	//Cpart += alpha/beta*(-1.0/y -1.0/yp + 1.0/z);
	//Cpart -= 2.0*gamma*beta/(x*xp*xpp*xpp);
	//Cpart += gamma/beta * (-(1.0/(y*y) + 1.0/(yp*yp)) + 1.0/(z*z));
      }
    /*    double Cpart = betainv*(xpp*xpp - (y* il_scaled(l+1,y) /il_scaled(l,y) +
	  yp*il_scaled(l+1,yp)/il_scaled(l,yp)));*/
    double dUpart = -(dU1+dU2);
    //HACK
    //    dUpart = 0.0;
    double PartB = 0.5*(Cpart + dUpart);
    sign = 2.0*((double)(PartB>0.0)) - 1.0;
    if (PartB == 0.0)
      return (-1.0e90);
    Log += log (fabs(PartB));
    return Log;
  }


  
  inline void SetScaleExp()
  {
//     const int NumSamples = 200;
//     double x = (*(Ul->grid))(row);
//     double xp = (*(Ul->grid))(col);
//     double center = 0.5*(x+xp);
//     double sigma = sqrt(0.5*Ul->beta);
//     double low =  -9.0*sigma;
//     double high = 9.0*sigma;
//     if ((low+center) < 1.0e-10)
//       low = 1.0e-10 - center;
//     //    double Max = Exp(1.0e-10-center);
//     // if (Exp(low) > Max)
//     double Max = Exp(low);
//     double deltaMax = low;
//     int imax = 0;
//     for (int i=0; i<NumSamples; i++)
//       {
// 	double delta = low + (high-low)*i/(NumSamples-1);
// 	double val = Exp(delta);
// 	if (val > Max)
// 	  {
// 	    Max = val;
// 	    deltaMax = delta;
// 	    imax = i;
// 	  }
//       }
//     // Now, zero in more with binary search
//     double a = low + (high-low)*(imax-1)/(NumSamples-1);
//     if ((a+center) < 1.0e-10)
//       a = 1.0e-10 - center;
//     double b = low + (high-low)*(imax+1)/(NumSamples-1);
//     if ((b+center) < 1.0e-10)
//       b = 1.0e-10 - center;
//     double mid, midb, mida;
//     int done = 0;
//     while (!done)
//       {
// 	mid = 0.5 *(a+b);
// 	mida = 0.5 *(a+mid);
// 	midb = 0.5 *(b+mid);
	
// 	double valmid = Exp(mid);
// 	double vala = Exp(mida);
// 	double valb = Exp(midb);
   
// 	if (vala > valb)
// 	  b = midb;
// 	else
// 	  a = mida;
// 	if ((b-a)/sigma < 1.0e-10)
// 	  done = 1;
// 	if (valmid > Max)
// 	  Max = valmid;
//       }
    

//     ScaleExp = -Max;
    // HACK
    ScaleExp = 0.0;
  }

  inline double operator()(double delta)
  {
    double expn = Exp(delta);
    Ncalls++;
    return (sign*exp(expn+ScaleExp));
  }
  dUSquaringIntegrand()
  {
    Ncalls = 0;
  }
};




class dUSquaringIntegrand_Hermite
{
public:
  int row, col;
  U_l *Ul;
  double ScaleExp;
  double sign;

  int Ncalls;
  
  inline double Exp (double delta)
  {
    double beta = Ul->beta;
    double lambda = Ul->lambda;
    int l = Ul->lchannel();
    double x = (*(Ul->grid))(row);
    double xp = (*(Ul->grid))(col);
    double xbar = 0.5*(x+xp);
    double xpp = xbar + delta;
    double Log = LogFreeIntegrand_Hermite(l, x,xp,delta,lambda,beta);
    double U1  = Ul->U (row,xpp);
    double U2  = Ul->U (col,xpp);
    double dU1 = Ul->dU(row,xpp);
    double dU2 = Ul->dU(col,xpp);
    Log -= (U1 + U2);

    double lbetainv = 1.0/(lambda*beta);
    double y = 0.5*x*xpp*lbetainv;
    double yp = 0.5*xp*xpp*lbetainv;
    double z = 0.25*x*xp*lbetainv;
    double Cpart;
//     if (z < 20.0)
//       Cpart = betainv *(0.25*betainv*(x*x+xp*xp+4.0*xpp*xpp) 
// 			-(1.5 + (double)l)
// 			- y*il_scaled(l+1,y)/il_scaled(l,y)
// 			- yp*il_scaled(l+1,yp)/il_scaled(l,yp)
// 			+ z*il_scaled(l+1,z)/il_scaled(l,z));
    if ((((double)l/z) > 0.005)||(z<20.0))//(z < 5000.0)
      Cpart = 0.5*lbetainv *(0.125*lbetainv*(x*x+xp*xp+4.0*xpp*xpp) 
			-(1.5 + (double)l)
			- y*exp(log_il_scaled(l+1,y)-log_il_scaled(l,y))
			- yp*exp(log_il_scaled(l+1,yp)-log_il_scaled(l,yp))
			+ z*exp(log_il_scaled(l+1,z)-log_il_scaled(l,z)));

    else
      {
	/// He we have used 6 terms in an aymptotic expansion for i_l
	/// to compute the free-particle part of the integrand.
	double sl = (double)l;
	double alpha = 0.5 * sl*(sl+1.0);
	double gamma = -0.125*(sl-2.0)*sl*(sl+1.0)*(sl+3.0);
	double epsilon = -0.5*sl*(sl+1.0)*(sl*sl+sl-3.0);
	Cpart = 0.5*lbetainv*(delta*delta*0.5*lbetainv - 0.5);
	Cpart += 2.0*alpha *delta/(x*xp*xpp);
	Cpart += 8.0*alpha*lambda*beta/(x*x*xp*xp*xpp*xpp)*
	  (delta*(xpp+xbar) + 0.5*x*xp);
	Cpart += 0.5*lbetainv*gamma*(-1.0/(y*y*y) - 
				     1.0/(yp*yp*yp) + 1.0/(z*z*z));
	Cpart += 0.5*lbetainv*epsilon*(-1.0/(y*y*y*y) - 1.0/(yp*yp*yp*yp) 
				       + 1.0/(z*z*z*z));
      }
    double dUpart = -(dU1+dU2);
    //HACK
    //dUpart = 0.0;
    double PartB = 0.5*(Cpart + dUpart);
    sign = 2.0*((double)(PartB>0.0)) - 1.0;
    if (PartB == 0.0)
      return (-1.0e90);
    Log += log (fabs(PartB));
    return Log;
  }


  
  inline void SetScaleExp()
  {
    ScaleExp = 0.0;
  }

  inline double operator()(double delta)
  {
    double x = (*(Ul->grid))(row);
    double xp = (*(Ul->grid))(col);
    double xbar = 0.5*(x+xp);
    double xpp = xbar + delta;
    if (xpp < 0.0)
      return (0.0);

    double expn = Exp(delta);
    Ncalls++;
    return (sign*exp(expn+ScaleExp));
  }
  dUSquaringIntegrand_Hermite()
  {
    Ncalls = 0;
  }
};





class dUSquaringIntegrandFD
{
public:
  int row, col;
  U_l *Ulplus, *Ulminus, *Ul;
  double delta_beta;
  double ScaleExp;

  inline double Val (double delta)
  {
    double beta = 0.5*(Ulplus->beta + Ulminus->beta);
    double lambda = Ul->lambda;
    beta = Ul->beta;
    //double deltabeta = Ulplus->beta - Ulminus->beta;
    double lbetainv = 1.0/(lambda*beta);
    int l = Ul->lchannel();
    double x = (*(Ul->grid))(row);
    double xp = (*(Ul->grid))(col);
    double xbar = 0.5*(x+xp);
    double xpp = xbar + delta;

    double rhoplus = Ulplus->rho_l(x,xpp);
    double rhominus = Ulminus->rho_l(x,xpp);
    double drho = (rhoplus-rhominus)/(2.0*delta_beta);
    double rho = Ul->rho_l(xp,xpp);
    double log_rho0 = Log_rho0_l(l, x,xp,lambda,2.0*beta);
    double rho0 = exp(log_rho0);

    return (rho * drho/ rho0);
  }

  inline double operator()(double delta)
  {
    return (Val(delta));
  }
};


#endif
