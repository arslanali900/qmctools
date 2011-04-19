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

#include "FreeParticles.h"
#include "../SpecialFunctions/SpecialFunctions.h"



//////////////////////////////////////////////////////////////////////
// LogFreeIntegrand:  computes the logorithm of the quantity:       //
// \frac{\rho_l^0(x,x'';\beta) \rho_l^0(x'',x';\beta)}              //
//      {\rho_l^0(x,x';2\beta),                                     //
// where delta = xpp - xbar, and xbar = 0.5*(x+xp).                 //
//////////////////////////////////////////////////////////////////////

double LogFreeIntegrand(int l, double x, double xp, double delta, 
			double lambda, double beta)
{
  double lbetainv = 1.0/(lambda*beta);
  double xbar = 0.5*(x+xp);
  double xpp = xbar + delta;
  double Log;
  double y = 0.5*x*xpp*lbetainv;
  double y2 = 2.0*y;
  double yp = 0.5*xp*xpp*lbetainv;
  double yp2 = 2.0*yp;
  double z = x*xp*0.25*lbetainv;
  double z2 = 2.0*z;
  double fl = (double)l;

  Log = log(4.0*M_PI*xpp*xpp) -1.5*log(2.0*M_PI*lambda*beta);
  Log += (log_il_scaled(l, 0.5*x*xpp*lbetainv) + 
	  log_il_scaled(l, 0.5*xp*xpp*lbetainv) -
	  log_il_scaled(l,0.25*x*xp*lbetainv));
  Log -= (delta*delta)/(2.0*lambda*beta);
  return (Log);
}




double LogFreeIntegrand_Hermite(int l, double x, double xp, double delta, 
				double lambda, double beta)
{
  double lbetainv = 1.0/(lambda*beta);
  double xbar = 0.5*(x+xp);
  double xpp = xbar + delta;
  double Log;
  double y = 0.5*x*xpp*lbetainv;
  double y2 = 2.0*y;
  double yp = 0.5*xp*xpp*lbetainv;
  double yp2 = 2.0*yp;
  double z = x*xp*0.25*lbetainv;
  double z2 = 2.0*z;
  double fl = (double)l;

  Log = log(4.0*M_PI*xpp*xpp) -1.5*log(2.0*M_PI*lambda*beta);
  Log += (log_il_scaled(l, 0.5*x*xpp*lbetainv) + 
	  log_il_scaled(l, 0.5*xp*xpp*lbetainv) -
	  log_il_scaled(l,0.25*x*xp*lbetainv));
  return (Log);
}




////////////////////////////////////////////////////////////
// This is the free-particle density matrix with the      //
// l(l+1)/(2*r^2) term in the kinetic energy operator.    //
////////////////////////////////////////////////////////////
double Log_rho0_l (int l, double x, double xp, double lambda, double beta)
{
  double logDMvalue;
  double lbetainv = 1.0/(lambda*beta);

  logDMvalue = -1.5*log(4.0*M_PI*lambda*beta) + log(4.0*M_PI*x*xp);
  logDMvalue += -0.25*lbetainv * (x-xp)*(x-xp);
  logDMvalue += log_il_scaled(l, 0.5*x*xp*lbetainv);
  return (logDMvalue);
}

  

double drho0_dx_l(int l, double x, double xp, double lambda, double beta)
{
  double logDeriv;
  double lbeta = lambda*beta;
  double lbetainv = 1.0/(lambda *beta);
  double z = 0.5*x*xp*lbetainv;
  double deriv = 4.0*M_PI *xp / sqrt(64.0*M_PI*M_PI*M_PI*lbeta*lbeta*lbeta);
  deriv *= exp(-(x-xp)*(x-xp)*0.25*lbetainv);
  deriv *= (exp(log_il_scaled(l,z))*(1.0-0.5*x*x*lbetainv) 
	    + z*dil_dz_scaled(l, z));
  return deriv;
}


double drho0_dx_l_FD(int l, double x, double xp, double lambda, double beta)
{
  double delta = 1.0e-6;
  double plus = exp(Log_rho0_l(l, x+delta, xp, lambda, beta));
  double minus = exp(Log_rho0_l(l, x-delta, xp, lambda, beta));
  return (((plus-minus)/(2.0*delta)));
}
		

double d2rho0_dx2_l(int l, double x, double xp, double lambda, double beta)
{
  double deriv2;
  double lbeta = lambda *beta;
  double lbetainv = 1.0/beta;
  double z = 0.5*x*xp*lbetainv;
  deriv2 = 2.0*M_PI*xp*lbetainv /sqrt(64.0*M_PI*M_PI*M_PI*lbeta*lbeta*lbeta);
  deriv2 *= exp(-(x-xp)*(x-xp)*0.25*lbetainv);
  deriv2 *= ((0.5*x*x*x*lbetainv - 3.0*x)*exp(log_il_scaled(l,z)) +
	     2.0*xp*(1.0-0.5*x*x*lbetainv)*dil_dz_scaled(l,z) +
	     xp*z*d2il_dz2_scaled(l,z));
  return (deriv2);
}


double d2rho0_dx2_l_FD(int l, double x, double xp, double lambda, double beta)
{
  double delta = 2.0e-4;
  double plus = exp(Log_rho0_l(l,x+delta, xp, lambda, beta));
  double zero = exp(Log_rho0_l(l,x,       xp, lambda, beta));
  double minus = exp(Log_rho0_l(l,x-delta, xp, lambda, beta));
  return ((plus - 2.0*zero + minus)/(delta*delta));
}



double drho0_dbeta_l_scaled(int l, double x, double xp, 
			    double lambda, double beta)
{
  double lbeta = lambda * beta;
  double lbetainv = 1.0/(lambda*beta);
  double z = 0.5*x*xp*lbetainv;
  double drho = 2.0*M_PI*x*xp*lbetainv/
    sqrt(64.0*M_PI*M_PI*M_PI*lbeta*lbeta*lbeta);
  //drho *= exp(-(x-xp)*(x-xp)*0.25*lbetainv);
  drho *= (((x*x+xp*xp)*0.25*lbetainv - 1.5)*exp(log_il_scaled(l, z))
	   -z*dil_dz_scaled(l,z));

  return (drho);
}



     

double drho0_dbeta_l(int l, double x, double xp, double lambda, double beta)
{
  double lbeta = lambda * beta;
  double lbetainv = 1.0/(lambda*beta);
  double z = 0.5*x*xp*lbetainv;
  double drho = 2.0*M_PI*x*xp*lbetainv/
    sqrt(64.0*M_PI*M_PI*M_PI*lbeta*lbeta*lbeta);
  drho *= exp(-(x-xp)*(x-xp)*0.25*lbetainv);
  drho *= (((x*x+xp*xp)*0.25*lbetainv - 1.5)*exp(log_il_scaled(l, z))
	   -z*dil_dz_scaled(l,z));

  return (drho);
}


double drho0_dbeta_l_FD (int l, double x, double xp, 
			 double lambda, double beta)
{
  double delta = 1.0e-6;
  double rhoplus = exp(Log_rho0_l(l, x, xp, lambda, beta+delta));
  double rhominus = exp(Log_rho0_l(l, x, xp, lambda, beta-delta));
  return ((rhoplus - rhominus)/(2.0*delta));
}



void Test_rho0_Derivs()
{
  double xp = 1.235;
  double beta = 1.7321;
  double lambda = 0.5;
  int l = 0;
  for (double x=0.01; x<13.0; x+=0.01)
    {
      double analytic1 = drho0_dx_l(l, x, xp, lambda, beta);
      double FD1       = drho0_dx_l_FD(l, x, xp, lambda, beta);
      double error1 = fabs(analytic1-FD1) / analytic1;
      double analytic2 = d2rho0_dx2_l(l, x, xp, lambda, beta);
      double FD2       = d2rho0_dx2_l_FD(l, x, xp, lambda, beta);
      double error2 = fabs(analytic2-FD2) / analytic2;
      double dbeta_analytic = drho0_dbeta_l   (l, x, xp, lambda, beta);
      double dbeta_FD       = drho0_dbeta_l_FD(l, x, xp, lambda, beta);
      double dbeta_error = fabs (dbeta_analytic - dbeta_FD) / dbeta_analytic;
      fprintf (stderr, "%1.2f %1.8e %1.8e %1.6e\n", x, dbeta_analytic, 
	       dbeta_FD, dbeta_error); 
    }
}





double Log_rho0 (double r, double rp, double costheta, 
		 double lambda, double beta)
{
  double distsq = r*r + rp*rp - 2.0*r*rp*costheta;
  
  double rho0 = -1.5*log(4.0*M_PI*lambda*beta) -
    distsq/(4.0*lambda*beta);

  return (rho0);
}



double rho0 (double r, double rp, double costheta, double lambda,
	     double beta)
{
  double dist2 = r*r + rp*rp - 2.0*r*rp*costheta;
  return(pow(4.0*M_PI*lambda*beta, -1.5) *
    exp(-dist2 / (4.0*lambda*beta)));
}

double drho0_dbeta(double r, double rp, double costheta, double lambda,
		   double beta)
{
  double dist2 = r*r + rp*rp - 2.0*r*rp*costheta;
  return ((dist2/(4.0*lambda*beta) -1.5)/beta *
	  rho0(r, rp, costheta, lambda, beta));
}


double drho0_dbeta_scaled(double r, double rp, double costheta, 
			  double lambda, double beta)
{
  double dist2 = r*r + rp*rp - 2.0*r*rp*costheta;
  return ((dist2/(4.0*lambda*beta) -1.5)/beta);
}
