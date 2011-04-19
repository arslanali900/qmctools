#ifndef DEBYE_MODEL_H
#define DEBYE_MODEL_H

#include "../Blitz.h"
#include "gsl/gsl_sf_debye.h"
#include <vector>

class DebyeModel
{
private:
  double Theta, ThetaInv;
  double kB; // In Hartrees per Kelvin
  double N;
  double Chi2_F (double theta, Array<double,1> &Fvals,
		 Array<double,1> &Tvals);
  double MinTheta_F (Array<double,1> &Fvals, Array<double,1> &Tvals,
		     double startTheta, double endTheta);
  double Chi2_U (double theta, Array<double,1> &Fvals,
		 Array<double,1> &Tvals);
  double MinTheta_U (Array<double,1> &Fvals, Array<double,1> &Tvals,
		     double startTheta, double endTheta);
public:
  inline void SetN    (int n);
  inline void SetTheta(double theta);
  inline double F(double T);
  inline double U(double T);
  inline double dF_dTheta        (double T);
  inline double dF_dT            (double T);
  inline double dU_dTheta        (double T);
  inline double d2F_dTheta_dT    (double T);
  inline double d2F_dTheta_dT_FD (double T);
  inline double d2F_dTheta2      (double T);
  // Heat capacity
  inline double C_V (double T);
  double OptTheta_F (Array<double,1> &Fvals, Array<double,1> &Tvals);
  double OptTheta_U (Array<double,1> &Uvals, Array<double,1> &Tvals);
  inline double CalcTheta (double F, double T);
  DebyeModel() : kB(3.16681526543384e-06), N(2)
  {
  }
};


  
class DebyeFreeEnergy
{
private:
  vector<double> V, Theta, F0, U0;
  DebyeModel Debye;
  Array<double,1> ThetaCoefs;
public:
  void AddVolume_F (double V, vector<double> &F,
		    vector<double> &T);
  void AddVolume_U (double V, vector<double> &U,
		    vector<double> &T);
  inline double Theta_V (double V);
  inline double dTheta_dV (double V);
  inline double d2Theta_dV2 (double V);
  
  void FitTheta_V(int numCoefs);
  // Helmholtz free energy
  double F(double V, double T);
  double dF_dT (double V, double T);
  double dF_dT_FD (double V, double T);
  double dF_dTheta (double V, double T);
  double d2F_dTheta_dT (double V, double T);
  double d2F_dTheta_dT_FD (double V, double T);
  double d2F_dTheta2 (double V, double T);
  // Internal energy
  double U(double V, double T);
  
  // Thermal pressure
  double P(double V, double T);
  double P_FD(double V, double T);
  double dP_dV (double V, double T);
  double dP_dT (double V, double T);
  double dP_dT_FD (double V, double T);

  // Bulk modulus
  double K_T (double V, double T);
  double K_T_FD (double V, double T);

  // Heat capacity
  double C_V (double V, double T);

  DebyeFreeEnergy()
  {
    Debye.SetN(2);
  }
};

inline double 
DebyeFreeEnergy::Theta_V(double V)
{
  double V2i = 1.0;
  double tsum = 0.0;
  for (int i=0; i<ThetaCoefs.size(); i++) {
    tsum += ThetaCoefs(i) * V2i;
    V2i *= V;
  }
  return tsum;
}

inline double 
DebyeFreeEnergy::dTheta_dV(double V)
{
  double V2im1 = 1.0;
  double tsum = 0.0;
  for (int i=1; i<ThetaCoefs.size(); i++) {
    tsum += ThetaCoefs(i) * (double)i * V2im1;
    V2im1 *= V;
  }
  return tsum;
}

inline double
DebyeFreeEnergy::d2Theta_dV2 (double V)
{
  double V2im2 = 1.0;
  double tsum = 0.0;
  for (int i=2; i<ThetaCoefs.size(); i++) {
    tsum += ThetaCoefs(i) * (double)(i*(i-1))*V2im2;
    V2im2 *= V;
  }
  return tsum;
}


inline double
DebyeModel::CalcTheta(double Fval, double T)
{
  double thetaLow = 0.0; 
  double thetaHigh = 10000.0;

  while ((thetaHigh - thetaLow) > 1.0e-10) {
    double theta = 0.5*(thetaHigh + thetaLow);
    SetTheta (theta);
    double Ftry = F(T) - F(0.0);
    //cerr << "Ftry = " << Ftry << "  Fval = " << Fval << endl;
    if (Ftry < Fval) 
      thetaLow = theta;
    else
      thetaHigh = theta;
  }
  return 0.5*(thetaHigh + thetaLow);
}


inline void
DebyeModel::SetN(int n)
{
  N = (double)n;
}

inline void
DebyeModel::SetTheta (double theta)
{
  Theta = theta; 
  ThetaInv = 1.0/theta;
}

inline double
DebyeModel::F(double T)
{
  double Tstar = T * ThetaInv;
  if (Tstar == 0.0)
    Tstar = 1.0e-10;
  
  double D = gsl_sf_debye_3(1.0/Tstar);
  double f = N*kB*T*(3.0*log1p(-exp(-1.0/Tstar)) - D  /*-9.0/(8.0*Tstar) */);
  // cerr << "theta = " << Theta << "  T = " << T << "  F = " << f << endl;
//   fprintf (stderr, "theta = %1.5f   T = %1.5f   F = %1.8e  D = %1.8e\n", 
// 	   Theta, T, F);
  return f;
}

inline double
DebyeModel::U(double T)
{
  double x = Theta/T;
  if (isinf(x) || x < 0.0)
    x = 1.0e10;
  return 3.0 * N * kB * T * gsl_sf_debye_3(x);
}

inline double
DebyeModel::C_V(double T)
{
  double x = (T == 0) ? 1.0e10 : Theta/T;
  return 3.0*N*kB *
    (4.0*gsl_sf_debye_3(x) - 3.0*x/expm1(x));
}

inline double
DebyeModel::dF_dTheta (double T)
{
  double x = T * ThetaInv;
  double D3 = gsl_sf_debye_3(1.0/x);
  
  return ThetaInv * 3.0*N*kB*T*D3;
}


inline double
DebyeModel::dU_dTheta (double T)
{
  double eps = 1.0e-6;
  double temp = Theta;
  Theta = temp + eps;
  double Uplus = U(T);
  Theta = temp - eps;
  double Uminus = U(T);
  Theta = temp;
  return ((Uplus - Uminus)/(2.0*eps));
}


inline double
DebyeModel::dF_dT (double T)
{
  double x = T * ThetaInv;
  return N*kB*(3.0*log1p(-exp(-Theta/T)) - 4.0*gsl_sf_debye_3(1.0/x));
}

inline double
DebyeModel::d2F_dTheta_dT (double T)
{
  double x = Theta / T;
  return 3.0*N*kB*(4.0*gsl_sf_debye_3(x)/Theta - 3.0/(T*(expm1(x))));
}

inline double
DebyeModel::d2F_dTheta_dT_FD (double T)
{
  double eps = 1.0e-8;
  return (dF_dTheta(T+eps)-dF_dTheta(T-eps))/(2.0*eps);
}


inline double
DebyeModel::d2F_dTheta2 (double T)
{
  return 9.0 * N * kB/ (Theta * expm1(Theta/T))
    - 12.0*N*kB*T/(Theta*Theta) * gsl_sf_debye_3(Theta/T);
}


// Model fit for Debye temperature as a function of T
class DebyeTemp
{
private:
  double T0, alpha, c;
public:
  inline double operator()(double T);
  inline TinyVector<double,3> Grad(double T);
  inline void SetParams(TinyVector<double,3> params);
  inline TinyVector<double,3> GetParams();

  DebyeTemp() : T0 (2000.0), alpha(0.005), c(-100.0)
  { }

};

inline double
DebyeTemp::operator()(double T)
{
  return T0 + c * exp(-alpha * T);
}

inline TinyVector<double,3> 
DebyeTemp::Grad(double T)
{
  TinyVector<double,3> grad;
  grad[0] = 1.0;
  grad[1] = -c * T * exp(-alpha*T);
  grad[2] = exp(-alpha * T);
  return grad;
}

inline void
DebyeTemp::SetParams(TinyVector<double,3> params)
{
  T0    = params[0];
  alpha = params[1];
  c     = params[2];
}


inline TinyVector<double,3>
DebyeTemp::GetParams()
{
  return TinyVector<double,3>(T0, alpha, c);
}

#endif
