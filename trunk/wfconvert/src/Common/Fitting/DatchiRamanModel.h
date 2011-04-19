#ifndef DATCHI_RAMAN_MODEL_H
#define DATCHI_RAMAN_MODEL_H

#include "../Blitz.h"

class DatchiRamanModel
{
private:
  // This is fixed
  double B0p;
  // These can vary
  double a, b, c, d, e, gamma;
public:
  inline void SetB0p (double b0p) { B0p = b0p; }

  inline void SetParams (TinyVector<double,6> params)
  {
    a = params[0];    b = params[1];
    c = params[2];    d = params[3];
    e = params[4];    gamma = params[5];
  }

  inline void GetParams (TinyVector<double,6> &params)
  {
    params[0] = a;    params[1] = b;
    params[2] = c;    params[3] = c;
    params[4] = e;    params[5] = gamma;
  }

  
  inline TinyVector<double,6> 
  Grad_FD(TinyVector<double,2> PT)
  {
    TinyVector<double,6> save, grad, plus, minus;
    GetParams(save);
    for (int i=0; i<6; i++) {
      double eps = max(1.0e-6 * fabs(save[i]), 1.0e-10);
      plus = minus = save;
      plus[i]  = save[i] + eps;
      minus[i] = save[i] - eps;
      SetParams(plus);
      double nu_plus  = (*this)(PT);
      SetParams(minus);
      double nu_minus = (*this)(PT);
      grad[i] = (nu_plus - nu_minus) / (2.0*eps);
    }
    SetParams(save);
    return grad;
  }

  inline double operator() (TinyVector<double,2> PT)
  {
    double P = PT[0];
    double T = PT[1];
    double nu0 = a*T + b*T*T + c;
    double B0 = d + e*T;
    double x = 1.0 + B0p * P / B0;
    return nu0 * pow(x, gamma/B0p);
  }

  inline TinyVector<double,6> 
  Grad(TinyVector<double,2> PT)
  {
    // return Grad_FD(PT);
    // grad = 6 [    978.673    979651  0.977695 0.0230118   23.0348  -19.1968 ]

    double P = PT[0];
    double T = PT[1];
    TinyVector<double,6> grad;
    
    double nu0 = a*T + b*T*T + c;
    double B0 = d + e*T;
    double x = 1.0 + B0p * P / B0;

    grad[0] = T   * pow(x, gamma/B0p);
    grad[1] = T*T * pow(x, gamma/B0p);
    grad[2] =       pow(x, gamma/B0p);
    double dnu_dB0 = - nu0 * gamma/B0p * pow(x, gamma/B0p-1.0) * B0p*P/(B0*B0);
    grad[3] = dnu_dB0;
    grad[4] = T * dnu_dB0;
    grad[5] = nu0 * log(x)/B0p * pow(x, gamma/B0p);
    // cerr << " grad   = " << grad << endl;
    // cerr << " gradFD = " << Grad_FD(PT) << endl;
    return grad;
  }
};



class DatchiRamanModelNoT
{
private:
  // This is fixed
  double B0p;
  // These can vary
  double c, d, gamma;
public:
  inline void SetB0p (double b0p) { B0p = b0p; }

  inline void SetParams (TinyVector<double,3> params)
  {
    c     = params[0];    d = params[1];
    gamma = params[2];
  }

  inline void GetParams (TinyVector<double,3> &params)
  {
    params[0] = c;    params[1] = d;
    params[2] = gamma;
  }

  
  inline TinyVector<double,3> 
  Grad_FD(double P)
  {
    TinyVector<double,3> save, grad, plus, minus;
    GetParams(save);
    for (int i=0; i<3; i++) {
      double eps = max(1.0e-6 * fabs(save[i]), 1.0e-10);
      plus = minus = save;
      plus[i]  = save[i] + eps;
      minus[i] = save[i] - eps;
      SetParams(plus);
      double nu_plus  = (*this)(P);
      SetParams(minus);
      double nu_minus = (*this)(P);
      grad[i] = (nu_plus - nu_minus) / (2.0*eps);
    }
    SetParams(save);
    return grad;
  }

  inline double operator() (double P)
  {
    double nu0 = c;
    double B0 = d;
    double x = 1.0 + B0p * P / B0;
    return nu0 * pow(x, gamma/B0p);
  }

  inline TinyVector<double,3> 
  Grad(double P)
  {
    TinyVector<double,3> grad;
    
    double nu0 = c;
    double B0 = d;
    double x = 1.0 + B0p * P / B0;

    // fprintf (stderr, "x = %1.5e  nu0 = %1.5e  B0 = %1.5e  gamma=%1.5e\n",
    // 	     x, nu0, B0, gamma);

    grad[0] =       pow(x, gamma/B0p);
    double dnu_dB0 = - nu0 * gamma/B0p * pow(x, gamma/B0p-1.0) * B0p*P/(B0*B0);
    grad[1] = dnu_dB0;
    grad[2] = nu0 * log(x)/B0p * pow(x, gamma/B0p);
    // cerr << " grad   = " << grad << endl;
    // cerr << " gradFD = " << Grad_FD(P) << endl;
    return grad;
  }
};



#endif
