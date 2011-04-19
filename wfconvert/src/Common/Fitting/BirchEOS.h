#ifndef BIRCH_EOS_H
#define BIRCH_EOS_H

#include <blitz/tinyvec.h>

using namespace blitz;

template<int N>
class BirchEOSClass
{
private:
  TinyVector<double,N> bn;
  double V0, E0;
  double third;
  double Atomic2GPa;
public:
  inline double operator()(double V);
  inline TinyVector<double,N+2> Grad(double V);
  inline TinyVector<double,N+2> GradFD(double V);
  inline double Pressure(double V);
  inline double PressureFD(double V);
  inline void SetParams (TinyVector<double,N+2> params);
  inline TinyVector<double,N+2> GetParams ();
  inline double GetB0();
  inline double GetB(double V);
  inline double GetBFD(double V);
  inline double GetB0FD();
  inline double GetB0p();
  inline double GetB0pFD();
  inline double GetBpFD(double V);
  inline double K_T (double V);

  BirchEOSClass() : third(1.0/3.0), Atomic2GPa(29421.01)
  { 
  }
};

template<int N>
inline void
BirchEOSClass<N>::SetParams (TinyVector<double,N+2> p)
{
  V0  = p[0];
  E0  = p[1];
  for (int i=0; i<N; i++)
    bn[i] = p[i+2];
}

template<int N>
inline TinyVector<double,N+2>
BirchEOSClass<N>::GetParams ()
{
  TinyVector<double,N+2> p;
  p[0] = V0;
  p[1] = E0;
  for (int i=0; i<N; i++)
    p[i+2] = bn[i];
  return p;
}


template<int N>
inline double
BirchEOSClass<N>::operator()(double V)
{
  double x = cbrt(V0/V);
  double t = x*x-1.0;
  double t_np2 = t*t;
  double E = E0;
  for (int n=0; n<N; n++) {
    E += bn[n] * t_np2;
    t_np2 *= t;
  }
  return E;
}


template<int N>
inline double
BirchEOSClass<N>::Pressure(double V)
{
  double x = cbrt(V0/V);
  double t = x*x-1;
  double tprime = -2.0/(3.0*V)*(t + 1.0);
  double P = 0.0;
  double t_np1 = t;
  for (int n=0; n<N; n++) {
    P += (double)(n+2)*bn[n]*t_np1;
    t_np1 *= t;
  }
  P *= tprime;

  return -Atomic2GPa * P;
}

template<int N>
inline double
BirchEOSClass<N>::K_T (double V)
{
  double eps = 1.0e-7;
  return -V*(Pressure(V+eps)-Pressure(V-eps))/(2.0*eps);
}

template<int N>
inline double
BirchEOSClass<N>::PressureFD(double V)
{
  return Atomic2GPa*((*this)(V+1.0e-5)-(*this)(V-1.0e-5))/-2.0e-5;
}


template<int N>
inline double
BirchEOSClass<N>::GetB0()
{
  return Atomic2GPa*8.0*bn[0]/(9.0*V0);
}

// template<int N>
// inline double
// BirchEOSClass<N>::GetB(double V)
// {
//   double x = cbrt(V/V0);
//   double t = x*x - 1.0;
//   double prefact = 2.0/(9.0*V)*(t + 1.0);
//   double sum = 0.0;
//   double t_n = 1.0;
//   for (int n=0; n<N; n++) {
//     sum += (double)(n+2)*((double)(2*n+2) + (double)(2*n+7)*t)*bn[n]*t_n;
//     t_n *= t;
//   }
//   return Atomic2GPa*prefact*sum;
// }

template<int N>
inline double
BirchEOSClass<N>::GetB(double V)
{
  double x = cbrt(V0/V);
  double t = x*x - 1.0;
  double prefact =  Atomic2GPa * 2.0/(9.0*V) * (t+1.0);
  double sum = 0.0;
  double t2n = 1.0;
  for (int n=0; n<N; n++) {
    sum += (double)(n+2)*t2n*bn[n]*(2.0*(n+1) + (double)(2*n+7)*t);
    t2n *= t;
  }
  return prefact*sum;
}

template<int N>
inline double
BirchEOSClass<N>::GetBFD(double V)
{
  double eps = 1.0e-7;
  double Pplus = Pressure(V+eps);
  double Pminus = Pressure(V-eps);
  return -V*(Pplus-Pminus)/(2.0*eps);
}


template<int N>
inline double
BirchEOSClass<N>::GetB0FD()
{
  double eps = 1.0e-7;
  double Pplus = Pressure(V0+eps);
  double Pminus = Pressure(V0-eps);
  return -V0 * (Pplus-Pminus)/(2.0*eps);
}

template<int N>
inline double
BirchEOSClass<N>::GetB0p()
{
  return 4.0 + Atomic2GPa*(16.0*bn[1])/(9.0*V0*GetB0());
}

template<int N>
inline double
BirchEOSClass<N>::GetB0pFD()
{
  double eps = 1.0e-7;
  double Bplus = GetB(V0+eps);
  double Bminus = GetB(V0-eps);
  double dB_dV = (Bplus - Bminus) / (2.0*eps);
  double B0 = GetB(V0);
  return -(V0/B0*dB_dV);
}

template<int N>
inline double
BirchEOSClass<N>::GetBpFD(double V)
{
  double eps = 1.0e-7;
  double Bplus = GetB(V+eps);
  double Bminus = GetB(V-eps);
  double dB_dV = (Bplus - Bminus) / (2.0*eps);
  double B = GetB(V);
  return -(V/B*dB_dV);
}


template<int N>
inline TinyVector<double,N+2>
BirchEOSClass<N>::Grad (double V)
{
  TinyVector<double,N+2> grad;
  grad[1] = 1.0;
  double x = cbrt(V0/V);
  double t = x*x - 1.0;

  // Find grad[0]
  grad[0] = 0.0;
  double dt_dV0 = 2.0/3.0*x*x/V0;
  double t_np1 = t;
  for (int n=0; n<N; n++) {
    grad[0] += bn[n] * (double)(n+2) * t_np1;
    t_np1 *= t;
  }
  grad[0] *= dt_dV0;
  
  // Find the rest
  double t_np2 = t*t;
  for (int n=0; n<N; n++) {
    grad[n+2] = t_np2;
    t_np2 *= t;
  }
  return grad;
}

template<int N>
inline TinyVector<double,N+2>
BirchEOSClass<N>::GradFD(double V)
{
  TinyVector<double,N+2> p, p_plus, p_minus, grad;
  double Eplus, Eminus;
  p = GetParams();
  for (int i=0; i<N+2; i++) {
    p_plus = p; p_minus = p;
    p_plus[i]  += 1.0e-7;
    p_minus[i] -= 1.0e-7;
    SetParams(p_plus);   Eplus  = (*this)(V);
    SetParams(p_minus);  Eminus = (*this)(V);
    grad[i] = (Eplus-Eminus)/(2.0e-7);
  }
  return grad;
}

#endif
