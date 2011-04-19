#ifndef SHIFTED_CHEBYSHEV_H
#define SHIFTED_CHEBYSHEV_H
#include <vector>
#include <cmath>

class ShiftedChebyshevQuad
{
private:
  std::vector<double> Abcissas, Weights;
public:
  inline void SetRule (int num);
  inline double x(int i) { return Abcissas[i]; }
  inline double w(int i) { return Weights[i];  }
  ShiftedChebyshevQuad()        {                     }
  ShiftedChebyshevQuad(int n)   { SetRule(n);         }
};

void
ShiftedChebyshevQuad::SetRule (int n)
{
  Abcissas.resize(n);
  Weights.resize(n);
  double ninv = 1.0/(double)n;
  double norm = 0.5*M_PI/(double)n;
  double w = M_PI/(double)n;
  for (int i=0; i<n; i++) {
    double k = (double)(i+1);
    Abcissas[i] = 0.5*(std::cos((2.0*k-1.0)*norm) + 1.0);
    Weights[i]  = w;
  }
}




template<typename T>
class ShiftedChebyshevPoly
{
private:
  std::vector<T> Coefs;
public:
  inline int size() { return Coefs.size(); }
  inline T operator()(double x);
  inline void operator() (double x, T &P, T& dP, T& d2P);
  inline void operator() (double x, std::vector<double> &Tn);
  void operator()(double x, std::vector<double> &Tn,
		  std::vector<double> &dTn, std::vector<double> &d2Tn);
  inline void resize(int N)  {  Coefs.resize(N); }
  const T& operator[](int n) const { return Coefs[n]; }
  T&       operator[](int n)       { return Coefs[n]; }

  ShiftedChebyshevPoly() {}
  ShiftedChebyshevPoly(int N) { resize(N); }
};

template<typename T> inline T
ShiftedChebyshevPoly<T>::operator()(double x)
{
  x = 2.0*x - 1.0;
  T sum = T();
  sum = Coefs[0];
  double Tnm1 = 1.0;
  double Tn   = x;
  for (int n=1; n<size(); n++) {
    sum += Coefs[n] * Tn;
    double Tnp1 = 2.0*x*Tn - Tnm1;
    Tnm1 = Tn;
    Tn = Tnp1;
  }
  return sum;
}


template<typename T> inline void
ShiftedChebyshevPoly<T>::operator()(double x,  std::vector<double> &Tn)
{
  x = 2.0*x - 1.0;
  Tn[0] = 1.0;
  Tn[1] = x;
  for (int n=1; n<(Tn.size()-1); n++) 
    Tn[n+1] = 2.0*x*Tn[n] - Tn[n-1];
}


template<typename T> inline void
ShiftedChebyshevPoly<T>::operator()(double x,  
			     std::vector<double> &Tn,
			     std::vector<double> &dTn,
			     std::vector<double> &d2Tn)
{
  x = 2.0*x - 1.0;
  Tn[0] = 1.0; dTn[0] = 0.0; d2Tn[0] = 0.0;
  Tn[1] = x;   dTn[1] = 1.0; d2Tn[0] = 0.0;
  for (int n=1; n<(Tn.size()-1); n++) {
    Tn[n+1]   = 2.0*x*Tn[n] - Tn[n-1];
    dTn[n+1]  = 2.0*(Tn[n] + x*dTn[n]) - dTn[n-1];
    d2Tn[n+1] = 4.0*dTn[n] + 2.0*x*d2Tn[n] - d2Tn[n-1];
  }
  for (int n=1; n<Tn.size(); n++) {
    dTn[n]  *= 2.0;
    d2Tn[n] *= 4.0;
  }
    
}


template<typename T> inline void
ShiftedChebyshevPoly<T>::operator()(double x, T& P, T& dP, T& d2P)
{
  x = 2.0*x - 1.0;
  P = dP = d2P = T();
  P = Coefs[0];

  double Tnm1 = 1.0;
  double Tn   = x;
  double dTn=1.0, dTnm1=0.0, d2Tn=0.0, d2Tnm1=0.0;
  for (int n=1; n<size(); n++) {
    P   +=     Coefs[n] *   Tn;
    dP  += 2.0*Coefs[n] *  dTn;
    d2P += 4.0*Coefs[n] * d2Tn;
    double   Tnp1 = 2.0*x*Tn - Tnm1;
    double  dTnp1 = 2.0*(Tn + x*dTn) - dTnm1;
    double d2Tnp1 = 2.0*(2.0*dTn + x*d2Tn) - d2Tnm1;
    Tnm1   =   Tn;    Tn   =   Tnp1;
    dTnm1  =  dTn;    dTn  =  dTnp1;
    d2Tnm1 = d2Tn;    d2Tn = d2Tnp1;
  }
}

#endif
