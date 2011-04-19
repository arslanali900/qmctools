#ifndef RAMAN_H
#define RAMAN_H

#include "../Splines/CubicNUBspline.h"
#include "../Integration/RungeKutta.h"
#include "../IO/IO.h"
#include "../Random/Random.h"

#include <vector>

class Vfunc
{
public:
  bool UseHarmonic;
  string fname;
  virtual double operator()(double x)=0;
  virtual void Read (IOSectionClass &in)=0;
  virtual double GetGridStart() = 0;
  virtual double GetGridEnd()   = 0;
  virtual void Perturb() {}
  virtual void Reset()   {}

  Vfunc() : UseHarmonic(false)
  { }
};

class VFitClass : public Vfunc
{
public:
  // Data from which the fit is derived.
  Array<double,1> Grid, E, Sigma;
  Array<double,2> F;
  void DoFit (Array<double,1> &Etry);

  TinyVector<double,4> Vcoefs;
  double xstart, xend;
  double Vstart, dVstart, d2Vstart, Vend, dVend, d2Vend;

  double operator()(double x);
  double Deriv     (double x);
  double Deriv2    (double x);
  void Read (IOSectionClass &in);
  double GetGridStart();
  double GetGridEnd();
};

class VSplineClass : public Vfunc
{
public:
  GeneralGrid SplineGrid;
  CubicNUBspline<GeneralGrid> Vspline;
  double Vmin, xmin, d2Vmin;
  double Vstart, dVstart, d2Vstart, Vend, dVend, d2Vend;
  double operator()(double x);
  void Read (IOSectionClass &in); 
  double GetGridStart();
  double GetGridEnd();

};

class PhononClass
{
private:
  GeneralGrid SplineGrid;
  LinearGrid IntegrationGrid;
  Array<Vec2,1> u_du;
  CubicNUBspline<GeneralGrid> Vspline;
  // lambda is hbar^2/(2*ReducedMass)
  double ReducedMass, lambda, lambdaInv;
  double Etrial;
  inline double V(double x);
  // This returns the cusp value at the midpoint.
  double Integrate();
  int CountNodes();
  double Vstart, dVstart, d2Vstart, Vend, dVend, d2Vend;
  double Vmin, xmin, d2Vmin;
  void Write();
  bool UseFit;
  RandomClass &Random;
public:
  Vfunc *Vfunction;
  inline Vec2 operator()(double x, Vec2 u);
  void Read (IOSectionClass &in);

  // Returns the eigenenergy
  double Solve (int desiredNodes);
  TinyVector<double,2> SolveError (int desiredNoes);
  PhononClass(RandomClass &rand) : Random(rand)
  {
  }
};

class RamanSpectrum
{
private:
  vector<double> Energies;
  const double kB, h, c;
  double T, Sigma;
  double V;
public:
  inline void SetTemp(double temp)     { T = temp;              }
  inline void SetVolume (double vol)   { V = vol;               }
  inline double GetVolume()            { return V;              }
  inline void SetWidth (double sigma)  { Sigma = sigma;         }
  inline void AddEnergy(double E)      { Energies.push_back(E); } 
  inline void SetEnergies (vector<double> E) { Energies = E;    }
  inline double Spectrum (double nu);
  // Returns intensity of a peak relative to the 0->1 transition
  inline double Intensity (int n);
  inline double MeanFrequency();
  
  RamanSpectrum() :
    kB(3.1668152e-06), c(2.9979246e+10), h(1.5198298e-16)
  {  }
};


// class DatchiRamanModel
// {
// private:
//   // This is fixed
//   double B0p;
//   // These can vary
//   double a, b, c, d, e, gamma;
// public:
//   inline void SetB0p (double b0p) { B0p = b0p; }

//   inline double operator() (TinyVector<double,2> PT)
//   {
//     double P = PT[0];
//     double T = PT[1];
//     double nu0 = a*T + b*T*T + c;
//     double B0 = d + e*T;
//     double x = 1.0 + B0p * P / B0;
//     return nu0 * pow(x, gamma/B0p);
//   }

//   inline void SetParams (TinyVector<double,6> params)
//   {
//     a = params[0];    b = params[1];
//     c = params[2];    d = params[3];
//     e = params[4];    gamma = params[5];
//   }

//   inline void GetParams (TinyVector<double,6> &params)
//   {
//     params[0] = a;    params[1] = b;
//     params[2] = c;    params[3] = c;
//     params[4] = e;    params[5] = gamma;
//   }

  
//   inline TinyVector<double,6> 
//   Grad_FD(TinyVector<double,2> PT)
//   {
//     TinyVector<double,6> save, grad, plus, minus;
//     GetParams(save);
//     for (int i=0; i<6; i++) {
//       double eps = max(1.0e-6 * fabs(save[i]), 1.0e-10);
//       plus = minus = save;
//       plus[i]  = save[i] + eps;
//       minus[i] = save[i] - eps;
//       SetParams(plus);
//       double nu_plus  = (*this)(PT);
//       SetParams(minus);
//       double nu_minus = (*this)(PT);
//       grad[i] = (nu_plus - nu_minus) / (2.0*eps);
//     }
//     SetParams(save);
//     return grad;
//   }

//   inline TinyVector<double,6> 
//   Grad(TinyVector<double,2> PT)
//   {
//     double P = PT[0];
//     double T = PT[1];
//     TinyVector<double,6> grad;
    
//     double nu0 = a*T + b*T*T + c;
//     double B0 = d + e*T;
//     double x = 1.0 + B0p * P / B0;

//     grad[0] = T   * pow(x, gamma/B0);
//     grad[1] = T*T * pow(x, gamma/B0);
//     grad[2] =       pow(x, gamma/B0);
//     double dnu_dB0 = - nu0 * gamma * pow(x, gamma/B0p-1.0) * B0p*P/(B0*B0);
//     grad[3] = dnu_dB0;
//     grad[4] = T * dnu_dB0;
//     grad[5] = nu0 * log(gamma/B0p)/B0p * pow(x, gamma/B0p);
//     return grad;
//   }
// };


class RamanModel
{
private:
  RamanSpectrum Spectrum;
  const double kB, h, c;
  const int NumEnergies, NumCoefs;
  RandomClass Random;
  CommunicatorClass Comm;
  // Polynomial fit for E_i(1/V)
  // First index is the energy level, i.  The second index is the
  // exponent of (1/V).
  Array<double,2> EnergyCoefs;
  Array<double,1> E0Coefs;
public:
  // First index is the volume.  The second index is the phonon energy
  // level, n.
  vector<Array<double,1> > Energies;
  vector<Array<double,1> > Sigma;
  vector<double> Volume;



  void AddVolume (string fname);
  double MeanFrequency(double V, double T);
  void Read(string fname);
  void DoFits();

  inline double E(int i, double V)
  {
    double Vinv = 1.0/V;
    double basis = 1.0;
    double e = 0.0;
    for (int n=0; n<NumCoefs; n++) {
      e += E0Coefs(n) * basis;
      e += EnergyCoefs(i,n) * basis;
      basis *= Vinv;
    }
    return e;
  }

  inline double Frequency (int i, double V)
  {
    const double h = 1.5198298e-16;
    const double c = 2.9979246e+10;
    double delta = E(i + 1, V) - E(i,V);
    return delta/(h*c);
  }


  RamanModel() : NumEnergies(11), NumCoefs(3), Random(Comm),
		 kB(3.1668152e-06), c(2.9979246e+10), h(1.5198298e-16) 
  {
    Random.Init();
  }
};


inline double
RamanSpectrum::MeanFrequency()
{
  double Z = 0.0;
  double sum = 0.0;
  // The initial state is n and the final state is n-1
  // Thus, the frequency is (E_n - E_{n-1})/(h*c)
  // The matrix element is proportial to sqrt(n), so
  // we multiply by n for Fermi's Golden Rule.  This could
  // be done more precisely by integrating <psi_{n+1}|x|psi_n>.
  for (int n=1; n<Energies.size(); n++) {
    double dE = (Energies[n] - Energies[0]);
    // HACK HACK HACK
    double a = exp(-dE/(kB * T));
    double factor = (double)n;
    Z += factor * a;
    double nu0 = (Energies[n] - Energies[n-1])/(h*c);
    sum += nu0 * factor * a;
  }
  return sum/Z;

  // HACK HACK HACK
  //  return (Energies[1] - Energies[0])/(h*c);
}


inline double
RamanSpectrum::Intensity(int n)
{
  double Z = 0.0;
  for (int i=0; i<Energies.size(); i++) {
    double dE = (Energies[i] - Energies[0]);
    double a = exp(-dE/(kB * T));
    double factor = (double)(i+1);
    Z += factor * a;
  }
  return (double)n*exp((Energies[0]-Energies[n-1])/(kB*T))/Z;
}


inline double
RamanSpectrum::Spectrum(double nu)
{
  double spec = 0.0;
  for (int n=1; n<Energies.size(); n++) {
    double nu0 = (Energies[n] - Energies[n-1])/(h*c);
    double I = Intensity (n);
    spec += I*sqrt(M_PI/(Sigma*Sigma))*exp(-(nu-nu0)*(nu-nu0)/(2.0*Sigma*Sigma));
  }
  return spec;
}



inline double
PhononClass::V(double x)
{
  return (*Vfunction)(x);

  // if (UseHarmonic) 
  //   return Vmin + 0.5*d2Vmin*(x-xmin)*(x-xmin);
  // else if (x < SplineGrid.Start) {
  //   double dx = x - SplineGrid.Start;
  //   return Vstart + dVstart * dx + 0.5*d2Vstart*dx*dx;
  // }
  // else if (x >= SplineGrid.End) {
  //   double dx = x - SplineGrid.End;
  //   return Vend + dVend*dx + 0.5*d2Vend*dx*dx;
  // }
  // else
  //   return Vspline(x);
}

inline Vec2
PhononClass::operator()(double x, Vec2 u)
{
  double du, d2u;

  du  = u[1];
  d2u = (V(x)-Etrial)*lambdaInv * u[0];
  return Vec2(du, d2u);
}

#endif
