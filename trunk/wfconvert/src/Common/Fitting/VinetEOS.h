#ifndef LI_EOS_H
#define LI_EOS_H

#include <blitz/tinyvec.h>

using namespace blitz;

class VinetEOSClass
{
private:
  double V0, Ec, B0, E0;
  double third;

public:
  inline double operator()(double V);
  inline TinyVector<double,3> Grad(double V);
  inline TinyVector<double,3> GradFD(double V);
  inline double Pressure(double V);
  inline double PressureFD(double V);
  inline void SetParams (TinyVector<double,3> params);
  inline TinyVector<double,3> GetParams ();
  inline double GetB0p();
  VinetEOSClass() : third(1.0/3.0)
  { 
  }
};

inline void
VinetEOSClass::SetParams (TinyVector<double,3> p)
{
  V0  = p[0];
  Ec  = p[1];
  B0  = p[2];
}

inline TinyVector<double,3>
VinetEOSClass::GetParams ()
{
  TinyVector<double,3> p;
  p[0] = V0;
  p[1] = Ec;
  p[2] = B0;
  return p;
}



inline double
VinetEOSClass::operator()(double V)
{
  double x = cbrt(V/V0);
  double eta = sqrt(9.0*B0*V0/Ec);
  double a = eta*(x-1.0);
  
  double E = -Ec*(1.0+a)*exp(-a);
  return E;
}


inline double
VinetEOSClass::Pressure(double V)
{
  double x = cbrt(V/V0);
  double eta = sqrt(9.0*B0*V0/Ec);
  double a = eta*(x-1.0);
  double dx_dV = third * cbrt(1.0/(V0*V*V));
  double da_dV = eta * dx_dV;
  
  double dE_dV = -Ec*da_dV*exp(-a) + Ec*(1.0 + a)*exp(-a)*da_dV;
  double E = -Ec*(1.0 + a)*exp(-a);
  return -29421.01*dE_dV;
}


inline double
VinetEOSClass::PressureFD(double V)
{
  return 29421.01*((*this)(V+1.0e-5)-(*this)(V-1.0e-5))/-2.0e-5;
}

inline double
VinetEOSClass::GetB0p()
{
  double eta = sqrt(9.0*B0*V0/Ec);
  return 1.0 + (2.0/3.0)*eta;
}


inline TinyVector<double,3>
VinetEOSClass::VinetEOSClass::Grad (double V)
{
  double x = cbrt(V/V0);
  double eta = sqrt(9.0*B0*V0/Ec);
  double a = eta*(x-1.0);
  
  double dx_dV0      = -third *cbrt(V/(V0*V0*V0*V0));
  double deta_dB0    =  0.5*sqrt(9.0*V0/(Ec*B0));
  double deta_dV0    =  0.5*sqrt(9.0*B0/(Ec*V0));
  double deta_dEc    = -0.5*sqrt(9.0*B0*V0/(Ec*Ec*Ec));
  double da_dV0      =  deta_dV0*(x-1.0) + eta * dx_dV0;
  double da_dB0      =  deta_dB0*(x-1.0);
  double da_dEc      =  deta_dEc*(x-1.0);

  double dE_dV0  = -Ec*da_dV0*exp(-a) +Ec*(1.0+a)*exp(-a)*da_dV0; 
  double dE_dEc  = -(1.0+a)*exp(-a) 
    -Ec*da_dEc*exp(-a) +Ec*(1.0+a)*exp(-a)*da_dEc;
  double dE_dB0  =  -Ec*da_dB0*exp(-a) + Ec*(1.0+a)*exp(-a)*da_dB0;

  return TinyVector<double,3> (dE_dV0, dE_dEc, dE_dB0);
}

inline TinyVector<double,3>
VinetEOSClass::GradFD(double V)
{
  TinyVector<double,3> p, p_plus, p_minus, grad;
  double Eplus, Eminus;
  p = GetParams();
  for (int i=0; i<3; i++) {
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
