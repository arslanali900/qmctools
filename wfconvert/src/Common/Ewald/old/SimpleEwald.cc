#include "SimpleEwald.h"

void SimpleEwaldClass::BreakUp(double cutoff)
{
  CutOff = cutoff;
}


double SimpleEwaldClass::Vshort(double r)
{
  assert (CutOff != 0.0);
  assert (Pot != NULL);
  double alpha = 3.5/CutOff;
  double Vlong;
  if (r <= 0.0)
    Vlong = 2.0/sqrt(M_PI)*Z*alpha;
  else 
    Vlong = Z/r*erf(alpha*r);
  return Pot->V(r)-Vlong;
}


double SimpleEwaldClass::Vlong_k (double k)
{
  double alpha = 3.5/CutOff;
  if (k <= 0.0)
    k = 1.0e-30;
  return 4.0/(k*k)*exp(-k*k/(4.0*alpha*alpha));
}
