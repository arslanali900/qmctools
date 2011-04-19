#include "PhononFreeEnergy.h"
#include "Fitting.h"

void
PhononFreeEnergy::EvalBasis(double V, double T,
			    Array<double,1> basis)
{
  double V2i = 1.0;
  assert (basis.size() == Vorder*Torder);
  for (int i=0; i<Vorder; i++) {
    double T2j = 1.0;
    for (int j=0; j<Torder; j++) {
      int index = (i*Torder)+j;
      basis (index) = V2i*T2j;
      T2j *= T;
    }
    V2i *= V;
  }
}


void
PhononFreeEnergy::FitF(int vorder, int torder,
		       Array<double,1> &F,
		       Array<double,1> &V,
		       Array<double,1> &T)
{
  Vorder = vorder;
  Torder = torder;
  FCoefs.resize(Vorder*Torder);
  Btmp.resize(Vorder*Torder);
  Array<double,1> errors(Vorder*Torder);
  Array<double,1> sigma(F.size());
  Array<double,2> basis(F.size(),Vorder*Torder);
  for (int i=0; i<F.size(); i++) {
//     fprintf (stderr, "F(i) = %1.8f T(i) = %1.8f  V(i) = %1.8f\n", 
// 	     F(i), T(i), V(i));
    sigma(i) = 1.0e-5;
    EvalBasis(V(i), T(i), basis(i,Range::all()));
//     if (T(i) != 0.0)
//       F(i) += 3.0 * kB*T(i)*log(T(i));
    // cerr << "basis(" << i << " = " << basis(i,Range::all()) << endl;
  }
  //LinFitSVD (F, sigma, basis, FCoefs, errors, 1.0e-50);
  LinFitLU (F, sigma, basis, FCoefs, errors);
  // cerr << "FCoefs = " << FCoefs << endl;
}

double 
PhononFreeEnergy::F_VT(double V, double T)
{
  EvalBasis (V, T, Btmp);
  double F = 0.0;
  for (int i=0; i<FCoefs.size(); i++)
    F += FCoefs(i)*Btmp(i);

//   F -= 3.0 * kB*T*log(T);
  return F;
}

double 
PhononFreeEnergy::P_VT(double V, double T)
{
  double V2im1 = 1.0;
  double T2j = 1.0;
  double P = 0.0;
  for (int i=1; i<Vorder; i++) {
    T2j = 1.0;
    for (int j=0; j<Torder; j++) {
      int index = (i*Torder) + j;
      P += T2j * (double)i * V2im1 * FCoefs(index);
      T2j *= T;
    }
    V2im1 *= V;
  }
  return P;
}

double
PhononFreeEnergy::P_VT_FD(double V, double T)
{
  double eps = 1.0e-6;
  return (1.0/(2.0*eps)*(F_VT(V+eps,T) - F_VT(V-eps,T)));
}

double
PhononFreeEnergy::dP_dT (double V, double T)
{
  double V2im1 = 1.0;
  double dP = 0.0;
  for (int i=1; i<Vorder; i++) {
    double T2jm1 = 1.0;
    for (int j=1; j<Torder; j++) {
      int index = (i*Torder)+j;
      dP += FCoefs(index)*(double)(i*j)*V2im1*T2jm1;
      T2jm1 *= T;
    }
    V2im1 *= V;
  }
  return dP;
}

double
PhononFreeEnergy::dP_dT_FD (double V, double T)
{
  double eps = 1.0e-6;
  return (1.0/(2.0*eps)*(P_VT(V,T+eps) - P_VT(V,T-eps)));
}

double
PhononFreeEnergy::dP_dV (double V, double T)
{
  double V2im2 = 1.0;
  double dP = 0.0;
  for (int i=2; i<Vorder; i++) {
    double T2j = 1.0;
    for (int j=0; j<Torder; j++) {
      int index = (i*Torder)+j;
      dP += FCoefs(index)*(double)(i*(i-1))*V2im2+T2j;
      T2j *= T;
    }
    V2im2 *= i;
  }
  return dP;
}

double
PhononFreeEnergy::dP_dV_FD (double V, double T)
{
  double eps = 1.0e-6;
  return (1.0/(2.0*eps)*(P_VT(V+eps,T) - P_VT(V-eps,T)));
}
