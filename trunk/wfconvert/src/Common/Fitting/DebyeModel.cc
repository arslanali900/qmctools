#include "DebyeModel.h"
#include "Fitting.h"

double
DebyeModel::Chi2_F (double theta, 
		  Array<double,1> &Fvals,
		  Array<double,1> &Tvals)
{
  SetTheta(theta);
  double chi2 = 0.0;
  for (int i=0; i<Fvals.size(); i++) {
    double diff = F(Tvals(i)) - Fvals(i);
    chi2 += diff *diff;
  }
  return chi2;
}

double
DebyeModel::Chi2_U (double theta, 
		  Array<double,1> &Uvals,
		  Array<double,1> &Tvals)
{
  SetTheta(theta);
  double chi2 = 0.0;
  for (int i=0; i<Uvals.size(); i++) {
    double diff = U(Tvals(i)) - Uvals(i);
    chi2 += diff *diff;
  }
  return chi2;
}


double
DebyeModel::MinTheta_F(Array<double,1> &Fvals, Array<double,1> &Tvals,
		      double startTheta, double endTheta)
{
  double minTheta = startTheta;
  double minChi2 = Chi2_F(minTheta, Fvals, Tvals);
  
  double dTheta = (endTheta - startTheta)/999.0;
  for (double theta = startTheta; theta <= (1.00000001*endTheta); 
       theta += dTheta) {
    double chi2 = Chi2_F(theta, Fvals, Tvals);
    // cerr << "Theta = " << theta << "  Chi2 = " << chi2 << endl;
    if (chi2 < minChi2) {
      minChi2 = chi2;
      minTheta = theta;
    }
  }
  // cerr << "minChi2 = " << minChi2 << endl;
  return minTheta;
}

double
DebyeModel::MinTheta_U(Array<double,1> &Uvals, Array<double,1> &Tvals,
		       double startTheta, double endTheta)
{
  double minTheta = startTheta;
  double minChi2 = Chi2_U(minTheta, Uvals, Tvals);
  
  double dTheta = (endTheta - startTheta)/999.0;
  for (double theta = startTheta; theta <= (1.00000001*endTheta); 
       theta += dTheta) {
    double chi2 = Chi2_U(theta, Uvals, Tvals);
    if (chi2 < minChi2) {
      minChi2 = chi2;
      minTheta = theta;
    }
  }
  return minTheta;
}

double
DebyeModel::OptTheta_F (Array<double,1> &Fvals, Array<double,1> &Tvals)
{
  double range = 7998.0;
  double mid = 4000.0;

  for (int i=0; i< 8; i++) {
    double start = mid - 0.5 * range;
    double end   = mid + 0.5 * range;
    mid = MinTheta_F (Fvals, Tvals, start, end);
    range = 0.1 * range;
  }
  return mid;
}

double
DebyeModel::OptTheta_U (Array<double,1> &Uvals, Array<double,1> &Tvals)
{
  double range = 7998.0;
  double mid = 4000.0;

  for (int i=0; i< 8; i++) {
    double start = mid - 0.5 * range;
    double end   = mid + 0.5 * range;
    mid = MinTheta_U (Uvals, Tvals, start, end);
    range = 0.1 * range;
  }
  return mid;
}

void
DebyeFreeEnergy::AddVolume_F (double Vval, vector<double> &Fvec,
			      vector<double> &Tvec)
{
  DebyeModel debye;
  debye.SetN(2);
  Array<double,1> F(Fvec.size()), T(Tvec.size());
  for (int i=0; i<Fvec.size(); i++) {
    F(i) = Fvec[i] - Fvec[0];
    T(i) = Tvec[i];
  }
  double theta = debye.OptTheta_F (F, T);
  V.push_back(Vval);
  F0.push_back(Fvec[0]);
  Theta.push_back(theta);
}


void
DebyeFreeEnergy::AddVolume_U (double Vval, vector<double> &Uvec,
			      vector<double> &Tvec)
{
  DebyeModel debye;
  debye.SetN(2);
  Array<double,1> U(Uvec.size()), T(Tvec.size());
  for (int i=0; i<Uvec.size(); i++) {
    U(i) = Uvec[i] - Uvec[0];
    T(i) = Tvec[i];
  }
  double theta = debye.OptTheta_U (U, T);
  V.push_back(Vval);
  U0.push_back(Uvec[0]);
  Theta.push_back(theta);
}


void
DebyeFreeEnergy::FitTheta_V (int numCoefs)
{
  int N = V.size();
  Array<double,2> basis(N, numCoefs);
  for (int i=0; i<N; i++) {
    basis(i,0) = 1.0;
    for (int j=1; j<numCoefs; j++)
      basis(i,j) = V[i]*basis(i,j-1);
  }
  ThetaCoefs.resize(numCoefs);
  Array<double,1> sigma(N), theta(N), errors(numCoefs);
  for (int i=0; i<N; i++) {
    sigma(i) = 1.0e-7;
    theta(i) = Theta[i];
  }
  LinFitLU(theta, sigma, basis, ThetaCoefs, errors);
}


double
DebyeFreeEnergy::F(double Vval, double T)
{
  Debye.SetTheta (Theta_V(Vval));
  return Debye.F(T);
}

double
DebyeFreeEnergy::U(double Vval, double T)
{
  Debye.SetTheta (Theta_V(Vval));
  return Debye.U(T);
}



double
DebyeFreeEnergy::P(double Vval, double T)
{
  Debye.SetTheta (Theta_V(Vval));
  return -dTheta_dV(Vval) * Debye.dF_dTheta(T);
}
  
double
DebyeFreeEnergy::P_FD (double Vval, double Tval)
{
  double eps = 1.0e-6;
  return ((F(Vval-eps,Tval)-F(Vval+eps,Tval))/(2.0*eps));
}  

double
DebyeFreeEnergy::dF_dT (double Vval, double Tval)
{
  Debye.SetTheta (Theta_V(Vval));
  return Debye.dF_dT (Tval);
}

double
DebyeFreeEnergy::d2F_dTheta_dT (double Vval, double Tval)
{
  Debye.SetTheta (Theta_V(Vval));
  return Debye.d2F_dTheta_dT (Tval);
}

double
DebyeFreeEnergy::d2F_dTheta_dT_FD (double Vval, double Tval)
{
  Debye.SetTheta (Theta_V(Vval));
  return Debye.d2F_dTheta_dT_FD (Tval);
}


double
DebyeFreeEnergy::dF_dT_FD (double Vval, double Tval)
{
  double eps = 1.0e-8;
  return ((F(Vval, Tval+eps) - F(Vval, Tval-eps))/(2.0*eps));
}

double
DebyeFreeEnergy::dF_dTheta (double Vval, double Tval)
{
  Debye.SetTheta (Theta_V(Vval));
  return Debye.dF_dTheta(Tval);
}


double
DebyeFreeEnergy::d2F_dTheta2 (double Vval, double Tval)
{
  Debye.SetTheta (Theta_V(Vval));
  return Debye.d2F_dTheta2(Tval);
//   double eps = 1.0e-7;
//   double theta = Theta_V(Vval);
//   Debye.SetTheta (theta + eps);
//   double plus = Debye.dF_dTheta(Tval);
//   Debye.SetTheta (theta - eps);
//   double minus = Debye.dF_dTheta (Tval);
//   return (plus-minus)/(2.0*eps);
}


double
DebyeFreeEnergy::K_T (double Vval, double Tval)
{
  return Vval * (dF_dTheta(Vval,Tval) * d2Theta_dV2(Vval) +
		 dTheta_dV(Vval) * dTheta_dV(Vval) * d2F_dTheta2(Vval, Tval));
}

double
DebyeFreeEnergy::K_T_FD(double Vval, double Tval)
{
  double eps = 1.0e-7;
  
  return -Vval *(P(Vval+eps,Tval)-P(Vval-eps,Tval))/(2.0*eps);
}

double
DebyeFreeEnergy::C_V (double Vval, double Tval)
{
  Debye.SetTheta(Theta_V(Vval));
  return Debye.C_V(Tval);
}


double
DebyeFreeEnergy::dP_dT (double Vval, double Tval)
{
  Debye.SetTheta(Theta_V(Vval));
  return -dTheta_dV(Vval) * Debye.d2F_dTheta_dT (Tval);
}

double
DebyeFreeEnergy::dP_dT_FD (double Vval, double Tval)
{
  double eps = 1.0e-7;
  return (P(Vval,Tval+eps)-P(Vval,Tval-eps))/(2.0*eps);
}

