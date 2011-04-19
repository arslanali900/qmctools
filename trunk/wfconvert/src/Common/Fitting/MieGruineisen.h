#ifndef MIE_GRUINEISEN_H
#define MIE_GRUINEISEN_H

#include "DebyeModel.h"

class MieGruineisen
{
private:
  DebyeModel Debye;
  double Theta0, Gamma0, q, V0, F0;
public:
  // This is used in the nonlinear fit.
  inline double operator()(TinyVector<double,2> V_T) 
  { return F(V_T[0], V_T[1]); }

  inline double P(double V, double T);
  inline double F(double V, double T); 
  inline TinyVector<double,4> GradF(double V, double T);
  inline TinyVector<double,4> Grad(TinyVector<double,2> V_T)
  { return GradF(V_T[0], V_T[1]); }

  inline void SetParams (TinyVector<double,4> params);
  inline TinyVector<double,4> GetParams();

  inline void SetV0 (double v0) { V0 = v0; }
  MieGruineisen() 
  {
    Debye.SetN(2);
  }
};

inline void
MieGruineisen::SetParams(TinyVector<double,4> params)
{
  Theta0 = params[0];
  Gamma0 = params[1];
  q      = params[2];
  F0     = params[3];
  Debye.SetTheta(Theta0);
}

inline TinyVector<double,4>
MieGruineisen::GetParams()
{
  return TinyVector<double,4>(Theta0, Gamma0, q, F0);
}



inline double
MieGruineisen::P(double V, double T)
{
  double gamma = Gamma0*pow(V/V0, q);
  return gamma*Debye.U(T)/V;
}


inline double
MieGruineisen::F(double V, double T)
{
//   fprintf (stderr, "Theta0 = %1.8f\n", Theta0);
//   fprintf (stderr, "Gamma0 = %1.8f\n", Gamma0);
//   fprintf (stderr, "q      = %1.8f\n",      q);
//   fprintf (stderr, "F0     = %1.8f\n", F0);
//   fprintf (stderr, "V0     = %1.8f\n", V0);

  double u = Debye.U(T);
//   fprintf (stderr, "U = %1.8e\n", u);

  double f = -Debye.U(T)*(Gamma0*pow(V/V0, q)/q + F0);
//   fprintf (stderr, "F= %1.8e\n", f);
  return f;
}

inline TinyVector<double,4>
MieGruineisen::GradF(double V, double T)
{
  TinyVector<double,4> grad, save, temp;
  double eps = 1.0e-6;
  save = GetParams();
  for (int i=0; i<4; i++) {
    temp = save;
    temp[i] = save[i] + eps;
    SetParams(temp);
    double Fplus = F(V,T);
    temp[i] = save[i] - eps;
    SetParams(temp);
    double Fminus = F(V,T);
    SetParams(save);
    grad[i] = (Fplus-Fminus)/(2.0*eps);
  }
  return grad;




  TinyVector<double,4> G;
  double DU = -Debye.U(T);
  double powVq = pow(V/V0,q)/q;

  G[0] = -Debye.dU_dTheta(T)*(Gamma0*powVq + F0);
  G[1] = -DU * powVq;
  G[2] = -DU*Gamma0*(log(V/V0)*powVq - powVq/q);
  G[3] = -DU;
  return G;
}


#endif
