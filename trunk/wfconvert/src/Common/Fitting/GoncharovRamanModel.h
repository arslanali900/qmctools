#ifndef GONCHAROV_RAMAN_MODEL_H
#define GONCHAROV_RAMAN_MODEL_H

#include "../Blitz.h"

class GoncharovRamanModel
{
private:
  double b0, b1, b2, R0, R1, R2, nu0, nu1, nu2;
public:
  inline void SetParams (TinyVector<double,9> params)
  {
    b0  = params[0];  b1  = params[1];  b2  = params[2];
    R0  = params[3];  R1  = params[4];  R2  = params[5];
    nu0 = params[6];  nu1 = params[7];  nu2 = params[8];
  }

  inline void GetParams (TinyVector<double,9> &params)
  {
    params[0] =  b0;  params[1] =  b1;  params[2] =  b2;
    params[3] =  R0;  params[4] =  R1;  params[5] =  R2;
    params[6] = nu0;  params[7] = nu1;  params[8] = nu2;
  }

  
  inline TinyVector<double,9> 
  Grad_FD(TinyVector<double,2> PT)
  {
    TinyVector<double,9> save, grad, plus, minus;
    GetParams(save);
    for (int i=0; i<9; i++) {
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

    double b  =  b0 +  b1*T +  b2*T*T;
    double R  =  R0 +  R1*T +  R2*T*T;
    double nubar = nu0 + nu1*T + nu2*T*T;
    return nubar * pow(b*P/R + 1.0, 1.0/b);
  }

  inline TinyVector<double,9> 
  Grad(TinyVector<double,2> PT)
  {
    TinyVector<double,9> G, GFD;
    double P = PT[0];
    double T = PT[1];

    double b  =  b0 +  b1*T +  b2*T*T;
    double R  =  R0 +  R1*T +  R2*T*T;
    double nubar = nu0 + nu1*T + nu2*T*T;
    double nu = nubar * pow(b*P/R + 1.0, 1.0/b);

    // cerr << "b  = " << b << endl;
    // cerr << "R  = " << R << endl;
    // cerr << "nu = " << nu << endl;

    double dnu_db     = nu * (-1.0/(b*b) * log(b*P/R + 1.0) + (P/b)/(b*P+R));
    double dnu_dR     = -nubar*P/(R*R) * pow(b*P/R + 1.0, (1.0/b-1.0));
    double dnu_dnubar = pow(b*P/R + 1.0, 1.0/b);
    
    G[0] = dnu_db;
    G[1] = dnu_db*T;
    G[2] = dnu_db*T*T;

    G[3] = dnu_dR;
    G[4] = dnu_dR*T;
    G[5] = dnu_dR*T*T;

    G[6] = dnu_dnubar;
    G[7] = dnu_dnubar*T;
    G[8] = dnu_dnubar*T*T;

    GFD = Grad_FD(PT);
    // fprintf (stderr, "GFD = %16.12e %16.12e %12.12e\n", GFD[0], GFD[1], GFD[2]);
    // fprintf (stderr, "G   = %16.12e %16.12e %12.12e\n",   G[0],   G[1],   G[2]);

    // fprintf (stderr, "GFD = %16.12e %16.12e %12.12e\n", GFD[3], GFD[4], GFD[5]);
    // fprintf (stderr, "G   = %16.12e %16.12e %12.12e\n",   G[3],   G[4],   G[5]);

    // fprintf (stderr, "GFD = %16.12e %16.12e %12.12e\n", GFD[6], GFD[7], GFD[8]);
    // fprintf (stderr, "G   = %16.12e %16.12e %12.12e\n",   G[6],   G[7],   G[8]);

    return G;
  }
};




class GoncharovRamanModel2
{
private:
  double b0, b1, b2, b3, R0, R1, R2, R3, nu0, nu1, nu2, nu3;
public:
  inline void SetParams (TinyVector<double,12> params)
  {
    b0  = params[0];  b1  = params[1];   b2  = params[2];  b3  = params[3];
    R0  = params[4];  R1  = params[5];   R2  = params[6];  R3  = params[7];
    nu0 = params[8];  nu1 = params[9];   nu2 = params[10]; nu3 = params[11];
  }

  inline void GetParams (TinyVector<double,12> &params)
  {
    params[0] =  b0;  params[1] =  b1;  params[2] =  b2;  params[3]  = b3;
    params[4] =  R0;  params[5] =  R1;  params[6] =  R2;  params[7]  = R3;
    params[8] = nu0;  params[9] = nu1;  params[10] = nu2; params[11] = nu3;
  }

  
  inline TinyVector<double,12> 
  Grad_FD(TinyVector<double,2> PT)
  {
    TinyVector<double,12> save, grad, plus, minus;
    GetParams(save);
    for (int i=0; i<10; i++) {
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

    double b  =  b0 +  b1*T +  b2*T*T + b3*T*T*T;
    double R  =  R0 +  R1*T +  R2*T*T + R3*T*T*T;
    double nubar = nu0 + nu1*T + nu2*T*T + nu3*T*T*T;
    return nubar * pow(b*P/R + 1.0, 1.0/b);
  }

  inline TinyVector<double,12> 
  Grad(TinyVector<double,2> PT)
  {
    TinyVector<double,12> G, GFD;
    double P = PT[0];
    double T = PT[1];

    double b  =  b0 +  b1*T +  b2*T*T + b3*T*T*T;
    double R  =  R0 +  R1*T +  R2*T*T + R3*T*T*T;
    double nubar = nu0 + nu1*T + nu2*T*T + nu3*T*T*T;
    double nu = nubar * pow(b*P/R + 1.0, 1.0/b);

    if (b < 0.0)      cerr << "Negative b.\n";
    if (R < 0.0)      cerr << "Negative R.\n";
    if (nu < 0.0)      cerr << "Negative nu.\n";

    double dnu_db     = nu * (-1.0/(b*b) * log(b*P/R + 1.0) + (P/b)/(b*P+R));
    double dnu_dR     = -nubar*P/(R*R) * pow(b*P/R + 1.0, (1.0/b-1.0));
    double dnu_dnubar = pow(b*P/R + 1.0, 1.0/b);
    
    G[0] = dnu_db;
    G[1] = dnu_db*T;
    G[2] = dnu_db*T*T;
    G[3] = dnu_db*T*T*T;

    G[4] = dnu_dR;
    G[5] = dnu_dR*T;
    G[6] = dnu_dR*T*T;
    G[7] = dnu_dR*T*T*T;

    G[8]  = dnu_dnubar;
    G[9]  = dnu_dnubar*T;
    G[10] = dnu_dnubar*T*T;
    G[11] = dnu_dnubar*T*T*T;

    GFD = Grad_FD(PT);
    return G;
  }
};






#endif
