#include "Chebyshev.h"
#include <cstdio>
#include <cmath>

void TestChebyshev()
{
  int N = 6;
  ChebyshevPoly<double> cheb(N);

  FILE *fout   = fopen ("Chebyshev.dat", "w");
  FILE *dfout  = fopen ("dChebyshev.dat", "w");
  FILE *d2fout = fopen ("d2Chebyshev.dat", "w");
  for (double x=-1.0; x<=1.0; x+= 0.0001) {
    fprintf (fout, "%9.4f ", x);    
    fprintf (dfout, "%9.4f ", x);
    fprintf (d2fout, "%9.4f ", x);
    for (int i=0; i<N; i++) {
      for (int j=0; j<N; j++)
	cheb[j] = (double)(i==j);
      double P, dP, d2P, dP_FD, d2P_FD;
      cheb(x, P, dP, d2P);
      dP_FD = (cheb(x+1.0e-6) - cheb(x-1.0e-6))/2.0e-6;
      d2P_FD = (cheb(x+1.0e-6) + cheb(x-1.0e-6) - 2.0*cheb(x))  /1.0e-6;
      fprintf (fout, "%16.12e ", P);
      fprintf (dfout, "%16.12e %1.16e ", dP, dP_FD);
      fprintf (d2fout, "%16.12e %1.16e ", d2P, d2P_FD);
    }
    fprintf (fout, "\n");
    fprintf (dfout, "\n");
    fprintf (d2fout, "\n");
  }  
  fclose(fout);
  fclose(dfout);
  fclose(d2fout);
}


inline double
testfunc(double x)
{
  double k = 8.1943294112345;
  return cos(k*x);
}


void
TestFit ()
{
  double k = 8.1943294112345;
  int nquad = 100;
  ChebyshevQuad rule(nquad);
  int N = 21;
  ChebyshevPoly<double> cheb(N);
  std::vector<double> Tn(N);

  for (int n=0; n<N; n++)
    cheb[n] = 0.0;

  for (int iq=0; iq<nquad; iq++) {
    double x = rule.x(iq);
    //    fprintf (stderr, "x = %12.10f\n", x);
    cheb(x, Tn);
    for (int n=0; n<N; n++)
      cheb[n] += rule.w(iq) *  Tn[n] * testfunc(x);
  }
  cheb[0] /= M_PI;
  for (int n=1; n<N; n++)
    cheb[n] /= 0.5*M_PI;

  for (int i=0; i<N; i++)
    fprintf (stderr, "  cheb[%d] = %1.10e\n", i, cheb[i]);
  
  FILE *fout = fopen("testfit.dat", "w");
  for (double x=-1.0; x<=1.0; x+=0.001) {
    double f = testfunc(x);
    double c = cheb(x);
    fprintf (fout, "%8.3f  %20.18e %20.18e\n", x, f, c);
  }
  fclose(fout);
}




main()
{
  TestChebyshev();
  TestFit();
}
