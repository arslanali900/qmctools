#include "MieGruineisen.h"
#include "NonlinearFitTemp.h"

void
DoMGFit (string fname)
{
  FILE *fin = fopen (fname.c_str(), "r");
  assert (fin != NULL);
  vector<double> Vvec, Tvec, Fvec, Uvec;
  bool done = false;
  while (!done) {
    double v,t,f,u;
    int retval = fscanf (fin, "%lf %lf %lf %lf", &v, &t, &f, &u);
    if (retval == 4) {
      Vvec.push_back (v);
      Tvec.push_back (t);
      Fvec.push_back (f);
      Uvec.push_back (u);
    }
    else
      done = true;
  }
  int N = Vvec.size();

  Array<TinyVector<double,2>,1> VT(N);
  Array<double,1> F(N), sigma(N);
  double F0;
  for (int i=0; i<N; i++) {
    if (fabs(Tvec[i]) < 1.0e-6)
      F0 = Fvec[i];
    VT(i) = TinyVector<double,2> (Vvec[i],Tvec[i]);
    F(i) = Fvec[i]-F0;
    sigma(i) = 1.0e-6;
  }
  
  TinyVector<double,4> params (1700.0, 1.04, 4.0, 0.0001);
  MieGruineisen MG;
  MG.SetParams(params);
  MG.SetV0 (80.0);
  
  NonlinearFitClass<4,MieGruineisen,TinyVector<double,2> >
    Fitter(MG);

  Fitter.Fit(VT, F, sigma, params);

  fprintf (stderr, "Theta  = %1.8f\n", params[0]);
  fprintf (stderr, "Gamma0 = %1.8f\n", params[1]);
  fprintf (stderr, "q      = %1.8f\n", params[2]);
  fprintf (stderr, "F0     = %1.8f\n", params[3]);

  MG.SetParams(params);
  
  FILE *fout = fopen ("PhonE.dat", "w");
  for (int i=0; i<N; i++) {
    double Fmodel = MG.F(Vvec[i], Tvec[i]);
    fprintf (fout, "%5.1f  %8.1f  %16.12e %16.12e\n",
	     Vvec[i], Tvec[i], F(i), Fmodel);
  }
  fclose(fout);
  
}


main(int argc, char *argv[])
{
  DoMGFit (argv[1]);
}
