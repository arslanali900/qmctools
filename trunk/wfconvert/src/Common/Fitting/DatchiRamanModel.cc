#include "DatchiRamanModel.h"
#include "NonlinearFitTemp.h"
#include <iostream>
#include <fstream>
#include <strstream>
#include <cstdio>
#include <vector>

void FitRaman(string fname)
{
  ifstream fin;
  fin.open(fname.c_str());
  vector<double> Pvec, Tvec, nuvec, sigmavec;


  while (!fin.eof()) {
    char line[500];
    fin.getline(line,500);
    if (line[0] == '#') continue;

    istrstream sin(line);
    double V, T, P, nu, sigma;
    sin >> V >> T >> P >> nu >> sigma;
    Pvec.push_back(P);
    Tvec.push_back(T);
    nuvec.push_back(nu);
    sigmavec.push_back(sigma);
  }
  int N = Pvec.size();
  Array<double,1> nu(N), sigma(N);
  Array<TinyVector<double,2>,1> PT(N);
  for (int i=0; i<N; i++) {
    nu(i) = nuvec[i];
    sigma(i) = sigmavec[i];
    PT(i) = TinyVector<double,2>(Pvec[i],Tvec[i]);
  }

  DatchiRamanModel model, DatchiModel;
  model.SetB0p (3.86);
  DatchiModel.SetB0p(3.7);
  NonlinearFitClass<6, DatchiRamanModel, TinyVector<double,2> > fitter(model);
  TinyVector<double,6> params, Datchi_params;
  
  // params = 6 [ -0.0173179-1.06158e-05   1062.27   466.607-0.0380165   1.28938 ]
  params[0] = -0.02;
  params[1] = -1.82e-5;
  params[2] = 1050.3;
  params[3] = 361;
  params[4] = -0.02;
  params[5] = 1.388;

  Datchi_params[0] = -0.01;
  Datchi_params[1] = -1.42e-5;
  Datchi_params[2] = 1058.3;
  Datchi_params[3] = 381;
  Datchi_params[4] = -0.03;
  Datchi_params[5] = 1.188;


  DatchiModel.SetParams(Datchi_params);
  model.SetParams(params);
  fitter.Fit (PT, nu, sigma, params);
  model.SetParams(params);
  fprintf (stderr, "a     = %12.5f\n", params[0]);
  fprintf (stderr, "b     = %12.5e\n", params[1]);
  fprintf (stderr, "c     = %12.5f\n", params[2]);
  fprintf (stderr, "d     = %12.5f\n", params[3]);
  fprintf (stderr, "e     = %12.5f\n", params[4]);
  fprintf (stderr, "gamma = %12.5f\n", params[5]);

  FILE *fout = fopen ("DatchiRamanModel.dat", "w");
  
  for (double P=0; P<=1000.0; P+=1.0) {
    TinyVector<double,2> PT(P, 300.0);
    double my_nu = model (PT);
    PT[0] = P;
    double Datchi_nu  = DatchiModel(PT);
    fprintf (fout, "%10.4f  %12.4f %12.4f\n",
	     P, my_nu, Datchi_nu);
  }
  fclose (fout); 


  cerr << "params = " << params << endl;
}

main(int argc, char **argv)
{
  if (argc > 1)
    FitRaman(argv[1]);
}
