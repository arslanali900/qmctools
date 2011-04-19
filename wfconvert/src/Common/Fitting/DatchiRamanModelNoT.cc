#include "DatchiRamanModel.h"
#include "NonlinearFit.h"
#include <iostream>
#include <fstream>
#include <strstream>
#include <cstdio>
#include <vector>

void FitRaman(string fname)
{
  ifstream fin;
  fin.open(fname.c_str());
  vector<double> Pvec, nuvec, sigmavec;


  while (!fin.eof()) {
    char line[500];
    fin.getline(line,500);
    if (line[0] == '#') continue;

    istrstream sin(line);
    double V, T, P, nu, sigma = 0.1;
    sin >> P >> nu;
    Pvec.push_back(P);
    nuvec.push_back(nu);
    sigmavec.push_back(sigma);
  }
  int N = Pvec.size();
  Array<double,1> nu(N), sigma(N), P(N);
  Array<Vec2,1> PT(N);
  for (int i=0; i<N; i++) {
    nu(i) = nuvec[i];
    sigma(i) = sigmavec[i];
    P(i) = Pvec[i];
    PT(i) = TinyVector<double,2>(Pvec[i],300.0);
  }

  DatchiRamanModelNoT Goncharov;
  DatchiRamanModel DatchiModel;
  //Goncharov.SetB0p (3.86);  
  Goncharov.SetB0p (3.7);
  DatchiModel.SetB0p(3.7);
  NonlinearFitClass<3, DatchiRamanModelNoT> fitter(Goncharov);
  TinyVector<double,3> params;
  TinyVector<double,6> Datchi_params;
  
  params[0] = 1050.3;
  params[1] = 361;
  params[2] = 1.388;

  Datchi_params[0] = -0.01;
  Datchi_params[1] = -1.42e-5;
  Datchi_params[2] = 1058.3;
  Datchi_params[3] = 381;
  Datchi_params[4] = -0.03;
  Datchi_params[5] = 1.188;

  DatchiModel.SetParams(Datchi_params);
  Goncharov.SetParams(params);
  fitter.Fit (P, nu, sigma, params);
  Goncharov.SetParams(params);
  fprintf (stderr, "a     = %12.5f\n", params[0]);
  fprintf (stderr, "b     = %12.5e\n", params[1]);
  fprintf (stderr, "c     = %12.5f\n", params[2]);
  fprintf (stderr, "d     = %12.5f\n", params[3]);
  fprintf (stderr, "e     = %12.5f\n", params[4]);
  fprintf (stderr, "gamma = %12.5f\n", params[5]);

  FILE *fout = fopen ("GoncharovDatchi.dat", "w");
  
  for (double P=0; P<=1000.0; P+=1.0) {
    TinyVector<double,2> PT(P, 300.0);
    double my_nu = Goncharov (P);
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
