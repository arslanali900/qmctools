#include "VinetEOS.h"
#include "NonlinearFit.h"
#include "../IO/IO.h"

void
TestGrad()
{
  TinyVector<double,3> params(29.19, 4.42, 0.296471, 3.378);

  VinetEOSClass eos;
  eos.SetParams (params);
  TinyVector<double,3> grad   = eos.Grad(10.0);
  TinyVector<double,3> gradFD = eos.GradFD(10.0);

  cerr << "grad   = " << grad << endl;
  cerr << "gradFD = " << gradFD << endl;
}

void
TestVal()
{
  TinyVector<double,3> params(29.19, 4.42, 0.296471, 3.378);
  VinetEOSClass eos;
  eos.SetParams (params);
  
  for (double V=15.0; V<=65.0; V+=0.1)
    fprintf (stderr, "%18.14f %18.14f\n", V, eos(V));
}

void
TestFit()
{
  TinyVector<double,3> params(29.19, 4.42, 0.296471, 3.378);  
  TinyVector<double,3> params2(25.0, 4.0, 0.26, 3.8);
  VinetEOSClass eos;
  eos.SetParams (params);
  NonlinearFitClass<3,VinetEOSClass> fitter(eos);

  Array<double,1> V(5), E(5), sigma(5);
  V = 15.0, 25.0, 35.0, 50.0, 65.0;
  for (int i=0; i<E.size(); i++) {
    E(i) = eos(V(i));
    sigma(i) = 0.01;
  }

  fitter.Fit (V, E, sigma, params2);
  cerr << "params2 = " << params2 << endl;
}


using namespace IO;
void
TestFit2(string fname)
{
  TinyVector<double,3> params(80.0, 0.51, 0.0155), errors;  
  VinetEOSClass eos;
  eos.SetParams (params);
  NonlinearFitClass<3,VinetEOSClass> fitter(eos);
  Array<double,1> V, E, sigma;


  IOSectionClass in;
  assert(in.OpenFile (fname));
  assert(in.ReadVar("V0", params[0]));
  assert(in.ReadVar("Ec", params[1]));
  assert(in.ReadVar("B0", params[2]));
  params[2] /= 29421.01;
	
  assert(in.ReadVar("V", V));
  assert(in.ReadVar("E", E));
  assert(in.ReadVar("Sigma", sigma));
  cerr << "V = " << V << endl;
  cerr << "E = " << E << endl;
  cerr << "Sigma = " << sigma << endl;

  fprintf (stdout, "V0  = %12.8f\n", params[0]);
  fprintf (stdout, "Ec  = %12.8f\n", params[1]);
  fprintf (stdout, "B0  = %12.8f GPa\n", params[2] * 29421.01);
  eos.SetParams (params);
  fprintf (stdout, "B0p = %12.8f\n", eos.GetB0p());
  fitter.Fit (V, E, sigma, params);

  Array<double,2> &covar = fitter.GetCovariance();
  errors[0] = sqrt (covar(0,0));
  errors[1] = sqrt (covar(1,1));
  errors[2] = sqrt (covar(2,2));
  fprintf (stdout, "V0  = %12.8f +/- %12.8f\n", params[0], errors[0]);
  fprintf (stdout, "Ec  = %12.8f +/- %12.8f\n", params[1], errors[1]);
  fprintf (stdout, "B0  = %12.8f +/- %12.8f\n", 
	   params[2]* 29421.01, errors[2]* 29421.01);
  eos.SetParams(params);
  fprintf (stdout, "B0p = %12.8f\n", eos.GetB0p());
  fprintf (stdout, "Delta = %1.5f\n", 0.0);
  double V0 = params[0];
  double latConst = cbrt(4.0*V0);
  double latConstError = cbrt(4.0)/(3.0*cbrt(V0)*cbrt(V0)) * errors[0];
  fprintf (stdout, "Lattice const. = %1.6f +/- %1.6f bohr radii\n",
	   latConst, latConstError);
  fprintf (stdout, "               = %1.6f +/- %1.6f angstrom\n",  
	   0.52917721*latConst, 0.52917721*latConstError);
  string Ename = fname + ".dat";

  

  FILE *fout = fopen (Ename.c_str(), "w");
  double delta = (V(0) - V(V.size()-1))/1000.0;
  for (double v=V(V.size()-1); v<=V(0); v+=delta)
    fprintf (fout, "%20.16f %20.16f %20.16f %20.16f\n", v, eos(v), eos.Pressure(v), eos.PressureFD(v));
  fclose (fout);
}



main(int argc, char **argv)
{
  TestGrad();
  //TestVal();
  if (argc > 1)
    TestFit2(argv[1]);
  else
    cerr << "Usage:  eosfit fname.in\n";
}
