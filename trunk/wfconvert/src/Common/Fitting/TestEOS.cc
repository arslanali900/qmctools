#include "LiEOS.h"
#include "NonlinearFit.h"
#include "../IO/IO.h"

void
TestGrad()
{
  TinyVector<double,4> params(29.19, 4.42, 0.296471, 3.378);

  LiEOSClass eos;
  eos.SetParams (params);
  TinyVector<double,4> grad   = eos.Grad(10.0);
  TinyVector<double,4> gradFD = eos.GradFD(10.0);

  cerr << "grad   = " << grad << endl;
  cerr << "gradFD = " << gradFD << endl;
}

void
TestVal()
{
  TinyVector<double,4> params(29.19, 4.42, 0.296471, 3.378);
  LiEOSClass eos;
  eos.SetParams (params);
  
  for (double V=15.0; V<=65.0; V+=0.1)
    fprintf (stderr, "%18.14f %18.14f\n", V, eos(V));
}

void
TestFit()
{
  TinyVector<double,4> params(29.19, 4.42, 0.296471, 3.378);  
  TinyVector<double,4> params2(25.0, 4.0, 0.26, 3.8);
  LiEOSClass eos;
  eos.SetParams (params);
  NonlinearFitClass<4,LiEOSClass> fitter(eos);

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
  TinyVector<double,4> params(80.0, 0.51, 0.0155, 3.0);  
  LiEOSClass eos;
  eos.SetParams (params);
  NonlinearFitClass<4,LiEOSClass> fitter(eos);
  Array<double,1> V, E, sigma;


  IOSectionClass in;
  assert(in.OpenFile (fname));
  assert(in.ReadVar("V", V));
  assert(in.ReadVar("E", E));
  assert(in.ReadVar("Sigma", sigma));
  cerr << "V = " << V << endl;
  cerr << "E = " << E << endl;
  cerr << "Sigma = " << sigma << endl;

  fprintf (stdout, "V0  = %12.8f\n", params[0]);
  fprintf (stdout, "Ec  = %12.8f\n", params[1]);
  fprintf (stdout, "B0  = %12.8f\n", params[2]);
  fprintf (stdout, "B0p = %12.8f\n", params[3]);
  fitter.Fit (V, E, sigma, params);
  fprintf (stdout, "V0  = %12.8f\n", params[0]);
  fprintf (stdout, "Ec  = %12.8f\n", params[1]);
  fprintf (stdout, "B0  = %12.8f GPa\n", params[2] * 29421.01);
  fprintf (stdout, "B0p = %12.8f\n", params[3]);
  eos.SetParams(params);
  fprintf (stdout, "Delta = %1.5f\n", eos.GetDelta());
  fprintf (stdout, "Lattice const. = %1.6f bohr radii\n", cbrt(4.0*params[0]));
  fprintf (stdout, "               = %1.6f angstrom\n",  
	   0.52917721*cbrt(4.0*params[0]));
  string Ename = fname + ".dat";

  FILE *fout = fopen (Ename.c_str(), "w");
  double delta = (V(0) - V(V.size()-1))/1000.0;
  for (double v=V(V.size()-1); v<=V(0); v+=delta)
    fprintf (fout, "%20.16f %20.16f %20.16f %20.16f\n", v, eos(v), eos.Pressure(v), eos.PressureFD(v));
  fclose (fout);
}



main(int argc, char **argv)
{
  //TestGrad();
  //TestVal();
  if (argc > 1)
    TestFit2(argv[1]);
  else
    cerr << "Usage:  eosfit fname.in\n";
}
