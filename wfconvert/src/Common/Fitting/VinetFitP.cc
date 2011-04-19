#include "VinetEOS2.h"
#include "NonlinearFit.h"
#include "../IO/IO.h"
#include <vector>
#include "PhononFreeEnergy.h"
#include "DebyeModel.h"
#include "RamanFreq.h"
#include "../Random/Random.h"

void
TestGrad()
{
  TinyVector<double,3> params(80.0, 300.0, 3.0);
  VinetPressure eos;
  eos.SetParams (params);
  TinyVector<double,3> grad   = eos.Grad(10.0);
  TinyVector<double,3> gradFD = eos.GradFD(10.0);

  cerr << "grad   = " << grad << endl;
  cerr << "gradFD = " << gradFD << endl;
}

void
TestVal()
{
  TinyVector<double,3> params(80.0, -0.5, 0.296471, 3.378);
  VinetPressure eos;
  eos.SetParams (params);
  
  for (double V=15.0; V<=65.0; V+=0.1)
    fprintf (stderr, "%18.14f %18.14f\n", V, eos(V));
}

void
TestFit()
{
  TinyVector<double,3> params(80.0, -0.5, 0.296471, 3.378);
  TinyVector<double,3> params2(81.0, -0.6, 0.286471, 3.9);
  VinetPressure eos;
  eos.SetParams (params);
  NonlinearFitClass<3,VinetPressure> fitter(eos);

  Array<double,1> V(5), E(5), sigma(5);
  V = 15.0, 25.0, 35.0, 50.0, 65.0;
  for (int i=0; i<E.size(); i++) {
    E(i) = eos(V(i));
    sigma(i) = 0.01;
  }

  fitter.Fit (V, E, sigma, params2);
  cerr << "params2 = " << params2 << endl;
}




void 
EOSFit (string fname, VinetPressure &eos)
{
  TinyVector<double,3> params(80.0, 300.0, 3.0), errors;  
  eos.SetParams (params);
  NonlinearFitClass<3,VinetPressure> fitter(eos);
  Array<double,1> V, P, sigma;

  IOSectionClass in;
  assert(in.OpenFile (fname));
  assert(in.ReadVar("V0",  params[0]));
  assert(in.ReadVar("B0",  params[1]));
  assert(in.ReadVar("B0p", params[2]));
  assert(in.ReadVar("V", V));
  assert(in.ReadVar("P", P));
  assert(in.ReadVar("Sigma", sigma));
  cerr << "V = " << V << endl;
  cerr << "P = " << P << endl;
  cerr << "Sigma = " << sigma << endl;

  fprintf (stdout, "V0  = %12.8f\n", params[0]);
  fprintf (stdout, "Ec  = %12.8f\n", params[1]);
  fprintf (stdout, "B0  = %12.8f GPa\n", params[2] * 29421.01);
  eos.SetParams (params);
  fprintf (stdout, "B0p = %12.8f\n", params[3]);
  fitter.Fit (V, P, sigma, params);

  Array<double,2> &covar = fitter.GetCovariance();
  errors[0] = sqrt (covar(0,0));
  errors[1] = sqrt (covar(1,1));
  errors[2] = sqrt (covar(2,2));
  fprintf (stdout, "V0  = %12.8f +/- %12.8f\n", params[0], errors[0]);
  fprintf (stdout, "B0  = %12.8f +/- %12.8f\n", 
	   params[1], errors[1]);
  eos.SetParams(params);
  fprintf (stdout, "B0p = %12.8f +/- %12.8f\n", params[2], errors[2]);
  double V0 = params[0];
  double latConst = cbrt(4.0*V0);
  double latConstError = cbrt(4.0)/(3.0*cbrt(V0)*cbrt(V0)) * errors[0];
  fprintf (stdout, "Lattice const. = %1.6f +/- %1.6f bohr radii\n",
	   latConst, latConstError);
  fprintf (stdout, "               = %1.6f +/- %1.6f angstrom\n",  
	   0.52917721*latConst, 0.52917721*latConstError);
  string Ename = fname + ".dat";


  FILE *fout = fopen (Ename.c_str(), "w");
  double Vmin=V(0), Vmax=V(0), delta;
  for (int i=0; i<V.size(); i++) {
    Vmin = min(Vmin, V(i));
    Vmax = max(Vmax, V(i));
  }
  cerr << "Vmin = " << Vmin << "  Vmax = " << Vmax << endl;

  delta = (Vmax - Vmin) / 1000.0;
  for (double v=Vmin; v<=Vmax; v+=delta)
    fprintf (fout, "%20.16f %20.16f %20.16f %20.16f\n", v, eos(v), eos(v));
  fclose (fout);

  // // Write out Goncharovs residuals
  // //  FILE *fin = fopen ("GoncharovEOSData.png.dat", "r");
  // FILE *fin = fopen ("GoncharovEOS.Holzapfel.dat", "r");
  // double Vval, Pval, aval;
  // if (fin) {
  //   fout = fopen ("GoncharovResiduals.dat", "w");
  //   while (fscanf (fin, "%lf %lf\n", &Pval, &aval) == 2) {
  //     //      Vval *= 79.7535137356486;
  //     aval *= 1.8897261;  // Angstrom to bohr radius
  //     Vval = 0.25 * aval * aval * aval;
  //     double diff = Pval - eos.Pressure(Vval);
  //     fprintf (fout, "%1.12e %1.12e\n", Vval, diff);
  //   }
  //   fclose(fout);
  //   fclose(fin);
  // }

  // // Write out Datchi's residuals
  // fin = fopen ("DatchiEOS.dat", "r");
  // if (fin) {
  //   fout = fopen ("DatchiResiduals.dat", "w");
  //   double aVal, Tval;
  //   char dummyLine[100];
  //   fgets (dummyLine, 100, fin);
  //   fgets (dummyLine, 100, fin);
  //   while (fscanf (fin, "%lf %lf %lf %lf\n", &Pval, &Tval, &aVal, &Vval) == 4) {
  //     Vval = 1.8897261*1.8897261*1.8897261*aVal*aVal*aVal/4.0;
  //     double diff = Pval - eos.Pressure(Vval);
  //     fprintf (fout, "%1.12e %1.12e\n", Vval, diff);
  //   }
  //   fclose(fout);
  //   fclose (fin);
  // }

  // // Write out Knittles's residuals
  // fin = fopen ("Exp300K_P-V.dat", "r");
  // if (fin) {
  //   fout = fopen ("KnittleResiduals.dat", "w");
  //   while (fscanf (fin, "%lf %lf\n", &Vval, &Pval) == 2) {
  //     Vval *= 79.714027629;
  //     // Vval *= 79.7535137356486;
  //     double diff = Pval - eos.Pressure(Vval);
  //     fprintf (fout, "%1.12e %1.12e\n", Vval, diff);
  //   }
  //   fclose(fout);
  //   fclose (fin);
  // }
}




double FindV (double P, double T,
	      VinetPressure &EOS)
{
  double Vlo = 10.0;
  double Vhi = 150.0;
  double au2GPa = 29421.01;

  while ((Vhi - Vlo) > 1.0e-10) {
    double Vtry = 0.5*(Vlo + Vhi);
    double Ptry = EOS(Vtry);
    if (Ptry > P)
      Vlo = Vtry;
    else
      Vhi = Vtry;
  }
  return 0.5*(Vlo+Vhi);
}




main(int argc, char **argv)
{
  TestGrad();
  if (argc > 1) {
    VinetPressure EOS;
    EOSFit (argv[1], EOS);
  }
  else
    cerr << "Usage:  VinetFitP pressure.in\n";
}
