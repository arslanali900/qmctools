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
  TinyVector<double,4> params(80.0, -0.57, 0.011726, 3.378);
  VinetEOSClass eos;
  eos.SetParams (params);
  TinyVector<double,4> grad   = eos.Grad(10.0);
  TinyVector<double,4> gradFD = eos.GradFD(10.0);

  cerr << "grad   = " << grad << endl;
  cerr << "gradFD = " << gradFD << endl;
}

void
TestVal()
{
  TinyVector<double,4> params(80.0, -0.5, 0.296471, 3.378);
  VinetEOSClass eos;
  eos.SetParams (params);
  
  for (double V=15.0; V<=65.0; V+=0.1)
    fprintf (stderr, "%18.14f %18.14f\n", V, eos(V));
}

void
TestFit()
{
  TinyVector<double,4> params(80.0, -0.5, 0.296471, 3.378);
  TinyVector<double,4> params2(81.0, -0.6, 0.286471, 3.9);
  VinetEOSClass eos;
  eos.SetParams (params);
  NonlinearFitClass<4,VinetEOSClass> fitter(eos);

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
  TinyVector<double,4> params(80.0, 0.51, 0.0155), errors;  
  VinetEOSClass eos;
  eos.SetParams (params);
  NonlinearFitClass<4,VinetEOSClass> fitter(eos);
  Array<double,1> V, E, sigma;


  IOSectionClass in;
  assert(in.OpenFile (fname));
  assert(in.ReadVar("V0",  params[0]));
  assert(in.ReadVar("E0",  params[1]));
  assert(in.ReadVar("B0",  params[2]));
  assert(in.ReadVar("B0p", params[3]));
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
  fprintf (stdout, "B0p = %12.8f\n", params[3]);
  fitter.Fit (V, E, sigma, params);

  Array<double,2> &covar = fitter.GetCovariance();
  errors[0] = sqrt (covar(0,0));
  errors[1] = sqrt (covar(1,1));
  errors[2] = sqrt (covar(2,2));
  errors[3] = sqrt (covar(3,3));
  fprintf (stdout, "V0  = %12.8f +/- %12.8f\n", params[0], errors[0]);
  fprintf (stdout, "E0  = %12.8f +/- %12.8f\n", params[1], errors[1]);
  fprintf (stdout, "B0  = %12.8f +/- %12.8f\n", 
	   params[2]* 29421.01, errors[2]* 29421.01);
  eos.SetParams(params);
  fprintf (stdout, "B0p = %12.8f +/- %12.8f\n", params[3], errors[3]);
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

  // Write out Goncharovs residuals
  //  FILE *fin = fopen ("GoncharovEOSData.png.dat", "r");
  FILE *fin = fopen ("GoncharovEOS.Holzapfel.dat", "r");
  fout = fopen ("GoncharovResiduals.Holzapfel.dat", "w");
  double Vval, Pval;
  while (fscanf (fin, "%lf %lf\n", &Pval, &Vval) == 2) {
    // Vval *= 79.7535137356486;
    Vval *= 79.714027629;
    double diff = Pval - eos.Pressure(Vval);
    fprintf (fout, "%1.12e %1.12e\n", Vval, diff);
  }
  fclose(fout);
  fclose(fin);

  fin = fopen ("GoncharovEOS.png.dat", "r");
  fout = fopen ("GoncharovResiduals.dat", "w");
  while (fscanf (fin, "%lf %lf\n", &Pval, &Vval) == 2) {
    // Vval *= 79.7535137356486;
    Vval *= 79.714027629;
    double diff = Pval - eos.Pressure(Vval);
    fprintf (fout, "%1.12e %1.12e\n", Vval, diff);
  }
  fclose(fout);


}



void
TestFitErrors(string fname)
{
  TinyVector<double,4> params(80.0, 0.51, 0.0155), errors; 
  VinetEOSClass eos;
  eos.SetParams (params);
  NonlinearFitClass<4,VinetEOSClass> fitter(eos);
  Array<double,1> V, E, sigma;


  IOSectionClass in;
  assert(in.OpenFile (fname));
  assert(in.ReadVar("V0",  params[0]));
  assert(in.ReadVar("E0",  params[1]));
  assert(in.ReadVar("B0",  params[2]));
  assert(in.ReadVar("B0p", params[3]));
  params[2] /= 29421.01;
	
  assert(in.ReadVar("V", V));
  assert(in.ReadVar("E", E));
  assert(in.ReadVar("Sigma", sigma));
  eos.SetParams (params);
  fitter.Fit (V, E, sigma, params);

  int numSamples = 5000;
  int numVol     = 100;
  Array<double,1> Epert(E.size());
  CommunicatorClass local;
  RandomClass random(local);
  random.Init();
  Array<double,2> Pressures(numVol,numSamples);
  for (int i=0; i<numSamples; i++) {
    TinyVector<double,4> params2;
    params2 = params;
    VinetEOSClass eos2;
    for (int j=0; j<E.size(); j++) 
      Epert(j) = E(j) + random.LocalGaussian(sigma(j));
    eos2.SetParams (params2);
    fitter.Fit (V, Epert, sigma, params2);
    eos2.SetParams(params2);
    for (int vi=0; vi<numVol; vi++) {
      double v0 = V(0);
      double v1 = V(V.size()-1);
      double v = v0 + (v1-v0)*(double)(vi)/(double)(numVol-1);
      Pressures(vi, i) = eos2.Pressure (v);
    }
  }

  string errorName = fname + ".errors.dat";
  FILE *ferr = fopen (errorName.c_str(), "w");
  eos.SetParams (params);
  // Now compute error bars as a function of volume
  for (int vi=0; vi<numVol; vi++) {
    double v0 = V(0);
    double v1 = V(V.size()-1);
    double v = v0 + (v1-v0)*(double)vi/(double)(numVol-1);
    double Psum=0.0, P2sum=0.0;
    for (int j=0; j<numSamples; j++) {
      Psum += Pressures(vi,j);
      P2sum += Pressures(vi,j)*Pressures(vi,j);
    }
    double Pmean = Psum / (double)numSamples;
    double Pvar  = P2sum/ (double)numSamples - Pmean*Pmean;
    double Perr  = sqrt (Pvar);
    fprintf (ferr, "%1.8f %1.8f %1.12f %1.8f\n",
	     v, eos.Pressure(v), Pmean, Perr);
  }
      

  Array<double,2> &covar = fitter.GetCovariance();
  errors[0] = sqrt (covar(0,0));
  errors[1] = sqrt (covar(1,1));
  errors[2] = sqrt (covar(2,2));
  errors[3] = sqrt (covar(3,3));
  fprintf (stdout, "V0  = %12.8f +/- %12.8f\n", params[0], errors[0]);
  fprintf (stdout, "E0  = %12.8f +/- %12.8f\n", params[1], errors[1]);
  fprintf (stdout, "B0  = %12.8f +/- %12.8f\n", 
	   params[2]* 29421.01, errors[2]* 29421.01);
  eos.SetParams(params);
  fprintf (stdout, "B0p = %12.8f +/- %12.8f\n", params[3], errors[3]);
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
  double Vmin=V(0), Vmax=V(0), delta;
  for (int i=0; i<V.size(); i++) {
    Vmin = min(Vmin, V(i));
    Vmax = max(Vmax, V(i));
  }
  delta = (Vmax - Vmin) / 1000.0;
  for (double v=Vmin; v<=Vmax; v+=delta)
    fprintf (fout, "%20.16f %20.16f %20.16f %20.16f\n", v, eos(v), eos.Pressure(v), eos.PressureFD(v));
  fclose (fout);
}




void 
StaticFit (string fname, VinetEOSClass &eos)
{
  TinyVector<double,4> params(80.0, 0.51, 0.0155, 3.7), errors;  
  eos.SetParams (params);
  NonlinearFitClass<4,VinetEOSClass> fitter(eos);
  Array<double,1> V, E, sigma;

  IOSectionClass in;
  assert(in.OpenFile (fname));
  assert(in.ReadVar("V0",  params[0]));
  assert(in.ReadVar("E0",  params[1]));
  assert(in.ReadVar("B0",  params[2]));
  assert(in.ReadVar("B0p", params[3]));
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
  fprintf (stdout, "B0p = %12.8f\n", params[3]);
  fitter.Fit (V, E, sigma, params);

  Array<double,2> &covar = fitter.GetCovariance();
  errors[0] = sqrt (covar(0,0));
  errors[1] = sqrt (covar(1,1));
  errors[2] = sqrt (covar(2,2));
  errors[3] = sqrt (covar(3,3));
  fprintf (stdout, "V0  = %12.8f +/- %12.8f\n", params[0], errors[0]);
  fprintf (stdout, "E0  = %12.8f +/- %12.8f\n", params[1], errors[1]);
  fprintf (stdout, "B0  = %12.8f +/- %12.8f\n", 
	   params[2]* 29421.01, errors[2]* 29421.01);
  eos.SetParams(params);
  fprintf (stdout, "B0p = %12.8f +/- %12.8f\n", params[3], errors[3]);
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
  double Vmin=V(0), Vmax=V(0), delta;
  for (int i=0; i<V.size(); i++) {
    Vmin = min(Vmin, V(i));
    Vmax = max(Vmax, V(i));
  }
  cerr << "Vmin = " << Vmin << "  Vmax = " << Vmax << endl;

  delta = (Vmax - Vmin) / 1000.0;
  for (double v=Vmin; v<=Vmax; v+=delta)
    fprintf (fout, "%20.16f %20.16f %20.16f %20.16f\n", v, eos(v), eos.Pressure(v), eos.PressureFD(v));
  fclose (fout);

  // Write out Goncharovs residuals
  //  FILE *fin = fopen ("GoncharovEOSData.png.dat", "r");
  FILE *fin = fopen ("GoncharovEOS.Holzapfel.dat", "r");
  double Vval, Pval, aval;
  if (fin) {
    fout = fopen ("GoncharovResiduals.Holzapfel.dat", "w");
    while (fscanf (fin, "%lf %lf\n", &Pval, &aval) == 2) {
      //      Vval *= 79.7535137356486;
      aval *= 1.8897261;  // Angstrom to bohr radius
      Vval = 0.25 * aval * aval * aval;
      double diff = Pval - eos.Pressure(Vval);
      fprintf (fout, "%1.12e %1.12e\n", Vval, diff);
    }
    fclose(fout);
    fclose(fin);
  }


  fin = fopen ("GoncharovEOS.png.dat", "r");
  if (fin) {
    fout = fopen ("GoncharovResiduals.dat", "w");
    while (fscanf (fin, "%lf %lf\n", &Pval, &aval) == 2) {
      //      Vval *= 79.7535137356486;
      aval *= 1.8897261;  // Angstrom to bohr radius
      Vval = 0.25 * aval * aval * aval;
      double diff = Pval - eos.Pressure(Vval);
      fprintf (fout, "%1.12e %1.12e\n", Vval, diff);
    }
    fclose(fout);
    fclose(fin);
  }


  // Write out Datchi's residuals
  fin = fopen ("DatchiEOS.dat", "r");
  if (fin) {
    fout = fopen ("DatchiResiduals.dat", "w");
    double aVal, Tval;
    char dummyLine[100];
    fgets (dummyLine, 100, fin);
    fgets (dummyLine, 100, fin);
    while (fscanf (fin, "%lf %lf %lf %lf\n", &Pval, &Tval, &aVal, &Vval) == 4) {
      Vval = 1.8897261*1.8897261*1.8897261*aVal*aVal*aVal/4.0;
      double diff = Pval - eos.Pressure(Vval);
      fprintf (fout, "%1.12e %1.12e\n", Vval, diff);
    }
    fclose(fout);
    fclose (fin);
  }

  // Write out Knittles's residuals
  fin = fopen ("Exp300K_P-V.dat", "r");
  if (fin) {
    fout = fopen ("KnittleResiduals.dat", "w");
    while (fscanf (fin, "%lf %lf\n", &Vval, &Pval) == 2) {
      Vval *= 79.714027629;
      // Vval *= 79.7535137356486;
      double diff = Pval - eos.Pressure(Vval);
      fprintf (fout, "%1.12e %1.12e\n", Vval, diff);
    }
    fclose(fout);
    fclose (fin);
  }
}


void
ThermalFit (string fname, PhononFreeEnergy &phonons)
{
  FILE *fin = fopen (fname.c_str(), "r");
  std::vector<double> Vvec, Tvec, Fvec, Uvec;
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
  fclose (fin);
  int N = Vvec.size();
  Array<double,1> F(N), U(N), V(N), T(N);
  for (int i=0; i<N; i++) {
    V(i) = Vvec[i];
    T(i) = Tvec[i];
    F(i) = Fvec[i];
    U(i) = Uvec[i];
  }
  
  phonons.FitF(4, 6, F, V, T);
  
  FILE *fout = fopen ("300K_FofV.dat", "w");
  double Tval = 300.0;
  for (double Vval=40.0; Vval<95.0; Vval+=0.01) {
    double Fval = phonons.F_VT(Vval,Tval);
    fprintf (fout, "%1.8f %1.8e\n", Vval, Fval);
  }
  fclose(fout); 

  fout = fopen ("V0_FofT.dat", "w");
  double Vval = 79.7469;
  for (double Tval=0.0; Tval<3000.0; Tval+=1.0) {
   double Fval = phonons.F_VT(Vval,Tval);
    fprintf (fout, "%1.8f %1.8e\n", Tval, Fval);
  }
  fclose(fout);

}

void
DebyeFit (string fname, vector<DebyeModel> &fits)
{
  FILE *fin = fopen (fname.c_str(), "r");
  std::vector<double> Vvec, Tvec, Fvec, Uvec;
  double lastv = 0;
  bool done=false;
  int i = 0;
  while (!done) {
    double v,t,f,u;
    int retval = fscanf (fin, "%lf %lf %lf %lf", &v, &t, &f, &u);
    if (v != lastv || retval != 4) {
      if (Vvec.size() != 0) {
	DebyeModel debye;
	debye.SetN (2);
	Array<double,1> T(Tvec.size()), F(Fvec.size());
	for (int i=0; i<Tvec.size(); i++) {
	  F(i) = Fvec[i]-Fvec[0];  T(i) = Tvec[i];
	}
	double theta = debye.OptTheta_F (F, T);
	debye.SetTheta (theta);
	fprintf (stderr, "%8.5f %10.5f %1.8e\n",
		 cbrt(4.0*lastv), theta, Fvec[0]);
	
	fits.push_back(debye);
	char name[1000];
	snprintf (name, 1000, "Phonon_%1.2f.dat",
		  lastv);
	FILE *fout = fopen (name, "w");
	for (int i=0; i<Tvec.size(); i++) 
	  fprintf (fout, "%1.3f %1.9e %1.9e\n",
		   Tvec[i], Fvec[i]-Fvec[0], debye.F(Tvec[i]));
	fclose (fout);
      }
      Vvec.clear();
      Tvec.clear();
      Fvec.clear();
      Uvec.clear();
    }
    if (retval == 4) {
      Vvec.push_back (v);
      Tvec.push_back (t);
      Fvec.push_back (f);
      Uvec.push_back (u);
      lastv = v;
    }
    else
      done = true;
  }
  


}


void
DebyeFit (string fname, DebyeFreeEnergy &debye)
{
  FILE *fin = fopen (fname.c_str(), "r");
  std::vector<double> Vvec, Tvec, Fvec, Uvec;
  double lastv = 0;
  bool done=false;
  int i = 0;
  while (!done) {
    double v,t,f,u;
    int retval = fscanf (fin, "%lf %lf %lf %lf", &v, &t, &f, &u);
    if (v != lastv || retval != 4) {
      if (Vvec.size() != 0) 
	debye.AddVolume_F (lastv, Fvec, Tvec);
      //debye.AddVolume_U (lastv, Uvec, Tvec);
      Vvec.clear();
      Tvec.clear();
      Fvec.clear();
      Uvec.clear();
    }
    if (retval == 4) {
      Vvec.push_back (v);      Tvec.push_back (t);
      Fvec.push_back (f);      Uvec.push_back (u);
      lastv = v;
    }
    else
      done = true;
  }
  debye.FitTheta_V (4);
  fclose (fin);
  fin = fopen (fname.c_str(), "r");
  lastv = 0;
  done=false;
  i = 0;
  double F0 = 0.0;
  double U0 = 0.0;
  FILE *out = fopen ("ModelF.dat", "w");
  while (!done) {
    double v,t,f,u;
    int retval = fscanf (fin, "%lf %lf %lf %lf", &v, &t, &f, &u);
    if (retval != 4)
      done = true;
    else {
      if (v != lastv)
	F0 = f; U0 = u;
      fprintf (out, "%1.9f  %5.1f   %12.14f  %12.14f\n",
	       v, t, F0 + debye.F(v,t), U0 + debye.U(v,t));
      lastv = v;
    }
  }
  fclose(out);
  fclose (fin);



  FILE *fout = fopen ("Theta_Fdat", "w");
  for (double V=40.0; V<=95.0; V+=0.01) {
    double theta = debye.Theta_V(V);
    fprintf (fout, "%1.2f %1.10e\n", V, theta);
  }
  fclose (fout);
}


double FindV (double P, double T,
	      VinetEOSClass &staticEOS,
	      DebyeFreeEnergy &thermalEOS)
{
  double Vlo = 10.0;
  double Vhi = 150.0;
  double au2GPa = 29421.01;

  while ((Vhi - Vlo) > 1.0e-10) {
    double Vtry = 0.5*(Vlo + Vhi);
    double Ptry = staticEOS.Pressure(Vtry) +
      au2GPa * thermalEOS.P(Vtry,T);
    if (Ptry > P)
      Vlo = Vtry;
    else
      Vhi = Vtry;
  }
  return 0.5*(Vlo+Vhi);
}

void
CalcProperties (VinetEOSClass &staticEOS,
		DebyeFreeEnergy &thermalEOS)
{
  // Compute Thermal Pressure
  const double au2GPa = 29421.01;
  
  FILE *fout = fopen ("ThermalPressure.dat", "w");
  double V = 79.701;
  //double V = 60.0;
  for (double T=1.0; T<3000.0; T+=1.0) 
    fprintf (fout, "%5.1f %12.8e %12.8e\n", T, 
	     au2GPa * thermalEOS.P(V,T),
	     au2GPa * thermalEOS.P_FD(V,T));
  fclose (fout);

  // Test dF_dT
  fout = fopen ("dF_dT.dat", "w");
  for (double T=1.0; T<3000.0; T+=1.0) {
    fprintf (fout, "%5.1f %12.8e %12.8e\n", T,
	     thermalEOS.dF_dT (V,T), thermalEOS.dF_dT_FD(V,T));
  }
  fclose (fout);

  // Test d2F_dTheta_dT
  fout = fopen ("d2F_dTheta_dT.dat", "w");
  for (double T=1.0; T<3000.0; T+=1.0) {
    fprintf (fout, "%5.1f %12.8e %12.8e\n", T,
	     thermalEOS.d2F_dTheta_dT (V,T), 
	     thermalEOS.d2F_dTheta_dT_FD(V,T));
  }
  fclose (fout);

  // Test thermal bulk modulus
  fout = fopen ("K_T.dat", "w");
  for (double T=1.0; T<3000.0; T+=1.0) {
    fprintf (fout, "%5.1f %12.8e %12.8e\n", T,
	     au2GPa*thermalEOS.K_T(V,T), 
	     au2GPa*thermalEOS.K_T_FD(V,T));
  }
  fclose (fout);

  // Test dP_dT
  fout = fopen ("dP_dT.dat", "w");
  for (double T=1.0; T<3000.0; T+=1.0) {
    fprintf (fout, "%5.1f %12.8e %12.8e\n", T,
	     au2GPa*thermalEOS.dP_dT(V,T), 
	     au2GPa*thermalEOS.dP_dT_FD(V,T));
  }
  fclose (fout);

  // Thermal expansion coefficient
  fout = fopen ("alpha.dat", "w");
  for (double T=1.0; T<3000.0; T+=1.0) {
    double K_T = au2GPa*thermalEOS.K_T(V,T) 
      + staticEOS.K_T(V);
    double dP_dT = au2GPa * thermalEOS.dP_dT(V,T);
    double alpha = dP_dT/K_T;
    double C_V = thermalEOS.C_V (V, T);
    fprintf (fout, "%5.1f %12.8e %12.8f %12.8f \n", T, alpha, K_T);
  }
  fclose (fout);
  // Gruineissen parameter
  fout = fopen ("gamma.dat", "w");
  for (double T=1.0; T<3000.0; T+=1.0) {
    double K_T = thermalEOS.K_T(V,T) 
      + staticEOS.K_T(V)/au2GPa;
    double dP_dT = thermalEOS.dP_dT(V,T);
    double alpha = dP_dT/K_T;
    double C_V = thermalEOS.C_V(V,T);
    double gamma = V * dP_dT / C_V;
//     double theta = thermalEOS.Theta_V(V);
//     double gamma = -(V/theta)*thermalEOS.dTheta_dV(V);
    fprintf (fout, "%5.1f %12.8e\n", T, gamma);
  }
  fclose (fout);

  // Heat capacity
  fout = fopen ("C_V_from_F.dat", "w");
  for (double T=1.0; T<3000.0; T+=1.0) {
    double C_V = thermalEOS.C_V(V,T);
    fprintf (fout, "%5.1f %12.8e\n", T, C_V);
  }
  fclose (fout);

  // TO phonon mode as a function of P,T
  fout = fopen ("0P_Raman.dat", "w");
  RamanFreq raman;
  for (double T=300.0; T<1250.0; T+=1.0) {
    double V = FindV (0.0, T, staticEOS, thermalEOS);
    double freq = raman(V);
    fprintf (fout, "%5.1f %12.8f\n", T, freq);
  }
  fclose (fout);

  // TO phonon mode as a function of P,T
  vector<int> Tlist;
  Tlist.push_back(300);
  Tlist.push_back(373);
  Tlist.push_back(473);
  Tlist.push_back(573);
  Tlist.push_back(673);
  Tlist.push_back(723);
  for (int ti=0; ti<Tlist.size(); ti++) {
    char fname[100];
    snprintf (fname, 100, "%dK_Raman.dat", Tlist[ti]);
    double T = (double)Tlist[ti];
    double P = 0.0;
    fout = fopen (fname, "w");
    double deltaV = FindV(0.0, 300, staticEOS, thermalEOS) - 
      (6.8317379*6.8317379*6.8317379)/4.0;
    for (double P=0.0; P<800.0; P+=1.0) {
      double V = FindV (P, T, staticEOS, thermalEOS);
      double freq = raman(V);
      fprintf (fout, "%5.1f %12.8f\n", P, freq);
    }
    fclose (fout);
  }
}



main(int argc, char **argv)
{
  TestGrad();
  //TestVal();
  if (argc > 2) {
    VinetEOSClass staticEOS;
    StaticFit (argv[1], staticEOS);
    PhononFreeEnergy thermalEOS;
    ThermalFit (argv[2], thermalEOS);
    //vector<DebyeModel> models;
    //DebyeFit (argv[2], models);
    DebyeFreeEnergy debyeEOS;
    DebyeFit (argv[2], debyeEOS);
    CalcProperties (staticEOS, debyeEOS);
  }
  else if (argc > 1) {
    VinetEOSClass staticEOS;
    StaticFit (argv[1], staticEOS);
    TestFitErrors(argv[1]);
  }

  else
    cerr << "Usage:  VinetFit2 static.in phonon.dat\n";
}
