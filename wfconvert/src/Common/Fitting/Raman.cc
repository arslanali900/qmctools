#include "Raman.h"
#include "../Fitting/Fitting.h"
#include "../Random/Random.h"
#include <vector>

double
VFitClass::GetGridStart()
{
  return xstart;
}

double
VFitClass::GetGridEnd()
{
  return xend;
}

double
VFitClass::operator()(double x)
{
  if (UseHarmonic) 
    return Vcoefs[0] + x*x*Vcoefs[1];
  else {
    // if (x > xend) 
    //   return Vend + (x-xend)*dVend + 0.5*(x-xend)*(x-xend)*d2Vend;
    // else if (x < xstart)
    //   return Vstart + (x-xstart)*dVstart + 0.5*(x-xstart)*(x-xstart)*d2Vstart;
    // else
      return Vcoefs[0] + x*x*Vcoefs[1] + x*x*x*Vcoefs[2] + x*x*x*x*Vcoefs[3];
  }
}

double
VFitClass::Deriv(double x)
{
  return 2.0*x*Vcoefs[1] + 3.0*x*x*Vcoefs[2] + 4.0*x*x*x*Vcoefs[3];
}

double
VFitClass::Deriv2(double x)
{
  return 2.0*Vcoefs[1] + 6.0*x*Vcoefs[2] + 12.0*x*x*Vcoefs[3];
}

void
VFitClass::DoFit (Array<double,1> &Etry)
{
  int N = Grid.size();
  Array<double,1> coefs(4), dcoefs(4);
  Array<double,2> F(N,4);
  for (int i=0; i<N; i++) {
    double x = Grid(i);
    F(i,0) = 1.0;
    F(i,1) = x*x;
    F(i,2) = x*x*x;
    F(i,3) = x*x*x*x;
  }
  LinFitSVD (Etry, Sigma, F, coefs, dcoefs, 1.0e-14);
  for (int i=0; i<4; i++)
    Vcoefs[i] = coefs(i);
}

void
VFitClass::Read (IOSectionClass &in)
{
  Array<double,1> coefs, dcoefs;
  assert (in.ReadVar("GridPoints", Grid));
  assert (in.ReadVar("E", E));
  assert (in.ReadVar("error", Sigma));
  int N = Grid.size();
  assert (E.size() == N);
  assert (Sigma.size() == N);
  // coefs.resize(4);
  // dcoefs.resize(4);
  // for (int i=0; i<N; i++) {
  //   double x = Grid(i);
  //   F(i,0) = 1.0;
  //   F(i,1) = x*x;
  //   F(i,2) = x*x*x;
  //   F(i,3) = x*x*x*x;
  // }
  // LinFitSVD (E, Sigma, F, coefs, dcoefs, 1.0e-14);
  // for (int i=0; i<4; i++)
  //   Vcoefs[i] = coefs(i);
  DoFit (E);
  cerr << "coefs = " << Vcoefs << endl;

  xstart   = 1.0*Grid(0);
  xend     = 1.0*Grid(N-1);
  Vstart   = (*this)(xstart+1.0e-12);
  dVstart  = Deriv(xstart);
  // HACK HACK HACK
  //d2Vstart = fabs(Deriv2(xstart));
  d2Vstart = Deriv2(xstart);
  Vend     = (*this)(xend-1.0e-12);
  dVend    = Deriv(xend);
  d2Vend   = Deriv2(xend);

  fname = in.GetFileName();
  FILE *fout = fopen ((fname+".fit.dat").c_str(), "w");
  for (double x=2.2*xstart; x<=2.2*xend; x+=0.001)
    fprintf (fout, "%1.12f %1.12f\n", x, (*this)(x));
  fclose (fout);
}


double
VSplineClass::GetGridStart()
{
  return SplineGrid.Start;
}

double
VSplineClass::GetGridEnd()
{
  return SplineGrid.End;
}


double
VSplineClass::operator()(double x)
{
  if (UseHarmonic) 
    return Vmin + 0.5*d2Vmin*(x-xmin)*(x-xmin);
  else if (x < SplineGrid.Start) {
    double dx = x - SplineGrid.Start;
    return Vstart + dVstart * dx + 0.5*d2Vstart*dx*dx;
  }
  else if (x >= SplineGrid.End) {
    double dx = x - SplineGrid.End;
    return Vend + dVend*dx + 0.5*d2Vend*dx*dx;
  }
  else
    return Vspline(x);
}

void
VSplineClass::Read(IOSectionClass &in)
{
  // Read potential data and initialize spline
  Array<double,1> gridPoints, Edata;
  assert (in.ReadVar("GridPoints", gridPoints));
  assert (in.ReadVar("E"         , Edata));
  assert (Edata.size() == gridPoints.size());
  double d2E_left, d2E_right;
  d2E_left = ((Edata(2)-Edata(1))/(gridPoints(2)-gridPoints(1)) -
	      (Edata(1)-Edata(0))/(gridPoints(1)-gridPoints(0))) /
    (gridPoints(1) - gridPoints(0));
  int n = gridPoints.size();
  d2E_right = ((Edata(n-1)-Edata(n-2))/(gridPoints(n-1)-gridPoints(n-2)) -
	       (Edata(n-2)-Edata(n-3))/(gridPoints(n-2)-gridPoints(n-3)))/
    (gridPoints(n-1) - gridPoints(n-2));
  BoundaryCondition<double> lBC(FIXED_SECOND, d2E_left);
  BoundaryCondition<double> rBC(FIXED_SECOND, d2E_right);
  SplineGrid.Init (gridPoints);
  // Initialize potential spline
  Vspline.Init (SplineGrid, Edata, lBC, rBC);

  Vstart   = Vspline(SplineGrid.Start);
  dVstart  = Vspline.Deriv(SplineGrid.Start);
  d2Vstart = Vspline.Deriv2(SplineGrid.Start);

  Vend   = Vspline(SplineGrid.End);
  dVend  = Vspline.Deriv(SplineGrid.End);
  d2Vend = Vspline.Deriv2(SplineGrid.End);

  xmin = 0.0;
  Vmin = 1.0e200;
  for (int i=0; i<10000; i++) {
    double t = (double)i/9999.0;
    double x = (1.0-t)*SplineGrid.Start + t*SplineGrid.End;
    if (Vspline(x) < Vmin) {
      Vmin = Vspline(x);
      xmin = x;
    }
  }
  d2Vmin = Vspline.Deriv2(xmin);
}


double
PhononClass::Integrate()
{
  RungeKutta<PhononClass, Vec2> integrator(*this);

  int N = IntegrationGrid.NumPoints;
  int Nmatch = N/3;
  if (u_du.size() != N)
    u_du.resize(N);

  u_du(0)   = Vec2 (1.0e-6, 1.0);
  u_du(N-1) = Vec2 (1.0e-6, -1.0);
  integrator.Integrate (IntegrationGrid, 0, Nmatch, u_du);
  Vec2 umid = u_du(Nmatch);
  integrator.Integrate (IntegrationGrid, N-1, Nmatch, u_du);

  // Scale right half to match left half
  double factor = umid[0]/u_du(Nmatch)[0];
  for (int i=Nmatch; i<N; i++)
    u_du(i) *= factor;

  // Now, normalize
  double nrm = 0.0;
  for (int i=0; i<N; i++)
    nrm += u_du(i)[0]*u_du(i)[0];
  nrm = 1.0/sqrt(nrm);
  for (int i=0; i<N; i++)
    u_du(i) *= nrm;
  
  double log_deriv_left = umid[1]/umid[0];
  double log_deriv_right = u_du(Nmatch)[1] / u_du(Nmatch)[0];
  
  return (u_du(Nmatch)[1]/u_du(Nmatch)[0] - umid[1]/umid[0]);
}

void
PhononClass::Write()
{
  FILE *fout = fopen ("u.dat", "w");
  int N = IntegrationGrid.NumPoints;
  bool useHarm = Vfunction->UseHarmonic;
  for (int i=0; i<N; i++) {
    double x = IntegrationGrid(i);
    Vfunction->UseHarmonic = true;
    double Vharm = V(x);
    Vfunction->UseHarmonic = false;
    double Vanharm = V(x);
    fprintf (fout, "%1.16e %1.16e %1.16e %1.16e %1.16e\n",
	     x, u_du(i)[0], u_du(i)[1], Vharm, Vanharm);
  }
  fclose (fout);
  Vfunction->UseHarmonic = useHarm;
}  

void
PhononClass::Read(IOSectionClass &in)
{
  Array<double,1> dummy;
  UseFit = in.ReadVar("error", dummy);
  if (UseFit) 
    Vfunction = new VFitClass;
  else
    Vfunction = new VSplineClass;
  Vfunction->Read(in);

  assert (in.ReadVar("ReducedMass", ReducedMass));
  lambda = 0.5/ReducedMass;
  lambdaInv = 1.0/lambda;
  int numPoints = 2000;
  in.ReadVar ("IntegrationPoints", numPoints);
  u_du.resize(numPoints);
  
  double gridStart = Vfunction->GetGridStart();
  double gridEnd   = Vfunction->GetGridEnd();

  double center = 0.5*(gridStart + gridEnd);
  double width = 1.10*(gridEnd - gridStart);
  double start = center - width;
  double end   = center + width;

  IntegrationGrid.Init (start, end, numPoints);
  if (UseFit) {
    const double h = 1.5198298e-16;
    const double c = 2.9979246e+10;
    VFitClass &vfit (*dynamic_cast<VFitClass*>(Vfunction));
    int N = vfit.E.size();
    Array<double,1> Etry(N);
    const int numFits = 100;
    double fsum = 0.0, f2sum = 0.0;
    for (int i=0; i<numFits; i++) {
      for (int j=0; j<N; j++)
	Etry(j) = vfit.E(j) + Random.LocalGaussian(vfit.Sigma(j));
      vfit.DoFit (Etry);
      double E0 = Solve(0);
      double E1 = Solve(1);
      double f = (E1-E0)/(h*c);
      //      cerr << "frequency = " << f << endl;
      fsum += f;
      f2sum += f*f;
    }
    double sigma_f = 
      sqrt(f2sum/(double)numFits - fsum*fsum/(double)(numFits*numFits));
    vfit.DoFit (vfit.E);
    double E0 = Solve(0);
    double E1 = Solve(1);
    double f = (E1-E0)/(h*c);
    fprintf (stderr, "Lowest Raman freq. = %1.3f +/- %1.3f\n",
	     f, sigma_f);
  }
  

  Vfunction->UseHarmonic = false;
}

int
PhononClass::CountNodes()
{
  int numNodes = 0;
  for (int i=0; i<u_du.size()-1; i++) 
    if (u_du(i)[0]*u_du(i+1)[0] < 0.0) {
      numNodes ++;
      //cerr << "Node found at i = " << i << endl;
    }
  return numNodes;
}

double
PhononClass::Solve(int desiredNodes)
{
  double Ehi, Elo;
  Ehi = Elo = V(IntegrationGrid(0));
  for (int i=1; i<IntegrationGrid.NumPoints; i++) {
    Elo = min (V(IntegrationGrid(i)), Elo);
    Ehi = max (V(IntegrationGrid(i)), Ehi);
  }  

  while ((Ehi - Elo) > 1.0e-14) {
    Etrial = 0.5*(Ehi + Elo);
    double cusp = Integrate();
    int numNodes = CountNodes();
    if (numNodes > desiredNodes)
      Ehi = Etrial;
    else if (numNodes < desiredNodes)
      Elo = Etrial;
    else if (cusp < 0.0)
      Elo = Etrial;
    else
      Ehi = Etrial;
  }
  //if (desiredNodes==0) Write();
  return 0.5 * (Ehi + Elo);
}

TinyVector<double,2>
PhononClass::SolveError(int desiredNodes)
{
  double E, sigma = 0.0;
  if (UseFit) {
    VFitClass &vfit (*dynamic_cast<VFitClass*>(Vfunction));
    int N = vfit.E.size();
    Array<double,1> Etry(N);
    const int numFits = 100;
    double esum = 0.0, e2sum = 0.0;
    for (int i=0; i<numFits; i++) {
      for (int j=0; j<N; j++)
	Etry(j) = vfit.E(j) + Random.LocalGaussian(vfit.Sigma(j));
      vfit.DoFit (Etry);
      double en = Solve(desiredNodes);
      double en1 = Solve(desiredNodes+1);
      double e = en1-en;
      esum += e;
      e2sum += e*e;
    }
    sigma = 
      sqrt(e2sum/(double)numFits - esum*esum/(double)(numFits*numFits));
    vfit.DoFit (vfit.E);
  }
  E = Solve(desiredNodes);
  ostringstream outname;
  outname << Vfunction->fname << "_" << desiredNodes << ".dat";
  FILE *fout = fopen (outname.str().c_str(), "w");
  for (int i=0; i<IntegrationGrid.NumPoints; i++)
    fprintf (fout, "%1.8e %1.8e %1.8e %1.8e\n",
	     IntegrationGrid(i), u_du(i)[0], u_du(i)[1], E);
  fclose(fout);


  return TinyVector<double,2> (E, sigma);
}


void
RamanModel::AddVolume(string fname)
{
  IOSectionClass in;
  assert (in.OpenFile(fname));
  PhononClass phonon (Random);
  double V;
  in.ReadVar("V", V);
  Volume.push_back(V);
  fprintf (stderr, "Adding volume %1.6f to phonon database.\n", V);
  phonon.Read(in);

  string freqName = fname+".freq.dat";
  FILE *fout = fopen (freqName.c_str(), "w");

  Array<double,1> E(NumEnergies), sigma(NumEnergies);
  TinyVector<double,2> E_sigma;
  for (int i=0; i<NumEnergies; i++) {
    E_sigma = phonon.SolveError(i);
    E(i)     = E_sigma[0];
    sigma(i) = E_sigma[1];
  }
  for (int i=1; i<NumEnergies; i++) 
    fprintf (fout, "%12.5f %12.5f %12.5f\n",
	     (E(i)-E(0))/kB, (E(i) - E(i-1))/(h*c), 
	     sigma(i)/(h*c));
  fclose(fout);

  Energies.push_back(E);
  Sigma.push_back(sigma);
}

void
RamanModel::Read(string fname)
{
  if (fname.find(".files") < fname.size()) {
    cerr << "\nReading data for Raman model:\n";
    FILE *fin = fopen (fname.c_str(), "r");
    assert (fin != NULL);
    char name[1000];
    while (fgets(name, 1000, fin) != NULL) {
      string filename = name;
      int N = filename.size();
      if (filename[N-1] == '\n')
	filename.erase(N-1, 1);
      AddVolume (filename);  
    }
    DoFits();
  }
}

void
RamanModel::DoFits()
{
  EnergyCoefs.resize(NumEnergies, NumCoefs);
  E0Coefs.resize(NumCoefs);

  Array<double,2> basis(Volume.size(), NumCoefs);
  for (int i=0; i<Volume.size(); i++) {
    double v = Volume[i];
    basis(i,0) = 1.0;
    for (int j=1; j<NumCoefs; j++)
      basis(i,j) = (1.0/v)*basis(i,j-1);
  }

  Array<double,1> coefs(NumCoefs), errors(NumCoefs);
  Array<double,1> e(Volume.size()), sigma(Volume.size());
  for (int i=0; i<NumEnergies; i++) {
    for (int j=0; j<Volume.size(); j++) {
      if (i > 0)
	e(j) = Energies[j](i) - Energies[j](0);
      else
	e(j) = Energies[j](0);
      sigma(j) = Sigma[j](i);
    }
    LinFitSVD (e, sigma, basis, coefs, errors, 1.0e-14);
    for (int n=0; n<NumCoefs; n++) {
      if (i > 0)
	EnergyCoefs(i, n) = coefs(n);
      else {
	EnergyCoefs(i, 0) = 0.0;
	E0Coefs(n) = coefs(n);
      }
      fprintf (stderr, "%10.5e ", coefs(n));
    }
    fprintf (stderr, "\n");
  }
  
  IOSectionClass out;
  assert (out.NewFile ("RamanModel.h5"));
  out.WriteVar ("EnergyCoefs", EnergyCoefs);
  out.WriteVar ("E0Coefs", E0Coefs);
  out.CloseFile();

  FILE *fout = fopen ("ModelFrequencies.dat", "w");
  assert (fout != NULL);
  for (double V=40.0; V<95.0; V+=0.1) {
    fprintf (fout, "%5.1f ", V);
    for (int i=0; i<NumEnergies-1; i++) 
      fprintf (fout, "%12.6f ", Frequency(i, V));
    fprintf (fout, "\n");
  }
  fclose (fout);

}


double
RamanModel::MeanFrequency (double V, double T)
{
  vector<double> e(NumEnergies);
  for (int i=0; i<NumEnergies; i++)
    e[i] = E(i,V);
  Spectrum.SetEnergies(e);
  Spectrum.SetTemp(T);
  return Spectrum.MeanFrequency();
}




/*
main(int argc, char **argv)
{
  PhononClass phonon;
  RamanSpectrum spectrum;
  IOSectionClass in;
  if ((argc < 2) || !in.OpenFile (argv[1]))  {
    cerr << "Usage:  phonon myfile.in\n";
    exit(-1);
  }
  
  phonon.Read(in);
  vector<double> harmVec, anharmVec;
  for (int i=0; i<15; i++) {
    phonon.Vfunction->UseHarmonic = false;
    double E = phonon.Solve(i);
    anharmVec.push_back(E);
    spectrum.AddEnergy(E);
    phonon.Vfunction->UseHarmonic = true;
    harmVec.push_back(phonon.Solve(i));
  }
  const double h = 1.5198298e-16;
  const double c = 2.9979246e+10;
  const double kB = 3.1668152e-06;
  double T1 = 300.0;
  double T2 = 723.0;
  
  spectrum.SetTemp(T2);

  FILE *maxout = fopen ("maxloc.dat", "w");
  for (double width=1.0; width<25.0; width+=4.0) {
    spectrum.SetWidth(width);
    char fname[100];
    snprintf (fname, 100, "RamanSpectrum_%1.1f.dat", width);
    FILE *fout = fopen (fname, "w");
    double maxloc = 900.0;
    double maxval = 0.0;
    for (double nu= 900.0; nu<1060.0; nu += 0.01) {
      double s = spectrum.Spectrum(nu);
      if (s > maxval) {
	maxloc = nu;
	maxval = s;
      }
      fprintf (fout, "%1.8f %1.12e\n", nu, spectrum.Spectrum(nu));
    }
    fprintf (maxout, "%1.3f %1.3f %1.8f\n", width, maxloc, maxval);
    fclose (fout);
  }
  fclose(maxout);

  if (false) {
    // Create animation of temperature-dependent frequency shift.
    double Tstart =  300.0;
    double Tend   = 2000.0;
    double widthStart = 3.0;
    double widthEnd = 15.0;
    maxout = fopen ("maxAnim.dat", "w");
    for (int i=0; i<=200; i++) {
      double x = (double)i/200.0;
      double T = (1.0-x)*Tstart + x*Tend;
      double width = (1.0-x)*widthStart + x * widthEnd;
      spectrum.SetTemp(T);
      spectrum.SetWidth(width);
      char fname[100];
      snprintf (fname, 100, "AnimFrame%d.dat", i);
      FILE *fout = fopen (fname, "w");
      
      double maxloc = 900.0;
      double maxval = 0.0;    
      for (double nu= 900.0; nu<1060.0; nu += 0.05) {
	double s = spectrum.Spectrum(nu);
	if (s > maxval) {
	  maxloc = nu;
	  maxval = s;
	}
	fprintf (fout, "%1.8f %1.12e\n", nu, spectrum.Spectrum(nu));
      }
      fprintf (maxout, "%1.3f %1.3f %1.8f %1.8f\n", width, maxloc, maxval, T);
      fclose (fout);
    }
  }

  fprintf (stderr, "Raman frequencies (1/cm):\n");
  fprintf (stderr, " Excitation   Harmonic        Anharmonic   Relative intensity   Norm. intensity\n");
  for (int i=1; i<harmVec.size(); i++) {
    double harmFreq   = (  harmVec[i] -   harmVec[i-1])/(h*c);
    double anharmFreq = (anharmVec[i] - anharmVec[i-1])/(h*c);
    double R = (double)i * exp(-(anharmVec[i-1]-anharmVec[0])/(kB*T1));
    double I = spectrum.Intensity(i);
    fprintf (stderr, "     %2d   %12.3f     %12.3f          %6.4f             %6.4f\n",
	     i, harmFreq, anharmFreq, R, I);
  }
  double Z = 0.0, dZ_dT=0.0;
  cerr << "kBT = " << kB*T1 << endl;
  double N=0.0, dN_dT = 0.0;
  for (int i=0; i<anharmVec.size()-1; i++) {
    double dE = anharmVec[i] - anharmVec[0];
    double a = exp(-dE/(kB*T1));
    double anharmFreq = (anharmVec[i+1] - anharmVec[i])/(h*c);
    double n = (double)(i+1);
    //double n = 1.0;
    N     += n * a * anharmFreq;
    dN_dT += n * a * anharmFreq*dE/(kB*T1*T1);
    dZ_dT += n * a * dE/(kB*T1*T1);
    Z     += n * a;
  }
  double avgFreq = N / Z;
  double davgFreq_dT = dN_dT/Z - N/(Z*Z)*dZ_dT;
  
  N = Z = 0.0;
  for (int i=0; i<anharmVec.size()-1; i++) {
    double dE = anharmVec[i] - anharmVec[0];
    double a = exp(-dE/(kB*T2));
    double anharmFreq = (anharmVec[i+1] - anharmVec[i])/(h*c);
    double n = (double)(i+1);
    //double n = 1.0;
    N += n * a * anharmFreq;
    Z += n * a;
  }
  double avgFreq2 = N/Z;
  double davgFreq_dT2 = avgFreq2 - avgFreq;
  
  fprintf (stderr, "Average frequency at 300K = %12.5f\n", avgFreq);
  fprintf (stderr, "Derivative of avg freq = %12.5f\n", 
	   davgFreq_dT);
  fprintf (stderr, "Average frequency at 723K = %12.5f\n", avgFreq2);
}
*/

