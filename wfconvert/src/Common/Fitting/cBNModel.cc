#include "cBNModel.h"
#include "DatchiRamanModel.h"

double
cBNModel::FindV(double P, double T)
{
  double Vlo = 10.0;
  double Vhi = 150.0;
  double au2GPa = 29421.01;

  while ((Vhi - Vlo) > 1.0e-12) {
    double Vtry = 0.5*(Vlo + Vhi);
    double Ptry = Static.Pressure(Vtry) +
      Phonon.P(Vtry,T);
    if (Ptry > P)
      Vlo = Vtry;
    else
      Vhi = Vtry;
  }
  return 0.5*(Vlo+Vhi);
}



void
cBNModel::SetStatic (string fname)
{
  TinyVector<double,4> params(80.0, 0.51, 0.0155, 3.7), errors;  
  Static.SetParams (params);
  NonlinearFitClass<4,VinetEOSClass> fitter(Static);
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
  Static.SetParams (params);
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
  Static.SetParams(params);
  fprintf (stdout, "B0p = %12.8f +/- %12.8f\n", params[3], errors[3]);
  fprintf (stdout, "Delta = %1.5f\n", 0.0);
  double V0 = params[0];
  double latConst = cbrt(4.0*V0);
  double latConstError = cbrt(4.0)/(3.0*cbrt(V0)*cbrt(V0)) * errors[0];
  fprintf (stdout, "Lattice const. = %1.6f +/- %1.6f bohr radii\n",
	   latConst, latConstError);
  fprintf (stdout, "               = %1.6f +/- %1.6f angstrom\n",  
	   0.52917721*latConst, 0.52917721*latConstError);
}



void
cBNModel::SetPhonon (string fname)
{
  FILE *fin = fopen (fname.c_str(), "r");
  assert (fin != NULL);
  std::vector<double> F_T, U_T;
  std::vector<std::vector<double> > F_VT, U_VT;
  std::vector<double> V, T;
  bool done = false;

  double lastv = 0.0;
  bool first = true;

  double v,t,f,u;
  int numFit = 0;
  while (!done) {
    int retval = fscanf (fin, "%lf %lf %lf %lf", &v, &t, &f, &u);
    if (retval == 4) {
      if (v != lastv) {
	if (first) {
	  first = false;
	}
	else {
	  F_VT.push_back (F_T);
	  U_VT.push_back (U_T);
	  V.push_back(lastv);
	  F_T.clear();
	  U_T.clear();
	  T.clear();
	  numFit = 0;
	}
      }
      if (t <= 10000.0 && t >= 0.0) {
	T.push_back (t);
	F_T.push_back (f);
	U_T.push_back (u);
	if (t > 500.0)
	  numFit++;
      }
      lastv = v;
    }
    else {
      F_VT.push_back (F_T);
      U_VT.push_back (U_T);
      V.push_back(lastv);
      done = true;
    }
  }
  F_VT.push_back(F_T);
  U_VT.push_back(U_T);
  fclose (fin);
  
  int numV = V.size();
  int numT = T.size();
  fprintf (stderr, "numT = %d\n", numT);
  double startV = V[0];
  double endV = V[V.size()-1];
  LinearGrid *Vgrid = new LinearGrid (startV, endV, numV);
  double startT = T[0];
  double endT   = T[T.size()-1];
  LinearGrid *Tgrid = new LinearGrid (startT, endT, numT);

  Array<double,2> F (numV, numT);
  Array<double,2> U (numV, numT);
  for (int iV=0; iV<numV; iV++) {
    //fprintf (stderr, "V = %1.2f\n", (*Vgrid)(iV));
    for (int iT=0; iT<numT; iT++) {
      //    fprintf (stderr, "T = %1.2f\n", (*Tgrid)(iV));
      F(iV, iT) = F_VT[iV][iT];
      U(iV, iT) = U_VT[iV][iT];
    }
  }
  Phonon.SetF(Vgrid, Tgrid, F);
  Phonon.SetU(Vgrid, Tgrid, U);

  fin = fopen (fname.c_str(), "r");
  lastv = 0;
  done=false;
  double F0 = 0.0;
  double U0 = 0.0;
  FILE *out = fopen ("ModelF2.dat", "w");
  while (!done) {
    double v,t,f,u;
    int retval = fscanf (fin, "%lf %lf %lf %lf", &v, &t, &f, &u);
    if (retval != 4)
      done = true;
    else {
      if (v != lastv)
	F0 = f; U0 = u;
      fprintf (out, "%1.9f  %5.1f   %12.14f  %12.14f\n",
	       v, t, Phonon.F(v,t), Phonon.U(v,t));
      lastv = v;
    }
  }
  fclose(out);
  fclose (fin);

  out = fopen ("TestF.dat", "w");
  for (double V=40.0; V<95.0; V+=0.1) {
    for (double T=0.0; T<2500.0; T+= 3.32342341)
      fprintf (out, "%1.12e ", Phonon.F(V,T));
    fprintf (out, "\n");
  }
  fclose(out);


  // Setup Debye model
  Debye.SetN(2);
  FILE *dout = fopen ("DebyeModel.dat", "w");
  FILE *pout = fopen ("DebyeTempParams.dat", "w");
  for (int iV=0; iV<numV; iV++) {
    Array<double,1> F_T(numT), T2(numT);
    Array<double,1> fitT(numFit), fitTheta(numFit), fitSigma(numFit);

    for (int iT=0; iT<numT; iT++) {
      F_T(iT) = F_VT[iV][iT];
      T2(iT)   = T[iT];
    }
    double theta = Debye.OptTheta_F(F_T, T2);
    cerr << "theta = " << theta << endl;
    Debye.SetTheta(theta);

    DebyeTemp tempModel;
    NonlinearFitClass<3,DebyeTemp> fitter(tempModel);
    TinyVector<double,3> params(2000.0,0.001, -10.0);
    vector<TinyVector<double,3> > DebyeTempParams;

    int iFit = 0;
    for (int iT=0; iT<numT; iT++) {
      if (T[iT] > 500.0) {
	double th = Debye.CalcTheta(F_VT[iV][iT], T[iT]);
	fitT(iFit) = T[iT];
	fitTheta(iFit) = th;
	fitSigma(iFit) = 0.01;
	iFit++;
      }
    }
    if (V[iV] < 79.9) {
      fitter.Fit (fitT, fitTheta, fitSigma, params);
      fprintf (pout, "%10.5f %14.8e  %14.8e  %14.8e\n",
	       V[iV], params[0], params[1], params[2]);
      cerr << "params = " << params << endl;
      tempModel.SetParams(params);
      DebyeTempParams.push_back(params);
    }
    // // Now fit theta parameters vs. V
    // Array<double,1> Varray(DebyeTempParams.size());
    // Array<double,1> fitTheta, 



    for (int iT=0; iT<numT; iT++) {
      double th = Debye.CalcTheta(F_VT[iV][iT], T[iT]);
      //double th = 0.0;
      Debye.SetTheta(theta);
      double modelTheta = tempModel(T[iT]);
      
      double v = V[iV];
      // double theta0 = (4.510998e+03      + -5.723577e+01*v +
      // 		       3.103155e-01*v*v  + -6.082507e-04*v*v*v);
      // double alpha  = (3.344318e-03      + -1.533690e-04*v +
      // 		       3.712041e-06*v*v  + -2.682415e-08*v*v*v);
      // double beta   = ( -5.133928e+03    +  1.715907e+02*v +
      // 		       -1.937482e+00*v*v + 7.446629e-03*v*v*v);
      double theta0 = (4.836656e+03      + -6.929704e+01*v +
		       4.634278e-01*v*v  + -1.273468e-03*v*v*v);
      double alpha  = (2.598608e-03      + -1.200504e-04*v +
		       2.789128e-06*v*v  + -2.008994e-08*v*v*v);
      double beta   =(-4.869563e+03      + 1.763443e+02*v +
		      -2.186053e+00*v*v  + 9.303181e-03*v*v*v);


      double modelTheta2 = theta0 + beta*exp(-fabs(alpha)*T[iT]);

      fprintf (dout, "%12.5f %12.5f  %14.8e  %14.8e  %14.8e  %14.8e ",
	       V[iV], T[iT], Debye.F(T[iT])-Debye.F(0.0), 
	       Debye.U(T[iT])-Debye.U(0.0), th,
	       modelTheta2);
      Debye.SetTheta(modelTheta2);
      fprintf (dout, "%14.8e %14.8e\n", Debye.F(T[iT])-Debye.F(0.0),
	       Debye.U(T[iT])-Debye.U(0.0));
    }
  }
  fclose(dout);
  fclose(pout);
}

void
cBNModel::SetRaman(string fname)
{
  Raman.Read (fname);
}

double
cBNModel::MeanFrequency_PT(double P, double T)
{
  double V = FindV(P, T);
  return Raman.MeanFrequency (V,T);
}


void
ThermalExpansion (cBNModel &model)
{
  FILE *fout = fopen ("ThermalExpansion.dat", "w");
  double P = 0.0;
  fprintf (fout, "# T      V(P=0)\n");
  vector<double> Plist;
  Plist.push_back (0.0);
  Plist.push_back (10.0);
  Plist.push_back (20.0);
  Plist.push_back (50.0);
  Plist.push_back (60.0);
  Plist.push_back (80.0);
  Plist.push_back (100.0);

  for (double T=0.0; T< 1400.0; T+=1.0) {
    fprintf (fout, "%6.1f ", T);
    for (int ip=0; ip<Plist.size(); ip++) {
      double P = Plist[ip];
      double V = model.FindV(P,T);
      double Vplus = model.FindV(P,T+0.5);
      double Vminus = model.FindV(P, T-0.5);
      double alpha = (Vplus-Vminus)/V;
      fprintf (fout, "%16.12f  %16.12e ", V, alpha);
    }
    fprintf (fout, "\n");
  }
  fclose(fout);
}


void 
CalcProperties(cBNModel &model)
{
  FILE *fout = fopen ("MeanFrequency_PT.dat", "w");
  vector<double> Tlist;
  Tlist.push_back(300.0);
  Tlist.push_back(373.0);
  Tlist.push_back(473.0);
  Tlist.push_back(573.0);
  Tlist.push_back(673.0);
  Tlist.push_back(723.0);

  fprintf (fout, "# Pressure  nu(300K)    nu(373K)   nu(473K)   nu(573K)   nu(673K)   nu(723K)\n");
  for (double P=0; P < 900.0; P+=1.0) {
    fprintf (fout, "%12.1f ", P);
    for (int ti=0; ti<Tlist.size(); ti++)
      fprintf (fout, "%11.6f ", model.MeanFrequency_PT(P, Tlist[ti]));
    fprintf (fout, "\n");
  }
  fclose(fout);

  fout = fopen ("MeanFrequency_T.dat", "w");
  for (double T=0.0; T<10000.0; T+= 1.0) {
    fprintf (fout, "%12.1f ", T);
    for (double P=0; P<=900.0; P+=20.0) 
      fprintf (fout, "%11.6f ", model.MeanFrequency_PT(P, T));
    fprintf (fout, "\n");
  }
}



void
cBNModel::Write_nuPT_Table()
{
  FILE *fout = fopen ("nu_PT.dat", "w");
  fprintf (fout, "# Volume      Temperature   Pressure      nu_TO    sigam_nu\n");

  const double h = 1.5198298e-16;
  const double c = 2.9979246e+10;
  const double hc = h*c;

  double T_P0 = 2000.0;
  double T_P900 = 10000.0;
  // double T_P0   = 2000.0;
  // double T_P900 = 2000.0;

  RamanSpectrum spectrum;

  FILE *f0 = fopen ("nu_T0.dat", "w");
  fprintf (f0, "# Volume      Temperature   Pressure      nu_TO    sigam_nu\n");

  for (int iV=0; iV < Raman.Volume.size(); iV++) {
    double V = Raman.Volume[iV];
    vector<double> E(Raman.Energies.size());
    for (int i=0; i<E.size(); i++)  
      E[i] = Raman.Energies[iV](i);
    spectrum.SetEnergies(E);
    spectrum.SetVolume(V);
    double Ps = Static.Pressure(V);
    double x = Ps/900.0;
    double Tmax = T_P0 + x*(T_P900 - T_P0);
    

    spectrum.SetTemp (10.0);
    double P0 = Static.Pressure(V) + Phonon.P(V,10.0);
    double nu0 = spectrum.MeanFrequency();
    fprintf (f0, "%12.5f  %12.5f  %12.5f   %12.5f  %12.5f %12.5f\n",
	     V, 0.0, P0, nu0, Raman.Sigma[iV](0)/hc, Phonon.P(V,10.0));
    for (int iT=0; iT<100; iT++) {
      double x = (double)iT/99.0;
      double T = 10.0*(1.0-x) + x*Tmax;
      spectrum.SetTemp (T);
      double P = Static.Pressure(V) + Phonon.P(V,T);
      double nu = spectrum.MeanFrequency();
      fprintf (fout, "%12.5f  %12.5f  %12.5f   %12.5f  %12.5f %12.5f\n",
    	       V, T, P, nu, Raman.Sigma[iV](0)/hc, Phonon.P(V,T));
    }
    // for (double T=10.0; T <= Tmax; T+=10.0) {
    //   spectrum.SetTemp (T);
    //   double P = Static.Pressure(V) + Phonon.P(V,T);
    //   double nu = spectrum.MeanFrequency();
    //   fprintf (fout, "%12.5f  %12.5f  %12.5f   %12.5f  %12.5f %12.5f\n",
    // 	       V, T, P, nu, Raman.Sigma[iV](0)/hc, Phonon.P(V,T));
    // }
  }
  fclose (fout);
  fclose (f0);
}


void 
cBNModel::SetAll (string staticName, 
		  string phononName, 
		  string ramanName)
{
  const int Nsamples=1000;
  const double h = 1.5198298e-16;
  const double c = 2.9979246e+10;
  const double hc = h*c;

  SetStatic (staticName);
  // Set phonon
  SetPhonon (phononName);
  // Set Raman
  SetRaman (ramanName);

  CommunicatorClass comm;
  RandomClass rand(comm);
  rand.Init();
  TinyVector<double,4> params(80.0, 0.51, 0.0155, 3.7), errors;  
  Static.SetParams (params);
  Array<double,1> V, E, sigma, Epert;
  
  IOSectionClass stin;
  assert(stin.OpenFile (staticName));
  assert(stin.ReadVar("V0",  params[0]));
  assert(stin.ReadVar("E0",  params[1]));
  assert(stin.ReadVar("B0",  params[2]));
  assert(stin.ReadVar("B0p", params[3]));
  params[2] /= 29421.01;
	
  assert(stin.ReadVar("V", V));
  assert(stin.ReadVar("E", E));
  Epert.resize(E.size());
  assert(stin.ReadVar("Sigma", sigma));
  Static.SetParams (params);

  RamanSpectrum spec;
  double T = 300.0;
  spec.SetTemp(T);

  const double Pmax = 1000.0;
  const int numP = 1001;
  Array<double,1> nuSum(numP), nu2Sum(numP), nuError(numP);
  for (int i=0; i<numP; i++)
    nuSum(i) = nu2Sum(i) = nuError(i) = 0.0;

  for (int i=0; i<Nsamples; i++) {
    DatchiRamanModelNoT datchiModel;
    TinyVector<double,3> datchiParams(1062.3, 466.6, 1.289);
    datchiModel.SetParams(datchiParams);
    datchiModel.SetB0p (3.87);
    NonlinearFitClass<3,DatchiRamanModelNoT> datchiFit(datchiModel);
    Array<double,1> Pfit(Raman.Volume.size());
    Array<double,1> nufit(Raman.Volume.size());
    Array<double,1> sigmanu (Raman.Volume.size());


    NonlinearFitClass<4,VinetEOSClass> fitter(Static);
    for (int j=0; j<E.size(); j++) Epert(j) = E(j) + rand.LocalGaussian(sigma(j));
    fitter.Fit (V, E, sigma, params);
    Static.SetParams(params);
    double B0p = params[3];
    datchiModel.SetB0p (B0p);
    // Loop over Raman volumes
    for (int iV=0; iV<Raman.Volume.size(); iV++) {
      double V = Raman.Volume[iV];
      // Compute frequency for this volume.
      vector<double> E(Raman.Energies.size());
      for (int i=0; i<E.size(); i++)  
      E[i] = Raman.Energies[iV](i);
      spec.SetEnergies(E);
      spec.SetVolume(V);
      sigmanu(iV) = Raman.Sigma[iV](0)/hc;
      nufit(iV) = spec.MeanFrequency() + rand.LocalGaussian(sigmanu(iV));

      // Compute pressure for this volume.
      Pfit(iV) = Static.Pressure(V) + Phonon.P(V, 300.0);
    }
    cerr << "Pfit    = " << Pfit << endl;
    cerr << "nufit   = " << nufit << endl;
    cerr << "sigmanu = " << sigmanu << endl;

    // Now fit to the Datchi form
    datchiFit.Fit(Pfit, nufit, sigmanu, datchiParams);
    datchiModel.SetParams(datchiParams);
    cerr << "Params = " << datchiParams << endl;
    for (int iP=0; iP< numP; iP++) {
      double P = (double)iP/(double)(numP-1) * Pmax;
      double nu = datchiModel(P);
      nuSum(iP) += nu;
      nu2Sum(iP) += nu*nu;
    }
  }
  FILE *fout = fopen ("MyDatchiModel.dat", "w");
  fprintf (fout, "# Pressure (GPa)   nu (1/cm)       error\n");
  for (int iP=0; iP< numP; iP++) {
    double P = (double)iP/(double)(numP-1) * Pmax;
    double nuAvg = nuSum(iP)/(double)Nsamples;
    double nuErr = sqrt((nu2Sum(iP)/(double)Nsamples - nuAvg*nuAvg));
    fprintf (fout, "%16.12e %16.12e %16.12e\n", P, nuAvg, nuErr);
  }
  fclose(fout);
}




main(int argc, char **argv)
{
  cBNModel model;
  // if (argc > 3) {
  //   model.SetAll (argv[1], argv[2], argv[3]);
    
  // }
  if (argc > 3) {
    model.SetStatic (argv[1]);
    model.SetPhonon (argv[2]);
    model.SetRaman  (argv[3]);
    FILE *fout = fopen ("300K_Model.dat", "w");
    for (double V=40.0; V<=95.0; V+=0.01)
      fprintf (fout, "%5.2f %10.6f\n", V, model.P(V, 300.0));
    fclose(fout);
  }

  model.Write_nuPT_Table();
  ThermalExpansion (model);
  CalcProperties(model);
}
