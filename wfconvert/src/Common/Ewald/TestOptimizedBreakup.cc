/////////////////////////////////////////////////////////////
// Copyright (C) 2003-2006 Bryan Clark and Kenneth Esler   //
//                                                         //
// This program is free software; you can redistribute it  //
// and/or modify it under the terms of the GNU General     //
// Public License as published by the Free Software        //
// Foundation; either version 2 of the License, or         //
// (at your option) any later version.  This program is    //
// distributed in the hope that it will be useful, but     //
// WITHOUT ANY WARRANTY; without even the implied warranty //
// of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. //  
// See the GNU General Public License for more details.    //
// For more information, please see the PIMC++ Home Page:  //
//           http://pathintegrals.info                     //
/////////////////////////////////////////////////////////////

#include "OptimizedBreakup2.h"

#include <cstdio>
#include <complex>

void WriteBasis(BasisClass &basis)
{
  FILE *fout = fopen ("LPQHI.dat", "w");
  for (double r=0.0; r<=basis.Get_rc(); r+=0.001) {
    fprintf (fout, "%1.12e ", r);
    for (int n=0; n<basis.NumElements(); n++)
      fprintf (fout, "%1.12e ", basis.h(n, r));
    fprintf (fout, "\n");
  }
  fclose (fout);
}

void  TestLPQHI()
{
  LPQHI_BasisClass basis;

  basis.SetNumKnots(6);
  basis.Set_rc(5.0);

  WriteBasis (basis);
}

void Test_kDeriv()
{
  LPQHI_BasisClass basis;
  basis.SetNumKnots(6);
  basis.Set_rc(5.0);
  double epsilon=1.0e-6;
  double k = 0.347;

  for (int i=0; i<5; i++)
    for (int n=0; n<5; n++) {
      complex<double> FD = 
	(basis.Dminus(i,k+epsilon,n)-basis.Dminus(i,k-epsilon,n))/(2.0*epsilon);
      complex<double> an = basis.dDminus_dk(i,k,n);
      fprintf (stderr, "%1.12e + %1.12ei      %1.12e + %1.12ei\n",
	       real(FD), imag(FD), real(an), imag(an));
    }
}

void TestBasisDeriv()
{
  LPQHI_BasisClass basis;
  basis.SetNumKnots(3);
  basis.Set_rc(2.0);
  basis.SetBox(TinyVector<double,3>(2.0, 2.0, 2.0));
  double epsilon=1.0e-5;
  double k = 1.0;
  
  
  for (int i=0; i<basis.NumElements(); i++) {
    double FD = (basis.c(i,k+epsilon)-basis.c(i,k-epsilon))/(2.0*epsilon);
    double an = basis.dc_dk(i,k);
    fprintf (stderr, "%1.12e  %1.12e\n", FD, an);
  }
}







void TestCoulomb()
{
  LPQHI_BasisClass basis;
  TinyVector<double,3> box;
  box = 10.0, 10.0, 10.0;
  basis.SetBox (box);
  double Omega = box[0]*box[1]*box[2];
  basis.SetNumKnots(20);
  basis.Set_rc(5.0);

//   double k = 1;
//   for (int n=0; n<basis.NumElements(); n++) {
//     double c = basis.c(n,k);
//     double c_numerical = basis.c_numerical(n,k);
//     fprintf (stderr, "n = %d\n", n);
//     fprintf (stderr, "c           = %1.12e\n", c);
//     fprintf (stderr, "c_numerical = %1.12e\n", c_numerical);
//   }

  OptimizedBreakupClass breakup(basis);
  breakup.SetkVecs (2.0, 25.0, 1000.0);
  double delta = basis.GetDelta();

  Array<double,1> Vk(breakup.kpoints.size());
  for (int i=0; i<breakup.kpoints.size(); i++) {
    double k = breakup.kpoints(i)[0];
    double k0 = 2.0*M_PI;
    Vk(i) = 4.0*M_PI/(Omega*(k*k)) * cos(k*basis.Get_rc());
  }
  int N = basis.NumElements();
  Array<double,1> t(N);
  Array<bool,1> adjust(N);
  adjust = true;
  t = 0.0;
  adjust(N-3) = false; t(N-3) = -0.2;
  adjust(N-2) = false; t(N-2) = 1.0/25.0*delta;
  adjust(N-1) = false; t(N-1) = -1/125.0*delta*delta;
  adjust(1)   = false; t(1) = 0.0;

  double chi2 = breakup.DoBreakup(Vk, t, adjust);
  cerr << "chi-squared = " << chi2 << endl;
  cerr << "t = " << t << endl;

  FILE *fout = fopen ("Vlong.dat", "w");
  
  for (double r=0.0; r<10.0; r+=0.001) {
    double v = 0.0;
    for (int n=0; n<t.rows(); n++)
      v += t(n)*basis.h(n,r);
    fprintf (fout, "%1.12e %1.12e\n", r, v);

  }
  fclose(fout);

}

void TestFeO()
{
  Mat3 a;
  a(0,0)=0.5; a(0,1)=0.5; a(0,2)=1.0;  
  a(1,0)=0.5; a(1,1)=1.0; a(1,2)=0.5;
  a(2,0)=1.0; a(2,1)=0.5; a(2,2)=0.5;
  a = 2.0*8.703 * a;

  LatticeClass lattice(a);

  vector<Vec3> Olist, Felist;
  Olist.push_back (Vec3(0.125, 0.125, 0.125));
  Olist.push_back (Vec3(0.375, 0.375, 0.375));
  Felist.push_back(Vec3(0.000, 0.000, 0.000));
  Felist.push_back(Vec3(0.250, 0.250, 0.250));
  Olist.push_back (Vec3(0.125, 0.125, 0.625));
  Olist.push_back (Vec3(0.375, 0.375, 0.875));
  Felist.push_back(Vec3(0.000, 0.000, 0.500));
  Felist.push_back(Vec3(0.250, 0.250, 0.750));
  Olist.push_back (Vec3(0.125, 0.625, 0.125));
  Olist.push_back (Vec3(0.375, 0.875, 0.375));
  Felist.push_back(Vec3(0.000, 0.500, 0.000));
  Felist.push_back(Vec3(0.250, 0.750, 0.250));
  Olist.push_back (Vec3(0.125, 0.625, 0.625));
  Olist.push_back (Vec3(0.375, 0.875, 0.875));
  Felist.push_back(Vec3(0.000, 0.500, 0.500));
  Felist.push_back(Vec3(0.250, 0.750, 0.750));
  Olist.push_back (Vec3(0.625, 0.125, 0.125));
  Olist.push_back (Vec3(0.875, 0.375, 0.375));
  Felist.push_back(Vec3(0.500, 0.000, 0.000));
  Felist.push_back(Vec3(0.750, 0.250, 0.250));
  Olist.push_back (Vec3(0.625, 0.125, 0.625));
  Olist.push_back (Vec3(0.875, 0.375, 0.875));
  Felist.push_back(Vec3(0.500, 0.000, 0.500));
  Felist.push_back(Vec3(0.750, 0.250, 0.750));
  Olist.push_back (Vec3(0.625, 0.625, 0.125));
  Olist.push_back (Vec3(0.875, 0.875, 0.375));
  Felist.push_back(Vec3(0.500, 0.500, 0.000));
  Felist.push_back(Vec3(0.750, 0.750, 0.250));
  Olist.push_back (Vec3(0.625, 0.625, 0.625));
  Olist.push_back (Vec3(0.875, 0.875, 0.875));
  Felist.push_back(Vec3(0.500, 0.500, 0.500));
  Felist.push_back(Vec3(0.750, 0.750, 0.750));

  vector<Vec3> OCart, FeCart;
  for (int i=0; i<Olist.size(); i++) 
    OCart.push_back(lattice.Reduced2Cart(Olist[i]));
  for (int i=0; i<Felist.size(); i++) 
    FeCart.push_back(lattice.Reduced2Cart(Felist[i]));

  double rcut = lattice.rMax();
  cerr << "rcut = " << rcut << endl;
  double Omega = lattice.Volume();
  LPQHI_BasisClass basis;
  basis.SetLattice(a);
  basis.Set_rc(rcut);
  basis.SetNumKnots(40);
  OptimizedBreakupClass breakup(basis);
  double nonDimCut = 15.0;
  double kcut = nonDimCut / rcut;
  breakup.SetkVecs (kcut, 30.0, 1000.0);
  Array<double,1> Xk(breakup.kpoints.size());
  for (int i=0; i<breakup.kpoints.size(); i++) {
    double k = breakup.kpoints(i)[0];
    Xk(i) = -4.0*M_PI/(Omega*(k*k)) * cos(k*basis.Get_rc());
  }
  int N = basis.NumElements();
  Array<double,1> t(N);
  Array<bool,1> adjust(N);
  adjust = true;
  t = 0.0;
  double delta = basis.GetDelta();
  t(N-3) = 1.0/rcut;                           adjust(N-3) = false;
  t(N-2) = -1.0/(rcut*rcut)*delta;             adjust(N-2) = false;
  t(N-1) = 2.0/(rcut*rcut*rcut)*delta*delta;   adjust(N-1) = false;
  
  double chi2 = breakup.DoBreakup(Xk, t, adjust);
  fprintf (stderr, "%1.16e %1.16e\n", rcut*kcut, chi2);

  double V = 0.0;
  FILE *fout = fopen ("Vlong.dat", "w");
  for (double r=0.001; r<rcut; r+=0.001) 
    fprintf (fout, "%5.3f %16.12e\n", r, breakup.Vlong_r(t, r));
  fclose(fout);


  //////////////////////
  // Short-range part //
  //////////////////////
  double Vshort = 0.0;
  // O-O contribution
  for (int i=0; i<OCart.size(); i++) {
    for (int j=i+1; j<OCart.size(); j++) {
      Vec3 disp = lattice.MinImage(OCart[i]-OCart[j]);
      double dist = sqrt(dot(disp,disp));
      if (dist < rcut)
	Vshort += (6.0*6.0)*(1.0/dist - breakup.Vlong_r(t, dist));
    }
  }
  // O-O contribution
  for (int i=0; i<FeCart.size(); i++) {
    for (int j=i+1; j<FeCart.size(); j++) {
      Vec3 disp = lattice.MinImage(FeCart[i]-FeCart[j]);
      double dist = sqrt(dot(disp,disp));
      if (dist < rcut)
	Vshort += (16.0*16.0)*(1.0/dist - breakup.Vlong_r(t, dist));
    }
  }
  // Fe-O contribution
  for (int i=0; i<OCart.size(); i++) {
    for (int j=0; j<FeCart.size(); j++) {
      Vec3 disp = lattice.MinImage(OCart[i]-FeCart[j]);
      double dist = sqrt(dot(disp,disp));
      if (dist < rcut)
	Vshort += (6.0*16.0)*(1.0/dist - breakup.Vlong_r(t, dist));
    }
  }
  cerr << "Short-range contribution = " << Vshort << endl;
  
  // Long-range contribution
  Vec3 kvec (0.0, 0.0, 0.0);
  vector<Vec3> Gvecs = lattice.GenGvecs(kvec, 0.5*kcut*kcut);
  cerr << "Using " << Gvecs.size() << " G-vectors.\n";
  vector<double> Vlong_k(Gvecs.size());
  for (int i=0; i<Gvecs.size(); i++) {
    double k = sqrt(dot(Gvecs[i], Gvecs[i]));
    if (k > 1.0e-12) {
      double xk = -4.0*M_PI/(Omega*(k*k)) * cos(k*basis.Get_rc());
      Vlong_k[i] = breakup.Vlong_k(t, k) - xk;
    }
    else 
      Vlong_k[i] = breakup.Vlong_k(t, 0.0);
  }
		     
  vector<complex<double> > rhok_O, rhok_Fe;;
  // Compute rho_k
  for (int gi=0; gi<Gvecs.size(); gi++) {
    complex<double> rhoO(0.0, 0.0);
    complex<double> rhoFe(0.0, 0.0);
    for (int i=0; i<OCart.size(); i++) {
      double phase = dot(OCart[i], Gvecs[gi]);
      rhoO += complex<double>(cos(phase), sin(phase));
    }
    rhok_O.push_back(rhoO);
    for (int i=0; i<FeCart.size(); i++) {
      double phase = dot(FeCart[i], Gvecs[gi]);
      rhoFe+= complex<double>(cos(phase), sin(phase));
    }
    rhok_Fe.push_back(rhoFe);
  }
    
  // Now add long-range contribution
  // O-O
  double Vlong = 0.0;
  for (int ki=0; ki<rhok_O.size(); ki++) {
    double k = sqrt(dot(Gvecs[ki], Gvecs[ki]));
    if (fabs(k) > 1.0e-10) {
      Vlong += 0.5* 6.0* 6.0*norm(rhok_O[ki])*Vlong_k[ki];
      Vlong += 0.5*16.0*16.0*norm(rhok_Fe[ki])*Vlong_k[ki];
      Vlong += 6.0*16.0*real(rhok_O[ki]*conj(rhok_Fe[ki]))*Vlong_k[ki];
    }
  }
  double ZO = 6.0;
  double ZFe = 16.0;
  //Vlong += 0.5*16.0*16.0*(ZO*ZO+ZFe*ZFe)*breakup.Vlong_k(t,0.0);
  //Vlong += 16.0*16.0*(ZO*ZFe)*breakup.Vlong_k(t,0.0);
  Vlong -= 0.5 * (ZO*ZO + ZFe*ZFe)*16.0*breakup.Vlong_r(t,0.0);
  cerr << "Vlong_r0 = " << breakup.Vlong_r(t, 0.0) << endl;

  double Vshort0 = 2.0*M_PI/Omega * rcut*rcut -
    breakup.Vlong_k(t, 0.0);
  cerr << "Vshort0 = " << Vshort0 << endl;
  cerr << "Vlong0  = " << breakup.Vlong_k(t, 0.0) << endl;
  Vlong -= 0.5*16.0*16.0*(ZFe*ZFe + ZO*ZO)*Vshort0;
  Vlong -= 1.0*16.0*16.0*(ZFe*ZO)         *Vshort0;


  cerr << "Long-range contribution = " << Vlong << endl;
  fprintf (stderr, "Total = %1.12f\n", (Vlong+Vshort)/8.0);
}



void CoulombError()
{
  LPQHI_BasisClass basis;
  TinyVector<double,3> box;
  //  box = 10.0, 10.0, 10.0;
  box = 5.0, 5.0, 5.0;
  basis.SetBox (box);
  double Omega = box[0]*box[1]*box[2];
  basis.SetNumKnots(20);
  double rc = 0.5*min(box[0],min(box[1],box[2]));
  basis.Set_rc(rc);

  FILE *fout = fopen ("chi.dat", "w");
  assert (fout != NULL);
  for (double kcut=0.6; kcut <7; kcut += 0.2) {
    OptimizedBreakupClass breakup(basis);
    breakup.SetkVecs (kcut, 25.0, 1000.0);
    double delta = basis.GetDelta();

    Array<double,1> Xk(breakup.kpoints.size());
    for (int i=0; i<breakup.kpoints.size(); i++) {
      double k = breakup.kpoints(i)[0];
      double k0 = 2.0*M_PI;
      Xk(i) = 4.0*M_PI/(Omega*(k*k)) * cos(k*basis.Get_rc());
    }
    int N = basis.NumElements();
    Array<double,1> t(N);
    Array<bool,1> adjust(N);
    adjust = true;
    t = 0.0;

    double chi2 = breakup.DoBreakup(Xk, t, adjust);
    fprintf (fout, "%1.16e %1.16e\n", rc*kcut, chi2);
    fflush (fout);
//     cerr << "chi-squared = " << chi2 << endl;
//     cerr << "t = " << t << endl;
    
//     FILE *fout = fopen ("Vlong.dat", "w");
//     for (double r=0.0; r<10.0; r+=0.001) {
//       double v = 0.0;
//       for (int n=0; n<t.rows(); n++)
// 	v += t(n)*basis.h(n,r);
//       fprintf (fout, "%1.12e %1.12e\n", r, v);

//     }
//     fclose(fout);
  }
  fclose (fout);
}

double 
GetChi2 (double kc, double L)
{
  double rc = 0.5 * L;
  TinyVector<double,3> box (L,L,L);
  LPQHI_BasisClass basis;
  double Omega = box[0]*box[1]*box[2];
  basis.SetNumKnots(20);
  basis.Set_rc(rc);
  
  OptimizedBreakupClass breakup(basis);
  breakup.SetkVecs (kc, 25.0, 1000.0);
  double delta = basis.GetDelta();
  
  Array<double,1> Xk(breakup.kpoints.size());
  for (int i=0; i<breakup.kpoints.size(); i++) {
    double k = breakup.kpoints(i)[0];
    double k0 = 2.0*M_PI;
    Xk(i) = 4.0*M_PI/(Omega*(k*k)) * cos(k*basis.Get_rc());
  }
  int N = basis.NumElements();
  Array<double,1> t(N);
  Array<bool,1> adjust(N);
  adjust = true;
  t = 0.0;
  
  double chi2 = breakup.DoBreakup(Xk, t, adjust);
  return chi2;
}



void CoulombError2()
{

  FILE *fout = fopen ("/home/esler/Thesis/Ewald/chi2.dat", "w");
  double kc = 3.48;
  assert (fout != NULL);
  for (double L=0.2; L <10.0; L += 0.5) {
    LPQHI_BasisClass basis;
    TinyVector<double,3> box;
    //  box = 10.0, 10.0, 10.0;
    box = L, L, L;
    basis.SetBox (box);
    double Omega = box[0]*box[1]*box[2];
    basis.SetNumKnots(20);
    double rc = 0.5*min(box[0],min(box[1],box[2]));
    basis.Set_rc(rc);

    OptimizedBreakupClass breakup(basis);
    breakup.SetkVecs (kc, 25.0, 1000.0);
    double delta = basis.GetDelta();

    Array<double,1> Xk(breakup.kpoints.size());
    for (int i=0; i<breakup.kpoints.size(); i++) {
      double k = breakup.kpoints(i)[0];
      double k0 = 2.0*M_PI;
      Xk(i) = 4.0*M_PI/(Omega*(k*k)) * cos(k*basis.Get_rc());
    }
    int N = basis.NumElements();
    Array<double,1> t(N);
    Array<bool,1> adjust(N);
    adjust = true;
    t = 0.0;

    double chi2 = breakup.DoBreakup(Xk, t, adjust);
    double chi = sqrt(chi2);
    fprintf (fout, "%1.16e %1.16e\n", rc*kc, log(L*chi));
    fflush (fout);
//     cerr << "chi-squared = " << chi2 << endl;
//     cerr << "t = " << t << endl;
    
//     FILE *fout = fopen ("Vlong.dat", "w");
//     for (double r=0.0; r<10.0; r+=0.001) {
//       double v = 0.0;
//       for (int n=0; n<t.rows(); n++)
// 	v += t(n)*basis.h(n,r);
//       fprintf (fout, "%1.12e %1.12e\n", r, v);

//     }
//     fclose(fout);
  }
  fclose (fout);
}




main()
{
  //  TestBasisDeriv();
//   TestLPQHI();
//   TestCoulomb();
  // CoulombError();
  // TestCoulomb();
  TestFeO();
}
