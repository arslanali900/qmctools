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

#ifndef OPTIMIZED_BREAKUP_H
#define OPTIMIZED_BREAKUP_H

#include <blitz/array.h>
#include "LatticeClass.h"

using namespace blitz;

class BasisClass
{
public:
  //protected:
  double r_c;
  LatticeClass Lattice;
  // Unit cell volume
  double Omega;
public:
  /// Set the cutoff radius
  virtual void Set_rc(double rc) = 0;
  inline double Get_rc() { return r_c; }
  inline void SetBox (TinyVector<double,3> box);
  inline void SetLattice (TinyMatrix<double,3,3> lattice);
  inline LatticeClass& GetLattice();
  /// Returns the number of basis elements
  virtual int NumElements() = 0;
  /// Returns the basis element n evaluated in real space at r
  virtual double h(int n, double r) = 0;
  /// Returns the basis element n evaluated in k space at k
  virtual double c(int n, double k) = 0;
  virtual double c0(int n);
  virtual double dc_dk(int n, double k) = 0;
  double c_numerical (int n, double k);
  /// This returns the coefficent of the nth basis function
  //virtual double  Get_t(int n) const     = 0;
  /// This sets the coefficent of the nth basis function
  //virtual void    Set_t(int n, double t) = 0;
  /// This returns the linear combination of the basis functions with
  /// coefficients t_n
  //virtual double f (double r) = 0;
  BasisClass() : r_c (0.0)
  { /* do nothing */ }
};


class OptimizedBreakupClass
{
private:
  BasisClass &Basis;
  void Addk(double k, double degeneracy=1.0);
public:
  // First element is |k|, second is degeneracy of the point.
  Array<TinyVector<double,2>,1> kpoints;
  void SetkVecs(double kc, double kcont, double kMax);
  /// kc is the k-space cutoff for the Ewald sum.  
  /// kMax is largest k we use in determining the error in the breakup.  
  /// t is the set of coefficients of the breakup.
  /// inFit is a boolean array telling whether t_n should be optimized
  /// or left at its initial value.  Returns chi-squared for the breakup.
  double DoBreakup (const Array<double,1> &Vk, Array<double,1> &t, 
		    const Array<bool,1> &adjust);
  /// Same as above, but we assume that all t's are adjusted.
  double DoBreakup (const Array<double,1> &Vk, Array<double,1> &t);
  inline double Vlong_r (Array<double,1> &t, double r);
  inline double Vlong_k (Array<double,1> &t, double k);
  OptimizedBreakupClass (BasisClass &basis) : Basis(basis)
  { /* Do nothing */ }
};


inline double
OptimizedBreakupClass::Vlong_r (Array<double,1> &t, double r)
{
  int N = Basis.NumElements();
  double V = 0.0;
  for (int i=0; i<N; i++)
    V += t(i) * Basis.h(i,r);
  return V;
}

inline double
OptimizedBreakupClass::Vlong_k (Array<double,1> &t, double k)
{
  int N = Basis.NumElements();
  double V = 0.0;
  for (int i=0; i<N; i++)
    V += t(i) * Basis.c(i,k);
  return V;
}


inline void BasisClass::SetBox (TinyVector<double,3> box)
{
  TinyMatrix<double,3,3> a;
  a(0,0)=box[0]; a(0,1)=0.0;    a(0,2)=0.0;
  a(1,0)=0.0;    a(1,1)=box[1]; a(1,2)=0.0;
  a(2,0)=0.0;    a(2,1)=0.0;    a(2,2)=box[2];
  
  SetLattice (a);
}

inline void
BasisClass::SetLattice(TinyMatrix<double,3,3> lattice)
{
  Lattice.Set(lattice);
  Omega = Lattice.Volume();
}

inline LatticeClass&
BasisClass::GetLattice ()
{
  return Lattice;
}



/// Locally Piecewise Quintic Hermite Interpolant
class LPQHI_BasisClass : public BasisClass
{
  /// public is HACK
  //private:
public:
  int NumKnots;
  double delta, deltaInv;
  TinyMatrix<double,3,6> S;
  /// The following are helpers to calculate the Fourier tranform of
  /// the basis functions
  inline complex<double> Eplus(int i, double k, int n);
  inline complex<double> Eminus(int i, double k, int n);
  inline double Dplus(int i, double k, int n);
  inline double Dminus(int i, double k, int n);
  inline complex<double> dEplus_dk(int i, double k, int n);
  inline complex<double> dEminus_dk(int i, double k, int n);
  inline double dDplus_dk(int i, double k, int n);
  inline double dDminus_dk(int i, double k, int n);
  Array<double,1> tvec;
public:
  inline double GetDelta() { return delta; }
  // n must be at least 2;
  void SetNumKnots(int n);
  void Set_rc(double rc);
  int NumElements();
  double h(int n, double r);
  double c(int n, double k);
  double dc_dk(int n, double k);
  LPQHI_BasisClass() : NumKnots(0), delta(0.0) 
  { 
    S = 
      1.0,   0.0,   0.0, -10.0,  15.0,  -6.0,
      0.0,   1.0,   0.0,  -6.0,   8.0,  -3.0,
      0.0,   0.0,   0.5,  -1.5,   1.5,  -0.5;
  }
};

inline complex<double> LPQHI_BasisClass::Eplus(int i, double k, int n)
{
  complex<double> eye(0.0, 1.0);

  if (n == 0) {
    complex<double> e1(cos(k*delta)-1.0, sin(k*delta));
    complex<double> e2(cos(k*delta*i),   sin(k*delta*i));
    return (-(eye/k)*e1*e2);
  }
  else {
    complex<double> t1, t2;
    double sign = 1.0;
    t1 = complex<double>(cos(k*delta*(i+1)),sin(k*delta*(i+1)));
    t2=-(double)n/delta*Eplus(i,k,n-1);;
    return (-(eye/k)*(t1+t2));
  }
}

inline complex<double> LPQHI_BasisClass::dEplus_dk(int i, double k, int n)
{
  complex<double> eye(0.0, 1.0);

  if (n == 0) {
    complex<double>  e1(cos(k*delta)-1.0,        sin(k*delta));
    complex<double> de1(-delta*sin(k*delta),     delta*cos(k*delta));
    complex<double>  e2(cos(k*delta*i),          sin(k*delta*i));
    complex<double> de2(-delta*i*sin(k*delta*i), delta*i*cos(k*delta*i));
    return ((eye/(k*k))*e1*e2-(eye/k)*(de1*e2+e1*de2));
  }
  else {
    complex<double> t1, t2, dt1, dt2;
    double sign = 1.0;
    t1  = complex<double>(cos(k*delta*(i+1)),sin(k*delta*(i+1)));
    dt1 = complex<double>(-delta*(i+1)*sin(k*delta*(i+1)),
			  delta*(i+1)*cos(k*delta*(i+1)));
    t2  = -(double)n/delta*Eplus(i,k,n-1);
    dt2 = -(double)n/delta*dEplus_dk(i,k,n-1);
    return ((eye/(k*k))*(t1+t2) -(eye/k)*(dt1+dt2));
  }
}

inline complex<double> LPQHI_BasisClass::Eminus(int i, double k, int n)
{
  complex<double> eye(0.0, 1.0);

  if (n == 0) {
    complex<double> e1(cos(k*delta)-1.0, -sin(k*delta));
    complex<double> e2(cos(k*delta*i),    sin(k*delta*i));
    return (-(eye/k)*e1*e2);
  }
  else {
    complex<double> t1, t2;
    double sign = (n & 1) ? -1.0 : 1.0;
    t1 = sign*
      complex<double> (cos(k*delta*(i-1)),sin(k*delta*(i-1)));
    t2=-(double)n/delta*Eminus(i,k,n-1);
    return (-(eye/k)*(t1+t2));
  }
}

inline complex<double> LPQHI_BasisClass::dEminus_dk(int i, double k, int n)
{
  complex<double> eye(0.0, 1.0);

  if (n == 0) {
    complex<double>  e1(cos(k*delta)-1.0,        -sin(k*delta));
    complex<double> de1(-delta*sin(k*delta),     -delta*cos(k*delta));
    complex<double>  e2(cos(k*delta*i),          sin(k*delta*i));
    complex<double> de2(-delta*i*sin(k*delta*i), delta*i*cos(k*delta*i));
    return ((eye/(k*k))*e1*e2-(eye/k)*(de1*e2+e1*de2));
  }
  else {
    complex<double> t1, t2, dt1, dt2;
    double sign = (n & 1) ? -1.0 : 1.0;
    t1  = sign * complex<double>(cos(k*delta*(i-1)),sin(k*delta*(i-1)));
    dt1 = sign * complex<double>(-delta*(i-1)*sin(k*delta*(i-1)),
			  delta*(i-1)*cos(k*delta*(i-1)));
    t2  = -(double)n/delta*Eminus(i,k,n-1);
    dt2 = -(double)n/delta*dEminus_dk(i,k,n-1);
    return ((eye/(k*k))*(t1+t2) -(eye/k)*(dt1+dt2));
  }
}



inline double LPQHI_BasisClass::Dplus(int i, double k, int n)
{
  complex<double> eye(0.0, 1.0); 
  complex<double> Z1 = Eplus(i,k,n+1);
  complex<double> Z2 = Eplus(i,k,n);
  return 4.0*M_PI/(k*Omega)*(delta* Z1.imag() + i*delta*Z2.imag());
}

inline double LPQHI_BasisClass::dDplus_dk(int i, double k, int n)
{
  complex<double> eye(0.0, 1.0); 
  complex<double> Z1 = Eplus(i,k,n+1);
  complex<double> Z2 = Eplus(i,k,n);
  complex<double> dZ1 = dEplus_dk(i,k,n+1);
  complex<double> dZ2 = dEplus_dk(i,k,n);
  return -4.0*M_PI/(k*k*Omega)*(delta* Z1.imag() + i*delta*Z2.imag()) +
    4.0*M_PI/(k*Omega)*(delta*dZ1.imag() + i*delta*dZ2.imag());
}

inline double LPQHI_BasisClass::Dminus(int i, double k, int n)
{
  complex<double> eye(0.0, 1.0); 
  complex<double> Z1 = Eminus(i,k,n+1);
  complex<double> Z2 = Eminus(i,k,n);
  return -4.0*M_PI/(k*Omega)*(delta* Z1.imag() + i*delta*Z2.imag());
}

inline double LPQHI_BasisClass::dDminus_dk(int i, double k, int n)
{
  complex<double> eye(0.0, 1.0); 
  complex<double> Z1  = Eminus(i,k,n+1);
  complex<double> dZ1 = dEminus_dk(i,k,n+1);
  complex<double> Z2  = Eminus(i,k,n);
  complex<double> dZ2 = dEminus_dk(i,k,n);
  return 4.0*M_PI/(k*k*Omega)*(delta* Z1.imag() + i*delta*Z2.imag())
    -4.0*M_PI/(k*Omega)*(delta* dZ1.imag() + i*delta*dZ2.imag());
}


#endif
