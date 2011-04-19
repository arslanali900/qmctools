#ifndef APW_CLASS_H
#define APW_CLASS_H

#include <complex>
#include <Common/Blitz.h>
#include "CubicSpline.h"

// One instance of this per (spin, atom, k-point)
class APWMatch
{
private:
  // Indices are (lm, order, G-vector index)
  Array<complex<double>,3> MatchingCoefs;
public:
  inline int getOrder() { return MatchingCoefs.extent(1); }
  inline int getNumG()  { return MatchingCoefs.extent(2); }

  inline void resize(int lmax, int order, int numG) 
  { MatchingCoefs.resize((lmax+1)*(lmax+1), order,numG); }

  inline complex<double> operator()(int l, int m, int ord, int iG) const
  { return MatchingCoefs(l*(l+1)+m, ord, iG); }

  inline complex<double>& operator()(int l, int m, int ord, int iG)
  { return MatchingCoefs(l*(l+1)+m, ord, iG); }
};



class APWClass 
{
public:
  int lMax, NumAtoms;

 // APW coefficients.  Indices are (spin,ik,iband)(ig).
  Array<Array<complex<double>,1>,3> Phi;
  // Indices are (atom, l)[order]
  Array<vector<CubSpline>,2> u;
  // Indices are (atom,l)(order,ir);
  Array<Array<double,2>,2> u_data;
  // This stores the derivative at the muffin-tin boundary
  // Indices are (atom,l)(order)
  Array<Array<double,1>,2> du_dr_data;
  // Indices are (atom)(ir)
  Array<Array<double,1>,1> r_data;
  
 
  // Indices are (spin, atom, k-point)
  Array<APWMatch,3> MatchCoefs;

  // For a given (spin, atom, ik, band), do the sum over 
  // plane-wave coefficients, giving u_lm(r).  The indices of the 
  // output, u_lm_r, are (l*(l+1)+m, ir)
  void evaluate (int spin, int atom, int ik, int band,
		 Array<double,1> &r, Array<complex<double>,2> &u_lm_r);

  // Do the same as above, but evaluate on the grid in r_data;
  void evaluate (int spin, int atom, int ik, int band,
		 Array<complex<double>,2> &u_lm_r,
		 Array<complex<double>,1> &du_lm_dr);

  APWClass() : lMax(0)
  {

  }
};



class LocalOrbitalClass
{
public:
  // Local orbital coefficients
  // Indices are (spin, ik, band, iloc)
  Array<complex<double>,4> Phi;

  // Index is iloc
  Array<int,1> AtomMap, lMap, mMap, iloMap;

  // Indices are (atom, ilo)
  Array<vector<CubSpline>,1> v;
  // Indices are (atom)(ilo,ir)
  Array<Array<double,2>,1> v_data;
  
  void add (int spin, int atom, int ik, int band,
	    Array<double,1> &r, Array<complex<double>,2> &u_lm_r);
  void add (int spin, int atom, int ik, int band,
	    Array<complex<double>,2> &u_lm_r);
};


class CoreStateClass
{
private:
  CubSpline f0Spline, g0Spline;

public:
  void SetupSplines();

  inline double f(double r) { return f0Spline(r); }
  inline double g(double r) { return g0Spline(r); }
  inline double dfdr(double r) { return f0Spline.Deriv(r); }
  inline double dgdr(double r) { return g0Spline.Deriv(r); }

  Array<double,1> r, g0, f0;
  int atom, n, l, k;
  double eigenvalue;
  
  CoreStateClass() { }

  CoreStateClass(const CoreStateClass &state)
  {
    r.resize(state.r.size()); r = state.r;
    g0.resize(state.g0.size()); g0 = state.g0;
    f0.resize(state.f0.size()); f0 = state.f0;
    atom = state.atom;
    n    = state.n;
    l    = state.l;
    k    = state.k;
    eigenvalue = state.eigenvalue;
  }
};


#endif
