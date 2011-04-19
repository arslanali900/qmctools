#include "APWClass.h"

void
APWClass::evaluate (int spin, int atom, int ik, int band,
		    Array<double,1> &r, Array<complex<double>,2> &u_lm_r)
{
  u_lm_r = complex<double>();

  for (int l=0; l<=lMax; l++) {
    int order = MatchCoefs(spin,atom,ik).getOrder();
    int numG   = MatchCoefs(spin,atom,ik).getNumG();
    for (int iord=0; iord<order; iord++) {
      for (int ir=0; ir<r.size(); ir++) {
	double uval = u(atom, l)[iord](r(ir));
	for (int m=-l; m<=l; m++) {
	  int lm = l*(l+1)+m;
	  for (int iG=0; iG<numG; iG++) {
	    complex<double> z = MatchCoefs(spin, atom, ik)(l, m, iord, iG);
	    u_lm_r(lm, ir) += uval* Phi(spin,ik,band)(iG)
	      * MatchCoefs(spin, atom, ik)(l, m, iord, iG);
	  }
	}
      }
    }
  }
}


void
APWClass::evaluate (int spin, int atom, int ik, int band,
		    Array<complex<double>,2> &u_lm_r,
		    Array<complex<double>,1> &du_lm_dr)
{
  int num_r = r_data(atom).size();
  int num_lm = (lMax+1)*(lMax+1);
  u_lm_r.resize  (num_lm, num_r);
  du_lm_dr.resize(num_lm);

  u_lm_r   = complex<double>();
  du_lm_dr = complex<double>();

  for (int l=0; l<=lMax; l++) {
    int order = MatchCoefs(spin,atom,ik).getOrder();
    int numG   = MatchCoefs(spin,atom,ik).getNumG();
    for (int iord=0; iord<order; iord++) {
      for (int m=-l; m<=l; m++) {
	int lm = l*(l+1)+m;
	complex<double> *restrict match = 
	  &(MatchCoefs(spin,atom,ik)(l, m, iord, 0));
	for (int iG=0; iG<numG; iG++) {
	  //complex<double> z = MatchCoefs(spin, atom, ik)(l, m, iord, iG);
	  complex<double> phiMatch = Phi(spin,ik,band)(iG) * match[iG];
	  double *restrict uval = &(u_data(atom, l)(iord,0));
	  for (int ir=0; ir<num_r; ir++) {
	    u_lm_r(lm, ir) += uval[ir]* phiMatch;//Phi(spin,ik,band)(iG) * match[iG];
	    //* MatchCoefs(spin, atom, ik)(l, m, iord, iG);
	  }
	}
      }
    }
  }
  

  for (int l=0; l<=lMax; l++) {
    int order = MatchCoefs(spin,atom,ik).getOrder();
    int numG   = MatchCoefs(spin,atom,ik).getNumG();
    for (int iord=0; iord<order; iord++) {
      double duval = du_dr_data(atom, l)(iord);
      for (int m=-l; m<=l; m++) {
	int lm = l*(l+1)+m;
	complex<double> *match = &(MatchCoefs(spin,atom,ik)(l, m, iord, 0));
	for (int iG=0; iG<numG; iG++) {
	  complex<double> z = MatchCoefs(spin, atom, ik)(l, m, iord, iG);
	  du_lm_dr(lm) += duval * Phi(spin,ik,band)(iG) * match[iG];
	  //* MatchCoefs(spin, atom, ik)(l, m, iord, iG);
	}
      }
    }
  }
}



// Adds on the local orbital contribution
void
LocalOrbitalClass::add (int spin, int atom, int ik, int band,
			Array<double,1> &r, Array<complex<double>,2> &u_lm_r)
{
  int numLocal = AtomMap.size();
  for (int iloc=0; iloc<numLocal; iloc++) {
    if (AtomMap(iloc) == atom) {
      int l = lMap(iloc);
      int m = mMap(iloc);
      int ilo = iloMap(iloc);
      int lm = l*(l+1)+m;
      for (int ir=0; ir<r.size(); ir++ ) 
	u_lm_r(lm,ir) += Phi(spin,ik,band,iloc) * v(atom)[ilo](r(ir));
    }
  }
}

void
LocalOrbitalClass::add (int spin, int atom, int ik, int band,
			Array<complex<double>,2> &u_lm_r)
{
  int numLocal = AtomMap.size();
  int num_r = u_lm_r.extent(1);
  for (int iloc=0; iloc<numLocal; iloc++) {
    if (AtomMap(iloc) == atom) {
      int l = lMap(iloc);
      int m = mMap(iloc);
      int ilo = iloMap(iloc);
      int lm = l*(l+1)+m;
      for (int ir=0; ir<num_r; ir++) 
	u_lm_r(lm,ir) += Phi(spin,ik,band,iloc) * v_data(atom)(ilo,ir);
    }
  }
}


void
CoreStateClass::SetupSplines()
{
  int N = r.size();
  vector<double> rvec(N), fvec(N), gvec(N);
  for (int i=0; i<N; i++) {
    rvec[i] = r(i);
    fvec[i] = f0(i);
    gvec[i] = g0(i);
  }
  SimpleGrid rgrid;
  rgrid.Init(rvec);
  f0Spline.Init (rgrid, fvec);
  g0Spline.Init (rgrid, gvec);
  

}
