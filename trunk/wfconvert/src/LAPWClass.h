#ifndef LAPW_CLASS_H
#define LAPW_CLASS_H

#include "LatticeClass.h"
#include "CubicSpline.h"
#include "Ylm.h"
#include <complex>
#include "MuffinTin.h"

class LAPWClass
{
private:
  double Volume, PWNorm;
  // Index is atom number
  Array<int,1> AtomTypes, AtomSpecies, APWorder;
  Array<Vec3,1> AtomPos;
  Array<double,1> AtomRadii;
  int lMax;
  // Indices are atom, lm, order, ik, ig
  Array<complex<double>,5> MatchCoefs;
  // Indices are (ik, iband, iG)
  Array<complex<double>,3> APWcoefs;
  // Indices are (ik, iband, ilocal)
  Array<complex<double>,3> LocalCoefs;
  // Index is lm_index(l,m)
  // Indices are (ik, iG)
  Array<Vec3,2> GVecs;
  Array<Vec3,1> kVecs;
  // Index is ik
  Array<int,1> NumG, NumBands;
  // APW radial functions.  Indices are atom, l, apworder
  Array<CubSpline,3> u;
  ////////////////////////
  // Local orbital info //
  ////////////////////////
  // Local radial functions.  Indices are atom, ilo
  Array<CubSpline,2> v;
  // Index is local coefficient number
  Array<int,1> LocalAtom, LocalSpecies, Local_l, Local_m, Local_ilo;

  ///////////////////////
  // Core orbital info //
  ///////////////////////
  // First index is atom, second is state number
  vector<vector<CubSpline> > CoreSplines;
  Array<int,1> CoreAtom, Core_l;
  
  complex<double> evalCore (int atom, int orb);
  complex<double> evalMuffinTin (int ik, int iband, int iatom, Vec3 r,
				 Vec3 L);
  complex<double> evalInterstitial (int ik, int iband, Vec3 r);

  inline int lm_index(int l, int m);
public:
  LatticeClass Lattice;

  void Read (string fname);
  void SetupMuffinTins(vector<MuffinTinClass> &tins);
  inline complex<double> operator()(int ik, int iband, Vec3 r,
				    bool &inTin);
};


inline int
LAPWClass::lm_index(int l, int m)
{
  return l*(l+1) + m;
}



complex<double>
LAPWClass::operator()(int ik, int iband, Vec3 r,
		      bool &inTin)
{
  for (int ia=0; ia<AtomPos.size(); ia++) {
    Vec3 disp = r-AtomPos(ia);
    Vec3 min_disp = Lattice.MinImage(disp);
    Vec3 L = disp-min_disp;
    double dist2 = dot(min_disp,min_disp);
    if (dist2 < AtomRadii(ia)*AtomRadii(ia)) {
      inTin = true;
      return evalMuffinTin (ik, iband, ia, min_disp,L);
    }
  }
  inTin = false;
  return evalInterstitial (ik, iband, r);
}

#endif
