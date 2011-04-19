#ifndef FIT_SK_H
#define FIT_SK_H

#include <Common/PlaneWavePHDFT/GVecs.h>
#include <Common/Ewald/OptimizedBreakup.h>
#include <vector>


class HermiteBasisClass
{
private:
  TinyMatrix<double,3,6> S;
  int NumKnots;
  double delta, deltaInv, k_c;
public:
  inline void SetNumKnots (int num);
  inline void SetCutoff (double kcut);
  inline double operator() (int i, double k);
  // returns the definite integral of the nth basis function
  // from 0 to k;
  inline double Integral (int n, double k);
  inline int NumElements() 
  { return 3 * NumKnots; }
  HermiteBasisClass() : NumKnots(0), delta(0.0)
  {
    S = 
      1.0,   0.0,   0.0, -10.0,  15.0,  -6.0,
      0.0,   1.0,   0.0,  -6.0,   8.0,  -3.0,
      0.0,   0.0,   0.5,  -1.5,   1.5,  -0.5;
  }
};


template<typename T>
class RecipVal
{
public:
  double Gmag;
  double Degeneracy;
  T ValSum;
};

template<typename T>
class RecipSet
{
public:
  int ParticleType1, ParticleType2;
  vector<RecipVal<T> > Vals;
  void AddVal (double Gmag, T val);
};

class SkClass
{
private:
  vector<RecipSet<complex<double> > > RhoSets;
  HermiteBasisClass Basis;
  // First index is the set number.  Second index is basis func
  // number.
  Array<double,2> FitCoefs;
  double BoxVol;
  double kcut, width;
  double BlendFunc (double k);
public:
  vector<RecipSet<double> > SkSets;
  TinyMatrix<double,3,3> PrimLattice;
  TinyVector<int,3> Supercell;
  GVecsClass GVecs;
  // The first index is the set number.  The second index is the
  // G-vector number
  Array<double,2> Sk;
  Array<complex<double>,2> Rhok;
  int NumParticles;
  bool Read (string fname);
  void Fit (int numKnots, double k_c);
  template<int N> void FitNonlinear();
  void FitRon (int numStars);
  void FitRon2 (double kcut);
  double Sum (int set, int firstk, int lastk);
  double TotalSum (int firstk, int lastk);
  double Madelung();

};


inline double
HermiteBasisClass::operator()(int n, double k)
{
  int i=n/3;
  int alpha = n-3*i;
  double ka = delta*(i-1);
  double kb = delta*i;
  double kc = delta*(i+1);
  kc = min(k_c, kc);
  if ((k > ka) && (k <= kb)) {
    double sum = 0.0;
    double prod = 1.0;
    for (int j=0; j<=5; j++) {
      sum += (S(alpha,j) * prod);
      prod *= ((kb - k) * deltaInv);
    }
    for (int j=0; j<alpha; j++)
      sum *= -1.0;
    return (sum);
  }
  else if ((k > kb) && (k <= kc)) {
    double sum = 0.0;
    double prod = 1.0;
    for (int j=0; j<=5; j++) {
      sum += S(alpha,j) * prod;
      prod *= ((k-kb) * deltaInv);
    }
    return sum;
  }
  return 0.0;
}

inline double
HermiteBasisClass::Integral(int n, double k)
{
  int i=n/3;
  int alpha = n-3*i;
  double ka = delta*(i-1);
  double kb = delta*i;
  double kc = delta*(i+1);
  kc = min(k_c, kc);
  if ((k > ka) && (k <= kb)) {
    double sum = 0.0;
    double prod = (kb-k);
    double prodka = (kb-ka);
    for (int j=0; j<=5; j++) {
      double jInv = 1.0/(double)(j+1);
      sum -= jInv* S(alpha,j) * (prod-prodka);
      prod   *= ((kb - k)  * deltaInv);
      prodka *= ((kb - ka) * deltaInv);
    }
    for (int j=0; j<alpha; j++)
      sum *= -1.0;
    return (sum);
  }
  else if ((k > kb) && (k <= kc)) {
    double sum = 0.0;
    double prod = (k-kb);
    for (int j=0; j<=5; j++) {
      double jInv = 1.0/(double)(j+1);
      sum += jInv * S(alpha,j) * prod;
      prod *= ((k-kb) * deltaInv);
    }
    if (i > 0)
      sum += Integral (n, kb-1.0e-14);
    return sum;
  }
  else if (k > kc)
    return Integral (n, kc-1.0e-14);
  return 0.0;
}


inline void
HermiteBasisClass::SetNumKnots (int num)
{
  assert (num > 1);
  NumKnots = num;
  if (k_c != 0.0) {
    delta = k_c / (NumKnots - 1);
    deltaInv = 1.0/delta;
  }
}

inline void
HermiteBasisClass::SetCutoff (double kcut)
{
  k_c = kcut;
  if (NumKnots != 0) {
    delta = k_c / (NumKnots - 1);
    deltaInv = 1.0/delta;
  }
}


#endif
