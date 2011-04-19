#ifndef GAUSSIAN_H
#define GAUSSIAN_H

#include <vector>
#include <string>
#include <map>
using namespace std;
#include <Common/Blitz.h>
#include "ParserClass.h"
#include "config.h"

#define F77_DGEMM  F77_FUNC(dgemm,DGEMM)

extern "C" void 
F77_DGEMM (char *transA, char *transB, int *m, int *n, int *k,
	   double *alpha, const double *A, int *lda, const double *B, int *ldb,
	   double *beta,  double *C, int *ldc);

/* From GAMESS pltorb.code 

C     ----- NORMALIZATION FACTOR FOR A GAUSSIAN -----
C     GENERAL NORMALIZATION FORMULA FOR X**L,Y**M,Z**N
C     CARTESIAN GAUSSIAN IS: (FACT2=DOUBLE FACTORIAL)
C        A=ZETA**(L+M+N+1.5)
C        B=2**(2*L+2*M+2*N+1.5)
C        C=FACT2(L+L-1)*FACT2(M+M-1)*FACT2(N+N-1)*PITO32
C

Note:  A and B are applied below at the level of the Shell,
since they depend only on (L+M+N)=l, but C must be applied at
the level of Site, since they depend on L,M,N separately.

*/ 

class Shell
{
  vector<blitz::TinyVector<double,3> > Contraction;

public:
  int l;

  inline void resize(int size)
  { Contraction.resize(size); }

  inline void clear()
  { Contraction.clear(); }

  inline double operator()(double r2)
  {
    double val=0.0;
    vector<blitz::TinyVector<double,3> >::iterator iter;
    for (iter=Contraction.begin(); iter!=Contraction.end(); iter++)
      val += (*iter)[0] * exp(-(*iter)[1]*r2) * (*iter)[2];
    return val;
  };

  inline void push_back (blitz::TinyVector<double,2> prim)
  { 
    double c, zeta, norm;
    c = prim[0];
    zeta = prim[1];
    norm = sqrt(pow(zeta,1.5+l) * pow(2.0, 2.0*l+1.5));
    Contraction.push_back(TinyVector<double,3>(c,zeta,norm)); 
  }

  inline double nrm()
  {
    vector<blitz::TinyVector<double,3> >::iterator iter;
    double n = 0.0;
    for (iter=Contraction.begin(); iter!=Contraction.end(); iter++)
      n += (*iter)[0]*(*iter)[0];
    return n;
  }


  inline blitz::TinyVector<double,2> operator[](int i) const
  {    
    return TinyVector<double,2> (Contraction[i][0],
				 Contraction[i][1]);
  }

//   inline blitz::TinyVector<double,2>& operator[](int i)
//   {    return Contraction[i];  }

  inline int size() { return Contraction.size(); }
    

  bool ReadGAMESS (ParserClass &log);
};



class Site
{
private:
  blitz::TinyVector<double,3> pos;
  vector<Shell> Shells;
  vector<double> norm;
  inline double double_factorial(int n) {
    double dn = n;
    double f = 1.0;
    while (dn > 0.0) {
      f *= dn;
      dn -= 1.0;
    }
    return f;
  }
  
  double snorm;
  double pnorm;
  double dnorm[6], fnorm[10], gnorm[15];
    
public:
  Site() {
    snorm = sqrt(pow(M_PI,-1.5));
    pnorm = snorm;
    dnorm[0] = snorm/sqrt(3.0);
    dnorm[1] = snorm/sqrt(3.0);
    dnorm[2] = snorm/sqrt(3.0);
    dnorm[3] = snorm;
    dnorm[4] = snorm;
    dnorm[5] = snorm;

    fnorm[0] = snorm/sqrt(15.0);
    fnorm[1] = snorm/sqrt(15.0);
    fnorm[2] = snorm/sqrt(15.0);
    fnorm[3] = snorm/sqrt(3.0);
    fnorm[4] = snorm/sqrt(3.0);
    fnorm[5] = snorm/sqrt(3.0);
    fnorm[6] = snorm/sqrt(3.0);
    fnorm[7] = snorm/sqrt(3.0);
    fnorm[8] = snorm/sqrt(3.0);
    fnorm[9] = snorm;
  }

  inline void set_pos (TinyVector<double,3> r)
  { pos = r; }

  inline void push_back (Shell &shell)
  { 
    Shells.push_back(shell); 
    
    
  }

  inline void 
  operator()(blitz::TinyVector<double,3> r,
	     vector<double>::iterator &iter)
  {
    blitz::TinyVector<double,3> dr = r - pos;
    double x = dr[0];
    double y = dr[1];
    double z = dr[2];
    double r2 = dot(dr,dr);

    for (int ish=0; ish<Shells.size(); ish++) {
      Shell &sh = Shells[ish];
      double u = sh(r2);
      switch (sh.l) {
	case 0:
	  *iter++ = snorm*u;
	  break;
	case 1:
	  *iter++ = pnorm*x*u;
	  *iter++ = pnorm*y*u;
	  *iter++ = pnorm*z*u;
	  break;
	case 2:
	  *iter++ = dnorm[0]*x*x*u;
	  *iter++ = dnorm[1]*y*y*u;
	  *iter++ = dnorm[2]*z*z*u;
	  *iter++ = dnorm[3]*x*y*u;
	  *iter++ = dnorm[4]*x*z*u;
	  *iter++ = dnorm[5]*y*z*u;
	  break;
	case 3:
	  *iter++ = fnorm[0]*x*x*x*u;
	  *iter++ = fnorm[1]*y*y*y*u;
	  *iter++ = fnorm[2]*z*z*z*u;
	  *iter++ = fnorm[3]*x*x*y*u;
	  *iter++ = fnorm[4]*x*x*z*u;
	  *iter++ = fnorm[5]*y*y*x*u;
	  *iter++ = fnorm[6]*y*y*z*u;
	  *iter++ = fnorm[7]*z*z*x*u;
	  *iter++ = fnorm[8]*z*z*y*u;
	  *iter++ = fnorm[9]*x*y*z*u;
	  break;
        case 4:
	  // CHECK ORDERING and SIZE!
	  *iter++ = gnorm[ 0]*x*x*x*x*u;
	  *iter++ = gnorm[ 1]*y*y*y*y*u;
	  *iter++ = gnorm[ 2]*z*z*z*z*u;
	  *iter++ = gnorm[ 3]*x*x*x*y*u;
	  *iter++ = gnorm[ 4]*x*x*x*z*u;
	  *iter++ = gnorm[ 5]*x*x*y*y*u;
	  *iter++ = gnorm[ 6]*x*x*y*z*u;
	  *iter++ = gnorm[ 7]*x*x*z*z*u;
	  *iter++ = gnorm[ 8]*x*y*y*y*u;
	  *iter++ = gnorm[ 9]*x*y*y*z*u;
	  *iter++ = gnorm[10]*x*y*z*z*u;
	  *iter++ = gnorm[11]*x*z*z*z*u;
	  *iter++ = gnorm[12]*y*y*y*z*u;
	  *iter++ = gnorm[13]*y*y*z*z*u;
	  *iter++ = gnorm[14]*y*z*z*z*u;
	  break;

	default:
	  cerr << "Error:  not implemented above l=4.\n";
	  abort();
	}
    }
  }

  inline void 
  operator()(blitz::TinyVector<double,3> r, double* &iter)
  {
    blitz::TinyVector<double,3> dr = r - pos;
    double x = dr[0];
    double y = dr[1];
    double z = dr[2];
    double r2 = dot(dr,dr);

    for (int ish=0; ish<Shells.size(); ish++) {
      Shell &sh = Shells[ish];
      double u = sh(r2);
      switch (sh.l) {
	case 0:
	  *iter++ = snorm*u;
	  break;
	case 1:
	  *iter++ = pnorm*x*u;
	  *iter++ = pnorm*y*u;
	  *iter++ = pnorm*z*u;
	  break;
	case 2:
	  *iter++ = dnorm[0]*x*x*u;
	  *iter++ = dnorm[1]*y*y*u;
	  *iter++ = dnorm[2]*z*z*u;
	  *iter++ = dnorm[3]*x*y*u;
	  *iter++ = dnorm[4]*x*z*u;
	  *iter++ = dnorm[5]*y*z*u;
	  break;
	case 3:
	  *iter++ = fnorm[0]*x*x*x*u;
	  *iter++ = fnorm[1]*y*y*y*u;
	  *iter++ = fnorm[2]*z*z*z*u;
	  *iter++ = fnorm[3]*x*x*y*u;
	  *iter++ = fnorm[4]*x*x*z*u;
	  *iter++ = fnorm[5]*y*y*x*u;
	  *iter++ = fnorm[6]*y*y*z*u;
	  *iter++ = fnorm[7]*z*z*x*u;
	  *iter++ = fnorm[8]*z*z*y*u;
	  *iter++ = fnorm[9]*x*y*z*u;
	  break;
        case 4:
	  // CHECK ORDERING and SIZE!
	  *iter++ = gnorm[ 0]*x*x*x*x*u;
	  *iter++ = gnorm[ 1]*y*y*y*y*u;
	  *iter++ = gnorm[ 2]*z*z*z*z*u;
	  *iter++ = gnorm[ 3]*x*x*x*y*u;
	  *iter++ = gnorm[ 4]*x*x*x*z*u;
	  *iter++ = gnorm[ 5]*x*x*y*y*u;
	  *iter++ = gnorm[ 6]*x*x*y*z*u;
	  *iter++ = gnorm[ 7]*x*x*z*z*u;
	  *iter++ = gnorm[ 8]*x*y*y*y*u;
	  *iter++ = gnorm[ 9]*x*y*y*z*u;
	  *iter++ = gnorm[10]*x*y*z*z*u;
	  *iter++ = gnorm[11]*x*z*z*z*u;
	  *iter++ = gnorm[12]*y*y*y*z*u;
	  *iter++ = gnorm[13]*y*y*z*z*u;
	  *iter++ = gnorm[14]*y*z*z*z*u;
	  break;

	default:
	  cerr << "Error:  not implemented above l=4.\n";
	  abort();
	}
    }
  }

  inline int size()
  {
    int n_of_l[] = {1, 3, 6, 10, 15};
    int s = 0;
    for (int is=0; is<Shells.size(); is++) {
      int l = Shells[is].l;
      int n = n_of_l[l];
      s += n;
//       fprintf (stderr, "is=%d l=%d num=%d\n",
// 	       is, l, n);
    }
    return s;
  }

  inline int num_orbs()
  {
    int n_of_l[] = {1, 3, 5, 7, 9};
    int s = 0;
    for (int is=0; is<Shells.size(); is++) {
      int l = Shells[is].l;
      int n = n_of_l[l];
      s += n;
    }
    return s;
  }


};



class GaussianBasis
{
  vector<Site> Sites;
public:
  inline void operator()(blitz::TinyVector<double,3> r, 
			 vector<double> &basis)
  {
    vector<double>::iterator iter = basis.begin();
    for (int isite=0; isite<Sites.size(); isite++)
      Sites[isite](r, iter);
  }

  inline void operator()(vector<blitz::TinyVector<double,3> > &r,
			 Array<double,2> &basis)
  {
    for (int ir=0; ir<r.size(); ir++) {
      double *b = &(basis(ir,0));
      for (int isite=0; isite<Sites.size(); isite++)
	Sites[isite](r[ir], b);
    }
  }

  inline int size()
  {
    int s=0;
    for (int isite=0; isite<Sites.size(); isite++)
      s += Sites[isite].size();
    return s;
  }

  inline int num_orbs()
  {
    int s=0;
    for (int isite=0; isite<Sites.size(); isite++)
      s += Sites[isite].num_orbs();
    return s;
  }


  inline void push_back(Site &site)
  { Sites.push_back (site); }
};
	

struct Atom
{
  Vec3 pos;
  double charge;
  string name;
  int Z;
};

class GaussianOrbitalSet
{
private:
  vector<double> BasisVec;
  // First index is position, second indexes basis element
  Array<double,2> BasisArray;
  vector<double> Eigenvalues;
  vector<Atom> Atoms;
  string Functional;
  double TotalEnergy;
  map<int, string> ZToSymbolMap;
  map<string, int> SymbolToZMap;
  map<int,double> ZToMassMap;
  int NumElecs, NumSpins;
  TinyVector<int,2> SpinElecs;
public:
  GaussianOrbitalSet();
  // First index is orbital num.  Second index is basis element
  TinyVector<blitz::Array<double,2>,2> Coefs;
  GaussianBasis Basis;

  
  inline void set_num_orbs (int numorb) {
    int nb = Basis.size();
    BasisVec.resize(nb);
    Coefs[0].resize(numorb, nb);
    if (NumSpins == 2)
      Coefs[1].resize(numorb, nb);
    Eigenvalues.resize(numorb);
  }

  inline void
  operator()(int ispin,
	     blitz::TinyVector<double,3> r,
	     vector<double> &vals)
  {
    blitz::Array<double,2> &C = Coefs[ispin];
    Basis(r, BasisVec);
    // Call DGEMV for vals = C * BasisVec
    for (int iorb=0; iorb<C.extent(0); iorb++) {
      vals[iorb] = 0.0;
      for (int ib=0; ib<C.extent(1); ib++)
	vals[iorb] += C(iorb,ib)*BasisVec[ib];
    }      
  }

  inline void
  operator()(int ispin,
	     vector<blitz::TinyVector<double,3> > &r,
	     Array<double,2> &vals)
  {
    blitz::Array<double,2> &C = Coefs[ispin];
    int numr = r.size();
    int numorb = C.extent(0);
    int numbasis = C.extent(1);
      
    if (BasisArray.extent(0) < numr)
      BasisArray.resize(numr, numbasis);
    Basis(r, BasisArray);

    // Call DGEMV for vals = C * BasisVec
    char transA = 'T';
    char transB = 'N';
    double alpha = 1.0;
    double beta = 0.0;
    F77_DGEMM (&transA, &transB, &numorb, &numr, &numbasis, 
	       &alpha, C.data(), &numbasis, 
	       BasisArray.data(), &numbasis, &beta, vals.data(), &numorb);

//     for (int ir=0; ir<numr; ir++)
//       for (int iorb=0; iorb<Coefs.extent(0); iorb++) {
// 	vals(ir,iorb) = 0.0;
// 	for (int ib=0; ib<C.extent(1); ib++)
// 	  vals(ir,iorb) += C(iorb,ib) * BasisArray(ir,ib);
//       }
  }      



  bool ReadGAMESS (string logname);
  bool ReadGAMESS (string logname, string vecname);
  void Write_ESHDF (string name);
};


#endif
