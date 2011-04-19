#ifndef LATTICE_CLASS_H
#define LATTICE_CLASS_H

#include "../Blitz.h"
#include <vector>


class LatticeClass
{
private:
  Mat3 Lattice, LatticeInv;
  // RecipLattice = 2 pi Transpose(Inverse(Lattice))
  Mat3 RecipLattice;
public:
  // Returns i'th lattice vector
  inline Vec3 operator[](int i);
  inline double operator()(int i, int j);
  inline Vec3 a(int i);
  inline Vec3 b(int i);
  inline void Set (Mat3 lattice);
  inline Vec3 Cart2Reduced (Vec3 r);
  inline Vec3 Reduced2Cart (Vec3 r);
  // Returns the radius of the largest inscribed sphere
  inline double rMax();
  inline TinyVector<int,3> MinIndices(double kcut);
  inline vector<Vec3> GenGvecs(Vec3 k, double Ecut);
  inline Vec3 MinImage(Vec3 disp);
  inline double Volume();

  LatticeClass (Mat3 lattice)
  { Set (lattice); }
  LatticeClass()
  { /* do nothing */  }
};

inline Vec3
LatticeClass::a(int i)
{
  return Vec3 (Lattice(i,0), Lattice(i,1), Lattice(i,2));
}

inline Vec3
LatticeClass::b(int i)
{
  return Vec3 (RecipLattice(i,0), RecipLattice(i,1), RecipLattice(i,2));
}

inline Vec3
LatticeClass::operator[](int i)
{
  return Vec3 (Lattice(i,0), Lattice(i,1), Lattice(i,2));
}

inline double
LatticeClass::operator()(int i, int j)
{
  return Lattice(i,j);
}

inline void
LatticeClass::Set(Mat3 lattice)
{
  Lattice = lattice;
  LatticeInv = Inverse(lattice);
  RecipLattice = 2.0*M_PI*Transpose(LatticeInv);
}

inline Vec3
LatticeClass::Cart2Reduced (Vec3 r)
{
  return LatticeInv * r;
}

inline Vec3
LatticeClass::Reduced2Cart (Vec3 u)
{
  return u * Lattice;
}

inline double
LatticeClass::rMax()
{
  double r = 1.0e200;
  Vec3 c = 0.5*((*this)[0] + (*this)[1] + (*this)[2]);
  for (int i=0; i<3; i++) {
    // Vec3 v1 = (*this)[(i+1)%3];
    // Vec3 v2 = (*this)[(i+2)%3];
    
    // double beta1 = (dot(v2,v2)*dot(c,v1) - dot (v1,v2)*dot(c,v2))/
    //	(dot(v1,v1)*dot(v2,v2) - dot(v1,v2)*dot(v1,v2));
    // double beta2 = (dot(v1,v1)*dot(c,v2) - dot (v1,v2)*dot(c,v1))/
    //	(dot(v1,v1)*dot(v2,v2) - dot(v1,v2)*dot(v1,v2));
    
    // Vec3 p = beta1*v1 + beta2*v2;
    // double dist = sqrt(dot(c-p, c-p));
    
    Vec3 A = a(i);
    Vec3 B = a((i+1)%3);
    Vec3 C = a((i+2)%3);
    Vec3 BxC = cross(B,C);
    double dist = 0.5*fabs(dot(A,BxC))/sqrt(dot(BxC,BxC));

    r = min (r, dist);
  }
  return r;
}


// This returns the triple of integers defining the multiples of the
// reciprocal lattice vectors necessary to contain the G-vectors with
// a given reciprocal-space cutoff
TinyVector<int,3>
LatticeClass::MinIndices(double kcut)
{
  Vec3 I;
  for (int i=0; i<3; i++) {
    Vec3 v1 = b((i+1)%3);
    Vec3 v2 = b((i+2)%3);
    Vec3 perpDir = cross(v1,v2);
    perpDir = (1.0/sqrt(dot(perpDir, perpDir)))*perpDir;
    double perpDist = fabs(dot(b(i), perpDir));
    // We add one just do be safe
    I[i] = (int)ceil(kcut/perpDist)+1;
  }
  return I;
}

vector<Vec3>
LatticeClass::GenGvecs (Vec3 k, double Ecut)
{
  double kcut = sqrt(2.0*Ecut);
  TinyVector<int,3> minIndex = MinIndices(kcut);
  TinyVector<int,3> realMax (0,0,0);
  TinyVector<int,3> realMin (0,0,0);


  vector<Vec3> vecs;
  for (int i0=-minIndex[0]; i0<=minIndex[0]; i0++)
    for (int i1=-minIndex[1]; i1<=minIndex[1]; i1++)
      for (int i2=-minIndex[2]; i2<=minIndex[2]; i2++) {
	Vec3 G = (double)i0*b(0) + (double)i1*b(1) + (double)i2*b(2);
	double E = 0.5*dot(G+k,G+k);
	if (E<=Ecut) {
	  vecs.push_back(G);
	  realMax[0]=max(i0,realMax[0]);  realMin[0]=min(i0,realMin[0]);
	  realMax[1]=max(i1,realMax[1]);  realMin[1]=min(i1,realMin[1]);
	  realMax[2]=max(i2,realMax[2]);  realMin[2]=min(i2,realMin[2]);
	}
      }
//   cerr << "realMax = " << realMax << endl;
//   cerr << "realMin = " << realMin << endl;
  return vecs;
}

inline Vec3
LatticeClass::MinImage(Vec3 disp)
{
  Vec3 u = Cart2Reduced (disp);
  u[0] -= round(u[0]);
  u[1] -= round(u[1]);
  u[2] -= round(u[2]);
  return Reduced2Cart(u);
}


inline double
LatticeClass::Volume()
{
  return fabs(det(Lattice));
}

#endif
