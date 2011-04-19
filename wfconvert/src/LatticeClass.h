#ifndef LATTICE_CLASS_H
#define LATTICE_CLASS_H

#include <Common/Blitz.h>

class LatticeClass
{
private:
  // B = Transpose(Inverse(A))
  Mat3 A, B, BBt;
public:
  inline void SetDirect(Mat3 direct);
  inline void SetRecip (Mat3 recip);
  inline Mat3 GetDirect()         const;
  inline Mat3 GetRecip()          const;
  inline Mat3 GetInverse()        const; 
  inline double GetVolume()       const;
  inline Vec3 Twist2k(Vec3 twist) const;
  inline Vec3 k2Twist(Vec3 k)     const;
  inline Vec3 r2u (Vec3 r)        const;
  inline Vec3 u2r (Vec3 r)        const;
  inline Vec3 a(int i)            const;
  inline Vec3 b(int i)            const;
  inline Vec3 MinImage  (Vec3 r)  const;
  inline double MinDist (Vec3 r)  const;
  inline const Mat3& GetBBt() const { return BBt; }
  inline const Mat3& GetB() const { return B; }
};

void
LatticeClass::SetDirect(Mat3 direct)
{
  A = direct;
  B = Transpose(Inverse(direct));
  BBt = B * Transpose(B);
}

void
LatticeClass::SetRecip(Mat3 recip) 
{
  B = (1.0/(2.0*M_PI))*recip;
  A = Inverse(Transpose(B));
  BBt = B * Transpose(B);
}


Mat3
LatticeClass::GetDirect() const
{
  return A;
}

Mat3
LatticeClass::GetRecip() const
{
  return 2.0*M_PI*B;
}

Mat3
LatticeClass::GetInverse() const
{
  return Transpose(B);
}

double
LatticeClass::GetVolume() const
{
  return fabs(det(A));
}

inline Vec3
LatticeClass::Twist2k(Vec3 twist) const
{
  return 2.0*M_PI*Transpose(B)*twist;
}

inline Vec3
LatticeClass::k2Twist(Vec3 k) const
{
  return (1.0/(2.0*M_PI))*(A*k);
}

inline Vec3
LatticeClass::r2u (Vec3 r) const
{
  return B*r;
}

inline Vec3
LatticeClass::u2r(Vec3 u) const
{
  return u*A;
}

inline Vec3
LatticeClass::a(int i) const
{
  return Vec3(A(i,0), A(i,1), A(i,2));
}

inline Vec3
LatticeClass::b(int i) const
{
  return 2.0*M_PI*Vec3(B(i,0), B(i,1), B(i,2));
}

inline Vec3
LatticeClass::MinImage(Vec3 r) const
{
  Vec3 u, rp;
  Vec3 minDisp = r;
  double minDist = dot(r,r);
  
  Vec3 a0 = a(0);
  Vec3 a1 = a(1);
  Vec3 a2 = a(2);

  for (int ix=-1; ix<=1; ix++) 
    for (int iy=-1; iy<=1; iy++) 
      for (int iz=-1; iz<=1; iz++) {
	Vec3 tryDisp = r + (double)ix*a0 + (double)iy*a1 + (double)iz*a2;
	double tryDist = dot (tryDisp, tryDisp);
	if (tryDist < minDist) {
	  minDist = tryDist;
	  minDisp = tryDisp;
	}
      }
  return minDisp;

  u = r2u(r);
  u[0] -= round(u[0]);
  u[1] -= round(u[1]);
  u[2] -= round(u[2]);

  rp = u2r(u);
  return rp;
}

inline double
LatticeClass::MinDist(Vec3 r) const
{
  Vec3 u;
  u = r2u(r);
  u[0] -= round(u[0]);
  u[1] -= round(u[1]);
  u[2] -= round(u[2]);
  r = u2r(u);
  return sqrt(dot(r,r));
}



#endif
