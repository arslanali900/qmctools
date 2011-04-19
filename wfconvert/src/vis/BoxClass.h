#ifndef BOX_CLASS_H
#define BOX_CLASS_H

#include "OnePath.h"

class BoxClass
{
private:
  Vec3 Box, BoxInv;
  Mat3 Lattice, LatticeInv;
  inline void PutInBox(Vec3 &r, int dim);
  bool BreakSegment (Vec3 &r1, Vec3 &r2, Vec3 &wall1, Vec3 &wall2);
  bool BreakSegment (Vec3 &r1, Vec3 &r2, Vec3 &wall1, Vec3 &wall2, int dim);
public:
  inline void PutInBox(Vec3 &r);

  inline void Set (Vec3 box) 
  { 
    double lx = box[0]; double ly=box[1]; double lz=box[2];
    Mat3 lattice;
    lattice =  lx , 0.0, 0.0, 0.0, ly , 0.0, 0.0, 0.0, lz ;
    Set (lattice);
//     Box = box; 
//     BoxInv[0]=1.0/Box[0]; BoxInv[1]=1.0/Box[1]; BoxInv[2]=1.0/Box[2];
  }
  inline void Set (double lx, double ly, double lz)
  {
    Mat3 lattice;
    lattice =  lx , 0.0, 0.0, 0.0, ly , 0.0, 0.0, 0.0, lz ;
    Set (lattice);
//     Box[0] = lx; Box[1] = ly; Box[2] = lz;
//     BoxInv[0]=1.0/lx; BoxInv[1]=1.0/ly; BoxInv[2]=1.0/lz;
  }
  inline void Set (Mat3 lattice)
  {
    Box = Vec3 (lattice(0,0), lattice(1,1), lattice(2,2));
    BoxInv[0]=1.0/Box[0]; BoxInv[1]=1.0/Box[1]; BoxInv[2]=1.0/Box[2];

    Lattice = lattice;
    double vol = det (Lattice);
    double volInv = 1.0/vol;
    LatticeInv(0,0) = volInv * 
      (Lattice(1,1)*Lattice(2,2) - Lattice(2,1)*Lattice(1,2));
    LatticeInv(0,1) = -volInv * 
      (Lattice(1,0)*Lattice(2,2) - Lattice(1,2)*Lattice(2,0));
    LatticeInv(0,2) = volInv * 
      (Lattice(1,0)*Lattice(2,1) - Lattice(1,1)*Lattice(2,0));
    
    LatticeInv(1,0) = -volInv * 
      (Lattice(0,1)*Lattice(2,2) - Lattice(0,2)*Lattice(2,1));
    LatticeInv(1,1) = volInv * 
      (Lattice(0,0)*Lattice(2,2) - Lattice(0,2)*Lattice(2,0));
    LatticeInv(1,2) = -volInv * 
      (Lattice(0,0)*Lattice(2,1) - Lattice(0,1)*Lattice(2,0));
    
    LatticeInv(2,0) = volInv * 
      (Lattice(0,1)*Lattice(1,2) - Lattice(0,2)*Lattice(1,1));
    LatticeInv(2,1) = -volInv * 
      (Lattice(0,0)*Lattice(1,2) - Lattice(0,2)*Lattice(1,0));
    LatticeInv(2,2) = volInv * 
      (Lattice(0,0)*Lattice(1,1) - Lattice(0,1)*Lattice(1,0));
  }

  inline Mat3 GetLattice()
  { return Lattice; }
  
  inline Mat3 GetLatticeInv()
  { return LatticeInv; }

  inline Vec3 operator()(int i)
  {
    return Vec3 (Lattice(i,0), Lattice(i,1), Lattice(i,2));
  }

  inline double operator[](int i) const
  { return Lattice(i,i); }
  inline operator Vec3() const
  { return Box; }

  void PutPathsInBox (vector<OnePath*>& inList);
  void PutInBox (Array<double,1> r);
};

// inline void BoxClass::PutInBox (Vec3 &r)
// {
//   for (int i=0; i<3; i++) {
//     double n = -floor(r[i]*BoxInv[i]+0.5);
//     r[i] += n*Box[i];
//   }
// }

inline void BoxClass::PutInBox (Vec3 &r, int dim)
{
  double n = -floor(r[dim]*BoxInv[dim]+0.5);
  r[dim] += n*Box[dim];
}


inline void BoxClass::PutInBox (Vec3 &r)
{
  Vec3 i, rp;
  i[0] = LatticeInv(0,0)*r[0] + LatticeInv(0,1)*r[1] + LatticeInv(0,2)*r[2];
  i[1] = LatticeInv(1,0)*r[0] + LatticeInv(1,1)*r[1] + LatticeInv(1,2)*r[2];
  i[2] = LatticeInv(2,0)*r[0] + LatticeInv(2,1)*r[1] + LatticeInv(2,2)*r[2];
  i[0] -= round(i[0]);
  i[1] -= round(i[1]);
  i[2] -= round(i[2]);
  rp[0] = Lattice(0,0)*i[0] + Lattice(0,1)*i[1] + Lattice(0,2)*i[2];
  rp[1] = Lattice(1,0)*i[0] + Lattice(1,1)*i[1] + Lattice(1,2)*i[2];
  rp[2] = Lattice(2,0)*i[0] + Lattice(2,1)*i[1] + Lattice(2,2)*i[2];
  
  r = rp;
}



#endif
