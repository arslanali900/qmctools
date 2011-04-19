#ifndef ISOSURFACE_H
#define ISOSURFACE_H

#include <Common/Splines/MyTricubicSpline.h>
//#include <Common/Splines/TricubicBspline.h>
#include <GLObject.h>
#include <vector>
#include <Common/Blitz.h>


/// This class uses to Marching Cubes algorithm to construct an
/// isosurface from data tabulated on a 3D grid.  It inherits from
/// both MyTricubicSpline and GLObject.  You simply set up the
/// tricubic spline as usual and then call SetIsoval to compute set up
/// the isosurface with the Marching Cubes algorithm.
class Isosurface : public MyTricubicSpline, public GLObject
{
private:
  static int EdgeData[256][16];
  static int EdgeTable[12][7];
  /// Returns the number of real solutions
  int CubicFormula (double a, double b, double c, double d,
		     double &x1, double &x2, double &x3);
  Vec3 FindEdge(int ix, int iy, int iz, int edgeNum);
  bool UseCubicInterp;
  bool UseNormals;
  vector<double> Isovals;
  double Isoval;
  TinyVector<double,4> Color;
  vector<TinyVector<double,4> > Colors;
  double Alpha;
  void Set();
  Vec3 LatticeVecs[3];
  Mat3 Lattice, LatticeInv;
  Vec3 uCenter, uMin, uMax;
public:
  inline void SetLattice (Mat3 lattice);
  int NumTriangles(int i=0);
  inline void SetIsoval (double val) 
  { Isovals.resize(1);  Isovals[0] = val; Set();}
  inline void SetIsoval (vector<double> vals) 
  { Isovals.resize(vals.size()); Isovals = vals; Set(); }
  inline void SetCenter (Vec3 ucenter, Vec3 umin, Vec3 umax)
  { uCenter = ucenter;  uMin = umin;   uMax = umax; }
  void DrawPOV(FILE* out, string rotMatrix);
  void SetColor (Vec3 color);
  void SetColor (vector<TinyVector<double,3> > &colors);
  void SetAlpha (double alpha);
  inline void SetColor (double r, double g, double b) 
  { SetColor (TinyVector<double,3>(r,g,b)); }
  Isosurface() : UseCubicInterp(true), UseNormals(true), 
		 uCenter(0.0, 0.0, 0.0),
		 uMin(0.0, 0.0, 0.0), uMax(1.0, 1.0, 1.0)
  {
    Color[0] = 0.0; Color[1] = 0.8; Color[2] = 0.0; Color[3] = 0.3;
    Colors.push_back(Color);
    Alpha = 0.3;
    LatticeVecs[0] = 1.0, 0.0, 0.0;
    LatticeVecs[1] = 0.0, 1.0, 0.0;
    LatticeVecs[2] = 0.0, 0.0, 1.0;
  }
};

inline void
Isosurface::SetLattice (Mat3 lattice)
{
  Lattice = lattice;
  LatticeInv = Inverse(lattice);
  LatticeVecs[0] = Vec3 (lattice(0,0), lattice(0,1), lattice(0,2));
  LatticeVecs[1] = Vec3 (lattice(1,0), lattice(1,1), lattice(1,2));
  LatticeVecs[2] = Vec3 (lattice(2,0), lattice(2,1), lattice(2,2));
}

#endif
