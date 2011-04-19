#ifndef PLANE_OBJECT_H
#define PLANE_OBJECT_H

#include "GLObject.h"
#include <GL/glut.h>
#include <Common/Blitz.h>
#include "ColorMap.h"
#include <Common/Splines/MyTricubicSpline.h>

class PlaneObject : public GLObject
{
private:
  MyTricubicSpline &Spline;
  int Direction;
  double Position;
  ColorMap CMap;
  ColorMapType MapType;
  double MinVal, MaxVal;
  GLuint TextureNum;
  bool HaveTexture, BuiltTexture;
  bool IsInitialized;
  Mat3 Lattice;
  // Information for contours
  bool UseContours;
  static int EdgeTable[4][5];
  static int EdgeData[16][5];
  Array<double,2> ValData;
  Vec3 FindEdge (int ix, int iy, int edgeNum,
		 Vec3 u0, Vec3 s, Vec3 t,
		 double isoVal);
  Vec3 uMin, uMax, uCenter;
public:
  void Set();
  void Init();
  void DrawPOV (FILE *out, string rotMatrix);
  void SetPosition (int dir, double pos);
  void SetIsocontours (bool show);
  void SetLattice(Mat3 lattice);
  void SetColorMap (ColorMapType map);
  void SetCenter (Vec3 ucenter, Vec3 umin, Vec3 umax);
  PlaneObject (MyTricubicSpline &spline) : 
    Spline (spline), Direction(0), Position(0),
    HaveTexture(false), BuiltTexture(false),
    IsInitialized(false), MapType (BLUE_WHITE_RED)
  {
    Lattice = 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0;
  }
  ~PlaneObject();
};

#endif
