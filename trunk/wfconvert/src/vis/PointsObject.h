#ifndef POINTS_OBJECT_H
#define POINTS_OBJECT_H

#include "GLObject.h"
#include <Common/Blitz.h>

class PointsObject : public GLObject
{
protected:
  Vec4 Color;
  double Radius;
  bool OffScreen;
  Mat3 Lattice;
  void Set();
  int NumFacets;
public:
  void SetPoints (Array<Vec3,2> &points);
  void SetBox (Mat3 lattice);
  void SetBox (Vec3 box);
  void SetRadius (double radius);
  void SetColor (Vec3 color);
  void SetColor (Vec4 color);
  void DrawPOV (FILE *fout, string rotString);
  PointsObject(bool offScreen=false);
};


#endif
