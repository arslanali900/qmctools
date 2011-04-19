#ifndef SPHERE_OBJECT_H
#define SPHERE_OBJECT_H

#include "GLObject.h"
#include <Common/Blitz.h>

class SphereObject : public GLObject
{
protected:
  Vec4 Color;
  Vec3 Pos;
  double Radius;
  static int SphereListNum, OffScreenListNum;
  static bool SphereListCreated, OffScreenListCreated;
  bool OffScreen;
  Mat3 Lattice;
  void Set();
  int NumFacets;
public:
  void SetPos (Vec3 pos);
  void SetBox (Mat3 lattice);
  void SetBox (Vec3 box);
  void SetRadius (double radius);
  void SetColor (Vec3 color);
  void SetColor (Vec4 color);
  void SetNumFacets(int num);
  void DrawPOV (FILE *fout, string rotString);
  SphereObject(bool offScreen=false);
};


#endif
