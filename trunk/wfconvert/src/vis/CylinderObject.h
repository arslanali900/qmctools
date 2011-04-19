#ifndef CYLINDER_OBJECT_H
#define CYLINDER_OBJECT_H

#include "GLObject.h"
#include <Common/Blitz.h>

class CylinderObject : public GLObject
{
protected:
  Vec4 Color;
  Vec3 Pos1, Pos2;
  double Radius, Length;
  static int CylinderListNum, OffScreenListNum;
  static bool CylinderListCreated, OffScreenListCreated;
  bool OffScreen;
  Mat3 Lattice;
  void Set();
  int NumFacets;
public:
  void SetPos (Vec3 pos1, Vec3 pos2);
  void SetBox (Mat3 lattice);
  void SetBox (Vec3 box);
  void SetRadius (double radius);
  void SetColor (Vec3 color);
  void SetColor (Vec4 color);
  void SetNumFacets(int num);
  void DrawPOV (FILE *fout, string rotString);
  CylinderObject(bool offScreen=false);
};


#endif
