#ifndef DISK_OBJECT_H
#define DISK_OBJECT_H

#include "GLObject.h"
#include <Common/Blitz.h>

class DiskObject : public GLObject
{
protected:
  Vec3 Color, Pos, NormVec;
  int Axis;
  double Radius;
  static int DiskListNum, OffScreenListNum;
  static bool DiskListCreated, OffScreenListCreated;
  void Set();
  bool OffScreen;
public:
  void SetPos (Vec3 pos);
  void SetNormVec (Vec3 norm);
  void SetAxis (int axis);
  void SetRadius (double radius);
  void SetColor (Vec3 color);
  void DrawPOV (FILE *fout, string rotString);
  DiskObject(bool offScreen=false);
};


#endif
