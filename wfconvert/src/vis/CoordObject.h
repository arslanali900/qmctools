#ifndef COORD_OBJECT_H
#define COORD_OBJECT_H

#include "GLObject.h"
#include <Common/Blitz.h>

class CoordObject : public GLObject
{
  double Lx, Ly, Lz;
  void POVLine (FILE *fout, 
		double x1, double y1, double z1,
		double x2, double y2, double z2, double radius,
		string rotString);
public:
  void Set (Vec3 box);
  void Set (double lx, double ly, double lz);
  void DrawPOV (FILE *fout, string rotString);
  CoordObject() 
  {
    IsAxes = true;
  }
};

#endif
