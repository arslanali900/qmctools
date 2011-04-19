#ifndef BOX_OBJECT_H
#define BOX_OBJECT_H

#include "GLObject.h"
#include <Common/Blitz.h>

class BoxObject : public GLObject
{
  Vec3 Color;
  double Lx, Ly, Lz;
  Vec3 LatticeVecs[3];
  bool Visible;
  void POVLine (FILE *fout, 
		double x1, double y1, double z1,
		double x2, double y2, double z2, double radius,
		string rotString);
public:
  void Set (Vec3 box, bool useClip=false);
  void Set (double lx, double ly, double lz, bool useClip = false);
  void Set (Mat3 lattice, bool useClip=false);
  void Set (Mat3 lattice, bool visible, bool useClip);
  void SetColor (double red, double blue, double green);
  void DrawPOV (FILE *fout, string rotString);
  BoxObject() : Visible(true)
  {
    Color = Vec3(0.2, 0.2, 0.2);
  }
};

#endif
