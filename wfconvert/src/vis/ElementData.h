#ifndef ELEMENT_COLOR_H
#define ELEMENT_COLOR_H

#include <Common/Blitz.h>

class ElementData
{
private:
  static const int IntColors[109][3];
  static const double Radii[109];
  static const double MaxInv;
public:
  inline static Vec3 GetColor(int element) {
    return Vec3 ((double)IntColors[element-1][0]*MaxInv,
		 (double)IntColors[element-1][1]*MaxInv,
		 (double)IntColors[element-1][2]*MaxInv);
  }
  inline static double GetRadius (int element) {
    return 1.0*1.8897261 * Radii[element-1];
  }
};

#endif



  
