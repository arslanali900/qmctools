#ifndef COLOR_MAP_H
#define COLOR_MAP_H

#include <Common/Splines/CubicBspline.h>

typedef enum { AUTUMN, BONE, COLORCUBE, COOL, COPPER, FLAG,
               GRAY, HOT, HSV, JET, LINES, PINK, PRISM, SPRING, 
	       SUMMER, WHITE, WINTER, BLUE_WHITE_RED } ColorMapType;



class ColorMap
{
private:
  CubicBspline Splines[4];
  bool Initialized;
  double Min, Max;
  static bool MapsRead;
  void ReadMaps();
public:
  inline void operator()(double x, TinyVector<double,3> &val);
  inline void operator()(double x, TinyVector<double,4> &val);
  void Init (double min, double max, ColorMapType map=BLUE_WHITE_RED);
  ColorMap() : Initialized(false)
  {
    if (!MapsRead)
      ReadMaps();
  }
};

inline void
ColorMap::operator() (double x, TinyVector<double,3> &val)
{
  x = max (Min+1.0e-6, x);
  x = min (Max-1.0e-6, x);
  val[0] = max(0.0, min(Splines[0](x), 1.0));
  val[1] = max(0.0, min(Splines[1](x), 1.0));
  val[2] = max(0.0, min(Splines[2](x), 1.0));
}

inline void
ColorMap::operator() (double x, TinyVector<double,4> &val)
{
  if (Initialized) {
    x = max (Min+1.0e-6, x);
    x = min (Max-1.0e-6, x);
    val[0] = max(0.0, min(Splines[0](x), 1.0));
    val[1] = max(0.0, min(Splines[1](x), 1.0));
    val[2] = max(0.0, min(Splines[2](x), 1.0));
    val[3] = max(0.0, min(Splines[3](x), 1.0));
  }
  else
    val = 0.0;
}


#endif
