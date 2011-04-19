#ifndef ISOSURFACE_H
#define ISOSURFACE_H

#include "GLObject.h"
// #include <gts.h>
#include <Common/Splines/MyTricubicSpline.h>

extern "C" void 
GTS_IsoCartesianWrapper(gdouble **a, GtsCartesianGrid g, guint i,
			gpointer data);

class Isosurface 
{
private:
  MyTricubicSpline &Spline;
  double IsoVal;
  int InterpFactor;
  void IsoCartesian (gdouble **a, GtsCartesianGrid g, guint i);
  friend void GTS_IsoCartesianWrapper(gdouble **a, GtsCartesianGrid g, 
				      guint i, gpointer data);
  GtsSurface *Surface;
  void Construct();
public:
  void DrawGL();

  Isosurface (MyTricubicSpline &spline, double val,
	      int interpFactor=1) :
    Spline(spline), IsoVal(val), InterpFactor(interpFactor),
    Surface(NULL)
  {
    Construct();
  }
};



class IsoSurfaceObject : public GLObject
{
protected:
  

public:
  void SetSurface (MyTricubicSpline &spline, double val);
  
};

#endif
