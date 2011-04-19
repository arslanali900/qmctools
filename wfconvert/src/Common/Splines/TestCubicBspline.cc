#include "CubicBspline.h"
#include "FastCubicSpline.h"

void TestValue()
{
  Array<double,1> d(10), e(11);
  for (int i=0; i<d.size(); i++) {
    double x = 2.0*M_PI*i/(d.size());
    d(i) = sin(x);
    e(i) = d(i);
  }
  e(10) = e(0);

  CubicBspline bspline;
  bspline.Init (0.0, 2.0*M_PI, d, true, 
		BoundaryCondition<double>(FLAT), 
		BoundaryCondition<double>(FLAT));

  FastCubicSpline spline;
  spline.Init (0.0, 2.0*M_PI, e, true);

  for (double x=0.0; x<2.0*M_PI; x+=0.0001) 
    fprintf (stdout, "%20.16e %20.16e %20.16e\n", x, bspline.Deriv(x), spline(x));
}

main() 
{
  TestValue();

}
    
    
