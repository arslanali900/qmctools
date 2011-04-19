#include "CubicSpline.h"
#include <cmath>
#include <cstdio>

void TestCubicSpline()
{
  vector<double> data(100), x(100);
  for (int i=0; i<data.size(); i++) {
    x[i] = 2.0*M_PI*i/(double)(data.size()-1);
    data[i] = sin(x[i]);
  }

  SimpleGrid grid;
  grid.Init(x);

  CubSpline sp;
  sp.Init(grid, data);

  for (double t=0.0; t<=1.0; t+=0.001)
    fprintf (stdout, "%20.16e %20.16e %20.16e\n",
	     2.0*M_PI*t, sin (2.0*M_PI*t), sp(2.0*M_PI*t));
}



main()
{
  TestCubicSpline();
}
