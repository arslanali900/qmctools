#include "BsplineHelper.h"
#include "CubicNUBspline.h"
#include "Grid.h"

void TestBasis()
{
  LinearGrid grid(-1.0, 1.0, 11);
  NUBsplineBasis<LinearGrid> uniform;
  uniform.Init(&grid);

  FILE *fout = fopen ("uniformBasis.dat", "w");
  TinyVector<double,4> bfuncs;
  for (double x=0.0; x<=1.0; x+=0.0001) {
    int i = uniform((double)x, bfuncs);
    fprintf (fout, "%22.16e %22.16e %22.16e %22.16e %22.16e\n",
	     x, bfuncs[0], bfuncs[1], bfuncs[2], bfuncs[3]);
  }
  fclose(fout);
}

void TestNUBasis()
{
  Array<double,1> points(11);
  points = 0.0, 1.0, 1.4, 1.8, 6.0, 7.0, 9.0, 10.0, 12.8, 13.0, 15.0;
  GeneralGrid grid;
  grid.Init(points);
  NUBsplineBasis<GeneralGrid> nonuni;
  nonuni.Init (&grid);

  FILE *fout = fopen ("NUBasis.dat", "w");
  TinyVector<double,4> bfuncs, dbfuncs, d2bfuncs;
  for (double x=0.0; x<=15.0; x+=0.001) {
//   for (int i=0; i<grid.NumPoints; i++) {
//     double x = grid(i);
    nonuni(x, bfuncs, dbfuncs, d2bfuncs);
    fprintf (fout, "%22.16e %22.16e %22.16e %22.16e %22.16e "
	     "%22.16e %22.16e %22.16e %22.16e "
	     "%22.16e %22.16e %22.16e %22.16e\n",
	     x, bfuncs[0], bfuncs[1], bfuncs[2], bfuncs[3],
	     dbfuncs[0], dbfuncs[1], dbfuncs[2], dbfuncs[3],
	     d2bfuncs[0], d2bfuncs[1], d2bfuncs[2], d2bfuncs[3]);
  }
  fclose(fout);
}

void
TestNUBspline()
{
  Array<double,1> points(11);
  points = 0.0, 1.0, 1.4, 1.8, 6.0, 7.0, 9.0, 10.0, 12.8, 13.0, 15.0;
  // points = 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0;
  GeneralGrid grid;
  grid.Init(points);
  Array<double,1> data(10);
  data = sin (points(Range(0,9)));
  
//   BoundaryCondition<double> lBC(NATURAL);
//   BoundaryCondition<double> rBC(FLAT);
  BoundaryCondition<double> lBC(PERIODIC);
  BoundaryCondition<double> rBC(PERIODIC);
  CubicNUBspline<GeneralGrid> spline;
  spline.Init (grid, data, lBC, rBC);

  FILE *fout = fopen ("NUBspline.dat", "w");
  for (double x=0.0; x<=15.0; x+=0.001) {
    double y = spline (x);
    double dy = spline.Deriv(x);
    double d2y = spline.Deriv2(x);
    double ex = sin(x);
    fprintf (fout, "%20.16e %20.16e %20.16e %20.16e %20.16e\n", x, y, dy, d2y, ex);
  }
  fclose (fout);

}


main()
{
  TestBasis();
  TestNUBasis();
  TestNUBspline();
}
