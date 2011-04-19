#include "CubicSpline.h"
#include "../IO/IO.h"

inline double sqr (double x)
{ return x*x; }

double
Error (CubicSpline &sp1, CubicSpline &sp2,
       double start, double end)
{
  double err = 0.0;
  double w1 = 1.0/6.0;
  double w2 = 2.0/6.0;
  double w4 = 4.0/6.0;
  
  const int numPoints = 5001;
  double dx = (end - start)/(numPoints-1);
  
  err += w1*(sqr(sp1(start) - sp2(start)) +
	     sqr(sp1(end)   - sp2(end)));
  for (int i=1; i<numPoints; i+=2) {
    double x = start + i*dx;
    err += w4*sqr(sp1(x)-sp2(x));
  }
  for (int i=2; i<numPoints; i+=2) {
    double x = start + i*dx;
    err += w2*sqr(sp1(x)-sp2(x));
  }
  return dx * sqrt(err);
}


void TestUniform(double rmax)
{
  IOSectionClass in;
  assert (in.OpenFile ("FeAtom.h5"));
  assert (in.OpenSection("Grid"));
  Grid *denseGrid = ReadGrid(in);
  in.CloseSection();

  assert (in.OpenSection("RadialWF", 3));
  Array<double,1> uDense;
  assert (in.ReadVar ("u", uDense));
//   for (int i=0; i<denseGrid->NumPoints; i++)
//     uDense(i) /= (*denseGrid)(i);

  CubicSpline denseSpline;
  denseSpline.Init (denseGrid, uDense);

  FILE *fout = fopen ("uniform.dat", "w");
  for (int n=20; n<500; n++) {
    LinearGrid coarseGrid (0.0, rmax, n);
    Array<double,1> uCoarse (n);
    for (int i=0; i<n; i++) 
      uCoarse(i) = denseSpline(coarseGrid(i));
    CubicSpline coarseSpline;
    coarseSpline.Init(&coarseGrid, uCoarse);
    double err = Error (denseSpline, coarseSpline, 0.0, rmax);
    fprintf (fout, "%3d %22.16e\n", n, err);
  }

//   for (double r=0.0; r<rmax; r+= 0.0001) 
//     fprintf (fout, "%22.16e %22.16e\n", r, denseSpline(r));
  fclose (fout);

}


void TestNonuniform(double rmax)
{
  IOSectionClass in;
  assert (in.OpenFile ("FeAtom.h5"));
  assert (in.OpenSection("Grid"));
  Grid *denseGrid = ReadGrid(in);
  in.CloseSection();

  assert (in.OpenSection("RadialWF", 3));
  Array<double,1> uDense;
  assert (in.ReadVar ("u", uDense));
//   for (int i=0; i<denseGrid->NumPoints; i++)
//     uDense(i) /= (*denseGrid)(i);

  FILE *uout = fopen ("u.dat", "w");
  for (int i=0; i<uDense.size(); i++)
    fprintf (uout, "%22.16e %22.16e\n", (*denseGrid)(i), uDense(i));
  fclose (uout);

  CubicSpline denseSpline;
  denseSpline.Init (denseGrid, uDense);

  FILE *fout = fopen ("nonuniform.dat", "w");
  for (int n=20; n<500; n++) {
    double minError = 1e10;
    for (double ratio=2.0; ratio < 100.1; ratio++) {
      OptimalGrid2 coarseGrid (0.0, rmax, ratio, n);
      Array<double,1> uCoarse (n);
      for (int i=0; i<n; i++) 
	uCoarse(i) = denseSpline(coarseGrid(i));
      CubicSpline coarseSpline;
      coarseSpline.Init(&coarseGrid, uCoarse);
      double err = Error (denseSpline, coarseSpline, 0.0, rmax);
      minError = min (minError, err);
    }
    fprintf (fout, "%3d %22.16e\n", n, minError);
  }

//   for (double r=0.0; r<rmax; r+= 0.0001) 
//     fprintf (fout, "%22.16e %22.16e\n", r, denseSpline(r));
  fclose (fout);

}





main()
{
  TestUniform(5.0);
  TestNonuniform(5.0);
}
