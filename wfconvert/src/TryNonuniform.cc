#include <iostream>
#include <Common/Splines/TricubicBspline.h>
#include <Common/Splines/TricubicNUBspline.h>
#include <Common/IO/IO.h>
#include <cassert>

using namespace std;

using namespace IO;

void ReadData (TricubicBspline<complex<double> > &exactOrb)
{
  IOSectionClass in;
  //  assert (in.OpenFile("/home/kesler/BoronNitride/CASINO/ABINIT/BN64/BN64_1loc.h5"));
  assert (in.OpenFile("/home/kesler/MgO/MgO64/MgO64_oneloc.h5"));
  //  assert (in.OpenFile("/home/kesler/Silicon/Si64_loc.h5"));
  //  assert (in.OpenFile("/home/kesler/Aluminum/Al32_loc.h5"));
  //  assert (in.OpenFile("/home/kesler/Aluminum/Al108_oneloc.h5"));
  //assert (in.OpenFile("/home/kesler/FeO/rocksalt/FeO8_rocksalt_trunc.h5"));
  //assert (in.OpenFile("/home/kesler/FeO/rocksalt/FeO8_rocksalt_bigloc.h5"));
  //assert (in.OpenFile ("/home/kesler/Silicon/silicon64_loc.h5"));
  assert (in.OpenSection("parameters"));
  Array<double,2> lattice;
  assert (in.ReadVar("lattice", lattice));
  in.CloseSection (); // "parameters"

  Array<double,4> realVals;
  assert (in.OpenSection ("eigenstates"));
  assert (in.OpenSection ("twist"));
  assert (in.OpenSection ("band", 0));
  assert (in.ReadVar("eigenvector", realVals));
  Array<complex<double>,3> orb;

  orb.resize(realVals.extent(0)-1, realVals.extent(1)-1, realVals.extent(2)-1);
  for (int ix=0; ix<orb.extent(0); ix++)
    for (int iy=0; iy<orb.extent(1); iy++)
      for (int iz=0; iz<orb.extent(2); iz++)
	orb(ix,iy,iz) = complex<double>(realVals(ix,iy,iz,0),
					realVals(ix,iy,iz,1));

  exactOrb.Init (0.0, lattice(0,0), 0.0, lattice(1,1), 0.0, lattice(2,2),
		 orb, true, NATURAL, NATURAL, NATURAL);
  Array<double,1> xCurvature(orb.extent(0));
  xCurvature = 0.0;
  TinyVector<complex<double>,3> grad;
  TinyMatrix<complex<double>,3,3> secDerivs;
  int Nx=orb.extent(0); int Ny=orb.extent(1); int Nz=orb.extent(2);
  double deltax = lattice(0,0)/(double)(Nx);
  double deltay = lattice(1,1)/(double)(Ny);
  double deltaz = lattice(2,2)/(double)(Nz);
  complex<double> val, lapl;
  double prefactor = 1.0/(double)(Nx*Ny);
  for (int ix=0; ix<Nx; ix++) {
    double x = (double)ix * deltax;
    for (int iy=0; iy<Ny; iy++) {
      double y = (double)iy * deltay;
      for (int iz=0; iz<Nz; iz++) {
	double z = (double)iz * deltaz;
	exactOrb.Evaluate (x, y, z, val, grad, secDerivs);
	//exactOrb.Evaluate (x, y, z, val, grad, lapl);
	//	double curve = norm(secDerivs(0,0)/val);
	double curve = norm(secDerivs(0,0));
	xCurvature(ix) += prefactor * curve;
      }
    }
  }
  FILE *fout = fopen ("xcurvature.dat", "w");
  Array<double,1> xCurve2(Nx);
  xCurve2 = 0.0;
  // Average a few grid points
  for (int ix=0; ix<Nx; ix++) {
    for (int i=0; i<5; i++) {
      int j = (ix-2+i + Nx)%Nx;
      xCurve2(ix) += 0.2*xCurvature(j);
    }
  }
  for (int ix=0; ix<Nx; ix++) {
    double x = (double)ix * deltax;
    fprintf (fout, "%1.10e %1.10e\n", x, xCurve2(ix));
  }
  fclose (fout);

  cerr << "Read (" << orb.extent(0) << " x " << orb.extent(1) << " x " 
       << orb.extent(2) << ") size orbital.\n";
}


double 
UniformError (TricubicBspline<complex<double> > &exactSpline,
	      int n)
{
  TricubicBspline<complex<double> > smallSpline;
  double xi,xf,yi,yf,zi,zf;
  exactSpline.GetExtents(xi,xf,yi,yf,zi,zf);
  Array<complex<double>,3> smallVals (n,n,n);
  for (int ix=0; ix<n; ix++) {
    double x = xi + (xf-xi)*(double)ix/(double)(n-1);
    for (int iy=0; iy<n; iy++) {
      double y = yi + (yf-yi)*(double)iy/(double)(n-1);
      for (int iz=0; iz<n; iz++) {
	double z = zi + (zf-zi)*(double)iz/(double)(n-1);
	smallVals(ix,iy,iz) = exactSpline(x,y,z);
      }
    }
  }
  smallSpline.Init(xi,xf,yi,yf,zi,zf,smallVals, true,
		   NATURAL, NATURAL, NATURAL);

  double err = 0.0;
  int numSamples = 100000;
  double normSum = 0.0;
  for (int i=0; i<numSamples; i++) {
    double x = xi + (xf-xi)*drand48();
    double y = yi + (yf-yi)*drand48();
    double z = zi + (zf-zi)*drand48();
    complex<double> exVal = exactSpline (x,y,z);
    complex<double> smVal = smallSpline (x,y,z);
    err += norm (exVal - smVal);
    normSum += norm (exVal);
  }
  return (sqrt(err/normSum));
}


double 
UniformError2 (TricubicBspline<complex<double> > &exactSpline, int n)
{
  TricubicNUBspline<complex<double>,LinearGrid> smallSpline;
  double xi,xf,yi,yf,zi,zf;
  exactSpline.GetExtents(xi,xf,yi,yf,zi,zf);
  LinearGrid xGrid, yGrid, zGrid;
  xGrid.Init (xi, xf, n);
  yGrid.Init (yi, yf, n);
  zGrid.Init (zi, zf, n);

  Array<complex<double>,3> smallVals (n,n,n);
  for (int ix=0; ix<n; ix++) {
    double x = xGrid(ix)-1.0e-8;
    for (int iy=0; iy<n; iy++) {
      double y = yGrid(iy)-1.0e-8; 
      for (int iz=0; iz<n; iz++) {
	double z = zGrid(iz)-1.0e-8;
	smallVals(ix,iy,iz) = exactSpline(x,y,z);
      }
    }
  }
  smallSpline.Init(&xGrid, &yGrid, &zGrid, smallVals,
		   NATURAL, NATURAL, NATURAL);

  double err = 0.0;
  int numSamples = 100000;
  double normSum = 0.0;
  for (int i=0; i<numSamples; i++) {
    double x = xi + (xf-xi)*drand48();
    double y = yi + (yf-yi)*drand48();
    double z = zi + (zf-zi)*drand48();
    Vec3 r (x,y,z);
    complex<double> exVal = exactSpline (x,y,z);
    complex<double> smVal = smallSpline (r);
    err += norm (exVal - smVal);
    normSum += norm (exVal);
  }
  return (sqrt(err/normSum));
}


double 
NonuniformError (TricubicBspline<complex<double> > &exactSpline,
		 int n, double ratio)
{
  TricubicNUBspline<complex<double>,CenterGrid> smallSpline;
  double xi,xf,yi,yf,zi,zf;
  exactSpline.GetExtents(xi,xf,yi,yf,zi,zf);
  CenterGrid xGrid, yGrid, zGrid;
  xGrid.Init (xi, xf, ratio, n);
  yGrid.Init (yi, yf, ratio, n);
  zGrid.Init (zi, zf, ratio, n);

  Array<complex<double>,3> smallVals (n,n,n);
  for (int ix=0; ix<n; ix++) {
    double x = xGrid(ix)-1.0e-8;
    for (int iy=0; iy<n; iy++) {
      double y = yGrid(iy)-1.0e-8; 
      for (int iz=0; iz<n; iz++) {
	double z = zGrid(iz)-1.0e-8;
	smallVals(ix,iy,iz) = exactSpline(x,y,z);
      }
    }
  }
  smallSpline.Init(&xGrid, &yGrid, &zGrid, smallVals,
		   NATURAL, NATURAL, NATURAL);

  double err = 0.0;
  int numSamples = 100000;
  double normSum = 0.0;
  for (int i=0; i<numSamples; i++) {
    double x = xi + (xf-xi)*drand48();
    double y = yi + (yf-yi)*drand48();
    double z = zi + (zf-zi)*drand48();
    Vec3 r (x,y,z);
    complex<double> exVal = exactSpline (x,y,z);
    complex<double> smVal = smallSpline (r);
    err += norm (exVal - smVal);
    normSum += norm (exVal);
  }
  return (sqrt(err/normSum));
}

double
NonuniformError (TricubicBspline<complex<double> > &exactSpline, int n)
{
  double minError = NonuniformError(exactSpline, n, 1.001);
  double bestRatio = 1.001;
  for (double ratio = 1.1; ratio < 15.0; ratio += 0.05) {
    double err = NonuniformError (exactSpline, n, ratio);
    if (err < minError) {
      minError = err;
      bestRatio = ratio;
    }
  }
  cerr << "bestRatio = " << bestRatio << endl;
  return minError;
}




main()
{
  TricubicBspline<complex<double> > exactSpline;
  ReadData(exactSpline);

  for (int n = 10; n <= 80; n+= 2) {
    cerr << "n = " << n << endl;
    fprintf (stdout, "%d %15.10e %15.10e\n", n,
	     UniformError2 (exactSpline, n),
	     NonuniformError (exactSpline, n));
    fflush (stdout);
  }

//   cerr << "UniformError   (50) = " << UniformError (exactSpline, 50) << endl;
//   cerr << "NonuniformError(50) = " 
//        <<  NonuniformError (exactSpline, 50) << endl;

}
