#include "TricubicBspline.h"

void
TestValue()
{
  int nx=25, ny=31, nz=33;
  Array<double,3> data(nx,ny,nz);
  for (int ix=0; ix<nx; ix++) {
    double x = 2.0*M_PI*ix/(double)nx;
    for (int iy=0; iy<ny; iy++) {
      double y = 2.0*M_PI*iy/(double)ny;
      for (int iz=0; iz<nz; iz++) {
	double z = 2.0*M_PI*iz/(double)nz;
	data(ix,iy,iz) = sin(x)*cos(y)*sin(z);
      }
    }
  }
  double xi=0.0; double xf=2.0*M_PI;
  double yi=0.0; double yf=2.0*M_PI;
  double zi=0.0; double zf=2.0*M_PI;
  TricubicBspline<double> bspline;
  bspline.Init(xi, xf, yi, yf, zi, zf, data, true, 
	       PERIODIC, PERIODIC, PERIODIC);

  for (double t=0.0; t<=1.0; t+=0.001) {
    double y = 2.0*M_PI*t;
    double x = 1.1+0.3*M_PI*t;
    double z = 1.7+0.5*M_PI*t;
    double exact = sin(x)*cos(y)*sin(z);
    double sp    = bspline(x,y,z);
    fprintf (stdout, "%20.16e %20.16e %20.16e\n", t, sp, exact);
  }
}

void
TestComplexValue()
{
  int nx=25, ny=31, nz=33;
  Array<complex<double>,3> data(nx,ny,nz);
  for (int ix=0; ix<nx; ix++) {
    double x = 2.0*M_PI*ix/(double)nx;
    for (int iy=0; iy<ny; iy++) {
      double y = 2.0*M_PI*iy/(double)ny;
      for (int iz=0; iz<nz; iz++) {
	double z = 2.0*M_PI*iz/(double)nz;
	data(ix,iy,iz) = complex<double>(sin(x)*cos(y)*sin(z),
					 -cos(x)*sin(y)*cos(z));
      }
    }
  }
  double xi=0.0; double xf=2.0*M_PI;
  double yi=0.0; double yf=2.0*M_PI;
  double zi=0.0; double zf=2.0*M_PI;
  TricubicBspline<complex<double> > bspline;
  bspline.Init(xi, xf, yi, yf, zi, zf, data, true, 
	       PERIODIC, PERIODIC, PERIODIC);

  for (double t=0.0; t<=1.0; t+=0.001) {
    double y = 2.0*M_PI*t;
    double x = 1.1+0.3*M_PI*t;
    double z = 1.7+0.5*M_PI*t;
    complex<double> exact(sin(x)*cos(y)*sin(z), 
			  -cos(x)*sin(y)*cos(z));
    complex<double> sp    = bspline(x,y,z);
    fprintf (stdout, "%20.16e %20.16e %20.16e %20.16e %20.16e\n", 
	     t, sp.real(), sp.imag(), exact.real(), exact.imag());
  }
}

void
TestGrad()
{
  int nx=25, ny=31, nz=33;
  Array<double,3> data(nx,ny,nz);
  for (int ix=0; ix<nx; ix++) {
    double x = 2.0*M_PI*ix/(double)nx;
    for (int iy=0; iy<ny; iy++) {
      double y = 2.0*M_PI*iy/(double)ny;
      for (int iz=0; iz<nz; iz++) {
	double z = 2.0*M_PI*iz/(double)nz;
	data(ix,iy,iz) = sin(x)*cos(y)*sin(z);
      }
    }
  }
  double xi=0.0; double xf=2.0*M_PI;
  double yi=0.0; double yf=2.0*M_PI;
  double zi=0.0; double zf=2.0*M_PI;
  TricubicBspline<double> bspline;
  bspline.Init(xi, xf, yi, yf, zi, zf, data, true,
	       PERIODIC, PERIODIC, PERIODIC);

  for (double t=0.0; t<=1.0; t+=0.001) {
    double y = 2.0*M_PI*t;
    double x = 1.1+0.3*M_PI*t;
    double z = 1.7+0.5*M_PI*t;
    TinyVector<double,3> exact;
    exact[0] = cos(x)*cos(y)*sin(z);
    exact[1] =-sin(x)*sin(y)*sin(z);
    exact[2] = sin(x)*cos(y)*cos(z);
    TinyVector<double,3> sp    = bspline.Grad(x,y,z);
    fprintf (stdout, "%20.16e %20.16e %20.16e %20.16e %20.16e  %20.16e %20.16e\n", t, 
	     sp[0], sp[1], sp[2], exact[0], exact[1], exact[2]);
  }

}

class ExactFunc
{
private:
  Array<complex<double>,3> c;

public:
  double operator()(double x, double y, double z)
  {
    double val = 0.0;
    for (int i=0; i<c.extent(0); i++) {
      complex<double> e2inx(cos(x*i),sin(x*i));
      for (int j=0; j<c.extent(1); j++) {
	complex<double> e2iny(cos(y*j),sin(y*j));
	for (int k=0; k<c.extent(2); k++) {
	  complex<double> e2inz(cos(z*k),sin(z*k));
	  complex<double> e2inr = e2inx*e2iny*e2inz;
	  val += real(c(i,j,k)*e2inr);
	}
      }
    }
    return val;
  }

  void Evaluate (double x, double y, double z,
		 double &val, TinyVector<double,3> &grad,
		 double &laplacian)
  {
    val = (*this)(x,y,z);
    grad[0] = 0.0;  grad[1] = 0.0;  grad[2] = 0.0;
    laplacian = 0.0;
    complex<double> eye(0.0, 1.0);
    for (int i=0; i<c.extent(0); i++) {
      complex<double> e2inx(cos(x*i),sin(x*i));
      for (int j=0; j<c.extent(1); j++) {
	complex<double> e2iny(cos(y*j),sin(y*j));
	for (int k=0; k<c.extent(2); k++) {
	  complex<double> e2inz(cos(z*k),sin(z*k));
	  complex<double> e2inr = e2inx*e2iny*e2inz;
	  grad[0] += real(eye*(double)i * c(i,j,k)*e2inr);
	  grad[1] += real(eye*(double)j * c(i,j,k)*e2inr);
	  grad[2] += real(eye*(double)k * c(i,j,k)*e2inr);
	  laplacian +=
	    real(-(double)(i*i)*c(i,j,k)*e2inr)+
	    real(-(double)(j*j)*c(i,j,k)*e2inr)+
	    real(-(double)(k*k)*c(i,j,k)*e2inr);
	}
      }
    }
  }
  void Evaluate (double x, double y, double z,
		 double &val, TinyVector<double,3> &grad,
		 TinyMatrix<double,3,3> &secDerivs)
  {
    val = (*this)(x,y,z);
    grad[0] = 0.0;  grad[1] = 0.0;  grad[2] = 0.0;
    secDerivs = 0.0;
    complex<double> eye(0.0, 1.0);
    for (int i=0; i<c.extent(0); i++) {
      complex<double> e2inx(cos(x*i),sin(x*i));
      for (int j=0; j<c.extent(1); j++) {
	complex<double> e2iny(cos(y*j),sin(y*j));
	for (int k=0; k<c.extent(2); k++) {
	  complex<double> e2inz(cos(z*k),sin(z*k));
	  complex<double> e2inr = e2inx*e2iny*e2inz;
	  grad[0] += real(eye*(double)i * c(i,j,k)*e2inr);
	  grad[1] += real(eye*(double)j * c(i,j,k)*e2inr);
	  grad[2] += real(eye*(double)k * c(i,j,k)*e2inr);
	  secDerivs(0,0) +=  real(-(double)(i*i)*c(i,j,k)*e2inr);
	  secDerivs(1,1) +=  real(-(double)(j*j)*c(i,j,k)*e2inr);
	  secDerivs(2,2) +=  real(-(double)(k*k)*c(i,j,k)*e2inr);
	  secDerivs(0,1) +=  real(-(double)(i*j)*c(i,j,k)*e2inr);
	  secDerivs(0,2) +=  real(-(double)(i*k)*c(i,j,k)*e2inr);
	  secDerivs(1,2) +=  real(-(double)(j*k)*c(i,j,k)*e2inr);
	}
      }
    }
    secDerivs(1,0) = secDerivs(0,1);
    secDerivs(2,0) = secDerivs(0,2);
    secDerivs(2,1) = secDerivs(1,2);
  }

  ExactFunc(int nx, int ny, int nz)
  {
    c.resize(nx, ny, nz);
    for (int i=0; i<nx; i++)
      for (int j=0; j<ny; j++)
	for (int k=0; k<nz; k++)
	  c(i,j,k) = complex<double>(2.0*drand48()-1.0, 2.0*drand48()-1);
  }
};


void
TestAll()
{
  int Nplx=10, Nply=10, Nplz=10;
  ExactFunc ex(Nplx, Nply, Nplz);

  int Nspx=60, Nspy=60, Nspz=60;  
  Array<double,3> data(Nspx,Nspy,Nspz);
  for (int ix=0; ix<Nspx; ix++) {
    double x = 2.0*M_PI*ix/(double)Nspx;
    for (int iy=0; iy<Nspy; iy++) {
      double y = 2.0*M_PI*iy/(double)Nspy;
      for (int iz=0; iz<Nspz; iz++) {
	double z = 2.0*M_PI*iz/(double)Nspz;
	data (ix, iy, iz) = ex(x,y,z);
      }
    }
  }
  double xi=0.0; double xf=2.0*M_PI;
  double yi=0.0; double yf=2.0*M_PI;
  double zi=0.0; double zf=2.0*M_PI;
  TricubicBspline<double> interp, noInterp;
  interp.Init(xi, xf, yi, yf, zi, zf, data, true, 
	      PERIODIC, PERIODIC, PERIODIC);
  noInterp.Init(xi, xf, yi, yf, zi, zf, data, false, 
		PERIODIC, PERIODIC, PERIODIC);
  
  FILE *exOut = fopen ("exact.dat", "w");
  FILE *inOut = fopen ("interp.dat", "w");
  FILE *noinOut = fopen ("nointerp.dat", "w");

  for (double t=0.0; t<=1.0; t+=0.001) {
    double y = 2.0*M_PI*t;
    double x = 1.1+0.3*M_PI*t;
    double z = 1.7+0.5*M_PI*t;
    TinyVector<double,3> gradEx, gradInterp, gradNoInterp;
    double valEx, valInterp, valNoInterp, lapEx, lapInterp, lapNoInterp;
    ex.Evaluate (x, y, z, valEx, gradEx, lapEx);
    interp.Evaluate(x, y, z, valInterp, gradInterp, lapInterp);
    noInterp.Evaluate(x, y, z, valNoInterp, gradNoInterp, lapNoInterp);
    fprintf (exOut, "%20.16e %20.16e %20.16e %20.16e %20.16e  %20.16e\n",
	     t, valEx, gradEx[0], gradEx[1], gradEx[2], lapEx);
    fprintf (inOut, "%20.16e %20.16e %20.16e %20.16e %20.16e  %20.16e\n",
	     t, valInterp, gradInterp[0], gradInterp[1], gradInterp[2], lapInterp);
    fprintf (noinOut, "%20.16e %20.16e %20.16e %20.16e %20.16e  %20.16e\n",
	     t, valNoInterp, gradNoInterp[0], gradNoInterp[1], 
	     gradNoInterp[2], lapNoInterp);
  }
  fclose (exOut);
  fclose (inOut);
  fclose (noinOut);
}

void
TestAll2()
{
  int Nplx=10, Nply=10, Nplz=10;
  ExactFunc ex(Nplx, Nply, Nplz);

  int Nspx=60, Nspy=60, Nspz=60;  
  Array<double,3> data(Nspx,Nspy,Nspz);
  for (int ix=0; ix<Nspx; ix++) {
    double x = 2.0*M_PI*ix/(double)Nspx;
    for (int iy=0; iy<Nspy; iy++) {
      double y = 2.0*M_PI*iy/(double)Nspy;
      for (int iz=0; iz<Nspz; iz++) {
	double z = 2.0*M_PI*iz/(double)Nspz;
	data (ix, iy, iz) = ex(x,y,z);
      }
    }
  }
  double xi=0.0; double xf=2.0*M_PI;
  double yi=0.0; double yf=2.0*M_PI;
  double zi=0.0; double zf=2.0*M_PI;
  TricubicBspline<double> interp, noInterp;
  interp.Init(xi, xf, yi, yf, zi, zf, data, true, 
	      PERIODIC, PERIODIC, PERIODIC);
  noInterp.Init(xi, xf, yi, yf, zi, zf, data, false, 
		PERIODIC, PERIODIC, PERIODIC);
  
  FILE *exOut = fopen ("exact.dat", "w");
  FILE *inOut = fopen ("interp.dat", "w");
  FILE *noinOut = fopen ("nointerp.dat", "w");

  for (double t=0.0; t<=1.0; t+=0.001) {
    double y = 2.0*M_PI*t;
    double x = 1.1+0.3*M_PI*t;
    double z = 1.7+0.5*M_PI*t;
    TinyVector<double,3> gradEx, gradInterp, gradNoInterp;
    TinyMatrix<double,3,3> d2Ex, d2Interp, d2NoInterp;
    double valEx, valInterp, valNoInterp, lapEx, lapInterp, lapNoInterp;
    ex.Evaluate (x, y, z, valEx, gradEx, d2Ex);
    interp.Evaluate(x, y, z, valInterp, gradInterp, d2Interp);
    noInterp.Evaluate(x, y, z, valNoInterp, gradNoInterp, d2NoInterp);
    fprintf (exOut, "%20.16e %20.16e %20.16e %20.16e %20.16e %20.16e %20.16e %20.16e %20.16e %20.16e %20.16e %20.16e %20.16e %20.16e\n",
	     t, valEx, gradEx[0], gradEx[1], gradEx[2], 
	     d2Ex(0,0), d2Ex(0,1), d2Ex(0,2), 
	     d2Ex(1,0), d2Ex(1,1), d2Ex(1,2), 
	     d2Ex(2,0), d2Ex(2,1), d2Ex(2,2));
    fprintf (inOut, "%20.16e %20.16e %20.16e %20.16e %20.16e %20.16e %20.16e %20.16e %20.16e %20.16e %20.16e %20.16e %20.16e %20.16e\n",
	     t, valInterp, gradInterp[0], gradInterp[1], gradInterp[2], 
	     d2Interp(0,0), d2Interp(0,1), d2Interp(0,2), 
	     d2Interp(1,0), d2Interp(1,1), d2Interp(1,2), 
	     d2Interp(2,0), d2Interp(2,1), d2Interp(2,2));
	     
    fprintf (noinOut, "%20.16e %20.16e %20.16e %20.16e %20.16e %20.16e %20.16e %20.16e %20.16e %20.16e %20.16e %20.16e %20.16e %20.16e\n",
	     t, valNoInterp, gradNoInterp[0], gradNoInterp[1], gradNoInterp[2], 
	     d2NoInterp(0,0), d2NoInterp(0,1), d2NoInterp(0,2), 
	     d2NoInterp(1,0), d2NoInterp(1,1), d2NoInterp(1,2), 
	     d2NoInterp(2,0), d2NoInterp(2,1), d2NoInterp(2,2));

  }
  fclose (exOut);
  fclose (inOut);
  fclose (noinOut);
}


void
TestFlat()
{
  int Nplx=10, Nply=11, Nplz=12;
  ExactFunc ex(Nplx, Nply, Nplz);

  int Nspx=60, Nspy=61, Nspz=62;  
  Array<double,3> data(Nspx,Nspy,Nspz);
  for (int ix=0; ix<Nspx; ix++) {
    double x = 2.0*M_PI*ix/(double)Nspx;
    for (int iy=0; iy<Nspy; iy++) {
      double y = 2.0*M_PI*iy/(double)Nspy;
      for (int iz=0; iz<Nspz; iz++) {
	double z = 2.0*M_PI*iz/(double)Nspz;
	data (ix, iy, iz) = ex(x,y,z);
      }
    }
  }
  double xi=0.0; double xf=2.0*M_PI;
  double yi=0.0; double yf=2.0*M_PI;
  double zi=0.0; double zf=2.0*M_PI;
  TricubicBspline<double> spline, noInterp;
  spline.Init(xi, xf, yi, yf, zi, zf, data, true, 
	      NATURAL, NATURAL, NATURAL);
  
  FILE *fout = fopen ("flat.dat", "w");
  for (double t=0.0; t<=1.000001; t+=0.001) {
    double x = t*xf;
    double y = 0.2+0.1*t*yf;
    double z = 0.8+0.5*t*zf;
    TinyVector<double,3> r(x,y,z);
    TinyVector<double,3> grad;
    double val, lapl;

    spline.Evaluate(r, val, grad, lapl);
    fprintf (fout, "%20.16e %20.16e %20.16e %20.26e %20.26e %20.16e\n", t, val,
	     grad[0], grad[1], grad[2], lapl);

  }

  fclose(fout);

}


void
TestEvaluate()
{
  int nx=8, ny=8, nz=8;
  Array<double,3> data(nx,ny,nz);
  for (int ix=0; ix<nx; ix++) {
    double x = 2.0*M_PI*ix/(double)nx;
    for (int iy=0; iy<ny; iy++) {
      double y = 2.0*M_PI*iy/(double)ny;
      for (int iz=0; iz<nz; iz++) {
	double z = 2.0*M_PI*iz/(double)nz;
	data(ix,iy,iz) = sin(2.0*x)*cos(2.0*y)*sin(2.0*z)
	+1.31312*sin(1.0*x)*cos(1.0*y)*sin(1.0*z);
      }
    }
  }
  double xi=0.0; double xf=2.0*M_PI;
  double yi=0.0; double yf=2.0*M_PI;
  double zi=0.0; double zf=2.0*M_PI;
  TricubicBspline<double> interpSpline;
  TricubicBspline<double> nointerp;
  interpSpline.Init(xi, xf, yi, yf, zi, zf, data, true, 
		    PERIODIC, PERIODIC, PERIODIC);
  nointerp.Init(xi, xf, yi, yf, zi, zf, data, false, 
		PERIODIC, PERIODIC, PERIODIC);
  
  FILE *exOut = fopen ("exact.dat", "w");
  FILE *inOut = fopen ("interp.dat", "w");
  FILE *noinOut = fopen ("nointerp.dat", "w");

  for (double t=0.0; t<=1.0; t+=0.001) {
    double y = 2.0*M_PI*t;
    double x = 1.1+0.3*M_PI*t;
    double z = 1.7+0.5*M_PI*t;
    TinyVector<double,3> gradEx, gradSp;
    double valEx, valSp, lapEx, lapSp;
    valEx = sin(2.0*x)*cos(2.0*y)*sin(2.0*z);
    gradEx[0] = 2.0*cos(2.0*x)*cos(2.0*y)*sin(2.0*z);
    gradEx[1] = 2.0*-sin(2.0*x)*sin(2.0*y)*sin(2.0*z);
    gradEx[2] = 2.0*sin(2.0*x)*cos(2.0*y)*cos(2.0*z);
    lapEx = -12.0*valEx;
    fprintf (exOut, "%20.16e %20.16e %20.16e %20.16e %20.16e  %20.16e\n", 
	     t, valEx, gradEx[0], gradEx[1], gradEx[2], lapEx);
    interpSpline.Evaluate(x,y,z, valSp, gradSp, lapSp);
    fprintf (inOut, "%20.16e %20.16e %20.16e %20.16e %20.16e  %20.16e\n", 
	     t, valSp, gradSp[0], gradSp[1], gradSp[2], lapSp);
    nointerp.Evaluate(x,y,z, valSp, gradSp, lapSp);
    fprintf (noinOut, "%20.16e %20.16e %20.16e %20.16e %20.16e  %20.16e\n", 
	     t, valSp, gradSp[0], gradSp[1], gradSp[2], lapSp);
  }
  fclose (exOut);
  fclose (inOut);
  fclose (noinOut);

}





main()
{
  // TestAll2();
  TestFlat();
}
