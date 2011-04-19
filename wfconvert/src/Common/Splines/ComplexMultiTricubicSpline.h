/////////////////////////////////////////////////////////////
// Copyright (C) 2003-2006 Bryan Clark and Kenneth Esler   //
//                                                         //
// This program is free software; you can redistribute it  //
// and/or modify it under the terms of the GNU General     //
// Public License as published by the Free Software        //
// Foundation; either version 2 of the License, or         //
// (at your option) any later version.  This program is    //
// distributed in the hope that it will be useful, but     //
// WITHOUT ANY WARRANTY; without even the implied warranty //
// of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. //  
// See the GNU General Public License for more details.    //
// For more information, please see the PIMC++ Home Page:  //
//           http://pathintegrals.info                     //
/////////////////////////////////////////////////////////////

#ifndef COMPLEX_MULTI_TRICUBIC_SPLINE_H
#define COMPLEX_MULTI_TRICUBIC_SPLINE_H

#include "Common/config.h"

#ifndef NOUNDERSCORE
#define FORT(name) name ## _
#else
#define FORT(name) name
#endif

#define F77_Z3SPLINE F77_FUNC(z3spline,Z3SPLINE)
#define F77_Z3VALGRAD F77_FUNC(z3valgrad,Z3VALGRAD)

extern "C" void 
F77_Z3SPLINE (double *x, double *y, double *z,
	      double *x0, double *dx, int *nx,
	      double *y0, double *dy, int *ny,
	      double *z0, double *dz, int *nz,
	      void *F, int *num, void *vals);

extern "C" void 
F77_Z3VALGRAD (double *x, double *y, double *z,
	       double *x0, double *dx, int *nx,
	       double *y0, double *dy, int *ny,
	       double *z0, double *dz, int *nz,
	       void *F, int *num, void *vals, void *grads);

#include "Grid.h"
#include <cmath>
//#include <blitz/array.h>
//using namespace blitz;

/// Each point of F contains:
/// 0)  real [F(x,y,z)]
/// 1)  real [dF/dx]
/// 2)  real [dF/dy]
/// 3)  real [dF/dz]
/// 4)  real [d2F/dxdy]
/// 5)  real [d2F/dxdz]
/// 6)  real [d2F/dydz]
/// 7)  real [d3F/dxdydz]
/// 8)  imag [F(x,y,z)]
/// 9)  imag [dF/dx]
/// 10) imag [dF/dy]
/// 11) imag [dF/dz]
/// 12) imag [d2F/dxdy]
/// 13) imag [d2F/dxdz]
/// 14) imag [d2F/dydz]
/// 15) imag [d3F/dxdydz]


class ComplexMultiTricubicSpline
{
  inline double p1(double t)
  { return ((t-1.0)*(t-1.0)*(1.0+2.0*t)); }
  inline double p2(double t)
  { return (t*t*(3.0-2.0*t)); }
  inline double q1(double t)
  { return (t*(t-1.0)*(t-1.0)); }
  inline double q2(double t)
  { return (t*t*(t-1.0)); }
  inline double dp1(double t)
  { return (6.0*t*(t-1.0)); }
  inline double dq1(double t)
  { return ((t-1.0)*(3.0*t-1.0)); }
  inline double dp2(double t)
  { return (-dp1(t)); }
  inline double dq2 (double t)
  { return (3.0*t*t - 2.0*t); }
  inline double d2p1(double t)
  { return (12.0*t-6.0); }
  inline double d2q1 (double t)
  { return (6.0*t - 4.0); }
  inline double d2p2 (double t)
  { return (-d2p1(t)); }
  inline double d2q2 (double t)
  { return (6.0*t - 2.0); } 
  
  // The definite integrals of the hermite functions from 0 to 1
  static const double Int_p1, Int_p2, Int_q1, Int_q2;
  // The definite integrals of the second derivatives of the hermite
  // functions from 0 to 1 
  static const double Int_d2p1, Int_d2p2, Int_d2q1, Int_d2q2;

  // dim:     Dimension to calculate derivative w.r.t
  // source:  Function to differentiate
  // dest:    where to put result
  void UpdateX (int source, int dest, int i);
  void UpdateY (int source, int dest, int i);
  void UpdateZ (int source, int dest, int i);
  void UpdateXPeriodic (int source, int dest, int i);
  void UpdateYPeriodic (int source, int dest, int i);
  void UpdateZPeriodic (int source, int dest, int i);
  bool UpToDate, Periodic;
  double dx, dy, dz, dxInv, dyInv, dzInv;
public:
  Array<TinyVector<TinyVector<double,2>,8>,4> F;

  int Nx, Ny, Nz, N;
  //  double dx, dy,dz, dxInv, dyInv, dzInv;
  Grid *Xgrid, *Ygrid, *Zgrid;
  TinyVector<Grid*,3> Grids;
  void Update();
  inline complex<double> operator()(int ix, int iy, int iz, int i) const
  { return complex<double> (F(ix,iy,iz, i)[0][0], F(ix,iy,iz,i)[0][1]); }
  inline void Set (int ix, int iy, int iz, int i, complex<double> val)
    { UpToDate=false; F(ix, iy, iz, i)[0][0] = val.real(); F(ix,iy,iz,i)[0][1] = val.imag(); }
  inline void operator()(double x, double y, double z, 
			 Array<complex<double>,1> &vals);
  inline void d_dx      (double x, double y, double z, 
			 Array<complex<double>,1> &vals); 
  inline void d_dy      (double x, double y, double z, 
			 Array<complex<double>,1> &vals); 
  inline void d_dz      (double x, double y, double z, 
			 Array<complex<double>,1> &vals); 
  inline void d2_dx2    (double x, double y, double z, 
			 Array<complex<double>,1> &vals); 
  inline void d2_dy2    (double x, double y, double z, 
			 Array<complex<double>,1> &vals); 
  inline void d2_dz2    (double x, double y, double z, 
			 Array<complex<double>,1> &vals); 
  inline void d2_dxdy   (double x, double y, double z, 
			 Array<complex<double>,1> &vals); 
  inline void d2_dxdz   (double x, double y, double z, 
			 Array<complex<double>,1> &vals); 
  inline void d2_dydz   (double x, double y, double z, 
			 Array<complex<double>,1> &vals); 
  inline void Grad      (double x, double y, double z, 
			 Array<cVec3,1> &grads);
  inline void ValGrad   (double x, double y, double z, 
			 Array<complex<double>,1> &vals, 
			 Array<cVec3,  1> &grads);
  inline void FValGrad   (double x, double y, double z, 
			  Array<complex<double>,1> &vals, 
			  Array<cVec3,  1> &grads);
  inline void Laplacian (double x, double y, double z, 
			 Array<complex<double>,1> &vals);
  /// Returns the integral of the laplacian over the whole box.
  inline void IntegralOfLaplacian (Array<complex<double>,1> &vals);
  inline void KineticEnergy (Array<double,1> &ke);  
  inline void Norm (Array<double,1> &norm);
  inline void* Data () { return F.data(); }

  ComplexMultiTricubicSpline(Grid *xgrid, Grid *ygrid, Grid *zgrid, int n)
  {
    Xgrid = xgrid; Nx = xgrid->NumPoints;
    Ygrid = ygrid; Ny = ygrid->NumPoints;
    Zgrid = zgrid; Nz = zgrid->NumPoints;
    N = n;
    
    F.resize(Nx,Ny,Nz,N);
    UpToDate = false;
  }
  
  /// Copy constructor
  inline ComplexMultiTricubicSpline (const ComplexMultiTricubicSpline &a);

  /// Assigment operator -- necessary for array resizeAndPreserve
  inline ComplexMultiTricubicSpline& operator=(const ComplexMultiTricubicSpline &a);
  //inline ComplexMultiTricubicSpline & operator= (ComplexMultiTricubicSpline a);


  inline void Init (Grid *xgrid, Grid *ygrid, Grid *zgrid,
		    const Array<complex<double>,4> &init, bool periodic=false);

  ComplexMultiTricubicSpline(Grid *xgrid, Grid *ygrid, Grid *zgrid,
		      const Array<complex<double>,4> &init, bool periodic=false)
  {
    Init (xgrid, ygrid, zgrid, init, periodic);
  }
  ComplexMultiTricubicSpline() : UpToDate(false) 
  { /* Do nothing. */ }
};

inline 
ComplexMultiTricubicSpline::ComplexMultiTricubicSpline
(const ComplexMultiTricubicSpline &a)
{
  F.resize(a.F.shape());
  F = a.F;
  Nx = a.Nx;
  Ny = a.Ny;
  Nz = a.Nz;
  N  = a.N;
  Xgrid = a.Xgrid;
  Ygrid = a.Ygrid;
  Zgrid = a.Zgrid;
  UpToDate = a.UpToDate;
}


inline 
ComplexMultiTricubicSpline& ComplexMultiTricubicSpline::operator=
(const ComplexMultiTricubicSpline &a)
{
  F.resize(a.F.shape());
  F = a.F;
  Nx = a.Nx;
  Ny = a.Ny;
  Nz = a.Nz;
  N  = a.N;
  Xgrid = a.Xgrid;
  Ygrid = a.Ygrid;
  Zgrid = a.Zgrid;
  UpToDate = a.UpToDate;
  return (*this);
}

// inline ComplexMultiTricubicSpline& 
// ComplexMultiTricubicSpline::operator=(ComplexMultiTricubicSpline a)
// {
//   F.resize(a.F.shape());
//   F = a.F;
//   Nx = a.Nx;
//   Ny = a.Ny;
//   Nz = a.Nz;
//   N  = a.N;
//   Xgrid = a.Xgrid;
//   Ygrid = a.Ygrid;
//   Zgrid = a.Zgrid;
//   UpToDate = a.UpToDate;
//   return (*this);
// }



inline void 
ComplexMultiTricubicSpline::Init (Grid *xgrid, Grid *ygrid, Grid *zgrid,
				  const Array<complex<double>,4> &init,
				  bool periodic)
{
  Periodic = periodic;

  Xgrid = xgrid; Nx = xgrid->NumPoints;
  Ygrid = ygrid; Ny = ygrid->NumPoints;
  Zgrid = zgrid; Nz = zgrid->NumPoints;
  dx = (*xgrid)(1) - (*xgrid)(0);    dxInv = 1.0/dx;
  dy = (*ygrid)(1) - (*ygrid)(0);    dyInv = 1.0/dy;
  dz = (*zgrid)(1) - (*zgrid)(0);    dzInv = 1.0/dz;

  assert (init.extent(0) == Nx);
  assert (init.extent(1) == Ny);
  assert (init.extent(2) == Nz);
  N = init.extent(3);

  F.resize(Nx,Ny,Nz,N);
  for (int ix=0; ix<Nx; ix++)
    for (int iy=0; iy<Ny; iy++)
      for (int iz=0; iz<Nz; iz++)
	for (int i=0; i<N; i++) {
	  F(ix,iy,iz,i)[0][0] = init(ix,iy,iz,i).real();
	  F(ix,iy,iz,i)[0][1] = init(ix,iy,iz,i).imag();
	}
  UpToDate = false;
  Update();
}




inline void ComplexMultiTricubicSpline::operator()
  (double x, double y, double z, Array<complex<double>,1> &vals)
{
  if (!UpToDate)
    Update();

  int ix = Xgrid->ReverseMap(x);  
  int iy = Ygrid->ReverseMap(y);
  int iz = Zgrid->ReverseMap(z);


  ix = max(0,ix); ix = min(ix, Nx-2);
  iy = max(0,iy); iy = min(iy, Ny-2);
  iz = max(0,iz); iz = min(iz, Nz-2);

  double h = (*Xgrid)(ix+1) - (*Xgrid)(ix);
  double k = (*Ygrid)(iy+1) - (*Ygrid)(iy);
  double l = (*Zgrid)(iz+1) - (*Zgrid)(iz);
  double u = (x - (*Xgrid)(ix))/h;
  double v = (y - (*Ygrid)(iy))/k;
  double w = (z - (*Zgrid)(iz))/l;
  double a0 = p1(u);
  double a1 = p2(u);
  double a2 = h*q1(u);
  double a3 = h*q2(u);

  register double b0 = p1(v);
  register double b1 = p2(v);
  register double b2 = k*q1(v);
  register double b3 = k*q2(v);

  register double c0 = p1(w);
  register double c1 = p2(w);
  register double c2 = l*q1(w);
  register double c3 = l*q2(w);

  for (int i=0; i<N; i++) {
    //////////////////
    /// Real parts ///
    //////////////////
    TinyVector<double,2>& Y000 = F(ix,iy,iz,i)[0];      //   F
    TinyVector<double,2>& Y001 = F(ix,iy,iz+1,i)[0];    //   F
    TinyVector<double,2>& Y002 = F(ix,iy,iz,i)[3];      //  dF/dz
    TinyVector<double,2>& Y003 = F(ix,iy,iz+1,i)[3];    //  dF/dz
    TinyVector<double,2>& Y010 = F(ix,iy+1,iz,i)[0];    //   F
    TinyVector<double,2>& Y011 = F(ix,iy+1,iz+1,i)[0];  //   F
    TinyVector<double,2>& Y012 = F(ix,iy+1,iz,i)[3];    //  dF/dz
    TinyVector<double,2>& Y013 = F(ix,iy+1,iz+1,i)[3];  //  dF/dz
    TinyVector<double,2>& Y020 = F(ix,iy,iz,i)[2];      //  dF/dy
    TinyVector<double,2>& Y021 = F(ix,iy,iz+1,i)[2];    //  dF/dy
    TinyVector<double,2>& Y022 = F(ix,iy,iz,i)[6];      // d2F/dydz
    TinyVector<double,2>& Y023 = F(ix,iy,iz+1,i)[6];    // d2F/dydz
    TinyVector<double,2>& Y030 = F(ix,iy+1,iz,i)[2];    //  dF/dy
    TinyVector<double,2>& Y031 = F(ix,iy+1,iz+1,i)[2];  //  dF/dy
    TinyVector<double,2>& Y032 = F(ix,iy+1,iz,i)[6];    // d2F/dydz
    TinyVector<double,2>& Y033 = F(ix,iy+1,iz+1,i)[6];  // d2F/dydz
    
    TinyVector<double,2>& Y100 = F(ix+1,iy,iz,i)[0];      //   F
    TinyVector<double,2>& Y101 = F(ix+1,iy,iz+1,i)[0];    //   F
    TinyVector<double,2>& Y102 = F(ix+1,iy,iz,i)[3];      //  dF/dz
    TinyVector<double,2>& Y103 = F(ix+1,iy,iz+1,i)[3];    //  dF/dz
    TinyVector<double,2>& Y110 = F(ix+1,iy+1,iz,i)[0];    //   F
    TinyVector<double,2>& Y111 = F(ix+1,iy+1,iz+1,i)[0];  //   F
    TinyVector<double,2>& Y112 = F(ix+1,iy+1,iz,i)[3];    //  dF/dz
    TinyVector<double,2>& Y113 = F(ix+1,iy+1,iz+1,i)[3];  //  dF/dz
    TinyVector<double,2>& Y120 = F(ix+1,iy,iz,i)[2];      //  dF/dy
    TinyVector<double,2>& Y121 = F(ix+1,iy,iz+1,i)[2];    //  dF/dy
    TinyVector<double,2>& Y122 = F(ix+1,iy,iz,i)[6];      // d2F/dydz
    TinyVector<double,2>& Y123 = F(ix+1,iy,iz+1,i)[6];    // d2F/dydz
    TinyVector<double,2>& Y130 = F(ix+1,iy+1,iz,i)[2];    //  dF/dy
    TinyVector<double,2>& Y131 = F(ix+1,iy+1,iz+1,i)[2];  //  dF/dy
    TinyVector<double,2>& Y132 = F(ix+1,iy+1,iz,i)[6];    // d2F/dydz
    TinyVector<double,2>& Y133 = F(ix+1,iy+1,iz+1,i)[6];  // d2F/dydz
    
    TinyVector<double,2>& Y200 = F(ix,iy,iz,i)[1];      //  dF/dx
    TinyVector<double,2>& Y201 = F(ix,iy,iz+1,i)[1];    //  dF/dx
    TinyVector<double,2>& Y202 = F(ix,iy,iz,i)[5];      // d2F/dxdz
    TinyVector<double,2>& Y203 = F(ix,iy,iz+1,i)[5];    // d2F/dxdz
    TinyVector<double,2>& Y210 = F(ix,iy+1,iz,i)[1];    //  dF/dx
    TinyVector<double,2>& Y211 = F(ix,iy+1,iz+1,i)[1];  //  dF/dx
    TinyVector<double,2>& Y212 = F(ix,iy+1,iz,i)[5];    // d2F/dxdz
    TinyVector<double,2>& Y213 = F(ix,iy+1,iz+1,i)[5];  // d2F/dxdz
    TinyVector<double,2>& Y220 = F(ix,iy,iz,i)[4];      // d2F/dxdy
    TinyVector<double,2>& Y221 = F(ix,iy,iz+1,i)[4];    // d2F/dxdy
    TinyVector<double,2>& Y222 = F(ix,iy,iz,i)[7];      // d3F/dxdydz
    TinyVector<double,2>& Y223 = F(ix,iy,iz+1,i)[7];    // d3F/dxdydz
    TinyVector<double,2>& Y230 = F(ix,iy+1,iz,i)[4];    // d2F/dxdy
    TinyVector<double,2>& Y231 = F(ix,iy+1,iz+1,i)[4];  // d2F/dxdy
    TinyVector<double,2>& Y232 = F(ix,iy+1,iz,i)[7];    // d3F/dxdydz
    TinyVector<double,2>& Y233 = F(ix,iy+1,iz+1,i)[7];  // d3F/dxdydz
    
    TinyVector<double,2>& Y300 = F(ix+1,iy,iz,i)[1];      //  dF/dx
    TinyVector<double,2>& Y301 = F(ix+1,iy,iz+1,i)[1];    //  dF/dx
    TinyVector<double,2>& Y302 = F(ix+1,iy,iz,i)[5];      // d2F/dxdz
    TinyVector<double,2>& Y303 = F(ix+1,iy,iz+1,i)[5];    // d2F/dxdz
    TinyVector<double,2>& Y310 = F(ix+1,iy+1,iz,i)[1];    //  dF/dx
    TinyVector<double,2>& Y311 = F(ix+1,iy+1,iz+1,i)[1];  //  dF/dx
    TinyVector<double,2>& Y312 = F(ix+1,iy+1,iz,i)[5];    // d2F/dxdz
    TinyVector<double,2>& Y313 = F(ix+1,iy+1,iz+1,i)[5];  // d2F/dxdz
    TinyVector<double,2>& Y320 = F(ix+1,iy,iz,i)[4];      // d2F/dxdy
    TinyVector<double,2>& Y321 = F(ix+1,iy,iz+1,i)[4];    // d2F/dxdy
    TinyVector<double,2>& Y322 = F(ix+1,iy,iz,i)[7];      // d3F/dxdydz
    TinyVector<double,2>& Y323 = F(ix+1,iy,iz+1,i)[7];    // d3F/dxdydz
    TinyVector<double,2>& Y330 = F(ix+1,iy+1,iz,i)[4];    // d2F/dxdy
    TinyVector<double,2>& Y331 = F(ix+1,iy+1,iz+1,i)[4];  // d2F/dxdy
    TinyVector<double,2>& Y332 = F(ix+1,iy+1,iz,i)[7];    // d3F/dxdydz
    TinyVector<double,2>& Y333 = F(ix+1,iy+1,iz+1,i)[7];  // d3F/dxdydz
    
    TinyVector<double,2> Val;
    Val = 
      a0*
      (b0*(Y000*c0+Y001*c1+Y002*c2+Y003*c3) +
       b1*(Y010*c0+Y011*c1+Y012*c2+Y013*c3) +
       b2*(Y020*c0+Y021*c1+Y022*c2+Y023*c3) +
       b3*(Y030*c0+Y031*c1+Y032*c2+Y033*c3))+
      a1 *
      (b0*(Y100*c0+Y101*c1+Y102*c2+Y103*c3) +
       b1*(Y110*c0+Y111*c1+Y112*c2+Y113*c3) +
       b2*(Y120*c0+Y121*c1+Y122*c2+Y123*c3) +
       b3*(Y130*c0+Y131*c1+Y132*c2+Y133*c3))+
      a2 *
      (b0*(Y200*c0+Y201*c1+Y202*c2+Y203*c3) +
       b1*(Y210*c0+Y211*c1+Y212*c2+Y213*c3) +
       b2*(Y220*c0+Y221*c1+Y222*c2+Y223*c3) +
       b3*(Y230*c0+Y231*c1+Y232*c2+Y233*c3))+
      a3 *
      (b0*(Y300*c0+Y301*c1+Y302*c2+Y303*c3) +
       b1*(Y310*c0+Y311*c1+Y312*c2+Y313*c3) +
       b2*(Y320*c0+Y321*c1+Y322*c2+Y323*c3) +
       b3*(Y330*c0+Y331*c1+Y332*c2+Y333*c3));
    vals(i) = complex<double>(Val[0],Val[1]);
  }
}


inline void 
ComplexMultiTricubicSpline::d_dx (double x, double y, double z, 
				  Array<complex<double>,1> &vals)
{
  if (!UpToDate)
    Update();

  int ix = Xgrid->ReverseMap(x);  
  int iy = Ygrid->ReverseMap(y);
  int iz = Zgrid->ReverseMap(z);

  ix = max(0,ix); ix = min(ix, Nx-2);
  iy = max(0,iy); iy = min(iy, Ny-2);
  iz = max(0,iz); iz = min(iz, Nz-2);

  double h = (*Xgrid)(ix+1) - (*Xgrid)(ix);
  double hinv = 1.0/h;
  double k = (*Ygrid)(iy+1) - (*Ygrid)(iy);
  double l = (*Zgrid)(iz+1) - (*Zgrid)(iz);
  double u = (x - (*Xgrid)(ix))/h;
  double v = (y - (*Ygrid)(iy))/k;
  double w = (z - (*Zgrid)(iz))/l;

  double da0 = hinv*dp1(u);
  double da1 = hinv*dp2(u);
  double da2 = dq1(u);
  double da3 = dq2(u);

  register double b0 = p1(v);
  register double b1 = p2(v);
  register double b2 = k*q1(v);
  register double b3 = k*q2(v);

  register double c0 = p1(w);
  register double c1 = p2(w);
  register double c2 = l*q1(w);
  register double c3 = l*q2(w);
  
  for (int i=0; i<N; i++) {
    TinyVector<double,2>& Y000 = F(ix,iy,iz,i)[0];      //   F
    TinyVector<double,2>& Y001 = F(ix,iy,iz+1,i)[0];    //   F
    TinyVector<double,2>& Y002 = F(ix,iy,iz,i)[3];      //  dF/dz
    TinyVector<double,2>& Y003 = F(ix,iy,iz+1,i)[3];    //  dF/dz
    TinyVector<double,2>& Y010 = F(ix,iy+1,iz,i)[0];    //   F
    TinyVector<double,2>& Y011 = F(ix,iy+1,iz+1,i)[0];  //   F
    TinyVector<double,2>& Y012 = F(ix,iy+1,iz,i)[3];    //  dF/dz
    TinyVector<double,2>& Y013 = F(ix,iy+1,iz+1,i)[3];  //  dF/dz
    TinyVector<double,2>& Y020 = F(ix,iy,iz,i)[2];      //  dF/dy
    TinyVector<double,2>& Y021 = F(ix,iy,iz+1,i)[2];    //  dF/dy
    TinyVector<double,2>& Y022 = F(ix,iy,iz,i)[6];      // d2F/dydz
    TinyVector<double,2>& Y023 = F(ix,iy,iz+1,i)[6];    // d2F/dydz
    TinyVector<double,2>& Y030 = F(ix,iy+1,iz,i)[2];    //  dF/dy
    TinyVector<double,2>& Y031 = F(ix,iy+1,iz+1,i)[2];  //  dF/dy
    TinyVector<double,2>& Y032 = F(ix,iy+1,iz,i)[6];    // d2F/dydz
    TinyVector<double,2>& Y033 = F(ix,iy+1,iz+1,i)[6];  // d2F/dydz
    
    TinyVector<double,2>& Y100 = F(ix+1,iy,iz,i)[0];      //   F
    TinyVector<double,2>& Y101 = F(ix+1,iy,iz+1,i)[0];    //   F
    TinyVector<double,2>& Y102 = F(ix+1,iy,iz,i)[3];      //  dF/dz
    TinyVector<double,2>& Y103 = F(ix+1,iy,iz+1,i)[3];    //  dF/dz
    TinyVector<double,2>& Y110 = F(ix+1,iy+1,iz,i)[0];    //   F
    TinyVector<double,2>& Y111 = F(ix+1,iy+1,iz+1,i)[0];  //   F
    TinyVector<double,2>& Y112 = F(ix+1,iy+1,iz,i)[3];    //  dF/dz
    TinyVector<double,2>& Y113 = F(ix+1,iy+1,iz+1,i)[3];  //  dF/dz
    TinyVector<double,2>& Y120 = F(ix+1,iy,iz,i)[2];      //  dF/dy
    TinyVector<double,2>& Y121 = F(ix+1,iy,iz+1,i)[2];    //  dF/dy
    TinyVector<double,2>& Y122 = F(ix+1,iy,iz,i)[6];      // d2F/dydz
    TinyVector<double,2>& Y123 = F(ix+1,iy,iz+1,i)[6];    // d2F/dydz
    TinyVector<double,2>& Y130 = F(ix+1,iy+1,iz,i)[2];    //  dF/dy
    TinyVector<double,2>& Y131 = F(ix+1,iy+1,iz+1,i)[2];  //  dF/dy
    TinyVector<double,2>& Y132 = F(ix+1,iy+1,iz,i)[6];    // d2F/dydz
    TinyVector<double,2>& Y133 = F(ix+1,iy+1,iz+1,i)[6];  // d2F/dydz
    
    TinyVector<double,2>& Y200 = F(ix,iy,iz,i)[1];      //  dF/dx
    TinyVector<double,2>& Y201 = F(ix,iy,iz+1,i)[1];    //  dF/dx
    TinyVector<double,2>& Y202 = F(ix,iy,iz,i)[5];      // d2F/dxdz
    TinyVector<double,2>& Y203 = F(ix,iy,iz+1,i)[5];    // d2F/dxdz
    TinyVector<double,2>& Y210 = F(ix,iy+1,iz,i)[1];    //  dF/dx
    TinyVector<double,2>& Y211 = F(ix,iy+1,iz+1,i)[1];  //  dF/dx
    TinyVector<double,2>& Y212 = F(ix,iy+1,iz,i)[5];    // d2F/dxdz
    TinyVector<double,2>& Y213 = F(ix,iy+1,iz+1,i)[5];  // d2F/dxdz
    TinyVector<double,2>& Y220 = F(ix,iy,iz,i)[4];      // d2F/dxdy
    TinyVector<double,2>& Y221 = F(ix,iy,iz+1,i)[4];    // d2F/dxdy
    TinyVector<double,2>& Y222 = F(ix,iy,iz,i)[7];      // d3F/dxdydz
    TinyVector<double,2>& Y223 = F(ix,iy,iz+1,i)[7];    // d3F/dxdydz
    TinyVector<double,2>& Y230 = F(ix,iy+1,iz,i)[4];    // d2F/dxdy
    TinyVector<double,2>& Y231 = F(ix,iy+1,iz+1,i)[4];  // d2F/dxdy
    TinyVector<double,2>& Y232 = F(ix,iy+1,iz,i)[7];    // d3F/dxdydz
    TinyVector<double,2>& Y233 = F(ix,iy+1,iz+1,i)[7];  // d3F/dxdydz
    
    TinyVector<double,2>& Y300 = F(ix+1,iy,iz,i)[1];      //  dF/dx
    TinyVector<double,2>& Y301 = F(ix+1,iy,iz+1,i)[1];    //  dF/dx
    TinyVector<double,2>& Y302 = F(ix+1,iy,iz,i)[5];      // d2F/dxdz
    TinyVector<double,2>& Y303 = F(ix+1,iy,iz+1,i)[5];    // d2F/dxdz
    TinyVector<double,2>& Y310 = F(ix+1,iy+1,iz,i)[1];    //  dF/dx
    TinyVector<double,2>& Y311 = F(ix+1,iy+1,iz+1,i)[1];  //  dF/dx
    TinyVector<double,2>& Y312 = F(ix+1,iy+1,iz,i)[5];    // d2F/dxdz
    TinyVector<double,2>& Y313 = F(ix+1,iy+1,iz+1,i)[5];  // d2F/dxdz
    TinyVector<double,2>& Y320 = F(ix+1,iy,iz,i)[4];      // d2F/dxdy
    TinyVector<double,2>& Y321 = F(ix+1,iy,iz+1,i)[4];    // d2F/dxdy
    TinyVector<double,2>& Y322 = F(ix+1,iy,iz,i)[7];      // d3F/dxdydz
    TinyVector<double,2>& Y323 = F(ix+1,iy,iz+1,i)[7];    // d3F/dxdydz
    TinyVector<double,2>& Y330 = F(ix+1,iy+1,iz,i)[4];    // d2F/dxdy
    TinyVector<double,2>& Y331 = F(ix+1,iy+1,iz+1,i)[4];  // d2F/dxdy
    TinyVector<double,2>& Y332 = F(ix+1,iy+1,iz,i)[7];    // d3F/dxdydz
    TinyVector<double,2>& Y333 = F(ix+1,iy+1,iz+1,i)[7];  // d3F/dxdydz

    
    TinyVector<double,2> Val;
    Val = 
      da0*
      (b0*(Y000*c0+Y001*c1+Y002*c2+Y003*c3) +
       b1*(Y010*c0+Y011*c1+Y012*c2+Y013*c3) +
       b2*(Y020*c0+Y021*c1+Y022*c2+Y023*c3) +
       b3*(Y030*c0+Y031*c1+Y032*c2+Y033*c3))+
      da1 *
      (b0*(Y100*c0+Y101*c1+Y102*c2+Y103*c3) +
       b1*(Y110*c0+Y111*c1+Y112*c2+Y113*c3) +
       b2*(Y120*c0+Y121*c1+Y122*c2+Y123*c3) +
       b3*(Y130*c0+Y131*c1+Y132*c2+Y133*c3))+
      da2 *
      (b0*(Y200*c0+Y201*c1+Y202*c2+Y203*c3) +
       b1*(Y210*c0+Y211*c1+Y212*c2+Y213*c3) +
       b2*(Y220*c0+Y221*c1+Y222*c2+Y223*c3) +
       b3*(Y230*c0+Y231*c1+Y232*c2+Y233*c3))+
      da3 *
      (b0*(Y300*c0+Y301*c1+Y302*c2+Y303*c3) +
       b1*(Y310*c0+Y311*c1+Y312*c2+Y313*c3) +
       b2*(Y320*c0+Y321*c1+Y322*c2+Y323*c3) +
       b3*(Y330*c0+Y331*c1+Y332*c2+Y333*c3));

    vals(i) = complex<double>(Val[0],Val[1]);
  }
}



inline void
ComplexMultiTricubicSpline::d_dy (double x, double y, double z, 
				  Array<complex<double>,1> &vals)
{
  if (!UpToDate)
    Update();

  int ix = Xgrid->ReverseMap(x);  
  int iy = Ygrid->ReverseMap(y);
  int iz = Zgrid->ReverseMap(z);

  ix = max(0,ix); ix = min(ix, Nx-2);
  iy = max(0,iy); iy = min(iy, Ny-2);
  iz = max(0,iz); iz = min(iz, Nz-2);

  double h = (*Xgrid)(ix+1) - (*Xgrid)(ix);
  double k = (*Ygrid)(iy+1) - (*Ygrid)(iy);
  double kinv = 1.0/k;
  double l = (*Zgrid)(iz+1) - (*Zgrid)(iz);
  double u = (x - (*Xgrid)(ix))/h;
  double v = (y - (*Ygrid)(iy))/k;
  double w = (z - (*Zgrid)(iz))/l;

  double a0 = p1(u);
  double a1 = p2(u);
  double a2 = h*q1(u);
  double a3 = h*q2(u);

  register double db0 = kinv*dp1(v);
  register double db1 = kinv*dp2(v);
  register double db2 = dq1(v);
  register double db3 = dq2(v);

  register double c0 = p1(w);
  register double c1 = p2(w);
  register double c2 = l*q1(w);
  register double c3 = l*q2(w);

  for (int i=0; i<N; i++) {
    TinyVector<double,2>& Y000 = F(ix,iy,iz,i)[0];      //   F
    TinyVector<double,2>& Y001 = F(ix,iy,iz+1,i)[0];    //   F
    TinyVector<double,2>& Y002 = F(ix,iy,iz,i)[3];      //  dF/dz
    TinyVector<double,2>& Y003 = F(ix,iy,iz+1,i)[3];    //  dF/dz
    TinyVector<double,2>& Y010 = F(ix,iy+1,iz,i)[0];    //   F
    TinyVector<double,2>& Y011 = F(ix,iy+1,iz+1,i)[0];  //   F
    TinyVector<double,2>& Y012 = F(ix,iy+1,iz,i)[3];    //  dF/dz
    TinyVector<double,2>& Y013 = F(ix,iy+1,iz+1,i)[3];  //  dF/dz
    TinyVector<double,2>& Y020 = F(ix,iy,iz,i)[2];      //  dF/dy
    TinyVector<double,2>& Y021 = F(ix,iy,iz+1,i)[2];    //  dF/dy
    TinyVector<double,2>& Y022 = F(ix,iy,iz,i)[6];      // d2F/dydz
    TinyVector<double,2>& Y023 = F(ix,iy,iz+1,i)[6];    // d2F/dydz
    TinyVector<double,2>& Y030 = F(ix,iy+1,iz,i)[2];    //  dF/dy
    TinyVector<double,2>& Y031 = F(ix,iy+1,iz+1,i)[2];  //  dF/dy
    TinyVector<double,2>& Y032 = F(ix,iy+1,iz,i)[6];    // d2F/dydz
    TinyVector<double,2>& Y033 = F(ix,iy+1,iz+1,i)[6];  // d2F/dydz
    
    TinyVector<double,2>& Y100 = F(ix+1,iy,iz,i)[0];      //   F
    TinyVector<double,2>& Y101 = F(ix+1,iy,iz+1,i)[0];    //   F
    TinyVector<double,2>& Y102 = F(ix+1,iy,iz,i)[3];      //  dF/dz
    TinyVector<double,2>& Y103 = F(ix+1,iy,iz+1,i)[3];    //  dF/dz
    TinyVector<double,2>& Y110 = F(ix+1,iy+1,iz,i)[0];    //   F
    TinyVector<double,2>& Y111 = F(ix+1,iy+1,iz+1,i)[0];  //   F
    TinyVector<double,2>& Y112 = F(ix+1,iy+1,iz,i)[3];    //  dF/dz
    TinyVector<double,2>& Y113 = F(ix+1,iy+1,iz+1,i)[3];  //  dF/dz
    TinyVector<double,2>& Y120 = F(ix+1,iy,iz,i)[2];      //  dF/dy
    TinyVector<double,2>& Y121 = F(ix+1,iy,iz+1,i)[2];    //  dF/dy
    TinyVector<double,2>& Y122 = F(ix+1,iy,iz,i)[6];      // d2F/dydz
    TinyVector<double,2>& Y123 = F(ix+1,iy,iz+1,i)[6];    // d2F/dydz
    TinyVector<double,2>& Y130 = F(ix+1,iy+1,iz,i)[2];    //  dF/dy
    TinyVector<double,2>& Y131 = F(ix+1,iy+1,iz+1,i)[2];  //  dF/dy
    TinyVector<double,2>& Y132 = F(ix+1,iy+1,iz,i)[6];    // d2F/dydz
    TinyVector<double,2>& Y133 = F(ix+1,iy+1,iz+1,i)[6];  // d2F/dydz
    
    TinyVector<double,2>& Y200 = F(ix,iy,iz,i)[1];      //  dF/dx
    TinyVector<double,2>& Y201 = F(ix,iy,iz+1,i)[1];    //  dF/dx
    TinyVector<double,2>& Y202 = F(ix,iy,iz,i)[5];      // d2F/dxdz
    TinyVector<double,2>& Y203 = F(ix,iy,iz+1,i)[5];    // d2F/dxdz
    TinyVector<double,2>& Y210 = F(ix,iy+1,iz,i)[1];    //  dF/dx
    TinyVector<double,2>& Y211 = F(ix,iy+1,iz+1,i)[1];  //  dF/dx
    TinyVector<double,2>& Y212 = F(ix,iy+1,iz,i)[5];    // d2F/dxdz
    TinyVector<double,2>& Y213 = F(ix,iy+1,iz+1,i)[5];  // d2F/dxdz
    TinyVector<double,2>& Y220 = F(ix,iy,iz,i)[4];      // d2F/dxdy
    TinyVector<double,2>& Y221 = F(ix,iy,iz+1,i)[4];    // d2F/dxdy
    TinyVector<double,2>& Y222 = F(ix,iy,iz,i)[7];      // d3F/dxdydz
    TinyVector<double,2>& Y223 = F(ix,iy,iz+1,i)[7];    // d3F/dxdydz
    TinyVector<double,2>& Y230 = F(ix,iy+1,iz,i)[4];    // d2F/dxdy
    TinyVector<double,2>& Y231 = F(ix,iy+1,iz+1,i)[4];  // d2F/dxdy
    TinyVector<double,2>& Y232 = F(ix,iy+1,iz,i)[7];    // d3F/dxdydz
    TinyVector<double,2>& Y233 = F(ix,iy+1,iz+1,i)[7];  // d3F/dxdydz
    
    TinyVector<double,2>& Y300 = F(ix+1,iy,iz,i)[1];      //  dF/dx
    TinyVector<double,2>& Y301 = F(ix+1,iy,iz+1,i)[1];    //  dF/dx
    TinyVector<double,2>& Y302 = F(ix+1,iy,iz,i)[5];      // d2F/dxdz
    TinyVector<double,2>& Y303 = F(ix+1,iy,iz+1,i)[5];    // d2F/dxdz
    TinyVector<double,2>& Y310 = F(ix+1,iy+1,iz,i)[1];    //  dF/dx
    TinyVector<double,2>& Y311 = F(ix+1,iy+1,iz+1,i)[1];  //  dF/dx
    TinyVector<double,2>& Y312 = F(ix+1,iy+1,iz,i)[5];    // d2F/dxdz
    TinyVector<double,2>& Y313 = F(ix+1,iy+1,iz+1,i)[5];  // d2F/dxdz
    TinyVector<double,2>& Y320 = F(ix+1,iy,iz,i)[4];      // d2F/dxdy
    TinyVector<double,2>& Y321 = F(ix+1,iy,iz+1,i)[4];    // d2F/dxdy
    TinyVector<double,2>& Y322 = F(ix+1,iy,iz,i)[7];      // d3F/dxdydz
    TinyVector<double,2>& Y323 = F(ix+1,iy,iz+1,i)[7];    // d3F/dxdydz
    TinyVector<double,2>& Y330 = F(ix+1,iy+1,iz,i)[4];    // d2F/dxdy
    TinyVector<double,2>& Y331 = F(ix+1,iy+1,iz+1,i)[4];  // d2F/dxdy
    TinyVector<double,2>& Y332 = F(ix+1,iy+1,iz,i)[7];    // d3F/dxdydz
    TinyVector<double,2>& Y333 = F(ix+1,iy+1,iz+1,i)[7];  // d3F/dxdydz

    TinyVector<double,2> Val =
      a0*
      (db0*(Y000*c0+Y001*c1+Y002*c2+Y003*c3) +
       db1*(Y010*c0+Y011*c1+Y012*c2+Y013*c3) +
       db2*(Y020*c0+Y021*c1+Y022*c2+Y023*c3) +
       db3*(Y030*c0+Y031*c1+Y032*c2+Y033*c3))+
      a1 *
      (db0*(Y100*c0+Y101*c1+Y102*c2+Y103*c3) +
       db1*(Y110*c0+Y111*c1+Y112*c2+Y113*c3) +
       db2*(Y120*c0+Y121*c1+Y122*c2+Y123*c3) +
       db3*(Y130*c0+Y131*c1+Y132*c2+Y133*c3))+
      a2 *
      (db0*(Y200*c0+Y201*c1+Y202*c2+Y203*c3) +
       db1*(Y210*c0+Y211*c1+Y212*c2+Y213*c3) +
       db2*(Y220*c0+Y221*c1+Y222*c2+Y223*c3) +
       db3*(Y230*c0+Y231*c1+Y232*c2+Y233*c3))+
      a3 *
      (db0*(Y300*c0+Y301*c1+Y302*c2+Y303*c3) +
       db1*(Y310*c0+Y311*c1+Y312*c2+Y313*c3) +
       db2*(Y320*c0+Y321*c1+Y322*c2+Y323*c3) +
       db3*(Y330*c0+Y331*c1+Y332*c2+Y333*c3));

    vals(i) = complex<double> (Val[0], Val[1]);
  }
}

inline void 
ComplexMultiTricubicSpline::d_dz (double x, double y, double z, 
				  Array<complex<double>,1> &vals)
{
  if (!UpToDate)
    Update();

  int ix = Xgrid->ReverseMap(x);  
  int iy = Ygrid->ReverseMap(y);
  int iz = Zgrid->ReverseMap(z);

  ix = max(0,ix); ix = min(ix, Nx-2);
  iy = max(0,iy); iy = min(iy, Ny-2);
  iz = max(0,iz); iz = min(iz, Nz-2);

  double h = (*Xgrid)(ix+1) - (*Xgrid)(ix);
  double k = (*Ygrid)(iy+1) - (*Ygrid)(iy);
  double l = (*Zgrid)(iz+1) - (*Zgrid)(iz);
  double linv = 1.0/l;

  double u = (x - (*Xgrid)(ix))/h;
  double v = (y - (*Ygrid)(iy))/k;
  double w = (z - (*Zgrid)(iz))/l;

  double a0 = p1(u);
  double a1 = p2(u);
  double a2 = h*q1(u);
  double a3 = h*q2(u);

  register double b0 = p1(v);
  register double b1 = p2(v);
  register double b2 = k*q1(v);
  register double b3 = k*q2(v);

  register double dc0 = linv*dp1(w);
  register double dc1 = linv*dp2(w);
  register double dc2 = dq1(w);
  register double dc3 = dq2(w);
  

  for (int i=0; i<N; i++) {
    TinyVector<double,2>& Y000 = F(ix,iy,iz,i)[0];      //   F
    TinyVector<double,2>& Y001 = F(ix,iy,iz+1,i)[0];    //   F
    TinyVector<double,2>& Y002 = F(ix,iy,iz,i)[3];      //  dF/dz
    TinyVector<double,2>& Y003 = F(ix,iy,iz+1,i)[3];    //  dF/dz
    TinyVector<double,2>& Y010 = F(ix,iy+1,iz,i)[0];    //   F
    TinyVector<double,2>& Y011 = F(ix,iy+1,iz+1,i)[0];  //   F
    TinyVector<double,2>& Y012 = F(ix,iy+1,iz,i)[3];    //  dF/dz
    TinyVector<double,2>& Y013 = F(ix,iy+1,iz+1,i)[3];  //  dF/dz
    TinyVector<double,2>& Y020 = F(ix,iy,iz,i)[2];      //  dF/dy
    TinyVector<double,2>& Y021 = F(ix,iy,iz+1,i)[2];    //  dF/dy
    TinyVector<double,2>& Y022 = F(ix,iy,iz,i)[6];      // d2F/dydz
    TinyVector<double,2>& Y023 = F(ix,iy,iz+1,i)[6];    // d2F/dydz
    TinyVector<double,2>& Y030 = F(ix,iy+1,iz,i)[2];    //  dF/dy
    TinyVector<double,2>& Y031 = F(ix,iy+1,iz+1,i)[2];  //  dF/dy
    TinyVector<double,2>& Y032 = F(ix,iy+1,iz,i)[6];    // d2F/dydz
    TinyVector<double,2>& Y033 = F(ix,iy+1,iz+1,i)[6];  // d2F/dydz
    
    TinyVector<double,2>& Y100 = F(ix+1,iy,iz,i)[0];      //   F
    TinyVector<double,2>& Y101 = F(ix+1,iy,iz+1,i)[0];    //   F
    TinyVector<double,2>& Y102 = F(ix+1,iy,iz,i)[3];      //  dF/dz
    TinyVector<double,2>& Y103 = F(ix+1,iy,iz+1,i)[3];    //  dF/dz
    TinyVector<double,2>& Y110 = F(ix+1,iy+1,iz,i)[0];    //   F
    TinyVector<double,2>& Y111 = F(ix+1,iy+1,iz+1,i)[0];  //   F
    TinyVector<double,2>& Y112 = F(ix+1,iy+1,iz,i)[3];    //  dF/dz
    TinyVector<double,2>& Y113 = F(ix+1,iy+1,iz+1,i)[3];  //  dF/dz
    TinyVector<double,2>& Y120 = F(ix+1,iy,iz,i)[2];      //  dF/dy
    TinyVector<double,2>& Y121 = F(ix+1,iy,iz+1,i)[2];    //  dF/dy
    TinyVector<double,2>& Y122 = F(ix+1,iy,iz,i)[6];      // d2F/dydz
    TinyVector<double,2>& Y123 = F(ix+1,iy,iz+1,i)[6];    // d2F/dydz
    TinyVector<double,2>& Y130 = F(ix+1,iy+1,iz,i)[2];    //  dF/dy
    TinyVector<double,2>& Y131 = F(ix+1,iy+1,iz+1,i)[2];  //  dF/dy
    TinyVector<double,2>& Y132 = F(ix+1,iy+1,iz,i)[6];    // d2F/dydz
    TinyVector<double,2>& Y133 = F(ix+1,iy+1,iz+1,i)[6];  // d2F/dydz
    
    TinyVector<double,2>& Y200 = F(ix,iy,iz,i)[1];      //  dF/dx
    TinyVector<double,2>& Y201 = F(ix,iy,iz+1,i)[1];    //  dF/dx
    TinyVector<double,2>& Y202 = F(ix,iy,iz,i)[5];      // d2F/dxdz
    TinyVector<double,2>& Y203 = F(ix,iy,iz+1,i)[5];    // d2F/dxdz
    TinyVector<double,2>& Y210 = F(ix,iy+1,iz,i)[1];    //  dF/dx
    TinyVector<double,2>& Y211 = F(ix,iy+1,iz+1,i)[1];  //  dF/dx
    TinyVector<double,2>& Y212 = F(ix,iy+1,iz,i)[5];    // d2F/dxdz
    TinyVector<double,2>& Y213 = F(ix,iy+1,iz+1,i)[5];  // d2F/dxdz
    TinyVector<double,2>& Y220 = F(ix,iy,iz,i)[4];      // d2F/dxdy
    TinyVector<double,2>& Y221 = F(ix,iy,iz+1,i)[4];    // d2F/dxdy
    TinyVector<double,2>& Y222 = F(ix,iy,iz,i)[7];      // d3F/dxdydz
    TinyVector<double,2>& Y223 = F(ix,iy,iz+1,i)[7];    // d3F/dxdydz
    TinyVector<double,2>& Y230 = F(ix,iy+1,iz,i)[4];    // d2F/dxdy
    TinyVector<double,2>& Y231 = F(ix,iy+1,iz+1,i)[4];  // d2F/dxdy
    TinyVector<double,2>& Y232 = F(ix,iy+1,iz,i)[7];    // d3F/dxdydz
    TinyVector<double,2>& Y233 = F(ix,iy+1,iz+1,i)[7];  // d3F/dxdydz
    
    TinyVector<double,2>& Y300 = F(ix+1,iy,iz,i)[1];      //  dF/dx
    TinyVector<double,2>& Y301 = F(ix+1,iy,iz+1,i)[1];    //  dF/dx
    TinyVector<double,2>& Y302 = F(ix+1,iy,iz,i)[5];      // d2F/dxdz
    TinyVector<double,2>& Y303 = F(ix+1,iy,iz+1,i)[5];    // d2F/dxdz
    TinyVector<double,2>& Y310 = F(ix+1,iy+1,iz,i)[1];    //  dF/dx
    TinyVector<double,2>& Y311 = F(ix+1,iy+1,iz+1,i)[1];  //  dF/dx
    TinyVector<double,2>& Y312 = F(ix+1,iy+1,iz,i)[5];    // d2F/dxdz
    TinyVector<double,2>& Y313 = F(ix+1,iy+1,iz+1,i)[5];  // d2F/dxdz
    TinyVector<double,2>& Y320 = F(ix+1,iy,iz,i)[4];      // d2F/dxdy
    TinyVector<double,2>& Y321 = F(ix+1,iy,iz+1,i)[4];    // d2F/dxdy
    TinyVector<double,2>& Y322 = F(ix+1,iy,iz,i)[7];      // d3F/dxdydz
    TinyVector<double,2>& Y323 = F(ix+1,iy,iz+1,i)[7];    // d3F/dxdydz
    TinyVector<double,2>& Y330 = F(ix+1,iy+1,iz,i)[4];    // d2F/dxdy
    TinyVector<double,2>& Y331 = F(ix+1,iy+1,iz+1,i)[4];  // d2F/dxdy
    TinyVector<double,2>& Y332 = F(ix+1,iy+1,iz,i)[7];    // d3F/dxdydz
    TinyVector<double,2>& Y333 = F(ix+1,iy+1,iz+1,i)[7];  // d3F/dxdydz

    TinyVector<double,2> Val = 
      a0*
      (b0*(Y000*dc0+Y001*dc1+Y002*dc2+Y003*dc3) +
       b1*(Y010*dc0+Y011*dc1+Y012*dc2+Y013*dc3) +
       b2*(Y020*dc0+Y021*dc1+Y022*dc2+Y023*dc3) +
       b3*(Y030*dc0+Y031*dc1+Y032*dc2+Y033*dc3))+
      a1 *
      (b0*(Y100*dc0+Y101*dc1+Y102*dc2+Y103*dc3) +
       b1*(Y110*dc0+Y111*dc1+Y112*dc2+Y113*dc3) +
       b2*(Y120*dc0+Y121*dc1+Y122*dc2+Y123*dc3) +
       b3*(Y130*dc0+Y131*dc1+Y132*dc2+Y133*dc3))+
      a2 *
      (b0*(Y200*dc0+Y201*dc1+Y202*dc2+Y203*dc3) +
       b1*(Y210*dc0+Y211*dc1+Y212*dc2+Y213*dc3) +
       b2*(Y220*dc0+Y221*dc1+Y222*dc2+Y223*dc3) +
       b3*(Y230*dc0+Y231*dc1+Y232*dc2+Y233*dc3))+
      a3 *
      (b0*(Y300*dc0+Y301*dc1+Y302*dc2+Y303*dc3) +
       b1*(Y310*dc0+Y311*dc1+Y312*dc2+Y313*dc3) +
       b2*(Y320*dc0+Y321*dc1+Y322*dc2+Y323*dc3) +
       b3*(Y330*dc0+Y331*dc1+Y332*dc2+Y333*dc3));

    vals(i) = complex<double> (Val[0], Val[1]);
  }
}

inline void
ComplexMultiTricubicSpline::Grad (double x, double y, double z, 
				  Array<cVec3,1> &grads)
{
  if (!UpToDate)
    Update();

  int ix = Xgrid->ReverseMap(x);  
  int iy = Ygrid->ReverseMap(y);
  int iz = Zgrid->ReverseMap(z);

  ix = max(0,ix); ix = min(ix, Nx-2);
  iy = max(0,iy); iy = min(iy, Ny-2);
  iz = max(0,iz); iz = min(iz, Nz-2);

  double h = (*Xgrid)(ix+1) - (*Xgrid)(ix);
  double hinv = 1.0/h;
  double k = (*Ygrid)(iy+1) - (*Ygrid)(iy);
  double kinv = 1.0/k;
  double l = (*Zgrid)(iz+1) - (*Zgrid)(iz);
  double linv = 1.0/l;

  double u = (x - (*Xgrid)(ix))/h;
  double v = (y - (*Ygrid)(iy))/k;
  double w = (z - (*Zgrid)(iz))/l;

  double a0 = p1(u);
  double a1 = p2(u);
  double a2 = h*q1(u);
  double a3 = h*q2(u);
  double da0 = hinv*dp1(u);
  double da1 = hinv*dp2(u);
  double da2 = dq1(u);
  double da3 = dq2(u);

  register double b0 = p1(v);
  register double b1 = p2(v);
  register double b2 = k*q1(v);
  register double b3 = k*q2(v);
  register double db0 = kinv*dp1(v);
  register double db1 = kinv*dp2(v);
  register double db2 = dq1(v);
  register double db3 = dq2(v);

  register double c0 = p1(w);
  register double c1 = p2(w);
  register double c2 = l*q1(w);
  register double c3 = l*q2(w);
  register double dc0 = linv*dp1(w);
  register double dc1 = linv*dp2(w);
  register double dc2 = dq1(w);
  register double dc3 = dq2(w);
  
  for (int i=0; i<N; i++) {
    TinyVector<double,2>& Y000 = F(ix,iy,iz,i)[0];      //   F
    TinyVector<double,2>& Y001 = F(ix,iy,iz+1,i)[0];    //   F
    TinyVector<double,2>& Y002 = F(ix,iy,iz,i)[3];      //  dF/dz
    TinyVector<double,2>& Y003 = F(ix,iy,iz+1,i)[3];    //  dF/dz
    TinyVector<double,2>& Y010 = F(ix,iy+1,iz,i)[0];    //   F
    TinyVector<double,2>& Y011 = F(ix,iy+1,iz+1,i)[0];  //   F
    TinyVector<double,2>& Y012 = F(ix,iy+1,iz,i)[3];    //  dF/dz
    TinyVector<double,2>& Y013 = F(ix,iy+1,iz+1,i)[3];  //  dF/dz
    TinyVector<double,2>& Y020 = F(ix,iy,iz,i)[2];      //  dF/dy
    TinyVector<double,2>& Y021 = F(ix,iy,iz+1,i)[2];    //  dF/dy
    TinyVector<double,2>& Y022 = F(ix,iy,iz,i)[6];      // d2F/dydz
    TinyVector<double,2>& Y023 = F(ix,iy,iz+1,i)[6];    // d2F/dydz
    TinyVector<double,2>& Y030 = F(ix,iy+1,iz,i)[2];    //  dF/dy
    TinyVector<double,2>& Y031 = F(ix,iy+1,iz+1,i)[2];  //  dF/dy
    TinyVector<double,2>& Y032 = F(ix,iy+1,iz,i)[6];    // d2F/dydz
    TinyVector<double,2>& Y033 = F(ix,iy+1,iz+1,i)[6];  // d2F/dydz
    
    TinyVector<double,2>& Y100 = F(ix+1,iy,iz,i)[0];      //   F
    TinyVector<double,2>& Y101 = F(ix+1,iy,iz+1,i)[0];    //   F
    TinyVector<double,2>& Y102 = F(ix+1,iy,iz,i)[3];      //  dF/dz
    TinyVector<double,2>& Y103 = F(ix+1,iy,iz+1,i)[3];    //  dF/dz
    TinyVector<double,2>& Y110 = F(ix+1,iy+1,iz,i)[0];    //   F
    TinyVector<double,2>& Y111 = F(ix+1,iy+1,iz+1,i)[0];  //   F
    TinyVector<double,2>& Y112 = F(ix+1,iy+1,iz,i)[3];    //  dF/dz
    TinyVector<double,2>& Y113 = F(ix+1,iy+1,iz+1,i)[3];  //  dF/dz
    TinyVector<double,2>& Y120 = F(ix+1,iy,iz,i)[2];      //  dF/dy
    TinyVector<double,2>& Y121 = F(ix+1,iy,iz+1,i)[2];    //  dF/dy
    TinyVector<double,2>& Y122 = F(ix+1,iy,iz,i)[6];      // d2F/dydz
    TinyVector<double,2>& Y123 = F(ix+1,iy,iz+1,i)[6];    // d2F/dydz
    TinyVector<double,2>& Y130 = F(ix+1,iy+1,iz,i)[2];    //  dF/dy
    TinyVector<double,2>& Y131 = F(ix+1,iy+1,iz+1,i)[2];  //  dF/dy
    TinyVector<double,2>& Y132 = F(ix+1,iy+1,iz,i)[6];    // d2F/dydz
    TinyVector<double,2>& Y133 = F(ix+1,iy+1,iz+1,i)[6];  // d2F/dydz
    
    TinyVector<double,2>& Y200 = F(ix,iy,iz,i)[1];      //  dF/dx
    TinyVector<double,2>& Y201 = F(ix,iy,iz+1,i)[1];    //  dF/dx
    TinyVector<double,2>& Y202 = F(ix,iy,iz,i)[5];      // d2F/dxdz
    TinyVector<double,2>& Y203 = F(ix,iy,iz+1,i)[5];    // d2F/dxdz
    TinyVector<double,2>& Y210 = F(ix,iy+1,iz,i)[1];    //  dF/dx
    TinyVector<double,2>& Y211 = F(ix,iy+1,iz+1,i)[1];  //  dF/dx
    TinyVector<double,2>& Y212 = F(ix,iy+1,iz,i)[5];    // d2F/dxdz
    TinyVector<double,2>& Y213 = F(ix,iy+1,iz+1,i)[5];  // d2F/dxdz
    TinyVector<double,2>& Y220 = F(ix,iy,iz,i)[4];      // d2F/dxdy
    TinyVector<double,2>& Y221 = F(ix,iy,iz+1,i)[4];    // d2F/dxdy
    TinyVector<double,2>& Y222 = F(ix,iy,iz,i)[7];      // d3F/dxdydz
    TinyVector<double,2>& Y223 = F(ix,iy,iz+1,i)[7];    // d3F/dxdydz
    TinyVector<double,2>& Y230 = F(ix,iy+1,iz,i)[4];    // d2F/dxdy
    TinyVector<double,2>& Y231 = F(ix,iy+1,iz+1,i)[4];  // d2F/dxdy
    TinyVector<double,2>& Y232 = F(ix,iy+1,iz,i)[7];    // d3F/dxdydz
    TinyVector<double,2>& Y233 = F(ix,iy+1,iz+1,i)[7];  // d3F/dxdydz
    
    TinyVector<double,2>& Y300 = F(ix+1,iy,iz,i)[1];      //  dF/dx
    TinyVector<double,2>& Y301 = F(ix+1,iy,iz+1,i)[1];    //  dF/dx
    TinyVector<double,2>& Y302 = F(ix+1,iy,iz,i)[5];      // d2F/dxdz
    TinyVector<double,2>& Y303 = F(ix+1,iy,iz+1,i)[5];    // d2F/dxdz
    TinyVector<double,2>& Y310 = F(ix+1,iy+1,iz,i)[1];    //  dF/dx
    TinyVector<double,2>& Y311 = F(ix+1,iy+1,iz+1,i)[1];  //  dF/dx
    TinyVector<double,2>& Y312 = F(ix+1,iy+1,iz,i)[5];    // d2F/dxdz
    TinyVector<double,2>& Y313 = F(ix+1,iy+1,iz+1,i)[5];  // d2F/dxdz
    TinyVector<double,2>& Y320 = F(ix+1,iy,iz,i)[4];      // d2F/dxdy
    TinyVector<double,2>& Y321 = F(ix+1,iy,iz+1,i)[4];    // d2F/dxdy
    TinyVector<double,2>& Y322 = F(ix+1,iy,iz,i)[7];      // d3F/dxdydz
    TinyVector<double,2>& Y323 = F(ix+1,iy,iz+1,i)[7];    // d3F/dxdydz
    TinyVector<double,2>& Y330 = F(ix+1,iy+1,iz,i)[4];    // d2F/dxdy
    TinyVector<double,2>& Y331 = F(ix+1,iy+1,iz+1,i)[4];  // d2F/dxdy
    TinyVector<double,2>& Y332 = F(ix+1,iy+1,iz,i)[7];    // d3F/dxdydz
    TinyVector<double,2>& Y333 = F(ix+1,iy+1,iz+1,i)[7];  // d3F/dxdydz

    TinyVector<TinyVector<double,2>,3> Grad;
    Grad[0] = 
      da0*
      (b0*(Y000*c0+Y001*c1+Y002*c2+Y003*c3) +
       b1*(Y010*c0+Y011*c1+Y012*c2+Y013*c3) +
       b2*(Y020*c0+Y021*c1+Y022*c2+Y023*c3) +
       b3*(Y030*c0+Y031*c1+Y032*c2+Y033*c3))+
      da1 *
      (b0*(Y100*c0+Y101*c1+Y102*c2+Y103*c3) +
       b1*(Y110*c0+Y111*c1+Y112*c2+Y113*c3) +
       b2*(Y120*c0+Y121*c1+Y122*c2+Y123*c3) +
       b3*(Y130*c0+Y131*c1+Y132*c2+Y133*c3))+
      da2 *
      (b0*(Y200*c0+Y201*c1+Y202*c2+Y203*c3) +
       b1*(Y210*c0+Y211*c1+Y212*c2+Y213*c3) +
       b2*(Y220*c0+Y221*c1+Y222*c2+Y223*c3) +
       b3*(Y230*c0+Y231*c1+Y232*c2+Y233*c3))+
      da3 *
      (b0*(Y300*c0+Y301*c1+Y302*c2+Y303*c3) +
       b1*(Y310*c0+Y311*c1+Y312*c2+Y313*c3) +
       b2*(Y320*c0+Y321*c1+Y322*c2+Y323*c3) +
       b3*(Y330*c0+Y331*c1+Y332*c2+Y333*c3));

    Grad[1] = 
      a0*
      (db0*(Y000*c0+Y001*c1+Y002*c2+Y003*c3) +
       db1*(Y010*c0+Y011*c1+Y012*c2+Y013*c3) +
       db2*(Y020*c0+Y021*c1+Y022*c2+Y023*c3) +
       db3*(Y030*c0+Y031*c1+Y032*c2+Y033*c3))+
      a1 *
      (db0*(Y100*c0+Y101*c1+Y102*c2+Y103*c3) +
       db1*(Y110*c0+Y111*c1+Y112*c2+Y113*c3) +
       db2*(Y120*c0+Y121*c1+Y122*c2+Y123*c3) +
       db3*(Y130*c0+Y131*c1+Y132*c2+Y133*c3))+
      a2 *
      (db0*(Y200*c0+Y201*c1+Y202*c2+Y203*c3) +
       db1*(Y210*c0+Y211*c1+Y212*c2+Y213*c3) +
       db2*(Y220*c0+Y221*c1+Y222*c2+Y223*c3) +
       db3*(Y230*c0+Y231*c1+Y232*c2+Y233*c3))+
      a3 *
      (db0*(Y300*c0+Y301*c1+Y302*c2+Y303*c3) +
       db1*(Y310*c0+Y311*c1+Y312*c2+Y313*c3) +
       db2*(Y320*c0+Y321*c1+Y322*c2+Y323*c3) +
       db3*(Y330*c0+Y331*c1+Y332*c2+Y333*c3));

    Grad[2] = 
      a0*
      (b0*(Y000*dc0+Y001*dc1+Y002*dc2+Y003*dc3) +
       b1*(Y010*dc0+Y011*dc1+Y012*dc2+Y013*dc3) +
       b2*(Y020*dc0+Y021*dc1+Y022*dc2+Y023*dc3) +
       b3*(Y030*dc0+Y031*dc1+Y032*dc2+Y033*dc3))+
      a1 *
      (b0*(Y100*dc0+Y101*dc1+Y102*dc2+Y103*dc3) +
       b1*(Y110*dc0+Y111*dc1+Y112*dc2+Y113*dc3) +
       b2*(Y120*dc0+Y121*dc1+Y122*dc2+Y123*dc3) +
       b3*(Y130*dc0+Y131*dc1+Y132*dc2+Y133*dc3))+
      a2 *
      (b0*(Y200*dc0+Y201*dc1+Y202*dc2+Y203*dc3) +
       b1*(Y210*dc0+Y211*dc1+Y212*dc2+Y213*dc3) +
       b2*(Y220*dc0+Y221*dc1+Y222*dc2+Y223*dc3) +
       b3*(Y230*dc0+Y231*dc1+Y232*dc2+Y233*dc3))+
      a3 *
      (b0*(Y300*dc0+Y301*dc1+Y302*dc2+Y303*dc3) +
       b1*(Y310*dc0+Y311*dc1+Y312*dc2+Y313*dc3) +
       b2*(Y320*dc0+Y321*dc1+Y322*dc2+Y323*dc3) +
       b3*(Y330*dc0+Y331*dc1+Y332*dc2+Y333*dc3));

    grads(i)[0] = complex<double> (Grad[0][0], Grad[0][1]);
    grads(i)[1] = complex<double> (Grad[1][0], Grad[1][1]);
    grads(i)[2] = complex<double> (Grad[2][0], Grad[2][1]);
  }
}

inline void
ComplexMultiTricubicSpline::FValGrad(double x, double y, double z, 
				     Array<complex<double>,1> &vals, 
				     Array<cVec3,1> &grads)
{
  double x0=Xgrid->Start; double y0=Ygrid->Start; double z0=Zgrid->Start;
  double dx=(*Xgrid)(1)-(*Xgrid)(0);
  double dy=(*Ygrid)(1)-(*Ygrid)(0);
  double dz=(*Zgrid)(1)-(*Zgrid)(0);
  F77_Z3VALGRAD(&x,&y,&z,&x0,&dx,&Nx,&y0,&dy,&Ny,&z0,&dz,&Nz,
		F.data(), &N, vals.data(), grads.data());
  
}


/// Returns the value and computes the gradient
inline void
ComplexMultiTricubicSpline::ValGrad(double x, double y, double z, 
				    Array<complex<double>,1> &vals, 
				    Array<cVec3,1> &grads)
{
  if (!UpToDate)
    Update();

  int ix = Xgrid->ReverseMap(x);  
  int iy = Ygrid->ReverseMap(y);
  int iz = Zgrid->ReverseMap(z);

  ix = max(0,ix); ix = min(ix, Nx-2);
  iy = max(0,iy); iy = min(iy, Ny-2);
  iz = max(0,iz); iz = min(iz, Nz-2);



//   double h = (*Xgrid)(ix+1) - (*Xgrid)(ix);
//   double hinv = 1.0/h;
//   double k = (*Ygrid)(iy+1) - (*Ygrid)(iy);
//   double kinv = 1.0/k;
//   double l = (*Zgrid)(iz+1) - (*Zgrid)(iz);
//   double linv = 1.0/l;

  double u = (x - (*Xgrid)(ix))*dxInv;
  double v = (y - (*Ygrid)(iy))*dyInv;
  double w = (z - (*Zgrid)(iz))*dzInv;


  double a0 = p1(u); 
  double a1 = p2(u);
  double a2 = dx*q1(u);
  double a3 = dx*q2(u);
  double da0 = dxInv*dp1(u);
  double da1 = dxInv*dp2(u);
  double da2 = dq1(u);
  double da3 = dq2(u);

  register double b0 = p1(v);
  register double b1 = p2(v);
  register double b2 = dy*q1(v);
  register double b3 = dy*q2(v);
  register double db0 = dyInv*dp1(v);
  register double db1 = dyInv*dp2(v);
  register double db2 = dq1(v);
  register double db3 = dq2(v);

  register double c0 = p1(w);
  register double c1 = p2(w);
  register double c2 = dz*q1(w);
  register double c3 = dz*q2(w);
  register double dc0 = dzInv*dp1(w);
  register double dc1 = dzInv*dp2(w);
  register double dc2 = dq1(w);
  register double dc3 = dq2(w);

//   fprintf (stdout, "ix=%d iy=%d iz=%d\n", ix,iy,iz);
//   fprintf (stdout, "u=%1.15f v=%1.15f w=%1.15f\n", u,v,w);
//   fprintf (stdout, "a0=%1.6f a1=%1.6f a2=%1.6f a3=%1.6f\n", a0,a1,a2,a3);
//   fprintf (stdout, "b0=%1.6f b1=%1.6f b2=%1.6f b3=%1.6f\n", b0,b1,b2,b3);
//   fprintf (stdout, "c0=%1.6f c1=%1.6f c2=%1.6f c3=%1.6f\n", c0,c1,c2,c3);
//   fprintf (stdout, "F(9,9,9,9,7)=(%1.6f, %1.6f)\n", F(9,9,9,9)[7][0],
// 	   F(9,9,9,9)[7][1]);

  
  for (int i=0; i<N; i++) {
    TinyVector<double,2>& Y000 = F(ix,iy,iz,i)[0];      //   F
    TinyVector<double,2>& Y001 = F(ix,iy,iz+1,i)[0];    //   F
    TinyVector<double,2>& Y002 = F(ix,iy,iz,i)[3];      //  dF/dz
    TinyVector<double,2>& Y003 = F(ix,iy,iz+1,i)[3];    //  dF/dz
    TinyVector<double,2>& Y010 = F(ix,iy+1,iz,i)[0];    //   F
    TinyVector<double,2>& Y011 = F(ix,iy+1,iz+1,i)[0];  //   F
    TinyVector<double,2>& Y012 = F(ix,iy+1,iz,i)[3];    //  dF/dz
    TinyVector<double,2>& Y013 = F(ix,iy+1,iz+1,i)[3];  //  dF/dz
    TinyVector<double,2>& Y020 = F(ix,iy,iz,i)[2];      //  dF/dy
    TinyVector<double,2>& Y021 = F(ix,iy,iz+1,i)[2];    //  dF/dy
    TinyVector<double,2>& Y022 = F(ix,iy,iz,i)[6];      // d2F/dydz
    TinyVector<double,2>& Y023 = F(ix,iy,iz+1,i)[6];    // d2F/dydz
    TinyVector<double,2>& Y030 = F(ix,iy+1,iz,i)[2];    //  dF/dy
    TinyVector<double,2>& Y031 = F(ix,iy+1,iz+1,i)[2];  //  dF/dy
    TinyVector<double,2>& Y032 = F(ix,iy+1,iz,i)[6];    // d2F/dydz
    TinyVector<double,2>& Y033 = F(ix,iy+1,iz+1,i)[6];  // d2F/dydz
    
    TinyVector<double,2>& Y100 = F(ix+1,iy,iz,i)[0];      //   F
    TinyVector<double,2>& Y101 = F(ix+1,iy,iz+1,i)[0];    //   F
    TinyVector<double,2>& Y102 = F(ix+1,iy,iz,i)[3];      //  dF/dz
    TinyVector<double,2>& Y103 = F(ix+1,iy,iz+1,i)[3];    //  dF/dz
    TinyVector<double,2>& Y110 = F(ix+1,iy+1,iz,i)[0];    //   F
    TinyVector<double,2>& Y111 = F(ix+1,iy+1,iz+1,i)[0];  //   F
    TinyVector<double,2>& Y112 = F(ix+1,iy+1,iz,i)[3];    //  dF/dz
    TinyVector<double,2>& Y113 = F(ix+1,iy+1,iz+1,i)[3];  //  dF/dz
    TinyVector<double,2>& Y120 = F(ix+1,iy,iz,i)[2];      //  dF/dy
    TinyVector<double,2>& Y121 = F(ix+1,iy,iz+1,i)[2];    //  dF/dy
    TinyVector<double,2>& Y122 = F(ix+1,iy,iz,i)[6];      // d2F/dydz
    TinyVector<double,2>& Y123 = F(ix+1,iy,iz+1,i)[6];    // d2F/dydz
    TinyVector<double,2>& Y130 = F(ix+1,iy+1,iz,i)[2];    //  dF/dy
    TinyVector<double,2>& Y131 = F(ix+1,iy+1,iz+1,i)[2];  //  dF/dy
    TinyVector<double,2>& Y132 = F(ix+1,iy+1,iz,i)[6];    // d2F/dydz
    TinyVector<double,2>& Y133 = F(ix+1,iy+1,iz+1,i)[6];  // d2F/dydz
    
    TinyVector<double,2>& Y200 = F(ix,iy,iz,i)[1];      //  dF/dx
    TinyVector<double,2>& Y201 = F(ix,iy,iz+1,i)[1];    //  dF/dx
    TinyVector<double,2>& Y202 = F(ix,iy,iz,i)[5];      // d2F/dxdz
    TinyVector<double,2>& Y203 = F(ix,iy,iz+1,i)[5];    // d2F/dxdz
    TinyVector<double,2>& Y210 = F(ix,iy+1,iz,i)[1];    //  dF/dx
    TinyVector<double,2>& Y211 = F(ix,iy+1,iz+1,i)[1];  //  dF/dx
    TinyVector<double,2>& Y212 = F(ix,iy+1,iz,i)[5];    // d2F/dxdz
    TinyVector<double,2>& Y213 = F(ix,iy+1,iz+1,i)[5];  // d2F/dxdz
    TinyVector<double,2>& Y220 = F(ix,iy,iz,i)[4];      // d2F/dxdy
    TinyVector<double,2>& Y221 = F(ix,iy,iz+1,i)[4];    // d2F/dxdy
    TinyVector<double,2>& Y222 = F(ix,iy,iz,i)[7];      // d3F/dxdydz
    TinyVector<double,2>& Y223 = F(ix,iy,iz+1,i)[7];    // d3F/dxdydz
    TinyVector<double,2>& Y230 = F(ix,iy+1,iz,i)[4];    // d2F/dxdy
    TinyVector<double,2>& Y231 = F(ix,iy+1,iz+1,i)[4];  // d2F/dxdy
    TinyVector<double,2>& Y232 = F(ix,iy+1,iz,i)[7];    // d3F/dxdydz
    TinyVector<double,2>& Y233 = F(ix,iy+1,iz+1,i)[7];  // d3F/dxdydz
    
    TinyVector<double,2>& Y300 = F(ix+1,iy,iz,i)[1];      //  dF/dx
    TinyVector<double,2>& Y301 = F(ix+1,iy,iz+1,i)[1];    //  dF/dx
    TinyVector<double,2>& Y302 = F(ix+1,iy,iz,i)[5];      // d2F/dxdz
    TinyVector<double,2>& Y303 = F(ix+1,iy,iz+1,i)[5];    // d2F/dxdz
    TinyVector<double,2>& Y310 = F(ix+1,iy+1,iz,i)[1];    //  dF/dx
    TinyVector<double,2>& Y311 = F(ix+1,iy+1,iz+1,i)[1];  //  dF/dx
    TinyVector<double,2>& Y312 = F(ix+1,iy+1,iz,i)[5];    // d2F/dxdz
    TinyVector<double,2>& Y313 = F(ix+1,iy+1,iz+1,i)[5];  // d2F/dxdz
    TinyVector<double,2>& Y320 = F(ix+1,iy,iz,i)[4];      // d2F/dxdy
    TinyVector<double,2>& Y321 = F(ix+1,iy,iz+1,i)[4];    // d2F/dxdy
    TinyVector<double,2>& Y322 = F(ix+1,iy,iz,i)[7];      // d3F/dxdydz
    TinyVector<double,2>& Y323 = F(ix+1,iy,iz+1,i)[7];    // d3F/dxdydz
    TinyVector<double,2>& Y330 = F(ix+1,iy+1,iz,i)[4];    // d2F/dxdy
    TinyVector<double,2>& Y331 = F(ix+1,iy+1,iz+1,i)[4];  // d2F/dxdy
    TinyVector<double,2>& Y332 = F(ix+1,iy+1,iz,i)[7];    // d3F/dxdydz
    TinyVector<double,2>& Y333 = F(ix+1,iy+1,iz+1,i)[7];  // d3F/dxdydz
  
    TinyVector<double,2> Val;
    TinyVector<TinyVector<double,2>,3> Grad;
    Val = 
      a0*
      (b0*(Y000*c0+Y001*c1+Y002*c2+Y003*c3) +
       b1*(Y010*c0+Y011*c1+Y012*c2+Y013*c3) +
       b2*(Y020*c0+Y021*c1+Y022*c2+Y023*c3) +
       b3*(Y030*c0+Y031*c1+Y032*c2+Y033*c3))+
      a1 *
      (b0*(Y100*c0+Y101*c1+Y102*c2+Y103*c3) +
       b1*(Y110*c0+Y111*c1+Y112*c2+Y113*c3) +
       b2*(Y120*c0+Y121*c1+Y122*c2+Y123*c3) +
       b3*(Y130*c0+Y131*c1+Y132*c2+Y133*c3))+
      a2 *
      (b0*(Y200*c0+Y201*c1+Y202*c2+Y203*c3) +
       b1*(Y210*c0+Y211*c1+Y212*c2+Y213*c3) +
       b2*(Y220*c0+Y221*c1+Y222*c2+Y223*c3) +
       b3*(Y230*c0+Y231*c1+Y232*c2+Y233*c3))+
      a3 *
      (b0*(Y300*c0+Y301*c1+Y302*c2+Y303*c3) +
       b1*(Y310*c0+Y311*c1+Y312*c2+Y313*c3) +
       b2*(Y320*c0+Y321*c1+Y322*c2+Y323*c3) +
       b3*(Y330*c0+Y331*c1+Y332*c2+Y333*c3));

    vals(i) = complex<double> (Val[0], Val[1]);


    Grad[0] = 
      da0*
      (b0*(Y000*c0+Y001*c1+Y002*c2+Y003*c3) +
       b1*(Y010*c0+Y011*c1+Y012*c2+Y013*c3) +
       b2*(Y020*c0+Y021*c1+Y022*c2+Y023*c3) +
       b3*(Y030*c0+Y031*c1+Y032*c2+Y033*c3))+
      da1 *
      (b0*(Y100*c0+Y101*c1+Y102*c2+Y103*c3) +
       b1*(Y110*c0+Y111*c1+Y112*c2+Y113*c3) +
       b2*(Y120*c0+Y121*c1+Y122*c2+Y123*c3) +
       b3*(Y130*c0+Y131*c1+Y132*c2+Y133*c3))+
      da2 *
      (b0*(Y200*c0+Y201*c1+Y202*c2+Y203*c3) +
       b1*(Y210*c0+Y211*c1+Y212*c2+Y213*c3) +
       b2*(Y220*c0+Y221*c1+Y222*c2+Y223*c3) +
       b3*(Y230*c0+Y231*c1+Y232*c2+Y233*c3))+
      da3 *
      (b0*(Y300*c0+Y301*c1+Y302*c2+Y303*c3) +
       b1*(Y310*c0+Y311*c1+Y312*c2+Y313*c3) +
       b2*(Y320*c0+Y321*c1+Y322*c2+Y323*c3) +
       b3*(Y330*c0+Y331*c1+Y332*c2+Y333*c3));

    Grad[1] = 
      a0*
      (db0*(Y000*c0+Y001*c1+Y002*c2+Y003*c3) +
       db1*(Y010*c0+Y011*c1+Y012*c2+Y013*c3) +
       db2*(Y020*c0+Y021*c1+Y022*c2+Y023*c3) +
       db3*(Y030*c0+Y031*c1+Y032*c2+Y033*c3))+
      a1 *
      (db0*(Y100*c0+Y101*c1+Y102*c2+Y103*c3) +
       db1*(Y110*c0+Y111*c1+Y112*c2+Y113*c3) +
       db2*(Y120*c0+Y121*c1+Y122*c2+Y123*c3) +
       db3*(Y130*c0+Y131*c1+Y132*c2+Y133*c3))+
      a2 *
      (db0*(Y200*c0+Y201*c1+Y202*c2+Y203*c3) +
       db1*(Y210*c0+Y211*c1+Y212*c2+Y213*c3) +
       db2*(Y220*c0+Y221*c1+Y222*c2+Y223*c3) +
       db3*(Y230*c0+Y231*c1+Y232*c2+Y233*c3))+
      a3 *
      (db0*(Y300*c0+Y301*c1+Y302*c2+Y303*c3) +
       db1*(Y310*c0+Y311*c1+Y312*c2+Y313*c3) +
       db2*(Y320*c0+Y321*c1+Y322*c2+Y323*c3) +
       db3*(Y330*c0+Y331*c1+Y332*c2+Y333*c3));

    Grad[2] = 
      a0*
      (b0*(Y000*dc0+Y001*dc1+Y002*dc2+Y003*dc3) +
       b1*(Y010*dc0+Y011*dc1+Y012*dc2+Y013*dc3) +
       b2*(Y020*dc0+Y021*dc1+Y022*dc2+Y023*dc3) +
       b3*(Y030*dc0+Y031*dc1+Y032*dc2+Y033*dc3))+
      a1 *
      (b0*(Y100*dc0+Y101*dc1+Y102*dc2+Y103*dc3) +
       b1*(Y110*dc0+Y111*dc1+Y112*dc2+Y113*dc3) +
       b2*(Y120*dc0+Y121*dc1+Y122*dc2+Y123*dc3) +
       b3*(Y130*dc0+Y131*dc1+Y132*dc2+Y133*dc3))+
      a2 *
      (b0*(Y200*dc0+Y201*dc1+Y202*dc2+Y203*dc3) +
       b1*(Y210*dc0+Y211*dc1+Y212*dc2+Y213*dc3) +
       b2*(Y220*dc0+Y221*dc1+Y222*dc2+Y223*dc3) +
       b3*(Y230*dc0+Y231*dc1+Y232*dc2+Y233*dc3))+
      a3 *
      (b0*(Y300*dc0+Y301*dc1+Y302*dc2+Y303*dc3) +
       b1*(Y310*dc0+Y311*dc1+Y312*dc2+Y313*dc3) +
       b2*(Y320*dc0+Y321*dc1+Y322*dc2+Y323*dc3) +
       b3*(Y330*dc0+Y331*dc1+Y332*dc2+Y333*dc3));
    grads(i)[0] = complex<double>(Grad[0][0], Grad[0][1]);
    grads(i)[1] = complex<double>(Grad[1][0], Grad[1][1]);
    grads(i)[2] = complex<double>(Grad[2][0], Grad[2][1]);
  }
}


inline void 
ComplexMultiTricubicSpline::d2_dx2 (double x, double y, double z, 
				    Array<complex<double>,1> &vals)
{
  if (!UpToDate)
    Update();

  int ix = Xgrid->ReverseMap(x);  
  int iy = Ygrid->ReverseMap(y);
  int iz = Zgrid->ReverseMap(z);

  ix = max(0,ix); ix = min(ix, Nx-2);
  iy = max(0,iy); iy = min(iy, Ny-2);
  iz = max(0,iz); iz = min(iz, Nz-2);

  double h = (*Xgrid)(ix+1) - (*Xgrid)(ix);
  double hinv = 1.0/h;
  double k = (*Ygrid)(iy+1) - (*Ygrid)(iy);
  double l = (*Zgrid)(iz+1) - (*Zgrid)(iz);
  double u = (x - (*Xgrid)(ix))/h;
  double v = (y - (*Ygrid)(iy))/k;
  double w = (z - (*Zgrid)(iz))/l;

  double d2a0 = hinv*hinv*d2p1(u);
  double d2a1 = hinv*hinv*d2p2(u);
  double d2a2 = hinv*d2q1(u);
  double d2a3 = hinv*d2q2(u);

  register double b0 = p1(v);
  register double b1 = p2(v);
  register double b2 = k*q1(v);
  register double b3 = k*q2(v);

  register double c0 = p1(w);
  register double c1 = p2(w);
  register double c2 = l*q1(w);
  register double c3 = l*q2(w);
  
  for (int i=0; i<N; i++) {
    TinyVector<double,2>& Y000 = F(ix,iy,iz,i)[0];      //   F
    TinyVector<double,2>& Y001 = F(ix,iy,iz+1,i)[0];    //   F
    TinyVector<double,2>& Y002 = F(ix,iy,iz,i)[3];      //  dF/dz
    TinyVector<double,2>& Y003 = F(ix,iy,iz+1,i)[3];    //  dF/dz
    TinyVector<double,2>& Y010 = F(ix,iy+1,iz,i)[0];    //   F
    TinyVector<double,2>& Y011 = F(ix,iy+1,iz+1,i)[0];  //   F
    TinyVector<double,2>& Y012 = F(ix,iy+1,iz,i)[3];    //  dF/dz
    TinyVector<double,2>& Y013 = F(ix,iy+1,iz+1,i)[3];  //  dF/dz
    TinyVector<double,2>& Y020 = F(ix,iy,iz,i)[2];      //  dF/dy
    TinyVector<double,2>& Y021 = F(ix,iy,iz+1,i)[2];    //  dF/dy
    TinyVector<double,2>& Y022 = F(ix,iy,iz,i)[6];      // d2F/dydz
    TinyVector<double,2>& Y023 = F(ix,iy,iz+1,i)[6];    // d2F/dydz
    TinyVector<double,2>& Y030 = F(ix,iy+1,iz,i)[2];    //  dF/dy
    TinyVector<double,2>& Y031 = F(ix,iy+1,iz+1,i)[2];  //  dF/dy
    TinyVector<double,2>& Y032 = F(ix,iy+1,iz,i)[6];    // d2F/dydz
    TinyVector<double,2>& Y033 = F(ix,iy+1,iz+1,i)[6];  // d2F/dydz
    
    TinyVector<double,2>& Y100 = F(ix+1,iy,iz,i)[0];      //   F
    TinyVector<double,2>& Y101 = F(ix+1,iy,iz+1,i)[0];    //   F
    TinyVector<double,2>& Y102 = F(ix+1,iy,iz,i)[3];      //  dF/dz
    TinyVector<double,2>& Y103 = F(ix+1,iy,iz+1,i)[3];    //  dF/dz
    TinyVector<double,2>& Y110 = F(ix+1,iy+1,iz,i)[0];    //   F
    TinyVector<double,2>& Y111 = F(ix+1,iy+1,iz+1,i)[0];  //   F
    TinyVector<double,2>& Y112 = F(ix+1,iy+1,iz,i)[3];    //  dF/dz
    TinyVector<double,2>& Y113 = F(ix+1,iy+1,iz+1,i)[3];  //  dF/dz
    TinyVector<double,2>& Y120 = F(ix+1,iy,iz,i)[2];      //  dF/dy
    TinyVector<double,2>& Y121 = F(ix+1,iy,iz+1,i)[2];    //  dF/dy
    TinyVector<double,2>& Y122 = F(ix+1,iy,iz,i)[6];      // d2F/dydz
    TinyVector<double,2>& Y123 = F(ix+1,iy,iz+1,i)[6];    // d2F/dydz
    TinyVector<double,2>& Y130 = F(ix+1,iy+1,iz,i)[2];    //  dF/dy
    TinyVector<double,2>& Y131 = F(ix+1,iy+1,iz+1,i)[2];  //  dF/dy
    TinyVector<double,2>& Y132 = F(ix+1,iy+1,iz,i)[6];    // d2F/dydz
    TinyVector<double,2>& Y133 = F(ix+1,iy+1,iz+1,i)[6];  // d2F/dydz
    
    TinyVector<double,2>& Y200 = F(ix,iy,iz,i)[1];      //  dF/dx
    TinyVector<double,2>& Y201 = F(ix,iy,iz+1,i)[1];    //  dF/dx
    TinyVector<double,2>& Y202 = F(ix,iy,iz,i)[5];      // d2F/dxdz
    TinyVector<double,2>& Y203 = F(ix,iy,iz+1,i)[5];    // d2F/dxdz
    TinyVector<double,2>& Y210 = F(ix,iy+1,iz,i)[1];    //  dF/dx
    TinyVector<double,2>& Y211 = F(ix,iy+1,iz+1,i)[1];  //  dF/dx
    TinyVector<double,2>& Y212 = F(ix,iy+1,iz,i)[5];    // d2F/dxdz
    TinyVector<double,2>& Y213 = F(ix,iy+1,iz+1,i)[5];  // d2F/dxdz
    TinyVector<double,2>& Y220 = F(ix,iy,iz,i)[4];      // d2F/dxdy
    TinyVector<double,2>& Y221 = F(ix,iy,iz+1,i)[4];    // d2F/dxdy
    TinyVector<double,2>& Y222 = F(ix,iy,iz,i)[7];      // d3F/dxdydz
    TinyVector<double,2>& Y223 = F(ix,iy,iz+1,i)[7];    // d3F/dxdydz
    TinyVector<double,2>& Y230 = F(ix,iy+1,iz,i)[4];    // d2F/dxdy
    TinyVector<double,2>& Y231 = F(ix,iy+1,iz+1,i)[4];  // d2F/dxdy
    TinyVector<double,2>& Y232 = F(ix,iy+1,iz,i)[7];    // d3F/dxdydz
    TinyVector<double,2>& Y233 = F(ix,iy+1,iz+1,i)[7];  // d3F/dxdydz
    
    TinyVector<double,2>& Y300 = F(ix+1,iy,iz,i)[1];      //  dF/dx
    TinyVector<double,2>& Y301 = F(ix+1,iy,iz+1,i)[1];    //  dF/dx
    TinyVector<double,2>& Y302 = F(ix+1,iy,iz,i)[5];      // d2F/dxdz
    TinyVector<double,2>& Y303 = F(ix+1,iy,iz+1,i)[5];    // d2F/dxdz
    TinyVector<double,2>& Y310 = F(ix+1,iy+1,iz,i)[1];    //  dF/dx
    TinyVector<double,2>& Y311 = F(ix+1,iy+1,iz+1,i)[1];  //  dF/dx
    TinyVector<double,2>& Y312 = F(ix+1,iy+1,iz,i)[5];    // d2F/dxdz
    TinyVector<double,2>& Y313 = F(ix+1,iy+1,iz+1,i)[5];  // d2F/dxdz
    TinyVector<double,2>& Y320 = F(ix+1,iy,iz,i)[4];      // d2F/dxdy
    TinyVector<double,2>& Y321 = F(ix+1,iy,iz+1,i)[4];    // d2F/dxdy
    TinyVector<double,2>& Y322 = F(ix+1,iy,iz,i)[7];      // d3F/dxdydz
    TinyVector<double,2>& Y323 = F(ix+1,iy,iz+1,i)[7];    // d3F/dxdydz
    TinyVector<double,2>& Y330 = F(ix+1,iy+1,iz,i)[4];    // d2F/dxdy
    TinyVector<double,2>& Y331 = F(ix+1,iy+1,iz+1,i)[4];  // d2F/dxdy
    TinyVector<double,2>& Y332 = F(ix+1,iy+1,iz,i)[7];    // d3F/dxdydz
    TinyVector<double,2>& Y333 = F(ix+1,iy+1,iz+1,i)[7];  // d3F/dxdydz
  
    TinyVector<double,2> Val =
      d2a0*
      (b0*(Y000*c0+Y001*c1+Y002*c2+Y003*c3) +
       b1*(Y010*c0+Y011*c1+Y012*c2+Y013*c3) +
       b2*(Y020*c0+Y021*c1+Y022*c2+Y023*c3) +
       b3*(Y030*c0+Y031*c1+Y032*c2+Y033*c3))+
      d2a1 *
      (b0*(Y100*c0+Y101*c1+Y102*c2+Y103*c3) +
       b1*(Y110*c0+Y111*c1+Y112*c2+Y113*c3) +
       b2*(Y120*c0+Y121*c1+Y122*c2+Y123*c3) +
       b3*(Y130*c0+Y131*c1+Y132*c2+Y133*c3))+
      d2a2 *
      (b0*(Y200*c0+Y201*c1+Y202*c2+Y203*c3) +
       b1*(Y210*c0+Y211*c1+Y212*c2+Y213*c3) +
       b2*(Y220*c0+Y221*c1+Y222*c2+Y223*c3) +
       b3*(Y230*c0+Y231*c1+Y232*c2+Y233*c3))+
      d2a3 *
      (b0*(Y300*c0+Y301*c1+Y302*c2+Y303*c3) +
       b1*(Y310*c0+Y311*c1+Y312*c2+Y313*c3) +
       b2*(Y320*c0+Y321*c1+Y322*c2+Y323*c3) +
       b3*(Y330*c0+Y331*c1+Y332*c2+Y333*c3));

    vals(i) = complex<double> (Val[0], Val[1]);
  }
}


inline void 
ComplexMultiTricubicSpline::d2_dy2 (double x, double y, double z, 
				    Array<complex<double>,1> &vals)
{
  if (!UpToDate)
    Update();

  int ix = Xgrid->ReverseMap(x);  
  int iy = Ygrid->ReverseMap(y);
  int iz = Zgrid->ReverseMap(z);

  ix = max(0,ix); ix = min(ix, Nx-2);
  iy = max(0,iy); iy = min(iy, Ny-2);
  iz = max(0,iz); iz = min(iz, Nz-2);

  double h = (*Xgrid)(ix+1) - (*Xgrid)(ix);
  double k = (*Ygrid)(iy+1) - (*Ygrid)(iy);
  double kinv = 1.0/k;
  double l = (*Zgrid)(iz+1) - (*Zgrid)(iz);
  double u = (x - (*Xgrid)(ix))/h;
  double v = (y - (*Ygrid)(iy))/k;
  double w = (z - (*Zgrid)(iz))/l;

  double a0 = p1(u);
  double a1 = p2(u);
  double a2 = h*q1(u);
  double a3 = h*q2(u);

  register double d2b0 = kinv*kinv*d2p1(v);
  register double d2b1 = kinv*kinv*d2p2(v);
  register double d2b2 = kinv*d2q1(v);
  register double d2b3 = kinv*d2q2(v);

  register double c0 = p1(w);
  register double c1 = p2(w);
  register double c2 = l*q1(w);
  register double c3 = l*q2(w);
  
  for (int i=0; i<N; i++) {
    TinyVector<double,2>& Y000 = F(ix,iy,iz,i)[0];      //   F
    TinyVector<double,2>& Y001 = F(ix,iy,iz+1,i)[0];    //   F
    TinyVector<double,2>& Y002 = F(ix,iy,iz,i)[3];      //  dF/dz
    TinyVector<double,2>& Y003 = F(ix,iy,iz+1,i)[3];    //  dF/dz
    TinyVector<double,2>& Y010 = F(ix,iy+1,iz,i)[0];    //   F
    TinyVector<double,2>& Y011 = F(ix,iy+1,iz+1,i)[0];  //   F
    TinyVector<double,2>& Y012 = F(ix,iy+1,iz,i)[3];    //  dF/dz
    TinyVector<double,2>& Y013 = F(ix,iy+1,iz+1,i)[3];  //  dF/dz
    TinyVector<double,2>& Y020 = F(ix,iy,iz,i)[2];      //  dF/dy
    TinyVector<double,2>& Y021 = F(ix,iy,iz+1,i)[2];    //  dF/dy
    TinyVector<double,2>& Y022 = F(ix,iy,iz,i)[6];      // d2F/dydz
    TinyVector<double,2>& Y023 = F(ix,iy,iz+1,i)[6];    // d2F/dydz
    TinyVector<double,2>& Y030 = F(ix,iy+1,iz,i)[2];    //  dF/dy
    TinyVector<double,2>& Y031 = F(ix,iy+1,iz+1,i)[2];  //  dF/dy
    TinyVector<double,2>& Y032 = F(ix,iy+1,iz,i)[6];    // d2F/dydz
    TinyVector<double,2>& Y033 = F(ix,iy+1,iz+1,i)[6];  // d2F/dydz
    
    TinyVector<double,2>& Y100 = F(ix+1,iy,iz,i)[0];      //   F
    TinyVector<double,2>& Y101 = F(ix+1,iy,iz+1,i)[0];    //   F
    TinyVector<double,2>& Y102 = F(ix+1,iy,iz,i)[3];      //  dF/dz
    TinyVector<double,2>& Y103 = F(ix+1,iy,iz+1,i)[3];    //  dF/dz
    TinyVector<double,2>& Y110 = F(ix+1,iy+1,iz,i)[0];    //   F
    TinyVector<double,2>& Y111 = F(ix+1,iy+1,iz+1,i)[0];  //   F
    TinyVector<double,2>& Y112 = F(ix+1,iy+1,iz,i)[3];    //  dF/dz
    TinyVector<double,2>& Y113 = F(ix+1,iy+1,iz+1,i)[3];  //  dF/dz
    TinyVector<double,2>& Y120 = F(ix+1,iy,iz,i)[2];      //  dF/dy
    TinyVector<double,2>& Y121 = F(ix+1,iy,iz+1,i)[2];    //  dF/dy
    TinyVector<double,2>& Y122 = F(ix+1,iy,iz,i)[6];      // d2F/dydz
    TinyVector<double,2>& Y123 = F(ix+1,iy,iz+1,i)[6];    // d2F/dydz
    TinyVector<double,2>& Y130 = F(ix+1,iy+1,iz,i)[2];    //  dF/dy
    TinyVector<double,2>& Y131 = F(ix+1,iy+1,iz+1,i)[2];  //  dF/dy
    TinyVector<double,2>& Y132 = F(ix+1,iy+1,iz,i)[6];    // d2F/dydz
    TinyVector<double,2>& Y133 = F(ix+1,iy+1,iz+1,i)[6];  // d2F/dydz
    
    TinyVector<double,2>& Y200 = F(ix,iy,iz,i)[1];      //  dF/dx
    TinyVector<double,2>& Y201 = F(ix,iy,iz+1,i)[1];    //  dF/dx
    TinyVector<double,2>& Y202 = F(ix,iy,iz,i)[5];      // d2F/dxdz
    TinyVector<double,2>& Y203 = F(ix,iy,iz+1,i)[5];    // d2F/dxdz
    TinyVector<double,2>& Y210 = F(ix,iy+1,iz,i)[1];    //  dF/dx
    TinyVector<double,2>& Y211 = F(ix,iy+1,iz+1,i)[1];  //  dF/dx
    TinyVector<double,2>& Y212 = F(ix,iy+1,iz,i)[5];    // d2F/dxdz
    TinyVector<double,2>& Y213 = F(ix,iy+1,iz+1,i)[5];  // d2F/dxdz
    TinyVector<double,2>& Y220 = F(ix,iy,iz,i)[4];      // d2F/dxdy
    TinyVector<double,2>& Y221 = F(ix,iy,iz+1,i)[4];    // d2F/dxdy
    TinyVector<double,2>& Y222 = F(ix,iy,iz,i)[7];      // d3F/dxdydz
    TinyVector<double,2>& Y223 = F(ix,iy,iz+1,i)[7];    // d3F/dxdydz
    TinyVector<double,2>& Y230 = F(ix,iy+1,iz,i)[4];    // d2F/dxdy
    TinyVector<double,2>& Y231 = F(ix,iy+1,iz+1,i)[4];  // d2F/dxdy
    TinyVector<double,2>& Y232 = F(ix,iy+1,iz,i)[7];    // d3F/dxdydz
    TinyVector<double,2>& Y233 = F(ix,iy+1,iz+1,i)[7];  // d3F/dxdydz
    
    TinyVector<double,2>& Y300 = F(ix+1,iy,iz,i)[1];      //  dF/dx
    TinyVector<double,2>& Y301 = F(ix+1,iy,iz+1,i)[1];    //  dF/dx
    TinyVector<double,2>& Y302 = F(ix+1,iy,iz,i)[5];      // d2F/dxdz
    TinyVector<double,2>& Y303 = F(ix+1,iy,iz+1,i)[5];    // d2F/dxdz
    TinyVector<double,2>& Y310 = F(ix+1,iy+1,iz,i)[1];    //  dF/dx
    TinyVector<double,2>& Y311 = F(ix+1,iy+1,iz+1,i)[1];  //  dF/dx
    TinyVector<double,2>& Y312 = F(ix+1,iy+1,iz,i)[5];    // d2F/dxdz
    TinyVector<double,2>& Y313 = F(ix+1,iy+1,iz+1,i)[5];  // d2F/dxdz
    TinyVector<double,2>& Y320 = F(ix+1,iy,iz,i)[4];      // d2F/dxdy
    TinyVector<double,2>& Y321 = F(ix+1,iy,iz+1,i)[4];    // d2F/dxdy
    TinyVector<double,2>& Y322 = F(ix+1,iy,iz,i)[7];      // d3F/dxdydz
    TinyVector<double,2>& Y323 = F(ix+1,iy,iz+1,i)[7];    // d3F/dxdydz
    TinyVector<double,2>& Y330 = F(ix+1,iy+1,iz,i)[4];    // d2F/dxdy
    TinyVector<double,2>& Y331 = F(ix+1,iy+1,iz+1,i)[4];  // d2F/dxdy
    TinyVector<double,2>& Y332 = F(ix+1,iy+1,iz,i)[7];    // d3F/dxdydz
    TinyVector<double,2>& Y333 = F(ix+1,iy+1,iz+1,i)[7];  // d3F/dxdydz


    TinyVector<double,2> Val = 
      a0*
      (d2b0*(Y000*c0+Y001*c1+Y002*c2+Y003*c3) +
       d2b1*(Y010*c0+Y011*c1+Y012*c2+Y013*c3) +
       d2b2*(Y020*c0+Y021*c1+Y022*c2+Y023*c3) +
       d2b3*(Y030*c0+Y031*c1+Y032*c2+Y033*c3))+
      a1 *
      (d2b0*(Y100*c0+Y101*c1+Y102*c2+Y103*c3) +
       d2b1*(Y110*c0+Y111*c1+Y112*c2+Y113*c3) +
       d2b2*(Y120*c0+Y121*c1+Y122*c2+Y123*c3) +
       d2b3*(Y130*c0+Y131*c1+Y132*c2+Y133*c3))+
      a2 *
      (d2b0*(Y200*c0+Y201*c1+Y202*c2+Y203*c3) +
       d2b1*(Y210*c0+Y211*c1+Y212*c2+Y213*c3) +
       d2b2*(Y220*c0+Y221*c1+Y222*c2+Y223*c3) +
       d2b3*(Y230*c0+Y231*c1+Y232*c2+Y233*c3))+
      a3 *
      (d2b0*(Y300*c0+Y301*c1+Y302*c2+Y303*c3) +
       d2b1*(Y310*c0+Y311*c1+Y312*c2+Y313*c3) +
       d2b2*(Y320*c0+Y321*c1+Y322*c2+Y323*c3) +
       d2b3*(Y330*c0+Y331*c1+Y332*c2+Y333*c3));

      vals(i) = complex<double> (Val[0], Val[1]);
  }
}


inline void 
ComplexMultiTricubicSpline::d2_dz2 (double x, double y, double z, 
				    Array<complex<double>,1> &vals)
{
  if (!UpToDate)
    Update();

  int ix = Xgrid->ReverseMap(x);  
  int iy = Ygrid->ReverseMap(y);
  int iz = Zgrid->ReverseMap(z);

  ix = max(0,ix); ix = min(ix, Nx-2);
  iy = max(0,iy); iy = min(iy, Ny-2);
  iz = max(0,iz); iz = min(iz, Nz-2);

  double h = (*Xgrid)(ix+1) - (*Xgrid)(ix);
  double k = (*Ygrid)(iy+1) - (*Ygrid)(iy);
  double l = (*Zgrid)(iz+1) - (*Zgrid)(iz);
  double linv = 1.0/l;

  double u = (x - (*Xgrid)(ix))/h;
  double v = (y - (*Ygrid)(iy))/k;
  double w = (z - (*Zgrid)(iz))/l;

  double a0 = p1(u);
  double a1 = p2(u);
  double a2 = h*q1(u);
  double a3 = h*q2(u);

  register double b0 = p1(v);
  register double b1 = p2(v);
  register double b2 = k*q1(v);
  register double b3 = k*q2(v);

  register double d2c0 = linv*linv*d2p1(w);
  register double d2c1 = linv*linv*d2p2(w);
  register double d2c2 = linv*d2q1(w);
  register double d2c3 = linv*d2q2(w);
  
  for (int i=0; i<N; i++) {
    TinyVector<double,2>& Y000 = F(ix,iy,iz,i)[0];      //   F
    TinyVector<double,2>& Y001 = F(ix,iy,iz+1,i)[0];    //   F
    TinyVector<double,2>& Y002 = F(ix,iy,iz,i)[3];      //  dF/dz
    TinyVector<double,2>& Y003 = F(ix,iy,iz+1,i)[3];    //  dF/dz
    TinyVector<double,2>& Y010 = F(ix,iy+1,iz,i)[0];    //   F
    TinyVector<double,2>& Y011 = F(ix,iy+1,iz+1,i)[0];  //   F
    TinyVector<double,2>& Y012 = F(ix,iy+1,iz,i)[3];    //  dF/dz
    TinyVector<double,2>& Y013 = F(ix,iy+1,iz+1,i)[3];  //  dF/dz
    TinyVector<double,2>& Y020 = F(ix,iy,iz,i)[2];      //  dF/dy
    TinyVector<double,2>& Y021 = F(ix,iy,iz+1,i)[2];    //  dF/dy
    TinyVector<double,2>& Y022 = F(ix,iy,iz,i)[6];      // d2F/dydz
    TinyVector<double,2>& Y023 = F(ix,iy,iz+1,i)[6];    // d2F/dydz
    TinyVector<double,2>& Y030 = F(ix,iy+1,iz,i)[2];    //  dF/dy
    TinyVector<double,2>& Y031 = F(ix,iy+1,iz+1,i)[2];  //  dF/dy
    TinyVector<double,2>& Y032 = F(ix,iy+1,iz,i)[6];    // d2F/dydz
    TinyVector<double,2>& Y033 = F(ix,iy+1,iz+1,i)[6];  // d2F/dydz
    
    TinyVector<double,2>& Y100 = F(ix+1,iy,iz,i)[0];      //   F
    TinyVector<double,2>& Y101 = F(ix+1,iy,iz+1,i)[0];    //   F
    TinyVector<double,2>& Y102 = F(ix+1,iy,iz,i)[3];      //  dF/dz
    TinyVector<double,2>& Y103 = F(ix+1,iy,iz+1,i)[3];    //  dF/dz
    TinyVector<double,2>& Y110 = F(ix+1,iy+1,iz,i)[0];    //   F
    TinyVector<double,2>& Y111 = F(ix+1,iy+1,iz+1,i)[0];  //   F
    TinyVector<double,2>& Y112 = F(ix+1,iy+1,iz,i)[3];    //  dF/dz
    TinyVector<double,2>& Y113 = F(ix+1,iy+1,iz+1,i)[3];  //  dF/dz
    TinyVector<double,2>& Y120 = F(ix+1,iy,iz,i)[2];      //  dF/dy
    TinyVector<double,2>& Y121 = F(ix+1,iy,iz+1,i)[2];    //  dF/dy
    TinyVector<double,2>& Y122 = F(ix+1,iy,iz,i)[6];      // d2F/dydz
    TinyVector<double,2>& Y123 = F(ix+1,iy,iz+1,i)[6];    // d2F/dydz
    TinyVector<double,2>& Y130 = F(ix+1,iy+1,iz,i)[2];    //  dF/dy
    TinyVector<double,2>& Y131 = F(ix+1,iy+1,iz+1,i)[2];  //  dF/dy
    TinyVector<double,2>& Y132 = F(ix+1,iy+1,iz,i)[6];    // d2F/dydz
    TinyVector<double,2>& Y133 = F(ix+1,iy+1,iz+1,i)[6];  // d2F/dydz
    
    TinyVector<double,2>& Y200 = F(ix,iy,iz,i)[1];      //  dF/dx
    TinyVector<double,2>& Y201 = F(ix,iy,iz+1,i)[1];    //  dF/dx
    TinyVector<double,2>& Y202 = F(ix,iy,iz,i)[5];      // d2F/dxdz
    TinyVector<double,2>& Y203 = F(ix,iy,iz+1,i)[5];    // d2F/dxdz
    TinyVector<double,2>& Y210 = F(ix,iy+1,iz,i)[1];    //  dF/dx
    TinyVector<double,2>& Y211 = F(ix,iy+1,iz+1,i)[1];  //  dF/dx
    TinyVector<double,2>& Y212 = F(ix,iy+1,iz,i)[5];    // d2F/dxdz
    TinyVector<double,2>& Y213 = F(ix,iy+1,iz+1,i)[5];  // d2F/dxdz
    TinyVector<double,2>& Y220 = F(ix,iy,iz,i)[4];      // d2F/dxdy
    TinyVector<double,2>& Y221 = F(ix,iy,iz+1,i)[4];    // d2F/dxdy
    TinyVector<double,2>& Y222 = F(ix,iy,iz,i)[7];      // d3F/dxdydz
    TinyVector<double,2>& Y223 = F(ix,iy,iz+1,i)[7];    // d3F/dxdydz
    TinyVector<double,2>& Y230 = F(ix,iy+1,iz,i)[4];    // d2F/dxdy
    TinyVector<double,2>& Y231 = F(ix,iy+1,iz+1,i)[4];  // d2F/dxdy
    TinyVector<double,2>& Y232 = F(ix,iy+1,iz,i)[7];    // d3F/dxdydz
    TinyVector<double,2>& Y233 = F(ix,iy+1,iz+1,i)[7];  // d3F/dxdydz
    
    TinyVector<double,2>& Y300 = F(ix+1,iy,iz,i)[1];      //  dF/dx
    TinyVector<double,2>& Y301 = F(ix+1,iy,iz+1,i)[1];    //  dF/dx
    TinyVector<double,2>& Y302 = F(ix+1,iy,iz,i)[5];      // d2F/dxdz
    TinyVector<double,2>& Y303 = F(ix+1,iy,iz+1,i)[5];    // d2F/dxdz
    TinyVector<double,2>& Y310 = F(ix+1,iy+1,iz,i)[1];    //  dF/dx
    TinyVector<double,2>& Y311 = F(ix+1,iy+1,iz+1,i)[1];  //  dF/dx
    TinyVector<double,2>& Y312 = F(ix+1,iy+1,iz,i)[5];    // d2F/dxdz
    TinyVector<double,2>& Y313 = F(ix+1,iy+1,iz+1,i)[5];  // d2F/dxdz
    TinyVector<double,2>& Y320 = F(ix+1,iy,iz,i)[4];      // d2F/dxdy
    TinyVector<double,2>& Y321 = F(ix+1,iy,iz+1,i)[4];    // d2F/dxdy
    TinyVector<double,2>& Y322 = F(ix+1,iy,iz,i)[7];      // d3F/dxdydz
    TinyVector<double,2>& Y323 = F(ix+1,iy,iz+1,i)[7];    // d3F/dxdydz
    TinyVector<double,2>& Y330 = F(ix+1,iy+1,iz,i)[4];    // d2F/dxdy
    TinyVector<double,2>& Y331 = F(ix+1,iy+1,iz+1,i)[4];  // d2F/dxdy
    TinyVector<double,2>& Y332 = F(ix+1,iy+1,iz,i)[7];    // d3F/dxdydz
    TinyVector<double,2>& Y333 = F(ix+1,iy+1,iz+1,i)[7];  // d3F/dxdydz
  
    TinyVector<double,2> Val = 
      a0*
      (b0*(Y000*d2c0+Y001*d2c1+Y002*d2c2+Y003*d2c3) +
       b1*(Y010*d2c0+Y011*d2c1+Y012*d2c2+Y013*d2c3) +
       b2*(Y020*d2c0+Y021*d2c1+Y022*d2c2+Y023*d2c3) +
       b3*(Y030*d2c0+Y031*d2c1+Y032*d2c2+Y033*d2c3))+
      a1 *
      (b0*(Y100*d2c0+Y101*d2c1+Y102*d2c2+Y103*d2c3) +
       b1*(Y110*d2c0+Y111*d2c1+Y112*d2c2+Y113*d2c3) +
       b2*(Y120*d2c0+Y121*d2c1+Y122*d2c2+Y123*d2c3) +
       b3*(Y130*d2c0+Y131*d2c1+Y132*d2c2+Y133*d2c3))+
      a2 *
      (b0*(Y200*d2c0+Y201*d2c1+Y202*d2c2+Y203*d2c3) +
       b1*(Y210*d2c0+Y211*d2c1+Y212*d2c2+Y213*d2c3) +
       b2*(Y220*d2c0+Y221*d2c1+Y222*d2c2+Y223*d2c3) +
       b3*(Y230*d2c0+Y231*d2c1+Y232*d2c2+Y233*d2c3))+
      a3 *
      (b0*(Y300*d2c0+Y301*d2c1+Y302*d2c2+Y303*d2c3) +
       b1*(Y310*d2c0+Y311*d2c1+Y312*d2c2+Y313*d2c3) +
       b2*(Y320*d2c0+Y321*d2c1+Y322*d2c2+Y323*d2c3) +
       b3*(Y330*d2c0+Y331*d2c1+Y332*d2c2+Y333*d2c3));

    vals(i) = complex<double>(Val[0], Val[1]);
  }
}


inline void 
ComplexMultiTricubicSpline::Laplacian (double x, double y, double z, 
				       Array<complex<double>,1> &vals)
{
  if (!UpToDate)
    Update();

  int ix = Xgrid->ReverseMap(x);  
  int iy = Ygrid->ReverseMap(y);
  int iz = Zgrid->ReverseMap(z);

  ix = max(0,ix); ix = min(ix, Nx-2);
  iy = max(0,iy); iy = min(iy, Ny-2);
  iz = max(0,iz); iz = min(iz, Nz-2);

  double h = (*Xgrid)(ix+1) - (*Xgrid)(ix);
  double hinv = 1.0/h;
  double k = (*Ygrid)(iy+1) - (*Ygrid)(iy);
  double kinv = 1.0/k;
  double l = (*Zgrid)(iz+1) - (*Zgrid)(iz);
  double linv = 1.0/l;

  double u = (x - (*Xgrid)(ix))/h;
  double v = (y - (*Ygrid)(iy))/k;
  double w = (z - (*Zgrid)(iz))/l;

  double a0 = p1(u);
  double a1 = p2(u);
  double a2 = h*q1(u);
  double a3 = h*q2(u);
  double d2a0 = hinv*hinv*d2p1(u);
  double d2a1 = hinv*hinv*d2p2(u);
  double d2a2 = hinv*d2q1(u);
  double d2a3 = hinv*d2q2(u);

  register double b0 = p1(v);
  register double b1 = p2(v);
  register double b2 = k*q1(v);
  register double b3 = k*q2(v);
  register double d2b0 = kinv*kinv*d2p1(v);
  register double d2b1 = kinv*kinv*d2p2(v);
  register double d2b2 = kinv*d2q1(v);
  register double d2b3 = kinv*d2q2(v);

  register double c0 = p1(w);
  register double c1 = p2(w);
  register double c2 = l*q1(w);
  register double c3 = l*q2(w);
  register double d2c0 = linv*linv*d2p1(w);
  register double d2c1 = linv*linv*d2p2(w);
  register double d2c2 = linv*d2q1(w);
  register double d2c3 = linv*d2q2(w);
  
  for (int i=0; i<N; i++) {
    TinyVector<double,2>& Y000 = F(ix,iy,iz,i)[0];      //   F
    TinyVector<double,2>& Y001 = F(ix,iy,iz+1,i)[0];    //   F
    TinyVector<double,2>& Y002 = F(ix,iy,iz,i)[3];      //  dF/dz
    TinyVector<double,2>& Y003 = F(ix,iy,iz+1,i)[3];    //  dF/dz
    TinyVector<double,2>& Y010 = F(ix,iy+1,iz,i)[0];    //   F
    TinyVector<double,2>& Y011 = F(ix,iy+1,iz+1,i)[0];  //   F
    TinyVector<double,2>& Y012 = F(ix,iy+1,iz,i)[3];    //  dF/dz
    TinyVector<double,2>& Y013 = F(ix,iy+1,iz+1,i)[3];  //  dF/dz
    TinyVector<double,2>& Y020 = F(ix,iy,iz,i)[2];      //  dF/dy
    TinyVector<double,2>& Y021 = F(ix,iy,iz+1,i)[2];    //  dF/dy
    TinyVector<double,2>& Y022 = F(ix,iy,iz,i)[6];      // d2F/dydz
    TinyVector<double,2>& Y023 = F(ix,iy,iz+1,i)[6];    // d2F/dydz
    TinyVector<double,2>& Y030 = F(ix,iy+1,iz,i)[2];    //  dF/dy
    TinyVector<double,2>& Y031 = F(ix,iy+1,iz+1,i)[2];  //  dF/dy
    TinyVector<double,2>& Y032 = F(ix,iy+1,iz,i)[6];    // d2F/dydz
    TinyVector<double,2>& Y033 = F(ix,iy+1,iz+1,i)[6];  // d2F/dydz
    
    TinyVector<double,2>& Y100 = F(ix+1,iy,iz,i)[0];      //   F
    TinyVector<double,2>& Y101 = F(ix+1,iy,iz+1,i)[0];    //   F
    TinyVector<double,2>& Y102 = F(ix+1,iy,iz,i)[3];      //  dF/dz
    TinyVector<double,2>& Y103 = F(ix+1,iy,iz+1,i)[3];    //  dF/dz
    TinyVector<double,2>& Y110 = F(ix+1,iy+1,iz,i)[0];    //   F
    TinyVector<double,2>& Y111 = F(ix+1,iy+1,iz+1,i)[0];  //   F
    TinyVector<double,2>& Y112 = F(ix+1,iy+1,iz,i)[3];    //  dF/dz
    TinyVector<double,2>& Y113 = F(ix+1,iy+1,iz+1,i)[3];  //  dF/dz
    TinyVector<double,2>& Y120 = F(ix+1,iy,iz,i)[2];      //  dF/dy
    TinyVector<double,2>& Y121 = F(ix+1,iy,iz+1,i)[2];    //  dF/dy
    TinyVector<double,2>& Y122 = F(ix+1,iy,iz,i)[6];      // d2F/dydz
    TinyVector<double,2>& Y123 = F(ix+1,iy,iz+1,i)[6];    // d2F/dydz
    TinyVector<double,2>& Y130 = F(ix+1,iy+1,iz,i)[2];    //  dF/dy
    TinyVector<double,2>& Y131 = F(ix+1,iy+1,iz+1,i)[2];  //  dF/dy
    TinyVector<double,2>& Y132 = F(ix+1,iy+1,iz,i)[6];    // d2F/dydz
    TinyVector<double,2>& Y133 = F(ix+1,iy+1,iz+1,i)[6];  // d2F/dydz
    
    TinyVector<double,2>& Y200 = F(ix,iy,iz,i)[1];      //  dF/dx
    TinyVector<double,2>& Y201 = F(ix,iy,iz+1,i)[1];    //  dF/dx
    TinyVector<double,2>& Y202 = F(ix,iy,iz,i)[5];      // d2F/dxdz
    TinyVector<double,2>& Y203 = F(ix,iy,iz+1,i)[5];    // d2F/dxdz
    TinyVector<double,2>& Y210 = F(ix,iy+1,iz,i)[1];    //  dF/dx
    TinyVector<double,2>& Y211 = F(ix,iy+1,iz+1,i)[1];  //  dF/dx
    TinyVector<double,2>& Y212 = F(ix,iy+1,iz,i)[5];    // d2F/dxdz
    TinyVector<double,2>& Y213 = F(ix,iy+1,iz+1,i)[5];  // d2F/dxdz
    TinyVector<double,2>& Y220 = F(ix,iy,iz,i)[4];      // d2F/dxdy
    TinyVector<double,2>& Y221 = F(ix,iy,iz+1,i)[4];    // d2F/dxdy
    TinyVector<double,2>& Y222 = F(ix,iy,iz,i)[7];      // d3F/dxdydz
    TinyVector<double,2>& Y223 = F(ix,iy,iz+1,i)[7];    // d3F/dxdydz
    TinyVector<double,2>& Y230 = F(ix,iy+1,iz,i)[4];    // d2F/dxdy
    TinyVector<double,2>& Y231 = F(ix,iy+1,iz+1,i)[4];  // d2F/dxdy
    TinyVector<double,2>& Y232 = F(ix,iy+1,iz,i)[7];    // d3F/dxdydz
    TinyVector<double,2>& Y233 = F(ix,iy+1,iz+1,i)[7];  // d3F/dxdydz
    
    TinyVector<double,2>& Y300 = F(ix+1,iy,iz,i)[1];      //  dF/dx
    TinyVector<double,2>& Y301 = F(ix+1,iy,iz+1,i)[1];    //  dF/dx
    TinyVector<double,2>& Y302 = F(ix+1,iy,iz,i)[5];      // d2F/dxdz
    TinyVector<double,2>& Y303 = F(ix+1,iy,iz+1,i)[5];    // d2F/dxdz
    TinyVector<double,2>& Y310 = F(ix+1,iy+1,iz,i)[1];    //  dF/dx
    TinyVector<double,2>& Y311 = F(ix+1,iy+1,iz+1,i)[1];  //  dF/dx
    TinyVector<double,2>& Y312 = F(ix+1,iy+1,iz,i)[5];    // d2F/dxdz
    TinyVector<double,2>& Y313 = F(ix+1,iy+1,iz+1,i)[5];  // d2F/dxdz
    TinyVector<double,2>& Y320 = F(ix+1,iy,iz,i)[4];      // d2F/dxdy
    TinyVector<double,2>& Y321 = F(ix+1,iy,iz+1,i)[4];    // d2F/dxdy
    TinyVector<double,2>& Y322 = F(ix+1,iy,iz,i)[7];      // d3F/dxdydz
    TinyVector<double,2>& Y323 = F(ix+1,iy,iz+1,i)[7];    // d3F/dxdydz
    TinyVector<double,2>& Y330 = F(ix+1,iy+1,iz,i)[4];    // d2F/dxdy
    TinyVector<double,2>& Y331 = F(ix+1,iy+1,iz+1,i)[4];  // d2F/dxdy
    TinyVector<double,2>& Y332 = F(ix+1,iy+1,iz,i)[7];    // d3F/dxdydz
    TinyVector<double,2>& Y333 = F(ix+1,iy+1,iz+1,i)[7];  // d3F/dxdydz

    TinyVector<double,2> Val =
      d2a0*
      (b0*(Y000*c0+Y001*c1+Y002*c2+Y003*c3) +
       b1*(Y010*c0+Y011*c1+Y012*c2+Y013*c3) +
       b2*(Y020*c0+Y021*c1+Y022*c2+Y023*c3) +
       b3*(Y030*c0+Y031*c1+Y032*c2+Y033*c3))+
      d2a1 *
      (b0*(Y100*c0+Y101*c1+Y102*c2+Y103*c3) +
       b1*(Y110*c0+Y111*c1+Y112*c2+Y113*c3) +
       b2*(Y120*c0+Y121*c1+Y122*c2+Y123*c3) +
       b3*(Y130*c0+Y131*c1+Y132*c2+Y133*c3))+
      d2a2 *
      (b0*(Y200*c0+Y201*c1+Y202*c2+Y203*c3) +
       b1*(Y210*c0+Y211*c1+Y212*c2+Y213*c3) +
       b2*(Y220*c0+Y221*c1+Y222*c2+Y223*c3) +
       b3*(Y230*c0+Y231*c1+Y232*c2+Y233*c3))+
      d2a3 *
      (b0*(Y300*c0+Y301*c1+Y302*c2+Y303*c3) +
       b1*(Y310*c0+Y311*c1+Y312*c2+Y313*c3) +
       b2*(Y320*c0+Y321*c1+Y322*c2+Y323*c3) +
       b3*(Y330*c0+Y331*c1+Y332*c2+Y333*c3));
  
    Val +=
      a0*
      (d2b0*(Y000*c0+Y001*c1+Y002*c2+Y003*c3) +
       d2b1*(Y010*c0+Y011*c1+Y012*c2+Y013*c3) +
       d2b2*(Y020*c0+Y021*c1+Y022*c2+Y023*c3) +
       d2b3*(Y030*c0+Y031*c1+Y032*c2+Y033*c3))+
      a1 *
      (d2b0*(Y100*c0+Y101*c1+Y102*c2+Y103*c3) +
       d2b1*(Y110*c0+Y111*c1+Y112*c2+Y113*c3) +
       d2b2*(Y120*c0+Y121*c1+Y122*c2+Y123*c3) +
       d2b3*(Y130*c0+Y131*c1+Y132*c2+Y133*c3))+
      a2 *
      (d2b0*(Y200*c0+Y201*c1+Y202*c2+Y203*c3) +
       d2b1*(Y210*c0+Y211*c1+Y212*c2+Y213*c3) +
       d2b2*(Y220*c0+Y221*c1+Y222*c2+Y223*c3) +
       d2b3*(Y230*c0+Y231*c1+Y232*c2+Y233*c3))+
      a3 *
      (d2b0*(Y300*c0+Y301*c1+Y302*c2+Y303*c3) +
       d2b1*(Y310*c0+Y311*c1+Y312*c2+Y313*c3) +
       d2b2*(Y320*c0+Y321*c1+Y322*c2+Y323*c3) +
       d2b3*(Y330*c0+Y331*c1+Y332*c2+Y333*c3));

    Val += 
      a0*
      (b0*(Y000*d2c0+Y001*d2c1+Y002*d2c2+Y003*d2c3) +
       b1*(Y010*d2c0+Y011*d2c1+Y012*d2c2+Y013*d2c3) +
       b2*(Y020*d2c0+Y021*d2c1+Y022*d2c2+Y023*d2c3) +
       b3*(Y030*d2c0+Y031*d2c1+Y032*d2c2+Y033*d2c3))+
      a1 *
      (b0*(Y100*d2c0+Y101*d2c1+Y102*d2c2+Y103*d2c3) +
       b1*(Y110*d2c0+Y111*d2c1+Y112*d2c2+Y113*d2c3) +
       b2*(Y120*d2c0+Y121*d2c1+Y122*d2c2+Y123*d2c3) +
       b3*(Y130*d2c0+Y131*d2c1+Y132*d2c2+Y133*d2c3))+
      a2 *
      (b0*(Y200*d2c0+Y201*d2c1+Y202*d2c2+Y203*d2c3) +
       b1*(Y210*d2c0+Y211*d2c1+Y212*d2c2+Y213*d2c3) +
       b2*(Y220*d2c0+Y221*d2c1+Y222*d2c2+Y223*d2c3) +
       b3*(Y230*d2c0+Y231*d2c1+Y232*d2c2+Y233*d2c3))+
      a3 *
      (b0*(Y300*d2c0+Y301*d2c1+Y302*d2c2+Y303*d2c3) +
       b1*(Y310*d2c0+Y311*d2c1+Y312*d2c2+Y313*d2c3) +
       b2*(Y320*d2c0+Y321*d2c1+Y322*d2c2+Y323*d2c3) +
       b3*(Y330*d2c0+Y331*d2c1+Y332*d2c2+Y333*d2c3));

    vals(i) = complex<double>(Val[0],Val[1]);
  }
}


inline void
ComplexMultiTricubicSpline::IntegralOfLaplacian(Array<complex<double>,1> &vals)
{
  vals = 0.0;
  for (int ix=0; ix<(Nx-1); ix++) {
    for (int iy=0; iy<(Ny-1); iy++) {
      for (int iz=0; iz<(Nz-1); iz++) {
	double h = (*Xgrid)(ix+1) - (*Xgrid)(ix);
	double hinv = 1.0/h;
	double k = (*Ygrid)(iy+1) - (*Ygrid)(iy);
	double kinv = 1.0/k;
	double l = (*Zgrid)(iz+1) - (*Zgrid)(iz);
	double linv = 1.0/l;
		
	double a0 = h*Int_p1;
	double a1 = h*Int_p2;
	double a2 = h*h*Int_q1;
	double a3 = h*h*Int_q2;
	double d2a0 = hinv*Int_d2p1;
	double d2a1 = hinv*Int_d2p2;
	double d2a2 = Int_d2q1;
	double d2a3 = Int_d2q2;
	
	register double b0 = k*Int_p1;
	register double b1 = k*Int_p2;
	register double b2 = k*k*Int_q1;
	register double b3 = k*k*Int_q2;
	register double d2b0 = kinv*Int_d2p1;
	register double d2b1 = kinv*Int_d2p2;
	register double d2b2 = Int_d2q1;
	register double d2b3 = Int_d2q2;
	
	register double c0 = l*Int_p1;
	register double c1 = l*Int_p2;
	register double c2 = l*l*Int_q1;
	register double c3 = l*l*Int_q2;
	register double d2c0 = linv*Int_d2p1;
	register double d2c1 = linv*Int_d2p2;
	register double d2c2 = Int_d2q1;
	register double d2c3 = Int_d2q2;
	
	for (int i=0; i<N; i++) {
	  TinyVector<double,2>& Y000 = F(ix,iy,iz,i)[0];      //   F
	  TinyVector<double,2>& Y001 = F(ix,iy,iz+1,i)[0];    //   F
	  TinyVector<double,2>& Y002 = F(ix,iy,iz,i)[3];      //  dF/dz
	  TinyVector<double,2>& Y003 = F(ix,iy,iz+1,i)[3];    //  dF/dz
	  TinyVector<double,2>& Y010 = F(ix,iy+1,iz,i)[0];    //   F
	  TinyVector<double,2>& Y011 = F(ix,iy+1,iz+1,i)[0];  //   F
	  TinyVector<double,2>& Y012 = F(ix,iy+1,iz,i)[3];    //  dF/dz
	  TinyVector<double,2>& Y013 = F(ix,iy+1,iz+1,i)[3];  //  dF/dz
	  TinyVector<double,2>& Y020 = F(ix,iy,iz,i)[2];      //  dF/dy
	  TinyVector<double,2>& Y021 = F(ix,iy,iz+1,i)[2];    //  dF/dy
	  TinyVector<double,2>& Y022 = F(ix,iy,iz,i)[6];      // d2F/dydz
	  TinyVector<double,2>& Y023 = F(ix,iy,iz+1,i)[6];    // d2F/dydz
	  TinyVector<double,2>& Y030 = F(ix,iy+1,iz,i)[2];    //  dF/dy
	  TinyVector<double,2>& Y031 = F(ix,iy+1,iz+1,i)[2];  //  dF/dy
	  TinyVector<double,2>& Y032 = F(ix,iy+1,iz,i)[6];    // d2F/dydz
	  TinyVector<double,2>& Y033 = F(ix,iy+1,iz+1,i)[6];  // d2F/dydz
	  
	  TinyVector<double,2>& Y100 = F(ix+1,iy,iz,i)[0];      //   F
	  TinyVector<double,2>& Y101 = F(ix+1,iy,iz+1,i)[0];    //   F
	  TinyVector<double,2>& Y102 = F(ix+1,iy,iz,i)[3];      //  dF/dz
	  TinyVector<double,2>& Y103 = F(ix+1,iy,iz+1,i)[3];    //  dF/dz
	  TinyVector<double,2>& Y110 = F(ix+1,iy+1,iz,i)[0];    //   F
	  TinyVector<double,2>& Y111 = F(ix+1,iy+1,iz+1,i)[0];  //   F
	  TinyVector<double,2>& Y112 = F(ix+1,iy+1,iz,i)[3];    //  dF/dz
	  TinyVector<double,2>& Y113 = F(ix+1,iy+1,iz+1,i)[3];  //  dF/dz
	  TinyVector<double,2>& Y120 = F(ix+1,iy,iz,i)[2];      //  dF/dy
	  TinyVector<double,2>& Y121 = F(ix+1,iy,iz+1,i)[2];    //  dF/dy
	  TinyVector<double,2>& Y122 = F(ix+1,iy,iz,i)[6];      // d2F/dydz
	  TinyVector<double,2>& Y123 = F(ix+1,iy,iz+1,i)[6];    // d2F/dydz
	  TinyVector<double,2>& Y130 = F(ix+1,iy+1,iz,i)[2];    //  dF/dy
	  TinyVector<double,2>& Y131 = F(ix+1,iy+1,iz+1,i)[2];  //  dF/dy
	  TinyVector<double,2>& Y132 = F(ix+1,iy+1,iz,i)[6];    // d2F/dydz
	  TinyVector<double,2>& Y133 = F(ix+1,iy+1,iz+1,i)[6];  // d2F/dydz
	  
	  TinyVector<double,2>& Y200 = F(ix,iy,iz,i)[1];      //  dF/dx
	  TinyVector<double,2>& Y201 = F(ix,iy,iz+1,i)[1];    //  dF/dx
	  TinyVector<double,2>& Y202 = F(ix,iy,iz,i)[5];      // d2F/dxdz
	  TinyVector<double,2>& Y203 = F(ix,iy,iz+1,i)[5];    // d2F/dxdz
	  TinyVector<double,2>& Y210 = F(ix,iy+1,iz,i)[1];    //  dF/dx
	  TinyVector<double,2>& Y211 = F(ix,iy+1,iz+1,i)[1];  //  dF/dx
	  TinyVector<double,2>& Y212 = F(ix,iy+1,iz,i)[5];    // d2F/dxdz
	  TinyVector<double,2>& Y213 = F(ix,iy+1,iz+1,i)[5];  // d2F/dxdz
	  TinyVector<double,2>& Y220 = F(ix,iy,iz,i)[4];      // d2F/dxdy
	  TinyVector<double,2>& Y221 = F(ix,iy,iz+1,i)[4];    // d2F/dxdy
	  TinyVector<double,2>& Y222 = F(ix,iy,iz,i)[7];      // d3F/dxdydz
	  TinyVector<double,2>& Y223 = F(ix,iy,iz+1,i)[7];    // d3F/dxdydz
	  TinyVector<double,2>& Y230 = F(ix,iy+1,iz,i)[4];    // d2F/dxdy
	  TinyVector<double,2>& Y231 = F(ix,iy+1,iz+1,i)[4];  // d2F/dxdy
	  TinyVector<double,2>& Y232 = F(ix,iy+1,iz,i)[7];    // d3F/dxdydz
	  TinyVector<double,2>& Y233 = F(ix,iy+1,iz+1,i)[7];  // d3F/dxdydz
	  
	  TinyVector<double,2>& Y300 = F(ix+1,iy,iz,i)[1];      //  dF/dx
	  TinyVector<double,2>& Y301 = F(ix+1,iy,iz+1,i)[1];    //  dF/dx
	  TinyVector<double,2>& Y302 = F(ix+1,iy,iz,i)[5];      // d2F/dxdz
	  TinyVector<double,2>& Y303 = F(ix+1,iy,iz+1,i)[5];    // d2F/dxdz
	  TinyVector<double,2>& Y310 = F(ix+1,iy+1,iz,i)[1];    //  dF/dx
	  TinyVector<double,2>& Y311 = F(ix+1,iy+1,iz+1,i)[1];  //  dF/dx
	  TinyVector<double,2>& Y312 = F(ix+1,iy+1,iz,i)[5];    // d2F/dxdz
	  TinyVector<double,2>& Y313 = F(ix+1,iy+1,iz+1,i)[5];  // d2F/dxdz
	  TinyVector<double,2>& Y320 = F(ix+1,iy,iz,i)[4];      // d2F/dxdy
	  TinyVector<double,2>& Y321 = F(ix+1,iy,iz+1,i)[4];    // d2F/dxdy
	  TinyVector<double,2>& Y322 = F(ix+1,iy,iz,i)[7];      // d3F/dxdydz
	  TinyVector<double,2>& Y323 = F(ix+1,iy,iz+1,i)[7];    // d3F/dxdydz
	  TinyVector<double,2>& Y330 = F(ix+1,iy+1,iz,i)[4];    // d2F/dxdy
	  TinyVector<double,2>& Y331 = F(ix+1,iy+1,iz+1,i)[4];  // d2F/dxdy
	  TinyVector<double,2>& Y332 = F(ix+1,iy+1,iz,i)[7];    // d3F/dxdydz
	  TinyVector<double,2>& Y333 = F(ix+1,iy+1,iz+1,i)[7];  // d3F/dxdydz
	  
	  TinyVector<double,2> Val =
	    d2a0*
	    (b0*(Y000*c0+Y001*c1+Y002*c2+Y003*c3) +
	     b1*(Y010*c0+Y011*c1+Y012*c2+Y013*c3) +
	     b2*(Y020*c0+Y021*c1+Y022*c2+Y023*c3) +
	     b3*(Y030*c0+Y031*c1+Y032*c2+Y033*c3))+
	    d2a1 *
	    (b0*(Y100*c0+Y101*c1+Y102*c2+Y103*c3) +
	     b1*(Y110*c0+Y111*c1+Y112*c2+Y113*c3) +
	     b2*(Y120*c0+Y121*c1+Y122*c2+Y123*c3) +
	     b3*(Y130*c0+Y131*c1+Y132*c2+Y133*c3))+
	    d2a2 *
	    (b0*(Y200*c0+Y201*c1+Y202*c2+Y203*c3) +
	     b1*(Y210*c0+Y211*c1+Y212*c2+Y213*c3) +
	     b2*(Y220*c0+Y221*c1+Y222*c2+Y223*c3) +
	     b3*(Y230*c0+Y231*c1+Y232*c2+Y233*c3))+
	    d2a3 *
	    (b0*(Y300*c0+Y301*c1+Y302*c2+Y303*c3) +
	     b1*(Y310*c0+Y311*c1+Y312*c2+Y313*c3) +
	     b2*(Y320*c0+Y321*c1+Y322*c2+Y323*c3) +
	     b3*(Y330*c0+Y331*c1+Y332*c2+Y333*c3));
	  
	  Val +=
	    a0*
	    (d2b0*(Y000*c0+Y001*c1+Y002*c2+Y003*c3) +
	     d2b1*(Y010*c0+Y011*c1+Y012*c2+Y013*c3) +
	     d2b2*(Y020*c0+Y021*c1+Y022*c2+Y023*c3) +
	     d2b3*(Y030*c0+Y031*c1+Y032*c2+Y033*c3))+
	    a1 *
	    (d2b0*(Y100*c0+Y101*c1+Y102*c2+Y103*c3) +
	     d2b1*(Y110*c0+Y111*c1+Y112*c2+Y113*c3) +
	     d2b2*(Y120*c0+Y121*c1+Y122*c2+Y123*c3) +
	     d2b3*(Y130*c0+Y131*c1+Y132*c2+Y133*c3))+
	    a2 *
	    (d2b0*(Y200*c0+Y201*c1+Y202*c2+Y203*c3) +
	     d2b1*(Y210*c0+Y211*c1+Y212*c2+Y213*c3) +
	     d2b2*(Y220*c0+Y221*c1+Y222*c2+Y223*c3) +
	     d2b3*(Y230*c0+Y231*c1+Y232*c2+Y233*c3))+
	    a3 *
	    (d2b0*(Y300*c0+Y301*c1+Y302*c2+Y303*c3) +
	     d2b1*(Y310*c0+Y311*c1+Y312*c2+Y313*c3) +
	     d2b2*(Y320*c0+Y321*c1+Y322*c2+Y323*c3) +
	     d2b3*(Y330*c0+Y331*c1+Y332*c2+Y333*c3));
	  
	  Val += 
	    a0*
	    (b0*(Y000*d2c0+Y001*d2c1+Y002*d2c2+Y003*d2c3) +
	     b1*(Y010*d2c0+Y011*d2c1+Y012*d2c2+Y013*d2c3) +
	     b2*(Y020*d2c0+Y021*d2c1+Y022*d2c2+Y023*d2c3) +
	     b3*(Y030*d2c0+Y031*d2c1+Y032*d2c2+Y033*d2c3))+
	    a1 *
	    (b0*(Y100*d2c0+Y101*d2c1+Y102*d2c2+Y103*d2c3) +
	     b1*(Y110*d2c0+Y111*d2c1+Y112*d2c2+Y113*d2c3) +
	     b2*(Y120*d2c0+Y121*d2c1+Y122*d2c2+Y123*d2c3) +
	     b3*(Y130*d2c0+Y131*d2c1+Y132*d2c2+Y133*d2c3))+
	    a2 *
	    (b0*(Y200*d2c0+Y201*d2c1+Y202*d2c2+Y203*d2c3) +
	     b1*(Y210*d2c0+Y211*d2c1+Y212*d2c2+Y213*d2c3) +
	     b2*(Y220*d2c0+Y221*d2c1+Y222*d2c2+Y223*d2c3) +
	     b3*(Y230*d2c0+Y231*d2c1+Y232*d2c2+Y233*d2c3))+
	    a3 *
	    (b0*(Y300*d2c0+Y301*d2c1+Y302*d2c2+Y303*d2c3) +
	     b1*(Y310*d2c0+Y311*d2c1+Y312*d2c2+Y313*d2c3) +
	     b2*(Y320*d2c0+Y321*d2c1+Y322*d2c2+Y323*d2c3) +
	     b3*(Y330*d2c0+Y331*d2c1+Y332*d2c2+Y333*d2c3));
	  
	  vals(i) += complex<double>(Val[0],Val[1]);
	}
      }
    }
  }
}

inline void 
ComplexMultiTricubicSpline::d2_dxdy (double x, double y, double z, 
				     Array<complex<double>,1> &vals)
{
  if (!UpToDate)
    Update();

  int ix = Xgrid->ReverseMap(x);  
  int iy = Ygrid->ReverseMap(y);
  int iz = Zgrid->ReverseMap(z);

  ix = max(0,ix); ix = min(ix, Nx-2);
  iy = max(0,iy); iy = min(iy, Ny-2);
  iz = max(0,iz); iz = min(iz, Nz-2);

  double h = (*Xgrid)(ix+1) - (*Xgrid)(ix);
  double hinv = 1.0/h;
  double k = (*Ygrid)(iy+1) - (*Ygrid)(iy);
  double kinv = 1.0/k;
  double l = (*Zgrid)(iz+1) - (*Zgrid)(iz);
  double u = (x - (*Xgrid)(ix))/h;
  double v = (y - (*Ygrid)(iy))/k;
  double w = (z - (*Zgrid)(iz))/l;

  double da0 = hinv*dp1(u);
  double da1 = hinv*dp2(u);
  double da2 = dq1(u);
  double da3 = dq2(u);

  register double db0 = kinv*dp1(v);
  register double db1 = kinv*dp2(v);
  register double db2 = dq1(v);
  register double db3 = dq2(v);

  register double c0 = p1(w);
  register double c1 = p2(w);
  register double c2 = l*q1(w);
  register double c3 = l*q2(w);
  
  for (int i=0; i<N; i++) {
    TinyVector<double,2>& Y000 = F(ix,iy,iz,i)[0];      //   F
    TinyVector<double,2>& Y001 = F(ix,iy,iz+1,i)[0];    //   F
    TinyVector<double,2>& Y002 = F(ix,iy,iz,i)[3];      //  dF/dz
    TinyVector<double,2>& Y003 = F(ix,iy,iz+1,i)[3];    //  dF/dz
    TinyVector<double,2>& Y010 = F(ix,iy+1,iz,i)[0];    //   F
    TinyVector<double,2>& Y011 = F(ix,iy+1,iz+1,i)[0];  //   F
    TinyVector<double,2>& Y012 = F(ix,iy+1,iz,i)[3];    //  dF/dz
    TinyVector<double,2>& Y013 = F(ix,iy+1,iz+1,i)[3];  //  dF/dz
    TinyVector<double,2>& Y020 = F(ix,iy,iz,i)[2];      //  dF/dy
    TinyVector<double,2>& Y021 = F(ix,iy,iz+1,i)[2];    //  dF/dy
    TinyVector<double,2>& Y022 = F(ix,iy,iz,i)[6];      // d2F/dydz
    TinyVector<double,2>& Y023 = F(ix,iy,iz+1,i)[6];    // d2F/dydz
    TinyVector<double,2>& Y030 = F(ix,iy+1,iz,i)[2];    //  dF/dy
    TinyVector<double,2>& Y031 = F(ix,iy+1,iz+1,i)[2];  //  dF/dy
    TinyVector<double,2>& Y032 = F(ix,iy+1,iz,i)[6];    // d2F/dydz
    TinyVector<double,2>& Y033 = F(ix,iy+1,iz+1,i)[6];  // d2F/dydz
    
    TinyVector<double,2>& Y100 = F(ix+1,iy,iz,i)[0];      //   F
    TinyVector<double,2>& Y101 = F(ix+1,iy,iz+1,i)[0];    //   F
    TinyVector<double,2>& Y102 = F(ix+1,iy,iz,i)[3];      //  dF/dz
    TinyVector<double,2>& Y103 = F(ix+1,iy,iz+1,i)[3];    //  dF/dz
    TinyVector<double,2>& Y110 = F(ix+1,iy+1,iz,i)[0];    //   F
    TinyVector<double,2>& Y111 = F(ix+1,iy+1,iz+1,i)[0];  //   F
    TinyVector<double,2>& Y112 = F(ix+1,iy+1,iz,i)[3];    //  dF/dz
    TinyVector<double,2>& Y113 = F(ix+1,iy+1,iz+1,i)[3];  //  dF/dz
    TinyVector<double,2>& Y120 = F(ix+1,iy,iz,i)[2];      //  dF/dy
    TinyVector<double,2>& Y121 = F(ix+1,iy,iz+1,i)[2];    //  dF/dy
    TinyVector<double,2>& Y122 = F(ix+1,iy,iz,i)[6];      // d2F/dydz
    TinyVector<double,2>& Y123 = F(ix+1,iy,iz+1,i)[6];    // d2F/dydz
    TinyVector<double,2>& Y130 = F(ix+1,iy+1,iz,i)[2];    //  dF/dy
    TinyVector<double,2>& Y131 = F(ix+1,iy+1,iz+1,i)[2];  //  dF/dy
    TinyVector<double,2>& Y132 = F(ix+1,iy+1,iz,i)[6];    // d2F/dydz
    TinyVector<double,2>& Y133 = F(ix+1,iy+1,iz+1,i)[6];  // d2F/dydz
    
    TinyVector<double,2>& Y200 = F(ix,iy,iz,i)[1];      //  dF/dx
    TinyVector<double,2>& Y201 = F(ix,iy,iz+1,i)[1];    //  dF/dx
    TinyVector<double,2>& Y202 = F(ix,iy,iz,i)[5];      // d2F/dxdz
    TinyVector<double,2>& Y203 = F(ix,iy,iz+1,i)[5];    // d2F/dxdz
    TinyVector<double,2>& Y210 = F(ix,iy+1,iz,i)[1];    //  dF/dx
    TinyVector<double,2>& Y211 = F(ix,iy+1,iz+1,i)[1];  //  dF/dx
    TinyVector<double,2>& Y212 = F(ix,iy+1,iz,i)[5];    // d2F/dxdz
    TinyVector<double,2>& Y213 = F(ix,iy+1,iz+1,i)[5];  // d2F/dxdz
    TinyVector<double,2>& Y220 = F(ix,iy,iz,i)[4];      // d2F/dxdy
    TinyVector<double,2>& Y221 = F(ix,iy,iz+1,i)[4];    // d2F/dxdy
    TinyVector<double,2>& Y222 = F(ix,iy,iz,i)[7];      // d3F/dxdydz
    TinyVector<double,2>& Y223 = F(ix,iy,iz+1,i)[7];    // d3F/dxdydz
    TinyVector<double,2>& Y230 = F(ix,iy+1,iz,i)[4];    // d2F/dxdy
    TinyVector<double,2>& Y231 = F(ix,iy+1,iz+1,i)[4];  // d2F/dxdy
    TinyVector<double,2>& Y232 = F(ix,iy+1,iz,i)[7];    // d3F/dxdydz
    TinyVector<double,2>& Y233 = F(ix,iy+1,iz+1,i)[7];  // d3F/dxdydz
    
    TinyVector<double,2>& Y300 = F(ix+1,iy,iz,i)[1];      //  dF/dx
    TinyVector<double,2>& Y301 = F(ix+1,iy,iz+1,i)[1];    //  dF/dx
    TinyVector<double,2>& Y302 = F(ix+1,iy,iz,i)[5];      // d2F/dxdz
    TinyVector<double,2>& Y303 = F(ix+1,iy,iz+1,i)[5];    // d2F/dxdz
    TinyVector<double,2>& Y310 = F(ix+1,iy+1,iz,i)[1];    //  dF/dx
    TinyVector<double,2>& Y311 = F(ix+1,iy+1,iz+1,i)[1];  //  dF/dx
    TinyVector<double,2>& Y312 = F(ix+1,iy+1,iz,i)[5];    // d2F/dxdz
    TinyVector<double,2>& Y313 = F(ix+1,iy+1,iz+1,i)[5];  // d2F/dxdz
    TinyVector<double,2>& Y320 = F(ix+1,iy,iz,i)[4];      // d2F/dxdy
    TinyVector<double,2>& Y321 = F(ix+1,iy,iz+1,i)[4];    // d2F/dxdy
    TinyVector<double,2>& Y322 = F(ix+1,iy,iz,i)[7];      // d3F/dxdydz
    TinyVector<double,2>& Y323 = F(ix+1,iy,iz+1,i)[7];    // d3F/dxdydz
    TinyVector<double,2>& Y330 = F(ix+1,iy+1,iz,i)[4];    // d2F/dxdy
    TinyVector<double,2>& Y331 = F(ix+1,iy+1,iz+1,i)[4];  // d2F/dxdy
    TinyVector<double,2>& Y332 = F(ix+1,iy+1,iz,i)[7];    // d3F/dxdydz
    TinyVector<double,2>& Y333 = F(ix+1,iy+1,iz+1,i)[7];  // d3F/dxdydz
    
    TinyVector<double,2> Val =
      da0*
      (db0*(Y000*c0+Y001*c1+Y002*c2+Y003*c3) +
       db1*(Y010*c0+Y011*c1+Y012*c2+Y013*c3) +
       db2*(Y020*c0+Y021*c1+Y022*c2+Y023*c3) +
       db3*(Y030*c0+Y031*c1+Y032*c2+Y033*c3))+
      da1 *
      (db0*(Y100*c0+Y101*c1+Y102*c2+Y103*c3) +
       db1*(Y110*c0+Y111*c1+Y112*c2+Y113*c3) +
       db2*(Y120*c0+Y121*c1+Y122*c2+Y123*c3) +
       db3*(Y130*c0+Y131*c1+Y132*c2+Y133*c3))+
      da2 *
      (db0*(Y200*c0+Y201*c1+Y202*c2+Y203*c3) +
       db1*(Y210*c0+Y211*c1+Y212*c2+Y213*c3) +
       db2*(Y220*c0+Y221*c1+Y222*c2+Y223*c3) +
       db3*(Y230*c0+Y231*c1+Y232*c2+Y233*c3))+
      da3 *
      (db0*(Y300*c0+Y301*c1+Y302*c2+Y303*c3) +
       db1*(Y310*c0+Y311*c1+Y312*c2+Y313*c3) +
       db2*(Y320*c0+Y321*c1+Y322*c2+Y323*c3) +
       db3*(Y330*c0+Y331*c1+Y332*c2+Y333*c3));

    vals(i) = complex<double>(Val[0], Val[1]);
  }
}

inline void 
ComplexMultiTricubicSpline::d2_dxdz (double x, double y, double z, 
				     Array<complex<double>,1> &vals)
{
  if (!UpToDate)
    Update();

  int ix = Xgrid->ReverseMap(x);  
  int iy = Ygrid->ReverseMap(y);
  int iz = Zgrid->ReverseMap(z);

  ix = max(0,ix); ix = min(ix, Nx-2);
  iy = max(0,iy); iy = min(iy, Ny-2);
  iz = max(0,iz); iz = min(iz, Nz-2);

  double h = (*Xgrid)(ix+1) - (*Xgrid)(ix);
  double hinv = 1.0/h;
  double k = (*Ygrid)(iy+1) - (*Ygrid)(iy);
  double l = (*Zgrid)(iz+1) - (*Zgrid)(iz);
  double linv = 1.0/l;
  double u = (x - (*Xgrid)(ix))/h;
  double v = (y - (*Ygrid)(iy))/k;
  double w = (z - (*Zgrid)(iz))/l;

  double da0 = hinv*dp1(u);
  double da1 = hinv*dp2(u);
  double da2 = dq1(u);
  double da3 = dq2(u);

  register double b0 = p1(v);
  register double b1 = p2(v);
  register double b2 = k*q1(v);
  register double b3 = k*q2(v);

  register double dc0 = linv*dp1(w);
  register double dc1 = linv*dp2(w);
  register double dc2 = dq1(w);
  register double dc3 = dq2(w);


  for (int i=0; i<N; i++) {
    TinyVector<double,2>& Y000 = F(ix,iy,iz,i)[0];      //   F
    TinyVector<double,2>& Y001 = F(ix,iy,iz+1,i)[0];    //   F
    TinyVector<double,2>& Y002 = F(ix,iy,iz,i)[3];      //  dF/dz
    TinyVector<double,2>& Y003 = F(ix,iy,iz+1,i)[3];    //  dF/dz
    TinyVector<double,2>& Y010 = F(ix,iy+1,iz,i)[0];    //   F
    TinyVector<double,2>& Y011 = F(ix,iy+1,iz+1,i)[0];  //   F
    TinyVector<double,2>& Y012 = F(ix,iy+1,iz,i)[3];    //  dF/dz
    TinyVector<double,2>& Y013 = F(ix,iy+1,iz+1,i)[3];  //  dF/dz
    TinyVector<double,2>& Y020 = F(ix,iy,iz,i)[2];      //  dF/dy
    TinyVector<double,2>& Y021 = F(ix,iy,iz+1,i)[2];    //  dF/dy
    TinyVector<double,2>& Y022 = F(ix,iy,iz,i)[6];      // d2F/dydz
    TinyVector<double,2>& Y023 = F(ix,iy,iz+1,i)[6];    // d2F/dydz
    TinyVector<double,2>& Y030 = F(ix,iy+1,iz,i)[2];    //  dF/dy
    TinyVector<double,2>& Y031 = F(ix,iy+1,iz+1,i)[2];  //  dF/dy
    TinyVector<double,2>& Y032 = F(ix,iy+1,iz,i)[6];    // d2F/dydz
    TinyVector<double,2>& Y033 = F(ix,iy+1,iz+1,i)[6];  // d2F/dydz
    
    TinyVector<double,2>& Y100 = F(ix+1,iy,iz,i)[0];      //   F
    TinyVector<double,2>& Y101 = F(ix+1,iy,iz+1,i)[0];    //   F
    TinyVector<double,2>& Y102 = F(ix+1,iy,iz,i)[3];      //  dF/dz
    TinyVector<double,2>& Y103 = F(ix+1,iy,iz+1,i)[3];    //  dF/dz
    TinyVector<double,2>& Y110 = F(ix+1,iy+1,iz,i)[0];    //   F
    TinyVector<double,2>& Y111 = F(ix+1,iy+1,iz+1,i)[0];  //   F
    TinyVector<double,2>& Y112 = F(ix+1,iy+1,iz,i)[3];    //  dF/dz
    TinyVector<double,2>& Y113 = F(ix+1,iy+1,iz+1,i)[3];  //  dF/dz
    TinyVector<double,2>& Y120 = F(ix+1,iy,iz,i)[2];      //  dF/dy
    TinyVector<double,2>& Y121 = F(ix+1,iy,iz+1,i)[2];    //  dF/dy
    TinyVector<double,2>& Y122 = F(ix+1,iy,iz,i)[6];      // d2F/dydz
    TinyVector<double,2>& Y123 = F(ix+1,iy,iz+1,i)[6];    // d2F/dydz
    TinyVector<double,2>& Y130 = F(ix+1,iy+1,iz,i)[2];    //  dF/dy
    TinyVector<double,2>& Y131 = F(ix+1,iy+1,iz+1,i)[2];  //  dF/dy
    TinyVector<double,2>& Y132 = F(ix+1,iy+1,iz,i)[6];    // d2F/dydz
    TinyVector<double,2>& Y133 = F(ix+1,iy+1,iz+1,i)[6];  // d2F/dydz
    
    TinyVector<double,2>& Y200 = F(ix,iy,iz,i)[1];      //  dF/dx
    TinyVector<double,2>& Y201 = F(ix,iy,iz+1,i)[1];    //  dF/dx
    TinyVector<double,2>& Y202 = F(ix,iy,iz,i)[5];      // d2F/dxdz
    TinyVector<double,2>& Y203 = F(ix,iy,iz+1,i)[5];    // d2F/dxdz
    TinyVector<double,2>& Y210 = F(ix,iy+1,iz,i)[1];    //  dF/dx
    TinyVector<double,2>& Y211 = F(ix,iy+1,iz+1,i)[1];  //  dF/dx
    TinyVector<double,2>& Y212 = F(ix,iy+1,iz,i)[5];    // d2F/dxdz
    TinyVector<double,2>& Y213 = F(ix,iy+1,iz+1,i)[5];  // d2F/dxdz
    TinyVector<double,2>& Y220 = F(ix,iy,iz,i)[4];      // d2F/dxdy
    TinyVector<double,2>& Y221 = F(ix,iy,iz+1,i)[4];    // d2F/dxdy
    TinyVector<double,2>& Y222 = F(ix,iy,iz,i)[7];      // d3F/dxdydz
    TinyVector<double,2>& Y223 = F(ix,iy,iz+1,i)[7];    // d3F/dxdydz
    TinyVector<double,2>& Y230 = F(ix,iy+1,iz,i)[4];    // d2F/dxdy
    TinyVector<double,2>& Y231 = F(ix,iy+1,iz+1,i)[4];  // d2F/dxdy
    TinyVector<double,2>& Y232 = F(ix,iy+1,iz,i)[7];    // d3F/dxdydz
    TinyVector<double,2>& Y233 = F(ix,iy+1,iz+1,i)[7];  // d3F/dxdydz
    
    TinyVector<double,2>& Y300 = F(ix+1,iy,iz,i)[1];      //  dF/dx
    TinyVector<double,2>& Y301 = F(ix+1,iy,iz+1,i)[1];    //  dF/dx
    TinyVector<double,2>& Y302 = F(ix+1,iy,iz,i)[5];      // d2F/dxdz
    TinyVector<double,2>& Y303 = F(ix+1,iy,iz+1,i)[5];    // d2F/dxdz
    TinyVector<double,2>& Y310 = F(ix+1,iy+1,iz,i)[1];    //  dF/dx
    TinyVector<double,2>& Y311 = F(ix+1,iy+1,iz+1,i)[1];  //  dF/dx
    TinyVector<double,2>& Y312 = F(ix+1,iy+1,iz,i)[5];    // d2F/dxdz
    TinyVector<double,2>& Y313 = F(ix+1,iy+1,iz+1,i)[5];  // d2F/dxdz
    TinyVector<double,2>& Y320 = F(ix+1,iy,iz,i)[4];      // d2F/dxdy
    TinyVector<double,2>& Y321 = F(ix+1,iy,iz+1,i)[4];    // d2F/dxdy
    TinyVector<double,2>& Y322 = F(ix+1,iy,iz,i)[7];      // d3F/dxdydz
    TinyVector<double,2>& Y323 = F(ix+1,iy,iz+1,i)[7];    // d3F/dxdydz
    TinyVector<double,2>& Y330 = F(ix+1,iy+1,iz,i)[4];    // d2F/dxdy
    TinyVector<double,2>& Y331 = F(ix+1,iy+1,iz+1,i)[4];  // d2F/dxdy
    TinyVector<double,2>& Y332 = F(ix+1,iy+1,iz,i)[7];    // d3F/dxdydz
    TinyVector<double,2>& Y333 = F(ix+1,iy+1,iz+1,i)[7];  // d3F/dxdydz

    TinyVector<double,2> Val = 
      da0*
      (b0*(Y000*dc0+Y001*dc1+Y002*dc2+Y003*dc3) +
       b1*(Y010*dc0+Y011*dc1+Y012*dc2+Y013*dc3) +
       b2*(Y020*dc0+Y021*dc1+Y022*dc2+Y023*dc3) +
       b3*(Y030*dc0+Y031*dc1+Y032*dc2+Y033*dc3))+
      da1 *
      (b0*(Y100*dc0+Y101*dc1+Y102*dc2+Y103*dc3) +
       b1*(Y110*dc0+Y111*dc1+Y112*dc2+Y113*dc3) +
       b2*(Y120*dc0+Y121*dc1+Y122*dc2+Y123*dc3) +
       b3*(Y130*dc0+Y131*dc1+Y132*dc2+Y133*dc3))+
      da2 *
      (b0*(Y200*dc0+Y201*dc1+Y202*dc2+Y203*dc3) +
       b1*(Y210*dc0+Y211*dc1+Y212*dc2+Y213*dc3) +
       b2*(Y220*dc0+Y221*dc1+Y222*dc2+Y223*dc3) +
       b3*(Y230*dc0+Y231*dc1+Y232*dc2+Y233*dc3))+
      da3 *
      (b0*(Y300*dc0+Y301*dc1+Y302*dc2+Y303*dc3) +
       b1*(Y310*dc0+Y311*dc1+Y312*dc2+Y313*dc3) +
       b2*(Y320*dc0+Y321*dc1+Y322*dc2+Y323*dc3) +
       b3*(Y330*dc0+Y331*dc1+Y332*dc2+Y333*dc3));
    vals(i) = complex<double>(Val[0], Val[1]);
  }
}


inline void 
ComplexMultiTricubicSpline::d2_dydz (double x, double y, double z, 
				     Array<complex<double>,1> &vals)
{
  if (!UpToDate)
    Update();

  int ix = Xgrid->ReverseMap(x);  
  int iy = Ygrid->ReverseMap(y);
  int iz = Zgrid->ReverseMap(z);

  ix = max(0,ix); ix = min(ix, Nx-2);
  iy = max(0,iy); iy = min(iy, Ny-2);
  iz = max(0,iz); iz = min(iz, Nz-2);

  double h = (*Xgrid)(ix+1) - (*Xgrid)(ix);
  double k = (*Ygrid)(iy+1) - (*Ygrid)(iy);
  double kinv = 1.0/k;
  double l = (*Zgrid)(iz+1) - (*Zgrid)(iz);
  double linv = 1.0/l;
  double u = (x - (*Xgrid)(ix))/h;
  double v = (y - (*Ygrid)(iy))/k;
  double w = (z - (*Zgrid)(iz))/l;

  double a0 = p1(u);
  double a1 = p2(u);
  double a2 = h*q1(u);
  double a3 = h*q2(u);

  register double db0 = kinv*dp1(v);
  register double db1 = kinv*dp2(v);
  register double db2 = dq1(v);
  register double db3 = dq2(v);

  register double dc0 = linv*dp1(w);
  register double dc1 = linv*dp2(w);
  register double dc2 = dq1(w);
  register double dc3 = dq2(w);

  for (int i=0; i<N; i++) {
    TinyVector<double,2>& Y000 = F(ix,iy,iz,i)[0];      //   F
    TinyVector<double,2>& Y001 = F(ix,iy,iz+1,i)[0];    //   F
    TinyVector<double,2>& Y002 = F(ix,iy,iz,i)[3];      //  dF/dz
    TinyVector<double,2>& Y003 = F(ix,iy,iz+1,i)[3];    //  dF/dz
    TinyVector<double,2>& Y010 = F(ix,iy+1,iz,i)[0];    //   F
    TinyVector<double,2>& Y011 = F(ix,iy+1,iz+1,i)[0];  //   F
    TinyVector<double,2>& Y012 = F(ix,iy+1,iz,i)[3];    //  dF/dz
    TinyVector<double,2>& Y013 = F(ix,iy+1,iz+1,i)[3];  //  dF/dz
    TinyVector<double,2>& Y020 = F(ix,iy,iz,i)[2];      //  dF/dy
    TinyVector<double,2>& Y021 = F(ix,iy,iz+1,i)[2];    //  dF/dy
    TinyVector<double,2>& Y022 = F(ix,iy,iz,i)[6];      // d2F/dydz
    TinyVector<double,2>& Y023 = F(ix,iy,iz+1,i)[6];    // d2F/dydz
    TinyVector<double,2>& Y030 = F(ix,iy+1,iz,i)[2];    //  dF/dy
    TinyVector<double,2>& Y031 = F(ix,iy+1,iz+1,i)[2];  //  dF/dy
    TinyVector<double,2>& Y032 = F(ix,iy+1,iz,i)[6];    // d2F/dydz
    TinyVector<double,2>& Y033 = F(ix,iy+1,iz+1,i)[6];  // d2F/dydz
    
    TinyVector<double,2>& Y100 = F(ix+1,iy,iz,i)[0];      //   F
    TinyVector<double,2>& Y101 = F(ix+1,iy,iz+1,i)[0];    //   F
    TinyVector<double,2>& Y102 = F(ix+1,iy,iz,i)[3];      //  dF/dz
    TinyVector<double,2>& Y103 = F(ix+1,iy,iz+1,i)[3];    //  dF/dz
    TinyVector<double,2>& Y110 = F(ix+1,iy+1,iz,i)[0];    //   F
    TinyVector<double,2>& Y111 = F(ix+1,iy+1,iz+1,i)[0];  //   F
    TinyVector<double,2>& Y112 = F(ix+1,iy+1,iz,i)[3];    //  dF/dz
    TinyVector<double,2>& Y113 = F(ix+1,iy+1,iz+1,i)[3];  //  dF/dz
    TinyVector<double,2>& Y120 = F(ix+1,iy,iz,i)[2];      //  dF/dy
    TinyVector<double,2>& Y121 = F(ix+1,iy,iz+1,i)[2];    //  dF/dy
    TinyVector<double,2>& Y122 = F(ix+1,iy,iz,i)[6];      // d2F/dydz
    TinyVector<double,2>& Y123 = F(ix+1,iy,iz+1,i)[6];    // d2F/dydz
    TinyVector<double,2>& Y130 = F(ix+1,iy+1,iz,i)[2];    //  dF/dy
    TinyVector<double,2>& Y131 = F(ix+1,iy+1,iz+1,i)[2];  //  dF/dy
    TinyVector<double,2>& Y132 = F(ix+1,iy+1,iz,i)[6];    // d2F/dydz
    TinyVector<double,2>& Y133 = F(ix+1,iy+1,iz+1,i)[6];  // d2F/dydz
    
    TinyVector<double,2>& Y200 = F(ix,iy,iz,i)[1];      //  dF/dx
    TinyVector<double,2>& Y201 = F(ix,iy,iz+1,i)[1];    //  dF/dx
    TinyVector<double,2>& Y202 = F(ix,iy,iz,i)[5];      // d2F/dxdz
    TinyVector<double,2>& Y203 = F(ix,iy,iz+1,i)[5];    // d2F/dxdz
    TinyVector<double,2>& Y210 = F(ix,iy+1,iz,i)[1];    //  dF/dx
    TinyVector<double,2>& Y211 = F(ix,iy+1,iz+1,i)[1];  //  dF/dx
    TinyVector<double,2>& Y212 = F(ix,iy+1,iz,i)[5];    // d2F/dxdz
    TinyVector<double,2>& Y213 = F(ix,iy+1,iz+1,i)[5];  // d2F/dxdz
    TinyVector<double,2>& Y220 = F(ix,iy,iz,i)[4];      // d2F/dxdy
    TinyVector<double,2>& Y221 = F(ix,iy,iz+1,i)[4];    // d2F/dxdy
    TinyVector<double,2>& Y222 = F(ix,iy,iz,i)[7];      // d3F/dxdydz
    TinyVector<double,2>& Y223 = F(ix,iy,iz+1,i)[7];    // d3F/dxdydz
    TinyVector<double,2>& Y230 = F(ix,iy+1,iz,i)[4];    // d2F/dxdy
    TinyVector<double,2>& Y231 = F(ix,iy+1,iz+1,i)[4];  // d2F/dxdy
    TinyVector<double,2>& Y232 = F(ix,iy+1,iz,i)[7];    // d3F/dxdydz
    TinyVector<double,2>& Y233 = F(ix,iy+1,iz+1,i)[7];  // d3F/dxdydz
    
    TinyVector<double,2>& Y300 = F(ix+1,iy,iz,i)[1];      //  dF/dx
    TinyVector<double,2>& Y301 = F(ix+1,iy,iz+1,i)[1];    //  dF/dx
    TinyVector<double,2>& Y302 = F(ix+1,iy,iz,i)[5];      // d2F/dxdz
    TinyVector<double,2>& Y303 = F(ix+1,iy,iz+1,i)[5];    // d2F/dxdz
    TinyVector<double,2>& Y310 = F(ix+1,iy+1,iz,i)[1];    //  dF/dx
    TinyVector<double,2>& Y311 = F(ix+1,iy+1,iz+1,i)[1];  //  dF/dx
    TinyVector<double,2>& Y312 = F(ix+1,iy+1,iz,i)[5];    // d2F/dxdz
    TinyVector<double,2>& Y313 = F(ix+1,iy+1,iz+1,i)[5];  // d2F/dxdz
    TinyVector<double,2>& Y320 = F(ix+1,iy,iz,i)[4];      // d2F/dxdy
    TinyVector<double,2>& Y321 = F(ix+1,iy,iz+1,i)[4];    // d2F/dxdy
    TinyVector<double,2>& Y322 = F(ix+1,iy,iz,i)[7];      // d3F/dxdydz
    TinyVector<double,2>& Y323 = F(ix+1,iy,iz+1,i)[7];    // d3F/dxdydz
    TinyVector<double,2>& Y330 = F(ix+1,iy+1,iz,i)[4];    // d2F/dxdy
    TinyVector<double,2>& Y331 = F(ix+1,iy+1,iz+1,i)[4];  // d2F/dxdy
    TinyVector<double,2>& Y332 = F(ix+1,iy+1,iz,i)[7];    // d3F/dxdydz
    TinyVector<double,2>& Y333 = F(ix+1,iy+1,iz+1,i)[7];  // d3F/dxdydz

    TinyVector<double,2> Val = 
      a0*
      (db0*(Y000*dc0+Y001*dc1+Y002*dc2+Y003*dc3) +
       db1*(Y010*dc0+Y011*dc1+Y012*dc2+Y013*dc3) +
       db2*(Y020*dc0+Y021*dc1+Y022*dc2+Y023*dc3) +
       db3*(Y030*dc0+Y031*dc1+Y032*dc2+Y033*dc3))+
      a1 *
      (db0*(Y100*dc0+Y101*dc1+Y102*dc2+Y103*dc3) +
       db1*(Y110*dc0+Y111*dc1+Y112*dc2+Y113*dc3) +
       db2*(Y120*dc0+Y121*dc1+Y122*dc2+Y123*dc3) +
       db3*(Y130*dc0+Y131*dc1+Y132*dc2+Y133*dc3))+
      a2 *
      (db0*(Y200*dc0+Y201*dc1+Y202*dc2+Y203*dc3) +
       db1*(Y210*dc0+Y211*dc1+Y212*dc2+Y213*dc3) +
       db2*(Y220*dc0+Y221*dc1+Y222*dc2+Y223*dc3) +
       db3*(Y230*dc0+Y231*dc1+Y232*dc2+Y233*dc3))+
      a3 *
      (db0*(Y300*dc0+Y301*dc1+Y302*dc2+Y303*dc3) +
       db1*(Y310*dc0+Y311*dc1+Y312*dc2+Y313*dc3) +
       db2*(Y320*dc0+Y321*dc1+Y322*dc2+Y323*dc3) +
       db3*(Y330*dc0+Y331*dc1+Y332*dc2+Y333*dc3));
    vals(i) = complex<double>(Val[0], Val[1]);
  }
}

template<typename T, int M, int N>
TinyMatrix<T,M,N> operator*(T val, const TinyMatrix<T,M,N>& mat)
{
  TinyMatrix<T,M,N> prod;
  for (int i=0; i<M; i++)
    for (int j=0; j<N; j++)
      prod(i,j) = val*mat(i,j);
  return prod;
}

inline void
ComplexMultiTricubicSpline::KineticEnergy(Array<double,1> &KE)
{
  if (!UpToDate)
    Update();
  KE = 0.0;
  
  Array<complex<double>,3> Y(4,4,4);

  for (int ix=0; ix<(Nx-1); ix++) {
    for (int iy=0; iy<(Ny-1); iy++) {
      for (int iz=0; iz<(Nz-1); iz++) {
	TinyMatrix<double,4,4> IntMatX, IntMatY, IntMatZ, d2IntMatX, d2IntMatY, d2IntMatZ;
	double h = (*Xgrid)(ix+1) - (*Xgrid)(ix);
	double hinv = 1.0/h;
	double k = (*Ygrid)(iy+1) - (*Ygrid)(iy);
	double kinv = 1.0/k;
	double l = (*Zgrid)(iz+1) - (*Zgrid)(iz);
	double linv = 1.0/l;

	IntMatX =
	  0.37142857142857144126  ,     0.12857142857142855874,
	  0.05238095238095241690*h,    -0.03095238095238106446*h,
	  0.12857142857142855874,       0.37142857142857144126,   
	  0.03095238095238095344*h,    -0.05238095238095230588*h,
	  0.05238095238095241690*h,     0.03095238095238095344*h,   
	  0.00952380952380943446*h*h,  -0.00714285714285717299*h*h,
         -0.03095238095238106446*h,    -0.05238095238095230588*h,  
         -0.00714285714285717299*h*h,   0.00952380952380954549*h*h;
	IntMatX = h * IntMatX;
	d2IntMatX = 
         -1.200000000000000,   1.200000000000000,  -1.100000000000000*h,-0.100000000000000*h,
	  1.200000000000000,  -1.200000000000000,   0.100000000000000*h, 1.100000000000000*h,
	 -0.100000000000000*h, 0.100000000000000*h,-0.133333333333333*h*h, 0.033333333333333*h*h,
         -0.100000000000000*h, 0.100000000000000*h, 0.033333333333333*h*h,-0.133333333333333*h*h;
	d2IntMatX = hinv * d2IntMatX;

	IntMatY =
	  0.37142857142857144126  ,     0.12857142857142855874,
	  0.05238095238095241690*k,    -0.03095238095238106446*k,
	  0.12857142857142855874,       0.37142857142857144126,   
	  0.03095238095238095344*k,    -0.05238095238095230588*k,
	  0.05238095238095241690*k,     0.03095238095238095344*k,   
	  0.00952380952380943446*k*k,  -0.00714285714285717299*k*k,
         -0.03095238095238106446*k,    -0.05238095238095230588*k,  
         -0.00714285714285717299*k*k,   0.00952380952380954549*k*k;
	IntMatY = k * IntMatY;
	d2IntMatY = 
         -1.200000000000000,   1.200000000000000,  -1.100000000000000*k,-0.100000000000000*k,
	  1.200000000000000,  -1.200000000000000,   0.100000000000000*k, 1.100000000000000*k,
	 -0.100000000000000*k, 0.100000000000000*k,-0.133333333333333*k*k, 0.033333333333333*k*k,
         -0.100000000000000*k, 0.100000000000000*k, 0.033333333333333*k*k,-0.133333333333333*k*k;
	d2IntMatY = kinv * d2IntMatY;

	IntMatZ =
	  0.37142857142857144126  ,     0.12857142857142855874,
	  0.05238095238095241690*l,    -0.03095238095238106446*l,
	  0.12857142857142855874,       0.37142857142857144126,   
	  0.03095238095238095344*l,    -0.05238095238095230588*l,
	  0.05238095238095241690*l,     0.03095238095238095344*l,   
	  0.00952380952380943446*l*l,  -0.00714285714285717299*l*l,
         -0.03095238095238106446*l,    -0.05238095238095230588*l,  
         -0.00714285714285717299*l*l,   0.00952380952380954549*l*l;
	IntMatZ = l * IntMatZ;
	d2IntMatZ = 
         -1.200000000000000,   1.200000000000000,  -1.100000000000000*l,-0.100000000000000*l,
	  1.200000000000000,  -1.200000000000000,   0.100000000000000*l, 1.100000000000000*l,
	 -0.100000000000000*l, 0.100000000000000*l,-0.133333333333333*l*l, 0.033333333333333*l*l,
         -0.100000000000000*l, 0.100000000000000*l, 0.033333333333333*l*l,-0.133333333333333*l*l;
	d2IntMatZ = linv * d2IntMatZ;

		
	double a0 = h*Int_p1;
	double a1 = h*Int_p2;
	double a2 = h*h*Int_q1;
	double a3 = h*h*Int_q2;
	double d2a0 = hinv*Int_d2p1;
	double d2a1 = hinv*Int_d2p2;
	double d2a2 = Int_d2q1;
	double d2a3 = Int_d2q2;
	
	register double b0 = k*Int_p1;
	register double b1 = k*Int_p2;
	register double b2 = k*k*Int_q1;
	register double b3 = k*k*Int_q2;
	register double d2b0 = kinv*Int_d2p1;
	register double d2b1 = kinv*Int_d2p2;
	register double d2b2 = Int_d2q1;
	register double d2b3 = Int_d2q2;
	
	register double c0 = l*Int_p1;
	register double c1 = l*Int_p2;
	register double c2 = l*l*Int_q1;
	register double c3 = l*l*Int_q2;
	register double d2c0 = linv*Int_d2p1;
	register double d2c1 = linv*Int_d2p2;
	register double d2c2 = Int_d2q1;
	register double d2c3 = Int_d2q2;
	
	for (int i=0; i<N; i++) {
	  Y(0,0,0) = complex<double> (F(ix,iy,iz,i)[0][0]      , F(ix,iy,iz,i)[0][1]      );  //   F
	  Y(0,0,1) = complex<double> (F(ix,iy,iz+1,i)[0][0]    , F(ix,iy,iz+1,i)[0][1]    );  //   F
	  Y(0,0,2) = complex<double> (F(ix,iy,iz,i)[3][0]      , F(ix,iy,iz,i)[3][1]      );  //  dF/dz
	  Y(0,0,3) = complex<double> (F(ix,iy,iz+1,i)[3][0]    , F(ix,iy,iz+1,i)[3][1]    );  //  dF/dz
	  Y(0,1,0) = complex<double> (F(ix,iy+1,iz,i)[0][0]    , F(ix,iy+1,iz,i)[0][1]    );  //   F
	  Y(0,1,1) = complex<double> (F(ix,iy+1,iz+1,i)[0][0]  , F(ix,iy+1,iz+1,i)[0][1]  );  //   F
	  Y(0,1,2) = complex<double> (F(ix,iy+1,iz,i)[3][0]    , F(ix,iy+1,iz,i)[3][1]    );  //  dF/dz
	  Y(0,1,3) = complex<double> (F(ix,iy+1,iz+1,i)[3][0]  , F(ix,iy+1,iz+1,i)[3][1]  );  //  dF/dz
	  Y(0,2,0) = complex<double> (F(ix,iy,iz,i)[2][0]      , F(ix,iy,iz,i)[2][1]      );  //  dF/dy
	  Y(0,2,1) = complex<double> (F(ix,iy,iz+1,i)[2][0]    , F(ix,iy,iz+1,i)[2][1]    );  //  dF/dy
	  Y(0,2,2) = complex<double> (F(ix,iy,iz,i)[6][0]      , F(ix,iy,iz,i)[6][1]      );  // d2F/dydz
	  Y(0,2,3) = complex<double> (F(ix,iy,iz+1,i)[6][0]    , F(ix,iy,iz+1,i)[6][1]    );  // d2F/dydz
	  Y(0,3,0) = complex<double> (F(ix,iy+1,iz,i)[2][0]    , F(ix,iy+1,iz,i)[2][1]    );  //  dF/dy
	  Y(0,3,1) = complex<double> (F(ix,iy+1,iz+1,i)[2][0]  , F(ix,iy+1,iz+1,i)[2][1]  );  //  dF/dy
	  Y(0,3,2) = complex<double> (F(ix,iy+1,iz,i)[6][0]    , F(ix,iy+1,iz,i)[6][1]    );  // d2F/dydz
	  Y(0,3,3) = complex<double> (F(ix,iy+1,iz+1,i)[6][0]  , F(ix,iy+1,iz+1,i)[6][1]  );  // d2F/dydz
	 						       						       
	  Y(1,0,0) = complex<double> (F(ix+1,iy,iz,i)[0][0]    , F(ix+1,iy,iz,i)[0][1]    );  //   F
	  Y(1,0,1) = complex<double> (F(ix+1,iy,iz+1,i)[0][0]  , F(ix+1,iy,iz+1,i)[0][1]  );  //   F
	  Y(1,0,2) = complex<double> (F(ix+1,iy,iz,i)[3][0]    , F(ix+1,iy,iz,i)[3][1]    );  //  dF/dz
	  Y(1,0,3) = complex<double> (F(ix+1,iy,iz+1,i)[3][0]  , F(ix+1,iy,iz+1,i)[3][1]  );  //  dF/dz
	  Y(1,1,0) = complex<double> (F(ix+1,iy+1,iz,i)[0][0]  , F(ix+1,iy+1,iz,i)[0][1]  );  //   F
	  Y(1,1,1) = complex<double> (F(ix+1,iy+1,iz+1,i)[0][0], F(ix+1,iy+1,iz+1,i)[0][1]);  //   , F
	  Y(1,1,2) = complex<double> (F(ix+1,iy+1,iz,i)[3][0]  , F(ix+1,iy+1,iz,i)[3][1]  );  //  dF/dz
	  Y(1,1,3) = complex<double> (F(ix+1,iy+1,iz+1,i)[3][0], F(ix+1,iy+1,iz+1,i)[3][1]);  //  dF/dz
	  Y(1,2,0) = complex<double> (F(ix+1,iy,iz,i)[2][0]    , F(ix+1,iy,iz,i)[2][1]    );  //  dF/dy
	  Y(1,2,1) = complex<double> (F(ix+1,iy,iz+1,i)[2][0]  , F(ix+1,iy,iz+1,i)[2][1]  );  //  dF/dy
	  Y(1,2,2) = complex<double> (F(ix+1,iy,iz,i)[6][0]    , F(ix+1,iy,iz,i)[6][1]    );  // d2F/dydz
	  Y(1,2,3) = complex<double> (F(ix+1,iy,iz+1,i)[6][0]  , F(ix+1,iy,iz+1,i)[6][1]  );  // d2F/dydz
	  Y(1,3,0) = complex<double> (F(ix+1,iy+1,iz,i)[2][0]  , F(ix+1,iy+1,iz,i)[2][1]  );  //  dF/dy
	  Y(1,3,1) = complex<double> (F(ix+1,iy+1,iz+1,i)[2][0], F(ix+1,iy+1,iz+1,i)[2][1]);  //  dF/dy
	  Y(1,3,2) = complex<double> (F(ix+1,iy+1,iz,i)[6][0]  , F(ix+1,iy+1,iz,i)[6][1]  );  // d2F/dydz
	  Y(1,3,3) = complex<double> (F(ix+1,iy+1,iz+1,i)[6][0], F(ix+1,iy+1,iz+1,i)[6][1]);  // d2F/dydz
	 						       						       
	  Y(2,0,0) = complex<double> (F(ix,iy,iz,i)[1][0]      , F(ix,iy,iz,i)[1][1]      );  //  dF/dx
	  Y(2,0,1) = complex<double> (F(ix,iy,iz+1,i)[1][0]    , F(ix,iy,iz+1,i)[1][1]    );  //  dF/dx
	  Y(2,0,2) = complex<double> (F(ix,iy,iz,i)[5][0]      , F(ix,iy,iz,i)[5][1]      );  // d2F/dxdz
	  Y(2,0,3) = complex<double> (F(ix,iy,iz+1,i)[5][0]    , F(ix,iy,iz+1,i)[5][1]    );  // d2F/dxdz
	  Y(2,1,0) = complex<double> (F(ix,iy+1,iz,i)[1][0]    , F(ix,iy+1,iz,i)[1][1]    );  //  dF/dx
	  Y(2,1,1) = complex<double> (F(ix,iy+1,iz+1,i)[1][0]  , F(ix,iy+1,iz+1,i)[1][1]  );  //  dF/dx
	  Y(2,1,2) = complex<double> (F(ix,iy+1,iz,i)[5][0]    , F(ix,iy+1,iz,i)[5][1]    );  // d2F/dxdz
	  Y(2,1,3) = complex<double> (F(ix,iy+1,iz+1,i)[5][0]  , F(ix,iy+1,iz+1,i)[5][1]  );  // d2F/dxdz
	  Y(2,2,0) = complex<double> (F(ix,iy,iz,i)[4][0]      , F(ix,iy,iz,i)[4][1]      );  // d2F/dxdy
	  Y(2,2,1) = complex<double> (F(ix,iy,iz+1,i)[4][0]    , F(ix,iy,iz+1,i)[4][1]    );  // d2F/dxdy
	  Y(2,2,2) = complex<double> (F(ix,iy,iz,i)[7][0]      , F(ix,iy,iz,i)[7][1]      );  // d3F/dxdydz
	  Y(2,2,3) = complex<double> (F(ix,iy,iz+1,i)[7][0]    , F(ix,iy,iz+1,i)[7][1]    );  // d3F/dxdydz
	  Y(2,3,0) = complex<double> (F(ix,iy+1,iz,i)[4][0]    , F(ix,iy+1,iz,i)[4][1]    );  // d2F/dxdy
	  Y(2,3,1) = complex<double> (F(ix,iy+1,iz+1,i)[4][0]  , F(ix,iy+1,iz+1,i)[4][1]  );  // d2F/dxdy
	  Y(2,3,2) = complex<double> (F(ix,iy+1,iz,i)[7][0]    , F(ix,iy+1,iz,i)[7][1]    );  // d3F/dxdydz
	  Y(2,3,3) = complex<double> (F(ix,iy+1,iz+1,i)[7][0]  , F(ix,iy+1,iz+1,i)[7][1]  );  // d3F/dxdydz
	 						       						       
	  Y(3,0,0) = complex<double> (F(ix+1,iy,iz,i)[1][0]    , F(ix+1,iy,iz,i)[1][1]    );  //  dF/dx
	  Y(3,0,1) = complex<double> (F(ix+1,iy,iz+1,i)[1][0]  , F(ix+1,iy,iz+1,i)[1][1]  );  //  dF/dx
	  Y(3,0,2) = complex<double> (F(ix+1,iy,iz,i)[5][0]    , F(ix+1,iy,iz,i)[5][1]    );  // d2F/dxdz
	  Y(3,0,3) = complex<double> (F(ix+1,iy,iz+1,i)[5][0]  , F(ix+1,iy,iz+1,i)[5][1]  );  // d2F/dxdz
	  Y(3,1,0) = complex<double> (F(ix+1,iy+1,iz,i)[1][0]  , F(ix+1,iy+1,iz,i)[1][1]  );  //  dF/dx
	  Y(3,1,1) = complex<double> (F(ix+1,iy+1,iz+1,i)[1][0], F(ix+1,iy+1,iz+1,i)[1][1]);  //  dF/dx
	  Y(3,1,2) = complex<double> (F(ix+1,iy+1,iz,i)[5][0]  , F(ix+1,iy+1,iz,i)[5][1]  );  // d2F/dxdz
	  Y(3,1,3) = complex<double> (F(ix+1,iy+1,iz+1,i)[5][0], F(ix+1,iy+1,iz+1,i)[5][1]);  // d2F/dxdz
	  Y(3,2,0) = complex<double> (F(ix+1,iy,iz,i)[4][0]    , F(ix+1,iy,iz,i)[4][1]    );  // d2F/dxdy
	  Y(3,2,1) = complex<double> (F(ix+1,iy,iz+1,i)[4][0]  , F(ix+1,iy,iz+1,i)[4][1]  );  // d2F/dxdy
	  Y(3,2,2) = complex<double> (F(ix+1,iy,iz,i)[7][0]    , F(ix+1,iy,iz,i)[7][1]    );  // d3F/dxdydz
	  Y(3,2,3) = complex<double> (F(ix+1,iy,iz+1,i)[7][0]  , F(ix+1,iy,iz+1,i)[7][1]  );  // d3F/dxdydz
	  Y(3,3,0) = complex<double> (F(ix+1,iy+1,iz,i)[4][0]  , F(ix+1,iy+1,iz,i)[4][1]  );  // d2F/dxdy
	  Y(3,3,1) = complex<double> (F(ix+1,iy+1,iz+1,i)[4][0], F(ix+1,iy+1,iz+1,i)[4][1]);  // d2F/dxdy
	  Y(3,3,2) = complex<double> (F(ix+1,iy+1,iz,i)[7][0]  , F(ix+1,iy+1,iz,i)[7][1]  );  // d3F/dxdydz
	  Y(3,3,3) = complex<double> (F(ix+1,iy+1,iz+1,i)[7][0], F(ix+1,iy+1,iz+1,i)[7][1]);  // d3F/dxdydz
	  
	  for (int i1=0; i1<4; i1++)
	    for (int i2=0; i2<4; i2++)
	      for (int j1=0; j1<4; j1++)
		for (int j2=0; j2<4; j2++)
		  for (int k1=0; k1<4; k1++)
		    for (int k2=0; k2<4; k2++) {
		      KE(i) += real(d2IntMatX(i1,i2) *  IntMatY(j1,j2) *  IntMatZ(k1,k2)*conj(Y(i1,j1,k1))*Y(i2,j2,k2));
		      KE(i) += real(  IntMatX(i1,i2) *d2IntMatY(j1,j2) *  IntMatZ(k1,k2)*conj(Y(i1,j1,k1))*Y(i2,j2,k2));
		      KE(i) += real(  IntMatX(i1,i2) *  IntMatY(j1,j2) *d2IntMatZ(k1,k2)*conj(Y(i1,j1,k1))*Y(i2,j2,k2));
		    }
		      
// 	  TinyVector<double,2> Val =
// 	    d2a0*
// 	    (b0*(Y000*c0+Y001*c1+Y002*c2+Y003*c3) +
// 	     b1*(Y010*c0+Y011*c1+Y012*c2+Y013*c3) +
// 	     b2*(Y020*c0+Y021*c1+Y022*c2+Y023*c3) +
// 	     b3*(Y030*c0+Y031*c1+Y032*c2+Y033*c3))+
// 	    d2a1 *
// 	    (b0*(Y100*c0+Y101*c1+Y102*c2+Y103*c3) +
// 	     b1*(Y110*c0+Y111*c1+Y112*c2+Y113*c3) +
// 	     b2*(Y120*c0+Y121*c1+Y122*c2+Y123*c3) +
// 	     b3*(Y130*c0+Y131*c1+Y132*c2+Y133*c3))+
// 	    d2a2 *
// 	    (b0*(Y200*c0+Y201*c1+Y202*c2+Y203*c3) +
// 	     b1*(Y210*c0+Y211*c1+Y212*c2+Y213*c3) +
// 	     b2*(Y220*c0+Y221*c1+Y222*c2+Y223*c3) +
// 	     b3*(Y230*c0+Y231*c1+Y232*c2+Y233*c3))+
// 	    d2a3 *
// 	    (b0*(Y300*c0+Y301*c1+Y302*c2+Y303*c3) +
// 	     b1*(Y310*c0+Y311*c1+Y312*c2+Y313*c3) +
// 	     b2*(Y320*c0+Y321*c1+Y322*c2+Y323*c3) +
// 	     b3*(Y330*c0+Y331*c1+Y332*c2+Y333*c3));
	  
// 	  Val +=
// 	    a0*
// 	    (d2b0*(Y000*c0+Y001*c1+Y002*c2+Y003*c3) +
// 	     d2b1*(Y010*c0+Y011*c1+Y012*c2+Y013*c3) +
// 	     d2b2*(Y020*c0+Y021*c1+Y022*c2+Y023*c3) +
// 	     d2b3*(Y030*c0+Y031*c1+Y032*c2+Y033*c3))+
// 	    a1 *
// 	    (d2b0*(Y100*c0+Y101*c1+Y102*c2+Y103*c3) +
// 	     d2b1*(Y110*c0+Y111*c1+Y112*c2+Y113*c3) +
// 	     d2b2*(Y120*c0+Y121*c1+Y122*c2+Y123*c3) +
// 	     d2b3*(Y130*c0+Y131*c1+Y132*c2+Y133*c3))+
// 	    a2 *
// 	    (d2b0*(Y200*c0+Y201*c1+Y202*c2+Y203*c3) +
// 	     d2b1*(Y210*c0+Y211*c1+Y212*c2+Y213*c3) +
// 	     d2b2*(Y220*c0+Y221*c1+Y222*c2+Y223*c3) +
// 	     d2b3*(Y230*c0+Y231*c1+Y232*c2+Y233*c3))+
// 	    a3 *
// 	    (d2b0*(Y300*c0+Y301*c1+Y302*c2+Y303*c3) +
// 	     d2b1*(Y310*c0+Y311*c1+Y312*c2+Y313*c3) +
// 	     d2b2*(Y320*c0+Y321*c1+Y322*c2+Y323*c3) +
// 	     d2b3*(Y330*c0+Y331*c1+Y332*c2+Y333*c3));
	  
// 	  Val += 
// 	    a0*
// 	    (b0*(Y000*d2c0+Y001*d2c1+Y002*d2c2+Y003*d2c3) +
// 	     b1*(Y010*d2c0+Y011*d2c1+Y012*d2c2+Y013*d2c3) +
// 	     b2*(Y020*d2c0+Y021*d2c1+Y022*d2c2+Y023*d2c3) +
// 	     b3*(Y030*d2c0+Y031*d2c1+Y032*d2c2+Y033*d2c3))+
// 	    a1 *
// 	    (b0*(Y100*d2c0+Y101*d2c1+Y102*d2c2+Y103*d2c3) +
// 	     b1*(Y110*d2c0+Y111*d2c1+Y112*d2c2+Y113*d2c3) +
// 	     b2*(Y120*d2c0+Y121*d2c1+Y122*d2c2+Y123*d2c3) +
// 	     b3*(Y130*d2c0+Y131*d2c1+Y132*d2c2+Y133*d2c3))+
// 	    a2 *
// 	    (b0*(Y200*d2c0+Y201*d2c1+Y202*d2c2+Y203*d2c3) +
// 	     b1*(Y210*d2c0+Y211*d2c1+Y212*d2c2+Y213*d2c3) +
// 	     b2*(Y220*d2c0+Y221*d2c1+Y222*d2c2+Y223*d2c3) +
// 	     b3*(Y230*d2c0+Y231*d2c1+Y232*d2c2+Y233*d2c3))+
// 	    a3 *
// 	    (b0*(Y300*d2c0+Y301*d2c1+Y302*d2c2+Y303*d2c3) +
// 	     b1*(Y310*d2c0+Y311*d2c1+Y312*d2c2+Y313*d2c3) +
// 	     b2*(Y320*d2c0+Y321*d2c1+Y322*d2c2+Y323*d2c3) +
// 	     b3*(Y330*d2c0+Y331*d2c1+Y332*d2c2+Y333*d2c3));
	  
// 	  vals(i) += complex<double>(Val[0],Val[1]);
	}
      }
    }
  }
}

inline void
ComplexMultiTricubicSpline::Norm(Array<double,1> &nrm)
{
  if (!UpToDate)
    Update();
  nrm = 0.0;
  
  Array<complex<double>,3> Y(4,4,4);

  for (int ix=0; ix<(Nx-1); ix++) {
    for (int iy=0; iy<(Ny-1); iy++) {
      for (int iz=0; iz<(Nz-1); iz++) {
	TinyMatrix<double,4,4> IntMatX, IntMatY, IntMatZ;
	double h = (*Xgrid)(ix+1) - (*Xgrid)(ix);
	double hinv = 1.0/h;
	double k = (*Ygrid)(iy+1) - (*Ygrid)(iy);
	double kinv = 1.0/k;
	double l = (*Zgrid)(iz+1) - (*Zgrid)(iz);
	double linv = 1.0/l;

	IntMatX =
	  0.37142857142857144126  ,     0.12857142857142855874,
	  0.05238095238095241690*h,    -0.03095238095238106446*h,
	  0.12857142857142855874,       0.37142857142857144126,   
	  0.03095238095238095344*h,    -0.05238095238095230588*h,
	  0.05238095238095241690*h,     0.03095238095238095344*h,   
	  0.00952380952380943446*h*h,  -0.00714285714285717299*h*h,
         -0.03095238095238106446*h,    -0.05238095238095230588*h,  
         -0.00714285714285717299*h*h,   0.00952380952380954549*h*h;
	IntMatX = h * IntMatX;

	IntMatY =
	  0.37142857142857144126  ,     0.12857142857142855874,
	  0.05238095238095241690*k,    -0.03095238095238106446*k,
	  0.12857142857142855874,       0.37142857142857144126,   
	  0.03095238095238095344*k,    -0.05238095238095230588*k,
	  0.05238095238095241690*k,     0.03095238095238095344*k,   
	  0.00952380952380943446*k*k,  -0.00714285714285717299*k*k,
         -0.03095238095238106446*k,    -0.05238095238095230588*k,  
         -0.00714285714285717299*k*k,   0.00952380952380954549*k*k;
	IntMatY = k * IntMatY;

	IntMatZ =
	  0.37142857142857144126  ,     0.12857142857142855874,
	  0.05238095238095241690*l,    -0.03095238095238106446*l,
	  0.12857142857142855874,       0.37142857142857144126,   
	  0.03095238095238095344*l,    -0.05238095238095230588*l,
	  0.05238095238095241690*l,     0.03095238095238095344*l,   
	  0.00952380952380943446*l*l,  -0.00714285714285717299*l*l,
         -0.03095238095238106446*l,    -0.05238095238095230588*l,  
         -0.00714285714285717299*l*l,   0.00952380952380954549*l*l;
	IntMatZ = l * IntMatZ;
	
	for (int i=0; i<N; i++) {
	  Y(0,0,0) = complex<double> (F(ix,iy,iz,i)[0][0]      , F(ix,iy,iz,i)[0][1]      );  //   F
	  Y(0,0,1) = complex<double> (F(ix,iy,iz+1,i)[0][0]    , F(ix,iy,iz+1,i)[0][1]    );  //   F
	  Y(0,0,2) = complex<double> (F(ix,iy,iz,i)[3][0]      , F(ix,iy,iz,i)[3][1]      );  //  dF/dz
	  Y(0,0,3) = complex<double> (F(ix,iy,iz+1,i)[3][0]    , F(ix,iy,iz+1,i)[3][1]    );  //  dF/dz
	  Y(0,1,0) = complex<double> (F(ix,iy+1,iz,i)[0][0]    , F(ix,iy+1,iz,i)[0][1]    );  //   F
	  Y(0,1,1) = complex<double> (F(ix,iy+1,iz+1,i)[0][0]  , F(ix,iy+1,iz+1,i)[0][1]  );  //   F
	  Y(0,1,2) = complex<double> (F(ix,iy+1,iz,i)[3][0]    , F(ix,iy+1,iz,i)[3][1]    );  //  dF/dz
	  Y(0,1,3) = complex<double> (F(ix,iy+1,iz+1,i)[3][0]  , F(ix,iy+1,iz+1,i)[3][1]  );  //  dF/dz
	  Y(0,2,0) = complex<double> (F(ix,iy,iz,i)[2][0]      , F(ix,iy,iz,i)[2][1]      );  //  dF/dy
	  Y(0,2,1) = complex<double> (F(ix,iy,iz+1,i)[2][0]    , F(ix,iy,iz+1,i)[2][1]    );  //  dF/dy
	  Y(0,2,2) = complex<double> (F(ix,iy,iz,i)[6][0]      , F(ix,iy,iz,i)[6][1]      );  // d2F/dydz
	  Y(0,2,3) = complex<double> (F(ix,iy,iz+1,i)[6][0]    , F(ix,iy,iz+1,i)[6][1]    );  // d2F/dydz
	  Y(0,3,0) = complex<double> (F(ix,iy+1,iz,i)[2][0]    , F(ix,iy+1,iz,i)[2][1]    );  //  dF/dy
	  Y(0,3,1) = complex<double> (F(ix,iy+1,iz+1,i)[2][0]  , F(ix,iy+1,iz+1,i)[2][1]  );  //  dF/dy
	  Y(0,3,2) = complex<double> (F(ix,iy+1,iz,i)[6][0]    , F(ix,iy+1,iz,i)[6][1]    );  // d2F/dydz
	  Y(0,3,3) = complex<double> (F(ix,iy+1,iz+1,i)[6][0]  , F(ix,iy+1,iz+1,i)[6][1]  );  // d2F/dydz
	 						       						       
	  Y(1,0,0) = complex<double> (F(ix+1,iy,iz,i)[0][0]    , F(ix+1,iy,iz,i)[0][1]    );  //   F
	  Y(1,0,1) = complex<double> (F(ix+1,iy,iz+1,i)[0][0]  , F(ix+1,iy,iz+1,i)[0][1]  );  //   F
	  Y(1,0,2) = complex<double> (F(ix+1,iy,iz,i)[3][0]    , F(ix+1,iy,iz,i)[3][1]    );  //  dF/dz
	  Y(1,0,3) = complex<double> (F(ix+1,iy,iz+1,i)[3][0]  , F(ix+1,iy,iz+1,i)[3][1]  );  //  dF/dz
	  Y(1,1,0) = complex<double> (F(ix+1,iy+1,iz,i)[0][0]  , F(ix+1,iy+1,iz,i)[0][1]  );  //   F
	  Y(1,1,1) = complex<double> (F(ix+1,iy+1,iz+1,i)[0][0], F(ix+1,iy+1,iz+1,i)[0][1]);  //   , F
	  Y(1,1,2) = complex<double> (F(ix+1,iy+1,iz,i)[3][0]  , F(ix+1,iy+1,iz,i)[3][1]  );  //  dF/dz
	  Y(1,1,3) = complex<double> (F(ix+1,iy+1,iz+1,i)[3][0], F(ix+1,iy+1,iz+1,i)[3][1]);  //  dF/dz
	  Y(1,2,0) = complex<double> (F(ix+1,iy,iz,i)[2][0]    , F(ix+1,iy,iz,i)[2][1]    );  //  dF/dy
	  Y(1,2,1) = complex<double> (F(ix+1,iy,iz+1,i)[2][0]  , F(ix+1,iy,iz+1,i)[2][1]  );  //  dF/dy
	  Y(1,2,2) = complex<double> (F(ix+1,iy,iz,i)[6][0]    , F(ix+1,iy,iz,i)[6][1]    );  // d2F/dydz
	  Y(1,2,3) = complex<double> (F(ix+1,iy,iz+1,i)[6][0]  , F(ix+1,iy,iz+1,i)[6][1]  );  // d2F/dydz
	  Y(1,3,0) = complex<double> (F(ix+1,iy+1,iz,i)[2][0]  , F(ix+1,iy+1,iz,i)[2][1]  );  //  dF/dy
	  Y(1,3,1) = complex<double> (F(ix+1,iy+1,iz+1,i)[2][0], F(ix+1,iy+1,iz+1,i)[2][1]);  //  dF/dy
	  Y(1,3,2) = complex<double> (F(ix+1,iy+1,iz,i)[6][0]  , F(ix+1,iy+1,iz,i)[6][1]  );  // d2F/dydz
	  Y(1,3,3) = complex<double> (F(ix+1,iy+1,iz+1,i)[6][0], F(ix+1,iy+1,iz+1,i)[6][1]);  // d2F/dydz
	 						       						       
	  Y(2,0,0) = complex<double> (F(ix,iy,iz,i)[1][0]      , F(ix,iy,iz,i)[1][1]      );  //  dF/dx
	  Y(2,0,1) = complex<double> (F(ix,iy,iz+1,i)[1][0]    , F(ix,iy,iz+1,i)[1][1]    );  //  dF/dx
	  Y(2,0,2) = complex<double> (F(ix,iy,iz,i)[5][0]      , F(ix,iy,iz,i)[5][1]      );  // d2F/dxdz
	  Y(2,0,3) = complex<double> (F(ix,iy,iz+1,i)[5][0]    , F(ix,iy,iz+1,i)[5][1]    );  // d2F/dxdz
	  Y(2,1,0) = complex<double> (F(ix,iy+1,iz,i)[1][0]    , F(ix,iy+1,iz,i)[1][1]    );  //  dF/dx
	  Y(2,1,1) = complex<double> (F(ix,iy+1,iz+1,i)[1][0]  , F(ix,iy+1,iz+1,i)[1][1]  );  //  dF/dx
	  Y(2,1,2) = complex<double> (F(ix,iy+1,iz,i)[5][0]    , F(ix,iy+1,iz,i)[5][1]    );  // d2F/dxdz
	  Y(2,1,3) = complex<double> (F(ix,iy+1,iz+1,i)[5][0]  , F(ix,iy+1,iz+1,i)[5][1]  );  // d2F/dxdz
	  Y(2,2,0) = complex<double> (F(ix,iy,iz,i)[4][0]      , F(ix,iy,iz,i)[4][1]      );  // d2F/dxdy
	  Y(2,2,1) = complex<double> (F(ix,iy,iz+1,i)[4][0]    , F(ix,iy,iz+1,i)[4][1]    );  // d2F/dxdy
	  Y(2,2,2) = complex<double> (F(ix,iy,iz,i)[7][0]      , F(ix,iy,iz,i)[7][1]      );  // d3F/dxdydz
	  Y(2,2,3) = complex<double> (F(ix,iy,iz+1,i)[7][0]    , F(ix,iy,iz+1,i)[7][1]    );  // d3F/dxdydz
	  Y(2,3,0) = complex<double> (F(ix,iy+1,iz,i)[4][0]    , F(ix,iy+1,iz,i)[4][1]    );  // d2F/dxdy
	  Y(2,3,1) = complex<double> (F(ix,iy+1,iz+1,i)[4][0]  , F(ix,iy+1,iz+1,i)[4][1]  );  // d2F/dxdy
	  Y(2,3,2) = complex<double> (F(ix,iy+1,iz,i)[7][0]    , F(ix,iy+1,iz,i)[7][1]    );  // d3F/dxdydz
	  Y(2,3,3) = complex<double> (F(ix,iy+1,iz+1,i)[7][0]  , F(ix,iy+1,iz+1,i)[7][1]  );  // d3F/dxdydz
	 						       						       
	  Y(3,0,0) = complex<double> (F(ix+1,iy,iz,i)[1][0]    , F(ix+1,iy,iz,i)[1][1]    );  //  dF/dx
	  Y(3,0,1) = complex<double> (F(ix+1,iy,iz+1,i)[1][0]  , F(ix+1,iy,iz+1,i)[1][1]  );  //  dF/dx
	  Y(3,0,2) = complex<double> (F(ix+1,iy,iz,i)[5][0]    , F(ix+1,iy,iz,i)[5][1]    );  // d2F/dxdz
	  Y(3,0,3) = complex<double> (F(ix+1,iy,iz+1,i)[5][0]  , F(ix+1,iy,iz+1,i)[5][1]  );  // d2F/dxdz
	  Y(3,1,0) = complex<double> (F(ix+1,iy+1,iz,i)[1][0]  , F(ix+1,iy+1,iz,i)[1][1]  );  //  dF/dx
	  Y(3,1,1) = complex<double> (F(ix+1,iy+1,iz+1,i)[1][0], F(ix+1,iy+1,iz+1,i)[1][1]);  //  dF/dx
	  Y(3,1,2) = complex<double> (F(ix+1,iy+1,iz,i)[5][0]  , F(ix+1,iy+1,iz,i)[5][1]  );  // d2F/dxdz
	  Y(3,1,3) = complex<double> (F(ix+1,iy+1,iz+1,i)[5][0], F(ix+1,iy+1,iz+1,i)[5][1]);  // d2F/dxdz
	  Y(3,2,0) = complex<double> (F(ix+1,iy,iz,i)[4][0]    , F(ix+1,iy,iz,i)[4][1]    );  // d2F/dxdy
	  Y(3,2,1) = complex<double> (F(ix+1,iy,iz+1,i)[4][0]  , F(ix+1,iy,iz+1,i)[4][1]  );  // d2F/dxdy
	  Y(3,2,2) = complex<double> (F(ix+1,iy,iz,i)[7][0]    , F(ix+1,iy,iz,i)[7][1]    );  // d3F/dxdydz
	  Y(3,2,3) = complex<double> (F(ix+1,iy,iz+1,i)[7][0]  , F(ix+1,iy,iz+1,i)[7][1]  );  // d3F/dxdydz
	  Y(3,3,0) = complex<double> (F(ix+1,iy+1,iz,i)[4][0]  , F(ix+1,iy+1,iz,i)[4][1]  );  // d2F/dxdy
	  Y(3,3,1) = complex<double> (F(ix+1,iy+1,iz+1,i)[4][0], F(ix+1,iy+1,iz+1,i)[4][1]);  // d2F/dxdy
	  Y(3,3,2) = complex<double> (F(ix+1,iy+1,iz,i)[7][0]  , F(ix+1,iy+1,iz,i)[7][1]  );  // d3F/dxdydz
	  Y(3,3,3) = complex<double> (F(ix+1,iy+1,iz+1,i)[7][0], F(ix+1,iy+1,iz+1,i)[7][1]);  // d3F/dxdydz
	  
	  for (int i1=0; i1<4; i1++)
	    for (int i2=0; i2<4; i2++)
	      for (int j1=0; j1<4; j1++)
		for (int j2=0; j2<4; j2++)
		  for (int k1=0; k1<4; k1++)
		    for (int k2=0; k2<4; k2++) 
		      nrm(i) += real(IntMatX(i1,i2) *  IntMatY(j1,j2) *  IntMatZ(k1,k2)*conj(Y(i1,j1,k1))*Y(i2,j2,k2));
		      
// 	  TinyVector<double,2> Val =
// 	    d2a0*
// 	    (b0*(Y000*c0+Y001*c1+Y002*c2+Y003*c3) +
// 	     b1*(Y010*c0+Y011*c1+Y012*c2+Y013*c3) +
// 	     b2*(Y020*c0+Y021*c1+Y022*c2+Y023*c3) +
// 	     b3*(Y030*c0+Y031*c1+Y032*c2+Y033*c3))+
// 	    d2a1 *
// 	    (b0*(Y100*c0+Y101*c1+Y102*c2+Y103*c3) +
// 	     b1*(Y110*c0+Y111*c1+Y112*c2+Y113*c3) +
// 	     b2*(Y120*c0+Y121*c1+Y122*c2+Y123*c3) +
// 	     b3*(Y130*c0+Y131*c1+Y132*c2+Y133*c3))+
// 	    d2a2 *
// 	    (b0*(Y200*c0+Y201*c1+Y202*c2+Y203*c3) +
// 	     b1*(Y210*c0+Y211*c1+Y212*c2+Y213*c3) +
// 	     b2*(Y220*c0+Y221*c1+Y222*c2+Y223*c3) +
// 	     b3*(Y230*c0+Y231*c1+Y232*c2+Y233*c3))+
// 	    d2a3 *
// 	    (b0*(Y300*c0+Y301*c1+Y302*c2+Y303*c3) +
// 	     b1*(Y310*c0+Y311*c1+Y312*c2+Y313*c3) +
// 	     b2*(Y320*c0+Y321*c1+Y322*c2+Y323*c3) +
// 	     b3*(Y330*c0+Y331*c1+Y332*c2+Y333*c3));
	  
// 	  Val +=
// 	    a0*
// 	    (d2b0*(Y000*c0+Y001*c1+Y002*c2+Y003*c3) +
// 	     d2b1*(Y010*c0+Y011*c1+Y012*c2+Y013*c3) +
// 	     d2b2*(Y020*c0+Y021*c1+Y022*c2+Y023*c3) +
// 	     d2b3*(Y030*c0+Y031*c1+Y032*c2+Y033*c3))+
// 	    a1 *
// 	    (d2b0*(Y100*c0+Y101*c1+Y102*c2+Y103*c3) +
// 	     d2b1*(Y110*c0+Y111*c1+Y112*c2+Y113*c3) +
// 	     d2b2*(Y120*c0+Y121*c1+Y122*c2+Y123*c3) +
// 	     d2b3*(Y130*c0+Y131*c1+Y132*c2+Y133*c3))+
// 	    a2 *
// 	    (d2b0*(Y200*c0+Y201*c1+Y202*c2+Y203*c3) +
// 	     d2b1*(Y210*c0+Y211*c1+Y212*c2+Y213*c3) +
// 	     d2b2*(Y220*c0+Y221*c1+Y222*c2+Y223*c3) +
// 	     d2b3*(Y230*c0+Y231*c1+Y232*c2+Y233*c3))+
// 	    a3 *
// 	    (d2b0*(Y300*c0+Y301*c1+Y302*c2+Y303*c3) +
// 	     d2b1*(Y310*c0+Y311*c1+Y312*c2+Y313*c3) +
// 	     d2b2*(Y320*c0+Y321*c1+Y322*c2+Y323*c3) +
// 	     d2b3*(Y330*c0+Y331*c1+Y332*c2+Y333*c3));
	  
// 	  Val += 
// 	    a0*
// 	    (b0*(Y000*d2c0+Y001*d2c1+Y002*d2c2+Y003*d2c3) +
// 	     b1*(Y010*d2c0+Y011*d2c1+Y012*d2c2+Y013*d2c3) +
// 	     b2*(Y020*d2c0+Y021*d2c1+Y022*d2c2+Y023*d2c3) +
// 	     b3*(Y030*d2c0+Y031*d2c1+Y032*d2c2+Y033*d2c3))+
// 	    a1 *
// 	    (b0*(Y100*d2c0+Y101*d2c1+Y102*d2c2+Y103*d2c3) +
// 	     b1*(Y110*d2c0+Y111*d2c1+Y112*d2c2+Y113*d2c3) +
// 	     b2*(Y120*d2c0+Y121*d2c1+Y122*d2c2+Y123*d2c3) +
// 	     b3*(Y130*d2c0+Y131*d2c1+Y132*d2c2+Y133*d2c3))+
// 	    a2 *
// 	    (b0*(Y200*d2c0+Y201*d2c1+Y202*d2c2+Y203*d2c3) +
// 	     b1*(Y210*d2c0+Y211*d2c1+Y212*d2c2+Y213*d2c3) +
// 	     b2*(Y220*d2c0+Y221*d2c1+Y222*d2c2+Y223*d2c3) +
// 	     b3*(Y230*d2c0+Y231*d2c1+Y232*d2c2+Y233*d2c3))+
// 	    a3 *
// 	    (b0*(Y300*d2c0+Y301*d2c1+Y302*d2c2+Y303*d2c3) +
// 	     b1*(Y310*d2c0+Y311*d2c1+Y312*d2c2+Y313*d2c3) +
// 	     b2*(Y320*d2c0+Y321*d2c1+Y322*d2c2+Y323*d2c3) +
// 	     b3*(Y330*d2c0+Y331*d2c1+Y332*d2c2+Y333*d2c3));
	  
// 	  vals(i) += complex<double>(Val[0],Val[1]);
	}
      }
    }
  }
}

/// This replicates the first points of the 3D array to the last
/// points, making the function periodic
void MakePeriodic(Array<complex<double>,4> &A);

#endif
