#include "TricubicBspline.h"
#include "../PlaneWavePHDFT/FFTBox.h"

inline double conj(double x)
{ return x; }

inline double real(double x)
{ return x; }

template<typename T>
double norm (T &func, Vec3 box,
	     int nx, int ny, int nz)
{
  double sixth = 1.0/6.0;
  Array<double,1> wx(nx), wy(ny), wz(nz);
 
  // x weights
  wx(0) = sixth;
  wx(nx-1) = sixth;
  for (int ix=1; ix<nx-1; ix++)
    if ((ix & 1) == 1)
      wx(ix) = 4.0*sixth;
    else
      wx(ix) = 2.0*sixth;

  // y weights
  wy(0) = sixth;
  wy(ny-1) = sixth;
  for (int iy=1; iy<ny-1; iy++)
    if ((iy & 1) == 1)
      wy(iy) = 4.0*sixth;
    else
      wy(iy) = 2.0*sixth;

  // z weights
  wz(0) = sixth;
  wz(nz-1) = sixth;
  for (int iz=1; iz<nz-1; iz++)
    if ((iz & 1) == 1)
      wz(iz) = 4.0*sixth;
    else
      wz(iz) = 2.0*sixth;

  double nrm = 0.0;
  Vec3 r;
  for (int ix=0; ix<(nx); ix++) {
    r[0] = box[0]*(double)ix/(double)(nx-1);
    for (int iy=0; iy<(ny); iy++) {
      r[1] = box[1]*(double)iy/(double)(ny-1);
      for (int iz=0; iz<(nz); iz++) {
	r[2] = box[2]*(double)iz/(double)(nz-1);
	double weight = wx(ix)*wy(iy)*wz(iz);
	nrm += weight*real(conj(func(r))*func(r));
      }
    }
  }
  nrm *= box[0]/(nx-1);
  nrm *= box[1]/(ny-1);
  nrm *= box[2]/(nz-1);
  return nrm;
}



template<typename T1, typename T2>
double InnerProd (T1 &func1, T2 &func2, Vec3 box,
		  int nx, int ny, int nz)
{
  double sixth = 1.0/6.0;
  Array<double,1> wx(nx), wy(ny), wz(nz);
 
  // x weights
  wx(0) = sixth;
  wx(nx-1) = sixth;
  for (int ix=1; ix<nx-1; ix++)
    if ((ix & 1) == 1)
      wx(ix) = 4.0*sixth;
    else
      wx(ix) = 2.0*sixth;

  // y weights
  wy(0) = sixth;
  wy(ny-1) = sixth;
  for (int iy=1; iy<ny-1; iy++)
    if ((iy & 1) == 1)
      wy(iy) = 4.0*sixth;
    else
      wy(iy) = 2.0*sixth;

  // z weights
  wz(0) = sixth;
  wz(nz-1) = sixth;
  for (int iz=1; iz<nz-1; iz++)
    if ((iz & 1) == 1)
      wz(iz) = 4.0*sixth;
    else
      wz(iz) = 2.0*sixth;

  double nrm = 0.0;
  Vec3 r;
  for (int ix=0; ix<(nx); ix++) {
    r[0] = box[0]*(double)ix/(double)(nx-1);
    for (int iy=0; iy<(ny); iy++) {
      r[1] = box[1]*(double)iy/(double)(ny-1);
      for (int iz=0; iz<(nz); iz++) {
	r[2] = box[2]*(double)iz/(double)(nz-1);
	double weight = wx(ix)*wy(iy)*wz(iz);
	nrm += weight*real(conj(func1(r))*func2(r));
      }
    }
  }
  nrm *= box[0]/(nx-1);
  nrm *= box[1]/(ny-1);
  nrm *= box[2]/(nz-1);
  return nrm;
}


void MakeBlipCoefs (FFTBox &fft,  const zVec &c,
		    Array<complex<double>,3> &coefs)
{
  zVec gc(c.size());
  int nx, ny, nz;
  fft.GetDims(nx,ny,nz);
  Vec3 box = fft.GVecs.GetBox();
  Vec3 a;
  a[0] = box[0]/nx;
  a[1] = box[1]/ny;
  a[2] = box[2]/nz;
  
  for (int i=0; i<c.size(); i++) {
    Vec3 G = fft.GVecs(i);
    G[0] *= 1.0*a[0];
    G[1] *= 1.0*a[1];
    G[2] *= 1.0*a[2];
    double gamma = 1.0;
    if (fabs(G[0]) > 1.0e-14)
      gamma *= (3.0/(G[0]*G[0]*G[0]*G[0])*(3.0 - 4.0*cos(G[0]) + cos(2.0*G[0])));
    else
      gamma *= 1.5;
    if (fabs(G[1]) > 1.0e-14)
      gamma *= (3.0/(G[1]*G[1]*G[1]*G[1])*(3.0 - 4.0*cos(G[1]) + cos(2.0*G[1])));
    else
      gamma *= 1.5;
    if (fabs(G[2]) > 1.0e-14)
      gamma *= (3.0/(G[2]*G[2]*G[2]*G[2])*(3.0 - 4.0*cos(G[2]) + cos(2.0*G[2])));
    else
      gamma *= 1.5;

    gc(i) = 3.375*c(i)/gamma;
  }
  fft.PutkVec (gc);
  fft.k2r();
  coefs.resize(nx,ny,nz);
  for (int ix=0; ix<nx; ix++)
    for (int iy=0; iy<ny; iy++)
      for (int iz=0; iz<nz; iz++)
	coefs(ix,iy,iz) = fft.rBox(ix,iy,iz);
}

void
MakeInterpCoefs (FFTBox &fft, const zVec &c,
		 Array<complex<double>,3> &coefs)
{
  fft.kBox = complex<FFT_FLOAT>();
  fft.PutkVec (c);
  fft.k2r();
  Array<complex<double>,3> data;
  int nx, ny, nz;
  fft.GetDims(nx,ny,nz);
  data.resize(nx,ny,nz);
  for (int ix=0; ix<nx; ix++)
    for (int iy=0; iy<ny; iy++)
      for (int iz=0; iz<nz; iz++)
	data(ix,iy,iz) = fft.rBox(ix,iy,iz);
  TricubicBspline<complex<double> > spline;
  Vec3 box = fft.GVecs.GetBox();
  spline.Init (0.0, box[0], 0.0, box[1], 0.0, box[2],
	       data, true, PERIODIC, PERIODIC, PERIODIC);
  coefs.resize(nx,ny,nz);
  for (int ix=0; ix<nx; ix++)
    for (int iy=0; iy<ny; iy++)
      for (int iz=0; iz<nz; iz++)
	coefs(ix,iy,iz) = spline.GetControlPoint (ix,iy,iz);

}

complex<double> eval(double x, double y, double z, 
		     zVec &coefs, GVecsClass &gvecs)
{
  complex<double> val = complex<double>();
  for (int i=0; i<coefs.size(); i++) {
    Vec3 G = gvecs(i);
    double phase = -(x*G[0]+ y *G[1] + z*G[2]);
    val += coefs(i)*complex<double> (cos(phase), sin(phase));
  }
  return val;
}



class PeriodicFunc
{
private:
  double NormFactor;

public:
  Array<complex<double>,1> c;

  GVecsClass &GVecs;
  FFTBox FFT;
  complex<double> operator()(TinyVector<double,3> r)
  {
    complex<double> val = complex<double>();
    for (int i=0; i<GVecs.size(); i++) {
      double phase = -dot(r, GVecs(i));
      complex<double> e2iGr = complex<double>(cos(phase), sin(phase));
      val += c(i) * e2iGr;
    }
    return val;
  }

//   void Norms (double &valExp, double gradExp, double &laplExp)
//   {
//     valExp = gradExp = laplExp = 0.0;
//     for (int i=0; i<c.size(); i++) {
//       TinyVector<double,3> G = GVecs(i);
//       valExp  += norm(c(i));
//       gradExp += norm (G[0]*c(i)) + norm(G[1]*c(i)) + norm(G[2]*c(i));
//       laplExp += dot (G,G)*norm(c(i));
//     }
//   }

  void Evaluate (TinyVector<double,3> r, complex<double> &val, 
		 TinyVector<complex<double>,3> &grad, complex<double> &laplacian)
  {
    val = (*this)(r);
    grad = complex<double>();
    laplacian = complex<double>();
    complex<double> eye(0.0, 1.0);
    for (int i=0; i<c.size(); i++) {
      Vec3 G = GVecs(i);
      double phase = -dot(r,G);
      complex<double> e2iGr = complex<double>(cos(phase), sin(phase));
      grad += -(eye * c(i)*G*e2iGr);
      laplacian += -dot(G,G)*c(i)*e2iGr;
    }
  }
  void Evaluate (TinyVector<double,3> r, 
		 complex<double> &val, TinyVector<complex<double>,3> &grad,
		 TinyMatrix<complex<double>,3,3> &secDerivs)
  {
    val = (*this)(r);
    grad = complex<double>();
    secDerivs = complex<double>();
    complex<double> eye(0.0, 1.0);
    for (int i=0; i<c.size(); i++) {
      Vec3 G = GVecs(i);
      double phase = -dot(r,G);
      complex<double> e2iGr = complex<double>(cos(phase), sin(phase));
      grad += NormFactor*-(eye * c(i)*G*e2iGr);
      secDerivs(0,0) +=  -G[0]*G[0]*c(i)*e2iGr;
      secDerivs(0,1) +=  -G[0]*G[1]*c(i)*e2iGr;
      secDerivs(0,2) +=  -G[0]*G[2]*c(i)*e2iGr;
      secDerivs(1,0) +=  -G[1]*G[0]*c(i)*e2iGr;
      secDerivs(1,1) +=  -G[1]*G[1]*c(i)*e2iGr;
      secDerivs(1,2) +=  -G[1]*G[2]*c(i)*e2iGr;
      secDerivs(2,0) +=  -G[2]*G[0]*c(i)*e2iGr;
      secDerivs(2,1) +=  -G[2]*G[1]*c(i)*e2iGr;
      secDerivs(2,2) +=  -G[2]*G[2]*c(i)*e2iGr;
    }
  }

  PeriodicFunc(GVecsClass &gvecs) : GVecs(gvecs), FFT(gvecs)
  {
    FFT.Setup();
    c.resize(gvecs.size());
    int nx, ny, nz;
    FFT.GetDims(nx, ny, nz);
    NormFactor = sqrt(1.0/(gvecs.GetBoxVol()));    
    double nrm;
    for (int i=0; i<gvecs.size(); i++) {
      c(i) = complex<double>(2.0*drand48()-1.0, 2.0*drand48()-1);
      nrm += norm (c(i));
    }
    c /= sqrt(nrm);
  }
};



class LocalFunc
{
private:
  TinyVector<double,3> Box, BoxInv;
  Array<TinyVector<double,3>,1> Center;
  Array<double,1> Prefactor, Alpha;
  double a;
  
  inline double clamp(double t)
  {
    double u = 1.0-t;
    double e1 = exp(-a*a*a*t*t*t);
    double e2 = exp(-a*a*a*u*u*u);
    return (1.0-e1)*(1.0-e2);
  }
  inline double dclamp(double t)
  {
   double u = 1.0-t;
   double e1 = exp(-a*a*a*t*t*t);
   double e2 = exp(-a*a*a*u*u*u);

   return (+3.0*a*a*a*t*t*e1*(1.0-e2) +
	   -3.0*a*a*a*u*u*e2*(1.0-e1));
  }
  inline double d2clamp(double t)
  {
    double u = 1.0-t;
    double a3 = a*a*a;
    double e1 = exp(-a3*t*t*t);
    double e2 = exp(-a3*u*u*u);
    double de1 = -3.0*a3*t*t*e1;
    double de2 = +3.0*a3*u*u*e2;
    double d2e1 = (-6.0*a3*t + 3.0*a3*a3*t*t*t*t)*e1;
    double d2e2 = (-6.0*a3*u + 3.0*a3*a3*u*u*u*u)*e2;
    return (d2e1*(e2-1.0)+d2e2*(e1-1.0) + 2.0*de1*de2);
//     return (+6.0*a3*t*e1*(1.0-e2) 
// 	    -9.0*a3*a3*t*t*t*t*e1*(1.0-e2)
// 	    -18.0*a3*a3*t*t*u*u*e1*e2
// 	    -6.0*a3*u*e2*(1.0-e1)
// 	    +9.0*a3*a3*u*u*u*u*e2*(1.0-e1));
  }


public:
  double operator() (TinyVector<double,3> r)
  {
    double val = 0.0;
    r[0]*=BoxInv[0];  r[1]*=BoxInv[1];  r[2]*=BoxInv[2];
    for (int i=0; i<Center.size(); i++) {
      TinyVector<double,3> dr = r - Center(i);
      val += Prefactor(i)*exp(-Alpha(i)*dot(dr,dr));
    }
    val *= clamp(r[0])*clamp(r[1])*clamp(r[2]);
    return val;
  }
  void Evaluate (TinyVector<double,3> r, double &val, 
		 TinyVector<double,3> &grad, double &laplacian)
  {
    TinyVector<double,3> t = r;
    t[0]*=BoxInv[0];  t[1]*=BoxInv[1];  t[2]*=BoxInv[2];
    val = 0.0;
    grad = 0.0;
    laplacian = 0.0;
    for (int i=0; i<Center.size(); i++) {
      TinyVector<double,3> dt = t - Center(i);
      TinyVector<double,3> dr = r - TinyVector<double,3>
	(Center(i)[0]*Box[0], Center(i)[1]*Box[1], Center(i)[2]*Box[2]);
      double expFact = exp(-Alpha(i)*dot(dt,dt));
      double expx = exp(-Alpha(i)*dt[0]*dt[0]);
      double expy = exp(-Alpha(i)*dt[1]*dt[1]);
      double expz = exp(-Alpha(i)*dt[2]*dt[2]);

      double phix = expx*clamp(t[0]);
      double phiy = expy*clamp(t[1]);
      double phiz = expz*clamp(t[2]);

      val += Prefactor(i) * phix*phiy*phiz;

      double dexpx = -2.0*Alpha(i)*dt[0]*BoxInv[0]*expx;
      double dexpy = -2.0*Alpha(i)*dt[1]*BoxInv[1]*expy;
      double dexpz = -2.0*Alpha(i)*dt[2]*BoxInv[2]*expz;
      
//       double dphix = -2.0*Alpha(i)*dt[0]*BoxInv[0]*phix +
// 	BoxInv[0]*dclamp(t[0])*exp(-Alpha(i)*dt[0]*dt[0]);
//       double dphiy = -2.0*Alpha(i)*dt[1]*BoxInv[1]*phiy +
// 	BoxInv[1]*dclamp(t[1])*exp(-Alpha(i)*dt[1]*dt[1]);
//       double dphix = -2.0*Alpha(i)*dt[2]*BoxInv[2]*phiy +
// 	BoxInv[2]*dclamp(t[2])*exp(-Alpha(i)*dt[2]*dt[2]);

      double dphix = BoxInv[0]*dclamp(t[0])*expx + clamp(t[0])*dexpx;
      double dphiy = BoxInv[1]*dclamp(t[1])*expy + clamp(t[1])*dexpy;
      double dphiz = BoxInv[2]*dclamp(t[2])*expz + clamp(t[2])*dexpz;

      grad[0] += Prefactor(i) *(dphix *  phiy *  phiz);
      grad[1] += Prefactor(i) *( phix * dphiy *  phiz);
      grad[2] += Prefactor(i) *( phix *  phiy * dphiz);

      double d2expx = expx*Alpha(i)*BoxInv[0]*BoxInv[0]*
	(-2.0 + 4.0*Alpha(i)*dt[0]*dt[0]);
      double d2expy = expy*Alpha(i)*BoxInv[1]*BoxInv[1]*
	(-2.0 + 4.0*Alpha(i)*dt[1]*dt[1]);
      double d2expz = expz*Alpha(i)*BoxInv[2]*BoxInv[2]*
	(-2.0 + 4.0*Alpha(i)*dt[2]*dt[2]);
	
      double d2phix = BoxInv[0]*BoxInv[0]*d2clamp(t[0])*expx + clamp(t[0])*d2expx
	+2.0*BoxInv[0]*dclamp(t[0])*dexpx;
      double d2phiy = BoxInv[1]*BoxInv[1]*d2clamp(t[1])*expy + clamp(t[1])*d2expy
	+2.0*BoxInv[1]*dclamp(t[1])*dexpy;
      double d2phiz = BoxInv[2]*BoxInv[2]*d2clamp(t[2])*expz + clamp(t[2])*d2expz
	+2.0*BoxInv[2]*dclamp(t[2])*dexpz;

      laplacian += Prefactor(i)*(d2phix *   phiy *   phiz+
				 phix   * d2phiy *   phiz +
				 phix   *   phiy * d2phiz  );


//       val += Prefactor(i)*expFact;
//       grad -= 2.0*Prefactor(i)*Alpha(i)*expFact*TinyVector<double,3>
// 	(dt[0]*BoxInv[0], dt[1]*BoxInv[1], dt[2]*BoxInv[2]);
//       laplacian += Prefactor(i)*Alpha(i)*expFact*
// 	(-2.0*dot(BoxInv,BoxInv) + 4.0*Alpha(i)*
// 	 (BoxInv[0]*BoxInv[0]*dt[0]*dt[0] +
// 	  BoxInv[1]*BoxInv[1]*dt[1]*dt[1] +
// 	  BoxInv[2]*BoxInv[2]*dt[2]*dt[2]));
    }
  }	 
					  

  

  LocalFunc (TinyVector<double,3> box, int num)
  {
    Box = box;
    a = 100.0;
    BoxInv[0] = 1.0/Box[0];
    BoxInv[1] = 1.0/Box[1];
    BoxInv[2] = 1.0/Box[2];
    Center.resize(num);
    Prefactor.resize(num);
    Alpha.resize(num);
    for (int i=0; i<num; i++) {
      Center(i)[0] = 0.1+0.8*drand48();
      Center(i)[1] = 0.1+0.8*drand48();
      Center(i)[2] = 0.1+0.8*drand48();
      Prefactor(i) = 2.0*drand48()-1.0;
      Alpha(i) = 20.0 + 1000.0*drand48();
    }
    double nrm = norm (*this, box, 200, 200, 200);
    Prefactor *= (1.0/sqrt(nrm));
  }
};


void
MakeBlipCoefs (Array<double,3> &data,
	       Array<double,3> &coefs)
{
  int nx = data.extent(0);
  int ny = data.extent(1);
  int nz = data.extent(2);
  int nx2 = (nx-1)/2;
  int ny2 = (ny-1)/2;
  int nz2 = (nz-1)/2;
  FFT3D fft;
  fft.resize(nx, ny, nz);
  for (int ix=0; ix<nx; ix++)
    for (int iy=0; iy<ny; iy++)
      for (int iz=0; iz<nz; iz++)
	fft.rBox(ix,iy,iz) = fft_complex(data(ix,iy,iz));
  fft.r2k();

  double fracx = 2.0*M_PI/(double)nx;
  double fracy = 2.0*M_PI/(double)ny;
  double fracz = 2.0*M_PI/(double)nz;
  // Now we have the fourier transform -- multiply coefs apprpriately
  double gammax, gammay, gammaz;
  for (int ix=0; ix<nx; ix++) {
    double fx = fracx*(double)((ix+nx2)%nx - nx2);
    if (fabs(fx) > 1.0e-12)
      gammax = (3.0/(fx*fx*fx*fx)*(3.0 - 4.0*cos(fx) + cos(2.0*fx)))/1.5;
    else
      gammax = 1.0;
    for (int iy=0; iy<ny; iy++) {
      double fy = fracy*(double)((iy+ny2)%ny - ny2);
      if (fabs(fy) > 1.0e-12)
	gammay = (3.0/(fy*fy*fy*fy)*(3.0 - 4.0*cos(fy) + cos(2.0*fy)))/1.5;
      else
	gammay = 1.0;
      for (int iz=0; iz<nz; iz++) {
	double fz = fracz*(double)((iz+nz2)%nz - nz2);
	if (fabs(fz) > 1.0e-12)
	  gammaz = (3.0/(fz*fz*fz*fz)*(3.0 - 4.0*cos(fz) + cos(2.0*fz)))/1.5;
	else
	  gammaz = 1.0;

	fft.kBox(ix,iy,iz) *= (1.0/(gammax*gammay*gammaz));
      }
    }
  }
  fft.k2r();
  for (int ix=0; ix<nx; ix++)  
    for (int iy=0; iy<ny; iy++)
      for (int iz=0; iz<nz; iz++)
	coefs(ix,iy,iz) = fft.kBox(ix,iy,iz).real();
}

inline double sqr(double x)
{
  return x*x;
}




void
TestLocal(LocalFunc &local, Vec3 &box, 
	  int nx, int ny, int nz)
{
  Vec3 r;

  Array<double,3> spdata(nx,ny,nz);
  for (int ix=0; ix<nx; ix++) {
    r[0] = box[0]*(double)ix/(double)(nx-1);
    for (int iy=0; iy<ny; iy++) {
      r[1] = box[1]*(double)iy/(double)(ny-1);
      for (int iz=0; iz<nz; iz++) {
	r[2] = box[2]*(double)iz/(double)(nz-1);
	spdata(ix,iy,iz) = local(r);
      }
    }
  }
  Array<double,3> blipdata(nx-1, ny-1, nz-1), 
    blipcoefs (nx-1, ny-1, nz-1);
  blipdata = spdata(Range(0,nx-2),Range(0,ny-2), Range(0,nz-2));
  MakeBlipCoefs(blipdata, blipcoefs);

  TricubicBspline<double> interpSpline, blipSpline;
  //  spline.Init (0.0, box[0], 0.0, box[1], 0.0, box[2], 
  //	       spdata, true, PERIODIC, PERIODIC, PERIODIC);
  interpSpline.Init (0.0, box[0], 0.0, box[1], 0.0, box[2], 
		     spdata, true, FLAT, FLAT, FLAT);
  //   interpSpline.Init (0.0, box[0], 0.0, box[1], 0.0, box[2], 
  // 		     spdata, true, NATURAL, NATURAL, NATURAL);
  blipSpline.Init (0.0, box[0], 0.0, box[1], 0.0, box[2], 
		   blipcoefs, false, PERIODIC, PERIODIC, PERIODIC);


  // Compute norms
//   double blipNorm=0.0, spNorm=0.0, exNorm=0.0;
//   double sixth = 1.0/6.0;
//   Array<double,1> wx(2*nx), wy(2*ny), wz(2*nz);
 
//   // x weights
//   wx(0) = sixth;
//   wx(2*nx-1) = sixth;
//   for (int ix=1; ix<2*nx-1; ix++)
//     if ((ix & 1) == 1)
//       wx(ix) = 4.0*sixth;
//     else
//       wx(ix) = 2.0*sixth;

//   // y weights
//   wy(0) = sixth;
//   wy(2*ny-1) = sixth;
//   for (int iy=1; iy<2*ny-1; iy++)
//     if ((iy & 1) == 1)
//       wy(iy) = 4.0*sixth;
//     else
//       wy(iy) = 2.0*sixth;

//   // z weights
//   wz(0) = sixth;
//   wz(2*nz-1) = sixth;
//   for (int iz=1; iz<2*nz-1; iz++)
//     if ((iz & 1) == 1)
//       wz(iz) = 4.0*sixth;
//     else
//       wz(iz) = 2.0*sixth;

//   for (int ix=0; ix<(2*nx); ix++) {
//     r[0] = box[0]*(double)ix/(double)(2*nx-1);
//     for (int iy=0; iy<(2*ny); iy++) {
//       r[1] = box[1]*(double)iy/(double)(2*ny-1);
//       for (int iz=0; iz<(2*nz); iz++) {
// 	r[2] = box[2]*(double)iz/(double)(2*nz-1);
// 	double weight = wx(ix)*wy(iy)*wz(iz);
// 	exNorm += weight*local(r)*local(r);
// 	spNorm += weight*interpSpline(r)*interpSpline(r);
// 	blipNorm += weight * blipSpline(r)*blipSpline(r);
//       }
//     }
//   }

  double spNorm   = norm (interpSpline, box, 2*nx+1, 2*ny+1, 2*nz+1);
  double blipNorm = norm (blipSpline,   box, 2*nx+1, 2*ny+1, 2*nz+1);

  double interpEx = InnerProd (local, interpSpline, 
			       box, 2*nx+1, 2*ny+1, 2*nz+1);
  double blipEx   = InnerProd (local, blipSpline, 
			       box, 2*nx+1, 2*ny+1, 2*nz+1);
  double interpAlpha = 1.0 - interpEx / sqrt(spNorm);
  double blipAlpha  =  1.0 -  blipEx / sqrt(blipNorm);
  cerr << "interp alpha = " << interpAlpha  << endl;
  cerr << "blip   alpha = " << blipAlpha << endl;

		 
  // renormalize
  spdata    *= (1.0/sqrt(spNorm));
  blipcoefs *= (1.0/sqrt(blipNorm));
  interpSpline.Init (0.0, box[0], 0.0, box[1], 0.0, box[2], 
		     spdata, true, 
		     FLAT, FLAT, FLAT);
  blipSpline.Init (0.0, box[0], 0.0, box[1], 0.0, box[2], 
		   blipcoefs, false, 
		   PERIODIC, PERIODIC, PERIODIC);
  fprintf (stderr, "%20.16f %20.16f\n", 
	   spNorm, blipNorm);
//   exNorm   = norm (local,        box, 2*nx+1, 2*ny+1, 2*nz+1);
//   spNorm   = norm (interpSpline, box, 2*nx+1, 2*ny+1, 2*nz+1);
//   blipNorm = norm (blipSpline,   box, 2*nx+1, 2*ny+1, 2*nz+1);
//   fprintf (stderr, "%20.16f %20.16f %20.16f\n", 
// 	   exNorm, spNorm, blipNorm);



  // Compare errors
  double bValError =0.0, sValError =0.0;
  double bGradError=0.0, sGradError=0.0;
  double bLaplError=0.0, sLaplError=0.0;
  double xValNorm = 0.0, xGradNorm=0.0, xLaplNorm=0.0;
  for (int i=0; i<10000; i++) {
    TinyVector<double,3> xGrad, bGrad, sGrad;
    double xLapl, bLapl, sLapl;
    double xVal, bVal, sVal;
    TinyVector<double,3> r (box[0]*drand48(),
			    box[1]*drand48(),
			    box[2]*drand48());
			      
    local.Evaluate (r, xVal, xGrad, xLapl);
    blipSpline.Evaluate (r, bVal, bGrad, bLapl);
    interpSpline.Evaluate (r, sVal, sGrad, sLapl);
    sValError  += sqr(sVal - xVal);
    bValError  += sqr(bVal - xVal);
    sGradError += sqr(sGrad[0]-xGrad[0]) + sqr(sGrad[1]-xGrad[1]) + sqr(sGrad[2]-xGrad[2]);
    bGradError += sqr(bGrad[0]-xGrad[0]) + sqr(bGrad[1]-xGrad[1]) + sqr(bGrad[2]-xGrad[2]);
    sLaplError += sqr(sLapl-xLapl);
    bLaplError += sqr(bLapl-xLapl);
    xValNorm   += sqr(xVal);
    xGradNorm  += sqr(xGrad[0])+sqr(xGrad[1]) + sqr(xGrad[2]);
    xLaplNorm  += sqr (xLapl);
  } 
  sValError = sqrt(sValError/xValNorm);
  bValError = sqrt(bValError/xValNorm);
  sGradError = sqrt(sGradError/xGradNorm);
  bGradError = sqrt(bGradError/xGradNorm);
  sLaplError = sqrt(sLaplError/xLaplNorm);
  bLaplError = sqrt(bLaplError/xLaplNorm);
  double spacing = box[0] / blipcoefs.extent(0);
  fprintf (stdout, "%20.16e %20.16e %20.16e %20.16e %20.16e %20.16e %20.16e %20.16e %20.16e\n", 
	   spacing, sValError, bValError,
	   sGradError, bGradError,
	   sLaplError, bLaplError,
	   interpAlpha, blipAlpha);
  fflush (stdout);


//   r[1] = 2.8;  r[2] = 5.7;
//   for (r[0]=0.0; r[0]<=box[0]; r[0]+=0.001)
//     fprintf (stdout, "%20.16e %20.16e %20.16e %20.16e\n", 
// 	     r[0], local(r), interpSpline(r), blipSpline(r));
}


void TestLocal()
{
  Vec3 box(10.0, 10.0, 10.0);

  LocalFunc local (box, 100);
  for (int n=21; n<300; n+=20)
    TestLocal (local, box, n,n,n);
}



void
TestOnePeriodic(PeriodicFunc &pfunc)
{
  Vec3 box = pfunc.GVecs.GetBox();
  int nx, ny, nz;
  pfunc.FFT.GetDims(nx, ny, nz);

  // Make blip function
  Array<complex<double>,3> blipCoefs;
  MakeBlipCoefs (pfunc.FFT, pfunc.c, blipCoefs);
  TricubicBspline<complex<double> > blipSpline;
  blipSpline.Init (0.0, box[0], 0.0, box[1], 0.0, box[2],
		   blipCoefs, false, 
		   PERIODIC, PERIODIC, PERIODIC);

  // Make interpolating B-spline
  Array<complex<double>,3> interpCoefs;
  MakeInterpCoefs (pfunc.FFT, pfunc.c, interpCoefs);
  TricubicBspline<complex<double> > interpSpline;
  interpSpline.Init (0.0, box[0], 0.0, box[1], 0.0, box[2],
		     interpCoefs, false, 
		     PERIODIC, PERIODIC, PERIODIC);
  
//   double interpEx   = InnerProd (pfunc, interpSpline, box,
// 			       2*nx+1, 2*ny+1, 2*nz+1);
//   double blipEx     = InnerProd (pfunc, blipSpline, box,
// 			        2*nx+1, 2*ny+1, 2*nz+1);
//   double interpNorm = norm (interpSpline, box,
// 			    2*nx+1, 2*ny+1, 2*nz+1);
//   interpCoefs *= (1.0/sqrt(interpNorm));
//   interpSpline.Init (0.0, box[0], 0.0, box[1], 0.0, box[2],
// 		     interpCoefs, false, 
// 		     PERIODIC, PERIODIC, PERIODIC);

//   double blipNorm   = norm (blipSpline, box,
// 			    2*nx+1, 2*ny+1, 2*nz+1);
//   blipCoefs *= (1.0/sqrt(blipNorm));
//   blipSpline.Init (0.0, box[0], 0.0, box[1], 0.0, box[2],
// 		   blipCoefs, false, 
// 		   PERIODIC, PERIODIC, PERIODIC);
//   double exNorm     = norm (pfunc, box, 2*nx+1, 2*nx+1, 2*nx+1);

//   double sAlpha = 1.0 - interpEx/sqrt(interpNorm*exNorm);
//   double bAlpha = 1.0 - blipEx/sqrt(blipNorm*exNorm);

  // Compare errors
  double bValError =0.0, sValError =0.0;
  double bGradError=0.0, sGradError=0.0;
  double bLaplError=0.0, sLaplError=0.0;
  double xValNorm  =0.0, xGradNorm=0.0, xLaplNorm=0.0;
  double sAlpha    =1.0, bAlpha = 1.0;
  double SpSpVal   =0.0, SpExVal = 0.0;
  double BlBlVal   =0.0, BlExVal = 0.0;
  double ExExVal   =0.0;
  for (int i=0; i<10000; i++) {
    TinyVector<complex<double>,3> xGrad, bGrad, sGrad;
    complex<double> xLapl, bLapl, sLapl;
    complex<double> xVal, bVal, sVal;
    TinyVector<double,3> r (box[0]*drand48(),
			    box[1]*drand48(),
			    box[2]*drand48());
			      
    pfunc.Evaluate (r, xVal, xGrad, xLapl);
    blipSpline.Evaluate (r, bVal, bGrad, bLapl);
    interpSpline.Evaluate (r, sVal, sGrad, sLapl);
    sValError  += norm(sVal - xVal);
    bValError  += norm(bVal - xVal);
    sGradError += norm(sGrad[0]-xGrad[0]) + norm(sGrad[1]-xGrad[1]) + norm(sGrad[2]-xGrad[2]);
    bGradError += norm(bGrad[0]-xGrad[0]) + norm(bGrad[1]-xGrad[1]) + norm(bGrad[2]-xGrad[2]);
    sLaplError += norm (sLapl-xLapl);
    bLaplError += norm (bLapl-xLapl);
    xValNorm   += norm (xVal);
    xGradNorm  += norm(xGrad[0])+norm(xGrad[1]) + norm(xGrad[2]);
    xLaplNorm  += norm (xLapl);
    ExExVal += norm (xVal);
    SpSpVal += norm (sVal);
    BlBlVal += norm (bVal);
    SpExVal += real(conj(sVal)*xVal);
    BlExVal += real(conj(bVal)*xVal);
  } 
  sAlpha = 1.0-SpExVal / sqrt(SpSpVal*ExExVal);
  bAlpha = 1.0-BlExVal / sqrt(BlBlVal*ExExVal);
  sValError = sqrt(sValError/xValNorm);
  bValError = sqrt(bValError/xValNorm);
  sGradError = sqrt(sGradError/xGradNorm);
  bGradError = sqrt(bGradError/xGradNorm);
  sLaplError = sqrt(sLaplError/xLaplNorm);
  bLaplError = sqrt(bLaplError/xLaplNorm);
  double spacing = box[0] / interpCoefs.extent(0);
  fprintf (stdout, "%20.16e %20.16e %20.16e %20.16e %20.16e %20.16e %20.16e %20.16e %20.16e\n", 
	   spacing, sValError, bValError,
	   sGradError, bGradError,
	   sLaplError, bLaplError,
	   sAlpha, bAlpha);
  fflush (stdout);
}

void
TestPeriodic()
{
  Array<complex<double>,1> coefs;
  const double kcut = 3.0*M_PI;
  Vec3 box (10.0, 11.0, 12.0);
  GVecsClass GVecs;
  Vec3 k (0.0, 0.0, 0.0);
  GVecs.Set(box, k, kcut, 2.0);
  PeriodicFunc pfunc1(GVecs);
  coefs.resize(pfunc1.c.size());
  coefs = pfunc1.c;
  for (double factor=1.0; factor<=4.0; factor+=0.2) {
    GVecs.Set(box, k, kcut, factor);
    PeriodicFunc pfunc(GVecs);
    pfunc.c = coefs;
    TestOnePeriodic(pfunc);
  }

}




void
TestBlips()
{
  const double kcut = 3.0*M_PI;
  Vec3 box (10.0, 10.0, 10.0);
  GVecsClass GVecs;
  Vec3 k (0.0, 0.0, 0.0);
  
  GVecs.Set(box, k, kcut, 4.0);
  FFTBox fft(GVecs);
  fft.Setup();
  zVec c(GVecs.size());
  for (int i=0; i<GVecs.size(); i++)
    c(i) = complex<double> (-1.0+2.0*drand48(),-1.0+2.0*drand48());

  Array<complex<double>,3> blipCoefs, interpCoefs;
  int nx, ny, nz;
  fft.GetDims(nx, ny, nz);
  blipCoefs.resize(nx, ny, nz);
  interpCoefs.resize(nx, ny, nz);
  MakeBlipCoefs(fft, c, blipCoefs);
   MakeInterpCoefs(fft, c, interpCoefs);
//   fft.PutkVec(c);
//   fft.k2r();
//   interpCoefs = fft.rBox;
  //blipCoefs   = fft.rBox;
//   for (int ix=0; ix<interpCoefs.extent(0); ix++)
//     for (int iy=0; iy<interpCoefs.extent(1); iy++)
//       for (int iz=0; iz<interpCoefs.extent(2); iz++) {
// 	double x = box[0]*(double)ix/interpCoefs.extent(0);
// 	double y = box[1]*(double)iy/interpCoefs.extent(1);
// 	double z = box[2]*(double)iz/interpCoefs.extent(2);
// 	interpCoefs(ix,iy,iz) = eval(x,y,z,c,GVecs);
//       }
  TricubicBspline<complex<double> > blipSpline, interpSpline;
  blipSpline.Init (0.0, box[0], 0.0, box[1], 0.0, box[2],
		   blipCoefs, false, PERIODIC, PERIODIC, PERIODIC);
  interpSpline.Init (0.0, box[1], 0.0, box[1], 0.0, box[2],
		     interpCoefs, false, PERIODIC, PERIODIC, PERIODIC);
//   for (int ix=0; ix<blipCoefs.extent(0); ix++)
//     for (int iy=0; iy<blipCoefs.extent(1); iy++)
//       for (int iz=0; iz<blipCoefs.extent(2); iz++)
// 	fprintf (stdout, "%22.16e %22.16e\n", blipCoefs(ix,iy,iz).real(),
// 		 interpCoefs(ix,iy,iz).real());
  double y= 2.3;
  double z = 3.8;
  for (double x=0.0; x<=10.0; x+=0.001)
    fprintf (stdout, "%22.16e %22.16e %22.16e %22.16e\n", x, interpSpline(x,y,z).real(),
	     blipSpline(x,y,z).real(), eval(x,y,z,c,GVecs).real());

  double blValError =0.0, spValError =0.0;
  double blGradError=0.0, spGradError=0.0;
  double blLaplError=0.0, spLaplError=0.0;
  for (int i=0; i<10000; i++) {
    double x = drand48()*box[0];
    double y = drand48()*box[1];
    double z = drand48()*box[2];
  }

	     
}


main()
{
  TestPeriodic();
  //TestLocal();
}
