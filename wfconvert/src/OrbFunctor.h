#ifndef ORB_FUNCTOR_H
#define ORB_FUNCTOR_H

#include <Common/PlaneWavePHDFT/GVecs.h>


class OrbFunctor
{
public:
  virtual complex<double> operator()(Vec3 r) = 0;
  virtual void operator()(Vec3 r, complex<double> &val,
			  TinyVector<complex<double>,3> &grad,
			  complex<double> &lapl) = 0;
  virtual void operator()(vector<Vec3> &r,
			  vector<complex<double> > &val)
  {
    for (int ir=0; ir<r.size(); ir++)
      val[ir] = (*this)(r[ir]);
  }

  virtual void operator()(vector<Vec3> &r,
			  vector<complex<double> > &val,
			  vector<complex<double> > &lapl)
  {
    TinyVector<complex<double>,3> grad;
    for (int ir=0; ir<r.size(); ir++)
      (*this)(r[ir], val[ir], grad, lapl[ir]);
  }

  virtual void operator()(vector<Vec3> &r,
			  vector<complex<double> > &val,
			  vector<TinyVector<complex<double>,3> > &grad,
			  vector<complex<double> > &lapl)
  {
    for (int ir=0; ir<r.size(); ir++)
      (*this)(r[ir], val[ir], grad[ir], lapl[ir]);
  }
};


#include "Common/config.h"
#include "PlaneWaveSum.h"

#ifdef HAVE_CUDA
  #include <cuda_runtime_api.h>
#endif


class PWOrb_z : public OrbFunctor
{
public:
  zVec &Coefs;
  GVecsClass &GVecs;
  complex<double> Prefactor;
  Vec3 k;
#ifdef HAVE_CUDA
  double *GVecs_cuda, *r_cuda;
  complex<double> *Coefs_cuda, *val_cuda, *grad_cuda,
    *lapl_cuda;
  int num_r;
#endif
public:
  complex<double> operator()(Vec3 r)
  {
    complex<double> sum;
    for (int iG=0; iG<GVecs.size(); iG++) {
      double phase = -dot(GVecs(iG)+k,r);
      double s,c;
      sincos(phase, &s, &c);
      sum += Coefs(iG)*complex<double>(c,s);
    }
    
    return Prefactor*sum;
  }

#ifndef HAVE_CUDA
  void operator()(vector<Vec3> &r,
                  vector<complex<double> > &val) 
  {
    for (int i=0; i<val.size(); i++)
      val[i] = complex<double>();
    for (int iG=0; iG<GVecs.size(); iG++) {
      for (int ir=0; ir<r.size(); ir++) {
	double phase = -dot(GVecs(iG)+k,r[ir]);
	double s,c;
	sincos(phase, &s, &c);
	val[ir] += Coefs(iG)*complex<double>(c,s);
      }
    }
  }
#endif

  void operator()(Vec3 r, complex<double> &val,
		  TinyVector<complex<double>,3> &grad,
		  complex<double> &lapl)
  {
    val = lapl = complex<double>();
    grad = TinyVector<complex<double>,3>();
    for (int iG=0; iG<GVecs.size(); iG++) {
      double phase = -dot(GVecs(iG)+k,r);
      double s,c;
      sincos(phase, &s, &c);
      complex<double> z(c,s);
      val +=  Coefs(iG)*z;
      grad -= Coefs(iG)*complex<double>(-s,c)*(GVecs(iG)+k);
      lapl -= Coefs(iG)*dot(GVecs(iG)+k,GVecs(iG)+k)*z;
    }
    val  *= Prefactor;
    grad *= Prefactor;
    lapl *= Prefactor;
  }
  

#ifdef HAVE_CUDA
  void alloc_cuda(int num)
  {
    if (val_cuda) {
      cudaFree (r_cuda);
      cudaFree (val_cuda);
      cudaFree (grad_cuda);
      cudaFree (lapl_cuda);
    }
    num_r = num;
    cudaMalloc((void**)&r_cuda,    3*num*sizeof(double));
    cudaMalloc((void**)&val_cuda,    num*sizeof(complex<double>));
    cudaMalloc((void**)&grad_cuda, 3*num*sizeof(complex<double>));
    cudaMalloc((void**)&lapl_cuda,   num*sizeof(complex<double>));
  }


  void operator()(vector<Vec3> &r,
		  vector<complex<double> > &val)
  {
    if (num_r < r.size())  alloc_cuda(r.size());
    cudaMemcpy (r_cuda, &r[0][0], 3*r.size()*sizeof(double),
    		cudaMemcpyHostToDevice);
    cudaThreadSynchronize();
    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess) {
      fprintf (stderr, "CUDA error in memcopy:\n  %s\n",
	       cudaGetErrorString(err));
      abort();
    }
    plane_wave_sum (r_cuda, GVecs_cuda, k[0], k[1], k[2],
		    Coefs_cuda,
    		    val_cuda, r.size(), GVecs.size());
    cudaMemcpy (&val[0], val_cuda, val.size()*sizeof(complex<double>),
    		cudaMemcpyDeviceToHost);
    for (int i=0; i<val.size(); i++) 
      val[i] *= Prefactor;

    // vector<complex<double> > gpu_val = val;
    // OrbFunctor::operator()(r,val);
    // for (int i=0; i<val.size(); i++)
    //   fprintf (stderr, "%12.5e %12.5ei  %12.5e %12.5ei\n",
    // 	       val[i].real(), val[i].imag(), 
    // 	       gpu_val[i].real(), gpu_val[i].imag());
  
  }

  void operator()(vector<Vec3> &r,
		  vector<complex<double> > &val,
		  vector<complex<double> > &lapl)
  {
    if (num_r != r.size())  alloc_cuda(r.size());
    cudaMemcpy (r_cuda, &r[0][0], 3*r.size()*sizeof(double),
    		cudaMemcpyHostToDevice);
    cudaThreadSynchronize();
    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess) {
      fprintf (stderr, "CUDA error in memcopy:\n  %s\n",
	       cudaGetErrorString(err));
      abort();
    }
    plane_wave_sum (r_cuda, GVecs_cuda, k[0], k[1], k[2],
		    Coefs_cuda,
    		    val_cuda, lapl_cuda, r.size(), GVecs.size());
    cudaMemcpy (&val[0], val_cuda, val.size()*sizeof(complex<double>),
    		cudaMemcpyDeviceToHost);
    cudaMemcpy (&lapl[0], lapl_cuda, val.size()*sizeof(complex<double>),
    		cudaMemcpyDeviceToHost);
    for (int i=0; i<val.size(); i++) {
      val[i]  *= Prefactor;
      lapl[i] *= Prefactor;
    }
  }

#else
  void operator()(vector<Vec3> &r,
		  vector<complex<double> > &val,
		  vector<complex<double> > &lapl)
  {
    TinyVector<complex<double>,3> grad;
    for (int ir=0; ir<r.size(); ir++)
      (*this)(r[ir], val[ir], grad, lapl[ir]);
  }
#endif



  PWOrb_z (zVec &coefs, GVecsClass &gvecs, 
	   Vec3 kpoint, complex<double> pref) :
    Coefs(coefs), GVecs(gvecs), Prefactor(pref), k(kpoint)
#ifdef HAVE_CUDA
    ,val_cuda(NULL), grad_cuda(NULL), lapl_cuda(NULL), r_cuda(NULL),
    num_r(0)
#endif
  {
#ifdef HAVE_CUDA
    cudaMalloc((void**)&GVecs_cuda, gvecs.size()*sizeof(Vec3));
    cudaMalloc((void**)&Coefs_cuda, gvecs.size()*sizeof(complex<double>));
    cudaMemcpy(GVecs_cuda, &(GVecs(0)[0]), 3*sizeof(double)*gvecs.size(),
	       cudaMemcpyHostToDevice);
    cudaMemcpy(Coefs_cuda, &Coefs(0), 2*sizeof(double)*gvecs.size(),
	       cudaMemcpyHostToDevice);
#endif
  }

#ifdef HAVE_CUDA
  ~PWOrb_z() 
  {
    cudaFree (GVecs_cuda);
    cudaFree (Coefs_cuda);
    if (r_cuda)       cudaFree (r_cuda);
    if (val_cuda)     cudaFree (val_cuda);
    if (grad_cuda)    cudaFree (grad_cuda);
    if (lapl_cuda)    cudaFree (lapl_cuda);
  }
#endif

};

class SplineOrb_z : public OrbFunctor
{
private:
  UBspline_3d_z *Spline;
  LatticeClass &Lattice;
  Vec3 k;
public:
  complex<double> operator()(Vec3 r)
  {
    Vec3 u = Lattice.r2u(r);
    for (int i=0; i<3; i++)
      u[i] -= floor(u[i]);
    complex<double> val;
    eval_UBspline_3d_z (Spline, u[0], u[1], u[2], &val);
    double s,c, phase;
    phase = -dot (r, k);
    sincos(phase, &s,&c);
    val *= complex<double>(c,s);
    return val;
  }

  void operator()(Vec3 r, complex<double> &val,
		  TinyVector<complex<double>,3> &grad,
		  complex<double> &lapl)
  {
    Vec3 u = Lattice.r2u(r);
    for (int i=0; i<3; i++)
      u[i] -= floor(u[i]);
    TinyVector<complex<double>,3> g;
    TinyMatrix<complex<double>,3,3> hess;
    eval_UBspline_3d_z_vgh (Spline, u[0], u[1], u[2], &val,
			    &g[0], &hess(0,0));
    grad = complex<double>();
    lapl = complex<double>();
    for (int i=0; i<3; i++)
      for (int j=0; j<3; j++) {
	grad[i] += Lattice.GetB()(i,j) * g[j];
	lapl += hess (i,j) * Lattice.GetBBt()(i,j);
      }
    double s,c, phase;
    phase = -dot (r, k);
    TinyVector<complex<double>,3> ck;
    ck[0] = k[0];
    ck[1] = k[1];
    ck[2] = k[2];
    sincos(phase, &s,&c);
    complex<double> emikr(c,s);
    complex<double> eye(0.0,1.0);
    complex<double> k_dot_grad = k[0]*grad[0] + k[1]*grad[1] + k[2]*grad[2];
    lapl = emikr*(lapl - dot(k,k)*val - 2.0*eye*k_dot_grad);
    grad = emikr*(grad - eye*val*k);
    val *= emikr;
  }

  SplineOrb_z(LatticeClass &lattice, UBspline_3d_z *spline,
	      Vec3 kpoint)
    : Spline(spline), Lattice(lattice), k(kpoint)
  {
    // nothing else for now.
  }
};
  



#endif
