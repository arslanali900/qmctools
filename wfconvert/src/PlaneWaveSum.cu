#include "PlaneWaveSum.h"

template<int BS, typename T>
__global__ void
plane_wave_sum_kernel_z (T *r, T *g,
			 T kx, T ky, T kz,
			 T *coefs, T *vals, 
			 int numr, int numg)
{
  __shared__ T r_shared[BS][3], g_shared[BS][3];
  __shared__ T coefs_shared[BS][2];
  T my_r[3];

  for (int i=0; i<3; i++) {
    int off = (3*blockIdx.x+i)*BS + threadIdx.x;
    if (off < 3*numr)
      r_shared[0][i*BS+threadIdx.x] = r[off];
  }
  __syncthreads();

  my_r[0] = r_shared[threadIdx.x][0];
  my_r[1] = r_shared[threadIdx.x][1];
  my_r[2] = r_shared[threadIdx.x][2];
  
  int numblocks = (numg+BS-1)/BS;
  T myval_real=T(), myval_imag=T();
  for (int iblock=0; iblock<numblocks; iblock++) {
    for (int i=0; i<3; i++) {
      int off = (3*iblock+i)*BS + threadIdx.x;
      if (off < 3*numg)
  	g_shared[0][i*BS+threadIdx.x] = g[off];
    }
    g_shared[threadIdx.x][0] += kx;
    g_shared[threadIdx.x][1] += ky;
    g_shared[threadIdx.x][2] += kz;
    int off = 2*iblock*BS + threadIdx.x;
    if (off < 2*numg)      coefs_shared[0][threadIdx.x] = coefs[off];
    off += BS;
    if (off < 2*numg)      coefs_shared[0][threadIdx.x+BS] = coefs[off];

    __syncthreads();
    int end = min (numg-iblock*BS, BS);
    for (int j=0; j<end; j++) {
      T phase = -(my_r[0]*g_shared[j][0] +
  		  my_r[1]*g_shared[j][1] +
  		  my_r[2]*g_shared[j][2]);
      T s,c;
      sincos(phase, &s, &c);
      myval_real += (c*coefs_shared[j][0] -
  		     s*coefs_shared[j][1]);
      myval_imag += (c*coefs_shared[j][1] +
  		     s*coefs_shared[j][0]);
    }
    __syncthreads();
  }
  /// Use coallesced writes
  __shared__ T vals_shared[2*BS];
  vals_shared[2*threadIdx.x+0] = myval_real;
  vals_shared[2*threadIdx.x+1] = myval_imag;
  __syncthreads();
  int off = 2*blockIdx.x*BS + threadIdx.x;
  if (off < 2*numr) 
    vals[off] = vals_shared[threadIdx.x];
  if (off+BS < 2*numr)
    vals[off+BS] = vals_shared[threadIdx.x+BS];
}

void
plane_wave_sum (double r[], double g[], 
		double kx, double ky, double kz,
		std::complex<double> coefs[],
		std::complex<double> vals[], 
		int numr, int numg)
{
  const int BS=32;
  dim3 dimBlock(BS);
  dim3 dimGrid((numr+BS-1)/BS);
  
  plane_wave_sum_kernel_z<BS,double><<<dimGrid,dimBlock>>> 
    (r, g, kx, ky, kz, (double*)coefs, (double*)vals, numr, numg);

  cudaThreadSynchronize();
  cudaError_t err = cudaGetLastError();
  if (err != cudaSuccess) {
    fprintf (stderr, "CUDA error in plane_wave_sum:\n  %s\n",
	     cudaGetErrorString(err));
    abort();
  }
}


template<int BS, typename T>
__global__ void
plane_wave_sum_kernel_z (T *r, T *g,
			 T kx, T ky, T kz,
			 T *coefs, T *vals, T *lapl,
			 int numr, int numg)
{
  __shared__ T r_shared[BS][3], g_shared[BS][3], g2_shared[BS];
  __shared__ T coefs_shared[BS][2];
  T my_r[3];

  for (int i=0; i<3; i++) {
    int off = (3*blockIdx.x+i)*BS + threadIdx.x;
    if (off < 3*numr)
      r_shared[0][i*BS+threadIdx.x] = r[off];
  }
  __syncthreads();

  my_r[0] = r_shared[threadIdx.x][0];
  my_r[1] = r_shared[threadIdx.x][1];
  my_r[2] = r_shared[threadIdx.x][2];
  
  int numblocks = (numg+BS-1)/BS;
  T myval_real=T(), myval_imag=T(), mylapl_real=T(), mylapl_imag=T();
  for (int iblock=0; iblock<numblocks; iblock++) {
    for (int i=0; i<3; i++) {
      int off = (3*iblock+i)*BS + threadIdx.x;
      if (off < 3*numg)
  	g_shared[0][i*BS+threadIdx.x] = g[off];
    }
    __syncthreads();
    g_shared[threadIdx.x][0] += kx;
    g_shared[threadIdx.x][1] += ky;
    g_shared[threadIdx.x][2] += kz;

    g2_shared[threadIdx.x] = (g_shared[threadIdx.x][0] * g_shared[threadIdx.x][0] +
			      g_shared[threadIdx.x][1] * g_shared[threadIdx.x][1] +
			      g_shared[threadIdx.x][2] * g_shared[threadIdx.x][2] );

    int off = 2*iblock*BS + threadIdx.x;
    if (off < 2*numg)      coefs_shared[0][threadIdx.x] = coefs[off];
    off += BS;
    if (off < 2*numg)      coefs_shared[0][threadIdx.x+BS] = coefs[off];

    __syncthreads();
    int end = min (numg-iblock*BS, BS);
    for (int j=0; j<end; j++) {
      T phase = -(my_r[0]*g_shared[j][0] +
  		  my_r[1]*g_shared[j][1] +
  		  my_r[2]*g_shared[j][2]);
      T s,c;
      sincos(phase, &s, &c);
      T re = (c*coefs_shared[j][0] -
	      s*coefs_shared[j][1]);
      T im = (c*coefs_shared[j][1] +
	      s*coefs_shared[j][0]);
      myval_real += re;
      myval_imag += im;
      mylapl_real -= g2_shared[j] * re;
      mylapl_imag -= g2_shared[j] * im;
    }
    __syncthreads();
  }
  /// Use coallesced writes for vals
  __shared__ T vals_shared[2*BS];
  vals_shared[2*threadIdx.x+0] = myval_real;
  vals_shared[2*threadIdx.x+1] = myval_imag;
  __syncthreads();
  int off = 2*blockIdx.x*BS + threadIdx.x;
  if (off < 2*numr) 
    vals[off] = vals_shared[threadIdx.x];
  if (off+BS < 2*numr)
    vals[off+BS] = vals_shared[threadIdx.x+BS];

  /// Use coallesced writes for lapl
  vals_shared[2*threadIdx.x+0] = mylapl_real;
  vals_shared[2*threadIdx.x+1] = mylapl_imag;
  __syncthreads();
  off = 2*blockIdx.x*BS + threadIdx.x;
  if (off < 2*numr) 
    lapl[off] = vals_shared[threadIdx.x];
  if (off+BS < 2*numr)
    lapl[off+BS] = vals_shared[threadIdx.x+BS];
}


void
plane_wave_sum (double r[], double g[], 
		double kx, double ky, double kz,
		std::complex<double> coefs[],
		std::complex<double> vals[],
		std::complex<double> lapl[], 
		int numr, int numg)
{
  const int BS=32;
  dim3 dimBlock(BS);
  dim3 dimGrid((numr+BS-1)/BS);
  
  plane_wave_sum_kernel_z<BS,double><<<dimGrid,dimBlock>>> 
    (r, g, kx, ky, kz, (double*)coefs, (double*)vals, 
     (double*)lapl, numr, numg);

  cudaThreadSynchronize();
  cudaError_t err = cudaGetLastError();
  if (err != cudaSuccess) {
    fprintf (stderr, "CUDA error in plane_wave_sum:\n  %s\n",
	     cudaGetErrorString(err));
    abort();
  }
}
