#include <cuda_runtime_api.h>
#include <vector_types.h>
#include <iostream>
//#include <gsl/gsl_sf_bessel.h>
#include <omp.h>
using namespace std;

template <typename T> struct vec3         { typedef float   Type; }; // dummy 
template <>           struct vec3<float>  { typedef float3  Type; };
template <>           struct vec3<double> { typedef double3 Type; };


extern int gsl_sf_bessel_jl_steed_array(const int lmax, const double val, double* jlVec);

template<typename tp> __device__
void bessel_steed_array(const int lmax, const tp x, tp* jl_x) {
// use steed/barnett method to calculate regular array of spherical bessel
// function from 0 to l.  This follows the gsl implementation, but is meant
// to be used on the gpus
   if (lmax < 0 || x < 0.0) {
      int j;
      for (j = 0; j<=lmax; j++) jl_x[j] = 0.0;
   } else if (x == 0.0) {
      int j;
      for (j = 1; j<=lmax; j++) jl_x[j] = 0.0;
      jl_x[0] = 1.0;
   } else if (x < 0.00024) {
      tp inv_fact = 1.0;
      tp x_l = 1.0;
      int l;
      for(l=0; l<=lmax; l++) {
	 jl_x[l]  = x_l * inv_fact;
	 jl_x[l] *= 1.0 - 0.5*x*x/(2.0*l+3.0);
	 inv_fact /= 2.0*l+3.0;
	 x_l *= x;
      }
   } else {
      // Steed/Barnett algorithm [Comp. Phys. Comm. 21, 297 (1981)]
      tp x_inv = 1.0/x;
      tp W = 2.0*x_inv;
      tp F = 1.0;
      tp FP = (lmax+1.0) * x_inv;
      tp B = 2.0*FP + x_inv;
      tp D = 1.0/B;
      tp del = -D;
      
      FP += del;
      
      /* continued fraction */
      do {
	 B += W;
	 D = 1.0/(B-D);
	 del *= (B*D - 1.);
	 FP += del;
	 if(D < 0.0) F = -F;
      } while(fabs(del) >= fabs(FP) * 1.19209e-07);
      
      FP *= F;
      
      if(lmax > 0) {
	 /* downward recursion */
	 tp XP2 = FP;
	 tp PL = lmax * x_inv;
	 int L  = lmax;
	 int LP;
	 jl_x[lmax] = F;
	 for(LP = 1; LP<=lmax; LP++) {
	    jl_x[L-1] = PL * jl_x[L] + XP2;
	    FP = PL*jl_x[L-1] - jl_x[L];
	    XP2 = FP;
	    PL -= x_inv;
	    --L;
	 }
	 F = jl_x[0];
      }  
      
      /* normalization */
      W = x_inv / hypot(FP, F);
      jl_x[0] = W*F;
      if(lmax > 0) {
	int L;
	for(L=1; L<=lmax; L++) {
          jl_x[L] *= W;
	}
      }
   }
}

// This is returning what is needed by wfconv, 
// NOTE!!!! This is not what would normally be returned,  
// this is the complex conjugate of the result multiplied by 4*Pi* (-i)^l
template<typename T,int lmax> __global__
void CalcYlmComplexCuda (T* Ghats, T *ylmptr, int numG) {
  const T fourPiInv = 0.0795774715459477;
  
  int rid = blockIdx.x*blockDim.x + threadIdx.x;
  
  if (rid < numG) {
    T costheta, sintheta, cottheta;
    T cosphi, sinphi;
    costheta = Ghats[3*rid+2];
    sintheta = sqrt(1.0-costheta*costheta);
    if (sintheta > 1.0e-6) {
      cottheta = costheta/sintheta;
      cosphi=Ghats[3*rid+0]/sintheta;
      sinphi=Ghats[3*rid+1]/sintheta; 
    }
    else {
      sintheta = 1.0e-8;
      cottheta = costheta * 1.0e8;
      cosphi = 1.0;
      sinphi = 0.0;
    }
    
    T e2iphi_r = cosphi;
    T e2iphi_i = sinphi;
    
    T lsign = 1.0;
    T dl = 0.0;
    
    T XlmVec[2*lmax+1];
    T dXlmVec[2*lmax+1];
    
    for (int l=0; l<=lmax; l++) {
      XlmVec[2*l]  = lsign;  
      dXlmVec[2*l] = dl * cottheta * XlmVec[2*l];
      XlmVec[0]    = lsign*XlmVec[2*l];
      dXlmVec[0]   = lsign*dXlmVec[2*l];
      T dm = dl;
      T msign = lsign;
      for (int m=l; m>0; m--) {
	T tmp = sqrt((dl+dm)*(dl-dm+1.0));
	XlmVec[l+m-1]  = -(dXlmVec[l+m] + dm*cottheta*XlmVec[l+m])/ tmp;
	dXlmVec[l+m-1] = (dm-1.0)*cottheta*XlmVec[l+m-1] + XlmVec[l+m]*tmp;
	// Copy to negative m
	XlmVec[l-(m-1)]  = -msign* XlmVec[l+m-1];
	dXlmVec[l-(m-1)] = -msign*dXlmVec[l+m-1];
	msign *= -1.0;
	dm -= 1.0;
      }
      T sum = 0.0;
      for (int m=-l; m<=l; m++) 
	sum += XlmVec[l+m]*XlmVec[l+m];
      // Now, renormalize the Ylms for this l
      T norm = sqrt((2.0*dl+1.0)*fourPiInv / sum);
      for (int m=-l; m<=l; m++) {
	XlmVec[l+m]  *= norm;
	dXlmVec[l+m] *= norm;
      }
      
      // Multiply by azimuthal phase and store in YlmVec
      T e2imphi_r = 1.0;
      T e2imphi_i = 0.0;
      for (int m=0; m<=l; m++) {
	ylmptr[2*(rid+numG*(l*(l+1)+m))] = XlmVec[l+m]*e2imphi_r;
	ylmptr[2*(rid+numG*(l*(l+1)+m))+1] = XlmVec[l+m]*e2imphi_i;
	ylmptr[2*(rid+numG*(l*(l+1)-m))] = XlmVec[l-m]*e2imphi_r;
	ylmptr[2*(rid+numG*(l*(l+1)-m))+1] = XlmVec[l-m]*-e2imphi_i;
	//ylmptr[maxindex*rid+2*(l*(l+1)+m)] = XlmVec[l+m]*e2imphi_r;
	//ylmptr[maxindex*rid+2*(l*(l+1)+m)+1] = XlmVec[l+m]*e2imphi_i;
	//ylmptr[maxindex*rid+2*(l*(l+1)-m)] = XlmVec[l-m]*e2imphi_r;
	//ylmptr[maxindex*rid+2*(l*(l+1)-m)+1] = XlmVec[l-m]*-e2imphi_i;
	
	T real = e2imphi_r*e2iphi_r - e2imphi_i*e2iphi_i;
	T imag = e2imphi_r*e2iphi_i + e2imphi_i*e2iphi_r;
	e2imphi_r = real;
	e2imphi_i = imag;
      } 
      
      dl += 1.0;
      lsign *= -1.0;
    }
    
    // This is the part that takes the complex conjugate and multiplies the result by 4pi*(-i)^l
    T minusi2l_r = 1.0;
    T minusi2l_i = 0.0;
    for (int l = 0, lm=0; l <= lmax; l++) {
      for (int m=-l; m<=l; m++, lm++) {
	const int ylmindex = 2*(rid+numG*lm);
	const T rbit = ylmptr[ylmindex] * minusi2l_r + ylmptr[ylmindex+1] * minusi2l_i;
	const T ibit = ylmptr[ylmindex] * minusi2l_i - ylmptr[ylmindex+1] * minusi2l_r;
	ylmptr[ylmindex] = rbit / fourPiInv;
	ylmptr[ylmindex+1] = ibit / fourPiInv;
	//T rbit = ylmptr[maxindex*rid+2*lm] * minusi2l_r + ylmptr[maxindex*rid+2*lm+1] * minusi2l_i;
	//T ibit = ylmptr[maxindex*rid+2*lm] * minusi2l_i - ylmptr[maxindex*rid+2*lm+1] * minusi2l_r;
	//ylmptr[maxindex*rid+2*lm] = rbit / fourPiInv;
	//ylmptr[maxindex*rid+2*lm+1] = ibit / fourPiInv;
      }
      T tmp1 = minusi2l_i;
      T tmp2 = -minusi2l_r;
      minusi2l_r = tmp1;
      minusi2l_i = tmp2;
    }
  }
}


// basically this is just here to contain the template stuff
// memory should be allocated elsewhere
template<typename T>
void CalcYlmComplexCudaLayer(T* dev_rhats, T* dev_ylms, int lmax, int N) {
  
  const int BS = 32;
  int Nblocks = (N+BS-1)/BS;
  dim3 dimGrid(Nblocks);
  dim3 dimBlock(BS);
  
  if (lmax == 0)
    return;
  else if (lmax == 1)
    CalcYlmComplexCuda<T,1><<<dimGrid,dimBlock>>>(dev_rhats, dev_ylms, N);
  else if (lmax == 2)
    CalcYlmComplexCuda<T,2><<<dimGrid,dimBlock>>>(dev_rhats, dev_ylms, N);
  else if (lmax == 3)
    CalcYlmComplexCuda<T,3><<<dimGrid,dimBlock>>>(dev_rhats, dev_ylms, N);
  else if (lmax == 4)
    CalcYlmComplexCuda<T,4><<<dimGrid,dimBlock>>>(dev_rhats, dev_ylms, N);
  else if (lmax == 5)
    CalcYlmComplexCuda<T,5><<<dimGrid,dimBlock>>>(dev_rhats, dev_ylms, N);
  else if (lmax == 6)
    CalcYlmComplexCuda<T,6><<<dimGrid,dimBlock>>>(dev_rhats, dev_ylms, N);
  else if (lmax == 7)
    CalcYlmComplexCuda<T,7><<<dimGrid,dimBlock>>>(dev_rhats, dev_ylms, N);
  else if (lmax == 8)
    CalcYlmComplexCuda<T,8><<<dimGrid,dimBlock>>>(dev_rhats, dev_ylms, N);
  else if (lmax == 9)
    CalcYlmComplexCuda<T,9><<<dimGrid,dimBlock>>>(dev_rhats, dev_ylms, N);
  else if (lmax == 10)
    CalcYlmComplexCuda<T,10><<<dimGrid,dimBlock>>>(dev_rhats, dev_ylms, N);
  
  cudaThreadSynchronize();
  cudaError_t err = cudaGetLastError();
  if (err != cudaSuccess) {
    fprintf (stderr, "CUDA error in CalcYlmComplexCuda:\n  %s\n",
	     cudaGetErrorString(err));
    abort();
  }
}


// On input, Ghats contains the g-vectors themselves, on output it has the ghats
template<typename T, int lmax> __global__ void
setupComputeKernel(const typename vec3<T>::Type k, T* pos, 
                   T* Ghats, T* phase_shifts, T* jlArray_spline, 
                   T* jlArray_poly, int numG, int numats, int SplinePoints, int PolyPoints, 
                   T dr, T poly_dr) {
  
  T jlVec[lmax+1];

  int gn = blockIdx.x* blockDim.x + threadIdx.x;
  if (gn < numG) {

  //  for (int gn=0; gn < numG; gn++) {
    
    // Compute Phase Shift
    Ghats[3*gn] += k.x;
    Ghats[3*gn+1] += k.y;
    Ghats[3*gn+2] += k.z;

    for (int iat = 0; iat < numats; iat++) {
      T val = -Ghats[3*gn]*pos[3*iat] - Ghats[3*gn+1]*pos[3*iat+1] - Ghats[3*gn+2]*pos[3*iat+2];
      phase_shifts[gn*numats+iat] = val;
    }

    // Compute appropriate Ghat
    T Gmag = sqrt(Ghats[3*gn]*Ghats[3*gn] + Ghats[3*gn+1]*Ghats[3*gn+1] + Ghats[3*gn+2]*Ghats[3*gn+2]);
    if (Gmag > 0.00001) {
      Ghats[gn*3] = (1.0/Gmag)*Ghats[gn*3];
      Ghats[gn*3+1] = (1.0/Gmag)*Ghats[gn*3+1];
      Ghats[gn*3+2] = (1.0/Gmag)*Ghats[gn*3+2];
    } else {
      Ghats[gn*3] = 0.0;
      Ghats[gn*3+1] = 0.0;
      Ghats[gn*3+2] = 1.0;
    }
    
    // Compute bessel functions
    for (int ir=0; ir<SplinePoints; ir++) {
      T r = dr*ir;
      bessel_steed_array(lmax, Gmag*r, jlVec);
      for (int l = 0; l <=lmax; l++) {
	jlArray_spline[ir+SplinePoints*(gn+numG*l)] = jlVec[l];
	//jlArray_spline[gn+numG*(ir+SplinePoints*l)] = jlVec[l];
	//jlArray_spline[gn+numG*(l+(lmax+1)*ir)] = jlVec[l];
	//jlArray_spline[l+(lmax+1)*(ir+SplinePoints*gn)] = jlVec[l];
      }
    }
    for (int ir = 0; ir<PolyPoints; ir++) {
      T r = poly_dr * ir;
      bessel_steed_array(lmax, Gmag*r, jlVec);
      for (int l = 0; l <=lmax; l++) {
	jlArray_poly[ir+PolyPoints*(gn+numG*l)] = jlVec[l];
	//jlArray_poly[gn+numG*(ir+PolyPoints*l)] = jlVec[l];
	//jlArray_poly[gn+numG*(l+(lmax+1)*ir)] = jlVec[l];
	//jlArray_poly[l+(lmax+1)*(ir+PolyPoints*gn)] = jlVec[l];
      }
    }
  }
}

template<typename T> void 
setupComputeKernelLayer(const typename vec3<T>::Type k, T* pos,
			T* Ghats, T* phase_shifts, T* jlArray_spline,
			T* jlArray_poly, int numG, int numats, int SplinePoints, int lmax, int PolyPoints,
			T dr, T poly_dr) {
  const int BS = 32;
  int Nblocks = (numG+BS-1)/BS;
  dim3 dimGrid(Nblocks);
  dim3 dimBlock(BS);
  
  if (lmax == 0)
    return;
  else if (lmax == 1)
    setupComputeKernel<T,1><<<dimGrid,dimBlock>>>(k,pos,Ghats,phase_shifts,jlArray_spline,jlArray_poly,numG,numats,SplinePoints,PolyPoints,dr,poly_dr);
  else if (lmax == 2)
    setupComputeKernel<T,2><<<dimGrid,dimBlock>>>(k,pos,Ghats,phase_shifts,jlArray_spline,jlArray_poly,numG,numats,SplinePoints,PolyPoints,dr,poly_dr);
  else if (lmax == 3)
    setupComputeKernel<T,3><<<dimGrid,dimBlock>>>(k,pos,Ghats,phase_shifts,jlArray_spline,jlArray_poly,numG,numats,SplinePoints,PolyPoints,dr,poly_dr);
  else if (lmax == 4)
    setupComputeKernel<T,4><<<dimGrid,dimBlock>>>(k,pos,Ghats,phase_shifts,jlArray_spline,jlArray_poly,numG,numats,SplinePoints,PolyPoints,dr,poly_dr);
  else if (lmax == 5)
    setupComputeKernel<T,5><<<dimGrid,dimBlock>>>(k,pos,Ghats,phase_shifts,jlArray_spline,jlArray_poly,numG,numats,SplinePoints,PolyPoints,dr,poly_dr);
  else if (lmax == 6)
    setupComputeKernel<T,6><<<dimGrid,dimBlock>>>(k,pos,Ghats,phase_shifts,jlArray_spline,jlArray_poly,numG,numats,SplinePoints,PolyPoints,dr,poly_dr);
  else if (lmax == 7)
    setupComputeKernel<T,7><<<dimGrid,dimBlock>>>(k,pos,Ghats,phase_shifts,jlArray_spline,jlArray_poly,numG,numats,SplinePoints,PolyPoints,dr,poly_dr);
  else if (lmax == 8)
    setupComputeKernel<T,8><<<dimGrid,dimBlock>>>(k,pos,Ghats,phase_shifts,jlArray_spline,jlArray_poly,numG,numats,SplinePoints,PolyPoints,dr,poly_dr);
  else if (lmax == 9)
    setupComputeKernel<T,9><<<dimGrid,dimBlock>>>(k,pos,Ghats,phase_shifts,jlArray_spline,jlArray_poly,numG,numats,SplinePoints,PolyPoints,dr,poly_dr);
  else if (lmax == 10)
    setupComputeKernel<T,10><<<dimGrid,dimBlock>>>(k,pos,Ghats,phase_shifts,jlArray_spline,jlArray_poly,numG,numats,SplinePoints,PolyPoints,dr,poly_dr);

  cudaThreadSynchronize();
  cudaError_t err = cudaGetLastError();
  if (err != cudaSuccess) {
    fprintf (stderr, "CUDA error in setupComputeKernel:\n  %s\n",
	     cudaGetErrorString(err));
    abort();
  }
} 


template<typename T, int BS> __global__ void
projectAnalyticPolyCudaKernelNew(const double* const orbitalCoefs, const double* const phase_shifts, 
				 const double* const jlArray_poly, const double* const YlmPtr, T* polyData,
				 int numats, int numSpins, int num_bands, int block_size, int numG, int PolyPoints, 
				 int lmax) {
  const int rblm = blockIdx.x * blockDim.x + threadIdx.x;
  const int lmmax = (lmax+1)*(lmax+1);
  const int lm = rblm / (PolyPoints*num_bands);
  const int rb = rblm % (PolyPoints*num_bands);
  const int ib = rb / PolyPoints;
  const int ir = rb % PolyPoints;
  const int l = static_cast<int>(floor(0.1+sqrt(static_cast<double>(lm))));

    
  // Number of atoms we can do sum for at once
  // Use half of the shared memory for this
  const int SimultaneousAts = 16;

  // It looks like if BS = 32 that this can be 16 and still fit in registers
  // Which is a must for performance!!!
  double reregisterbuf[SimultaneousAts];
  double imregisterbuf[SimultaneousAts];

  
  __shared__ float phase_shifts_re[BS][SimultaneousAts];
  __shared__ float phase_shifts_im[BS][SimultaneousAts];
  
  //__shared__ double phase_shifts_re[BS][SimultaneousAts];
  //__shared__ double phase_shifts_im[BS][SimultaneousAts];

  __shared__ double ylm_buffer_re[BS];
  __shared__ double ylm_buffer_im[BS];  
  __shared__ double orb_buffer_re[BS];
  __shared__ double orb_buffer_im[BS];
  
  if (ir < PolyPoints) {
    for (int is = 0; is < numSpins; is++) {
      int NAtBlocks = (numats+SimultaneousAts-1) / SimultaneousAts;
      for (int atblock = 0; atblock < NAtBlocks; atblock++) {
	const int offset = atblock*SimultaneousAts;
	const int end = min(numats, (atblock+1)*SimultaneousAts)-offset;
	for (int i = 0; i < end; i++) {
	  reregisterbuf[i] = 0;
	  imregisterbuf[i] = 0;
	}
	
	const int NSharedBlocks = (numG+BS-1)/BS;
	for (int sblock = 0; sblock < NSharedBlocks; sblock++) {
	  const int sboffset = sblock*BS;
	  const int sbend = min(numG, (sblock+1)*BS)-sboffset;
	  
	  if (threadIdx.x < sbend) {
 	    ylm_buffer_re[threadIdx.x] = YlmPtr[2*((sboffset+threadIdx.x)+numG*lm)];
	    ylm_buffer_im[threadIdx.x] = YlmPtr[2*((sboffset+threadIdx.x)+numG*lm)+1];
	    orb_buffer_re[threadIdx.x] = orbitalCoefs[2*((sboffset+threadIdx.x)+block_size*(ib+num_bands*is))];
	    orb_buffer_im[threadIdx.x] = orbitalCoefs[2*((sboffset+threadIdx.x)+block_size*(ib+num_bands*is))+1];
	    //Each thread now has one G vector, this is a loop over atoms
	    for (int i = 0; i < end; i++) {
	      sincos(phase_shifts[(sboffset+threadIdx.x)*numats+offset+i], &(phase_shifts_im[threadIdx.x][i]), &(phase_shifts_re[threadIdx.x][i]));
	    }
	  }
	  __syncthreads();
	  for (int iGIndex = 0; iGIndex < sbend; iGIndex++) {
	    const int iG = sboffset+iGIndex;	    
	    const int jlapindex = ir+PolyPoints*(iG+numG*l);
	    
	    const double jloylmre = jlArray_poly[jlapindex] *
	      (orb_buffer_re[iGIndex]*ylm_buffer_re[iGIndex] - orb_buffer_im[iGIndex]*ylm_buffer_im[iGIndex]);
	    const double jloylmim = jlArray_poly[jlapindex] *
	      (orb_buffer_re[iGIndex]*ylm_buffer_im[iGIndex] + orb_buffer_im[iGIndex]*ylm_buffer_re[iGIndex]);

	    for (int ti = 0; ti < end; ti++) {
	      
	      reregisterbuf[ti] += jloylmre*phase_shifts_re[iGIndex][ti] - jloylmim*phase_shifts_im[iGIndex][ti];
	      imregisterbuf[ti] += jloylmre*phase_shifts_im[iGIndex][ti] + jloylmim*phase_shifts_re[iGIndex][ti];
	    }
	    __syncthreads();
	  }
	  __syncthreads();
	}
	__syncthreads();
	
	for (int ti = 0; ti < end; ti++) {
	  int iat = offset+ti;
	  const int pdIndex = 2*(lm+lmmax*(ir+PolyPoints*(ib+num_bands*(is+numSpins*iat))));
	  polyData[pdIndex] += reregisterbuf[ti];
	  polyData[pdIndex+1] += imregisterbuf[ti];	  
	}
	__syncthreads();
	
      }
    }
  }
}

/*
void setupCPUComputeKernelDebug(const typename vec3<double>::Type k, double* gpupositions,
				double* gpughats, double* gpuPhaseShifts, double* gpujlArray_spline,
				double* gpujlArray_poly, double* gpuylmptr, int numG, int numats, 
				int SplinePoints, int lmax, int PolyPoints, double dr, double poly_dr) {

  const int numlm = (lmax+1)*(lmax+1);
  double* pos = new double[3*numats];
  double* Ghats = new double[3*numG];
  double* phase_shifts = new double[numG*numats];
  double* jlArray_spline = new double[numG*SplinePoints*(lmax+1)];
  double* jlArray_poly = new double[numG*PolyPoints*(lmax+1)];
  double* ylmptr = new double[2*numG*numlm];
  
  // now that memory is allocated, initialize that memory with the values stored on the GPU
  cudaMemcpy(pos, gpupositions, 3*numats*sizeof(double), cudaMemcpyDeviceToHost);
  cudaMemcpy(Ghats, gpughats, 3*numG*sizeof(double), cudaMemcpyDeviceToHost);
  cudaMemcpy(phase_shifts, gpuPhaseShifts, numG*numats*sizeof(double), cudaMemcpyDeviceToHost);
  cudaMemcpy(jlArray_spline, gpujlArray_spline, numG*SplinePoints*(lmax+1)*sizeof(double), cudaMemcpyDeviceToHost);
  cudaMemcpy(jlArray_poly, gpujlArray_poly, numG*PolyPoints*(lmax+1)*sizeof(double), cudaMemcpyDeviceToHost);
  cudaMemcpy(ylmptr, gpuylmptr, 2*numG*numlm*sizeof(double), cudaMemcpyDeviceToHost);

  // now do the appropriate calculations
  double jlVec[lmax+1];
  for (int gn = 0; gn < numG; gn++) {
    //    cout << "gn = " << gn << endl;
    // Compute Phase Shift
    Ghats[3*gn] += k.x;
    Ghats[3*gn+1] += k.y;
    Ghats[3*gn+2] += k.z;
    
    for (int iat = 0; iat < numats; iat++) {
      double val = -Ghats[3*gn]*pos[3*iat] - Ghats[3*gn+1]*pos[3*iat+1] - Ghats[3*gn+2]*pos[3*iat+2];
      //double val = -Ghats[3*gn]*pos.x -Ghats[3*gn+1]*pos.y -Ghats[3*gn+2]*pos.z;
      phase_shifts[gn*numats+iat] = val;
      //phase_shifts[gn*2] = cos(val);
      //phase_shifts[gn*2+1] = sin(val);
    }
    
    // Compute appropriate Ghat
    double Gmag = std::sqrt(Ghats[3*gn]*Ghats[3*gn] + Ghats[3*gn+1]*Ghats[3*gn+1] + Ghats[3*gn+2]*Ghats[3*gn+2]);
    if (Gmag > 0.0) {
      Ghats[gn*3] = (1.0/Gmag)*Ghats[gn*3];
      Ghats[gn*3+1] = (1.0/Gmag)*Ghats[gn*3+1];
      Ghats[gn*3+2] = (1.0/Gmag)*Ghats[gn*3+2];
    } else {
      Ghats[gn*3] = 0.0;
      Ghats[gn*3+1] = 0.0;
      Ghats[gn*3+2] = 1.0;
    }
    
    // Compute bessel functions
    for (int ir=0; ir<SplinePoints; ir++) {
      double r = dr*ir;
      gsl_sf_bessel_jl_steed_array(lmax, Gmag*r, jlVec);
      for (int l = 0; l <=lmax; l++) {
	//jlArray_spline[l+(lmax+1)*(ir+SplinePoints*gn)] = jlVec[l];
	//jlArray_spline[gn+numG*(l+(lmax+1)*ir)] = jlVec[l];
	jlArray_spline[ir+SplinePoints*(gn+numG*l)] = jlVec[l];
      }
    }
    for (int ir = 0; ir<PolyPoints; ir++) {
      double r = poly_dr * ir;
      gsl_sf_bessel_jl_steed_array(lmax, Gmag*r, jlVec);
      for (int l = 0; l <=lmax; l++) {
	//jlArray_poly[l+(lmax+1)*(ir+PolyPoints*gn)] = jlVec[l];
	//jlArray_poly[gn+numG*(l+(lmax+1)*ir)] = jlVec[l];
	jlArray_poly[ir+PolyPoints*(gn+numG*l)] = jlVec[l];
      }
    }
  }
  
  cout << "Finished first loop over iG" << endl;
  
  for (int gn = 0; gn < numG; gn++) {
    const double fourPiInv = 0.0795774715459477;
    //int maxindex = (lmax+1)*(lmax+1)*2;
    int rid = gn;
      
    double costheta, sintheta, cottheta;
    double cosphi, sinphi;
    costheta = Ghats[3*rid+2];
    sintheta = std::sqrt(1.0-costheta*costheta);
    if (sintheta > 1.0e-6) {
      cottheta = costheta/sintheta;
      cosphi=Ghats[3*rid+0]/sintheta;
      sinphi=Ghats[3*rid+1]/sintheta; 
    }
    else {
      sintheta = 1.0e-8;
      cottheta = costheta * 1.0e8;
      cosphi = 1.0;
      sinphi = 0.0;
    }

    double e2iphi_r = cosphi;
    double e2iphi_i = sinphi;
  
    double lsign = 1.0;
    double dl = 0.0;
  
    double XlmVec[2*lmax+1];
    double dXlmVec[2*lmax+1];

    for (int l=0; l<=lmax; l++) {
      XlmVec[2*l]  = lsign;  
      dXlmVec[2*l] = dl * cottheta * XlmVec[2*l];
      XlmVec[0]    = lsign*XlmVec[2*l];
      dXlmVec[0]   = lsign*dXlmVec[2*l];
      double dm = dl;
      double msign = lsign;
      for (int m=l; m>0; m--) {
	double tmp = std::sqrt((dl+dm)*(dl-dm+1.0));
	XlmVec[l+m-1]  = -(dXlmVec[l+m] + dm*cottheta*XlmVec[l+m])/ tmp;
	dXlmVec[l+m-1] = (dm-1.0)*cottheta*XlmVec[l+m-1] + XlmVec[l+m]*tmp;
	// Copy to negative m
	XlmVec[l-(m-1)]  = -msign* XlmVec[l+m-1];
	dXlmVec[l-(m-1)] = -msign*dXlmVec[l+m-1];
	msign *= -1.0;
	dm -= 1.0;
      }
      double sum = 0.0;
      for (int m=-l; m<=l; m++) 
	sum += XlmVec[l+m]*XlmVec[l+m];
      // Now, renormalize the Ylms for this l
      double norm = std::sqrt((2.0*dl+1.0)*fourPiInv / sum);
      for (int m=-l; m<=l; m++) {
	XlmVec[l+m]  *= norm;
	dXlmVec[l+m] *= norm;
      }
      
      // Multiply by azimuthal phase and store in YlmVec
      double e2imphi_r = 1.0;
      double e2imphi_i = 0.0;
      for (int m=0; m<=l; m++) {
	ylmptr[2*(rid+numG*(l*(l+1)+m))] = XlmVec[l+m]*e2imphi_r;
	ylmptr[2*(rid+numG*(l*(l+1)+m))+1] = XlmVec[l+m]*e2imphi_i;
	ylmptr[2*(rid+numG*(l*(l+1)-m))] = XlmVec[l-m]*e2imphi_r;
	ylmptr[2*(rid+numG*(l*(l+1)-m))+1] = XlmVec[l-m]*-e2imphi_i;
	//ylmptr[maxindex*rid+2*(l*(l+1)+m)] = XlmVec[l+m]*e2imphi_r;
	//ylmptr[maxindex*rid+2*(l*(l+1)+m)+1] = XlmVec[l+m]*e2imphi_i;
	//ylmptr[maxindex*rid+2*(l*(l+1)-m)] = XlmVec[l-m]*e2imphi_r;
	//ylmptr[maxindex*rid+2*(l*(l+1)-m)+1] = XlmVec[l-m]*-e2imphi_i;
	
	double real = e2imphi_r*e2iphi_r - e2imphi_i*e2iphi_i;
	double imag = e2imphi_r*e2iphi_i + e2imphi_i*e2iphi_r;
	e2imphi_r = real;
	e2imphi_i = imag;
      } 
      
      dl += 1.0;
      lsign *= -1.0;
    }
    
    // This is the part that takes the complex conjugate and multiplies the result by 4pi*(-i)^l    
    double minusi2l_r = 1.0;
    double minusi2l_i = 0.0;
    for (int l = 0, lm=0; l <= lmax; l++) {
      for (int m=-l; m<=l; m++, lm++) {	
	const int ylmindex = 2*(rid+numG*lm);
	const double rbit = ylmptr[ylmindex] * minusi2l_r + ylmptr[ylmindex+1] * minusi2l_i;
	const double ibit = ylmptr[ylmindex] * minusi2l_i - ylmptr[ylmindex+1] * minusi2l_r;
	ylmptr[ylmindex] = rbit / fourPiInv;
	ylmptr[ylmindex+1] = ibit / fourPiInv;
	//T rbit = ylmptr[maxindex*rid+2*lm] * minusi2l_r + ylmptr[maxindex*rid+2*lm+1] * minusi2l_i;
	//T ibit = ylmptr[maxindex*rid+2*lm] * minusi2l_i - ylmptr[maxindex*rid+2*lm+1] * minusi2l_r;
	//ylmptr[maxindex*rid+2*lm] = rbit / fourPiInv;
	//ylmptr[maxindex*rid+2*lm+1] = ibit / fourPiInv;
      }
      double tmp1 = minusi2l_i;
      double tmp2 = -minusi2l_r;
      minusi2l_r = tmp1;
      minusi2l_i = tmp2;
    }
  }
  


  // now that everything has been calculated, can copy it back to the gpus and deallocate all of the memory
  cudaMemcpy(gpupositions, pos, 3*numats*sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(gpughats, Ghats, 3*numG*sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(gpuPhaseShifts, phase_shifts, numG*numats*sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(gpujlArray_spline, jlArray_spline, numG*SplinePoints*(lmax+1)*sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(gpujlArray_poly, jlArray_poly, numG*PolyPoints*(lmax+1)*sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(gpuylmptr, ylmptr, 2*numG*numlm*sizeof(double), cudaMemcpyHostToDevice);

  delete[] pos;
  delete[] Ghats;
  delete[] phase_shifts;
  delete[] jlArray_spline;
  delete[] jlArray_poly;
  delete[] ylmptr;
}    
*/

void projectAnalyticCudaDriver (const vec3<double>::Type k, double* pos, const double* gvecs, 
				const double* orbitalCoefs, double* splineData, double* polyData,
				int num_ks, int knum, int numG, int lmax, int numats, int SplinePoints, int PolyPoints, 
				int numSpins, int num_bands, double dr, double poly_dr) {


  int deviceCount;
  int devices[omp_get_num_threads()];
  cudaGetDeviceCount(&deviceCount);
  int num_appropriate=0, device=0;
  for (device = 0; device < deviceCount; ++device) {
    cudaDeviceProp deviceProp;
    cudaGetDeviceProperties(&deviceProp, device);
    if (((deviceProp.major >= 1) && (deviceProp.minor >= 3)) ||
        deviceProp.major >= 2) {
      devices[num_appropriate] = device;
      num_appropriate++;
    }
  }

  #pragma omp parallel for
  for (int i = 0; i < omp_get_num_threads(); i++) {
    cudaSetDevice(devices[i]);

    double* gpughats;
    double* gpuPhaseShifts;
    double* gpujlArray_spline;
    double* gpujlArray_poly;
    double* gpuylmptr;
    double* gpuOrbitalCoefs;
    double* gpuSplineData;
    double* gpuPolyData;  
    double* gpupositions;

    const int numlm = (lmax+1)*(lmax+1);
    static const int G_block_size = 16384; // In the future it might make sense to change this based on the
                                           // number of bands so that the GPU memory might be better used
    const int num_G_blocks = (numG+G_block_size-1) / G_block_size;

    cudaMalloc((void**) &gpupositions, 3*numats*sizeof(double));
    cudaMalloc((void**) &gpughats, 3*G_block_size*sizeof(double));
    //cudaMalloc((void**) &gpuPhaseShifts, 2*G_block_size*sizeof(double));
    cudaMalloc((void**) &gpuPhaseShifts, G_block_size*numats*sizeof(double));
    cudaMalloc((void**) &gpujlArray_spline, G_block_size*SplinePoints*(lmax+1)*sizeof(double));
    cudaMalloc((void**) &gpujlArray_poly, G_block_size*PolyPoints*(lmax+1)*sizeof(double));
    cudaMalloc((void**) &gpuylmptr, 2*G_block_size*numlm*sizeof(double));
    cudaMalloc((void**) &gpuOrbitalCoefs, numSpins*G_block_size*num_bands*2*sizeof(double));
    cudaMalloc((void**) &gpuSplineData, 2*num_bands*numSpins*SplinePoints*numlm*numats*sizeof(double));
    cudaMalloc((void**) &gpuPolyData, 2*num_bands*numSpins*PolyPoints*numlm*numats*sizeof(double));
    
    // These arrays are passed in as 0, so this serves to initialize the arrays on the GPUs
    cudaMemcpy(gpuSplineData, splineData, 2*num_bands*numSpins*SplinePoints*numlm*numats*sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(gpuPolyData, polyData, 2*num_bands*numSpins*PolyPoints*numlm*numats*sizeof(double), cudaMemcpyHostToDevice);  
    
    // Copy the positions in to the gpu
    cudaMemcpy(gpupositions, pos, 3*numats*sizeof(double), cudaMemcpyHostToDevice);
    
    // Space to use for copying orbital data in to the GPUs
    double* floatorbs = new double[numSpins*G_block_size*num_bands*2];
    
    for (int Gblock=i; Gblock < num_G_blocks; Gblock+= omp_get_num_threads()) {
      cout << "Starting Gblock " << Gblock+1 << "/" << num_G_blocks << " on thread " << i << endl;
      const int offset = Gblock*G_block_size;
      const int end = min(numG, (Gblock+1)*G_block_size) - offset;
      
      cudaMemcpy(gpughats, &(gvecs[3*offset]), 3*end*sizeof(double), cudaMemcpyHostToDevice);


      
      //cout << "About to start setupComputeKernel on thread " << i << endl;
      setupComputeKernelLayer<double>(k, gpupositions, gpughats, gpuPhaseShifts, gpujlArray_spline, gpujlArray_poly, 
				      end, numats, SplinePoints, lmax, PolyPoints, dr, poly_dr);
      //cout << "Finised setupComputeKernel on thread " << i << endl;
      //cout << "About to call calcYlmCuda on thread" << i << endl;
      CalcYlmComplexCudaLayer<double>(gpughats, gpuylmptr, lmax, end);
      //cout << "Finished calling calcYlmCuda on thread " << i << endl;
      

      /*
      cout << "About to start setupCPUComputeKernelDebug on thread " << i << endl;
      setupCPUComputeKernelDebug(k, gpupositions, gpughats, gpuPhaseShifts, gpujlArray_spline, 
				 gpujlArray_poly, gpuylmptr, end, numats, SplinePoints, lmax, 
				 PolyPoints, dr, poly_dr);
      cout << "Finished setupCPUComputeKernelDebug on thread " << i << endl;
      cout << "About to call calcYlmCuda on thread" << i << endl;
      CalcYlmComplexCudaLayer<double>(gpughats, gpuylmptr, lmax, end);
      cout << "Finished calling calcYlmCuda on thread " << i << endl;      
      */

      // Need to create an array that has the proper subset of orbitalCoefficients for this block
      for (int ib = 0; ib < num_bands; ib++) {
	for (int sigma = 0; sigma < numSpins; sigma++) {
	  for (int j = 0; j < 2*end; j++) {
	    floatorbs[j+2*G_block_size*(ib+num_bands*sigma)] = orbitalCoefs[2*offset+j+2*numG*(ib+num_bands*sigma)];
	  }
	}
      }
      
      cudaMemcpy(gpuOrbitalCoefs, floatorbs, numSpins*G_block_size*num_bands*2*sizeof(double), cudaMemcpyHostToDevice);
      
      const int BS=32;
      int numrblm = (lmax+1)*(lmax+1)*num_bands*SplinePoints;
      int Nblocks = (numrblm+BS-1)/BS;
      dim3 dimGrid(Nblocks);
      dim3 dimBlock(BS);
      
      cout << "About to call projectAnalyticSplineCudaKernel on thread " << i << endl;
      projectAnalyticPolyCudaKernelNew<double,BS><<<dimGrid,dimBlock>>>(gpuOrbitalCoefs, gpuPhaseShifts, gpujlArray_spline,
									gpuylmptr, gpuSplineData, numats, numSpins, num_bands, 
									G_block_size, end, SplinePoints, lmax);
      cudaThreadSynchronize();
      cudaError_t err = cudaGetLastError();
      if (err != cudaSuccess) {
	fprintf (stderr, "CUDA error in projectAnalyticSplineCudaKernel:\n  %s\n",
		 cudaGetErrorString(err));
	abort();
      }
      //cout << "Finished call to projectAnalyticSplineCudaKernel on thread " << i << endl;
      
      const int BS2=32;
      numrblm = (lmax+1)*(lmax+1)*num_bands*PolyPoints;
      Nblocks = (numrblm+BS2-1)/BS2;
      dim3 dimGrid2(Nblocks);
      dim3 dimBlock2(BS2);
      
      cout << "About to call projectAnalyticPolyCudaKernel on thread " << i << endl;
      projectAnalyticPolyCudaKernelNew<double,BS2><<<dimGrid2,dimBlock2>>>(gpuOrbitalCoefs, gpuPhaseShifts, gpujlArray_poly,
									   gpuylmptr, gpuPolyData, numats, numSpins, num_bands, 
									   G_block_size, end, PolyPoints, lmax);
      cudaThreadSynchronize();
      err = cudaGetLastError();
      if (err != cudaSuccess) {
	fprintf (stderr, "CUDA error in projectAnalyticPolyCudaKernel:\n  %s\n",
		 cudaGetErrorString(err));
	abort();
      }
      //cout << "Finished call to projectAnalyticPolyCudaKernel on thread " << i << endl;
    }
    
    #pragma omp critical
    {
      cout << "About to start copying data out on thread " << i << endl;
      double* polydatabuffer = new double[2*num_bands*numSpins*PolyPoints*numlm];
      double* splinedatabuffer = new double[2*num_bands*numSpins*SplinePoints*numlm];
      
      for (int ii = 0; ii < numats; ii++) {
	const int singleAtomSBufSize = 2*num_bands*numSpins*SplinePoints*numlm;
	cudaMemcpy(splinedatabuffer, gpuSplineData+ii*singleAtomSBufSize, singleAtomSBufSize*sizeof(double), cudaMemcpyDeviceToHost);
	for (int j = 0; j < singleAtomSBufSize; j++) {
	  splineData[j+ii*singleAtomSBufSize] += splinedatabuffer[j];
	}
	const int singleAtomPBufSize = 2*num_bands*numSpins*PolyPoints*numlm;
	cudaMemcpy(polydatabuffer, gpuPolyData+ii*singleAtomPBufSize, singleAtomPBufSize*sizeof(double), cudaMemcpyDeviceToHost);
	for (int j = 0; j < singleAtomPBufSize; j++) {
	  polyData[j+ii*singleAtomPBufSize] += polydatabuffer[j];
	}
      }
      
      delete[] polydatabuffer;
      delete[] splinedatabuffer;
      cout << "Finished copying data out on thread " << i << endl;
    }
    
    cudaFree(gpughats);
    cudaFree(gpuPhaseShifts);
    cudaFree(gpujlArray_spline);
    cudaFree(gpujlArray_poly);
    cudaFree(gpuylmptr);
    cudaFree(gpuOrbitalCoefs);
    cudaFree(gpuSplineData);
    cudaFree(gpuPolyData);
    if (floatorbs) {
      delete[] floatorbs;
    }
    cudaThreadExit();
  }
}
