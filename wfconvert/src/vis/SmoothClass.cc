#include "SmoothClass.h"
#include <complex>


OnePath* SmoothClass::SmoothClosedPath (OnePath &oldPath)
{
  OnePath &newPath = *(new OnePath);
  newPath.Closed = true;
  newPath.Winding = oldPath.Winding;
  // The -1 ignores the path close to get proper periodicity
  int N = oldPath.Path.size()-1;
  double Ninv = 1.0/(double)N;
  int Numk = (int) floor (SmoothLevel+0.0001)+1;
  double NumkInv = 1.0/Numk;
  // First, compute fourier coefficients
  vector<TinyVector<complex<double>,3> > Fk(Numk);
  for (int k=0; k<Numk; k++) 
    for (int dim=0; dim<3; dim++) {
      Fk[k][dim] = complex<double>(0.0, 0.0);
      for (int n=0; n<N; n++) {
	double real, imag, phase;
	phase = 2.0*M_PI*n*k*Ninv;
	real = oldPath.Path[n][dim] * cos(phase);
	imag = oldPath.Path[n][dim] * sin(phase);
	Fk[k][dim] += Ninv*complex<double>(real, imag);
      }
    }
  
  // Figure out the appropriate scale factor by computing
  // power in real and fourier spectrum;
  Vec3 ravg(0.0, 0.0, 0.0);
  double realPower = 0.0;
  for (int i=0; i<N; i++)
    ravg += Ninv * oldPath.Path[i];
  for (int i=0; i<N; i++)
    realPower += Ninv*dot (oldPath.Path[i]-ravg, oldPath.Path[i]-ravg);
  
  double kPower = 0.0;
  for (int k=1; k<Numk; k++)
    for (int dim=0; dim<3; dim++)
      kPower += real(Fk[k][dim]*conj(Fk[k][dim]));
  
  double scale = sqrt(2.0*realPower/kPower);
  
  // Now, create new path, taking real part of inverse transform. 
  N = 15*Numk+1;
  Ninv = 1.0/(N-1);
  newPath.Path.resize(N);
  for (int n=0; n<N; n++) 
    for (int dim=0; dim<3; dim++) {
      complex<double> r(0.0, 0.0);
      for (int k=0; k<Numk; k++) {
	double phase, real, imag;
	phase = -2.0*M_PI*n*k*Ninv;
	real = cos(phase);
	imag = sin(phase);
	complex<double> z(real, imag);
	if (k != 0)
	  r += scale * Fk[k][dim]*z;
	else
	  r += Fk[k][dim]*z;
      }
      newPath.Path[n][dim] = r.real();
    }
  return (&newPath);
}


OnePath* SmoothClass::SmoothOpenPath2 (OnePath &oldPath)
{
  OnePath &newPath = *(new OnePath);
  newPath.Winding = oldPath.Winding;
  newPath.Closed = false;
  // The -1 ignores the path close to get proper periodicity
  int N = oldPath.Path.size()-1;
  double Ninv = 1.0/(double)N;
  int Numk = (int) floor (SmoothLevel+0.00001)+1;
  double NumkInv = 1.0/Numk;
  // First, compute fourier coefficients
  vector<TinyVector<complex<double>,3> > Fk(Numk);
  Vec3 diff = oldPath.Path[N] - oldPath.Path[0];

  for (int n=0; n<N; n++) 
    oldPath.Path[n] -= Ninv*n*diff;

  for (int k=0; k<Numk; k++) 
    for (int dim=0; dim<3; dim++) {
      Fk[k][dim] = complex<double>(0.0, 0.0);
      for (int n=0; n<N; n++) {
	double frac = n*Ninv;
	double real, imag, phase;
	phase = 2.0*M_PI*n*k*Ninv;
	real = oldPath.Path[n][dim] * cos(phase);
	imag = oldPath.Path[n][dim] * sin(phase);
	Fk[k][dim] += Ninv*complex<double>(real, imag);
      }
    }
  
  // Figure out the appropriate scale factor by computing
  // power in real and fourier spectrum;
  Vec3 ravg(0.0, 0.0, 0.0);
  double realPower = 0.0;
  for (int i=0; i<N; i++)
    ravg += Ninv * oldPath.Path[i];
  for (int i=0; i<N; i++)
    realPower += Ninv*dot (oldPath.Path[i]-ravg, oldPath.Path[i]-ravg);
  
  double kPower = 0.0;
  for (int k=1; k<Numk; k++)
    for (int dim=0; dim<3; dim++)
      kPower += real(Fk[k][dim]*conj(Fk[k][dim]));
  
  double scale = sqrt(2.0*realPower/kPower);
  
  // Now, create new path, taking real part of inverse transform. 
  N = 20*Numk+1;
  Ninv = 1.0/(N-1);
  newPath.Path.resize(N);
  for (int n=0; n<N; n++) {
    double frac = Ninv * n;
    for (int dim=0; dim<3; dim++) {
      complex<double> r(0.0, 0.0);
      for (int k=0; k<Numk; k++) {
	double phase, real, imag;
	phase = -2.0*M_PI*n*k*Ninv;
	real = cos(phase);
	imag = sin(phase);
	complex<double> z(real, imag);
	if (k != 0)
	  r += scale * Fk[k][dim]*z;
	else
	  r += Fk[k][dim]*z;
      }
      newPath.Path[n][dim] = r.real();
    }
  }

  for (int n=0; n<N; n++) 
    newPath.Path[n] += Ninv*n*diff;

  return (&newPath);
}



OnePath* SmoothClass::SmoothOpenPath (OnePath &oldPath)
{
  OnePath &newPath = *(new OnePath);
  newPath.Closed = false;
  newPath.Winding = oldPath.Winding;
  int N = oldPath.Path.size();
  double Ninv = 1.0/(double)(2*N-1);
  int Numk = (int) floor (SmoothLevel+0.00001) + 1;
  double NumkInv = 1.0/Numk;
  // First, compute fourier coefficients
  vector<TinyVector<complex<double>,3> > Fk(Numk);
  for (int k=0; k<Numk; k++) 
    for (int dim=0; dim<3; dim++) {
      Fk[k][dim] = complex<double>(0.0, 0.0);
      for (int n=-N+1; n<N; n++) {
	double real, imag, phase;
	phase = 2.0*M_PI*n*k*Ninv;
	real = oldPath.Path[abs(n)][dim] * cos(phase);
	imag = oldPath.Path[abs(n)][dim] * sin(phase);
	Fk[k][dim] += Ninv*complex<double>(real, imag);
      }
    }
  
  // Figure out the appropriate scale factor by computing
  // power in real and fourier spectrum;
  Vec3 ravg(0.0, 0.0, 0.0);
  double realPower = 0.0;
  for (int i=0; i<N; i++)
    ravg += Ninv * oldPath.Path[i];
  for (int i=0; i<N; i++)
    realPower += Ninv*dot (oldPath.Path[i]-ravg, oldPath.Path[i]-ravg);
  
  double kPower = 0.0;
  for (int k=1; k<Numk; k++)
    for (int dim=0; dim<3; dim++)
      kPower += 0.5*real(Fk[k][dim]*conj(Fk[k][dim]));
  
  double scale = sqrt(1.0*realPower/kPower);
  
  // Now, create new path, taking real part of inverse transform. 
  N = 20*Numk+1;
  Ninv = 1.0/(double)(2*N-1);
  newPath.Path.resize(N);
  for (int n=0; n<N; n++) 
    for (int dim=0; dim<3; dim++) {
      complex<double> r(0.0, 0.0);
      for (int k=0; k<Numk; k++) {
	double phase, real, imag;
	phase = -2.0*M_PI*n*k*Ninv;
	real = cos(phase);
	imag = sin(phase);
	complex<double> z(real, imag);
	if (k != 0)
	  r += scale * Fk[k][dim]*z;
	else
	  r += Fk[k][dim]*z;
      }
      newPath.Path[n][dim] = r.real();
    }

  return (&newPath);
}

void SmoothClass::SmoothPaths(vector<OnePath*> &inList)
{
  for (int pi=0; pi<inList.size(); pi++) {
    OnePath &oldPath = *(inList[pi]);
    OnePath *newPath;

    if (oldPath.Closed)
      newPath = SmoothClosedPath (oldPath);
    else
      newPath = SmoothOpenPath2 (oldPath);

    delete &oldPath;
    inList[pi] = newPath;
  }
}
