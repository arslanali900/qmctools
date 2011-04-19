#ifndef PLANE_WAVE_SUM_H
#define PLANE_WAVE_SUM_H

#include <complex>

void
plane_wave_sum (double *r, double g[], 
		double kx, double ky, double kz,
		std::complex<double> coefs[],
		std::complex<double> vals[], 
		int numr, int numg);

void
plane_wave_sum (double r[], double g[], 
		double kx, double ky, double kz,
		std::complex<double> coefs[],
		std::complex<double> vals[],
		std::complex<double> lapl[], 
		int numr, int numg);


void
plane_wave_sum (double r[], double g[], 
		double kx, double ky, double kz,
		std::complex<double> coefs[],
		std::complex<double> vals[],
		std::complex<double> grads[],
		std::complex<double> lapl[], 
		int numr, int numg);

#endif
