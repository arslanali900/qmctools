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

#include <gsl/gsl_math.h>
#include <gsl/gsl_specfunc.h>
#include <gsl/gsl_errno.h>
#include "SpecialFunctions.h"




/*double ModBesselScaled(int l, double x)
{
  gsl_sf_result result;
  int ErrorCode;
  double L = (double)l;
  double xinv = 1.0/x;
  double xinvsq = xinv*xinv;

  if (x > 1.0e2)
    {
      //fprintf (stderr, "New way.\n");
      double result = 0.5*xinv * exp(-0.5*xinv*(L*(L+1.0))) *
	(1.0 - 0.25*xinvsq*L*(L+1.0) +
	 L*(L+1.0)*(L-2.0)*(L+3.0)*xinv*xinvsq/24.0 +
	 L*(L+1.0)*(5.0*L*L + 5.0*L -12.0)*xinvsq*xinvsq/32.0);
      return (result);
    }
  // hack to catch small values
  if (x < 1.0e-10)
    {
      if (l==0)
	return (1.0);
      else
	return (0.0);
    }
    

  if (l==0)
    ErrorCode = gsl_sf_bessel_i0_scaled_e (x, &result);
  else if (l==1)
    ErrorCode = gsl_sf_bessel_i1_scaled_e (x, &result);
  else if (l==2)
    ErrorCode = gsl_sf_bessel_i2_scaled_e (x, &result);
  else
    ErrorCode = gsl_sf_bessel_il_scaled_e (l, x, &result);

  if (ErrorCode != GSL_SUCCESS)
    {
      fprintf (stderr,"Error in ModBessel.\n l = %d  x = %1.8e ErrCode = %d\n",
	       l, x, ErrorCode);
      exit(1);
    }
    return (//  exp(fabs(x)) *  
    result.val);
}
*/


double Log_il_Scaled (int l, double z)
{
  if (l==-1)
    {
      if (z > 25.0)
	{
	  return (- log(2.0*z));
	}
      else
	return (log (cosh(z)/z) - fabs(z));
    }
  else if (l==-2)
    {
      if (z > 30.0)
	{
	  return (log((1.0-(1.0/z))/(2.0*z)));
	}
      else
	{
	  double result = log (sinh (z)/z - cosh(z)/(z*z));
	  result -= fabs(z);
	  return (result);
	}
    }
  else
    {
      gsl_sf_result result;
      int ErrorCode;
      ErrorCode = gsl_sf_bessel_il_scaled_e (l, z, &result);
      if (ErrorCode != GSL_SUCCESS)
	{
	  fprintf (stderr, "Error in LogModBesselScaled.\n");
	  fprintf (stderr, "l = %d z = %1.8e ErrCode = %d\n",
		   l, z, ErrorCode);
	  exit(1);
	}
      return (log(result.val));
    }
}



double il(int l, double z)
{
  if (l == -1)
    {
      return (cosh(z)/z);
    }
  else if (l == -2)
    {
      return (sinh(z)/z - cosh(z)/(z*z));
    }
  else
    {
      gsl_sf_result result;
      int ErrorCode;
      ErrorCode = gsl_sf_bessel_il_scaled_e (l, z, &result);
      if (ErrorCode != GSL_SUCCESS)
	{
	  fprintf (stderr, "Error in il.\n");
	  fprintf (stderr, "l = %d z = %1.8e ErrCode = %d\n",
		   l, z, ErrorCode);
	  exit(1);
	}
      return (result.val*exp(fabs(z)));
    }
}

double jl(int l, double z)
{
  gsl_sf_result result;
  int ErrorCode;
  ErrorCode = gsl_sf_bessel_jl_e (l, z, &result);
  if (ErrorCode != GSL_SUCCESS)
    {
      fprintf (stderr, "Error in il.\n");
      fprintf (stderr, "l = %d z = %1.8e ErrCode = %d\n",
	       l, z, ErrorCode);
      exit(1);
    }
  return (result.val*exp(fabs(z)));
}



/*double il_scaled(int l, double z)
{
  return (sqrt(0.5*M_PI/z)*exp(LogRegModBesselScaled(l,z)));
}*/


const double ninv[] = 
  {0.000000000000000000e+00,  // This element is a dummy for 1.0/0.0;
   1.000000000000000000e+00,
   5.000000000000000000e-01,
   3.333333333333333333e-01,
   2.500000000000000000e-01,
   2.000000000000000000e-01,
   1.666666666666666667e-01,
   1.428571428571428571e-01,
   1.250000000000000000e-01,
   1.111111111111111111e-01,
   1.000000000000000000e-01,
   9.090909090909090910e-02,
   8.333333333333333333e-02,
   7.692307692307692308e-02,
   7.142857142857142857e-02,
   6.666666666666666667e-02,
   6.250000000000000000e-02,
   5.882352941176471e-02,
   5.555555555555555555e-02,
   5.263157894736842e-02,
   5.000000000000000e-02,
   4.761904761904762e-02,
   4.545454545454546e-02,
   4.347826086956522e-02,
   4.166666666666666e-02,
   4.000000000000000e-02,
   3.846153846153846e-02,
   3.703703703703703e-02,
   3.571428571428571e-02,
   3.448275862068965e-02,
   3.333333333333333e-02,
   3.225806451612903e-02,
   3.125000000000000e-02,
   3.030303030303030e-02,
   2.941176470588235e-02,
   2.857142857142857e-02,
   2.777777777777778e-02,
   2.702702702702703e-02,
   2.631578947368421e-02,
   2.564102564102564e-02,
   2.500000000000000e-02,
   2.439024390243903e-02,
   2.380952380952381e-02,
   2.325581395348837e-02,
   2.272727272727273e-02,
   2.222222222222222e-02,
   2.173913043478261e-02,
   2.127659574468085e-02,
   2.083333333333333e-02,
   2.040816326530612e-02,
   2.000000000000000e-02,
   1.960784313725490e-02,
   1.923076923076923e-02,
   1.886792452830189e-02,
   1.851851851851852e-02,
   1.818181818181818e-02,
   1.785714285714286e-02,
   1.754385964912281e-02,
   1.724137931034483e-02,
   1.694915254237288e-02,
   1.666666666666667e-02,
   1.639344262295082e-02,
   1.612903225806452e-02,
   1.587301587301587e-02,
   1.562500000000000e-02,
   1.538461538461539e-02,
   1.515151515151515e-02,
   1.492537313432836e-02,
   1.470588235294118e-02,
   1.449275362318841e-02,
   1.428571428571429e-02,
   1.408450704225352e-02,
   1.388888888888889e-02,
   1.369863013698630e-02,
   1.351351351351351e-02,
   1.333333333333333e-02,
   1.315789473684210e-02,
   1.298701298701299e-02,
   1.282051282051282e-02,
   1.265822784810127e-02,
   1.250000000000000e-02,
   1.234567901234568e-02,
   1.219512195121951e-02,
   1.204819277108434e-02,
   1.190476190476190e-02,
   1.176470588235294e-02,
   1.162790697674419e-02,
   1.149425287356322e-02,
   1.136363636363636e-02,
   1.123595505617977e-02,
   1.111111111111111e-02,
   1.098901098901099e-02,
   1.086956521739130e-02,
   1.075268817204301e-02,
   1.063829787234043e-02,
   1.052631578947368e-02,
   1.041666666666667e-02,
   1.030927835051546e-02,
   1.020408163265306e-02,
   1.010101010101010e-02,
   1.000000000000000e-02
  };

/*! Asymptotic expansion for il(z) for large z.  
    From Abramowitz and Stegun 9.7.1
   \[ i_l(z) \approx \frac{e^x}{2z} \left[1-\frac{\mu - 1}{8z} +
   \frac{(\mu -1)(\mu-9)}{2!(8z)^} -
   \frac{(\mu-1)(\mu-9)(\mu-25)}{3!(8z)^3} + \dots \] */
double il_scaled_series_large(int l, double z)
{
  double sl = (double)l;
  double nu = sl+0.5;
  const double mu = 4.0*nu*nu;
  double term = 1.0;
  double n = 0.0;
  double sum = 0.0;
  int nint = 0;
  double zinv = 1.0/z;

  while ((fabs(term) > (1e-20*(sum+term))) && (n<50.0))
    {
      sum += term;
      n+=1.0;
      nint++;
      double a = 2.0*n-1.0;
      term *= -(mu - a*a);
      term /= (8.0*z*n);
      //term *= 0.125*zinv*ninv[nint];
    }
  if (n>=50.0)
    cerr << "Warning:  Unconverged in il_scaled_series_large.\n";
  //cerr << "n = " << n << endl;
  return (0.5*zinv*sum);
}


/*! Asymptotic expansion for il(z) for small z.  
    From Abramowitz and Stegun 9.7.1
   \[ i_l(z) \approx \frac{z^l}{(2l+1)!! 
\left[ 1 + \frac{\frac{z^2}{2}}{1!(2l+3)} +
\frac{\left( \frac{z^2}{2}\right)^2}{2!(2l+3)(2l+5)+\dots} \right] \]
*/

double il_scaled_series_small(int l, double z)
{
  double sl=(double)l;
  double doubleFactorial = 1.0;
  for (double x=3.0; x<=(2.0*sl+1.0); x+=2.0)
    doubleFactorial *=x;

  double mu = 0.5*z*z;
  double sum = 0.0;
  double term = 1.0;
  double n = 1.0;
  while ((term > (1e-20*(term+sum))) && (n < 50.0))
    {
      sum += term;
      term *= mu /(n*(2.0*sl+2.0*n+1.0));
      n+=1.0;
    }
  return (pow(z,sl)/doubleFactorial * sum*exp(-z));
}



double log_il_scaled_series_small(int l, double z)
{
  double sl=(double)l;
  double doubleFactorial = 1.0;
  for (double x=3.0; x<=(2.0*sl+1.0); x+=2.0)
    doubleFactorial *=x;

  double mu = 0.5*z*z;
  double sum = 0.0;
  double term = 1.0;
  double n = 1.0;
  while ((term > (1e-20*(term+sum))) && (n < 50.0))
    {
      sum += term;
      term *= mu /(n*(2.0*sl+2.0*n+1.0));
      n+=1.0;
    }
  if (n>=50.0)
    cerr << "Warning:  Unconverged in log_il_scaled_series_small.\n";
  return (sl*log(z)+log(sum/doubleFactorial)-z);
}




double il_scaled(int l, double z)
{
  double zinv = 1.0/z;
  /* if (l == 0)
    return (0.5*(1.0-exp(-2.0*z))/z);
  else if (l==1)
    return (0.5*zinv*(1.0-zinv + (1.0+zinv)*exp(-2.0*z)));
    else */
  if (l == -1)
    return ((1.0+exp(-2.0*z))/(2.0*z));
  else if (l == -2)
    {
      double exp_neg_2z = exp (-2.0*z);
      double zinv = 1.0/z;
      return (0.5*(1.0-exp_neg_2z)*zinv - 0.5*(1.0+exp_neg_2z)*zinv*zinv);
    }
  else if (z > 1000.0)
    return il_scaled_series_large(l,z);
  else if (z == 0.0)
    return (0.0);
  else
    {
      gsl_sf_result result;
      int ErrorCode;
      ErrorCode = gsl_sf_bessel_il_scaled_e (l, z, &result);
      if (ErrorCode != GSL_SUCCESS)
	{
	  fprintf (stderr, "Error in il_scaled.\n");
	  fprintf (stderr, "l = %d z = %1.8e ErrCode = %d\n",
		   l, z, ErrorCode);
	  exit(1);
	}
      return (result.val);
    }
}



inline static double 
debye_u1(const double * tpow)
{
  return (3.0*tpow[1] - 5.0*tpow[3])/24.0;
}

inline static double 
debye_u2(const double * tpow)
{
  return (81.0*tpow[2] - 462.0*tpow[4] + 385.0*tpow[6])/1152.0;
}

inline
static double debye_u3(const double * tpow)
{
  return (30375.0*tpow[3] - 369603.0*tpow[5] + 765765.0*tpow[7] - 425425.0*tpow[9])/414720.0;
}

inline
static double debye_u4(const double * tpow)
{
  return (4465125.0*tpow[4] - 94121676.0*tpow[6] + 349922430.0*tpow[8] - 
          446185740.0*tpow[10] + 185910725.0*tpow[12])/39813120.0;
}

inline
static double debye_u5(const double * tpow)
{
  return (1519035525.0*tpow[5]     - 49286948607.0*tpow[7] + 
          284499769554.0*tpow[9]   - 614135872350.0*tpow[11] + 
          566098157625.0*tpow[13]  - 188699385875.0*tpow[15])/6688604160.0;
}

#if 0
inline
static double debye_u6(const double * tpow)
{
  return (2757049477875.0*tpow[6] - 127577298354750.0*tpow[8] + 
          1050760774457901.0*tpow[10] - 3369032068261860.0*tpow[12] + 
          5104696716244125.0*tpow[14] - 3685299006138750.0*tpow[16] + 
          1023694168371875.0*tpow[18])/4815794995200.0;
}
#endif


/* nu -> Inf; uniform in x > 0  [Abramowitz+Stegun, 9.7.7]
 *
 * error:
 *   The error has the form u_N(t)/nu^N  where  0 <= t <= 1.
 *   It is not hard to show that |u_N(t)| is small for such t.
 *   We have N=6 here, and |u_6(t)| < 0.025, so the error is clearly
 *   bounded by 0.025/nu^6. This gives the asymptotic bound on nu
 *   seen below as nu ~ 100. For general MACH_EPS it will be 
 *                     nu > 0.5 / MACH_EPS^(1/6)
 *   When t is small, the bound is even better because |u_N(t)| vanishes
 *   as t->0. In fact u_N(t) ~ C t^N as t->0, with C ~= 0.1.
 *   We write
 *                     err_N <= min(0.025, C(1/(1+(x/nu)^2))^3) / nu^6
 *   therefore
 *                     min(0.29/nu^2, 0.5/(nu^2+x^2)) < MACH_EPS^{1/3}
 *   and this is the general form.
 *
 * empirical error analysis, assuming 14 digit requirement:
 *   choose   x > 50.000 nu   ==>  nu >   3
 *   choose   x > 10.000 nu   ==>  nu >  15
 *   choose   x >  2.000 nu   ==>  nu >  50
 *   choose   x >  1.000 nu   ==>  nu >  75
 *   choose   x >  0.500 nu   ==>  nu >  80
 *   choose   x >  0.100 nu   ==>  nu >  83
 *
 * This makes sense. For x << nu, the error will be of the form u_N(1)/nu^N,
 * since the polynomial term will be evaluated near t=1, so the bound
 * on nu will become constant for small x. Furthermore, increasing x with
 * nu fixed will decrease the error.
 */
double log_il_scaled_asymp_unif_e(int l, const double x)
{
  double nu = (double)l + 0.5;
  double result;
  int i;
  double z = x/nu;
  double root_term = sqrt(1.0 + z*z);
  double pre = 1.0/sqrt(2.0*M_PI*nu * root_term);
  double eta = root_term + log(z/(1.0+root_term));
  double ex_arg = 
    ( z < 1.0/GSL_ROOT3_DBL_EPSILON ? nu*(-z + eta) : 
      -0.5*nu/z*(1.0 - 1.0/(12.0*z*z)) );
  gsl_sf_result ex_result;
  //int stat_ex = gsl_sf_exp_e(ex_arg, &ex_result);
  double t = 1.0/root_term;
  double sum;
  double tpow[16];
  tpow[0] = 1.0;
  for(i=1; i<16; i++) tpow[i] = t * tpow[i-1];
  sum = 1.0 + debye_u1(tpow)/nu + debye_u2(tpow)/(nu*nu) + 
    debye_u3(tpow)/(nu*nu*nu) + 
    debye_u4(tpow)/(nu*nu*nu*nu) + 
    debye_u5(tpow)/(nu*nu*nu*nu*nu);
  
  result = log(pre*sum) + ex_arg +0.5*log(M_PI/(2.0*x));
    //    result->val  = pre * ex_result.val * sum;
//   result->err  = pre * ex_result.val / (nu*nu*nu*nu*nu*nu);
//   result->err += pre * ex_result.err * fabs(sum);
//   result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
  return result;
}





double log_il_scaled(int l, double z)
{
  double sl = (double)l;
  double lp1 = sl+1.0;
  double approxlog = sl*log(z) - lp1*(log(2.0*lp1) -  1.0);
  /// Check for underflow;
  if (approxlog < -300.0)
    return (log_il_scaled_series_small(l,z));
  else if ((l> 83) && (z > 8.3))
    return (log_il_scaled_asymp_unif_e(l,z));
  else 
    {
      double il_val = il_scaled(l,z);
      if (il_val <= 0.0)
	return (-400.0);
      else
	return (log(il_val));
    }
}



void Test_il_series_large()
{
  FILE *fout;
  
  assert((fout=fopen("il_series_large.dat", "w"))!=NULL);
  
  int l = 100;

  for (int i=0; i<5000; i++)
    {
      double z = 300.0 + 20.0*i;
      double gslval = il_scaled(l, z);
      double seriesval = il_scaled_series_large(l, z);
      fprintf (fout, "%1.16e %1.16e %1.16e\n", z, gslval, seriesval);
    }
  fclose(fout);
}



void Test_il_series_small()
{
  FILE *fout;
  
  assert((fout=fopen("il_series_small.dat", "w"))!=NULL);
  
  int l = 100;

  for (int i=0; i<5000; i++)
    {
      double z = 1e-3*i;
      double gslval = (log_il_scaled(l, z));
      double seriesval = (log_il_scaled_series_small(l, z));
      fprintf (fout, "%1.16e %1.16e %1.16e\n", z, gslval, seriesval);
    }
  fclose(fout);
}




double dil_dz (int l, double z)
{
  //return (il(l-1,z) - (1.0+l)/z * il(l,z));
  return (il(l+1,z) + (double)l/z * il(l,z));
}

double dil_dz_scaled (int l, double z)
{
  //return (il_scaled(l-1,z) - (1.0+l)/z * il_scaled(l,z));
  //return (il_scaled(l+1,z) + (double)l/z * il_scaled(l,z));
  return (exp(log_il_scaled(l+1,z)) + 
	  (double)l/z * exp(log_il_scaled(l,z)));
}


double dil_dz_FD (int l, double z)
{
  double delta = 1.0e-6;
  if (z < delta)
    delta = 0.00001*z;
  return (1.0/(2.0*delta) * (il(l, z+delta) - il(l,z-delta)));
}

double d2il_dz2 (int l, double z)
{
  return (il(l-2,z) - (2.0*l+1.0)/z*il(l-1,z) + (1.0+l)*(2.0+l)/(z*z)*il(l,z));
}


double d2il_dz2_scaled (int l, double z)
{
//   return (il_scaled(l-2,z) - (2.0*l+1.0)/z*il_scaled(l-1,z) 
// 	  + (1.0+l)*(2.0+l)/(z*z)*il_scaled(l,z));
  return (exp(log_il_scaled(l-2,z)) - (2.0*l+1.0)/z*exp(log_il_scaled(l-1,z))
	  + (1.0+l)*(2.0+l)/(z*z)*exp(log_il_scaled(l,z)));
}


double d2il_dz2_FD (int l, double z)
{
  double delta = 1.0e-4;
  double il0, ilplus, ilminus;
  il0     = il(l, z      );
  ilplus  = il(l, z+delta);
  ilminus = il(l, z-delta);
  return (1.0/(delta*delta) * (ilplus - 2.0*il0 + ilminus));
}



void Test_il_Derivs()
{
  for (int l=0; l<10; l++)
    {
      for (double z=0.000001; z<50.0; z+=0.01)
	{
	  double dil = dil_dz(l, z);
	  double dil_FD = dil_dz_FD(l, z);

	  double d2il = d2il_dz2(l, z);
	  double d2il_FD = d2il_dz2_FD(l, z);

	  double err1 = fabs(dil_FD - dil) / dil;
	  double err2 = fabs(d2il_FD - d2il) / d2il;
	  fprintf (stderr, "%d %1.6f %1.6e %1.6e\n", l, z, err1, err2);
	}
    }
}


double LogRegModBesselScaled(int l, double z)
{
  gsl_sf_result result;
  int ErrorCode;
  double gsl_val, series_val;
  double f;
  // f will be used to make the transition from one representation to
  // the other continuous.
  double l_log10z = -log10(z) * l;
  double retval;

  if (l_log10z > 200.0)
    {
      // Compute series expression:  Abramowitz & Stegun 10.2.5
      double denom=1.0;
      series_val = 0.0;
      
      for (int i=3; i<=(2*l+1); i+=2)
	{
	  if (denom > 1.0e200)
	    {
	      series_val -= log(denom);
	      denom = 1.0;
	    }
	  denom *= (double)i;
	}
      series_val -= log(denom);
      series_val += ((double)l+0.5)*log(z); 
      series_val -= 0.5 * log (0.5 * M_PI);
      series_val += log (1.0 + 0.5*z*z/( 1.0+(2.0*l+3.0))
			 + 0.125*z*z*z*z/((2.0*l+3.0)*(2.0*l+5.0)));
      series_val -= fabs(z);
      retval = series_val;
    }
  else  // Compute gsl value
    {
      ErrorCode = gsl_sf_bessel_Inu_scaled_e ((double)l+0.5, z, &result);
      //ErrorCode = gsl_sf_bessel_il_scaled_e (l, z, &result);
      if (ErrorCode != GSL_SUCCESS)
	{
	  fprintf (stderr, "Error in LogModBesselScaled.\n");
	  fprintf (stderr, "l = %d z = %1.8e ErrCode = %d\n",
		   l, z, ErrorCode);
	  exit(1);
	}
      gsl_val = log(result.val);
      //gsl_val += 0.5*log(2.0*z/M_PI);
      // Interpolate between series val and gsl_val;
      //double sigma = 1.5;
      //double diff = (100.0 - l_log10z);
      //f = exp(-(diff*diff)/(sigma*sigma));
      //f = 1.0;
      //result = f*gsl_val + (1.0-f)*series_val;
      retval = gsl_val;
    }
  return (retval);
}


double RegModBesselFracScaled(double nu, double x)
{
  gsl_sf_result result;
  int ErrorCode;


  // hack to catch small values
  double exponent;
  // This means the value is < 10^-(200)
  if ((-log10(x) * nu) > 200.0)
    {
      if (fabs(nu-0.5) < 1e-5)
	return (1.0);
      else
	return (1.0e-200);
    }
    
  ErrorCode = gsl_sf_bessel_Inu_scaled_e (nu, x, &result);

  if (ErrorCode != GSL_SUCCESS)
    {
      fprintf (stderr,"Error in ModBesselFracScaled.\n nu = %1.8e x = %1.8e ErrCode = %d\n",
	       nu, x, ErrorCode);
      exit(1);
    }
  return (/* exp(fabs(x)) * */ result.val);
}

/*
double RegModBesselFrac(double nu, double x)
{
  gsl_sf_result result;
  int ErrorCode;

  ErrorCode = gsl_sf_bessel_Inu_e (nu, x, &result);

  if (ErrorCode != GSL_SUCCESS)
    {
      fprintf (stderr,"Error in ModBessel.\n l = %1.12e  x = %1.8e ErrCode = %d\n",
	       nu, x, ErrorCode);
      exit(1);
    }
  return (result.val);
}
*/
double Legendre (int l, double CosTheta)
{
  gsl_sf_result result;

  gsl_sf_legendre_Pl_e (l, CosTheta, &result);

  return (result.val);
}





