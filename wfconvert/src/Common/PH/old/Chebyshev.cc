#include "Chebyshev.h"

scalar Chebyshev (Array<scalar,1> a, scalar x)
{
  int m = a.rows()-1;
  scalar R_jplus2 = 0.0;
  scalar R_jplus1 = 0.0;
  scalar R_j;

  for (int j=m; j>0; j--)
    {
      R_j = 2.0*x*R_jplus1-R_jplus2+a(j);
      R_jplus2 = R_jplus1;
      R_jplus1 = R_j;
    }

  return (x*R_j-R_jplus2 + 0.5*a(0));
}



scalar Chebyshev2 (Array<scalar,1> a, scalar x)
{
  scalar sum = 0.5*a(0);
  scalar theta = acos(x);

  for (int j=1; j<a.rows(); j++)
    sum += a(j) * cos((scalar)j * theta);
  return (sum);
}


Array<scalar,1> ChebyshevDeriv (Array<scalar,1> a)
{
  int m = a.rows()-1;
  if (m<0)
    m=0;
  Array<scalar,1> deriv(m);
  

  if (m>0)
    deriv(m-1) = 2.0*(scalar)m * a(m);
  if (m>1)
    deriv(m-2) = 2.0*(scalar)(m-1) *a(m-1);

  for (int i=m-2; i>=1; i--)
    deriv(i-1) = deriv(i+1) + 2.0*(scalar)i*a(i);
  return (deriv);
}


/*main()
{
  Array<scalar,1> a(6);
  a = 1.0, 2.1, 1.4, 0.8, -2.3, 6.7;
  //a = 0.0, 1.0, 0.0, 0.0, 0.0, 0.0;
  Array<scalar,1> da(5);
  
  da = ChebyshevDeriv(a);
  const scalar delta = 1.0e-6;

  for (scalar x=-1.0; x<=1.0; x+=0.01)
    {
      scalar c1 = Chebyshev(a,x);
      scalar c2 = Chebyshev2(a,x);
      scalar cplus = Chebyshev(a,x+delta);
      scalar cminus = Chebyshev(a,x-delta);
      scalar FDderiv = (cplus-cminus)/(2.0*delta);
      scalar deriv = Chebyshev(da, x);
      //fprintf (stdout, "%21.16f %21.16f\n", c1, c2);
      fprintf (stdout, "%21.16f %21.16f\n", deriv, FDderiv);
    }
}*/
