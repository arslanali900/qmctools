#include "AtomicOrbital.h"
#include <einspline/multi_bspline.h>
#include <Common/Fitting/Fitting.h>
#include <gsl/gsl_sf_bessel.h>
#include <OrbitalClass.h>

inline int CubicRoots (double a, double b, double c, double d,
		       double &x1, double &x2, double &x3) {
   double A = b/a;
   double B = c/a;
   double C = d/a;
   double Q = (A*A - 3.0*B)/9.0;
   double R = (2.0*A*A*A - 9.0*A*B + 27.0*C)/54.0;
   //cerr << "Q = " << Q << " R = " << R << "\n";
   if ((R*R) < (Q*Q*Q)) {
      double theta = acos(R/sqrt(Q*Q*Q));
      double twosqrtQ = 2.0*sqrt(Q);
      double third = 1.0/3.0;
      double thirdA = third * A;
      x1 = -twosqrtQ*cos(third*theta) - thirdA;
      x2 = -twosqrtQ*cos(third*(theta + 2.0*M_PI)) - thirdA;
      x3 = -twosqrtQ*cos(third*(theta - 2.0*M_PI)) - thirdA;
      if (x2 < x3) 
	swap (x2, x3);
      if (x1 < x3)
	swap (x1, x3);
      if (x1 < x2)
	swap (x1, x2);
      
      assert (x1 > x2);
      assert (x2 > x3);
      assert (x3 > 0.0);
      
      cerr << "x1 = " << x1 << endl;
      cerr << "x2 = " << x2 << endl;
      cerr << "x3 = " << x3 << endl;
      
      return 3;
   }
   else {
      double D = -Q*Q*Q + R*R;
      double u = cbrt(-R + sqrt(D));
      double v = cbrt(-R - sqrt(D));
      double y1 = u+v;
      x1 = y1 - A/3.0;
      return 1;
   }
}


inline void
Ylm (TinyVector<double,3> rhat, int lMax,
     vector<complex<double> >& YlmVec) {
   YlmVec.resize(5*(lMax+1)*(lMax+1));
   
   if (fabs(rhat[0]) < 1.0e-10 &&
       fabs(rhat[1]) < 1.0e-10 &&
       fabs(rhat[2]- 1.0) < 1.0e-10)
     rhat = Vec3 (1.0e-6, 1.0e-6, 0.999999999999);
   
   if (fabs(rhat[0]) < 1.0e-10 &&
       fabs(rhat[1]) < 1.0e-10 &&
       fabs(rhat[2]+ 1.0) < 1.0e-10)
     rhat = Vec3 (1.0e-6, 1.0e-6, -0.999999999999);
   
   
   const double fourPiInv = 0.0795774715459477;
   
   double costheta = rhat[2];
   double sintheta = std::sqrt(1.0-costheta*costheta);
   double cottheta = costheta/sintheta;
    
   double cosphi, sinphi;
   cosphi=rhat[0]/sintheta;
   sinphi=rhat[1]/sintheta;
   
   complex<double> e2iphi(cosphi, sinphi);
   
   double lsign = 1.0;
   double dl = 0.0;
   double XlmVec[2*lMax+1], dXlmVec[2*lMax+1];
   for (int l=0; l<=lMax; l++) {
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
      complex<double> e2imphi (1.0, 0.0);
      for (int m=0; m<=l; m++) {
	 YlmVec[l*(l+1)+m]  =  XlmVec[l+m]*e2imphi;
	 YlmVec[l*(l+1)-m]  =  XlmVec[l-m]*conj(e2imphi);
	 e2imphi *= e2iphi;
      } 
      dl += 1.0;
      lsign *= -1.0;
   }
}


void 
QuadratureRule::Check() {
   vector<vector<complex<double> > > Ylms(rhat.size());
   FILE *fout = fopen ("QuadRule.dat", "w");
   for (int i=0; i<rhat.size(); i++) {
      fprintf (fout, "%16.12f %16.12f %16.12f\n", rhat[i][0], 
	       rhat[i][1], rhat[i][2]);
      Ylm(rhat[i], lexact, Ylms[i]);
   }
   fclose(fout);
   
   for (int l1=0; l1<=lexact; l1++) 
     for (int l2=0; l2 <= (lexact-l1); l2++) 
       for (int m1=-l1; m1<=l1; m1++)
	 for (int m2=-l2; m2<=l2; m2++) {
	    complex<double> sum(0.0, 0.0);
	    int lm1 = l1*(l1+1) + m1;
	    int lm2 = l2*(l2+1) + m2;
	    for (int k=0; k<rhat.size(); k++) {
	       complex<double> v1 = Ylms[k][lm1];
	       complex<double> v2 = Ylms[k][lm2];
	       sum += 4.0*M_PI*weight[k] * conj(v1)*v2;
	    }
	    double re = real (sum);
	    double im = imag (sum);
	    if ((l1==l2) && (m1==m2)) 
	      re -= 1.0;
	    if ((std::fabs(im) > 1.0e-10) || (std::fabs(re) > 1.0e-10)) {
	       cerr << "Broken spherical quadrature for " << rhat.size() << "-point rule.\n" << endl;
	       fprintf (stderr, "(l1,m1,l2,m2) = (%2d,%2d,%2d,%2d)  sum = (%20.16f %20.16f)\n",
			l1, m1, l2, m2, real(sum), imag(sum));
	       abort();
	    }
	    
	 }
}

void
QuadratureRule::MakeLebedev(int lexact) {
   if (lexact != 23) {
      cerr << "Rule for lexact = " << lexact << " not implemented.\n";
      abort();
   }
   double A[3];
   A[0] = pow(2.0,7.0)*73.0/5242545.0;
   A[1] = pow(2.0,14.0)*1663.0/((3.0*3.0*3.0)*11.0*11.0*1458821.0);
   A[2] = pow(3.0, 10.0)*1599797/(pow(2.0,9.0)*173.0*173.0*13.0*13.0*6545.0);
   // fprintf (stderr, "A[0] = %1.14e\n", A[0]);  
   // fprintf (stderr, "A[1] = %1.14e\n", A[1]);
   // fprintf (stderr, "A[2] = %1.14e\n", A[2]);
   double B[4];
   B[0] = 5.51877146727e-3;
   B[1] = 5.15823771181e-3;
   B[2] = 5.60870408259e-3;
   B[3] = 4.10677702817e-3;
   double m[4], l[4];
   m[0] = 0.777493219315;
   m[1] = 0.912509096867;
   m[2] = 0.314196994183;
   m[3] = 0.982972302707;
   l[0] = 1.0/sqrt(2.0) * sqrt(1.0-m[0]*m[0]);  //0.444593317871;
   l[1] = 1.0/sqrt(2.0) * sqrt(1.0-m[1]*m[1]);  //0.289246452758;
   l[2] = 1.0/sqrt(2.0) * sqrt(1.0-m[2]*m[2]);  //0.671297344270;
   l[3] = 1.0/sqrt(2.0) * sqrt(1.0-m[3]*m[3]);  //0.129933544765;
   // fprintf (stderr, "l[0] = %1.14f\n", l[0]);
   // fprintf (stderr, "l[1] = %1.14f\n", l[1]);
   // fprintf (stderr, "l[2] = %1.14f\n", l[2]);
   // fprintf (stderr, "l[3] = %1.14f\n", l[3]);
   double C1 = pow(38.0, 4.0)/(33.0*33.0*7.0*7.0*7.0*1105.0);
   double p1 = 0.938319218138;
   double q1 = 0.345770219761;
   double D1 = 19.0*19.0*pow(23.0, 6.0)/(7.0*7.0*7.0*pow(2.0,9.0)*3.0*3.0*6113965.0);
   double r1 = 0.836036015482;
   double u1 = 0.159041710538;
   double w1 = 0.525118572443;
   // fprintf (stderr, "C1 = %1.16e\n", C1);
   // fprintf (stderr, "D1 = %1.16e\n", D1);

   vector<Vec3> a[3], b[4], c, d;
   double p = 1.0/std::sqrt(2.0);
   double q = 1.0/std::sqrt(3.0);
   
   a[0].push_back(Vec3( 1.0, 0.0, 0.0));
   a[0].push_back(Vec3(-1.0, 0.0, 0.0));
   a[0].push_back(Vec3( 0.0, 1.0, 0.0));
   a[0].push_back(Vec3( 0.0,-1.0, 0.0));
   a[0].push_back(Vec3( 0.0, 0.0, 1.0));
   a[0].push_back(Vec3( 0.0, 0.0,-1.0));
   
   a[1].push_back(Vec3(   p,   p, 0.0));
   a[1].push_back(Vec3(   p,  -p, 0.0));
   a[1].push_back(Vec3(  -p,   p, 0.0));
   a[1].push_back(Vec3(  -p,  -p, 0.0));
   a[1].push_back(Vec3(   p, 0.0,   p));
   a[1].push_back(Vec3(   p, 0.0,  -p));
   a[1].push_back(Vec3(  -p, 0.0,   p));
   a[1].push_back(Vec3(  -p, 0.0,  -p));
   a[1].push_back(Vec3( 0.0,   p,   p));
   a[1].push_back(Vec3( 0.0,   p,  -p));
   a[1].push_back(Vec3( 0.0,  -p,   p));
   a[1].push_back(Vec3( 0.0,  -p,  -p));
      
   a[2].push_back(Vec3(   q,   q,   q));
   a[2].push_back(Vec3(   q,   q,  -q));
   a[2].push_back(Vec3(   q,  -q,   q));
   a[2].push_back(Vec3(   q,  -q,  -q));
   a[2].push_back(Vec3(  -q,   q,   q));
   a[2].push_back(Vec3(  -q,   q,  -q));
   a[2].push_back(Vec3(  -q,  -q,   q));
   a[2].push_back(Vec3(  -q,  -q,  -q));
   
   for (int i=0; i<4; i++) {
      b[i].push_back(Vec3(  l[i],  l[i],  m[i]));
      b[i].push_back(Vec3(  l[i],  l[i], -m[i]));
      b[i].push_back(Vec3(  l[i], -l[i],  m[i]));
      b[i].push_back(Vec3(  l[i], -l[i], -m[i]));
      b[i].push_back(Vec3( -l[i],  l[i],  m[i]));
      b[i].push_back(Vec3( -l[i],  l[i], -m[i]));
      b[i].push_back(Vec3( -l[i], -l[i],  m[i]));
      b[i].push_back(Vec3( -l[i], -l[i], -m[i]));
      
      b[i].push_back(Vec3(  l[i],  m[i],  l[i]));
      b[i].push_back(Vec3(  l[i],  m[i], -l[i]));
      b[i].push_back(Vec3(  l[i], -m[i],  l[i]));
      b[i].push_back(Vec3(  l[i], -m[i], -l[i]));
      b[i].push_back(Vec3( -l[i],  m[i],  l[i]));
      b[i].push_back(Vec3( -l[i],  m[i], -l[i]));
      b[i].push_back(Vec3( -l[i], -m[i],  l[i]));
      b[i].push_back(Vec3( -l[i], -m[i], -l[i]));
      
      b[i].push_back(Vec3(  m[i],  l[i],  l[i]));
      b[i].push_back(Vec3(  m[i],  l[i], -l[i]));
      b[i].push_back(Vec3(  m[i], -l[i],  l[i]));
      b[i].push_back(Vec3(  m[i], -l[i], -l[i]));
      b[i].push_back(Vec3( -m[i],  l[i],  l[i]));
      b[i].push_back(Vec3( -m[i],  l[i], -l[i]));
      b[i].push_back(Vec3( -m[i], -l[i],  l[i]));
      b[i].push_back(Vec3( -m[i], -l[i], -l[i]));
   }
   
   c.push_back( Vec3(  p1,  q1,  0.0));
   c.push_back( Vec3(  p1, -q1,  0.0));
   c.push_back( Vec3( -p1,  q1,  0.0));
   c.push_back( Vec3( -p1, -q1,  0.0));
   
   c.push_back( Vec3(  p1,  0.0,  q1));
   c.push_back( Vec3(  p1,  0.0, -q1));
   c.push_back( Vec3( -p1,  0.0,  q1));
   c.push_back( Vec3( -p1,  0.0, -q1));
   
   c.push_back( Vec3( 0.0,   p1,  q1));
   c.push_back( Vec3( 0.0,   p1, -q1));
   c.push_back( Vec3( 0.0,  -p1,  q1));
   c.push_back( Vec3( 0.0,  -p1, -q1));
   
   c.push_back( Vec3(  q1,  p1,  0.0));
   c.push_back( Vec3(  q1, -p1,  0.0));
   c.push_back( Vec3( -q1,  p1,  0.0));
   c.push_back( Vec3( -q1, -p1,  0.0));
   
   c.push_back( Vec3(  q1,  0.0,  p1));
   c.push_back( Vec3(  q1,  0.0, -p1));
   c.push_back( Vec3( -q1,  0.0,  p1));
   c.push_back( Vec3( -q1,  0.0, -p1));
   
   c.push_back( Vec3( 0.0,   q1,  p1));
   c.push_back( Vec3( 0.0,   q1, -p1));
   c.push_back( Vec3( 0.0,  -q1,  p1));
   c.push_back( Vec3( 0.0,  -q1, -p1));
   
   
   d.push_back(Vec3(  r1,  u1,  w1)); 
   d.push_back(Vec3(  r1,  u1, -w1)); 
   d.push_back(Vec3(  r1, -u1,  w1)); 
   d.push_back(Vec3(  r1, -u1, -w1)); 
   d.push_back(Vec3( -r1,  u1,  w1)); 
   d.push_back(Vec3( -r1,  u1, -w1)); 
   d.push_back(Vec3( -r1, -u1,  w1)); 
   d.push_back(Vec3( -r1, -u1, -w1)); 
   
   d.push_back(Vec3(  r1,  w1,  u1)); 
   d.push_back(Vec3(  r1,  w1, -u1)); 
   d.push_back(Vec3(  r1, -w1,  u1)); 
   d.push_back(Vec3(  r1, -w1, -u1)); 
   d.push_back(Vec3( -r1,  w1,  u1)); 
   d.push_back(Vec3( -r1,  w1, -u1)); 
   d.push_back(Vec3( -r1, -w1,  u1)); 
   d.push_back(Vec3( -r1, -w1, -u1)); 
   
   d.push_back(Vec3(  u1,  r1,  w1)); 
   d.push_back(Vec3(  u1,  r1, -w1)); 
   d.push_back(Vec3(  u1, -r1,  w1)); 
   d.push_back(Vec3(  u1, -r1, -w1)); 
   d.push_back(Vec3( -u1,  r1,  w1)); 
   d.push_back(Vec3( -u1,  r1, -w1)); 
   d.push_back(Vec3( -u1, -r1,  w1)); 
   d.push_back(Vec3( -u1, -r1, -w1)); 
   
   d.push_back(Vec3(  u1,  w1,  r1)); 
   d.push_back(Vec3(  u1,  w1, -r1)); 
   d.push_back(Vec3(  u1, -w1,  r1)); 
   d.push_back(Vec3(  u1, -w1, -r1)); 
   d.push_back(Vec3( -u1,  w1,  r1)); 
   d.push_back(Vec3( -u1,  w1, -r1)); 
   d.push_back(Vec3( -u1, -w1,  r1)); 
   d.push_back(Vec3( -u1, -w1, -r1)); 
   
   d.push_back(Vec3(  w1,  u1,  r1)); 
   d.push_back(Vec3(  w1,  u1, -r1)); 
   d.push_back(Vec3(  w1, -u1,  r1)); 
   d.push_back(Vec3(  w1, -u1, -r1)); 
   d.push_back(Vec3( -w1,  u1,  r1)); 
   d.push_back(Vec3( -w1,  u1, -r1)); 
   d.push_back(Vec3( -w1, -u1,  r1)); 
   d.push_back(Vec3( -w1, -u1, -r1)); 
   
   d.push_back(Vec3(  w1,  r1,  u1)); 
   d.push_back(Vec3(  w1,  r1, -u1)); 
   d.push_back(Vec3(  w1, -r1,  u1)); 
   d.push_back(Vec3(  w1, -r1, -u1)); 
   d.push_back(Vec3( -w1,  r1,  u1)); 
   d.push_back(Vec3( -w1,  r1, -u1)); 
   d.push_back(Vec3( -w1, -r1,  u1)); 
   d.push_back(Vec3( -w1, -r1, -u1)); 
   
   for (int k=0; k<3; k++)
     if (std::fabs(A[k]) > 1.0e-10) 
       for (int i=0; i<a[k].size(); i++)
	 addknot(a[k][i], A[k]);
   for (int k=0; k<4; k++)
     if (std::fabs(B[k]) > 1.0e-10) 
       for (int i=0; i<b[k].size(); i++)
	 addknot(b[k][i], B[k]);
   for (int i=0; i<c.size(); i++) 
     addknot(c[i], C1);
   
   for (int i=0; i<d.size(); i++) 
     addknot(d[i], D1);
   
   double wSum = 0.0;
   for (int i=0; i<weight.size(); i++)
     wSum += weight[i];
   for (int i=0; i<weight.size(); i++)
     weight[i] /= wSum;
   
   cerr << "Quadrature rule has " << rhat.size() << " points.\n";
}


void
QuadratureRule::SetRule (int rule) {
   int nk;
   double w;
   typedef enum {SINGLE, TETRA, OCTA, ICOSA} SymmType;
   SymmType symmetry;
   double A, B, C, D;
   A = B = C = D = 0.0;
   
   double p = 1.0/std::sqrt(2.0);
   double q = 1.0/std::sqrt(3.0);
   double r = 1.0/std::sqrt(11.0);
   double s = 3.0/std::sqrt(11.0);
   double u=0.0, v=0.0;
   
   // double tau1 = 651.0/ 391.0;
   // double tau2 = (986.0*tau1/31.0 - 53.0/17.0)/3.0;
   // double tau3 = ( 17.0*tau1/31.0 - 14.0/17.0)/3.0;
   // double t1, t2, t3;
   // CubicRoots (1.0, -tau1, tau2, -tau3, t1, t2, t3);
   // cerr << "t1 = " << t1 << endl;
   // assert (t1>t2  && t2>t3 && t3>0.0);
   // 
   // 
   switch (rule) {
    case 1:
      nk = 1;   symmetry = SINGLE; lexact = 0;
      A = 1.0;
      break;
    case 2:
      nk = 4;   symmetry = TETRA;  lexact = 2;
      A=0.25;
      break;
    case 3:
      nk = 6;   symmetry = OCTA;   lexact = 3;
      A=1.0/6.0;
      break;
    case 4:
      nk = 12;  symmetry = ICOSA;  lexact = 5;
      A = 1.0/12.0;
      B = 1.0/12.0;
      break;
    case 5:
    nk = 18;  symmetry = OCTA;   lexact = 5;
      A = 1.0/30.0; 
      B = 1.0/15.0;
      break;
    case 6:
      nk = 26;  symmetry = OCTA;   lexact = 7;
      A = 1.0  / 21.0;
      B = 4.0  / 105.0;
      C = 27.0 / 840.0;
    break;
    case 7:
      nk = 50;  symmetry = OCTA;   lexact = 11;
      A = 4.0/315.0;
      B = 64.0/2835.0;
      C = 27.0/1280.0;
      D = 14641.0/725760.0;
      break;
    case 8:
      nk = 194; symmetry = OCTA; lexact = 23;
      // fprintf (stderr, "m1 = %1.16f\n", sqrt(t1));
      // fprintf (stderr, "m2 = %1.16f\n", sqrt(t2));
      // fprintf (stderr, "m3 = %1.16f\n", sqrt(t3));
      break;
    default:
      cerr << "Unrecognized spherical quadrature rule " << rule << ".";
    abort();
   }
   
   if (lexact == 23)
     MakeLebedev (lexact);
   else {
      
      // First, build a_i, b_i, and c_i points
      vector<Vec3> a, b, c, d, e;
      
      if (symmetry == SINGLE) {
	 a.push_back (Vec3(1.0, 0.0, 0.0));
      }
      else if (symmetry == TETRA) {
	 a.push_back(Vec3( q, q, q));
	 a.push_back(Vec3( q,-q,-q));
	 a.push_back(Vec3(-q, q,-q));
	 a.push_back(Vec3(-q,-q, q));
      }
      else if (symmetry == OCTA) {
	 a.push_back(Vec3( 1.0, 0.0, 0.0));
	 a.push_back(Vec3(-1.0, 0.0, 0.0));
	 a.push_back(Vec3( 0.0, 1.0, 0.0));
	 a.push_back(Vec3( 0.0,-1.0, 0.0));
	 a.push_back(Vec3( 0.0, 0.0, 1.0));
	 a.push_back(Vec3( 0.0, 0.0,-1.0));
	 
	 b.push_back(Vec3(   p,   p, 0.0));
	 b.push_back(Vec3(   p,  -p, 0.0));
	 b.push_back(Vec3(  -p,   p, 0.0));
	 b.push_back(Vec3(  -p,  -p, 0.0));
	 b.push_back(Vec3(   p, 0.0,   p));
	 b.push_back(Vec3(   p, 0.0,  -p));
	 b.push_back(Vec3(  -p, 0.0,   p));
	 b.push_back(Vec3(  -p, 0.0,  -p));
	 b.push_back(Vec3( 0.0,   p,   p));
	 b.push_back(Vec3( 0.0,   p,  -p));
	 b.push_back(Vec3( 0.0,  -p,   p));
	 b.push_back(Vec3( 0.0,  -p,  -p));
	 
	 c.push_back(Vec3(   q,   q,   q));
	 c.push_back(Vec3(   q,   q,  -q));
	 c.push_back(Vec3(   q,  -q,   q));
	 c.push_back(Vec3(   q,  -q,  -q));
	 c.push_back(Vec3(  -q,   q,   q));
	 c.push_back(Vec3(  -q,   q,  -q));
	 c.push_back(Vec3(  -q,  -q,   q));
	 c.push_back(Vec3(  -q,  -q,  -q));
	 
	 d.push_back(Vec3(   r,   r,   s));
	 d.push_back(Vec3(   r,   r,  -s));
	 d.push_back(Vec3(   r,  -r,   s));
	 d.push_back(Vec3(   r,  -r,  -s));
	 d.push_back(Vec3(  -r,   r,   s));
	 d.push_back(Vec3(  -r,   r,  -s));
	 d.push_back(Vec3(  -r,  -r,   s));
	 d.push_back(Vec3(  -r,  -r,  -s));
	 
	 d.push_back(Vec3(   r,   s,   r));
	 d.push_back(Vec3(   r,   s,  -r));
	 d.push_back(Vec3(   r,  -s,   r));
	 d.push_back(Vec3(   r,  -s,  -r));
	 d.push_back(Vec3(  -r,   s,   r));
	 d.push_back(Vec3(  -r,   s,  -r));
	 d.push_back(Vec3(  -r,  -s,   r));
	 d.push_back(Vec3(  -r,  -s,  -r));
	 
	 d.push_back(Vec3(   s,   r,   r));
	 d.push_back(Vec3(   s,   r,  -r));
	 d.push_back(Vec3(   s,  -r,   r));
	 d.push_back(Vec3(   s,  -r,  -r));
	 d.push_back(Vec3(  -s,   r,   r));
	 d.push_back(Vec3(  -s,   r,  -r));
	 d.push_back(Vec3(  -s,  -r,   r));
	 d.push_back(Vec3(  -s,  -r,  -r));
	 
	 e.push_back(Vec3(   u,   v, 0.0));
	 e.push_back(Vec3(   u,  -v, 0.0));
	 e.push_back(Vec3(  -u,   v, 0.0));
	 e.push_back(Vec3(  -u,  -v, 0.0));
	 
	 e.push_back(Vec3(   u,  0.0,  v));
	 e.push_back(Vec3(   u,  0.0, -v));
	 e.push_back(Vec3(  -u,  0.0,  v));
	 e.push_back(Vec3(  -u,  0.0, -v));

	 e.push_back(Vec3(   0.0,  u,  v));
	 e.push_back(Vec3(   0.0,  u, -v));
	 e.push_back(Vec3(   0.0, -u,  v));
	 e.push_back(Vec3(   0.0, -u, -v));
	 
	 e.push_back(Vec3(   v,   u, 0.0));
	 e.push_back(Vec3(   v,  -u, 0.0));
	 e.push_back(Vec3(  -v,   u, 0.0));
	 e.push_back(Vec3(  -v,  -u, 0.0));
	 
	 e.push_back(Vec3(   v,  0.0,  u));
	 e.push_back(Vec3(   v,  0.0, -u));
	 e.push_back(Vec3(  -v,  0.0,  u));
	 e.push_back(Vec3(  -v,  0.0, -u));
	 
	 e.push_back(Vec3( 0.0,   v,   u));
	 e.push_back(Vec3( 0.0,   v,  -u));
	 e.push_back(Vec3( 0.0,  -v,   u));
	 e.push_back(Vec3( 0.0,  -v,  -u));
	 
	 
      }
      else if (symmetry == ICOSA) {
	 double t, p;  // theta and phi
      // a points
	 t = 0.0;  p=0.0;
	 a.push_back(Vec3(std::cos(t),std::sin(t)*std::cos(p),std::sin(t)*std::sin(p)));
	 t = M_PI;  p=0.0;
	 a.push_back(Vec3 (std::cos(t),std::sin(t)*std::cos(p),std::sin(t)*std::sin(p)));
      // b points
	 for (int k=0; k<5; k++) {
	    t = std::atan(2.0);          p = (double)(2*k+0)*M_PI/5.0;
	    b.push_back(Vec3 (std::cos(t), std::sin(t)*std::cos(p), std::sin(t)*std::sin(p)));
	    t = M_PI-std::atan(2.0);     p = (double)(2*k+1)*M_PI/5.0;
	    b.push_back(Vec3 (std::cos(t), std::sin(t)*std::cos(p), std::sin(t)*std::sin(p)));
	 }
	 // c points
	 double t1 = std::acos ((2.0+std::sqrt(5.0)) / std::sqrt(15.0+6.0*std::sqrt(5.0)));
	 double t2 = std::acos (      1.0            / std::sqrt(15.0+6.0*std::sqrt(5.0)));
	 for (int k=0; k<5; k++) {
	    t = t1; p = (double)(2*k+1)*M_PI/5.0;
	    c.push_back(Vec3 (std::cos(t), std::sin(t)*std::cos(p), std::sin(t)*std::sin(p)));
	    t = t2; p = (double)(2*k+1)*M_PI/5.0;
	    c.push_back(Vec3 (std::cos(t), std::sin(t)*std::cos(p), std::sin(t)*std::sin(p)));
	    t = M_PI - t1; p = (double)(2*k+0)*M_PI/5.0;
	    c.push_back(Vec3 (std::cos(t), std::sin(t)*std::cos(p), std::sin(t)*std::sin(p)));
	    t = M_PI - t2; p = (double)(2*k+0)*M_PI/5.0;
	    c.push_back(Vec3 (std::cos(t), std::sin(t)*std::cos(p), std::sin(t)*std::sin(p)));
	 }
      }
      // Now, construct rule
      if (std::fabs(A) > 1.0e-10) 
	for (int i=0; i<a.size(); i++)
	  addknot(a[i], A);
      if (std::fabs(B) > 1.0e-10) 
	for (int i=0; i<b.size(); i++)
	  addknot(b[i], B);
      if (std::fabs(C) > 1.0e-10) 
	for (int i=0; i<c.size(); i++)
	  addknot(c[i], C);
      if (std::fabs(D) > 1.0e-10) 
	for (int i=0; i<d.size(); i++)
	  addknot(d[i], D);
   }
   
   
   // Rotate rule so that it is not exactly aligned with the z-axis
   double alpha = 0.3823423409283498423424;
   double beta  = -1.98392480923809480293848;
   double gamma = 0.698390482380984238434208;
   Mat3 rx, ry, rz, rotmat;
   rx(0,0) =         1.0;  rx(0,1) =         0.0;   rx(0,2) =        0.0;
   rx(1,0) =         0.0;  rx(1,1) =  cos(alpha);   rx(1,2) = sin(alpha);
   rx(2,0) =         0.0;  rx(2,1) = -sin(alpha);   rx(2,2) = cos(alpha);
   
   ry(0,0) =   cos(beta);  ry(0,1) =         0.0;   ry(0,2) = -sin(beta);
   ry(1,0) =         0.0;  ry(1,1) =         1.0;   ry(1,2) =        0.0;
   ry(2,0) =   sin(beta);  ry(2,1) =         0.0;   ry(2,2) =  cos(beta); 
   
   rz(0,0) =  cos(gamma);  rz(0,1) =  sin(gamma);   rz(0,2) =        0.0;
   rz(1,0) = -sin(gamma);  rz(1,1) =  cos(gamma);   rz(1,2) =        0.0;
   rz(2,0) =         0.0;  rz(2,1) =         0.0;   rz(2,2) =        1.0;
   
   rotmat = rx * ry * rz;
   for (int k=0; k<rhat.size(); k++)
     if (dot (rhat[k], rhat[k]) < 0.99999)
       cerr << "rhat[" << k << "] = " << rhat[k] << endl;
   for (int k=0; k < rhat.size(); k++)
     rhat[k] = rotmat * rhat[k];
   
   // Finally, check the rule for correctness
   Check();
   
   double wSum = 0.0;
   for (int k=0; k < nk; k++) {
      Vec3 r = rhat[k];
      double nrm = dot(r,r);
      rhat[k] /= sqrt(nrm);
      nrm = dot(rhat[k], rhat[k]);
      assert (std::fabs(nrm-1.0) < 1.0e-14);
    wSum += weight[k];
    //cout << pp_nonloc->sgridxyz_m[k] << " " << pp_nonloc->sgridweight_m[k] << endl;
   }
   // fprintf (stderr, "wSum = %1.16f\n", wSum);
   assert (std::fabs(wSum - 1.0) < 1.0e-14);
}

#include <time.h>
void
AtomicOrbital::TimeYlm()
  {
    clock_t start, end;
    Vec3 rhat(0.1, 0.3, 0.8);
    start = clock();
    for (int i=0; i<1000000; i++)
      CalcYlm(rhat);
    end = clock();
    double time = (double)(end-start)/(double)CLOCKS_PER_SEC;
    fprintf (stderr, "Ylm evals per second = %1.4f\n",
	     1.0e6/time);
  }

  // Fast implementation
  // See Geophys. J. Int. (1998) 135,pp.307-309
  void
  AtomicOrbital::CalcYlm (TinyVector<double,3> rhat)
  {
    // if (fabs(rhat[0]) < 1.0e-10 &&
    // 	fabs(rhat[1]) < 1.0e-10 &&
    // 	fabs(rhat[2]- 1.0) < 1.0e-10)
    //   rhat = Vec3 (1.0e-6, 1.0e-6, 0.999999999999);

    // if (fabs(rhat[0]) < 1.0e-10 &&
    // 	fabs(rhat[1]) < 1.0e-10 &&
    // 	fabs(rhat[2]+ 1.0) < 1.0e-10)
    //   rhat = Vec3 (1.0e-6, 1.0e-6, -0.999999999999);


    const double fourPiInv = 0.0795774715459477;
    
    double costheta, sintheta, cottheta;
    double cosphi, sinphi;
    costheta = rhat[2];
    sintheta = std::sqrt(1.0-costheta*costheta);
    if (sintheta > 1.0e-6) {
      cottheta = costheta/sintheta;
      cosphi=rhat[0]/sintheta;
      sinphi=rhat[1]/sintheta; 
    }
    else {
      sintheta = 1.0e-8;
      cottheta = costheta * 1.0e8;
      cosphi = 1.0;
      sinphi = 0.0;
    }
    
    complex<double> e2iphi(cosphi, sinphi);
    
    
    double lsign = 1.0;
    double dl = 0.0;
    for (int l=0; l<=lMax; l++) {
      double XlmVec[2*l+1], dXlmVec[2*l+1];
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
      complex<double> e2imphi (1.0, 0.0);
      for (int m=0; m<=l; m++) {
	YlmVec[l*(l+1)+m]  =  XlmVec[l+m]*e2imphi;
	YlmVec[l*(l+1)-m]  =  XlmVec[l-m]*conj(e2imphi);
	dYlmVec[l*(l+1)+m] = dXlmVec[l+m]*e2imphi;
	dYlmVec[l*(l+1)-m] = dXlmVec[l-m]*conj(e2imphi);
	e2imphi *= e2iphi;
      } 
      
      dl += 1.0;
      lsign *= -1.0;
    }
  }



static int orbNum = 0;

void
AtomicOrbital::Project (OrbFunctor &phi,
			OrbFunctor &spline,
			Array<complex<double>,2> &radialFuncs)
{
  // HACK HACK HACK
  double dr = OuterRadius/(SplinePoints-1);
  double r = 0.0;
  int numquad = QuadRule.rhat.size();
  int numlm   = YlmVec.size();
  

  radialFuncs.resize(SplinePoints, numlm);

  for (int iq=0; iq<numquad; iq++) {
    CalcYlm (QuadRule.rhat[iq]);
    for (int lm=0; lm<numlm; lm++) 
      YlmArray(lm,iq) = YlmVec[lm];
  }
    
  vector<complex<double> > phival(numquad*SplinePoints);
  vector<Vec3> rvecs(numquad*SplinePoints);
// #pragma omp parallel 
// {
// #pragma omp for schedule(dynamic,1)


  int index=0;
  for (int ir=0; ir<SplinePoints; ir++) {
    double r = dr*ir;
    for (int iq=0; iq<numquad; iq++)
      rvecs[index++] = Pos + r * QuadRule.rhat[iq];
  }
  phi(rvecs, phival);

  double power[lMax+1];
  for (int l=0; l<=lMax; l++)
    power[l] = 0.0;
  for (int ir=0; ir<SplinePoints; ir++) {
    // double r = dr*ir;
    // for (int iq=0; iq < numquad; iq++) {
    //   Vec3 rvec = Pos + r * QuadRule.rhat[iq];
    //   phival[iq] = phi(rvec);
    // }
    //    double r2l = 1.0;
    int off = ir*numquad;
    double realPower = 0.0;
    double spherePower = 0.0;
    complex<double> mean = 0.0;
    for (int iq=0; iq < numquad; iq++)
      mean += QuadRule.weight[iq] * phival[off+iq];
    //    mean /= (double)numquad;
    for (int iq=0; iq < numquad; iq++)
      realPower += 4.0*M_PI*QuadRule.weight[iq] * norm(phival[off+iq]-mean);
    for (int l=0; l<=lMax; l++) {
      for (int m=-l; m<=l; m++) {
	int lm = l*(l+1)+m;
	complex<double> sum;
	for (int iq=0; iq < numquad; iq++) 
	  // sum += QuadRule.weight[iq] * phival[iq] * conj(YlmArray(lm,iq));
	  sum += QuadRule.weight[iq] * phival[off+iq] * conj(YlmArray(lm,iq));
	radialFuncs(ir,lm) = 4.0*M_PI *sum;
	if (lm != 0)
	  spherePower += norm(radialFuncs(ir,lm));
	power[l] += norm(radialFuncs(ir,lm));
      }
    }
    // double renorm = sqrt(realPower/spherePower);
    // fprintf (stderr, "Renormalizing by %1.18f\n", renorm);
    // for (int lm=1; lm< numlm; lm++) {
    //   radialFuncs(ir,lm) *= renorm;
    // }

    // renorm = sqrt(realPower/newspherePower);
    // fprintf (stderr, "New ratio =  by %1.18f\n", renorm);

    // fprintf (stderr, "realPower   = %18.16e\n", realPower);
    // fprintf (stderr, "spherePower = %18.16e\n", spherePower);
  }
  for (int l=0; l<=lMax; l++) 
    fprintf (stderr, "   Power[%d] = %10.6e\n", l, power[l]);
  fprintf (stderr, "\n");


  // for (int lm=0; lm<numlm; lm++)
  //   radialFuncs(0,lm) = ( 4.0*radialFuncs(1,lm)
  // 			 -6.0*radialFuncs(2,lm)
  // 			 +4.0*radialFuncs(3,lm)
  // 			 -1.0*radialFuncs(4,lm));

  
  vector<complex<double> > ulm(SplinePoints);
  for (int lm=0; lm<numlm; lm++) {
    for (int ir=0; ir < SplinePoints; ir++) 
      ulm[ir] = radialFuncs(ir,lm);
    set_multi_UBspline_1d_z (Spline, lm, &(ulm[0]));
  }
  
 
  // // Do Chebyshev fits
  // int Ncheb = 25, Nquad = 25;
  // ChebyshevCoefs.resize(Ncheb, numlm);
  // Chebyshev.resize(Ncheb);
  // Tn.resize(Ncheb); dTn.resize(Ncheb); d2Tn.resize(Ncheb);
  // ChebyshevCoefs = complex<double>();
  // ShiftedChebyshevQuad cquad(Nquad);
  // Array<complex<double>,2> chebFuncs(Nquad, numlm);
  // ShiftedChebyshevPoly<complex<double> > cheb(Ncheb);
  // for (int lm=0; lm<numlm; lm++)
  //   chebFuncs.resize(Ncheb);
  // for (int irad=0; irad<Nquad; irad++) {
  //   double x = cquad.x(irad);
  //   double r = OuterRadius*std::fabs(x);
  //   double w = cquad.w(irad);
  //   cheb(x, Tn);
  //   for (int isph=0; isph < numquad; isph++) {
  //     Vec3 rvec = Pos + r * QuadRule.rhat[isph];
  //     phival[isph] = phi(rvec);
  //   }
  //   for (int l=0; l<=lMax; l++) {
  //     for (int m=-l; m<=l; m++) {
  // 	int lm = l*(l+1)+m;
  // 	for (int n=0; n<Ncheb; n++) {
  // 	  double prefactor = (n==0) ? 1.0/M_PI : 2.0/M_PI;
  // 	  for (int isph=0; isph < numquad; isph++)
  // 	    ChebyshevCoefs(n,lm) += w * Tn[n] * prefactor * 
  // 	      4.0*M_PI * QuadRule.weight[isph] * phival[isph] * conj(YlmArray(lm, isph));
  // 	}
  //     }
  //   }
  // }
  
	
  // // Now, condition coefficients at the origin
  // cheb(0.0, Tn, dTn, d2Tn);

  // // For l >= 3
  // for (int lm=9; lm<numlm; lm++) {
  //   complex<double> d2u;
  //   for (int n=0; n<Ncheb; n++)
  //     d2u += d2Tn[n] * ChebyshevCoefs(n,lm);
  //   ChebyshevCoefs(2,lm) -= d2u/d2Tn[2];
  // }

  // // For l  >= 2
  // for (int lm=4; lm<numlm; lm++) {
  //   complex<double> du;
  //   for (int n=0; n<Ncheb; n++)
  //     du += dTn[n] * ChebyshevCoefs(n,lm);
  //   ChebyshevCoefs(1,lm) -= du/dTn[1];
  // }

  //   // For l  >= 1
  // for (int lm=1; lm<numlm; lm++) {
  //   complex<double> u;
  //   for (int n=0; n<Ncheb; n++)
  //     u += Tn[n] * ChebyshevCoefs(n,lm);
  //   ChebyshevCoefs(0,lm) -= u/Tn[0];
  //   } 
						  
  // Do polynomial fits
  Array<complex<double>,2> polyFuncs(SplinePoints, numlm);
  double pdr = 1.5*PolyRadius / (SplinePoints-1);
  //#pragma omp for schedule(dynamic,1)

  index=0;
  for (int ir=0; ir<SplinePoints; ir++) {
    double r = pdr*ir;
    for (int iq=0; iq<numquad; iq++)
      rvecs[index++] = Pos + r * QuadRule.rhat[iq];
  }
  phi(rvecs, phival);

  for (int ir=0; ir<SplinePoints; ir++) {
    // double r = ir*pdr;
    // for (int iq=0; iq < numquad; iq++) {
    //   Vec3 rvec = Pos + r * QuadRule.rhat[iq];
    //   phival[iq] = phi(rvec);
    // }
    int off = ir*numquad;
    for (int l=0; l<=lMax; l++) {
      for (int m=-l; m<=l; m++) {
	int lm = l*(l+1)+m;
	complex<double> sum;
	for (int iq=0; iq < numquad; iq++)
	  sum += QuadRule.weight[iq] * phival[off+iq] * conj(YlmArray(lm,iq));
	polyFuncs(ir,lm) = 4.0*M_PI*sum;
      }
    }
  }


  // Set up bases
  // Array<double,2> sBasis(SplinePoints, PolyOrder);
  // Array<double,2> non_sBasis(SplinePoints, PolyOrder-1);
  // for (int ir=0; ir<SplinePoints; ir++) {
  //   double r = pdr*ir;
  //   sBasis(ir,0) = 1.0;
  //   double r2n = r*r;
  //   for (int n=2; n <= PolyOrder; n++) {
  //     sBasis(ir, n-1) = r2n;
  //     non_sBasis(ir, n-2) = r2n;
  //     r2n *= r;
  //   }
  // }

  for (int lm=0; lm < numlm; lm++) {
    PolyCoefs(0,lm) = complex<double>();
    PolyCoefs(1,lm) = complex<double>();
  }

  // Do fits for l=m=0
  Array<double,1> y_re(SplinePoints), y_im(SplinePoints), sigma(SplinePoints);
  Array<double,1> coefs_re(PolyOrder), errors_re(PolyOrder),
    coefs_im(PolyOrder), errors_im(PolyOrder);

  for (int l=0; l<=lMax; l++) {
    coefs_re.resize(PolyOrder-l+1);    coefs_im.resize(PolyOrder-l+1);
    Array<double,2> Basis(SplinePoints, PolyOrder-l+1);
    for (int ir=0; ir<SplinePoints; ir++) {
      double r=pdr*ir;
      double r2lpn=1.0;
      for (int i=0; i<l; i++)	r2lpn *= r;
      for (int n=0; n<PolyOrder-l+1; n++) {
	Basis(ir,n) = r2lpn;
	r2lpn *= r;
      }
    }
    for (int m=-l; m<=l; m++) {
      int lm = l*(l+1)+m;
      for (int ir=0; ir<SplinePoints; ir++) {
	y_re(ir) = polyFuncs(ir,lm).real();
	y_im(ir) = polyFuncs(ir,lm).imag();
	sigma(ir) = 1.0e-40;
      }
      LinFitLU (y_re, sigma, Basis, coefs_re, errors_re);
      LinFitLU (y_im, sigma, Basis, coefs_im, errors_im);

      for (int n=0; n<l; n++)
	PolyCoefs(n,lm) = complex<double>();
      for (int n=l; n<PolyOrder; n++)
	PolyCoefs(n,lm) = complex<double>(coefs_re(n-l), coefs_im(n-l));
    }
  }

  /*
  // Fill y arrays
  for (int ir=0; ir<SplinePoints; ir++) {
    y_re(ir) = polyFuncs(ir,0).real();
    y_im(ir) = polyFuncs(ir,0).imag();
    sigma(ir) = 1.0e-16;
  }
      
  LinFitLU (y_re, sigma, sBasis, coefs_re, errors_re);
  LinFitLU (y_re, sigma, sBasis, coefs_im, errors_im);
  PolyCoefs(0,0) = complex<double> (coefs_re(0), coefs_im(0));
  for (int n=2; n<=PolyOrder; n++)
    PolyCoefs(n,0) = complex<double> (coefs_re(n-1), coefs_im(n-1));
  
  // Do fits for l!=0
  coefs_re.resize(PolyOrder-1);   coefs_im.resize(PolyOrder-1);
  errors_re.resize(PolyOrder-1);  errors_im.resize(PolyOrder-1);
  for (int lm=1; lm<numlm; lm++) {
    // Fill y arrays
    for (int ir=0; ir<SplinePoints; ir++) {
      y_re(ir) = polyFuncs(ir,lm).real();
      y_im(ir) = polyFuncs(ir,lm).imag();
      sigma(ir) = 1.0e-16;
    }
    // Do fits
    LinFitLU (y_re, sigma, non_sBasis, coefs_re, errors_re);
    LinFitLU (y_re, sigma, non_sBasis, coefs_im, errors_im);
    PolyCoefs(0,lm) = complex<double>();
    PolyCoefs(1,lm) = complex<double>();
    for (int n=2; n<=PolyOrder; n++)
      PolyCoefs(n,lm) = complex<double> (coefs_re(n-2), coefs_im(n-2));
  }
  */

  Vec3 rvec;
  rvec[1] = Pos[1]+0.000001;
  rvec[2] = Pos[2]+0.000001;

  char fname[500];
  snprintf (fname, 500, "test_orb_%d.dat", orbNum++);
  FILE *fout = fopen (fname, "w");
  fprintf (fout, "#  Pos = %10.5f %10.5f %10.5f\n",
  	   Pos[0], Pos[1], Pos[2]);

  
  rvecs.clear();
  // for (double x=Pos[0]-Radius-1.0; x <= Pos[0]+Radius+1.0; x+=0.02) {
  //   rvec[0] = x;

  Vec3 u(0.1, 0.2, 0.15);
  u = 1.0/sqrt(dot(u,u)) * u;
  for (double x=-1.1; x <= 1.1; x+=0.01) {
    rvec = Pos + x*Radius*u;
    rvecs.push_back(rvec);
  }
  int numr = rvecs.size();
  vector<complex<double> > ex_val(numr), ex_lapl(numr), sp_val(numr), sp_lapl(numr),
    at_val(numr), at_lapl(numr);
  phi(rvecs, ex_val, ex_lapl);
  spline(rvecs, sp_val, sp_lapl);
  for (int ir=0; ir<numr; ir++) {
    TinyVector<complex<double>,3> at_grad;
    eval(rvecs[ir], at_val[ir], at_grad, at_lapl[ir]);

    fprintf (fout, "%12.8f  "
	     "%16.12e %16.12e %16.12e %16.12e  " 
	     "%16.12e %16.12e %16.12e %16.12e   "
	     "%16.12e %16.12e %16.12e %16.12e\n",
	     rvecs[ir][0], 	     
	     ex_val[ir].real(), ex_val[ir].imag(), ex_lapl[ir].real(), ex_lapl[ir].imag(),
	     at_val[ir].real(), at_val[ir].imag(), at_lapl[ir].real(), at_lapl[ir].imag(),
	     sp_val[ir].real(), sp_val[ir].imag(), sp_lapl[ir].real(), sp_lapl[ir].imag());
  }
  fclose(fout);

  // for (double x=Pos[0]-Radius-1.0; x <= Pos[0]+Radius+1.0; x+=0.02) {
  //   rvec[0] = x;
  //   complex<double> ex_val, ex_lapl, at_val, at_lapl, poly_val, cheb_val, cheb_lapl;
  //   complex<double> spline_val, spline_lapl;
  //   TinyVector<complex<double>,3> ex_grad, at_grad, cheb_grad, spline_grad;
  //   phi(rvec, ex_val, ex_grad, ex_lapl);
  //   //    ex_lapl = lapl(rvec);
  //   spline (rvec, spline_val, spline_grad, spline_lapl);
  //   eval(rvec, at_val, at_grad, at_lapl);
  //   poly_val = eval_poly (rvec);
  //   //eval_cheb (rvec, cheb_val, cheb_grad, cheb_lapl);
  //   fprintf (fout, "%12.8f  %16.12e  %16.12e %16.12e %16.12e  %16.12e  %16.12e  %16.12e %16.12e %16.12e %16.12e %16.12e %16.12e\n", 
  // 	     rvec[0], ex_val.imag(), at_val.imag(),
  // 	     ex_lapl.imag(), at_lapl.imag(), poly_val.imag(),
  // 	     cheb_val.imag(), cheb_lapl.imag(),
  // 	     spline_val.imag(), spline_lapl.imag(),
  // 	     ex_grad[0].imag(), at_grad[0].imag(), spline_grad[0].imag());
  // }
  // fclose (fout);

}







void
AtomicOrbital::ProjectAnalytic (PWOrb_z &phi,
				OrbFunctor &spline,
				Array<complex<double>,2> &radialFuncs)
{
  //  phi.k = Vec3(0.2123,0.18343,0.082348);
  double dr = OuterRadius/(SplinePoints-1);
  int numquad = QuadRule.rhat.size();
  int numlm   = YlmVec.size();
  radialFuncs.resize(SplinePoints, numlm);

  Array<complex<double>,2> Ylm_G(phi.GVecs.size(), numlm);
  for (int iG=0; iG<phi.GVecs.size(); iG++) {
    Vec3 G = phi.GVecs(iG) + phi.k;
    double Gmag = sqrt(dot (G,G));
    if (Gmag > 0.0)
      G *= 1.0/Gmag;
    else
      G = Vec3(0.0, 0.0, 1.0);
    CalcYlm (G);
    for (int lm=0; lm<numlm; lm++) 
      Ylm_G(iG, lm) = YlmVec[lm];
  }

  // Convert to spherical harmonics with analytic sum
  complex<double> minus_i(0.0,  -1.0);
  vector<double> jl(lMax+1);
  vector<complex<double> > ulm_sum((lMax+1)*(lMax+1));

  for (int ir=0; ir<SplinePoints; ir++) {
    double r = dr*ir;
    for (int lm=0; lm<numlm; lm++)
      ulm_sum[lm] = 0.0;
    Vec3 rdummy(0.3,0.2,0.6);
    CalcYlm(rdummy/sqrt(dot(rdummy,rdummy)));
    for (int iG=0; iG<phi.GVecs.size(); iG++) {
      Vec3 G = phi.GVecs(iG) + phi.k;
      double re, im;
      sincos(-dot(rdummy,G), &im, &re);
      complex<double> z(re,im);
      
      double Gmag = sqrt(dot (G,G));
      double phase = -dot(G, Pos);
      sincos(phase, &im, &re);
      complex<double> phase_shift(re,im);
      double Gr = Gmag * r;
      gsl_sf_bessel_jl_steed_array(lMax, Gr, &(jl[0]));
      complex<double> i2l(1.0,0.0);
      complex<double> zsum(0.0, 0.0);
      for (int l=0; l<=lMax; l++) {
	for (int m=-l; m<=l; m++) {
	  int lm = l*(l+1)+m;
	  radialFuncs(ir,lm) += 4.0*M_PI*phi.Prefactor*
	    phi.Coefs(iG) * phase_shift * i2l * jl[l] * conj(Ylm_G(iG, lm));
	  // zsum += 4.0*M_PI*i2l*
	  //   gsl_sf_bessel_jl(l, Gmag*sqrt(dot(rdummy,rdummy))) *
	  //   conj(Ylm_G(iG,lm)) * YlmVec[lm];
	}
	i2l *= minus_i;
      }
      // fprintf (stderr, "z    = %10.6e + i*%10.6e\n", z.real(),    z.imag());
      // fprintf (stderr, "zsum = %10.6e + i*%10.6e\n", zsum.real(), zsum.imag());
    }
    for (int lm=0; lm<numlm; lm++) {
      //radialFuncs(ir,lm) = 4.0*M_PI * phi.Prefactor * ulm_sum[lm];
      // fprintf (stderr, "radialFuncs(%d,%d) = %12.8f + i*%12.8f\n",
      // 	       ir, lm,radialFuncs(ir,lm).real(),  radialFuncs(ir,lm).imag());
    }
    r += dr;
  }
  
  vector<complex<double> > ulm(SplinePoints);
  for (int lm=0; lm<numlm; lm++) {
    for (int ir=0; ir < SplinePoints; ir++) 
      ulm[ir] = radialFuncs(ir,lm);
    set_multi_UBspline_1d_z (Spline, lm, &(ulm[0]));
  }
   
  vector<Vec3> rvecs(numquad*SplinePoints);
  
  // Do polynomial fits
  for (int iq=0; iq<numquad; iq++) {
    CalcYlm (QuadRule.rhat[iq]);
    for (int lm=0; lm<numlm; lm++) 
      YlmArray(lm,iq) = YlmVec[lm];
  }

  vector<complex<double> > phival(numquad*SplinePoints);
  Array<complex<double>,2> polyFuncs(SplinePoints, numlm);
  double pdr = 1.5*PolyRadius / (SplinePoints-1);

  int index=0;
  for (int ir=0; ir<SplinePoints; ir++) {
    double r = pdr*ir;
    for (int iq=0; iq<numquad; iq++)
      rvecs[index++] = Pos + r * QuadRule.rhat[iq];
  }
  phi(rvecs, phival);

  for (int ir=0; ir<SplinePoints; ir++) {
    int off = ir*numquad;
    for (int l=0; l<=lMax; l++) {
      for (int m=-l; m<=l; m++) {
  	int lm = l*(l+1)+m;
  	complex<double> sum;
  	for (int iq=0; iq < numquad; iq++)
  	  sum += QuadRule.weight[iq] * phival[off+iq] * conj(YlmArray(lm,iq));
  	polyFuncs(ir,lm) = 4.0*M_PI*sum;
      }
    }
  }


  for (int lm=0; lm < numlm; lm++) {
    PolyCoefs(0,lm) = complex<double>();
    PolyCoefs(1,lm) = complex<double>();
  }

  // Do fits for l=m=0
  Array<double,1> y_re(SplinePoints), y_im(SplinePoints), sigma(SplinePoints);
  Array<double,1> coefs_re(PolyOrder), errors_re(PolyOrder),
    coefs_im(PolyOrder), errors_im(PolyOrder);

  for (int l=0; l<=lMax; l++) {
    coefs_re.resize(PolyOrder-l+1);    coefs_im.resize(PolyOrder-l+1);
    Array<double,2> Basis(SplinePoints, PolyOrder-l+1);
    for (int ir=0; ir<SplinePoints; ir++) {
      double r=pdr*ir;
      double r2lpn=1.0;
      for (int i=0; i<l; i++)	r2lpn *= r;
      for (int n=0; n<PolyOrder-l+1; n++) {
  	Basis(ir,n) = r2lpn;
  	r2lpn *= r;
      }
    }
    for (int m=-l; m<=l; m++) {
      int lm = l*(l+1)+m;
      for (int ir=0; ir<SplinePoints; ir++) {
  	y_re(ir) = polyFuncs(ir,lm).real();
  	y_im(ir) = polyFuncs(ir,lm).imag();
  	sigma(ir) = 1.0e-40;
      }
      LinFitLU (y_re, sigma, Basis, coefs_re, errors_re);
      LinFitLU (y_im, sigma, Basis, coefs_im, errors_im);

      for (int n=0; n<l; n++)
  	PolyCoefs(n,lm) = complex<double>();
      for (int n=l; n<PolyOrder; n++)
  	PolyCoefs(n,lm) = complex<double>(coefs_re(n-l), coefs_im(n-l));
    }
  }

  Vec3 rvec;
  rvec[1] = Pos[1]+0.000001;
  rvec[2] = Pos[2]+0.000001;

  char fname[500];
  snprintf (fname, 500, "test_orb_%d.dat", orbNum++);
  FILE *fout = fopen (fname, "w");
  fprintf (fout, "#  Pos = %10.5f %10.5f %10.5f\n",
  	   Pos[0], Pos[1], Pos[2]);

  
  rvecs.clear();
  // for (double x=Pos[0]-Radius-1.0; x <= Pos[0]+Radius+1.0; x+=0.02) {
  //   rvec[0] = x;

  //  Vec3 u(0.1, 0.2, 0.15);
  Vec3 u(0.4, 0.2, 0.15);
  //Vec3 u(0.0, 1.0, 0.0);
  u = 1.0/sqrt(dot(u,u)) * u;
  vector<double> xvec;
  for (double x=-1.1; x <= 1.1; x+=0.01) {
    rvec = Pos + x*Radius*u;
    rvecs.push_back(rvec);
    xvec.push_back(x);
  }
  int numr = rvecs.size();
  vector<complex<double> > ex_val(numr), ex_lapl(numr), sp_val(numr), sp_lapl(numr),
    at_val(numr), at_lapl(numr);
  phi(rvecs, ex_val, ex_lapl);
  spline(rvecs, sp_val, sp_lapl);
  for (int ir=0; ir<numr; ir++) {
    TinyVector<complex<double>,3> at_grad;
    eval(rvecs[ir], at_val[ir], at_grad, at_lapl[ir]);

    fprintf (fout, "%12.8f  "
	     "%16.12e %16.12e %16.12e %16.12e  " 
	     "%16.12e %16.12e %16.12e %16.12e   "
	     "%16.12e %16.12e %16.12e %16.12e\n",
	     xvec[ir],//norm(rvecs[ir])*sign[ir], 	     
	     ex_val[ir].real(), ex_val[ir].imag(), ex_lapl[ir].real(), ex_lapl[ir].imag(),
	     at_val[ir].real(), at_val[ir].imag(), at_lapl[ir].real(), at_lapl[ir].imag(),
	     sp_val[ir].real(), sp_val[ir].imag(), sp_lapl[ir].real(), sp_lapl[ir].imag());
  }
  fclose(fout);

  // for (double x=Pos[0]-Radius-1.0; x <= Pos[0]+Radius+1.0; x+=0.02) {
  //   rvec[0] = x;
  //   complex<double> ex_val, ex_lapl, at_val, at_lapl, poly_val, cheb_val, cheb_lapl;
  //   complex<double> spline_val, spline_lapl;
  //   TinyVector<complex<double>,3> ex_grad, at_grad, cheb_grad, spline_grad;
  //   phi(rvec, ex_val, ex_grad, ex_lapl);
  //   //    ex_lapl = lapl(rvec);
  //   spline (rvec, spline_val, spline_grad, spline_lapl);
  //   eval(rvec, at_val, at_grad, at_lapl);
  //   poly_val = eval_poly (rvec);
  //   //eval_cheb (rvec, cheb_val, cheb_grad, cheb_lapl);
  //   fprintf (fout, "%12.8f  %16.12e  %16.12e %16.12e %16.12e  %16.12e  %16.12e  %16.12e %16.12e %16.12e %16.12e %16.12e %16.12e\n", 
  // 	     rvec[0], ex_val.imag(), at_val.imag(),
  // 	     ex_lapl.imag(), at_lapl.imag(), poly_val.imag(),
  // 	     cheb_val.imag(), cheb_lapl.imag(),
  // 	     spline_val.imag(), spline_lapl.imag(),
  // 	     ex_grad[0].imag(), at_grad[0].imag(), spline_grad[0].imag());
  // }
  // fclose (fout);

}



complex<double> 
AtomicOrbital::eval (Vec3 r) 
{
  Vec3 diff = r - Pos;
  double rmag = sqrt (dot(diff,diff));
  if (rmag > Radius)
    return complex<double>();

  Vec3 rhat = (1.0/rmag) * diff;

  CalcYlm (rhat);

  complex<double> ulm[YlmVec.size()];

  eval_multi_UBspline_1d_z (Spline, rmag, &(ulm[0]));
  complex<double> z;
  double r2l = 1.0;
  int lm=0;
  for (int l=0; l<=lMax; l++) {
    for (int m=-l; m<=l; m++,lm++) 
      z += /*r2l * */ulm[lm] * YlmVec[lm];
    r2l *= rmag;
  }

  return z;
}

complex<double> 
AtomicOrbital::eval_poly (Vec3 r) 
{
  Vec3 diff = r - Pos;
  double rmag = sqrt (dot(diff,diff));
  if (rmag > Radius)
    return complex<double>();

  Vec3 rhat = (1.0/rmag) * diff;

  CalcYlm (rhat);

  complex<double> ulm[YlmVec.size()];
  for (int lm=0; lm<YlmVec.size(); lm++)
    ulm[lm] = complex<double>();
  double r2n=1.0;
  for (int n=0; n<=PolyOrder; n++) {
    for (int lm=0; lm<YlmVec.size(); lm++)
      ulm[lm] += r2n*PolyCoefs(n,lm);
    r2n *= rmag;
  }

  complex<double> z;
  for (int i=0; i<YlmVec.size(); i++) 
    z += ulm[i] * YlmVec[i];

  return z;
}



complex<double> 
AtomicOrbital::eval_cheb (Vec3 r) 
{
  Vec3 diff = r - Pos;
  double rmag = sqrt (dot(diff,diff));
  if (rmag > Radius)
    return complex<double>();

  double x = rmag / OuterRadius;
  Vec3 rhat = (1.0/rmag) * diff;

  CalcYlm (rhat);

  complex<double> ulm[YlmVec.size()];
  Chebyshev(x, Tn);

  for (int lm=0; lm<YlmVec.size(); lm++)
    ulm[lm] = complex<double>();
  double r2n=1.0;
  for (int n=0; n<=Tn.size(); n++) 
    for (int lm=0; lm<YlmVec.size(); lm++)
	ulm[lm] += Tn[n]*ChebyshevCoefs(n,lm);

  complex<double> z;
  for (int i=0; i<YlmVec.size(); i++) 
    z += ulm[i] * YlmVec[i];

  return z;
}

void
AtomicOrbital::eval_cheb (Vec3 r,
			  complex<double> &val,
			  TinyVector<complex<double>,3> &grad,
			  complex<double> &lapl)
{
  Vec3 diff = r - Pos;
  double rmag = sqrt (dot(diff,diff));
  val = lapl = complex<double>();
  grad = TinyVector<complex<double>,3>();
  if (rmag > Radius)
    return;

  if (rmag < 1.0e-4)
    rmag = 1.0e-4;

  double x = rmag / OuterRadius;
  Vec3 rhat, thetahat, phihat;
  rhat = (1.0/rmag) * diff;
  double costheta = rhat[2];
  double sintheta = std::sqrt(1.0-costheta*costheta);
  double cosphi = rhat[0]/sintheta;
  double sinphi = rhat[1]/sintheta;
  thetahat = TinyVector<double,3>(costheta*cosphi,
				  costheta*sinphi,
				  -sintheta);
  phihat   = TinyVector<double,3>(-sinphi,
				  cosphi,
				  0.0 );

  CalcYlm (rhat);

  complex<double> ulm[YlmVec.size()], 
    dulm[YlmVec.size()], d2ulm[YlmVec.size()];
  Chebyshev(x, Tn, dTn, d2Tn);
  for (int lm=0; lm<YlmVec.size(); lm++) 
    ulm[lm] = dulm[lm] = d2ulm[lm] = complex<double>();
  double gnorm = 1.0/OuterRadius;
  double lnorm = gnorm * gnorm;
  for (int n=0; n<=ChebyshevCoefs.extent(0); n++) {
    for (int lm=0; lm<YlmVec.size(); lm++) {
      ulm[lm]   +=   Tn[n] * ChebyshevCoefs(n,lm);
      dulm[lm]  +=  dTn[n] * ChebyshevCoefs(n,lm) * gnorm;
      d2ulm[lm] += d2Tn[n] * ChebyshevCoefs(n,lm) * lnorm;
    }
  }

  val = lapl = complex<double>();
  grad = TinyVector<complex<double>,3>();
  int lm =0;
  for (int l=0; l<=lMax; l++)
    for (int m=-l; m<=l; m++) {
      complex<double> im(0.0,(double)m);
      val += ulm[lm] * YlmVec[lm];

      grad += 
	(dulm[lm]                *     YlmVec[lm] * rhat     +
	 ulm[lm]/rmag            *    dYlmVec[lm] * thetahat +
	 ulm[lm]/(rmag*sintheta) * im *YlmVec[lm] * phihat);

      lapl += YlmVec[lm] * 
	(-(double)(l*(l+1))/(rmag*rmag) * ulm[lm]
	 + d2ulm[lm] + 2.0/rmag *dulm[lm]);
      lm++;
  }
}




void 
AtomicOrbital::eval (Vec3 r, LatticeClass &lattice,
		     complex<double> &val,
		     TinyVector<complex<double>,3> &grad,
		     complex<double> &lapl)
{
  Vec3 diff = r - Pos;
  Vec3 u = lattice.r2u(diff);
  for (int i=0; i<3; i++) u[i] -= round(u[i]);
  diff = lattice.u2r(u);
  r = Pos + diff;
  double rmag = sqrt (dot(diff,diff));
  val = lapl = complex<double>();
  grad = TinyVector<complex<double>,3>();
  if (rmag > Radius) 
    return;

  eval (r, val, grad, lapl);
}


void
AtomicOrbital::eval_poly (Vec3 r,
			  complex<double> &val,
			  TinyVector<complex<double>,3> &grad,
			  complex<double> &lapl)
{
  Vec3 diff = r - Pos;
  double rmag = sqrt (dot(diff,diff));
  val = lapl = complex<double>();
  grad = TinyVector<complex<double>,3>();
  if (rmag > Radius)
    return;
  if (rmag < 1.0e-4)
    rmag = 1.0e-4;

  Vec3 rhat, thetahat, phihat;
  rhat = (1.0/rmag) * diff;
  double costheta = rhat[2];
  double sintheta = std::sqrt(1.0-costheta*costheta);
  double cosphi = rhat[0]/sintheta;
  double sinphi = rhat[1]/sintheta;
  thetahat = TinyVector<double,3>(costheta*cosphi,
				  costheta*sinphi,
				  -sintheta);
  phihat   = TinyVector<double,3>(-sinphi,
				  cosphi,
				  0.0 );

  CalcYlm (rhat);

  complex<double> ulm[YlmVec.size()], 
    dulm[YlmVec.size()], d2ulm[YlmVec.size()];
  for (int lm=0; lm<YlmVec.size(); lm++) 
    ulm[lm] = dulm[lm] = d2ulm[lm] = complex<double>();
  double r2n=1.0; 
  double r2nm1;
  double r2nm2;
  for (int n=0; n<=PolyOrder; n++) {
    for (int lm=0; lm<YlmVec.size(); lm++) {
      ulm[lm]   +=                     r2n   * PolyCoefs(n,lm);
      dulm[lm]  += (double)n         * r2nm1 * PolyCoefs(n,lm);
      d2ulm[lm] += (double)(n*(n-1)) * r2nm2 * PolyCoefs(n,lm);
    }
    r2nm2 = r2nm1;
    r2nm1 = r2n;
    r2n *= rmag;
  }

  val = lapl = complex<double>();
  grad = TinyVector<complex<double>,3>();
  int lm =0;
  for (int l=0; l<=lMax; l++)
    for (int m=-l; m<=l; m++) {
      complex<double> im(0.0,(double)m);
      val += ulm[lm] * YlmVec[lm];

      grad += 
	(dulm[lm]                *     YlmVec[lm] * rhat     +
	 ulm[lm]/rmag            *    dYlmVec[lm] * thetahat +
	 ulm[lm]/(rmag*sintheta) * im *YlmVec[lm] * phihat);

      lapl += YlmVec[lm] * 
	(-(double)(l*(l+1))/(rmag*rmag) * ulm[lm]
	 + d2ulm[lm] + 2.0/rmag *dulm[lm]);
      lm++;
  }
}


void 
AtomicOrbital::eval (Vec3 r,
		     complex<double> &val,
		     TinyVector<complex<double>,3> &grad,
		     complex<double> &lapl)
{
  Vec3 diff = r - Pos;
  double rmag = sqrt (dot(diff,diff));
  val = lapl = complex<double>();
  grad = TinyVector<complex<double>,3>();
  if (rmag > Radius) 
    return;

  if (rmag < PolyRadius) {
    eval_poly (r, val, grad, lapl); 
    return;
  }

  Vec3 rhat, thetahat, phihat;
  rhat = (1.0/rmag) * diff;

  double costheta = rhat[2];
  double sintheta = std::sqrt(1.0-costheta*costheta);
  double cosphi = rhat[0]/sintheta;
  double sinphi = rhat[1]/sintheta;
  thetahat = TinyVector<double,3>(costheta*cosphi,
				  costheta*sinphi,
				  -sintheta);
  phihat   = TinyVector<double,3>(-sinphi,
				  cosphi,
				  0.0 );

  CalcYlm (rhat);

  int numlm = YlmVec.size();
  complex<double> ulm[numlm], dulm[numlm], d2ulm[numlm];

  eval_multi_UBspline_1d_z_vgl (Spline, rmag, 
				&(ulm[0]), &dulm[0], &d2ulm[0]);

  val = lapl = complex<double>();
  grad = TinyVector<complex<double>,3>();
  int lm =0;
  double r2l=1.0, r2lm1=0.0, r2lm2=0.0;

  for (int l=0; l<=lMax; l++) {
    double dl = (double)l;
    for (int m=-l; m<=l; m++) {
      complex<double> u, du, d2u;
      
      u   = r2l * ulm[lm];
      du  = dl*r2lm1*ulm[lm] + r2l*dulm[lm];
      d2u = ((double)(l*(l-1))*r2lm2*ulm[lm]
	     + 2.0*dl*r2lm1*dulm[lm]
	     + r2l*d2ulm[lm]);
      
      complex<double> im(0.0,(double)m);

      // val += u * YlmVec[lm];

      // grad += 
      // 	(du                *     YlmVec[lm] * rhat     +
      // 	 u/rmag            *    dYlmVec[lm] * thetahat +
      // 	 u/(rmag*sintheta) * im *YlmVec[lm] * phihat);

      // lapl += YlmVec[lm] * 
      // 	(-(double)(l*(l+1))/(rmag*rmag) * u
      // 	 + d2u + 2.0/rmag *du);

      val += ulm[lm] * YlmVec[lm];
      grad += 
      	(dulm[lm]                *     YlmVec[lm] * rhat     +
      	 ulm[lm]/rmag            *    dYlmVec[lm] * thetahat +
      	 ulm[lm]/(rmag*sintheta) * im *YlmVec[lm] * phihat);

      lapl += YlmVec[lm] * 
      	(-(double)(l*(l+1))/(rmag*rmag) * ulm[lm]
      	 + d2ulm[lm] + 2.0/rmag *dulm[lm]);

      lm++;
    }
    r2lm2 = r2lm1;
    r2lm1 = r2l;
    r2l *= rmag;
  }
}

void
AtomicOrbital::Write_ESHDF (IO::IOSectionClass &out)
{
  out.WriteVar ("position", Pos);
  out.WriteVar ("cutoff_radius", Radius);
  out.WriteVar ("spline_radius", OuterRadius);
  out.WriteVar ("spline_points", SplinePoints);
  out.WriteVar ("lmax", lMax);
  out.WriteVar ("polynomial_order", PolyOrder);
  out.WriteVar ("polynomial_radius", PolyRadius);
}



#ifdef HAVE_CUDA


void ProjectAnalyticCuda(int numats, int lMax, int SplinePoints, int PolyOrder, double PolyRadius, 
			 double OuterRadius, double* positions, 
			 TinyVector<Array<OrbitalClass*,2>,2>& orbitals) {
  int numSpins = 1 + (orbitals[1].size() > 0);
  int num_twists = orbitals[0].extent(0);
  int num_bands = orbitals[0].extent(1);
  int numlm = (lMax+1)*(lMax+1);

  for (int iat = 0; iat < numats; iat++) {
    for (int spin=0; spin< numSpins; spin++) {
      for (int ti=0; ti<num_twists; ti++) {
	for (int bi=0; bi<num_bands; bi++) {
	  orbitals[spin](ti,bi)->AtomicData[iat].SplineData.resize(SplinePoints, numlm);
	  orbitals[spin](ti,bi)->AtomicData[iat].SplineData = 0.0;
	  orbitals[spin](ti,bi)->AtomicData[iat].PolyCoefs.resize(PolyOrder+1, numlm);
	}
      }
    }
  }
  
  const int PolyPoints = 32;
  double poly_dr = 1.5*PolyRadius / (PolyPoints-1);   
  double dr = OuterRadius / (SplinePoints -1);
  
  double* PolyData = new double[2*num_bands*numSpins*PolyPoints*numlm*numats];
  for (int ik=0; ik<num_twists; ik++) {
    CellClass& cell = orbitals[0](ik,0)->GetCell();
    int numG = cell.GVecs.size();
    // Need to create array of g-vectors
    double* gvecs = new double[3*numG];
    for (int iG = 0; iG < numG; iG++) {
      gvecs[iG*3] = cell.GVecs(iG)[0];
      gvecs[iG*3+1] = cell.GVecs(iG)[1];
      gvecs[iG*3+2] = cell.GVecs(iG)[2];
    }    
    /*
    for (int iG = 0; iG < 20; iG++) {
      cout << "G vector " << iG << " = " << gvecs[iG*3] << "   " << gvecs[iG*3+1] << "   " << gvecs[iG*3+2] << endl;
    }
    */
    Vec3 k = orbitals[0](ik,0)->Getk();     
    
    vec3<double>::Type floatk;
    floatk.x = k[0];
    floatk.y = k[1];
    floatk.z = k[2];

    double* orbitalCoefs = new double[2*numSpins*num_bands*numG];
    for (int is = 0; is < numSpins; is++) {
      for (int ib = 0; ib < num_bands; ib++) {
	zVec &coefs = orbitals[is](ik, ib)->GetCoefs();
	for (int iG = 0; iG < numG; iG++) {
	  const int coefindex = 2*(iG+numG*(ib+num_bands*is));
	  orbitalCoefs[coefindex] = coefs(iG).real();
	  orbitalCoefs[coefindex+1] = coefs(iG).imag();
	}
      }
    }
    // Need to create data to hold everything for splineData and polyData
    double* SplineData = new double[2*num_bands*numSpins*SplinePoints*numlm*numats];
    for (int i = 0; i < 2*num_bands*numSpins*SplinePoints*numlm*numats; i++) {
      SplineData[i] = 0.0;
    }
    for (int i = 0; i < 2*num_bands*numSpins*PolyPoints*numlm*numats; i++) {
      PolyData[i] = 0.0;
    }
    
    projectAnalyticCudaDriver(floatk, positions, gvecs, orbitalCoefs, SplineData,
			      PolyData, ik, num_twists, numG, lMax, numats, SplinePoints, 
			      PolyPoints, numSpins, num_bands, dr, poly_dr);
    
    cout << "   About to copy SplineData" << endl;
    for (int iat = 0; iat < numats; iat++) {
      for (int is = 0; is < numSpins; is++) {
	for (int ib = 0; ib < num_bands; ib++) {
	  Array<complex<double>,2>& data = orbitals[is](ik,ib)->AtomicData[iat].SplineData;
	  for (int ir =0; ir < SplinePoints; ir++) {
	    for (int lm = 0; lm < numlm; lm++) {
	      const int index = 2*(lm+numlm*(ir+SplinePoints*(ib+num_bands*(is+numSpins*iat))));
	      data(ir,lm).real() = SplineData[index];
	      data(ir,lm).imag() = SplineData[index+1];
	    }
	  }
	}
      }	  
    }   
    cout << "   Finished copying SplineData" << endl;
    
    delete[] gvecs;
    delete[] orbitalCoefs;
    delete[] SplineData;
 
    // Now do Polynomial Fits
    for (int iat = 0; iat < numats; iat++) {
      Array<double,1> y_re(PolyPoints), y_im(PolyPoints), sigma(PolyPoints);
      Array<double,1> coefs_re(PolyOrder), errors_re(PolyOrder),
	coefs_im(PolyOrder), errors_im(PolyOrder);
      
      for (int spin=0; spin<numSpins; spin++) {
	for (int ib=0; ib<num_bands; ib++) {
	  Array<complex<double>,2> &coefs = orbitals[spin](ik, ib)->AtomicData[iat].PolyCoefs;
	  coefs.resize(PolyOrder+1, (lMax+1)*(lMax+1));
	  coefs = complex<double>();
	  
	  for (int l=0; l<=lMax; l++) {
	    coefs_re.resize(PolyOrder-l+1);    
	    coefs_im.resize(PolyOrder-l+1);
	    Array<double,2> Basis(PolyPoints, PolyOrder-l+1);
	    for (int ir=0; ir<PolyPoints; ir++) {
	      double r=poly_dr*ir;
	      double r2lpn=1.0;
	      for (int i=0; i<l; i++) {
		r2lpn *= r;
	      }
	      for (int n=0; n<PolyOrder-l+1; n++) {
		Basis(ir,n) = r2lpn;
		r2lpn *= r;
	      }
	    }
	    
	    for (int m=-l; m<=l; m++) {
	      int lm = l*(l+1)+m;
	      for (int ir=0; ir<PolyPoints; ir++) {
		const int index = 2*(lm+numlm*(ir+PolyPoints*(ib+num_bands*(spin+numSpins*iat))));
		y_re(ir) = PolyData[index];
		y_im(ir) = PolyData[index+1];
		sigma(ir) = 1.0e-40;
	      }
	      
	      LinFitLU (y_re, sigma, Basis, coefs_re, errors_re);
	      LinFitLU (y_im, sigma, Basis, coefs_im, errors_im);
	      for (int n=0; n<l; n++)
		coefs(n,lm) = complex<double>();
	      for (int n=l; n<PolyOrder; n++)
		coefs(n,lm) = complex<double>(coefs_re(n-l), coefs_im(n-l));
	    }
	  }
	}
      }
    }
  }
  delete[] PolyData;
}


void 
AtomicOrbital::ProjectAnalyticCudaDebug (int iat, TinyVector<Array<OrbitalClass*,2>,2> &orbitals)
{
  int numSpins   = 1 + (orbitals[1].size() > 0);
  int num_twists = orbitals[0].extent(0);
  int num_bands  = orbitals[0].extent(1);
  int numlm      = (lMax+1)*(lMax+1);

  // Resize data in all the orbitals
  for (int spin=0; spin<numSpins; spin++) {
    for (int ti=0; ti<num_twists; ti++)
      for (int bi=0; bi<num_bands; bi++) {
	orbitals[spin](ti,bi)->AtomicData[iat].SplineData.resize(SplinePoints, numlm);
	orbitals[spin](ti,bi)->AtomicData[iat].SplineData = 0.0;
	orbitals[spin](ti,bi)->AtomicData[iat].PolyCoefs.resize(PolyOrder+1, numlm);
      }
  }

  const int PolyPoints = 32;
  Array<complex<double>,5> poly_data(numSpins,num_twists,num_bands,PolyPoints, numlm);
  poly_data = complex<double>();
  double poly_dr = 1.5*PolyRadius / (PolyPoints-1);

  double dr = OuterRadius / (SplinePoints - 1);
  const int G_block_size = 65536;


  for (int ik=0; ik<num_twists; ik++) {
    Array<double,3> jlArray_spline(G_block_size, SplinePoints, lMax+1);
    Array<double,3> jlArray_poly  (G_block_size, PolyPoints, lMax+1);
    // We include the 4 pi * i^l term in this array
    Array<complex<double>,2> YlmArray(G_block_size, numlm);
    //      vector<double> jlVec(lMax+1);
    vector<complex<double> > phase_shift(G_block_size);
    
    Vec3 k = orbitals[0](ik,0)->Getk();
    cerr << "k = " << k << endl;
    CellClass &cell = orbitals[0](ik,0)->GetCell();
    int numG = cell.GVecs.size();
    int num_G_blocks = (numG+G_block_size-1) / G_block_size;
    
    typedef double Tp;
    typedef double Tp2;
    Tp2* orbitalCoefs = new Tp2[2*numSpins*num_bands*numG];
    for (int is = 0; is < numSpins; is++) {
      for (int ib = 0; ib < num_bands; ib++) {
	zVec &coefs = orbitals[is](ik, ib)->GetCoefs();
	for (int iG = 0; iG < numG; iG++) {
	  const int coefindex = 2*(iG+numG*(ib+num_bands*is));
	  orbitalCoefs[coefindex] = coefs(iG).real();
	  orbitalCoefs[coefindex+1] = coefs(iG).imag();
	}
      }
    }
    //cout << "Finished setting up array for orbitalCoefs" << endl;
    // Need to create data to hold everything for splineData and polyData
    Tp2* SplineData = new Tp2[2*num_bands*numSpins*SplinePoints*numlm];
    Tp2* PolyData = new Tp2[2*num_bands*numSpins*PolyPoints*numlm];
    for (int i = 0; i < 2*num_bands*numSpins*SplinePoints*numlm; i++) {
      SplineData[i] = 0.0;
    }
    for (int i = 0; i < 2*num_bands*numSpins*PolyPoints*numlm; i++) {
      PolyData[i] = 0.0;
    }
    

    Tp* gvecs = new Tp[3*numG];
    for (int iG = 0; iG < numG; iG++) {
      gvecs[iG*3] = cell.GVecs(iG)[0];
      gvecs[iG*3+1] = cell.GVecs(iG)[1];
      gvecs[iG*3+2] = cell.GVecs(iG)[2];
    }
    
    Tp* floatgvecs = new Tp[3*G_block_size];
    Tp* floatPhaseShifts = new Tp[2*G_block_size];
    Tp* floatjlArraySpline = new Tp[G_block_size*SplinePoints*(lMax+1)];
    Tp* floatjlArrayPoly = new Tp[G_block_size*PolyPoints*(lMax+1)];
    Tp* floatYlmarr = new Tp[2*G_block_size*numlm];
    
    
    for (int Gblock=0; Gblock < num_G_blocks; Gblock++) {
      cerr << "Gblock = " << Gblock << endl;
      int offset = Gblock * G_block_size;
      int end = min(numG, (Gblock+1)*G_block_size) - offset;
      
      
      vec3<Tp>::Type floatk;
      floatk.x = k[0];
      floatk.y = k[1];
      floatk.z = k[2];
      
      vec3<Tp>::Type floatpos;
      floatpos.x = Pos[0];
      floatpos.y = Pos[1];
      floatpos.z = Pos[2];
      
      for (int i = 0; i < end; i++) {
	//cout << "Gvec[" << i+offset << "] = " << gvecs[3*(i+offset)] << "   " << gvecs[3*(i+offset)+1] << "   " << gvecs[3*(i+offset)+2] << endl;
	floatgvecs[3*i] = gvecs[3*(i+offset)];
	floatgvecs[3*i+1] = gvecs[3*(i+offset)+1];
	floatgvecs[3*i+2] = gvecs[3*(i+offset)+2];
      }
      
      //cout << "About to call setupGPUComputeKernel" << endl;
      //setupGPUComputeKernelDebug(floatk, pos, floatgvecs, floatPhaseShifts, floatjlArraySpline,
      //		           floatjlArrayPoly, floatYlmarr, end, numats, SplinePoints, lMax,
      //		           PolyPoints, dr, poly_dr);
      //      cout << "Finished calling setupGPUComputeKernel" << endl;
      
      
      cout << "About to call setupCPUComputeKernel" << endl;
      setupCPUComputeKernel(floatk, floatpos, floatgvecs, floatPhaseShifts, floatjlArraySpline,
			    floatjlArrayPoly, floatYlmarr, end, SplinePoints, lMax,
			    PolyPoints, dr, poly_dr);
      cout << "finished calling setupCPUComputeKernel" << endl;

      //cout << "About to call setupCPUComputeKernelNew" << endl;
      //setupCPUComputeKernelNew(floatk, floatpos, floatgvecs, floatPhaseShifts, floatjlArraySpline,
      //			       floatjlArrayPoly, floatYlmarr, end, numats, SplinePoints, lMax,
      //			       PolyPoints, dr, poly_dr);
      //cout << "finished calling setupCPUComputeKernelNew" << endl;
      
      for (int iG = 0; iG < end; iG++) {
	phase_shift[iG] = complex<double> (floatPhaseShifts[2*iG], floatPhaseShifts[2*iG+1]);
	for (int ir=0; ir<SplinePoints; ir++) {
	  for (int l=0; l<=lMax; l++) {
	    //const int asloc = l+(lMax+1)*(ir+SplinePoints*(iG));
	    const int asloc = iG+end*(l+(lMax+1)*ir);
	    jlArray_spline(iG, ir, l) = floatjlArraySpline[asloc];
	  }
	}
	for (int ir=0; ir<PolyPoints; ir++) {
	  for (int l=0; l<=lMax; l++) {
	    //const int aploc = l+(lMax+1)*(ir+PolyPoints*(iG));
	    const int aploc = iG+end*(l+(lMax+1)*ir);
	    jlArray_poly(iG, ir, l) = floatjlArrayPoly[aploc];
	  }
	}
	for (int l=0,lm=0; l<=lMax; l++) {
	  for (int m=-l; m<=l; m++, lm++) {
	    const int ylmloc = 2*(iG+end*lm);
	    //const int ylmloc = (lMax+1)*(lMax+1)*2*iG+ 2*lm;
	    YlmArray(iG, lm) = complex<double>(floatYlmarr[ylmloc],floatYlmarr[ylmloc+1]); 
	  }
	}
      }
      


      Tp2* floatorbs = new Tp2[numSpins*G_block_size*num_bands*2];
      for (int ib = 0; ib < num_bands; ib++) {
	for (int sigma = 0; sigma < numSpins; sigma++) {
	  for (int j = 0; j < 2*end; j++) {
	    floatorbs[j+2*G_block_size*(ib+num_bands*sigma)] = orbitalCoefs[2*offset+j+2*numG*(ib+num_bands*sigma)];
	  }
	}
      }
      Tp2* float2PhaseShifts = new Tp2[2*G_block_size];
      Tp2* float2jlArraySpline = new Tp2[G_block_size*SplinePoints*(lMax+1)];
      Tp2* float2jlArrayPoly = new Tp2[G_block_size*PolyPoints*(lMax+1)];
      Tp2* float2Ylmarr = new Tp2[2*G_block_size*numlm];
      for (int i = 0; i < 2*G_block_size; i++) {
	float2PhaseShifts[i] = floatPhaseShifts[i];
      }
      for (int i = 0; i < G_block_size*SplinePoints*(lMax+1); i++) {
	float2jlArraySpline[i] = floatjlArraySpline[i];
      }
      for (int i = 0; i < G_block_size*PolyPoints*(lMax+1); i++) {
	float2jlArrayPoly[i] = floatjlArrayPoly[i];
      }
      for (int i = 0; i < 2*G_block_size*numlm; i++) {
	float2Ylmarr[i] = floatYlmarr[i];
      }

      
      cout << "About to call projectAnalyticKernel" << endl;
      
      projectAnalyticCPUKernel(floatorbs, float2PhaseShifts, float2jlArraySpline,
			       float2jlArrayPoly, float2Ylmarr, SplineData, PolyData,
			       numSpins, num_bands, G_block_size, end, SplinePoints,
			       PolyPoints, lMax);
      
      
      
      //projectAnalyticCUDAKernel(floatorbs, float2PhaseShifts, float2jlArraySpline,
      //				float2jlArrayPoly, float2Ylmarr, SplineData, PolyData,
      //				numSpins, num_bands, G_block_size, end, SplinePoints,
      //				PolyPoints, lMax);
      //      cout << "Finished calling projectAnalyticKernel" << endl;
      

      delete[] float2PhaseShifts;
      delete[] float2jlArraySpline;
      delete[] float2jlArrayPoly;
      delete[] float2Ylmarr;
    }
    //cout << "   Finished projectAnalyticCudaDriver" << endl;
    //cout << "   About to copy over poly_data" << endl;
    // Need to copy proper data back from polyData to poly_data 
    for (int is = 0; is < numSpins; is++) {
      for (int ib = 0; ib < num_bands; ib++) {
	for (int ir = 0; ir < PolyPoints; ir++) {
	  for (int lm = 0; lm < numlm; lm++) {
	    const int index = 2*(lm+numlm*(ir+PolyPoints*(ib+num_bands*is)));
	    poly_data(is,ik,ib,ir,lm).real() = PolyData[index];
	    poly_data(is,ik,ib,ir,lm).imag() = PolyData[index+1];
	  }
	}
      }
    }
    //cout << "   Finished copying poly_data" << endl;
    // splineData to  Array<complex<double>,2> &data = orbitals[spin](ik, ib)->AtomicData[iat].SplineData;
    //cout << "   About to copy SplineData" << endl;
    for (int is = 0; is < numSpins; is++) {
      for (int ib = 0; ib < num_bands; ib++) {
	Array<complex<double>,2>& data = orbitals[is](ik,ib)->AtomicData[iat].SplineData;
	for (int ir =0; ir < SplinePoints; ir++) {
	  for (int lm = 0; lm < numlm; lm++) {
	    const int index = 2*(lm+numlm*(ir+SplinePoints*(ib+num_bands*is)));
	    //	     cout << "(is,ib,ir,lm) = (" << is << "," << ib << "," << ir << "," << lm << "),  index = " << index << endl;
	    data(ir,lm).real() = SplineData[index];
	    data(ir,lm).imag() = SplineData[index+1];
	  }
	}
      }
    }	     
    //cout << "   Finished copying SplineData" << endl;
	

    // Now, do polynomial fits
    Array<double,1> y_re(PolyPoints), y_im(PolyPoints), sigma(PolyPoints);
    Array<double,1> coefs_re(PolyOrder), errors_re(PolyOrder),
      coefs_im(PolyOrder), errors_im(PolyOrder);
    
    for (int spin=0; spin<numSpins; spin++)
      for (int ib=0; ib<num_bands; ib++) {
	Array<complex<double>,2> &coefs = orbitals[spin](ik, ib)->AtomicData[iat].PolyCoefs;
	coefs.resize(PolyOrder+1, (lMax+1)*(lMax+1));
	coefs = complex<double>();
	
	for (int l=0; l<=lMax; l++) {
	  coefs_re.resize(PolyOrder-l+1);    coefs_im.resize(PolyOrder-l+1);
	  Array<double,2> Basis(PolyPoints, PolyOrder-l+1);
	  for (int ir=0; ir<PolyPoints; ir++) {
	    double r=poly_dr*ir;
	    double r2lpn=1.0;
	    for (int i=0; i<l; i++)	r2lpn *= r;
	    for (int n=0; n<PolyOrder-l+1; n++) {
	      Basis(ir,n) = r2lpn;
	      r2lpn *= r;
	    }
	  }
	  for (int m=-l; m<=l; m++) {
	    int lm = l*(l+1)+m;
	    for (int ir=0; ir<PolyPoints; ir++) {
	      y_re(ir) = poly_data(spin,ik,ib,ir,lm).real();
	      y_im(ir) = poly_data(spin,ik,ib,ir,lm).imag();
	      sigma(ir) = 1.0e-40;
	    }
	    LinFitLU (y_re, sigma, Basis, coefs_re, errors_re);
	    LinFitLU (y_im, sigma, Basis, coefs_im, errors_im);
	    
	    for (int n=0; n<l; n++)
	      coefs(n,lm) = complex<double>();
	    for (int n=l; n<PolyOrder; n++)
	      coefs(n,lm) = complex<double>(coefs_re(n-l), coefs_im(n-l));
	    
	  }
	}
      }
  }
}


void
projectAnalyticCPUKernel(const double* const orbitalCoefs, const double* const phase_shifts, 
                          const double* const jlArray_spline, const double* const jlArray_poly,
			  const double* const YlmPtr, double* splineData, double* polyData, 
                          int numSpins, int num_bands, int block_size, int numG, int SplinePoints, int PolyPoints, 
                          int lmax) {
  const int lmmax = (lmax+1)*(lmax+1);
  for (int is = 0; is < numSpins; is++) {
    for (int ib = 0; ib < num_bands; ib++) {
      for (int ir = 0; ir < SplinePoints; ir++) {
	for (int l=0, lm=0; l<=lmax; l++) {
	  for (int m=-l; m<=l; m++,lm++) {
	    const int sdIndex = 2*(lm+lmmax*(ir+SplinePoints*(ib+num_bands*is)));
	    // try kahan summation
	    double reregister = 0.0;
	    double recomp = 0.0;
	    double imregister = 0.0;
	    double imcomp = 0.0;
	    for (int iG = 0; iG < numG; iG++) {
	      const int orbindex = 2*(iG+block_size*(ib+num_bands*is));
	      const double creal = orbitalCoefs[orbindex]*phase_shifts[2*iG] - orbitalCoefs[orbindex+1]*phase_shifts[2*iG+1];
	      const double cimag = orbitalCoefs[orbindex]*phase_shifts[2*iG+1] + orbitalCoefs[orbindex+1]*phase_shifts[2*iG];
	      //const int jlasindex = l+(lmax+1)*(ir+SplinePoints*iG);
	      const int jlasindex = iG+numG*(l+(lmax+1)*ir);
	      const double cjlreal = creal*jlArray_spline[jlasindex];
	      const double cjlimag = cimag*jlArray_spline[jlasindex];
	      //const int YlmIndex = 2*(iG*lmmax+lm);
	      const int YlmIndex = 2*(iG+numG*lm);

	      const double yvalre = cjlreal*YlmPtr[YlmIndex]-cjlimag*YlmPtr[YlmIndex+1] - recomp;
	      const double tvalre = reregister + yvalre;
	      recomp = (tvalre - reregister) - yvalre;
	      reregister = tvalre;
	      //reregister += cjlreal*YlmPtr[YlmIndex]-cjlimag*YlmPtr[YlmIndex+1];
	      
	      const double yvalim = cjlreal*YlmPtr[YlmIndex+1]+cjlimag*YlmPtr[YlmIndex] - imcomp;
	      const double tvalim = imregister + yvalim;
	      imcomp = (tvalim - imregister) - yvalim;
	      imregister = tvalim;
	      //imregister += cjlreal*YlmPtr[YlmIndex+1]+cjlimag*YlmPtr[YlmIndex];
	    }
	    splineData[sdIndex] += reregister;
	    splineData[sdIndex+1] += imregister;
	  }
	}
      }
      for (int ir = 0; ir < PolyPoints; ir++) {
	for (int l=0, lm=0; l<=lmax; l++) {
	  for (int m=-l; m<=l; m++,lm++) {
	    const int polyIndex = 2*(lm+lmmax*(ir+PolyPoints*(ib+num_bands*is)));
	    // Try kahan summation
	    double reregister = 0.0;
	    double recomp = 0.0;
	    double imregister = 0.0;
	    double imcomp = 0.0;
	    for (int iG = 0; iG < numG; iG++) {
	      const int orbindex = 2*(iG+block_size*(ib+num_bands*is));
	      const double creal = orbitalCoefs[orbindex]*phase_shifts[2*iG] - orbitalCoefs[orbindex+1]*phase_shifts[2*iG+1];
	      const double cimag = orbitalCoefs[orbindex]*phase_shifts[2*iG+1] + orbitalCoefs[orbindex+1]*phase_shifts[2*iG];
	      //const int jlapindex = l+(lmax+1)*(ir+PolyPoints*iG);
	      const int jlapindex = iG+numG*(l+(lmax+1)*ir);
	      const double cjlreal = creal*jlArray_poly[jlapindex];
	      const double cjlimag = cimag*jlArray_poly[jlapindex];
	      //const int YlmIndex = 2*(iG*lmmax+lm);
	      const int YlmIndex = 2*(iG+numG*lm);
	      
	      const double yvalre = cjlreal*YlmPtr[YlmIndex]-cjlimag*YlmPtr[YlmIndex+1] - recomp;
	      const double tvalre = reregister + yvalre;
	      recomp = (tvalre - reregister) - yvalre;
	      reregister = tvalre;
	      //reregister += cjlreal*YlmPtr[YlmIndex]-cjlimag*YlmPtr[YlmIndex+1];
	      
	      const double yvalim = cjlreal*YlmPtr[YlmIndex+1]+cjlimag*YlmPtr[YlmIndex] - imcomp;
	      const double tvalim = imregister + yvalim;
	      imcomp = (tvalim - imregister) - yvalim;
	      imregister = tvalim;
	      //imregister += cjlreal*YlmPtr[YlmIndex+1]+cjlimag*YlmPtr[YlmIndex];
	    }
	    polyData[polyIndex] += reregister;
	    polyData[polyIndex+1] += imregister;
	  }
	}
      }    
    }
  }
};

void
setupCPUComputeKernel(const typename vec3<double>::Type k, const typename vec3<double>::Type pos,
		      double* Ghats, double* phase_shifts, double* jlArray_spline,
		      double* jlArray_poly, double* ylmptr,  int numG, int SplinePoints, int lmax, int PolyPoints,
		      double dr, double poly_dr) {
  double jlVec[lmax+1];
  for (int gn = 0; gn < numG; gn++) {
    //    cout << "gn = " << gn << endl;
    // Compute Phase Shift
    Ghats[3*gn] += k.x;
    Ghats[3*gn+1] += k.y;
    Ghats[3*gn+2] += k.z;
    
    double val = -Ghats[3*gn]*pos.x -Ghats[3*gn+1]*pos.y -Ghats[3*gn+2]*pos.z;
    phase_shifts[gn*2] = cos(val);
    phase_shifts[gn*2+1] = sin(val);
    
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
	jlArray_spline[gn+numG*(l+(lmax+1)*ir)] = jlVec[l];
      }
    }
    for (int ir = 0; ir<PolyPoints; ir++) {
      double r = poly_dr * ir;
      gsl_sf_bessel_jl_steed_array(lmax, Gmag*r, jlVec);
      for (int l = 0; l <=lmax; l++) {
	//jlArray_poly[l+(lmax+1)*(ir+PolyPoints*gn)] = jlVec[l];
	jlArray_poly[gn+numG*(l+(lmax+1)*ir)] = jlVec[l];
      }
    }
  }
  
  cout << "Finished first loop over iG" << endl;
  for (int gn = 0; gn < numG; gn++) {
    const double fourPiInv = 0.0795774715459477;
    int maxindex = (lmax+1)*(lmax+1)*2;
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
}    


#endif

void 
AtomicOrbital::ProjectAnalytic (int iat, TinyVector<Array<OrbitalClass*,2>,2> &orbitals)
{
  int numSpins   = 1 + (orbitals[1].size() > 0);
  int num_twists = orbitals[0].extent(0);
  int num_bands  = orbitals[0].extent(1);
  int numlm      = (lMax+1)*(lMax+1);

  // Resize data in all the orbitals
  for (int spin=0; spin<numSpins; spin++) {
    for (int ti=0; ti<num_twists; ti++)
      for (int bi=0; bi<num_bands; bi++) {
	orbitals[spin](ti,bi)->AtomicData[iat].SplineData.resize(SplinePoints, numlm);
	orbitals[spin](ti,bi)->AtomicData[iat].SplineData = 0.0;
	orbitals[spin](ti,bi)->AtomicData[iat].PolyCoefs.resize(PolyOrder+1, numlm);
      }
  }

  const int PolyPoints = 32;
  Array<complex<double>,5> poly_data(numSpins,num_twists,num_bands,PolyPoints, numlm);
  poly_data = complex<double>();
  double poly_dr = 1.5*PolyRadius / (PolyPoints-1);

  double dr = OuterRadius / (SplinePoints - 1);
  const int G_block_size = 65536;
  // Loop over k-points
//#pragma omp parallel
  {
//#pragma omp parallel for schedule(dynamic,1)
    for (int ik=0; ik<num_twists; ik++) {
      Array<double,3> jlArray_spline(G_block_size, SplinePoints, lMax+1);
      Array<double,3> jlArray_poly  (G_block_size, PolyPoints, lMax+1);
      // We include the 4 pi * i^l term in this array
      Array<complex<double>,2> YlmArray(G_block_size, numlm);
      //      vector<double> jlVec(lMax+1);
      vector<complex<double> > phase_shift(G_block_size);

      Vec3 k = orbitals[0](ik,0)->Getk();
      cerr << "k = " << k << endl;
      CellClass &cell = orbitals[0](ik,0)->GetCell();
      int numG = cell.GVecs.size();
      int num_G_blocks = (numG+G_block_size-1) / G_block_size;
      for (int Gblock=0; Gblock < num_G_blocks; Gblock++) {
	cerr << "Gblock = " << Gblock << endl;
	int offset = Gblock * G_block_size;
	int end = min(numG, (Gblock+1)*G_block_size) - offset;
	//#pragma omp parallel for schedule(dynamic,1)
	for (int iG=0; iG < end; iG++) {
	  vector<double> jlVec(lMax+1);
	  Vec3 G = cell.GVecs(iG+offset) + k;
	  double Gmag = sqrt(dot(G,G));
	  Vec3 Ghat;
	  if (Gmag > 0.0)
	    Ghat = (1.0/Gmag)*G;
	  else
	    Ghat = Vec3(0.0, 0.0, 1.0);
	  double s, c;
	  sincos (-dot(G,Pos), &s, &c);
	  phase_shift[iG] = complex<double> (c, s);
	  
	  // Compute bessel functions
	  for (int ir=0; ir<SplinePoints; ir++) {
	    double r = dr * ir;
	    gsl_sf_bessel_jl_steed_array(lMax, Gmag*r, &(jlVec[0]));
	    for (int l=0; l<=lMax; l++)
	      jlArray_spline(iG, ir, l) = jlVec[l];
	  }
	  for (int ir=0; ir<PolyPoints; ir++) {
	    double r = poly_dr * ir;
	    gsl_sf_bessel_jl_steed_array(lMax, Gmag*r, &(jlVec[0]));
	    for (int l=0; l<=lMax; l++)
	      jlArray_poly(iG, ir, l) = jlVec[l];
	  }
	
	  // Compute Ylms
	  CalcYlm(Ghat);
	  complex<double> minus_i2l(1.0,0.0);
	  complex<double> minus_i(0.0, -1.0);
	  for (int l=0,lm=0; l<=lMax; l++) {
	    for (int m=-l; m<=l; m++, lm++) 
	      YlmArray(iG, lm) = 4.0*M_PI*minus_i2l * conj(YlmVec[lm]);
	    minus_i2l *= minus_i;
	  }
	}
	 for (int spin=0; spin<numSpins; spin++) {
#pragma omp parallel for schedule(dynamic,1)
	    for (int ib=0; ib<num_bands; ib++) {
	       Array<complex<double>,2> &data = orbitals[spin](ik, ib)->AtomicData[iat].SplineData;
	       zVec &coefs = orbitals[spin](ik, ib)->GetCoefs();
	       for (int iG=0; iG<end; iG++) {
		  complex<double> c = coefs(iG+offset) * phase_shift[iG];
		  for (int ir=0; ir<SplinePoints; ir++)
		    for (int l=0,lm=0; l<=lMax; l++) {
		       complex<double> cjl = c*jlArray_spline(iG, ir, l);
		       for (int m=-l; m<=l; m++,lm++) 
			 data(ir,lm) += cjl * YlmArray(iG,lm);
		    }
		  for (int ir=0; ir<PolyPoints; ir++)
		    for (int l=0,lm=0; l<=lMax; l++) {
		       complex<double> cjl = c*jlArray_poly(iG, ir, l);
		       for (int m=-l; m<=l; m++,lm++) 
			 poly_data(spin, ik, ib, ir,lm) += cjl * YlmArray(iG,lm);
		    }
	       }
	    }
	 }
      }
       
       
       // Now, do polynomial fits
       Array<double,1> y_re(PolyPoints), y_im(PolyPoints), sigma(PolyPoints);
       Array<double,1> coefs_re(PolyOrder), errors_re(PolyOrder),
	 coefs_im(PolyOrder), errors_im(PolyOrder);
       
       for (int spin=0; spin<numSpins; spin++)
	 for (int ib=0; ib<num_bands; ib++) {
	    Array<complex<double>,2> &coefs = orbitals[spin](ik, ib)->AtomicData[iat].PolyCoefs;
	    coefs.resize(PolyOrder+1, (lMax+1)*(lMax+1));
	    coefs = complex<double>();
	    
	    for (int l=0; l<=lMax; l++) {
	       coefs_re.resize(PolyOrder-l+1);    coefs_im.resize(PolyOrder-l+1);
	       Array<double,2> Basis(PolyPoints, PolyOrder-l+1);
	       for (int ir=0; ir<PolyPoints; ir++) {
		  double r=poly_dr*ir;
		  double r2lpn=1.0;
		  for (int i=0; i<l; i++)	r2lpn *= r;
		  for (int n=0; n<PolyOrder-l+1; n++) {
		     Basis(ir,n) = r2lpn;
		     r2lpn *= r;
		  }
	       }
	       for (int m=-l; m<=l; m++) {
		  int lm = l*(l+1)+m;
		  for (int ir=0; ir<PolyPoints; ir++) {
		     y_re(ir) = poly_data(spin,ik,ib,ir,lm).real();
		     y_im(ir) = poly_data(spin,ik,ib,ir,lm).imag();
		     sigma(ir) = 1.0e-40;
		  }
		  LinFitLU (y_re, sigma, Basis, coefs_re, errors_re);
		  LinFitLU (y_im, sigma, Basis, coefs_im, errors_im);
		  
		  for (int n=0; n<l; n++)
		    coefs(n,lm) = complex<double>();
		  for (int n=l; n<PolyOrder; n++)
		    coefs(n,lm) = complex<double>(coefs_re(n-l), coefs_im(n-l));
	       
	       }
	    }
	 }
    }
     
     
     
  }
}
