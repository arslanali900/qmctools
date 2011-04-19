#ifndef PBC_H
#define PBC_H

#include "LatticeClass.h"

inline void
Minkowski_reduce (Vec3 a[3])
{
  bool found = true;
  int iter = 0;
  while (found && iter < 1000) {
    found = false;
    
    Vec3 temp;
    // sort the vectors
    if (dot(a[1],a[1]) < dot(a[0],a[0])) {
      temp = a[0];
      a[0] = a[1];
      a[1] = temp;
    }
    if (dot(a[2],a[2]) < dot(a[1],a[1])) {
      temp = a[1];
      a[1] = a[2];
      a[2] = temp;
    }
    if (dot(a[1],a[1]) < dot(a[0],a[0])) {
      temp = a[0];
      a[0] = a[1];
      a[1] = temp;
    }

    double l2[3];
    l2[0] = dot(a[0], a[0]);
    l2[1] = dot(a[1], a[1]);
    l2[2] = dot(a[2], a[2]);
    
    // Make all combinations
    for (int i=-1; i<=1; i++)
      for (int j=-1; j<=1; j++)
	for (int k=-1; k<=1; k++) {
	  Vec3 comb = (double)i*a[0] + (double)j*a[1] + (double)k*a[2];
	  double len2 = dot(comb,comb);
	  if (len2 < l2[2] && k!=0) {
	    found = true;
	    a[2] = comb;
	  }
	  else if (len2 < l2[1] && j!= 0) {
	    found = true;
	    a[1] = comb;
	  }
	  else if (len2 < l2[0] && i!=0) {
	    found = true;
	    a[0] = comb;
	  }
	  if (found)
	    break;
	}
    iter++;
  }
}


inline Vec3
MinImage1(Vec3 r, LatticeClass &lattice)
{
  Vec3 u = lattice.r2u(r);
  u[0] -= round(u[0]);
  u[1] -= round(u[1]);
  u[2] -= round(u[2]);

  Vec3 u2 = u;

  
  Vec3 rMin = r;
  double r2Min = dot(r,r);
  TinyVector<int,3> iMin(0,0,0);

  for (int i=-1; i<=1; i++) {
    u2[0] = u[0] + (double)i;
    for (int j=-1; j<=1; j++) {
      u2[1] = u[1] + (double)j;
      for (int k=-1; k<=1; k++) {
	u2[2] = u[2] + (double)k;
	Vec3 rtry = lattice.u2r(u2);
	double r2 = dot(rtry, rtry);
	//	cerr << "rtry = " << rtry  << "  r2 = " << r2 << endl;
	if (r2 < r2Min) {
	  r2Min = r2;
	  rMin = rtry;
	  iMin = TinyVector<int,3>(i,j,k);
	}
      }
    }
  }
  if ((iMin[0]!=0) && (iMin[1]!=0) && (iMin[2]!=0)) {
    cerr << "u = " << u << endl;
    cerr << "iMin = " << iMin << endl;
  }
  return rMin;
}


inline Vec3
MinImage2 (Vec3 r, LatticeClass &lattice)
{
  Vec3 u = lattice.r2u(r);
  u[0] -= round(u[0]);
  u[1] -= round(u[1]);
  u[2] -= round(u[2]);

  Vec3 utry[7], rtry[7];
  double r2[7];
  for (int i=0; i<7; i++)
    utry[i] = u;
  utry[1][0] += 1.0;  utry[2][0] -= 1.0;
  utry[3][1] += 1.0;  utry[4][1] -= 1.0;
  utry[5][2] += 1.0;  utry[6][2] -= 1.0;
  
  for (int i=0; i<7; i++) {
    rtry[i] = lattice.u2r(utry[i]);
    r2[i]   = dot(rtry[i], rtry[i]);
  }
  
  Vec3 rMin = rtry[0];
  double r2Min = r2[0];

  for (int i=1; i<7; i++) {
    if (r2[i] < r2Min) {
      r2Min = r2[i];
      rMin = rtry[i];
    }
  }
  return rMin;
}


inline Vec3
MinImage3 (Vec3 r, LatticeClass &lattice)
{
  Vec3 u = lattice.r2u(r);
  u[0] -= floor(u[0]);
  u[1] -= floor(u[1]);
  u[2] -= floor(u[2]);

  Vec3 utry[8], rtry[8];
  double r2[8];
  for (int i=0; i<8; i++)
    utry[i] = u;
  utry[1] -= Vec3(1.0, 0.0, 0.0);
  utry[2] -= Vec3(0.0, 1.0, 0.0);
  utry[3] -= Vec3(0.0, 0.0, 1.0);
  utry[4] -= Vec3(1.0, 1.0, 0.0);
  utry[5] -= Vec3(1.0, 0.0, 1.0);
  utry[6] -= Vec3(0.0, 1.0, 1.0);
  utry[7] -= Vec3(1.0, 1.0, 1.0);
    
  for (int i=0; i<8; i++) {
    rtry[i] = lattice.u2r(utry[i]);
    r2[i]   = dot(rtry[i], rtry[i]);
  }
  
  Vec3 rMin = rtry[0];
  double r2Min = dot(rtry[0],rtry[0]);

  for (int i=1; i<8; i++) {
    double r2 = dot(rtry[i], rtry[i]);
    if (r2 < r2Min) {
      r2Min = r2;
      rMin = rtry[i];
    }
  }
  return rMin;
}



#endif
