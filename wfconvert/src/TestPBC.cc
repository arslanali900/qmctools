#include "PBC.h"
#include <Common/Random/Random.h>

void
Test()
{
  CommunicatorClass comm;
  RandomClass rand(comm);
  rand.Init();
  const int numLattice = 1000;
  const int numVecs = 10000;

  for (int iL=0; iL<numLattice; iL++) {
    // construct a random lattice
    Vec3 a[3];
    Mat3 aMat;
    for (int i=0; i<3; i++)
      for (int j=0; j<3; j++)
	a[i][j] = 40.0*(rand.Local() - 0.5);
    for (int i=0; i<3; i++)
      for (int j=0; j<3; j++)
	aMat(i,j) = a[i][j];
    cerr << "det1 = " << det(aMat) << endl;
    Minkowski_reduce(a);
    for (int i=0; i<3; i++)
      for (int j=0; j<3; j++)
	aMat(i,j) = a[i][j];
    cerr << "det2 = " << det(aMat) << endl;
    LatticeClass L;
    L.SetDirect(aMat);
    for (int i=0; i<numVecs; i++) {
      Vec3 r;
      for (int j=0; j<3; j++)
	r[j] = 200.0*(rand.Local() - 0.5);
      Vec3 rMin1 = MinImage1(r, L);
      Vec3 rMin2 = MinImage2(r, L);
      Vec3 rMin3 = MinImage3(r, L);
      double r2_1 = dot(rMin1, rMin1);
      double r2_2 = dot(rMin2, rMin2);
      double r2_3 = dot(rMin3, rMin3);
      // if (fabs(r2_1 - r2_2) > 1.0e-10) 
      // 	cerr << "Error in rMin2.\n  rMin1 = " << rMin1 << "\n  rMin2 = "
      // 	     << rMin2 << endl;
      if (fabs(r2_1 - r2_3) > 1.0e-10) 
      	cerr << "Error in rMin3.\n  rMin1 = " << rMin1 << "\n  rMin3 = "
      	     << rMin3 << endl;
    }

  }

}


main()
{
  Test();
}
