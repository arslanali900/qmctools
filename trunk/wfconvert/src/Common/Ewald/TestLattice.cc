#include "LatticeClass.h"

void 
TestrMax()
{
  Mat3 a;
  a(0,0)=0.5;  a(0,1)=0.5; a(0,2)=1.0;
  a(1,0)=0.5;  a(1,1)=1.0; a(1,2)=0.5;
  a(2,0)=1.0;  a(2,1)=0.5; a(2,2)=0.5;
  a = 2.0*8.703 * a;

  LatticeClass lattice(a);

  fprintf (stderr, "rMax = %1.5f\n", lattice.rMax());
  double kcut = sqrt(2.0*200.0);
  TinyVector<int,3> minIndices = lattice.MinIndices(kcut);
  fprintf (stderr, "Min indices for Ecut=200 are (%d,%d,%d)\n",
	   minIndices[0], minIndices[1], minIndices[2]);

  fprintf (stderr, "Volume = %1.4f\n", lattice.Volume());

  Vec3 k(0.0, 0.0, 0.0);
  vector<Vec3> Gvecs = lattice.GenGvecs(k, 200.0);
  fprintf (stderr, "There are %d G-vectors.\n",
	   Gvecs.size());

}


main()
{
  TestrMax();
}
