#include "LatticeClass.h"
#include <vector>

using namespace std;


typedef TinyVector<int,3> Int3;

void
MakekMesh(Mat3 Aprim, Mat3 S, Int3 superMesh,
	  Vec3 superShift, vector<Vec3> &twistVecs)
{
  LatticeClass prim, super;
  prim.SetDirect(Aprim);
  super.SetDirect(S*Aprim);

  fprintf (stderr, "Super lattice vectors:\n");
  for (int i=0; i<3; i++)
    fprintf (stderr, "  [ %9.5f %9.5f %9.5f ]\n", 
	     super.GetDirect()(i,0), 
	     super.GetDirect()(i,1), 
	     super.GetDirect()(i,2));

  Int3 Nmax;
  for (int i=0; i<3; i++) {
    double maxDot = 0.0;
    for (int j=0; j<3; j++)
      maxDot = max(maxDot, (1.0/(2.0*M_PI))*fabs(dot(super.b(i), prim.a(j))));
    Nmax[i] = (int) ceil (1.0/maxDot);
  }

  vector<Vec3> superTwists;

  for (int ns0=0; ns0<superMesh[0]; ns0++) 
    for (int ns1=0; ns1<superMesh[1]; ns1++)
      for (int ns2=0; ns2<superMesh[2]; ns2++) {
	Vec3 shift;
	shift[0] = ((double)ns0+superShift[0])/(double)superMesh[0];
	shift[1] = ((double)ns1+superShift[1])/(double)superMesh[1];
	shift[2] = ((double)ns2+superShift[2])/(double)superMesh[2];

	for (int n0=-Nmax[0]; n0<=Nmax[0]; n0++)
	  for (int n1=-Nmax[1]; n1<=Nmax[1]; n1++)
	    for (int n2=-Nmax[2]; n2<=Nmax[2]; n2++) {
	      Vec3 G = (((double)n0+shift[0])*super.b(0) + 
			((double)n1+shift[1])*super.b(1) +
			((double)n2+shift[2])*super.b(2));
	      // Check if it's in the FBZ of primitive lattice
	      if ((dot(G,prim.a(0)) < 0.5*M_PI) &&
		  (dot(G,prim.a(1)) < 0.5*M_PI) &&
		  (dot(G,prim.a(2)) < 0.5*M_PI)) {
		bool found = false;
		Vec3 twist = prim.k2Twist (G);
		Vec3 superTwist = Vec3 ((double)n0, (double)n1, (double)n2);
		for (int j=0; j<3; j++)
		  twist[j] -= round (twist[j]);
		for (int i=0; i<twistVecs.size(); i++) {
		  Vec3 diff = twist - twistVecs[i];
		  for (int j=0; j<3; j++)
		    diff[j] -= round (diff[j]);
		  found = found || (dot(diff,diff) < 1.0e-10);
		}
		if (!found) {
		  twistVecs.push_back(twist);
		  superTwists.push_back(S*twist);
		}
	      }
	    }
      }
  cerr << "Found " << twistVecs.size() << " distict k-vectors.\n";


  fprintf (stderr, "Reduced k-vectors of primitive cell:\n");
  for (int i=0; i<twistVecs.size(); i++) 
    fprintf (stderr, "  %9.5f %9.5f %9.5f\n", 
	     twistVecs[i][0], twistVecs[i][1], twistVecs[i][2]);

  fprintf (stderr, "\nReduced k-vectors of supercell:\n");
  for (int i=0; i<twistVecs.size(); i++) 
    fprintf (stderr, "  %9.5f %9.5f %9.5f\n", 
	     superTwists[i][0], superTwists[i][1], superTwists[i][2]);

  fprintf (stderr, "\n");
  fprintf (stdout, "K_POINTS {crystal}\n");
  fprintf (stdout, "%d\n", twistVecs.size());
  for (int i=0; i<twistVecs.size(); i++) 
    fprintf (stderr, "  %9.5f %9.5f %9.5f  1.0\n", 
	     twistVecs[i][0], twistVecs[i][1], twistVecs[i][2]);

}


int
main(int argc, char **argv)
{
  Mat3 Aprim, S;
  cerr << "Enter primitive lattice vectors:\n";
  for (int i=0; i<3; i++) {
    cerr << "  Vector " << i+1 << ":  ";
    cin >> Aprim(i,0) >> Aprim(i,1) >> Aprim(i,2);
  }
  cerr << "\n";
  
  TinyMatrix<int,3,3> Sint;
  cerr << "Enter tiling matrix:\n";
  for (int i=0; i<3; i++) {
    cerr << "  Row " << i << ":  ";
    cin >> Sint(i,0) >> Sint(i,1) >> Sint(i,2);
    S(i,0) = (double)Sint(i,0);
    S(i,1) = (double)Sint(i,1);
    S(i,2) = (double)Sint(i,2);
  }
  cerr << "\n";
  
  
//   Mat3 Aprim;
//   Aprim =  0.5, 0.5, 1.0,
//            0.5, 1.0, 0.5,
//            1.0, 0.5, 0.5;

//   Mat3 S;
//   S = -1.0, -1.0,  3.0,
//       -1.0,  3.0, -1.0,
//        3.0, -1.0, -1.0;

  Vec3 superShift (0.0, 0.0, 0.0);
  Int3 superMesh(1,1,1);
  cerr << "Energy desired supercell mesh (3 integrers).\n";
  cin >> superMesh[0] >> superMesh[1] >> superMesh[2];
  cerr << "\nEnter supercell shift in reduced coordinates.\n";
  cin >> superShift[0] >> superShift[1] >> superShift[2];
  vector<Vec3> twistVecs;
  MakekMesh (Aprim, S, superMesh, superShift,
	     twistVecs);

}
