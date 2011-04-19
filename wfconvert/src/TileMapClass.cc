#include "TileMapClass.h"
#include "OrbitalClass.h"
#include <Common/MPI/Communication.h>
#include <vector>
#include <algorithm>

using namespace std;

struct OrbInfo
{
  double E;
  int ki, bi;

  OrbInfo (int ki_, int bi_, double E_) :
    E(E_), ki(ki_), bi(bi_)
  { }
};

struct Eless
{
  inline bool operator() (OrbInfo orb1, OrbInfo orb2)
  { return orb1.E < orb2.E; }
};

int
TileMapClass::GetNumPrimTwists()
{ return NumPrimTwists; }

int
TileMapClass::GetNumPrimBands()
{ return NumPrimBands; }

int
TileMapClass::GetNumSuperTwists()
{ return NumSuperTwists; }

int
TileMapClass::GetNumSuperBands()
{ return NumSuperBands; }

void
TileMapClass::Prim2Super (int primTwist, int primBand,
			  int &superTwist, int &superBand)
{
  OrbIndex prim(primTwist, primBand);
  if (Prim2SuperMap.find(prim) == Prim2SuperMap.end()) {
    cerr << "(twist, band) = (" << primTwist << ", " << primBand 
	 << ") not found in TwistMapClass::Prim2Super.\n";
    abort();
  }
  OrbIndex super = Prim2SuperMap[prim];
  superTwist = super[0];
  superBand  = super[1];
}
void
TileMapClass::Super2Prim (int superTwist, int superBand,
			  int &primTwist, int &primBand)
{
  OrbIndex super(superTwist, superBand);
  if (Super2PrimMap.find(super) == Super2PrimMap.end()) {
    cerr << "(twist, band) = (" << primTwist << ", " << primBand 
	 << ") not found in TwistMapClass::Prim2Super.\n";
    abort();
  }
  OrbIndex prim = Super2PrimMap[super];
  primTwist = prim[0];
  primBand  = prim[1];
}

inline Vec3 
IntPart (Vec3 twist)
{
  Vec3 intpart(round(twist[0]), round(twist[1]), round(twist[2]));
  Vec3 fracpart = twist - intpart;
  for (int i=0; i<3; i++) 
    if (fabs(fracpart[i] + 0.5) < 1.0e-5) {
      fracpart[i] += 1.0;
      intpart[i]  -= 1.0;
    }
  return intpart;
}

inline Vec3
FracPart (Vec3 twist)
{
  return twist - IntPart(twist);
}


void
TileMapClass::CreateMap(Array<OrbitalClass*,2> &primOrbs,
			LatticeClass &primLattice,
			LatticeClass &superLattice)
{
  NumPrimTwists = primOrbs.extent(0);
  NumPrimBands  = primOrbs.extent(1);

  // First, divide the primitive k-vectors into sets of supercell
  // k-vectors. 

  // A vector holding the distinct fractional parts of the supercell
  // twist vectors
  vector<Vec3> superFracs;
  // This holds to which supercell kpoint each primitive k-point belongs
  Array<int,1> superIndex(NumPrimTwists);
  for (int ki=0; ki < NumPrimTwists; ki++) {
    Vec3 k = primOrbs(ki,0)->Getk();
    Vec3 superTwist = superLattice.k2Twist(k);
    Vec3 frac = FracPart (superTwist);
    bool found = false;
    for (int i=0; i<superFracs.size(); i++) {
      Vec3 diff = superFracs[i] - frac;
      if ((dot (diff, diff) < 1.0e-10)) {
	found = true;
	superIndex(ki) = i;
      }
    }
    if (!found) {
      superIndex(ki) = superFracs.size();
      superFracs.push_back(frac);
    }
  }

  int numTwistsNeeded = (int)round (superLattice.GetVolume()/
				    primLattice.GetVolume());

  NumPrimBands   = primOrbs.extent(1);
  NumPrimTwists  = primOrbs.extent(0);
  NumSuperTwists = superFracs.size();
  NumSuperBands  = numTwistsNeeded * NumPrimBands;


  // For each supercell twist, create a list of primitive twists which
  // belong to it.
  Array<vector<int>,1> superSets(NumSuperTwists);
  for (int ki=0; ki < NumPrimTwists; ki++) 
    superSets(superIndex(ki)).push_back(ki);
  
  for (int si=0; si<NumSuperTwists; si++) {
    fprintf (stderr, "Super twist #%d:  [ %9.5f %9.5f %9.5f ]\n",
	     si, superFracs[si][0], superFracs[si][1], superFracs[si][2]);
    fprintf (stderr, "  Using k-points: ");
    for (int i=0; i<superSets(si).size(); i++) 
      fprintf (stderr, " %d", superSets(si)[i]);
    fprintf (stderr, "\n");
  }
  
  // Now check to see that each supercell twist has the right twists
  // to tile the primitive cell orbitals.
  for (int si=0; si<NumSuperTwists; si++) {
    // First make sure we have enough points
    if (superSets(si).size() != numTwistsNeeded) {
      fprintf (stderr, "Super twist %d should own %d k-points, but owns %d.\n",
	       si, numTwistsNeeded, superSets(si).size());
      abort();
    }
    // Now, make sure they are all distinct
    int N = superSets(si).size();
    for (int i=0; i<N; i++) {
      Vec3 ik = primOrbs(superSets(si)[i],0)->Getk();
      Vec3 iInt = IntPart(superLattice.k2Twist(ik));
      // Vec3 iInt = primkPoints(superSets(si)[i]).GetSuperTwistInt();
      for (int j=i+1; j<N; j++) {
	Vec3 jk = primOrbs(superSets(si)[j],0)->Getk();
	Vec3 jInt = IntPart(superLattice.k2Twist(jk));
	// Vec3 jInt = primkPoints(superSets(si)[j]).GetSuperTwistInt();
	if (dot(iInt-jInt, iInt-jInt) < 1.0e-6) {
	  cerr << "Identical k-points detected in super twist set "
	       << si << endl;
	  abort();
	}
      }
    }
  }
  ///////////////////////////////////////
  // Construct super twist orbital map //
  ///////////////////////////////////////
  NumPrimBands  = primOrbs.extent(1);
  NumPrimTwists = primOrbs.extent(0);
  cerr << "There are " << superSets.size() << " supersets.\n";
  for (int si=0; si < NumSuperTwists; si++) {
    // Create a vector of all the orbitals belonging to this
    // supertwist.  We will then sort them by energy
    vector<OrbInfo> myOrbitals;
    for (int i=0; i<superSets(si).size(); i++) {
      int ki = superSets(si)[i];
      for (int bi=0; bi<NumPrimBands; bi++) {
	OrbitalClass &orb = *primOrbs(ki,bi);
	double E = orb.GetEigVal();
	myOrbitals.push_back(OrbInfo(ki,bi,E));
      }
    }
    Eless comparison;
    sort (myOrbitals.begin(), myOrbitals.end(), comparison);
    for (int i=0; i<myOrbitals.size(); i++) {
      OrbIndex superIndex (si,i);
      OrbIndex primIndex (myOrbitals[i].ki, myOrbitals[i].bi);
      // fprintf (stderr, "i = %d Super index = (%d,%d)  Prim index = (%d,%d)\n",
      // 	       i, superIndex[0], superIndex[1], primIndex[0], primIndex[1]);
      map<OrbIndex, OrbIndex, OrbLess>::value_type s2p(superIndex, primIndex);
      map<OrbIndex, OrbIndex, OrbLess>::value_type p2s(primIndex, superIndex);
      Prim2SuperMap.insert(p2s);
      Super2PrimMap.insert(s2p);
      //Super2PrimMap[superIndex] =  primIndex;
      //Prim2SuperMap[ primIndex] = superIndex;
    }
  }
}


void
TileMapClass::Broadcast (CommunicatorClass &comm, int root)
{
  Mat3 Aprim  =  PrimLattice.GetDirect();
  Mat3 Asuper = SuperLattice.GetDirect();
  comm.Broadcast (root, Aprim);
  comm.Broadcast (root, Asuper);
  PrimLattice.SetDirect(Aprim);
  SuperLattice.SetDirect(Asuper);

  comm.Broadcast (root, NumPrimTwists);
  comm.Broadcast (root, NumPrimBands);
  comm.Broadcast (root, NumSuperTwists);
  comm.Broadcast (root, NumSuperBands);
  
  // Now, broadcast maps
  int mapSize = Prim2SuperMap.size();
  comm.Broadcast (root, mapSize);


  Array<OrbIndex,1> prim(mapSize), super(mapSize);
  map<OrbIndex, OrbIndex, OrbLess>::iterator iter;
  int i=0;
  for (iter=Prim2SuperMap.begin(); iter != Prim2SuperMap.end(); iter++) {
    prim(i)  = iter->first;
    super(i) = iter->second;
    i++;
  }
  comm.Broadcast(root, prim);
  comm.Broadcast(root, super);

  if (comm.MyProc() != root) 
    for (int i=0; i<mapSize; i++) {
      Prim2SuperMap[ prim(i)] = super(i);
      Super2PrimMap[super(i)] =  prim(i);
    }
}
