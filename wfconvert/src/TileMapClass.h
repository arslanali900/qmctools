#ifndef TILE_MAP_H
#define TILE_MAP_H

#include "LatticeClass.h"
#include <map>

typedef TinyVector<int,2> OrbIndex;
class CommunicatorClass;

struct OrbLess
{
  inline bool operator()(const OrbIndex &a, const OrbIndex &b)
  {
    if (a[0]<b[0])
      return true;
    if (a[0]>b[0])
      return false;
    if (a[1] >= b[1])
      return false;  
    return true;
  }
};

class OrbitalClass;
class kPointClass;

class TileMapClass
{
private:
  LatticeClass PrimLattice, SuperLattice;
  map<OrbIndex, OrbIndex, OrbLess> Prim2SuperMap, Super2PrimMap;
  int NumPrimTwists, NumPrimBands, NumSuperTwists, NumSuperBands;
public:
  void CreateMap (Array<OrbitalClass*,2> &primOrbs,
		  LatticeClass &primLattice,
		  LatticeClass &superLattice);
  int GetNumPrimTwists();
  int GetNumPrimBands();
  int GetNumSuperTwists();
  int GetNumSuperBands();
  void Prim2Super (int primTwist,   int primBand,
		   int &superTwist, int &superBand);
  void Super2Prim (int primTwist,   int primBand,
		   int &superTwist, int &superBand);
  void Broadcast (CommunicatorClass &comm, int root);
};


#endif

