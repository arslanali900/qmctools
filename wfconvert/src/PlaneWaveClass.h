#ifndef PLANE_WAVE_CLASS
#define PLANE_WAVE_CLASS

#include <Common/MPI/Communication.h>
#include <Common/IO/IO.h>
#include <Common/PlaneWavePHDFT/GVecs.h>
#include <Common/PlaneWavePHDFT/FFTBox.h>
#include "ParserClass.h"
#include "LatticeClass.h"

using namespace blitz;

class BitField3DClass
{
  int Nx, Ny, Nz;
  Array<unsigned int,1> Field;
public:

  inline bool operator()(int ix, int iy, int iz) {
    int num = iz + Nz*(iy + Ny*ix);
    int wordnum = num >> 5;
    int bitnum  = num & 0x001f;
    unsigned int mask = 1 << bitnum;
    return (Field(wordnum) & mask) ? true : false;
  }

  inline double GetDouble(int ix, int iy, int iz) {
    int num = iz + Nz*(iy + Ny*ix);
    int wordnum = num >> 5;
    int bitnum  = num & 0x001f;
    unsigned int mask = 1 << bitnum;
    return (Field(wordnum) & mask) ? 1.0 : 0.0;
  }
  
  inline void Set (int ix, int iy, int iz, bool val)
  {
    int num = iz + Nz*(iy + Ny*ix);
    int wordnum = num >> 5;
    int bitnum  = num & 0x001f;
    unsigned int mask = 1 << bitnum;
    if (val)
      Field(wordnum) |= mask;
    else 
      Field(wordnum) &= (mask ^ 0xffffffff);
  }
  inline void Init (int nx, int ny, int nz) 
  {
    Nx = nx; Ny = ny; Nz = nz;
    Field.resize((nx*ny*nz+31)/32);
    Field = 0;
  }
};
    

class CenterClass
{
public:
  Vec3 r;
  double Radius, SkinThickness;
  bool Spherical;
  int NumOrbitals;
  int NumUp, NumDown;
  BitField3DClass Bitfield;
  void Broadcast (CommunicatorClass &comm, int root=0);
  CenterClass (Vec3 center, double radius, 
	       int numOrb, bool spherical=true)
  {
    r = center;
    Radius = radius;
    Spherical = spherical;
    NumOrbitals = numOrb;
  }
  CenterClass() : 
    Spherical(true), r(0.0, 0.0, 0.0), Radius(0.0), NumOrbitals(1),
    SkinThickness(0.0)
  {

  }
};


class OrbitalClass;

class CellClass
{
private:
  bool IsSuper;
  bool FFTisSetup;
public:
  LatticeClass Lattice;
  GVecsClass GVecs;
  FFTBox FFT;
  Array<Vec3, 1> IonPos;
  Array<int,  1> AtomTypes;

  /////////////////////
  // Member functions//
  /////////////////////
  void SetLattice (Mat3 A, bool isSuper);
  void PutOrbital (OrbitalClass &orbital);
  void AddOrbital (OrbitalClass &orbital,
		   complex<double> prefactor);
  inline void SetupFFT()
  {  
    if (!FFTisSetup) {
      FFTisSetup = true;
      FFT.Setup();
    }
  }  

  void Broadcast (CommunicatorClass &comm, int root);

  // Constructor
  CellClass() : FFT(GVecs), FFTisSetup(false)
  {
  }
};



// This class stores the data necessary for each primitive k-point
// which is read in.  If we apply a tiling matrix, the k-point mesh
// will be folded into the bands of the superlattice k-point mesh.
class kPointClass
{
private:
  CellClass *PrimCell, *SuperCell;
  Array<Int3,1> PrimIndices, SuperIndices;
  // The twist vector w.r.t. the primitive lattice
  Vec3 PrimTwist;
  // The twist vector w.r.t. the super lattice
  // This is devided into an integer part and a factional part.
  Vec3 SuperTwistInt, SuperTwistFrac;
public:
  void SetCells (CellClass &primCell, CellClass &superCell,
		 Vec3 k);
  Int3 GetFFTBoxSize (Array<Vec3,1> &primGVecs, double fftFactor);
  void SetIndices (Array<Vec3,1> &primGvecs, Int3 fftBoxSize);
  void PutzVec (zVec &c, FFTBox &fft, bool useSuper);
  void AddzVec (zVec &c, complex<double> prefactor, 
		FFTBox &fft, bool useSuper);
  inline Vec3 GetPrimk()          { return PrimCell->Lattice.Twist2k (PrimTwist); }
  inline Vec3 GetPrimTwist()      { return PrimTwist;                             }
  inline Vec3 GetSuperTwistInt()  { return SuperTwistInt;                         }
  inline Vec3 GetSuperTwistFrac() { return SuperTwistFrac;                        }
  void Broadcast (CommunicatorClass &comm, int root);
};





class OrbitalClass
{
private:
  CellClass &PrimCell, &SuperCell;
  zVec Coefs;
  double Eigenvalue;
  // This stores the fraction of the norm^2 contained in the
  // localization radius if Center!=NULL
  double Norm2;
  int Spin;
  kPointClass &kPoint;
  int PrimTwistIndex, PrimBandIndex;
public:
  CenterClass *Center;

  bool Read (ParserClass &parser);
  inline kPointClass &GetkPoint()             { return kPoint; }
  inline Vec3 GetPrimk()                      { return kPoint.GetPrimk(); }
  inline Vec3 GetPrimTwist()                  { return kPoint.GetPrimTwist(); }
  inline Vec3 GetSuperTwistFrac()             { return kPoint.GetSuperTwistFrac(); }
  inline Vec3 GetSuperTwistInt()              { return kPoint.GetSuperTwistFrac(); }
  inline void SetEigVal (double eig)          { Eigenvalue = eig;  }
  inline double GetEigVal ()                  { return Eigenvalue; }
  inline void SetSpin (int spin)              { Spin = spin;       }
  inline int  GetSpin ()                      { return Spin;       }
  inline zVec& GetCoefs ()                    { return Coefs;      }
  inline void SetCenter (CenterClass &center) { Center = &center;  }
  inline void SetNorm2 (double norm2)         { Norm2 = norm2;     }
  inline double GetNorm2 ()                   { return Norm2;      }
  inline const CellClass& GetPrimCell()       { return PrimCell;   }
  inline const CellClass& GetSuperCell()      { return SuperCell;  }
  inline void SetPrimIndices (int twistIndex, int bandIndex)
  { PrimTwistIndex = twistIndex;    PrimBandIndex = bandIndex;     }
  inline int GetPrimTwistIndex()              { return PrimTwistIndex; }
  inline int GetPrimBandIndex()               { return PrimBandIndex;  }
  void Read  (IO::IOSectionClass &in);
  void Write (IO::IOSectionClass &out);
  void WriteSpline (IO::IOSectionClass &out, bool real, bool shift,
		    bool truncate, bool writeSuper);
  void WriteBWFN (FILE *fout, bool writeComplex, bool writeSuper);
  void WritePWFN (FILE *fout);
  void CheckKineticEnergy();
  
  OrbitalClass(CellClass &prim, CellClass &super, kPointClass &k) : 
    PrimCell (prim), SuperCell(super), kPoint(k), 
    Center(NULL), Spin(0), PrimTwistIndex (-1), PrimBandIndex(-1)
  { 
    // nothing for now.
  }
};



class TwistClass
{
private:
  CellClass *Cell;
  Vec3 Twist;
public:
  // The 2 corresponds to the up and down bands
  vector<OrbitalClass*> Bands[2];
  Vec3 GetTwist();
  Vec3 Getk();
  Vec3 GetTwistFrac();
  Vec3 GetTwistInt();
  void Set(CellClass &cell, Vec3 k);
};


class PlaneWaveSystem 
{
private:
  bool ReadRunInfo  (ParserClass &parser);
  bool ReadIonPos   (ParserClass &parser);
  void TileIonPos();
  bool ReadLattice  (ParserClass &parser);
  bool ReadGVecs    (ParserClass &parser);
  bool ReadOrbitals (ParserClass &parser);
  bool ReadCenters  (string fname);
  inline Vec3 MinImage (Vec3 dr);
  void SetBitfield (CenterClass &center);
  void LocalFunc (CenterClass center, Array<double,3> &func);
  void ThetaMatrix  (CenterClass &center,
		     Array<complex<double>,2> &theta);

  /////////////////////////////////////////////
  // k-point mesh related data and functions //
  /////////////////////////////////////////////
  // Stores the k-point mesh dimensions
  Int3 kPointMesh;
  void SetupSuperOrbitals();
  // Array of structures holding information necessary to "unfold" the
  // k-point mesh 
  Array<kPointClass,1> PrimkPoints;

  ///////////////////////
  // Cell-related data //
  ///////////////////////
  CellClass PrimCell, SuperCell;
  bool WriteSuper;


  /////////////////////////////////
  // Localization data and funcs //
  /////////////////////////////////
  double SkinThickness;
  CommunicatorClass Comm;
  // These store the number of up and down bands to use in localizing
  // creating localized orbitals.
  int UseUpBands, UseDownBands;
  // This stores the rotation matrix necessary to construct the
  // localized orbitals from the extended ones.  The indices are 
  // [spin](kpoint, newBand, oldBand)
  TinyVector<Array<complex<double>,3>,2> NewOrbCoefs;
  // The indices are [spin](kpoint, newBand)
  TinyVector<Array<double,2>,2> NewOrbNorms, NewEigVals;
  void LocalizeMPI (bool ortho, int spin);
  void DistributeBands();
  /// Returns which bands a processor is reponsible for
  void ProcBands(int proc, int numOcc, int &first, int &last);
  void ProcTasks (int proc, int numTasks, int &first, int &last);
  void Localize (bool ortho, int spin);
  // Returns true if converged.
  bool UpdateCenters(int spin);
  // Uses a polar decomposition to return the nearest set of
  // orthogonal orbitals to the present, nonorthogonal set.
  // void OrthogonalizeOrbitals(int spin);
  vector<CenterClass> Centers;
  /// Returns the nearest ion to pos in PBC
  int NearestIon (Vec3 pos);

  // CASINO export
  void WriteCasinoHeader (FILE *fout);
  void WriteCentres();

  // Original run data
  string GeneratedBy, Method, Pseudopotential, Functional;
  double TotalEnergy, KineticEnergy, LocalPotEnergy, NonlocalEnergy,
    eeEnergy, IonIonEnergy;

  // The matrix below transforms the primitive lattice vectors into
  // the superlattice vectors
  Mat3 TileMatrix;

  // This is the FFT box resolution enhancement factor.
  double FFTFactor;
  
  // Creates the a real-space orbital in FFT.rBox.  This function
  // applies the rotation in NewOrbCoefs as it is constructing the
  // real-space orbital.
  void MakeRealSpaceOrbital (int spin, int ki, int band);
  void WriteSpline (IO::IOSectionClass &out, int spin, int ki, int band);
  void MakeRealSpaceBlip    (int spin, int ki, int band);

public:
  Array<Vec3,1> GVecsArray;
  Array<Int3,1> GInts;
  // The outer tinyvector is for up and down
  // The first index is the kPoint index
  // The second index is the band index
  TinyVector<Array<OrbitalClass*,2>,2> Orbitals, SuperOrbitals;
  double ECut;
  int NumElectrons;
  bool Localized, OptimizeCenters, OptimizeRadii, 
    ShiftOrbitals, CheckKE, Spline, Real, Orthogonalize,
    Truncate, SpinPolarized;

  bool Read (string fname);
  bool Read (IO::IOSectionClass &in);
  void ReadABINIT_WFK (string fname);

  void Localize (string centersFile, double radius=3.0);
  /// foldFactor gives the unfolding factor.  E.g. if the k-point mesh
  /// was 4x4x4, and foldFactors was (2,2,2), we would unfolding into
  /// 2x2x2 supercells, each with a 2x2x2 k-point mesh.  This can
  /// allow combining localization with twist-averaging.
  /// void Unfold(int spin, Int3 foldFactor=Int3(0,0,0));

  /// Tile does basically the same thing, but does not actually unfold
  /// the data.  Rather, it relabels the bands and k-points and
  /// breaks each twist vector into two pieces:  a part which is a
  /// G-vector of the supercell and a remainder.
//   void Tile  (int spin, Int3 tileFactor);
//   void Tile  (int spin, Mat3 tileMatrix);
  void SetTileMatrix (Mat3 tileMatrix);
  bool Write (string fname);
  bool WritePWFN (string fname);
  bool WriteBWFN (string fname);
  bool WriteBWFN (string fname, Int3 split);
  void SetSkinThickness (double thickness);
  inline void SetFFTFactor (double factor) { FFTFactor = factor; }
  PlaneWaveSystem () : Localized(false), 
		       OptimizeCenters(false), OptimizeRadii(false), 
		       ShiftOrbitals(false), CheckKE(false), Spline(false), 
		       Real(false), Orthogonalize(false), Truncate(false),
		       FFTFactor(1.0), WriteSuper(false)
  {
    Comm.SetWorld();
    TileMatrix(0,0)=1.0; TileMatrix(0,1)=0.0; TileMatrix(0,2)=0.0;
    TileMatrix(1,0)=0.0; TileMatrix(1,1)=1.0; TileMatrix(1,2)=0.0;
    TileMatrix(2,0)=0.0; TileMatrix(2,1)=0.0; TileMatrix(2,2)=1.0;
    // do nothing
  }
};


inline Vec3
PlaneWaveSystem::MinImage (Vec3 r)
{
  CellClass &cell = WriteSuper ? SuperCell : PrimCell;

  Vec3 u, rp;
  Vec3 minDisp = r;
  double minDist = dot(r,r);
  
  Vec3 a0 = cell.Lattice.a(0);
  Vec3 a1 = cell.Lattice.a(1);
  Vec3 a2 = cell.Lattice.a(2);

  for (int ix=-1; ix<=1; ix++) 
    for (int iy=-1; iy<=1; iy++) 
      for (int iz=-1; iz<=1; iz++) {
	Vec3 tryDisp = r + (double)ix*a0 + (double)iy*a1 + (double)iz*a2;
	double tryDist = dot (tryDisp, tryDisp);
	if (tryDist < minDist) {
	  minDist = tryDist;
	  minDisp = tryDisp;
	}
      }
  return minDisp;

  u = cell.Lattice.r2u(r);
  u[0] -= round(u[0]);
  u[1] -= round(u[1]);
  u[2] -= round(u[2]);

  rp = cell.Lattice.u2r(u);
  return rp;
}

#endif
