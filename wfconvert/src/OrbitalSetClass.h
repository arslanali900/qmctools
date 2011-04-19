#ifndef ORBITAL_SET_CLASS_H
#define ORBITAL_SET_CLASS_H

#include "OrbitalClass.h"
#include "TileMapClass.h"
#include "CellClass.h"
#include "CenterClass.h"
#include "APWClass.h"
#include "AtomicOrbital.h"
#include <Common/MPI/Communication.h>
#include <vector>


class ParserClass;

class OrbitalSetClass 
{
private:
  ////////////////
  // Basic data //
  ////////////////
  CellClass PrimCell, SuperCell, FineCell;
  Array<Vec3,1> SuperGVecs;
  // The outer tinyvector is for up and down
  // The first index is the kPoint index and the second the band index
  TinyVector<Array<OrbitalClass*,2>,2> PrimOrbitals, SuperOrbitals, FineOrbitals;
  vector<TinyVector<Array<OrbitalClass*,2>,2> > PrimFirstOrder, SuperFirstOrder;
  DensityClass PrimDensity, SuperDensity;
  DensityClass PrimVHXC, SuperVHXC;
  bool WriteSuper;
  map<int, string> ZToSymbolMap;
  map<string, int> SymbolToZMap;
  map<int, double> ZToMassMap;
  map<string,double> UnitToHartreeMap;
  map<string,double> UnitToBohrMap;
  void SetupMaps();

  //////////////////////////////////////////////////
  // Functions for reading CASINO pwfn.data files //
  //////////////////////////////////////////////////
  bool ReadRunInfo  (ParserClass &parser);
  bool ReadIonPos   (ParserClass &parser);
  bool ReadLattice  (ParserClass &parser);
  bool ReadGVecs    (ParserClass &parser);
  bool ReadOrbitals (ParserClass &parser);
  bool ReadCenters  (string fname);
  bool ReadCenters  (IO::IOSectionClass &in);

  ////////////////////////////////////////////
  // Data for atomic orbital representation //
  ////////////////////////////////////////////
  vector<AtomicOrbital> AtomicOrbitals;
  
  //////////////////////////////////////////////////////
  // Functions and data for tiling the primitive cell //
  //////////////////////////////////////////////////////

  // The matrix below transforms the primitive lattice vectors into
  // the superlattice vectors
  Mat3 TileMatrix;
  TinyVector<TileMapClass,2> TileMap;
  void TileIonPos();
  // These arrays, one for each primitive k-point, store the indices
  // necessary to put an array of complex coefficients into the
  // SuperCell's FFT box.
  Array<Int3,1> PrimFFTIndices, FineFFTIndices;
  vector<Array<Int3,1> > SuperFFTIndices;
  vector<Vec3> SuperGshifts;
  void SetupFFTIndices();
  void CreateSuperOrbitals();
  void CreateFineOrbitals();

  /////////////////////////////////
  // Localization data and funcs //
  /////////////////////////////////
  CommunicatorClass Comm;
  // Stores which bands this processor is responsible for storing.
  // The two components are for up and down spin.
  int MyFirstBand[2], MyLastBand[2];

  double SkinThickness;
  // These store the number of up and down bands to use in localizing
  // creating localized orbitals.
  int UseUpBands, UseDownBands, NumExtendedBands;
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

  ///////////////////
  // CASINO export //
  ///////////////////
  void WriteCasinoHeader (FILE *fout);
  void WriteCentres();

  ///////////////////////
  // Original run data //
  ///////////////////////
  string GeneratedBy, Method, Pseudopotential, Functional;
  double TotalEnergy, KineticEnergy, LocalPotEnergy, NonlocalEnergy,
    eeEnergy, IonIonEnergy;

  // This is the FFT box resolution enhancement factor.
  double FFTFactor;

  //////////////////////
  // LAPW information //
  //////////////////////
  APWClass APW;
  LocalOrbitalClass LocalOrbitals;
  vector<CoreStateClass> CoreStates;
  const int APWRadialPoints, CoreRadialPoints;
  
  // Creates the a real-space orbital in FFT.rBox.  This function
  // applies the rotation in NewOrbCoefs as it is constructing the
  // real-space orbital.  Returns the scalar rotation factor.
  complex<double> MakeRealSpaceOrbital (int spin, int ki, int band, bool shiftOrbs);
  void MakeFineOrbital   (int spin, int ki, int band, bool shiftOrbs, complex<double> rotation);
  void MakeFineLaplacian (int spin, int ki, int band, bool shiftOrbs, complex<double> rotation);
  void MakeFirstOrderOrbital (int ideriv, int spin, int ki, int band, 
			      bool shiftOrbs, complex<double> rotation);
  void WriteSpline  (IO::IOSectionClass &out, int spin, int ki, int band);
  void WriteESHDFSpline  (IO::IOSectionClass &out, int spin, int ki, int band,
			 bool writeComplex);
  void WriteESHDFMultiRep (IO::IOSectionClass &out, int spin, int ki, int band,
			   bool writeComplex);

  void WriteBWFNOrb (FILE *fout, int spin, int ki, int band, bool writeComplex);
  void MakeRealSpaceBlip    (int spin, int ki, int band);
  void MakeRealSpaceInterp  (int spin, int ki, int band);
  void MakeBoxShifter (Vec3 rshift, Array<complex<double>,3> &shiftBox);
public:
  Array<Vec3,1> GVecsArray;
  Array<Int3,1> GInts;
  // The outer tinyvector is for up and down
  // The first index is the kPoint index
  // The second index is the band index
  double ECut;
  TinyVector<int,2> NumElectrons;
  bool Localized, OptimizeCenters, OptimizeRadii, 
    ShiftOrbitals, CheckKE, Spline, Real, Orthogonalize,
    Truncate, SpinPolarized, UseMultiRep, CompareHybrid,
    WriteInterp;


  bool Read (string fname);
  bool ReadPWFN (string fname);
  bool Read (IO::IOSectionClass &in);
  bool Read_ABINIT_WFK (string fname);
  bool Read_ABINIT_WF1 (string fname);
  bool Read_ABINIT_First_Order(string prefix);
  bool Read_ABINIT_First_Order_FD(string name);
  bool Read_ABINIT_DEN (string fname);
  bool Read_ABINIT_POT (string fname);
  bool Read_FPMD  (string fname);
  bool Read_LAPW  (string fname);
  bool Read_ESHDF (string fname);
  bool ReadMultiRep (string fname);

  void ComputeDensity();
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
  bool Write_qmcPACK (string fname);
  bool Write_ESHDF (string fname, int lappwr);
  bool WritePWFN (string fname);
  bool WriteBWFN (string fname, bool splitTwists);
  void SetSkinThickness (double thickness);
  inline void SetFFTFactor (double factor) { FFTFactor = factor; }
  OrbitalSetClass () : Localized(false), 
		       OptimizeCenters(false), OptimizeRadii(false), 
		       ShiftOrbitals(false), CheckKE(false), Spline(false), 
		       Real(false), Orthogonalize(false), Truncate(false),
		       FFTFactor(1.0), WriteSuper(true),
		       NumExtendedBands(0), APWRadialPoints(1000),
		       CoreRadialPoints(25000), CompareHybrid(false), 
		       WriteInterp(false)
  {
    Comm.SetWorld();
    TileMatrix(0,0)=1.0; TileMatrix(0,1)=0.0; TileMatrix(0,2)=0.0;
    TileMatrix(1,0)=0.0; TileMatrix(1,1)=1.0; TileMatrix(1,2)=0.0;
    TileMatrix(2,0)=0.0; TileMatrix(2,1)=0.0; TileMatrix(2,2)=1.0;
    SetupMaps();
  }
};


#endif
