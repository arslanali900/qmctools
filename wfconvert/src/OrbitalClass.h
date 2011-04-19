#ifndef ORBITAL_CLASS_H
#define ORBITAL_CLASS_H

#include <Common/Blitz.h>
#include <vector>
#include "CellClass.h"
typedef TinyVector<int,3> Int3;
typedef blitz::Array<complex<double>,1> zVec;

namespace IO {
  class IOSectionClass;
};

class CellClass;
class CenterClass;
class ParserClass;
class CommunicatorClass;

class AtomicDataClass
{
public:
  // First index is radial grid point
  // Second index is lm = l*(l+1) + m
  Array<complex<double>,2> SplineData;
  // First index is n in r^n
  // Second index is lm = l*(l+1) + m
  Array<complex<double>,2> PolyCoefs;
};

class OrbitalClass
{
private:
  CellClass *Cell;
  zVec* Coefs;
  Array<Int3,1> *FFTIndices;
  Vec3 Twist, Gshift;
  double Eigenvalue, Occupancy, KE;
  // This stores the fraction of the norm^2 contained in the
  // localization radius if Center!=NULL
  double Norm2;
  int Spin, TwistIndex, BandIndex;
  // int PrimTwistIndex;
  CenterClass *Center;
  int DerivAtom, DerivDir;
public:
  vector<AtomicDataClass> AtomicData;

  inline void SetCenter(CenterClass &center);
  inline CenterClass &GetCenter();
  inline void SetTwist (Vec3 twist);
  inline Vec3 GetTwist ();
  inline void SetGshift (Vec3 gshift);
  inline Vec3 GetGshift ();
  inline void   SetEigVal(double val);
  inline double GetEigVal();
  inline void   SetOccupancy(double val);
  inline double GetOccupancy();
  inline void SetSpin (int spin);
  inline int  GetSpin ();
  inline void SetLabels (int  spin, int  ik, int  iband);
  inline void GetLabels (int &spin, int &ik, int &iband);
  inline zVec& GetCoefs();
  inline void SetCoefs (zVec& coefs);
  inline void SetFFTIndices (Array<Int3,1> &fftIndices);
  inline Array<Int3,1>& GetFFTIndices ();
  inline CellClass& GetCell();
  inline void SetCell (CellClass &cell);
  inline void SetDeriv (int atom, int dir)   { DerivAtom=atom; DerivDir=dir;   }
  inline void GetDeriv (int &atom, int &dir) { atom=DerivAtom; dir = DerivDir; }

  void Setk (Vec3 k);
  Vec3 Getk ();
  void Read (IO::IOSectionClass &in);
  bool Read (ParserClass &parser, int numCoefs);
  void Write (IO::IOSectionClass &out);
  void WriteSpline (IO::IOSectionClass &out, bool real,
		    bool shift, bool truncate);
  void WriteBWFN (FILE *fout, bool real, bool shift,
			 bool truncate);
  void PutInFFTBox();
  void PutInFFTBox(zVec &coefs);
  void AddToFFTBox(zVec &coefs);
  void AddToFFTBox(complex<double> prefactor);
  void Broadcast (CommunicatorClass &comm, int root);
  void Send      (CommunicatorClass &comm, int toProc);
  void Receive   (CommunicatorClass &comm, int fromProc);

  OrbitalClass() :
    Cell(NULL), Coefs(NULL), FFTIndices(NULL), Twist(0.0, 0.0, 0.0),
    Eigenvalue(0.0), Norm2(0.0), Spin(0), Center(NULL),
    Gshift(0.0, 0.0, 0.0)
  { }
};

class DensityClass
{
private:
  CellClass Cell;
  zVec UpCoefs, DownCoefs;
  Array<Int3,1> GVecIndices, FFTIndices;
  Array<double,3> RhoUp, RhoDown;
  
public:
  inline void SetLattice (Mat3 lattice);
  inline const Array<Int3,1>& GetGVecIndices() { return GVecIndices; }
  void Write       (IO::IOSectionClass &out);
  void WriteSpline (IO::IOSectionClass &out);
  void WriteESHDF  (IO::IOSectionClass &out, string basename="density");
  void Set (const Array<double,3> &rhoTotal, double ecut);
  void Set (const Array<double,3> &rhoUp, 
	    const Array<double,3> &rhoDown, double ecut);
  void Set(const Array<Int3,1> &gvecs, const zVec &coefs);
  void Set(const Array<Int3,1> &gvecs,  
	   const zVec &upCoefs, const zVec &downCoefs);
  inline bool IsSet() { return RhoUp.size() > 0; }
  
};


inline void 
DensityClass::SetLattice(Mat3 lattice)
{
  Cell.SetLattice(lattice);
}


inline void
OrbitalClass::SetCenter(CenterClass &center)
{
  Center = &center;
}
  
inline CenterClass&
OrbitalClass::GetCenter()
{ 
  return *Center;
}

inline void 
OrbitalClass::SetTwist (Vec3 twist)
{ 
  Twist = twist;
}

inline Vec3 
OrbitalClass::GetTwist ()
{ 
  return Twist;
}

inline void
OrbitalClass::SetGshift (Vec3 gshift)
{
  Gshift = gshift;
}

inline Vec3
OrbitalClass::GetGshift ()
{
  return Gshift;
}

inline void 
OrbitalClass::SetEigVal(double val)
{ 
  Eigenvalue = val;
}

inline double 
OrbitalClass::GetEigVal()
{ 
  return Eigenvalue;
}

inline void 
OrbitalClass::SetOccupancy(double val)
{ 
  Occupancy = val;
}

inline double 
OrbitalClass::GetOccupancy()
{ 
  return Occupancy;
}

inline void 
OrbitalClass::SetSpin (int spin)
{ 
  Spin = spin;
}

inline int  
OrbitalClass::GetSpin ()
{ 
  return Spin;
}


inline void
OrbitalClass::SetLabels (int spin, int ik, int iband)
{
  Spin       = spin;
  TwistIndex = ik;
  BandIndex  = iband;
}

inline void
OrbitalClass::GetLabels (int &spin, int &ik, int &iband)
{
  spin  = Spin;
  ik    = TwistIndex;
  iband = BandIndex;
}

inline void 
OrbitalClass::SetCoefs (zVec& coefs)
{ 
  Coefs = &coefs;
}

inline zVec& 
OrbitalClass::GetCoefs()
{ 
  return *Coefs;
}

inline void 
OrbitalClass::SetFFTIndices (Array<Int3,1> &fftIndices)
{ 
  FFTIndices = &fftIndices;
}

inline Array<Int3,1>& 
OrbitalClass::GetFFTIndices ()
{ 
  return *FFTIndices;
}

inline void 
OrbitalClass::SetCell (CellClass &cell)
{
  Cell = &cell;
}

inline CellClass& 
OrbitalClass::GetCell()
{ 
  return *Cell;
}


#endif
