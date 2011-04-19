#include "ParseCommand.h"
#include "LatticeClass.h"
#include <Common/IO/IO.h>
#include <einspline/bspline.h>
#include <string>

using namespace IO;
using namespace std;

class OrbitalSplineClass
{
private:
  LatticeClass Lattice;
  Vec3 Twist;
  UBspline_3d_z *Spline;
public:
  void SetLattice(Mat3 A);
  inline LatticeClass& GetLattice() { return Lattice; }
  void SetTwist(Vec3 twist);
  void SetCoefs(Array<complex<double>,3> &coefs);
  complex<double> operator()(Vec3 r);
  OrbitalSplineClass() : Spline(NULL)
  {
  }
};


class TestTile
{
private:
  OrbitalSplineClass PrimSpline, SuperSpline;
  IOSectionClass PrimIO, SuperIO;
  int NumPrimTwists, NumPrimBands;
  int NumSuperTwists, NumSuperBands;
  double NormRatio;
public:
  void SetFiles (string primName, string superName);
  void ReadOrbs (int superkIndex, int superBandIndex);
  void CheckVals (int numVals);
};

void
OrbitalSplineClass::SetLattice (Mat3 A)
{
  Lattice.SetDirect(A);
}

void
OrbitalSplineClass::SetTwist (Vec3 twist)
{
  Twist = twist;
};

complex<double>
OrbitalSplineClass::operator()(Vec3 r) 
{
  Vec3 u = Lattice.r2u(r);
  Vec3 uf = u;
  for (int i=0; i<3; i++)
    uf[i] -= floor(uf[i]);
  
  complex<double> val;
  eval_UBspline_3d_z (Spline, uf[0], uf[1], uf[2], &val);
  
  double phase = -2.0*M_PI*dot(u, Twist);
  double s, c;
  sincos (phase, &s, &c);
  val *= complex<double>(c,s);
  return val;
}

void
OrbitalSplineClass::SetCoefs (Array<complex<double>,3> &coefs)
{
  BCtype_z xBC, yBC, zBC;
  xBC.lCode = xBC.rCode = PERIODIC;
  yBC.lCode = yBC.rCode = PERIODIC;
  zBC.lCode = zBC.rCode = PERIODIC;

  if (Spline != NULL) 
    destroy_Bspline (Spline);

  // Create spline
  Ugrid xGrid, yGrid, zGrid;
  xGrid.start=0.0;  xGrid.end=1.0;  xGrid.num=coefs.extent(0);
  yGrid.start=0.0;  yGrid.end=1.0;  yGrid.num=coefs.extent(1);
  zGrid.start=0.0;  zGrid.end=1.0;  zGrid.num=coefs.extent(2);
  Spline = create_UBspline_3d_z (xGrid, yGrid, zGrid, 
				 xBC, yBC, zBC, coefs.data());
}

void 
TestTile::SetFiles (string primName, string superName)
{
  assert (PrimIO.OpenFile (primName));
  assert (SuperIO.OpenFile (superName));
  Array<double,2> lattice;
  Mat3 A;

  assert (PrimIO.OpenSection("parameters"));
  assert (PrimIO.ReadVar ("lattice", lattice));
  assert (PrimIO.ReadVar ("num_twists", NumPrimTwists));
  assert (PrimIO.ReadVar ("num_bands", NumPrimBands));
  for (int i=0; i<3; i++) for (int j=0; j<3; j++)
    A(i,j) = lattice(i,j);
  PrimSpline.SetLattice (A);
  PrimIO.CloseSection();

  assert (SuperIO.OpenSection("parameters"));
  assert (SuperIO.ReadVar ("lattice", lattice));
  assert (SuperIO.ReadVar ("num_twists", NumSuperTwists));
  assert (SuperIO.ReadVar ("num_bands", NumSuperBands));
  for (int i=0; i<3; i++) for (int j=0; j<3; j++)
    A(i,j) = lattice(i,j);
  SuperSpline.SetLattice (A);
  SuperIO.CloseSection();
 
}

void 
TestTile::ReadOrbs (int superTwistIndex, int superBandIndex)
{
  Array<double,4> superVec, primVec;
  Array<complex<double>,3> superCoefs, primCoefs;
  Vec3 superTwist, primTwist;
  Array<double,1> twist;
  assert (SuperIO.OpenSection("eigenstates"));
  assert (SuperIO.OpenSection("twist", superTwistIndex));
  assert (SuperIO.ReadVar("twist_angle", twist));
  superTwist = Vec3 (twist(0), twist(1), twist(2));
  assert (SuperIO.OpenSection("band", superBandIndex));
  int primTwistIndex, primBandIndex;
  assert (SuperIO.ReadVar("prim_twist_index", primTwistIndex));
  assert (SuperIO.ReadVar("prim_band_index", primBandIndex));
  assert (SuperIO.ReadVar("eigenvector", superVec));
  SuperIO.CloseSection(); // "band"
  SuperIO.CloseSection(); // "twist"
  SuperIO.CloseSection(); // "eigenstates"


  assert (primTwistIndex < NumPrimTwists);
  assert (primBandIndex  < NumPrimBands);
  assert (PrimIO.OpenSection("eigenstates"));
  assert (PrimIO.OpenSection("twist", primTwistIndex));
  assert (PrimIO.ReadVar("twist_angle", twist));
  primTwist = Vec3 (twist(0), twist(1), twist(2));
  assert (PrimIO.OpenSection("band", primBandIndex));
  assert (PrimIO.ReadVar("eigenvector", primVec));
  PrimIO.CloseSection(); // "band"
  PrimIO.CloseSection(); // "twist"
  PrimIO.CloseSection(); // "eigenstates"

  superCoefs.resize(superVec.extent(0)-1, 
		    superVec.extent(1)-1,
		    superVec.extent(2)-1);
  for (int ix=0; ix<superCoefs.extent(0); ix++)
    for (int iy=0; iy<superCoefs.extent(1); iy++)
      for (int iz=0; iz<superCoefs.extent(2); iz++)
	superCoefs(ix,iy,iz) = complex<double>(superVec(ix,iy,iz,0),
					       superVec(ix,iy,iz,1));


  primCoefs.resize(primVec.extent(0)-1, 
		   primVec.extent(1)-1,
		   primVec.extent(2)-1);
  for (int ix=0; ix<primCoefs.extent(0); ix++)
    for (int iy=0; iy<primCoefs.extent(1); iy++)
      for (int iz=0; iz<primCoefs.extent(2); iz++)
	primCoefs(ix,iy,iz) = complex<double>(primVec(ix,iy,iz,0),
					      primVec(ix,iy,iz,1));
  BCtype_z xBC, yBC, zBC;
  xBC.lCode = xBC.rCode = PERIODIC;
  yBC.lCode = yBC.rCode = PERIODIC;
  zBC.lCode = zBC.rCode = PERIODIC;

  // Create splines
  PrimSpline.SetTwist  (primTwist);
  PrimSpline.SetCoefs  (primCoefs);
  SuperSpline.SetTwist (superTwist);
  SuperSpline.SetCoefs (superCoefs);
}


void 
TestTile::CheckVals (int numVals)
{
  double primVol, superVol;
  primVol  = fabs(det(PrimSpline.GetLattice().GetDirect()));
  superVol = fabs(det(SuperSpline.GetLattice().GetDirect()));
  double ratio = sqrt (superVol/primVol);
  complex<double> primVal, superVal;
  Vec3 uSuper, r;
  fprintf (stderr, "                 r                       prim val               super val:\n");
  for (int i=0; i<numVals; i++) {
    uSuper = Vec3 (drand48(), drand48(), drand48());
    r = SuperSpline.GetLattice().u2r(uSuper);
    primVal  =  PrimSpline (r);
    superVal = SuperSpline (r); 
    fprintf (stderr, " [%9.5f %9.5f %9.5f]  (%9.5f, %9.5f)  (%9.5f, %9.5f)\n", 
	     r[0], r[1], r[2], primVal.real(), primVal.imag(),
	     ratio*superVal.real(), ratio*superVal.imag());
  }
}


main(int argc, char **argv)
{
  list<ParamClass> argList;
  CommandLineParserClass parser (argList);
  bool success = parser.Parse (argc, argv);

  if (!success || parser.NumFiles() != 2) {
    cerr << "Usage.\n"
	 << "  TestTile PrimFile.h5 SuperFile.h5.\n";
    return -1;
  }
  
  TestTile tester;
  tester.SetFiles (parser.GetFile(0), parser.GetFile(1));
  tester.ReadOrbs (0, 179);
  tester.CheckVals(20);


}
