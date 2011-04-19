#include <Common/MPI/Communication.h>
#include "PlaneWaveClass.h"
#include "ParserClass.h"
#include "ParseCommand.h"
#include <Common/Splines/ComplexMultiTricubicSpline.h>
#include <algorithm>

void
CellClass::SetLattice (Mat3 A, bool isSuper)
{
  Lattice.SetDirect(A);
  IsSuper = isSuper;
}

void
CellClass::PutOrbital (OrbitalClass &orbital)
{
  zVec coefs = orbital.GetCoefs();
  kPointClass &kPoint = orbital.GetkPoint();
  kPoint.PutzVec (coefs, FFT, IsSuper);
}

void
CellClass::AddOrbital (OrbitalClass &orbital,
		       complex<double> prefactor)
{
  zVec coefs = orbital.GetCoefs();
  kPointClass &kPoint = orbital.GetkPoint();
  kPoint.AddzVec (coefs, prefactor, FFT, IsSuper);
}


void
TwistClass::Set (CellClass &cell, Vec3 k)
{
  Cell = &cell;
  Twist = cell.Lattice.k2Twist (k);
}

Vec3
TwistClass::Getk()
{
  return Cell->Lattice.Twist2k(Twist);
}

Vec3
TwistClass::GetTwist()
{
  return Twist;
}

Vec3
TwistClass::GetTwistFrac()
{
  Vec3 t = Twist;
  for (int i=0; i<3; i++)
    t[i] -= floor(t[i]);
  return t;
}

Vec3
TwistClass::GetTwistInt()
{
  return Vec3 (floor(Twist[0]), floor(Twist[1]), floor(Twist[2]));
}


void
kPointClass::SetCells (CellClass &primCell, 
		       CellClass &superCell, Vec3 k)
{
  PrimCell  =  &primCell;
  SuperCell = &superCell;
  PrimTwist = superCell.Lattice.k2Twist(k);
  Vec3 superTwist = superCell.Lattice.k2Twist(k);
  cerr << "superTwist = " << superTwist << endl;
  for (int i=0; i<3; i++) 
    SuperTwistInt[i]  = round (superTwist[i]);
  SuperTwistFrac = superTwist - SuperTwistInt;
}

Int3
kPointClass::GetFFTBoxSize (Array<Vec3,1> &primGVecs, double fftFactor)
{
  // This is the component of the k-vector of the primitive lattice
  // that is a G-vector of the superlattice.
  Vec3 Gshift = SuperCell->Lattice.Twist2k (SuperTwistInt);
  Array<Vec3,1> superGArray(primGVecs.size());
  for (int gi=0; gi<primGVecs.size(); gi++)
    superGArray(gi) = primGVecs(gi) + Gshift;
  
  Mat3 Asuper = SuperCell->Lattice.GetDirect();
  GVecsClass superGVecs;
  // Note:  this must be fixed to get the correct fft box size
  return (superGVecs.GetFFTBoxSize(Asuper, superGArray, fftFactor));
}

void
kPointClass::SetIndices (Array<Vec3,1> &primGVecs, Int3 fftBoxSize)
{
  PrimIndices.resize (PrimCell->GVecs.size());
  for (int i=0; i<PrimIndices.size(); i++) 
    PrimIndices(i) = PrimCell->GVecs.Index(i);

  // This is the component of the k-vector of the primitive lattice
  // that is a G-vector of the superlattice.
  Vec3 Gshift = SuperCell->Lattice.Twist2k (SuperTwistInt);
  Array<Vec3,1> superGArray(primGVecs.size());
  for (int gi=0; gi<primGVecs.size(); gi++)
    superGArray(gi) = primGVecs(gi) + Gshift;
  
  Mat3 Asuper = SuperCell->Lattice.GetDirect();
  GVecsClass superGVecs;
  superGVecs.Set(Asuper, superGArray, fftBoxSize);
  SuperIndices.resize(superGVecs.size());
  for (int i=0; i<superGVecs.size(); i++)
    SuperIndices(i) = superGVecs.Index(i);
}

void
kPointClass::PutzVec (zVec &c, FFTBox &fft, bool useSuper)
{
  Array<Int3,1> &indices = useSuper ? SuperIndices : PrimIndices;
  assert (c.size() == indices.size());
  fft.kBox = complex<double>();
  for (int gi=0; gi<c.size(); gi++)
    fft.kBox(indices(gi)) = c(gi);
}

void
kPointClass::AddzVec (zVec &c, complex<double> prefactor, 
		      FFTBox &fft, bool useSuper)
{
  Array<Int3,1> &indices = useSuper ? SuperIndices : PrimIndices;
  assert (c.size() == indices.size());
  for (int gi=0; gi<c.size(); gi++)
    fft.kBox(indices(gi)) += prefactor*c(gi);
}

void
PlaneWaveSystem::MakeRealSpaceOrbital (int spin, int ki, int band)
{
  CellClass &cell = WriteSuper ? SuperCell : PrimCell;
  FFTBox &FFT = cell.FFT;
  GVecsClass &GVecs = cell.GVecs;
  TinyVector<Array<OrbitalClass*,2>,2> &orbs = 
    WriteSuper ? SuperOrbitals : Orbitals;

  if (ShiftOrbitals) 
    ; // do something different
  else {
    if (band < NewOrbCoefs[spin].extent(1)) {
      FFT.kBox = complex<double>();
      for (int bi=0; bi<NewOrbCoefs[spin].extent(1); bi++)
	cell.AddOrbital (*orbs[spin](ki,bi), NewOrbCoefs[spin](ki,band,bi));
    }
    else 
      cell.PutOrbital (*orbs[spin](ki,band));
  }
   
  FFT.k2r();
}

void
PlaneWaveSystem::MakeRealSpaceBlip (int spin, int ki, int band)
{
  CellClass &cell = WriteSuper ? SuperCell : PrimCell;
  FFTBox &FFT = cell.FFT;
  GVecsClass &GVecs = cell.GVecs;
  TinyVector<Array<OrbitalClass*,2>,2> &orbs = 
    WriteSuper ? SuperOrbitals : Orbitals;

}


void
PlaneWaveSystem::WriteSpline (IO::IOSectionClass &out, int spin, int ki, int band)
{
  int nx, ny, nz;

  MakeRealSpaceOrbital (spin, ki, band);
  TinyVector<Array<OrbitalClass*,2>,2> &orbs = 
    WriteSuper ? SuperOrbitals : Orbitals;
  CellClass &cell = WriteSuper ? SuperCell : PrimCell;
  FFTBox &FFT = cell.FFT;
  FFT.GetDims (nx, ny, nz);

  double wfnorm = 0.0;
  double vol = fabs(det(cell.Lattice.GetDirect()));
  double elemvol = vol/(nx*ny*nz);
  double prefact = sqrt(1.0/(vol));
  OrbitalClass &orb = *orbs[spin](ki,band);

  double eigval = 0.0;
  if (band < NewOrbCoefs[spin].extent(1))
    for (int bi=0; bi<NewOrbCoefs[spin].extent(1); bi++)
      eigval += orbs[spin](ki,bi)->GetEigVal() * norm(NewOrbCoefs[spin](ki,band,bi));
  else
    eigval = orbs[spin](ki,band)->GetEigVal();

  if (Truncate && (orb.Center != NULL)) {
    cerr << "Please implement truncated orbitals in PlaneWaveSystem::WriteSpline.\n";
    abort();
  }
  else {
    Array<double,4> eigvec(nx+1, ny+1, nz+1, 2);
    for (int ix=0; ix<nx+1; ix++)
      for (int iy=0; iy<ny+1; iy++)
	for (int iz=0; iz<nz+1; iz++) {
	  wfnorm += elemvol*norm(prefact*(FFT.rBox(ix%nx,iy%ny,iz%nz)));
	  eigvec(ix, iy, iz, 0) = prefact*FFT.rBox(ix%nx,iy%ny,iz%nz).real();
	  eigvec(ix, iy, iz, 1) = prefact*FFT.rBox(ix%nx,iy%ny,iz%nz).imag();
	}
    
    // Finally write to file
    out.WriteVar ("eigenvector", eigvec);
  }
  out.WriteVar ("eigenvalue", eigval);
  out.WriteVar ("spin", spin);
  if (orb.Center != NULL) {
    Array<double,1> r(3);
    r(0) = orb.Center->r[0];  r(1)=orb.Center->r[1]; r(2)=orb.Center->r[2];
    out.WriteVar ("center", r);
    out.WriteVar ("radius", orb.Center->Radius);
  }
  // out.WriteVar ("prim_twist_index", PrimTwistIndex);
  // out.WriteVar ("prim_band_index",  PrimBandIndex);
}


void
PlaneWaveSystem::SetSkinThickness (double thickness)
{
  SkinThickness = thickness;
  for (int i=0; i<Centers.size(); i++)
    Centers[i].SkinThickness = thickness;
}

bool
PlaneWaveSystem::ReadIonPos(ParserClass &parser)
{
  bool success = parser.FindToken ("GEOMETRY");
  if (!success) return false;
  success = parser.FindToken ("cell");
  if (!success) return false;
  int numAtoms;
  success = parser.ReadInt (numAtoms);
  if (!success) return false;
  success = parser.FindToken ("(au)");
  if (!success) return false;

  PrimCell.IonPos.resize(numAtoms);
  PrimCell.AtomTypes.resize(numAtoms);
  for (int atom=0; atom<numAtoms; atom++) {
    success = parser.ReadInt(PrimCell.AtomTypes(atom));
    if (!success) return false;
    for (int i=0; i<3; i++) {
      success = parser.ReadDouble(PrimCell.IonPos(atom)[i]);
      if (!success) return false;
    }
  }
  fprintf (stderr,   "Atom Type           Position:\n");
  for (int atom=0; atom<numAtoms; atom++) {
    fprintf (stderr, "   %2d      %9.5f %9.5f %9.5f\n", PrimCell.AtomTypes(atom),
	     PrimCell.IonPos(atom)[0], PrimCell.IonPos(atom)[1], 
	     PrimCell.IonPos(atom)[2]);
  }
  fprintf (stderr, "\n");
  return true;
}

void
PlaneWaveSystem::TileIonPos()
{
  int numAtoms = PrimCell.IonPos.size();
  // Now, tile the supercell
  int numSuper = (int)round(fabs(det(TileMatrix)))*numAtoms;
  vector<Vec3> superPos;
  vector<int>  superTypes;
  Int3 iMax(5,5,5);
  for (int i0=-iMax[0]; i0<=iMax[0]; i0++)
    for (int i1=-iMax[1]; i1<=iMax[1]; i1++)
      for (int i2=-iMax[2]; i2<=iMax[2]; i2++) 
	for (int iat=0; iat<numAtoms; iat++) {
	  Vec3 r = PrimCell.IonPos(iat);
	  r += (double)i0*PrimCell.Lattice.a(0);
	  r += (double)i1*PrimCell.Lattice.a(1);
	  r += (double)i2*PrimCell.Lattice.a(2);
	  Vec3 uSuper = SuperCell.Lattice.r2u(r);
	  if ((uSuper[0] >= -1.0e-6) && (uSuper[0] < 0.9999) &&
	      (uSuper[1] >= -1.0e-6) && (uSuper[1] < 0.9999) &&
	      (uSuper[2] >= -1.0e-6) && (uSuper[2] < 0.9999)) {
// 	    fprintf (stderr, 
// 		     " [ %9.5f %9.5f %9.5f ]  r = [%9.5f %9.5f, %9.5f]\n",
// 		     uSuper[0], uSuper[1], uSuper[2],
// 		     r[0], r[1], r[2]);
	    superPos.push_back(r);
	    superTypes.push_back(PrimCell.AtomTypes(iat));
	  }
	}
  cerr << "There are " << superPos.size() << " atoms in the supercell.\n";

  SuperCell.IonPos.resize(superPos.size());
  SuperCell.AtomTypes.resize(superPos.size());
  for (int i=0; i<superPos.size(); i++) {
    SuperCell.IonPos(i)    =   superPos[i];
    SuperCell.AtomTypes(i) = superTypes[i];
  }
  int numExpected = (int)round(fabs(det(TileMatrix)))*numAtoms;
  if (numExpected != superPos.size()) {
    cerr << "Expected " << numExpected << " atoms in the supercell, but found "
	 << superPos.size() <<".  Aborting.\n";
    abort();
  }

}



bool
PlaneWaveSystem::ReadLattice(ParserClass &parser)
{
  bool success = parser.FindToken ("(au)");
  if (!success) return false;
  Mat3 Aprim;
  for (int i=0; i<3; i++)
    for (int j=0; j<3; j++) {
      success = parser.ReadDouble (Aprim(i,j));
      if (!success) return false;
    }
  fprintf (stderr, "Primitive lattice vectors:\n");
  for (int i=0; i<3; i++)
    fprintf (stderr, "  [ %9.6f %9.6f %9.6f ]\n", 
	    Aprim(i,0), Aprim(i,1), Aprim(i,2));

  PrimCell.SetLattice  (Aprim, false);
  Mat3 Asuper = TileMatrix*Aprim;
  SuperCell.SetLattice (Asuper, true);
  
  fprintf (stderr, "Superlattice vectors:\n");
  for (int i=0; i<3; i++)
    fprintf (stderr, "  [ %9.6f %9.6f %9.6f ]\n", 
	     Asuper(i,0), Asuper(i,1), Asuper(i,2));
  
  Mat3 recip = PrimCell.Lattice.GetRecip();
  fprintf (stderr, "Reciprocal lattice vectors:\n");
  for (int i=0; i<3; i++)
    fprintf (stderr, "  [ %9.6f %9.6f %9.6f ]\n", 
	    recip(i,0), recip(i,1), recip(i,2));
  fprintf (stderr, "\n");  
  return true;
}


bool
PlaneWaveSystem::ReadGVecs(ParserClass &parser)
{
  int numG;
  bool success = parser.FindToken("G-vectors");
  if (!success) return false;
  success = parser.ReadInt (numG);
  if (!success) return false;
  success = parser.FindToken ("(au)");
  if (!success) return false;

  GVecsArray.resize(numG);
  GInts.resize(numG);

  for (int iG=0; iG<numG; iG++) {
    for (int dim=0; dim<3; dim++) {
      success = parser.ReadDouble (GVecsArray(iG)[dim]);
      // HACK HACK HACK
      GVecsArray(iG)[dim] *= -1.0;
      if (!success) return false;
    }
  }
  perr << "Read " << numG << " G-vectors.\n";

  // Now compute integer values.
  double TwoPiInv = 1.0/(2.0*M_PI);
  Vec3 a0 = PrimCell.Lattice.a(0);
  Vec3 a1 = PrimCell.Lattice.a(1);
  Vec3 a2 = PrimCell.Lattice.a(2);
  Vec3 b0 = PrimCell.Lattice.b(0);
  Vec3 b1 = PrimCell.Lattice.b(1);
  Vec3 b2 = PrimCell.Lattice.b(2);
  for (int iG=0; iG<numG; iG++){
    GInts(iG)[0] = (int)round(TwoPiInv*dot(GVecsArray(iG), a0));
    GInts(iG)[1] = (int)round(TwoPiInv*dot(GVecsArray(iG), a1));
    GInts(iG)[2] = (int)round(TwoPiInv*dot(GVecsArray(iG), a2));
    Vec3 G = (double)GInts(iG)[0]*b0 + (double)GInts(iG)[1]*b1 +
      (double)GInts(iG)[2]*b2;
    Vec3 diff = GVecsArray(iG) - G;
    if (dot(diff,diff) > 1.0e-12) 
      perr << "Error:  mag = " << dot(diff,diff) << endl;
  }

  // Check this to see if the FFT box sizes are commensurate.
  PrimCell.GVecs.Set  ( PrimCell.Lattice.GetDirect(), GVecsArray, FFTFactor);
  SuperCell.GVecs.Set (SuperCell.Lattice.GetDirect(), GVecsArray, FFTFactor);
  return true;
}


bool
PlaneWaveSystem::ReadOrbitals(ParserClass &parser)
{
  if (!parser.FindToken("k-points")) return false;
  int numk, numUpBands, numDownBands, dummy;
  if (!parser.ReadInt (numk)) return false;
  if (!parser.FindToken("(au)")) return false;
  if (!parser.ReadInt (dummy)) return false;
  /// Number of up bands
  if (!parser.ReadInt (numUpBands)) return false;
  /// Number of down band -- not used if not spin-dependent
  if (!parser.ReadInt (numDownBands)) return false;

  Orbitals[0].resize(numk, numUpBands);
  if (numDownBands > 0)
    Orbitals[1].resize(numk,numDownBands);
  int numBands = numUpBands + numDownBands;
  PrimkPoints.resize(numk);
  for (int ik=0; ik<numk; ik++) {
    for (int iband=0; iband<numUpBands; iband++)
      Orbitals[0](ik,iband) = 
	new OrbitalClass(PrimCell, SuperCell, PrimkPoints(ik));
    for (int iband=0; iband<numDownBands; iband++)
      Orbitals[1](ik,iband) = 
	new OrbitalClass(PrimCell, SuperCell, PrimkPoints(ik));
  }

  parser.Reset();

  if (!parser.FindToken("k-points")) return false;

  Vec3 k;
  int kpoint, upbands, downbands, spin;
  double eigval;
  Int3 maxFFTBox (0,0,0);
  for (int ik=0; ik < numk; ik++) {
    if (!parser.FindToken("k-point"))          return false;
    if (!parser.FindToken("(au)"))             return false;
    if (!parser.ReadInt (kpoint))              return false;
    assert (kpoint == ik+1);
    if (!parser.ReadInt (upbands))             return false;
    assert (upbands == numUpBands); 
    if (!parser.ReadInt (downbands))           return false;
    assert (downbands == numDownBands);
    if (!parser.ReadDouble(k[0]))              return false;
    if (!parser.ReadDouble(k[1]))              return false;
    if (!parser.ReadDouble(k[2]))              return false;
    // HACK HACK HACK:  since we use the opposite FFT convention than
    // CASINO, k -> -k
    k = -1.0*k;
    PrimkPoints(ik).SetCells (PrimCell, SuperCell, k);
    Int3 fftBox = PrimkPoints(ik).GetFFTBoxSize(GVecsArray, FFTFactor);
    maxFFTBox[0] = max(fftBox[0], maxFFTBox[0]);
    maxFFTBox[1] = max(fftBox[1], maxFFTBox[1]);
    maxFFTBox[2] = max(fftBox[2], maxFFTBox[2]);
    int iband[2] = { 0, 0 };
    for (int band=0; band<numBands; band++) {
      if (!parser.FindToken ("Band"))          return false;
      if (!parser.FindToken ("\n"))            return false;
      if (!parser.ReadInt (dummy))             return false;
      if (!parser.ReadInt (spin))              return false;
      if (!parser.ReadDouble (eigval))         return false;
      if (!parser.FindToken ("\n"))            return false;
      if (!parser.FindToken ("coefficients"))  return false;
      // Make spin 0 or 1 instead 1 or 2
      spin--;
      Orbitals[spin] (ik, iband[spin])->SetEigVal (eigval);
      Orbitals[spin] (ik, iband[spin])->SetSpin (spin);
      Orbitals[spin] (ik, iband[spin])->SetPrimIndices (ik, band);
      if (!Orbitals[spin] (ik, iband[spin])->Read (parser)) return false;
      iband[spin]++;
    }
    // Gramm-Schmidt orthogonalize bands:
    for (int spin=0; spin<2; spin++) {
      int nBands = Orbitals[spin].extent(1);
      for (int m=0; m<nBands; m++) { 
	zVec &mCoefs = Orbitals[spin](ik, m)->GetCoefs();
	double nrm = 0.0;
	for (int i=0; i<mCoefs.size(); i++)
	  nrm += norm (mCoefs(i));
	mCoefs *= 1.0/sqrt(nrm);
	for (int n=m+1; n<nBands; n++) {
	  zVec &nCoefs = Orbitals[spin](ik, n)->GetCoefs();
	  complex<double> overlap(0.0, 0.0);
	  for (int i=0; i<mCoefs.size(); i++) 
	    overlap += conj(mCoefs(i))*nCoefs(i);
	  if (norm(overlap) > 1.0e-12) 
	    cerr << "Serious nonorthogonality for m=" << m << ", n=" << n 
		 << "  overlap = " << overlap << endl;
	  //	   complex<double> overlap = conjdot (mCoefs, nCoefs);
	  nCoefs -= overlap * mCoefs;
	}
      }
    }
  }
  cerr << "The max FFT box size is " << maxFFTBox << endl;
  cerr << "The supercell FFTBox is " << SuperCell.FFT.rBox.shape() << endl;

  // Now we know the maximum FFT box size that we need, so we can
  // setup the kPoint objects' Gvecs
  for (int ik=0; ik<numk; ik++)
    PrimkPoints(ik).SetIndices (GVecsArray, maxFFTBox);
  SuperCell.GVecs.Set (SuperCell.Lattice.GetDirect(), GVecsArray, maxFFTBox);
				

  perr << "Successfully read " << numk << " kpoints by " 
       << numBands << " bands.\n";
  return true;
}



bool
OrbitalClass::Read(ParserClass &parser)
{
  Coefs.resize(PrimCell.GVecs.size());
  double nrm = 0.0;
  for (int i=0; i<PrimCell.GVecs.size(); i++) {
    if (!parser.ReadComplex(Coefs(i))) return false;
    nrm += norm (Coefs(i));
  }
  Coefs *= sqrt(1.0/nrm);
  return true;
}


void
OrbitalClass::Read (IO::IOSectionClass &in)
{
  in.ReadVar ("eigenvalue", Eigenvalue);
  Array<double,2> eigvec;
  in.ReadVar ("eigenvector", eigvec);
  Coefs.resize(eigvec.extent(0));
  for (int i=0; i<eigvec.extent(0); i++)
    Coefs(i) = complex<double>(eigvec(i,0), eigvec(i,1));
}

void
OrbitalClass::Write(IO::IOSectionClass &out)
{
  out.WriteVar ("eigenvalue", Eigenvalue);
  Array<double,2> eigvec(PrimCell.GVecs.size(),2);
  for (int i=0; i<PrimCell.GVecs.size(); i++) {
    eigvec(i,0) = real(Coefs(i));
    eigvec(i,1) = imag(Coefs(i));
  }
  out.WriteVar ("eigenvector", eigvec);
  if (Center != NULL) {
    Array<double,1> r(3);
    r(0) = Center->r[0];  r(1)=Center->r[1]; r(2)=Center->r[2];
    out.WriteVar ("center", r);
    out.WriteVar ("radius", Center->Radius);
  }
  out.WriteVar ("prim_twist_index", PrimTwistIndex);
  out.WriteVar ("prim_band_index",  PrimBandIndex);
}

inline complex<float> operator*(double s, complex<float> c)
{
  return complex<float> (s*c.real(), s*c.imag());
}

void
OrbitalClass::CheckKineticEnergy()
{
  FFTBox &FFT = PrimCell.FFT;
  GVecsClass &GVecs = PrimCell.GVecs;

  // Calculate kinetic energy in reciprocal space:
  double ke = 0.0;
  for (int i=0; i<GVecs.size(); i++)
    ke += dot(GVecs(i),GVecs(i))*norm(Coefs(i));
  //  ke /= GVecs.GetBoxVol();

  int nx, ny, nz;
  FFT.GetDims(nx, ny, nz);
  Mat3 box = GVecs.GetLattice();
  LinearGrid xGrid(0.0, box(0,0), nx+1);
  LinearGrid yGrid(0.0, box(1,1), ny+1);
  LinearGrid zGrid(0.0, box(2,2), nz+1);
  Array<complex<double>,4> data(nx+1, ny+1, nz+1, 1);
  
  FFT.PutkVec (Coefs);
  FFT.k2r();
  double wfnorm = 0.0;
  double vol = GVecs.GetBoxVol();
  double elemvol = vol/(nx*ny*nz);
  double prefact = sqrt(1.0/(vol));
  for (int ix=0; ix<nx; ix++)
    for (int iy=0; iy<ny; iy++)
      for (int iz=0; iz<nz; iz++) {
	wfnorm += elemvol*norm(prefact*(FFT.rBox(ix,iy,iz)));
	data(ix, iy, iz, 0) = prefact*FFT.rBox(ix,iy,iz);
      }
  MakePeriodic(data);

  ComplexMultiTricubicSpline spline (&xGrid, &yGrid, &zGrid, data, true);
  Array<double,1> KE(1), norm(1);
  Array<complex<double>,1> val(1), laplacian(1);
  spline.Norm(norm);
  spline.KineticEnergy(KE);
  // Quick KE check:
  double pwKE, spKE;
  double pwNorm, spNorm;
  Array<complex<double>,1> dummy(1);
  double spSum = 0.0;
  double pwSum = 0.0;
  double pwspSum = 0.0;
  perr << "Quick start:\n";
  for (int i=0; i<1000; i++) {
    double x = drand48()*box(0,0);
    double y = drand48()*box(1,1);
    double z = drand48()*box(2,2);
    Vec3 r(x,y,z);
    complex<double> pwLap (0.0, 0.0), spLap (0.0, 0.0);
    spline.Laplacian(x,y,z,dummy);
    spLap = dummy(0);
    for (int j=0; j<GVecs.size(); j++) {
      double phase = dot(GVecs(j),r);
      complex<double> e2iGr(cos(phase),-sin(phase));
      pwLap -= 1.0/sqrt(vol)*dot(GVecs(j),GVecs(j))*Coefs(j)*e2iGr;
    }
    pwSum += real(conj(pwLap)*pwLap);
    spSum += real(conj(spLap)*spLap);
    pwspSum += real(conj(spLap)*pwLap);
  }
  perr << "1-alpha = " << 1-pwspSum/sqrt(pwSum*spSum) << endl;

  complex<double> quickKE = 0.0;
  for (int ix=0; ix<nx; ix++) {
    double x = 0.9*xGrid(ix)+0.1*xGrid(ix+1);
    for (int iy=0; iy<ny; iy++) {
      double y = 0.9*yGrid(iy)+0.1*yGrid(iy+1);
      for (int iz=0; iz<nz; iz++) {
	double z = 0.9*zGrid(iz)+0.1*zGrid(iz+1);
	spline.Laplacian(x,y,z,laplacian);
	spline(x,y,z,val);
	quickKE -= elemvol*conj(laplacian(0))*val(0);
      }
    }
  }
	

  fprintf (stderr, "Spline norm   = %14.10f\n", norm(0));
  fprintf (stderr, "Plane wave KE = %14.10f\n", ke);
  fprintf (stderr, "Spline KE     = %14.10f\n", -KE(0)/norm(0));
  fprintf (stderr, "Quick KE      = %14.10f\n", real(quickKE)/norm(0));
  fprintf (stderr, "Fractional Error = %1.4e\n", -KE(0)/norm(0)/ke - 1.0);
}

void
OrbitalClass::WriteSpline(IO::IOSectionClass &out, bool real,
			  bool shift, bool truncate,
			  bool writeSuper)
{
  CellClass &cell   = writeSuper ? SuperCell : PrimCell;

  GVecsClass &GVecs = cell.GVecs;
  FFTBox &FFT       = cell.FFT;

  Mat3 lattice = cell.Lattice.GetDirect();
  out.WriteVar ("eigenvalue", Eigenvalue);
  out.WriteVar ("spin", Spin);
  out.WriteVar ("prim_twist_index", PrimTwistIndex);
  out.WriteVar ("prim_band_index",  PrimBandIndex);
  int nx, ny, nz;
  FFT.GetDims(nx, ny, nz);
  if (Center == NULL || !shift) {
    Array<complex<double>,1> shiftedCoefs(Coefs.size());
    for (int i=0; i<Coefs.size(); i++) {
      Int3 in = GVecs.Index(i);
      shiftedCoefs(i) = Coefs(i);
    }
    kPoint.PutzVec (shiftedCoefs, FFT, writeSuper);
  }
  else {
    Array<complex<double>,1> shiftedCoefs(Coefs.size());
    for (int i=0; i<Coefs.size(); i++) {
      double phase = -dot (FFT.GVecs(i), Center->r);
      complex<double> e2iGr (cos(phase), sin(phase));
      shiftedCoefs(i) = e2iGr * Coefs(i);
      Int3 in = GVecs.Index(i);
      if ((in[0]+in[1]+in[2])%2==1)
	shiftedCoefs(i) *= -1.0;
    }
    kPoint.PutzVec (shiftedCoefs, FFT, writeSuper);
  }
  FFT.k2r();

  if (truncate && Center!=NULL) {
    int ixMin, ixMax, iyMin, iyMax, izMin, izMax;
    // First, find minimum and maximum grid points to include
    ixMin = ixMax = nx/2;
    iyMin = iyMax = ny/2;
    izMin = izMax = nz/2;
    Vec3 c = 0.5*(Vec3(lattice(0,0), lattice(0,1), lattice(0,2)) +
		  Vec3(lattice(1,0), lattice(1,1), lattice(1,2)) +
		  Vec3(lattice(2,0), lattice(2,1), lattice(2,2)));
    double nxInv = 1.0/(double)nx;
    double nyInv = 1.0/(double)ny;
    double nzInv = 1.0/(double)nz;  
    for (int ix=0; ix<nx; ix++) {
      double sx = nxInv*(double)ix;
      for (int iy=0; iy<ny; iy++) {
	double sy = nyInv*(double)iy;
	for (int iz=0; iz<nz; iz++) {
	  double sz = nzInv*(double)iz;
	  Vec3 r = (sx*Vec3(lattice(0,0), lattice(0,1), lattice(0,2))+
		    sy*Vec3(lattice(1,0), lattice(1,1), lattice(1,2))+
		    sz*Vec3(lattice(2,0), lattice(2,1), lattice(2,2)));
	  Vec3 diff = r - c;
	  if (dot(diff,diff) <= Center->Radius*Center->Radius) {
	    ixMin = min (ixMin, ix); ixMax = max(ixMax, ix);
	    iyMin = min (iyMin, iy); iyMax = max(iyMax, iy);
	    izMin = min (izMin, iz); izMax = max(izMax, iz);
	  }
	}
      }
    }
    perr << "ixMin = " << ixMin << "  ixMax = " << ixMax << endl;
    perr << "iyMin = " << iyMin << "  iyMax = " << iyMax << endl;
    perr << "izMin = " << izMin << "  izMax = " << izMax << endl;
    double x0 = (double)ixMin/(double)nx; double x1=(double)ixMax/(double)nx;
    double y0 = (double)iyMin/(double)ny; double y1=(double)iyMax/(double)ny;
    double z0 = (double)izMin/(double)ny; double z1=(double)izMax/(double)nz;
    LinearGrid xGrid(x0,x1,ixMax-ixMin+1);
    LinearGrid yGrid(y0,y1,iyMax-iyMin+1);
    LinearGrid zGrid(y0,z1,izMax-izMin+1);
    out.NewSection("Grid"); xGrid.Write(out); out.CloseSection();
    out.NewSection("Grid"); yGrid.Write(out); out.CloseSection();
    out.NewSection("Grid"); zGrid.Write(out); out.CloseSection();
    if (real) {
      Array<double,3> eigvec(ixMax-ixMin+1, iyMax-iyMin+1, izMax-izMin+1);
      for (int ix=ixMin; ix<=ixMax; ix++)
	for (int iy=iyMin; iy<=iyMax; iy++)
	  for (int iz=izMin; iz<=izMax; iz++)
	    eigvec(ix-ixMin, iy-iyMin, iz-izMin) = FFT.rBox (ix,iy,iz).real();
      out.WriteVar("eigenvector", eigvec);
    }
    else {
      Array<double,4> eigvec(ixMax-ixMin+1, iyMax-iyMin+1, izMax-izMin+1, 2);
      for (int ix=ixMin; ix<=ixMax; ix++)
	for (int iy=iyMin; iy<=iyMax; iy++)
	  for (int iz=izMin; iz<=izMax; iz++) {
	    eigvec(ix-ixMin, iy-iyMin, iz-izMin, 0) = FFT.rBox (ix,iy,iz).real();
	    eigvec(ix-ixMin, iy-iyMin, iz-izMin, 1) = FFT.rBox (ix,iy,iz).imag();
	  }
      out.WriteVar("eigenvector", eigvec);
    }
  }
  else {
    double wfnorm = 0.0;
    double vol = GVecs.GetBoxVol();
    double elemvol = vol/(nx*ny*nz);
    double prefact = sqrt(1.0/(vol));
    if (!real) {
      Array<double,4> eigvec(nx+1, ny+1, nz+1, 2);
      for (int ix=0; ix<nx+1; ix++)
	for (int iy=0; iy<ny+1; iy++)
	  for (int iz=0; iz<nz+1; iz++) {
	    wfnorm += elemvol*norm(prefact*(FFT.rBox(ix%nx,iy%ny,iz%nz)));
	    eigvec(ix, iy, iz, 0) = prefact*FFT.rBox(ix%nx,iy%ny,iz%nz).real();
	    eigvec(ix, iy, iz, 1) = prefact*FFT.rBox(ix%nx,iy%ny,iz%nz).imag();
	  }

      // Finally write to file
      out.WriteVar ("eigenvector", eigvec);
    }
    else {
      Array<double,3> eigvec(nx+1, ny+1, nz+1);
      for (int ix=0; ix<nx+1; ix++)
	for (int iy=0; iy<ny+1; iy++)
	  for (int iz=0; iz<nz+1; iz++) 
	    eigvec(ix, iy, iz) = FFT.rBox(ix%nx,iy%ny,iz%nz).real();

      // Finally write to file
      out.WriteVar ("eigenvector", eigvec);
    }
  }
  if (Center != NULL) {
    Array<double,1> r(3);
    r(0) = Center->r[0];  r(1)=Center->r[1]; r(2)=Center->r[2];
    out.WriteVar ("center", r);
    out.WriteVar ("radius", Center->Radius);
  }
}

bool
PlaneWaveSystem::ReadRunInfo (ParserClass &parser)
{
  if (!parser.FindToken ("Generated by:\n"))        return false;
  if (!parser.ReadLine(GeneratedBy))                return false;
  if (!parser.FindToken ("Method:"))                return false;
  if (!parser.ReadWord (Method))                    return false;
  if (!parser.FindToken("DFT Functional:\n"))       return false;
  if (!parser.ReadLine(Functional))                 return false;
  if (!parser.FindToken ("Pseudopotential"))        return false;
  if (!parser.FindToken ("\n"))                     return false;
  if (!parser.ReadLine (Pseudopotential))           return false;
  if (!parser.FindToken ("Plane wave cutoff (au)")) return false;
  if (!parser.ReadDouble (ECut))                    return false;
  if (!parser.FindToken ("Spin polarized:"))        return false;
  string tmp;
  if (!parser.ReadWord (tmp))                       return false;
  if ((tmp==".false.") || tmp=="F") 
    SpinPolarized = false;
  else if ((tmp==".true.") || (tmp=="T"))
    SpinPolarized = true;
  else {
    perr << "Unrecognized boolean string " << tmp << " in ReadRunInfo.\n";
    abort();
  }
  if (!parser.FindToken ("Total energy (au per primitive cell)"))
    return false;
  if (!parser.ReadDouble (TotalEnergy))             return false;
  if (!parser.FindToken ("Kinetic energy"))         return false;
  if (!parser.FindToken ("\n"))                     return false;
  if (!parser.ReadDouble (KineticEnergy))           return false;
  if (!parser.FindToken ("Local potential energy")) return false;
  if (!parser.FindToken ("\n"))                     return false;
  if (!parser.ReadDouble(LocalPotEnergy))           return false;
  if (!parser.FindToken ("Non"))                    return false;
  if (!parser.FindToken ("local potential energy")) return false;
  if (!parser.FindToken ("\n"))                     return false;
  if (!parser.ReadDouble(NonlocalEnergy))           return false;
  if (!parser.FindToken ("Electron"))               return false;
  if (!parser.FindToken ("electron energy"))        return false;
  if (!parser.FindToken("\n"))                      return false;
  if (!parser.ReadDouble (eeEnergy))                return false;
  if (!parser.FindToken ("Ion"))                    return false;
  if (!parser.FindToken ("ion energy"))             return false;
  if (!parser.FindToken ("\n"))                     return false;
  if (!parser.ReadDouble (IonIonEnergy))            return false;
  if (!parser.FindToken ("Number of electrons per primitive cell"))
    return false;
  if (!parser.ReadInt (NumElectrons))               return false;
  return true;
}
      
struct EigValLess
{
  inline bool operator() (OrbitalClass *orb1, OrbitalClass *orb2)
  {    return orb1->GetEigVal() < orb2->GetEigVal();  }
};



void
PlaneWaveSystem::SetupSuperOrbitals()
{
  // First, divide the primitive k-vectors into sets of supercell
  // k-vectors. 
  // A vector holding the unique fractional parts of the supercell
  // twist vectors
  vector<Vec3> superFracs;
  // This holds to which supercell kpoint each primitive k-point belongs
  Array<int,1> superIndex(PrimkPoints.size());
  for (int ki=0; ki < PrimkPoints.size(); ki++) {
    Vec3 frac = PrimkPoints(ki).GetSuperTwistFrac();
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
  cerr << "Detected " << superFracs.size() 
       << " distinct supercell twist vector"
       << ((superFracs.size()>1) ? "s.\n" : ".\n");

  int numSuperTwists = superFracs.size();

  // For each supercell twist, create a list of primitive twists which
  // belong to it.
  Array<vector<int>,1> superSets(numSuperTwists);
  for (int ki=0; ki < PrimkPoints.size(); ki++) 
    superSets(superIndex(ki)).push_back(ki);
  
  for (int si=0; si<numSuperTwists; si++) {
    fprintf (stderr, "Super twist #%d:  [ %9.5f %9.5f %9.5f ]\n",
	     si, superFracs[si][0], superFracs[si][1], superFracs[si][2]);
    fprintf (stderr, "  Using k-points: ");
    for (int i=0; i<superSets(si).size(); i++) 
      fprintf (stderr, " %d", superSets(si)[i]);
    fprintf (stderr, "\n");
  }
  
  // Now check to see that each supercell twist has the right twists
  // to tile the primitive cell orbitals.
  int numTwistsNeeded = (int)round (fabs(det(TileMatrix)));
  for (int si=0; si<numSuperTwists; si++) {
    // First make sure we have enough points
    if (superSets(si).size() != numTwistsNeeded) {
      fprintf (stderr, "Super twist %d should own %d k-points, but owns %d.\n",
	       si, numTwistsNeeded, superSets(si).size());
      abort();
    }
    // Now, make sure they are all distinct
    int N = superSets(si).size();
    for (int i=0; i<N; i++) {
      Vec3 iInt = PrimkPoints(superSets(si)[i]).GetSuperTwistInt();
      for (int j=i+1; j<N; j++) {
	Vec3 jInt = PrimkPoints(superSets(si)[j]).GetSuperTwistInt();
	if (dot(iInt-jInt, iInt-jInt) < 1.0e-6) {
	  cerr << "Identical k-points detected in super twist set "
	       << si << endl;
	  abort();
	}
      }
    }
  }

  ////////////////////////////////////
  // Construct super twist orbitals //
  ////////////////////////////////////
  int numPrimBands  = Orbitals[0].extent(1);
  int numSuperBands = numTwistsNeeded * numPrimBands;
  int maxSpin = SpinPolarized ? 1 : 0;
  for (int spin=0; spin<=maxSpin; spin++) 
    SuperOrbitals[spin].resize(numSuperTwists, numSuperBands);
  for (int si=0; si < numSuperTwists; si++) {
    for (int spin=0; spin<=maxSpin; spin++) {
      // Create a vector of all the orbitals belonging to this
      // supertwist.  We will then sort them by energy
      vector<OrbitalClass*> myOrbitals;
      for (int i=0; i<superSets(si).size(); i++) {
	int ki = superSets(si)[i];
	for (int bi=0; bi<numPrimBands; bi++)
	  myOrbitals.push_back(Orbitals[spin](ki,bi));
      }
      EigValLess comparison;
      sort (myOrbitals.begin(), myOrbitals.end(), comparison);
      for (int i=0; i<myOrbitals.size(); i++)
	SuperOrbitals[spin](si,i) = myOrbitals[i];
    }
  }


}


bool 
PlaneWaveSystem::Read (string fname)
{
  if (Comm.MyProc() != 0)
    return true;
  MemParserClass memParser;
  FileParserClass2 fileParser;
  streamsize fsize = memParser.FileSize(fname);
  if (fsize == -1)
    return false;
  
  ParserClass *parserPtr;
  if (fsize < (streamsize)(1<<28))
    parserPtr = &memParser;
  else
    parserPtr = &fileParser;

  ParserClass &parser = *parserPtr;

  bool success;

  if (!parser.OpenFile (fname))         return false;
//   // Read the plane-wave energy cutoff
//   if (!parser.FindToken("Plane wave"))  return false;
//   if (!parser.FindToken("(au)"))        return false;
//   if (!parser.ReadDouble(ECut))         return false;

//   // Read the number of electrons
//   if (!parser.FindToken ("Number of electrons per primitive cell"))
//     return false;
//   if (!parser.ReadInt (NumElectrons))   return false;
  if (!ReadRunInfo (parser)) {
    perr << "Error reading run info.\n";
    return false;
  }
  if (!ReadIonPos (parser))             return false;
  if (!ReadLattice (parser))            return false;
  TileIonPos();
  if (!ReadGVecs (parser))              return false;
  if (!ReadOrbitals(parser))            return false;
  parser.CloseFile();
  
  // Check the k-point mesh to see if it's valid
  SetupSuperOrbitals();
//   bool meshOK = CheckkMesh(0);
//   Int3 upMesh = kPointMesh;
//   if (SpinPolarized) {
//      meshOK = meshOK && CheckkMesh(1);
//     if (upMesh != kPointMesh) {
//       cerr << "Up and down k-point meshes to not match.\n";
//       abort();
//     }
//   }
//   if (!meshOK) {
//     cerr << "Error in k-point mesh in file " << fname << endl;
//     abort();
//   }
  
  return true;
}


bool
PlaneWaveSystem::Write (string fname)
{
  CellClass &cell = WriteSuper ? SuperCell : PrimCell;
  FFTBox &FFT = cell.FFT;
  GVecsClass &GVecs = cell.GVecs;
  TinyVector<Array<OrbitalClass*,2>,2> &orbs =
    WriteSuper ? SuperOrbitals : Orbitals;

  if (Comm.MyProc() != 0)
    return true;
  cell.SetupFFT();

  IO::IOSectionClass out;
  if (!out.NewFile (fname)) return false;
  Array<int,1> version(2);
  version(0) = 0;
  version(1) = 20;
  out.WriteVar("version", version);

  out.NewSection("parameters");
  Array<double,2> lattice(3,3);
  for (int i=0; i<3; i++)
    for (int j=0; j<3; j++)
      lattice(i,j) = cell.Lattice.GetDirect()(i,j);
  out.WriteVar("lattice", lattice);
  Array<double,2> reciprocal_lattice(3,3);
  for (int i=0; i<3; i++)
    for (int j=0; j<3; j++)
      reciprocal_lattice(i,j) = cell.Lattice.GetRecip()(i,j);
  out.WriteVar("reciprocal_lattice", reciprocal_lattice);
  int complex_coefficients = 1;
  out.WriteVar ("complex_coefficients", complex_coefficients);
  out.WriteVar ("num_bands",      orbs[0].extent(1)+orbs[1].extent(1));
  out.WriteVar ("num_up_bands",   orbs[0].extent(1));
  out.WriteVar ("num_down_bands", orbs[1].extent(1));
  out.WriteVar ("num_twists",     orbs[0].extent(0));
  int num_spins = SpinPolarized ? 2 : 1;
  out.WriteVar ("num_spins", num_spins);
  out.WriteVar ("maximum_ecut", ECut);
  int numElecs = NumElectrons;
  if (WriteSuper)
    numElecs *= (int)round(fabs(det(TileMatrix)));
  out.WriteVar ("num_electrons", numElecs);
  out.CloseSection(); // "parameters"

  out.NewSection("basis");
  if (Spline)
    out.WriteVar ("type", "spline");
  else {
    out.WriteVar ("type", "planewaves");
    out.WriteVar ("num_planewaves", GVecs.size());
    Array<double,2> gvecs (GVecs.size(), 3);
    Array<int,2> multipliers (GVecs.size(), 3);
    for (int i=0; i<GVecs.size(); i++)
      for (int j=0; j<3; j++) {
	gvecs(i,j) = GVecs(i)[j];
	multipliers(i,j) = GVecs.Multiplier(i)[j];
      }
    out.WriteVar ("planewaves", gvecs);
    out.WriteVar ("multipliers", multipliers);
  }
  out.CloseSection(); // "basis"

  out.NewSection("ions");
  Array<double,2> pos(cell.IonPos.size(),3);
  for (int i=0; i<cell.IonPos.size(); i++) {
    
    pos(i,0) = cell.IonPos(i)[0];
    pos(i,1) = cell.IonPos(i)[1];
    pos(i,2) = cell.IonPos(i)[2];
  }
  out.WriteVar ("pos", pos);
  out.WriteVar ("atom_types", cell.AtomTypes);
  out.CloseSection(); // "ions"
  

  out.NewSection ("eigenstates");
  for (int ik=0; ik<orbs[0].extent(0); ik++) {
    out.NewSection("twist");
    Vec3 t;
    if (WriteSuper)
      t = orbs[0](ik,0)->GetkPoint().GetSuperTwistFrac();
    else
      t = PrimkPoints(ik).GetPrimTwist();
    
    Array<double,1> twist(3);
    twist(0)=t[0]; twist(1)=t[1]; twist(2)=t[2];
    out.WriteVar ("twist_angle", twist);
    for (int spin=0; spin<2; spin++) {
      for (int band=0; band<orbs[spin].extent(1); band++) {
	out.NewSection("band");      

	if (WriteSuper) {
	  Vec3 ts = orbs[spin](ik,band)->GetkPoint().GetSuperTwistInt();
	  twist(0)=ts[0]; twist(1)=ts[1]; twist(2)=ts[2];
	  out.WriteVar ("super_twist_int", twist);
	}

	if (Spline) {
	  if (CheckKE)
	    orbs[spin](ik, band)->CheckKineticEnergy();
	  WriteSpline (out, spin, ik, band);
// 	  orbs[spin](ik, band)->WriteSpline 
// 	    (out, Real, ShiftOrbitals, Truncate, WriteSuper);
	}
	else
	  orbs[spin](ik, band)->Write (out);
	out.CloseSection(); // "band"
      }
    } // spin loop
    out.CloseSection(); // "twist"
  }
  out.CloseSection(); // "eigenstates"
  
  out.CloseFile();
  return true;
}



bool
PlaneWaveSystem::Read (IO::IOSectionClass &in) 
{
  Array<int,1> version;
  assert (in.ReadVar ("version", version));
  if (version(0) > 0 || version(1) > 11) {
    perr << "WF file version > program version.\n";
    abort();
  }
  
  assert (in.OpenSection("parameters"));
  Array<double,2> lattice;
  assert (in.ReadVar ("lattice", lattice));
  assert (lattice.extent(0) == 3);
  assert (lattice.extent(1) == 3);
  Mat3 primLattice, superLattice;
  for (int i=0; i<3; i++) 
    for (int j=0; j<3; j++)
      primLattice(i,j) = lattice(i,j);
  PrimCell.SetLattice(primLattice, false);
  SuperCell.SetLattice(TileMatrix*primLattice, true);
  
  assert (in.ReadVar ("reciprocal_lattice", lattice));
  assert (lattice.extent(0) == 3);
  assert (lattice.extent(1) == 3);
  for (int i=0; i<3; i++)
    for (int j=0; j<3; j++)
      assert (fabs(PrimCell.Lattice.GetRecip()(i,j) - lattice(i,j)) < 1.0e-6);
  int isComplex, numBands, numTwists, numSpins, numElecs;
  assert (in.ReadVar("complex_coefficients", isComplex));
  assert (in.ReadVar("num_bands", numBands));
  assert (in.ReadVar("num_twists", numTwists));
  assert (in.ReadVar("num_spins", numSpins));
  assert (in.ReadVar("num_electrons", numElecs));
  assert (in.ReadVar("maximum_ecut", ECut));
  in.CloseSection(); // "parameters"
  
  assert (in.OpenSection("basis"));

  in.CloseSection(); // "basis"
  string basisType;
  int numG;
  assert (in.ReadVar("type", basisType));
  assert (basisType == "planewaves");
  assert (in.ReadVar("num_planewaves", numG));
  Array<double,2> gVecsArray;
  Array<Vec3,1> gVecs;
  assert (in.ReadVar("planewaves", gVecsArray));
  assert (gVecsArray.extent(0) == numG);
  assert (gVecsArray.extent(1) == 3);
  gVecs.resize(numG);
  for (int i=0; i<numG; i++)
    gVecs(i) = Vec3 (gVecsArray(i,0), gVecsArray(i,1), gVecsArray(i,2));
  //  GVecs.Set (Lattice, gVecs, FFTFactor);
  
  return true;
}



void TestParse (int argc, char **argv)
{
  list<ParamClass> argList;
  
  argList.push_back(ParamClass("spline", false));
  argList.push_back(ParamClass("real",   false));

  CommandLineParserClass parser (argList);

  parser.Parse (argc, argv);

}


#include <time.h>

Mat3
string2Mat3 (string str)
{
  Mat3 mat;
  stringstream strin;
  strin << str;
  for (int i=0; i<3; i++)
    for (int j=0; j<3; j++) {
      int val;
      strin >> val;
      mat(i,j) = (double) val;
      if (strin.fail()) {
	cerr << "Error reading tilemat string.  Aborting.\n";
	abort();
      }
    } 
  return mat;
}

Int3 
string2Int3 (string str)
{
  Int3 vec;
  stringstream strin;
  strin << str;
  for (int i=0; i<3; i++) {
    int val;
    strin >> val;
    vec[i] = val;
    if (strin.fail()) {
      cerr << "Error reading tilemat string.  Aborting.\n";
      abort();
    }
  }
  return vec;
 }


void
PlaneWaveSystem::SetTileMatrix (Mat3 tileMatrix)
{
  TileMatrix = tileMatrix;
  WriteSuper = true;
}


main(int argc, char **argv) 
{
  COMM::Init (argc, argv);

  list<ParamClass> argList;

  argList.push_back(ParamClass("spline", false));
  argList.push_back(ParamClass("real",   false));
  argList.push_back(ParamClass("check",  false));
  argList.push_back(ParamClass("factor", true));
  argList.push_back(ParamClass("localize", true));
  argList.push_back(ParamClass("radius", true));
  argList.push_back(ParamClass("shift", false));
  argList.push_back(ParamClass("optimize-centers", false));
  argList.push_back(ParamClass("optimize-radii", false));
  argList.push_back(ParamClass("ortho", false));
  argList.push_back(ParamClass("truncate", false));
  argList.push_back(ParamClass("unfold",false));
  argList.push_back(ParamClass("bwfn", true));
  argList.push_back(ParamClass("split", 3));
  argList.push_back(ParamClass("pwfn", true));
  argList.push_back(ParamClass("qmcPACK", true));
  argList.push_back(ParamClass("skin-thickness", true));
  argList.push_back(ParamClass("fold-factor", 3));
  argList.push_back(ParamClass("tile", 3));
  argList.push_back(ParamClass("tilemat", 1));

  CommandLineParserClass parser (argList);
  bool success = parser.Parse (argc, argv);

  if (!success || parser.NumFiles() != 1) {
    perr << "Usage:  wfconvert [options...] infile\n"
	 << "Options:\n"
	 << "  --localize fname        construct localized linear combinations of\n"
	 << "        or random         the occupied bands\n"
	 << "  --radius x              Radius of localization\n"
	 << "  --skin-thickness x      Skin thickness for smooth truncation\n"
	 << "  --spline                Fourier transform orbitals to real space\n"
	 << "  --real                  take the real part of spline data\n"
         << "  --factor x              Multiply the number of real-space mesh points\n" 
	 << "                          by x in each direction.  x can be fractional\n"
	 << "  --shift                 shift localized and spline orbitals to the \n"
	 << "                          center of the simulation cell\n"
	 << "  --optimize-centers      optimize the localization centers\n" 
	 << "  --optimize-radii x      optimize the localization radii to contain x\n"
	 << "                          fraction of the norm\n"
	 << "  --ortho                 orthogonalize the local orbitals\n"
	 << "  --truncate              truncate localized orbitals beyond localization\n"
	 << "                          radius\n"
	 << "  --unfold                unfold k-point grid into a single k-point  supercell\n"
	 << "  --fold-factor nx ny nz  Partially unfold by factor nx x ny x nz.  Otherwise \n"
         << "                          unfold completely to a single k-point.\n"
	 << "  --tile nx ny nz         Tile the primitive cell by nx x ny x nz.\n"
         << "                          Do not really unfold, just update labeling and twists.\n"
         << "  --tilemat \"S11 S12 S13 S21 S22 S23 S31 S32 S33\"\n"
	 << "                          Tile the primitive cell by the integer matrix S.\n"
	 << "  --bwfn fname            write CASINO blip wave function file\n"
	 << "  --pwfn fname            write CASINO plane wave function file\n"
	 << "  --qmcPACK fname         write qmcPACK wave function file\n";     

    exit(1); 
  }
  
  string inName = parser.GetFile(0);
  double factor(1.0), radius(3.0);
  bool spline  = parser.Found("spline");
  bool real    = parser.Found("real");
  bool checkKE = parser.Found("check");
  bool shift   = parser.Found("shift");
  if (parser.Found("factor")) {
    factor = atof (parser.GetArg("factor").c_str());
    perr << "factor = " << factor << endl;
  }
  if (parser.Found("radius")) {
    radius = atof (parser.GetArg("radius").c_str());
    perr << "radius = " << radius << endl;
  }
  
  
  PlaneWaveSystem system;
  system.Localized       = parser.Found("localize");
  system.Spline          = parser.Found("spline");
  system.OptimizeCenters = parser.Found("optimize-centers");
  system.OptimizeRadii   = parser.Found("optimize-radii");
  system.ShiftOrbitals   = parser.Found("shift");
  system.Real            = parser.Found("real");
  system.CheckKE         = parser.Found("check");
  system.Orthogonalize   = parser.Found("ortho");
  system.Truncate        = parser.Found("truncate");
  bool unfold            = parser.Found("unfold");
  system.SetFFTFactor (factor);
  
  Int3 foldFactor (0,0,0), tileFactor(0,0,0), split(0,0,0);
  if (parser.Found("fold-factor")) {
    unfold = true;
    const char *arg;
    char *endptr;
    arg = parser.GetArg("fold-factor", 0).c_str();
    foldFactor[0] = (int)strtol (arg, &endptr, 10);
    assert (arg[0] != '\0' && (*endptr == '\0'));
    
    arg = parser.GetArg("fold-factor", 1).c_str();
    foldFactor[1] = (int)strtol (arg, &endptr, 10);
    assert (arg[0] != '\0' && (*endptr == '\0'));
    
    arg = parser.GetArg("fold-factor", 2).c_str();
    foldFactor[2] = (int)strtol (arg, &endptr, 10);
    assert (arg[0] != '\0' && (*endptr == '\0'));
    if (parser.Found("tile")) {
      cerr << "Cannot simultaneously tile and unfold.  Please choose one.\n";
      abort();
    }
  }
  if (parser.Found("tilemat")) 
    system.SetTileMatrix (string2Mat3 (parser.GetArg("tilemat",0)));
  else if (parser.Found("tile")) {
    unfold = true;

    Int3 tile = string2Int3 (parser.GetArg("tile",0));
    Mat3 tileMat;
    tileMat = 0.0;
    tileMat(0,0) = (double)tile[0];
    tileMat(1,1) = (double)tile[1];
    tileMat(2,2) = (double)tile[2];
    system.SetTileMatrix(tileMat);

    const char *arg;
    char *endptr;
    arg = parser.GetArg("tile", 0).c_str();
    tileFactor[0] = (int)strtol (arg, &endptr, 10);
    assert (arg[0] != '\0' && (*endptr == '\0'));
    
    arg = parser.GetArg("tile", 1).c_str();
    tileFactor[1] = (int)strtol (arg, &endptr, 10);
    assert (arg[0] != '\0' && (*endptr == '\0'));
    
    arg = parser.GetArg("tile", 2).c_str();
    tileFactor[2] = (int)strtol (arg, &endptr, 10);
    assert (arg[0] != '\0' && (*endptr == '\0'));
  }
  if (parser.Found("split")) {
    const char *arg;
    char *endptr;
    arg = parser.GetArg("split", 0).c_str();
    split[0] = (int)strtol (arg, &endptr, 10);
    assert (arg[0] != '\0' && (*endptr == '\0'));
    
    arg = parser.GetArg("split", 1).c_str();
    split[1] = (int)strtol (arg, &endptr, 10);
    assert (arg[0] != '\0' && (*endptr == '\0'));
    
    arg = parser.GetArg("split", 2).c_str();
    split[2] = (int)strtol (arg, &endptr, 10);
    assert (arg[0] != '\0' && (*endptr == '\0'));
  }

      
  if (system.ShiftOrbitals && !system.Localized) {
    perr << "Orbitals must be localized to be shifted.\n";
    abort();
  }
  if (system.OptimizeCenters && !system.Localized) {
    perr << "Cannot optimize localization centers without localizing.\n";
    abort();
  }
  if (system.Truncate && !system.ShiftOrbitals) {
    perr << "Orbitals must be shifted to be truncated.\n";
    abort();
  }
  

  system.SetFFTFactor (factor);
  if (!system.Read(inName)) {
    perr << "Could not read file " << inName << ".  Aborting.\n";
    abort();
  }

//   if (unfold) {
//     system.Unfold(0, foldFactor);
//     if (system.SpinPolarized)
//       system.Unfold(1, foldFactor);
//   }
//   if (parser.Found ("tile")) {
//     system.Tile (0, tileFactor);
//     if (system.SpinPolarized)
//       system.Tile (0, tileFactor);
//   }
      

  if (parser.Found("localize")) {
    string name = parser.GetArg("localize");
    clock_t start = clock();
    system.Localize(name, radius);
    clock_t end = clock();
    perr << "Total time = " << (end-start)/CLOCKS_PER_SEC << " seconds.\n";
  }

  if (parser.Found("skin-thickness")) {
    double thickness = atof (parser.GetArg("skin-thickness").c_str());
    perr << "skin thickness = " << thickness << endl;
    system.SetSkinThickness (thickness);
  }

  if (parser.Found ("bwfn")) {
    if (!parser.Found ("split"))
      system.WriteBWFN (parser.GetArg ("bwfn"));
    else
      system.WriteBWFN (parser.GetArg ("bwfn"), split);
  }
  if (parser.Found ("pwfn"))
    system.WritePWFN (parser.GetArg ("pwfn"));
  if (parser.Found("qmcPACK"))
    system.Write (parser.GetArg ("qmcPACK"));

  COMM::Finalize();
}
