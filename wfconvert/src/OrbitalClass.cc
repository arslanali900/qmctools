#include "OrbitalClass.h"
#include "CellClass.h"
#include "CenterClass.h"
#include "ParserClass.h"
#include <Common/IO/IO.h>
#include <Common/MPI/Communication.h>
#include <vector>
#include <algorithm>

using namespace IO;

void
OrbitalClass::Setk(Vec3 k)
{
  Twist = Cell->Lattice.k2Twist (k);
}

Vec3
OrbitalClass::Getk()
{
  return Cell->Lattice.Twist2k (Twist);
}

void
OrbitalClass::Read (IOSectionClass &in)
{

}

void
OrbitalClass::Write (IOSectionClass &out)
{


}

void
OrbitalClass::Broadcast(CommunicatorClass &comm, int root)
{
  comm.Broadcast (root, Twist);
  comm.Broadcast (root, Gshift);
  comm.Broadcast (root, Eigenvalue);
  comm.Broadcast (root, Norm2);
  comm.Broadcast (root, Spin);
  comm.Broadcast (root, TwistIndex);
  comm.Broadcast (root, BandIndex);
  if (Coefs == NULL)
    Coefs = new zVec;
  int numCoefs = Coefs->size();
  comm.Broadcast (root, numCoefs);
  if (Coefs->size() != numCoefs)
    Coefs->resize(numCoefs);
  comm.Broadcast (root, *Coefs);
}

void
OrbitalClass::Send (CommunicatorClass &comm, int toProc)
{
  comm.Send (toProc, Twist);
  comm.Send (toProc, Gshift);
  comm.Send (toProc, Eigenvalue);
  comm.Send (toProc, Norm2);
  comm.Send (toProc, Spin);
  comm.Send (toProc, TwistIndex);
  comm.Send (toProc, BandIndex);
  int numCoefs = Coefs->size();
  comm.Send (toProc, numCoefs);
  comm.Send (toProc, *Coefs);
}


void
OrbitalClass::Receive (CommunicatorClass &comm, int fromProc)
{
  comm.Receive (fromProc, Twist);
  comm.Receive (fromProc, Gshift);
  comm.Receive (fromProc, Eigenvalue);
  comm.Receive (fromProc, Norm2);
  comm.Receive (fromProc, Spin);
  comm.Receive (fromProc, TwistIndex);
  comm.Receive (fromProc, BandIndex);
  int numCoefs;
  comm.Receive (fromProc, numCoefs);
  if (Coefs == NULL)
    Coefs = new zVec;
  Coefs->resize(numCoefs);
  comm.Receive (fromProc, *Coefs);
}


void
OrbitalClass::WriteSpline (IOSectionClass &out, bool real,
			   bool shift, bool truncate)
{
  Cell->SetupFFT();
  FFTBox &FFT = Cell->FFT;
  

  out.WriteVar ("eigenvalue", Eigenvalue);
  out.WriteVar ("spin", Spin);
  out.WriteVar ("occupancy", Occupancy);
  
  if (Center == NULL || !shift) 
    PutInFFTBox();
  else {
    Array<complex<double>,1> shiftedCoefs(Coefs->size());
    for (int i=0; i<Coefs->size(); i++) {
      Vec3 G = FFT.GVecs(i);
      G += Gshift;
      double phase = -dot (FFT.GVecs(i), Center->r);
      complex<double> e2iGr (cos(phase), sin(phase));
      shiftedCoefs(i) = e2iGr * (*Coefs)(i);
      Int3 in = (*FFTIndices)(i);
      if ((in[0]+in[1]+in[2])%2==1)
	shiftedCoefs(i) *= -1.0;
    }
    PutInFFTBox(shiftedCoefs);
  }
  FFT.k2r();
  int nx, ny, nz;
  FFT.GetDims (nx, ny, nz);

  double wfnorm = 0.0;
  double vol = Cell->Lattice.GetVolume();
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

  if (Center != NULL) {
    Array<double,1> r(3);
    r(0) = Center->r[0];  r(1)=Center->r[1]; r(2)=Center->r[2];
    out.WriteVar ("center", r);
    out.WriteVar ("radius", Center->Radius);
  }
}


void
OrbitalClass::WriteBWFN (FILE *fout, bool real, bool shift,
			 bool truncate)
{


}


void
OrbitalClass::PutInFFTBox()
{
  FFTBox &FFT = Cell->FFT;
  FFT.kBox = complex<double>();
  for (int i=0; i<FFTIndices->size(); i++) 
    FFT.kBox((*FFTIndices)(i)) = (*Coefs)(i);
}

void
OrbitalClass::PutInFFTBox (zVec &coefs)
{
  FFTBox &FFT = Cell->FFT;
  FFT.kBox = complex<double>();
  for (int i=0; i<FFTIndices->size(); i++) 
    FFT.kBox((*FFTIndices)(i)) = coefs(i);
}

void
OrbitalClass::AddToFFTBox(complex<double> prefactor)
{
  FFTBox &FFT = Cell->FFT;
  for (int i=0; i<FFTIndices->size(); i++) 
    FFT.kBox((*FFTIndices)(i)) += prefactor*(*Coefs)(i);
}

void
OrbitalClass::AddToFFTBox(zVec &coefs)
{
  FFTBox &FFT = Cell->FFT;
  for (int i=0; i<FFTIndices->size(); i++) 
    FFT.kBox((*FFTIndices)(i)) += coefs(i);
}

bool
OrbitalClass::Read (ParserClass &parser, int numCoefs)
{
  double nrm = 0.0;
  Coefs = new Array<complex<double>,1> (numCoefs);
  // Just avoid pointer deferencing
  Array<complex<double>,1> &coefs = *Coefs;
  for (int i=0; i<numCoefs; i++) {
    if (!parser.ReadComplex(coefs(i))) return false;
    nrm += norm (coefs(i));
  }
  coefs *= sqrt(1.0/nrm);

  return true;
}

struct Vec3Less
{
  inline bool operator()(Vec3 a, Vec3 b)
  { return dot(a,a) < dot(b,b);   }
};


void
DensityClass::Set(const Array<double,3> &rhoTotal, double ecut)
{
  RhoUp.resize(rhoTotal.shape());
  RhoUp = rhoTotal;

  int n0 = RhoUp.extent(0);
  int n1 = RhoUp.extent(1);
  int n2 = RhoUp.extent(2);
  Vec3 b0 = Cell.Lattice.b(0);
  Vec3 b1 = Cell.Lattice.b(1);
  Vec3 b2 = Cell.Lattice.b(2);

  double Gcut = 2.0*sqrt(2.0*ecut);

  vector<Vec3> GVecs;
  // Now, setup the G-vector set
  for (int i0=-(n0/2+1); i0 < (n0/2+1); i0++)
    for (int i1=-(n1/2+1); i1 < (n1/2+1); i1++)
      for (int i2=-(n2/2+1); i2 < (n2/2+1); i2++) {
	Vec3 G = (double)i0 * b0 + (double)i1*b1 + (double)(i2)*b2;
	if (dot(G,G) < Gcut*Gcut) 
	  GVecs.push_back(G);
      }
  
  // Sort the G-vectors by magnitude
  Vec3Less comparator;
  sort (GVecs.begin(), GVecs.end(), comparator);

  int numG = GVecs.size();
  Array<Vec3,1> GVecsArray(numG);
  for (int iG=0; iG<numG; iG++)
    GVecsArray(iG) = GVecs[iG];

  // Setup the FFT box
  Cell.FFT.GVecs.Set (Cell.Lattice.GetDirect(), GVecsArray, RhoUp.shape());
  Cell.FFT.Setup();

  // Now compute Fouier transfroms
  UpCoefs.resize(numG);  
  for (int i0=0; i0<n0; i0++)
    for (int i1=0; i1<n1; i1++)
      for (int i2=0; i2<n2; i2++)
	Cell.FFT.rBox(i0,i1,i2) = rhoTotal(i0,i1,i2);
  Cell.FFT.r2k();
  Cell.FFT.GetkVec (UpCoefs);  
  // UpCoefs *= Cell.Lattice.GetVolume();
}

void
DensityClass::Set(const Array<double,3> &rhoUp,
		  const Array<double,3> &rhoDown,
		  double ecut)
{
  RhoUp.resize(rhoUp.shape());      RhoUp   = rhoUp;
  RhoDown.resize(rhoDown.shape());  RhoDown = rhoDown;

  int n0 = RhoUp.extent(0);
  int n1 = RhoUp.extent(1);
  int n2 = RhoUp.extent(2);
  Vec3 b0 = Cell.Lattice.b(0);
  Vec3 b1 = Cell.Lattice.b(1);
  Vec3 b2 = Cell.Lattice.b(2);

  double Gcut = 2.0*sqrt(2.0*ecut);

  vector<Vec3> GVecs;
  // Now, setup the G-vector set
  for (int i0=-(n0/2+1); i0 < (n0/2+1); i0++)
    for (int i1=-(n1/2+1); i1 < (n1/2+1); i1++)
      for (int i2=-(n2/2+1); i2 < (n2/2+1); i2++) {
	Vec3 G = (double)i0 * b0 + (double)i1*b1 + (double)(i2)*b2;
	if (dot(G,G) < Gcut*Gcut) 
	  GVecs.push_back(G);
      }
  
  Vec3Less comparator;
  sort (GVecs.begin(), GVecs.end(), comparator);

  int numG = GVecs.size();

  Array<Vec3,1> GVecsArray(numG);
  for (int iG=0; iG<numG; iG++)
    GVecsArray(iG) = GVecs[iG];

  // Setup the FFT box
  Cell.FFT.GVecs.Set (Cell.Lattice.GetDirect(), GVecsArray, rhoUp.shape());
  Cell.FFT.Setup();

  // Now compute Fouier transfroms
  UpCoefs.resize(GVecs.size());
  for (int i0=0; i0<n0; i0++)
    for (int i1=0; i1<n1; i1++)
      for (int i2=0; i2<n2; i2++)
	Cell.FFT.rBox(i0,i1,i2) = rhoUp(i0,i1,i2);
  Cell.FFT.r2k();
  Cell.FFT.GetkVec (UpCoefs);

  DownCoefs.resize(GVecs.size());
  for (int i0=0; i0<n0; i0++)
    for (int i1=0; i1<n1; i1++)
      for (int i2=0; i2<n2; i2++)
	Cell.FFT.rBox(i0,i1,i2) = rhoDown(i0,i1,i2);
  Cell.FFT.r2k();
  Cell.FFT.GetkVec (DownCoefs);
}


void
DensityClass::Set(const Array<Int3,1> &gvecs, const zVec &coefs)
{
  UpCoefs.resize(coefs.shape());        UpCoefs   =   coefs;

  Vec3 b0 = Cell.Lattice.b(0);
  Vec3 b1 = Cell.Lattice.b(1);
  Vec3 b2 = Cell.Lattice.b(2);

  Array<Vec3,1> GVecs(gvecs.size());
  Int3 imin(0,0,0), imax(0,0,0);
  for (int ig=0; ig<gvecs.size(); ig++) {
    // Flip sign since wfconv uses wrong sign convention internally
    Int3 g = -1*gvecs(ig);
    for (int j=0; j<3; j++) {
      imin[j] = min(imin[j],g[j]);
      imax[j] = max(imax[j],g[j]);
    }
    GVecs(ig) = (double)g[0]*b0 + (double)g[1]*b1 + (double)g[2]*b2;
  }
  Int3 mesh(imax - imin);
  int numG = GVecs.size();

  // Setup the FFT box
  Cell.FFT.GVecs.Set (Cell.Lattice.GetDirect(), GVecs, mesh);
  Cell.FFT.Setup();

  Cell.FFT.PutkVec(UpCoefs);
  Cell.FFT.k2r();
  RhoUp.resize(mesh);
  for (int i0=0; i0<mesh[0]; i0++)
    for (int i1=0; i1<mesh[1]; i1++)
      for (int i2=0; i2<mesh[2]; i2++)
	RhoUp(i0,i1,i2) = Cell.FFT.rBox(i0,i1,i2).real();
}


void
DensityClass::Set(const Array<Int3,1> &gvecs,
		  const zVec &upCoefs, const zVec &downCoefs)
{
  UpCoefs.resize(upCoefs.shape());        UpCoefs   =   upCoefs;
  DownCoefs.resize(downCoefs.shape());  DownCoefs   = downCoefs;

  Vec3 b0 = Cell.Lattice.b(0);
  Vec3 b1 = Cell.Lattice.b(1);
  Vec3 b2 = Cell.Lattice.b(2);

  Array<Vec3,1> GVecs(gvecs.size());
  Int3 imin(0,0,0), imax(0,0,0);
  for (int ig=0; ig<gvecs.size(); ig++) {
    // Flip sign since wfconv uses wrong sign convention internally
    Int3 g = -1*gvecs(ig);
    for (int j=0; j<3; j++) {
      imin[j] = min(imin[j],g[j]);
      imax[j] = max(imax[j],g[j]);
    }
    GVecs(ig) = (double)g[0]*b0 + (double)g[1]*b1 + (double)g[2]*b2;
  }
  Int3 mesh(imax - imin);
  int numG = GVecs.size();

  // Setup the FFT box
  Cell.FFT.GVecs.Set (Cell.Lattice.GetDirect(), GVecs, mesh);
  Cell.FFT.Setup();

  Cell.FFT.PutkVec(UpCoefs);
  Cell.FFT.k2r();
  RhoUp.resize(mesh);
  for (int i0=0; i0<mesh[0]; i0++)
    for (int i1=0; i1<mesh[1]; i1++)
      for (int i2=0; i2<mesh[2]; i2++)
	RhoUp(i0,i1,i2) = Cell.FFT.rBox(i0,i1,i2).real();
  
  Cell.FFT.PutkVec(DownCoefs);
  Cell.FFT.k2r();
  RhoDown.resize(mesh);
  for (int i0=0; i0<mesh[0]; i0++)
    for (int i1=0; i1<mesh[1]; i1++)
      for (int i2=0; i2<mesh[2]; i2++)
	RhoDown(i0,i1,i2) = Cell.FFT.rBox(i0,i1,i2).real();
}


void
DensityClass::Write(IOSectionClass &out) 
{
  int numSpins;
  if (RhoUp.size() && RhoDown.size())  numSpins = 2;
  else if (RhoUp.size())               numSpins = 1;
  else                                 numSpins = 0;

  Array<int,2> gvecs(Cell.FFT.GVecs.size(),3);
  double twoPiInv = 0.5/M_PI;
  for (int i=0; i<gvecs.extent(0); i++)
    for (int j=0; j<3; j++)
      gvecs(i,j) = (int)round(twoPiInv * dot(Cell.FFT.GVecs(i), Cell.Lattice.a(j)));


  if (numSpins) {
    out.NewSection("density");
    out.WriteVar ("reduced_gvecs", gvecs);
    Array<double,2> RealCoefs(UpCoefs.size(),2);
    for (int iG=0; iG<UpCoefs.size(); iG++) {
      RealCoefs(iG,0) = real(UpCoefs(iG));
      RealCoefs(iG,1) = imag(UpCoefs(iG));
    }

    if (numSpins == 1) {
      out.WriteVar ("rho_r", RhoUp);
      out.WriteVar ("rho_G", RealCoefs);
    }
    else if (numSpins == 2) {
      out.NewSection("spin");
      out.WriteVar ("rho_r", RhoUp);
      out.WriteVar ("rho_G", RealCoefs);
      out.CloseSection();
      for (int iG=0; iG<UpCoefs.size(); iG++) {
	RealCoefs(iG,0) = real(DownCoefs(iG));
	RealCoefs(iG,1) = imag(DownCoefs(iG));
      }
      out.NewSection("spin");
      out.WriteVar ("rho_r", RhoDown);
      out.WriteVar ("rho_G", RealCoefs);
      out.CloseSection();
    }
    out.CloseSection(); // "density"
  }
}


void
DensityClass::WriteESHDF(IOSectionClass &out, string basename) 
{
  //  cerr << "Writing \"" << basename << "\" section.\n";
  int numSpins;
  if (RhoUp.size() && RhoDown.size())  numSpins = 2;
  else if (RhoUp.size())               numSpins = 1;
  else                                 numSpins = 0;
  cerr << "numSpins =  " << numSpins << endl;


  Array<int,2> gvecs(Cell.FFT.GVecs.size(),3);
  double twoPiInv = 0.5/M_PI;
  // Note:  ESHDF uses opposite sign convention than wfconv does
  // internally, so we swap the sign of the gvecs.
  for (int i=0; i<gvecs.extent(0); i++)
    for (int j=0; j<3; j++)
      gvecs(i,j) = -(int)round(twoPiInv * dot(Cell.FFT.GVecs(i), Cell.Lattice.a(j)));


  if (numSpins) {
    out.NewSection(basename);
    out.WriteVar ("number_of_gvectors", gvecs.size());
    out.WriteVar ("gvectors", gvecs);

    Array<double,1> mesh(3);
    for (int i=0; i<3; i++) mesh(i) = RhoUp.extent(i);
    out.WriteVar ("mesh", mesh);
    Array<double,2> RealCoefs(UpCoefs.size(),2);
    for (int iG=0; iG<UpCoefs.size(); iG++) {
      RealCoefs(iG,0) = real(UpCoefs(iG));
      RealCoefs(iG,1) = imag(UpCoefs(iG));
    }
    out.SetUnderscores(true);
    out.NewSection("spin");
    out.SetUnderscores(false);
//     out.WriteVar ("density_r", RhoUp);
//     out.WriteVar ("density_g", RealCoefs);
    out.WriteVar(basename + "_r", RhoUp);
    out.WriteVar(basename + "_g", RealCoefs);
    out.CloseSection();
    if (numSpins > 1) {
      for (int iG=0; iG<UpCoefs.size(); iG++) {
	RealCoefs(iG,0) = real(DownCoefs(iG));
	RealCoefs(iG,1) = imag(DownCoefs(iG));
      }
      out.NewSection("spin");
//       out.WriteVar ("density_r", RhoDown);
//       out.WriteVar ("density_g", RealCoefs);
      out.WriteVar(basename + "_r", RhoDown);
      out.WriteVar(basename + "_g", RealCoefs);
      out.CloseSection();
    }
    out.CloseSection(); // basename
  }
}


