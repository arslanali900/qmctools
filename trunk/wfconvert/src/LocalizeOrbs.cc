#include "OrbitalSetClass.h"
#include "ParserClass.h"
#include "ParseCommand.h"
#include <Common/IO/IO.h>
#include <Common/MatrixOps/MatrixOps.h>
//#include <sys/sysinfo.h>
#include <sys/types.h>
#include <unistd.h>
#include <time.h>
#include <sstream>


long FreeMemory()
{
  long memFree;
  MemParserClass memParser;
  if (memParser.OpenFile ("/proc/meminfo")) {
    assert(memParser.FindToken ("MemFree:"));
    memParser.ReadLong (memFree);
    long buffers, inactive, cached;
    assert(memParser.FindToken("Buffers:"));
    memParser.ReadLong(buffers);
    assert(memParser.FindToken("Cached:"));
    memParser.ReadLong(cached);
    memFree += cached;
    assert(memParser.FindToken("Inactive:"));
    memParser.ReadLong(inactive);
    memFree += inactive/2;

    memFree *= 1024;
    memParser.CloseFile();
    return memFree;
  }
  else // Default to 256MB
    return 268435456;
}




////////////////////////////////////////////////////////////
//                Localization functions                  //
////////////////////////////////////////////////////////////

inline string 
Extension (string fileName)
{
  string extn;
  stack<char> bwExtn;
  int pos = fileName.length()-1;
  while ((pos >= 0) && fileName[pos]!='.') {
    bwExtn.push(fileName[pos]);
    pos--;
  }
  
  if (fileName[pos] == '.') 
    while (!bwExtn.empty()) {
      extn += bwExtn.top();
      bwExtn.pop();
    }
  else
    extn = "";
  return (extn);
}

bool 
OrbitalSetClass::ReadCenters(string fname) 
{
  if (Extension(fname) == "in") {
    IO::IOSectionClass in;
    assert (in.OpenFile (fname));
    return ReadCenters(in);
  }

  CellClass &cell = WriteSuper ? SuperCell : PrimCell;

  MemParserClass parser;
  int numCenters, locType;
  if (!parser.OpenFile (fname))                 return false;
  if (!parser.FindToken ("Number of centres"))  return false;
  if (!parser.ReadInt (numCenters))             return false;
  if (!parser.FindToken ("regions"))            return false;
  if (!parser.ReadInt (locType))                return false;
  if (!parser.FindToken ("(up and down)"))      return false;
  if (!parser.FindToken ("\n"))                 return false;
  if (!parser.ReadInt (UseUpBands))             return false;
  if (!parser.ReadInt (UseDownBands))           return false;
  if (!parser.FindToken ("\n"))                 return false;
  if (!parser.FindToken ("\n"))                 return false;
  Centers.resize(numCenters);
  for (int i=0; i<numCenters; i++) {
    if (!parser.ReadDouble(Centers[i].r[0]))    return false;
    if (!parser.ReadDouble(Centers[i].r[1]))    return false;
    if (!parser.ReadDouble(Centers[i].r[2]))    return false;
    if (!parser.ReadDouble(Centers[i].Radius))  return false;
    if (!parser.ReadInt   (Centers[i].NumUp))   return false;
    if (!parser.ReadInt   (Centers[i].NumDown)) return false;
    Centers[i].NumOrbitals = 
      max(Centers[i].NumUp, Centers[i].NumDown);
    Centers[i].SetupBitfield(cell);
  }
  if (parser.FindToken ("skin") && parser.FindToken ("\n"))
    parser.ReadDouble (SkinThickness);
  else
    SkinThickness = 0.0;
  for (int i=0; i<numCenters; i++)
    Centers[i].SkinThickness = SkinThickness;
    
  parser.CloseFile();
  return true;
}


bool
OrbitalSetClass::ReadCenters(IO::IOSectionClass &in)
{
  CellClass &cell = WriteSuper ? SuperCell : PrimCell;

  assert (in.ReadVar ("UseUpBands", UseUpBands));
  assert (in.ReadVar ("UseDownBands", UseDownBands));
  in.ReadVar ("NumExtendedBands", NumExtendedBands);
  int numSections = in.CountSections ("LocalCenter");
  // This is a new format, allowing identical orbitals to be cloned at
  // different sites without extra storage
  if (numSections > 0) {
    Centers.resize(numSections);
    int upTotal=0, downTotal=0;
    for (int i=0; i<numSections; i++) {
      in.OpenSection ("LocalCenter", i);
      Array<double,2> reducedCenters, reflections;
      assert (in.ReadVar ("ReducedCenters", reducedCenters));
      if (in.ReadVar ("Reflections", reflections))
	assert (reducedCenters.shape() == reflections.shape());
      Centers[i].IdenticalSites.resize(reducedCenters.extent(0));
      if (reflections.extent(0) > 0)
	Centers[i].Reflections.resize(reflections.extent(0));
      for (int j=0; j<reducedCenters.extent(0); j++) {
	Vec3 u(reducedCenters(j,0), reducedCenters(j,1), reducedCenters(j,2));
	Centers[i].IdenticalSites(j) = cell.Lattice.u2r (u);
	if (j < reflections.extent(0))
	  Centers[i].Reflections(j) = 
	    Vec3 (reflections(j,0), reflections(j,1), reflections(j,2));
      }
      Centers[i].r = Centers[i].IdenticalSites(0);
      assert (in.ReadVar ("Radius", Centers[i].Radius));
      if (!in.ReadVar ("SkinThickness", Centers[i].SkinThickness))
	Centers[i].SkinThickness = 0.0;
      assert (in.ReadVar ("NumUp", Centers[i].NumUp));
      assert (in.ReadVar ("NumDown", Centers[i].NumDown));
      Centers[i].NumOrbitals = max (Centers[i].NumUp, Centers[i].NumDown);
      upTotal   += Centers[i].NumUp   * reducedCenters.extent(0);
      downTotal += Centers[i].NumDown * reducedCenters.extent(0);
      in.CloseSection(); // "LocalCenter"
    }
    if (upTotal != UseUpBands) 
      perr << "Warning:  UseUpBands does not match the total number of "
	   << "localized up orbitals.\n"
	   << "UseUpBands = " << UseUpBands << "  total up orbitals = " 
	   << upTotal << ".\n";
    if (downTotal != UseDownBands) 
      perr << "Warning:  UseDownBands does not match the total number of "
	   << "localized down orbitals.\n"
	   << "UseDownBands = " << UseDownBands << "  total down orbitals = " 
	   << downTotal << ".\n";

  }
  else {  
    Array<double,2> reducedCenters;
    assert (in.ReadVar ("ReducedCenters", reducedCenters));
    Array<double,1> radii;
    assert (in.ReadVar ("Radii", radii));
    Array<int,2> numOrbitals;
    assert (in.ReadVar ("NumOrbitals", numOrbitals));
    // Check that sizes match
    int N = reducedCenters.extent(0);
    assert (reducedCenters.extent(1) == 3);
    assert (radii.size() == N);
    assert (numOrbitals.extent(0) == N);
    assert (numOrbitals.extent(1) == 2);
    
    double skinThickness (0.0);
    in.ReadVar("SkinThickness", skinThickness);
    // Now create centers
    Centers.resize (N);
    for (int i=0; i<N; i++) {
      Vec3 u(reducedCenters(i,0), reducedCenters(i,1), reducedCenters(i,2));
      Centers[i].r = cell.Lattice.u2r(u);
      Centers[i].IdenticalSites.resize(1);
      Centers[i].IdenticalSites(0) = Centers[i].r;
      Centers[i].Radius  = radii(i);
      Centers[i].NumUp   = numOrbitals(i,0);
      Centers[i].NumDown = numOrbitals(i,1);
      Centers[i].NumOrbitals = max (numOrbitals(i,0), numOrbitals(i,1));
      Centers[i].SkinThickness = skinThickness;
    }
  }
    return true;
}


// The following code roughly implements the nonorthonal localization
// methods of Reboredo and Williamson, PRB 71, 12115(R) (2005)
void
OrbitalSetClass::Localize (string centersFilename, double radius)
{
  CellClass &cell = WriteSuper ? SuperCell : PrimCell;
  cell.SetupFFT();
  FFTBox &FFT = cell.FFT;
  if (Comm.MyProc() == 0) {
    cell.SetupFFT();
    if (centersFilename == "random") {
      int numElecs = NumElectrons[0] + NumElectrons[1];
      if (WriteSuper)
	numElecs *= (int)round(fabs(det(TileMatrix)));      
      int numOcc = numElecs/2;
      UseUpBands = UseDownBands = numOcc;

      Centers.resize(numOcc);
      
      Vec3 r0 = cell.Lattice.a(0);
      Vec3 r1 = cell.Lattice.a(1);
      Vec3 r2 = cell.Lattice.a(2);
      for (int i=0; i<numOcc; i++) {
	Vec3 r = drand48()*r0 + drand48()*r1 + drand48()*r2;
	Centers[i] = CenterClass (r, radius, 1);
	Centers[i].NumUp = Centers[i].NumDown = 1;
	Centers[i].SetupBitfield (cell);
      }
    }
    else
      ReadCenters (centersFilename);
  }
  
  if (Comm.NumProcs() > 0) {
    DistributeBands();
    int numSpins = PrimOrbitals[1].size() > 0 ? 2 : 1;
    for (int spin=0; spin < numSpins; spin++) {
      if (OptimizeCenters)
	do {
	  LocalizeMPI(true, spin);
	} while (!UpdateCenters(spin));
      
      LocalizeMPI(Orthogonalize, spin);
    }
  }
  else {
    // Now do the localization  
    int numSpins = PrimOrbitals[1].size() > 0 ? 1 : 2;
    for (int spin=0; spin<numSpins; spin++) {
      if (OptimizeCenters)
	do {
	  Localize(true, spin);
	} while (!UpdateCenters(spin));
      
      Localize(Orthogonalize, spin);
    }
  }
  Localized = true;
}

bool 
OrbitalSetClass::UpdateCenters(int spin)
{
  CellClass &cell = WriteSuper ? SuperCell : PrimCell;
  TinyVector<Array<OrbitalClass*,2>,2> &orbs =
    WriteSuper ? SuperOrbitals : PrimOrbitals;
  FFTBox &FFT = cell.FFT;

  int nx, ny, nz;
  FFT.GetDims(nx,ny,nz);
  double nxInv=1.0/(double)nx;
  double nyInv=1.0/(double)ny;
  double nzInv=1.0/(double)nz;
  Vec3 r0 = cell.Lattice.a(0);
  Vec3 r1 = cell.Lattice.a(1);
  Vec3 r2 = cell.Lattice.a(2);
  int numk = orbs[spin].extent(0);
  bool converged = true;
  double drmax2 = 0.001*
    (nxInv*nxInv*dot(r0,r0) + nyInv*nyInv*dot(r1,r1) + nzInv*nzInv*dot(r2,r2));

  Array<Vec3,1> drAvg(Centers.size()), drSum(Centers.size());
  drAvg = Vec3(0.0, 0.0, 0.0);
  int firstCenter, lastCenter;
  ProcTasks (Comm.MyProc(), Centers.size(), firstCenter, lastCenter);
  int orb=0;
  for (int ci=0; ci<firstCenter; ci++)
    orb += Centers[ci].NumOrbitals;
  for (int ci=firstCenter; ci<=lastCenter; ci++) {
    CenterClass &center = Centers[ci];
    double nrm = 0.0;
    double prefact = 
      1.0/((double)center.NumOrbitals*(double)(nx*ny*nz));
    for (int iorb=0; iorb<center.NumOrbitals; iorb++) {
      int band = orb/numk;
      int k    = orb%numk;
      orbs[spin](k, band)->PutInFFTBox();
      FFT.k2r();
      double sx, sy, sz;
      for (int ix=0; ix<nx; ix++) {
	sx = (double)ix * nxInv;
	for (int iy=0; iy<ny; iy++) {
	  sy = (double)iy * nyInv;
	  for (int iz=0; iz<nz; iz++) {
	    sz = (double) iz * nzInv;
	    Vec3 r = sx*r0 + sy*r1 + sz*r2;
	    Vec3 diff = cell.Lattice.MinImage(r - center.r);
	    double rho = norm(FFT.rBox(ix,iy,iz));
	    nrm += prefact * rho;
	    drAvg(ci) += prefact*rho*diff;
	  }
	}
      }
      orb++;
    }
    drAvg(ci) *= 1.0/nrm;
  }
  // Communicate the results across all processors
  Comm.AllSum (drAvg, drSum);
  for (int ci=0; ci<Centers.size(); ci++) {
    if (dot(drSum(ci), drSum(ci)) > drmax2)
      converged = false;
    if (Comm.MyProc() == 0)
      fprintf (stderr, "drAvg = [ %11.6f %11.6f %11.6f ]\n", 
	       drSum(ci)[0], drSum(ci)[1], drSum(ci)[2]);
    Centers[ci].r += drSum(ci);
    Centers[ci].SetupBitfield(cell);
  }

  return converged;
}



inline complex<double> promote (complex<float> val)
{ return complex<double>(val.real(), val.imag()); }

inline complex<double> promote (complex<double> val)
{ return val; }


int
OrbitalSetClass::NearestIon(Vec3 pos)
{
  CellClass &cell = WriteSuper ? SuperCell : PrimCell;

  if (cell.IonPos.size() < 1) 
    return -1;
  Vec3 disp;
  double minDist;
  int minIon;

  disp    = cell.Lattice.MinImage (pos-cell.IonPos(0));
  minDist = sqrt (dot(disp, disp));
  minIon  = 0;
  for (int i=0; i<cell.IonPos.size(); i++) {
    disp = cell.Lattice.MinImage (pos-cell.IonPos(i));
    double dist = sqrt (dot (disp, disp));
    if (dist < minDist) {
      minIon = i;
      minDist = dist;
    }
  }
  return minIon;
}
    



void
OrbitalSetClass::Localize(bool ortho, int spin)
{
  CellClass &cell = WriteSuper ? SuperCell : PrimCell;
  TinyVector<Array<OrbitalClass*,2>,2> &orbs =
    WriteSuper ? SuperOrbitals : PrimOrbitals;
  FFTBox &FFT = cell.FFT;

  int N = orbs[spin](0,0)->GetCoefs().size();
  int numk = orbs[spin].extent(0);
  cerr << "numk = " << numk << endl;
  int numOcc = (spin==0) ? UseUpBands : UseDownBands;
  cerr << "Using " << numOcc << " bands for "
       << ((spin==0) ? "up" : "down") << "-spin electrons.\n";
  int numOrbs = 0;
  vector<CenterClass>::iterator iter;
  for (iter = Centers.begin(); iter != Centers.end(); iter++) 
    numOrbs += (spin==0) ? iter->NumUp : iter->NumDown;
  if (numOrbs != numOcc) {
    cerr << "Warning:  Number of localized orbitals should equal "
	 << "the number of occupied bands.\n";
    cerr << "NumOrbs = " << numOrbs << ",  numOcc = " << numOcc << endl;
    // exit(1);
  }

  // This tmp file stores the parts of the Phi matrices for each
  // center as we construct them.  This way, we don't have to store
  // all the matrices in RAM as we construct them.  The first index is
  // the center number.  The second two indices are for the 
  // numOcc x numOcc matrix of the band overlaps in the localization
  // region.  
  IO::IOSectionClass tmp;
  tmp.NewFile ("ThetaMats2.h5");
  Array<complex<double>,3> A(1, numOcc, numOcc);
  tmp.WriteVar("Theta", A);
  IO::IOVarBase &thetaVar = *(tmp.GetVarPtr("Theta"));
  thetaVar.Resize(Centers.size());

  // Determine the number of FFT boxes we can fit in ram
  int nx, ny, nz;
  FFT.GetDims(nx,ny,nz);
  double nrm = 1.0/(double)(nx*ny*nz);
  long fftSize = nx*ny*nz*sizeof(complex<double>);
  int numFFTs = max(1, (int)(FreeMemory() / fftSize));
  numFFTs = min (numFFTs, numOcc);
  cerr << "Using " << numFFTs << " cached FFT boxes.\n";

  // First, compute theta matrices and store them in the temporary
  // HDF5 file
  Array<complex<double>,4> phi_rows(numFFTs, nx, ny, nz);
  Array<complex<double>,2> thetaCols (Centers.size(), numOcc);
  thetaCols = complex<double>();
  int row = 0;
  Array<complex<double>,1> theta(Centers.size());

  cerr << "FFT.rBox.shape = " << FFT.rBox.shape() << endl;

  for (int ki=0; ki < numk; ki++) {
    while (row < numOcc) {
      cerr << "row = " << row << " of " << numOcc << endl;
      // How many rows to do in this iteration
      int todo = min (numFFTs, numOcc-row);
      for (int ri=0; ri<todo; ri++) {
	int band1 = (row+ri) / numk;
	int k1    = (row+ri) % numk;
	orbs[spin](k1,band1)->PutInFFTBox();
	FFT.k2r();
	for (int ix=0; ix<nx; ix++)
	  for (int iy=0; iy<ny; iy++)
	    for (int iz=0; iz<nz; iz++)
	      phi_rows(ri, ix, iy, iz) = conj (FFT.rBox(ix,iy,iz));
      }
      // Now loop through columns
      // First do part already in memory:
      for (int rowi=0; rowi<todo; rowi++)
	for (int coli=rowi; coli<todo; coli++) {
	  // cerr << " (rowi, coli) = (" << rowi << ", " << coli << ")\n";
	  int ci = 0;
	  theta = complex<double>();
	  for (iter = Centers.begin(); iter != Centers.end(); iter++) {
	    CenterClass &center = (*iter);
	    complex<double>* phi_row0 = &(phi_rows(rowi,0,0,0));
	    complex<double>* phi_col0 = &(phi_rows(coli,0,0,0));
	    for (int j=0; j<center.Offsets.size(); j++) {
	      int offset = center.Offsets[j];
	      theta(ci) += nrm * phi_row0[offset] * conj(phi_col0[offset]);
// 	      Int3 index = center.Indices[j];
// 	      int ix=index[0]; int iy=index[1]; int iz=index[2];
// 	      theta(ci) += nrm*(phi_rows(rowi, ix,iy,iz) * conj(phi_rows(coli, ix,iy,iz)));
	    }
// 	    for (int ix=0; ix<nx; ix++)
// 	      for (int iy=0; iy<ny; iy++)
// 		for (int iz=0; iz<nz; iz++) 
// 		  theta(ci) += nrm*center.Bitfield.GetDouble (ix,iy,iz)*
// 		    (phi_rows(rowi, ix,iy,iz) * conj(phi_rows(coli, ix,iy,iz)));
	    ci++;
	  }
	  // cerr << "theta = " << theta << endl;
	  thetaVar.Write(theta, Range::all(), row+rowi, row+coli);
	  for (int i=0; i<Centers.size(); i++)
	    theta(i) = conj(theta(i));
	  thetaVar.Write(theta, Range::all(), row+coli, row+rowi);
	  tmp.FlushFile();
	}
      // Now do the rest of the columns
      for (int col=row+todo; col < numOcc; col++) {
	// FFT orbital into real space
	int band = col / orbs[spin].extent(0);
	int k    = col % orbs[spin].extent(0);
	orbs[spin](k, band)->PutInFFTBox();
	FFT.k2r();
	for (int rowi=0; rowi < todo; rowi++) {
	  theta = complex<double>();
	  int ci=0;
	  for (iter = Centers.begin(); iter != Centers.end(); iter++) {
	    CenterClass &center = (*iter);
	    for (int j=0; j<center.Offsets.size(); j++) {
	      Int3 index = center.Indices[j];
	      int ix=index[0]; int iy=index[1]; int iz=index[2];
	      theta(ci) += nrm*(phi_rows(rowi,ix,iy,iz) * promote(FFT.rBox(ix,iy,iz)));
	    }
// 	    for (int ix=0; ix<nx; ix++)
// 	      for (int iy=0; iy<ny; iy++)
// 		for (int iz=0; iz<nz; iz++) 
// 		  theta(ci) += nrm*center.Bitfield.GetDouble (ix,iy,iz)*
// 		    (phi_rows(rowi,ix,iy,iz) * promote(FFT.rBox(ix,iy,iz)));
	    ci++;
	  }
	  thetaVar.Write(theta, Range::all(), row+rowi, col);
	  for (int i=0; i<Centers.size(); i++)
	    theta(i) = conj(theta(i));
	  thetaVar.Write(theta, Range::all(), col, row+rowi);
	}
      }
      row += numFFTs;
    }
    // Now we have the theta matrices stored.  We read them in one by
    // one and then diagonalize them.  The largest eigenvalues and
    // vectors are stored for each center.
    NewOrbCoefs[spin].resize(numk, numOcc, numOcc);
    NewOrbNorms[spin].resize(numk, numOcc);
    Array<complex<double>,2> thetaMat(numOcc, numOcc);
    Array<double,1> eigvals;
    Array<complex<double>,2> eigvecs;
    
    int orb = 0;
    int ci = 0;
    fprintf (stderr, "          Center             Norm     Nearest Ion\n");
    for (iter = Centers.begin(); iter != Centers.end(); iter++) {
      CenterClass &center = (*iter);
      thetaVar.Read(thetaMat, ci, Range::all(), Range::all());
      thetaMat = -1.0*thetaMat;
      SymmEigenPairs(thetaMat, center.NumOrbitals, eigvals, eigvecs);
      int numOrbs = (spin==0) ? center.NumUp : center.NumDown;
      for (int i=0; i<numOrbs; i++) {
	NewOrbNorms[spin](ki, orb) = -eigvals(i);
	NewOrbCoefs[spin](ki, orb,Range::all()) = eigvecs(i,Range::all());
	fprintf (stderr, "(%7.4f %7.4f %7.4f) %10.7f       %d\n", 
		 center.r[0], center.r[1], center.r[2],
		 NewOrbNorms[spin](ki,orb), NearestIon(center.r));
	// Store center in each orbital
	orbs[spin](ki, orb)->SetCenter(center);
	orb++;
      }
      ci++;
    }
    // Fill rest of coefs with identity
    for (int i=orb; i<numOcc; i++) 
      for (int j=0; j<numOcc; j++)
	NewOrbCoefs[spin](ki,i,j) = 0.0;
    for (int i=orb; i<numOcc; i++) 
      NewOrbCoefs[spin](ki, i,i) = 1.0;
    
    // Orthogonalize
    if (ortho) {
      // PolarOrthogonalize orthogonalizes the columns, so we transpose
      // first.  
      Array<complex<double>,2> orthoMat;
      orthoMat.reference(NewOrbCoefs[spin](ki,Range::all(),Range::all()));
      Transpose(orthoMat);
      PolarOrthogonalize (orthoMat);
      Transpose(orthoMat);
    } 
  }
  
  tmp.CloseFile();
}


void
OrbitalSetClass::DistributeBands()
{
  clock_t start, end;
  start = clock();
  Comm.Broadcast (0, WriteSuper);
  CellClass &cell = WriteSuper ? SuperCell : PrimCell;
  TinyVector<Array<OrbitalClass*,2>,2> &orbs =
    WriteSuper ? SuperOrbitals : PrimOrbitals;
  FFTBox &FFT = cell.FFT;

  // First, distribute all the essential data
  // Distribute the cells
  PrimCell.Broadcast  (Comm, 0);
  SuperCell.Broadcast (Comm, 0);

  Comm.Broadcast (0, TileMatrix);
  Comm.Broadcast (0, NumElectrons[0]);
  Comm.Broadcast (0, NumElectrons[1]);
  Comm.Broadcast (0, UseUpBands);
  Comm.Broadcast (0, UseDownBands);

  // Distribute the FFTIndices
  int numPrimIndices = PrimFFTIndices.size();
  Comm.Broadcast (0, numPrimIndices);
  int numk, numSuperIndices;
  numk = SuperFFTIndices.size();
  if (Comm.MyProc() == 0)
    numSuperIndices = SuperFFTIndices[0].size();
  Comm.Broadcast(0, numk);
  Comm.Broadcast(0, numSuperIndices);
  if (Comm.MyProc() != 0) {
    PrimFFTIndices.resize(numPrimIndices);
    SuperFFTIndices.resize(numk);
    for (int ik=0; ik<numk; ik++)
      SuperFFTIndices[ik].resize(numSuperIndices);
  }
  Comm.Broadcast (0, PrimFFTIndices);
  for (int ik=0; ik<numk; ik++) 
    Comm.Broadcast (0, SuperFFTIndices[ik]);

  // Send orbitals
  for (int spin=0; spin<2; spin++) {
    TileMap[spin].Broadcast (Comm, 0);
    int numk     = orbs[spin].extent(0);
    int numBands = orbs[spin].extent(1);
    Comm.Broadcast (0, numk);
    Comm.Broadcast (0, numBands);    

    if (numBands > 0) {
      if (Comm.MyProc() != 0)
	// Create orbitals if I'm not proc 0
	if ((orbs[spin].extent(0) != numk) ||
	    (orbs[spin].extent(1) != numBands)) {
	  orbs[spin].resize(numk, numBands);
	  for (int ik=0; ik<numk; ik++) 
	    for (int iband=0; iband<numBands; iband++) {
	      int primk, primBand;
	      if (WriteSuper)
		TileMap[spin].Super2Prim (ik, iband, primk, primBand);
	      else {
		primk    = ik;  
		primBand = iband; 
	      }
	      orbs[spin](ik,iband) = new OrbitalClass;
	      orbs[spin](ik,iband)->SetCell(cell);
	      if (WriteSuper) 
		orbs[spin](ik,iband)->SetFFTIndices (SuperFFTIndices[primk]);
	      else
		orbs[spin](ik,iband)->SetFFTIndices (PrimFFTIndices);
	    }
	}

      int num = (spin==0) ? UseUpBands : UseDownBands;
      ProcBands (Comm.MyProc(), num, MyFirstBand[spin], MyLastBand[spin]);

      // We now send each orbital only to the processor responsible
      // for it.
      if (Comm.MyProc() != 0) {
	for (int ik=0; ik<numk; ik++)
	  for (int iband=MyFirstBand[spin]; iband <= MyLastBand[spin];
	       iband++)
	    orbs[spin](ik,iband)->Receive(Comm, 0);
      }
      else {
	for (int proc=1; proc<Comm.NumProcs(); proc++) {
	  int first, last;
	  ProcBands (proc, num, first, last);
	  for (int ik=0; ik<numk; ik++)
	    for (int iband=first; iband <= last; iband++)
	      orbs[spin](ik,iband)->Send(Comm, proc);
	}
      }
       

//       // Now, broadcast eigenvalues and eigenvectors
//       for (int ik=0; ik<numk; ik++)
// 	for (int iband=0; iband<numBands; iband++) 
// 	  orbs[spin](ik,iband)->Broadcast (Comm, 0);
    }
  }
  // Send the centers
  int numCenters = Centers.size();
  Comm.Broadcast (0, numCenters);
  if (Comm.MyProc() != 0) 
    Centers.resize(numCenters);
  for (int i=0; i<numCenters; i++) 
    Centers[i].Broadcast(Comm, 0);
  end = clock();
  perr << "Time for DistributeBands() = " 
       << (double)(end-start)/(double)CLOCKS_PER_SEC << "s.\n";
}

inline complex<double> myconj (complex<double> z)
{
  return complex<double> (z.real(), -z.imag());
}


void
OrbitalSetClass::LocalizeMPI(bool ortho, int spin)
{
  CellClass &cell = WriteSuper ? SuperCell : PrimCell;
  TinyVector<Array<OrbitalClass*,2>,2> &orbs =
    WriteSuper ? SuperOrbitals : PrimOrbitals;
  FFTBox &FFT = cell.FFT;

  Int3 boxSize;
  cell.SetupFFT();
  FFT.GVecs.GetFFTBoxSize(boxSize[0], boxSize[1], boxSize[2]);
  for (int i=0; i<Centers.size(); i++)
    Centers[i].SetupBitfield(cell);

  //int N        = orbs[spin](0,0)->GetCoefs().size();
  int N = orbs[spin](0,MyFirstBand[spin])->GetCoefs().size();
  int numk     = orbs[spin].extent(0);
  int numBands = orbs[spin].extent(1);

  int numOcc = (spin==0) ? UseUpBands : UseDownBands;
  perr << "Using " << numOcc << " bands in localizing spin "
       << ((spin==0) ? "up" : "down") << " orbitals.\n";

  int numOrbs = 0;
  int elecs = 0;
  vector<CenterClass>::iterator iter;
  for (iter = Centers.begin(); iter != Centers.end(); iter++) {
    elecs += (iter->NumUp + iter->NumDown)*iter->IdenticalSites.size();
    numOrbs += ((spin == 0) ? iter->NumUp : iter->NumDown) * 
      iter->IdenticalSites.size();
  }
  if (elecs != NumElectrons[spin]) {
    perr << "Error:  Number of localized orbitals should equal "
	 << "the number of electrons.\n";
  }
  if (numOrbs != numOcc) 
    perr << "Warning:  Total number of localized orbitals does not equal"
	 << " the number of occupied bands.\n";
  if (numOrbs > numOcc) {
    perr << "Cannot create more independent localized orbitals than "
	 << "the number of bands we started with.\n";
    abort();
  }
  // Write out which orbitals are being used for localization
  if (Comm.MyProc() == 0) 
    for (int ki=0; ki<numk; ki++) {
      for (int primki=0; primki<TileMap[spin].GetNumPrimTwists(); primki++) {
	fprintf (stderr, "Occupied bands for primitive twist #%d:\n  ",
		 primki);
	for (int bi=0; bi<numOcc; bi++) {
	  int spinIndex, twistIndex, bandIndex;
	  int primTwistIndex, primBandIndex;
	  orbs[spin](ki,bi)->GetLabels (spinIndex, twistIndex, bandIndex);
	  //	  fprintf (stderr, "Super spin=%d ki=%d bi=%d    ", spinIndex, twistIndex, bandIndex);
	  TileMap[spin].Super2Prim (twistIndex, bandIndex, 
				    primTwistIndex, primBandIndex);
	  //fprintf (stderr, "Prim spin=%d ki=%d bi=%d\n", 
	  //         spinIndex, primTwistIndex, primBandIndex);
	  if (primTwistIndex == primki) 
	    fprintf (stderr, "%d ", primBandIndex);
	}
	fprintf (stderr, "\n");
      }
    }      

  // Allocate new orbital coefficients
  NewOrbCoefs[spin].resize(numk, numOcc, numOcc);
  NewOrbNorms[spin].resize(numk, numOcc);
  NewEigVals[spin].resize (numk, numOcc);
  // Initialize localized orbital matrices
  for (int ki=0; ki<numk; ki++)
    for (int bi=0; bi<numOcc; bi++) {
      for (int bj=0; bj<numOcc; bj++) 
	NewOrbCoefs[spin](ki,bi,bj) = 0.0;
      NewOrbCoefs[spin](ki,bi,bi) = 1.0;
      NewOrbNorms[spin](ki,bi) = 1.0;
      NewEigVals[spin](ki,bi)  = orbs[spin](ki,bi)->GetEigVal();
    }
  
  // Outer loop is over k-points
  for (int ki=0; ki<numk; ki++) {
    perr << "\nDoing MPI-based localization for k-point #" << ki << ":\n";
    // Allocate memory for storage of FFTs
    int nx, ny, nz;
    FFT.GetDims(nx,ny,nz);
    // First, store my processors FFTs
    int myFirstBand, myLastBand;
    ProcBands (Comm.MyProc(), numOcc, myFirstBand, myLastBand);
    Array<complex<double>,4> myBands(myLastBand-myFirstBand+1, nx, ny, nz);
    // I must also store the FFTs from other processors bands
    int first, last, maxBands;
    maxBands = 0;
    for (int proc=0; proc<Comm.NumProcs(); proc++) {
      ProcBands (proc, numOcc, first, last);
      int nBands = last-first+1;
      maxBands = max (nBands, maxBands);
    }
    Array<complex<double>,4> otherBands (maxBands, nx, ny, nz);
    
    /// FFT my bands and store them in myBands
    for (int band=myFirstBand; band<=myLastBand; band++) {
      double nrm = 0.0;
      orbs[spin](ki, band)->PutInFFTBox();
      fprintf (stderr, "rank = %d  FFT.size() = (%d,%d,%d)\n",
	       Comm.MyProc(), nx, ny, nz);
      FFT.k2r();
      for (int ix=0; ix<nx; ix++)
	for (int iy=0; iy<ny; iy++)
	  for (int iz=0; iz<nz; iz++) {
	    myBands(band-myFirstBand, ix, iy, iz) = 
	      promote(FFT.rBox(ix,iy,iz));
	    nrm += norm (FFT.rBox(ix,iy,iz));
	  }
    }
    Array<complex<double>,3> A(numOcc, numOcc, Centers.size()),
      Asum(numOcc, numOcc, Centers.size());
    A = complex<double>();
    // Compute my entries in A

    for (int b1=myFirstBand; b1<=myLastBand; b1++) {
      for (int b2=b1; b2<=myLastBand; b2++) {
	complex<double>* my1 = &(myBands(b1-myFirstBand,0,0,0));
	complex<double>* my2 = &(myBands(b2-myFirstBand,0,0,0));
	for (int ci=0; ci<Centers.size(); ci++) {
	  CenterClass &center = Centers[ci];
	  for (int j=0; j<center.Offsets.size(); j++) {
	    int offset = center.Offsets[j];
	    A(b1,b2,ci) += conj(my1[center.Offsets[j]])*my2[center.Offsets[j]];
	  }
	}
// 	for (int ix=0; ix<nx; ix++)
// 	  for (int iy=0; iy<ny; iy++)
// 	    for (int iz=0; iz<nz; iz++) {
// 	      complex<double> z = 
// 		conj(myBands(b1-myFirstBand, ix, iy, iz)) *
// 		myBands(b2-myFirstBand, ix, iy, iz);
// 	      for (int ci=0; ci<Centers.size(); ci++) 
// 		A(b1,b2,ci) += z*Centers[ci].Bitfield.GetDouble (ix,iy,iz);
// 	    }
	for (int ci=0; ci<Centers.size(); ci++)
	  A(b2,b1,ci) = conj(A(b1,b2,ci));
      }
    }
    // Now compute elements involving two different processors' bands.
    // Loop through other processors. 
    for (int offset=1; offset<Comm.NumProcs(); offset++) {
      /////////////////////////////
      // First, do communication //
      /////////////////////////////
      clock_t start = clock();
      int sendProc = (Comm.MyProc()+offset) % Comm.NumProcs();
      int recvProc = (Comm.MyProc()-offset + Comm.NumProcs()) % Comm.NumProcs();
      int otherFirst, otherLast;
      ProcBands(recvProc, numOcc, otherFirst, otherLast);
      int numSend = myLastBand-myFirstBand+1;
      int numRecv = otherLast-otherFirst+1;
      int numSendRecv = min (numSend, numRecv);
      for (int i=0; i<numSendRecv; i++)
	Comm.SendReceive 
	  (sendProc, myBands(i, Range::all(), Range::all(), Range::all()),
	   recvProc, otherBands (i, Range::all(), Range::all(), Range::all()));
      if (numSend > numRecv)
	for (int i=numSendRecv; i<numSend; i++)
	  Comm.Send (sendProc, myBands(i,Range::all(),Range::all(),Range::all()));
      else if (numRecv > numSend)
	for (int i=numSendRecv; i<numRecv; i++)
	  Comm.Receive (recvProc, otherBands(i,Range::all(),Range::all(),Range::all()));
      clock_t end = clock();
      perr << "Time for communication = " 
	   << (double)(end-start)/(double)(CLOCKS_PER_SEC) << "s.\n";
      start = end;
      //////////////////////////////////////
      // Now compute matrix elements of A //
      //////////////////////////////////////
      for (int b1=myFirstBand; b1<=myLastBand; b1++) 
	for (int b2=otherFirst; b2<=otherLast; b2++) {
// 	  for (int ix=0; ix<nx; ix++)
// 	    for (int iy=0; iy<ny; iy++)
// 	      for (int iz=0; iz<nz; iz++) {
// 		complex<double> z = 
// 		  conj(myBands(b1-myFirstBand, ix, iy, iz)) *
// 		  otherBands(b2-otherFirst, ix, iy, iz);
// 		if (b1 < b2)
// 		  for (int ci=0; ci<Centers.size()/2; ci++)
// 		    A(b1,b2,ci) += z*Centers[ci].Bitfield.GetDouble (ix,iy,iz);
// 		else
// 		  for (int ci=Centers.size()/2; ci<Centers.size(); ci++)
// 		    A(b1,b2,ci) += z*Centers[ci].Bitfield.GetDouble (ix,iy,iz);
// 	      }
// 	  for (int ci=0; ci<Centers.size(); ci++)
// 	    A(b2,b1,ci) = conj(A(b1,b2,ci));
	  int c1 = (b1 < b2) ?        0         : Centers.size()/2;
	  int c2 = (b1 < b2) ? Centers.size()/2 : Centers.size();
	  complex<double>* my0    = &(myBands(b1-myFirstBand,0,0,0));
	  complex<double>* other0 = &(otherBands(b2-otherFirst,0,0,0));
	  for (int ci=c1; ci<c2; ci++) {
	    CenterClass &center = Centers[ci];
	    for (int j=0; j<center.Offsets.size(); j++) {
	      int offset = center.Offsets[j];
	      A(b1,b2,ci) += conj (my0[offset])*other0[offset];
	    }
	  }
	  for (int ci=0; ci<Centers.size(); ci++)
 	    A(b2,b1,ci) = conj(A(b1,b2,ci));
	}
      end = clock();
      perr << "Time for A matrix compute = " 
	   << (double)(end-start)/(double)(CLOCKS_PER_SEC) << "s.\n";
    }
    double nrm = 1.0/(double)(nx*ny*nz);
    A *= (-nrm);
    Comm.AllSum (A, Asum);
    A = Asum;
    // for (int i=0; i<A.extent(0); i++) {
    //   for (int j=0; j<A.extent(1); j++)
    // 	fprintf (stderr, "%10.4e %10.4e  ", 
    // 		 A(i,j,0).real(),
    // 		 A(i,j,0).imag());
    //   fprintf (stderr, "\n");
    // }

    ////////////////////////////////////
    // Now diagonalize the A matrices //
    ////////////////////////////////////
    Array<complex<double>,2> eigvecs, Aone(numOcc, numOcc);
    Array<double,1> eigvals;
    
    if (Comm.MyProc()==0) {
      fprintf (stderr, "Spin = %s\n", (spin==0) ? "up" : "down");
      fprintf (stderr, "          Center             Norm     Nearest Ion\n");
      int orb = 0;
      for (int ci=0; ci<Centers.size(); ci++) {
	CenterClass &center = Centers[ci];
	fprintf (stderr, "NumUp = %d   NumDown = %d\n", 
		 center.NumUp, center.NumDown);
	Aone = A(Range::all(), Range::all(), ci);
	int nOrbs = (spin == 0) ? center.NumUp : center.NumDown;
	SymmEigenPairs(Aone, nOrbs, eigvals, eigvecs);
	for (int i=0; i<nOrbs; i++) {
	  NewOrbNorms[spin](ki,orb) = -eigvals(i);
	  NewOrbCoefs[spin](ki,orb,Range::all()) = eigvecs(i,Range::all());
	  fprintf (stderr, "(%7.4f %7.4f %7.4f) %10.7f       %d\n", 
		   center.r[0], center.r[1], center.r[2],
		   -eigvals(i), NearestIon(center.r));
	  // Store center in each orbital
	  orbs[spin](ki, orb)->SetCenter(center);
	  orb++;
	}
      }

      for (int i=0; i<numOcc; i++)
	fprintf (stderr, "NewOrbCoefs[0](0,%d) = %12.8e + %12.8ei\n", 
		 i, NewOrbCoefs[0](0,i).real(), NewOrbCoefs[0](0,i).imag());
      // Fill rest of coefs with identity.  This is for the nonlocalized
      // bands.
      for (int i=orb; i<numOcc; i++)
	for (int j=0; j<numOcc; j++)
	  NewOrbCoefs[spin](ki,i,j) = 0.0;
      for (int i=orb; i<numOcc; i++)
	NewOrbCoefs[spin](ki,i,i) = 1.0;
      

      // Orthogonalize
      if (ortho) {
	// PolarOrthogonalize orthogonalizes the columns, so we transpose
	// first.  
	Array<complex<double>,2> orthoMat;
	orthoMat.reference (NewOrbCoefs[spin](ki,Range::all(),Range::all()));
	Transpose(orthoMat);
	PolarOrthogonalize (orthoMat);
	Transpose(orthoMat);
      } 
    }
  }
  Comm.BarrierSync();
}


void
OrbitalSetClass::ProcBands(int proc, int numOcc, int &firstBand, int &lastBand)
{
  int numProcs = Comm.NumProcs();
  int first = 0;
  int last;
  for (int pr=0; pr<numProcs; pr++) {
    int bands = numOcc / numProcs + (pr < (numOcc%numProcs) ? 1 : 0);
    last = first + bands-1;
    if (pr == proc) {
      firstBand = first;
      lastBand  = last;
    }
    first = last+1;
  }
  
}


void
OrbitalSetClass::ProcTasks (int proc, int numTasks, int &firstTask, int &lastTask)
{
  int numProcs = Comm.NumProcs();
  int first = 0;
  int last;
  for (int pr=0; pr<numProcs; pr++) {
    int bands = numTasks / numProcs + (pr < (numTasks%numProcs) ? 1 : 0);
    last = first + bands-1;
    if (pr == proc) {
      firstTask = first;
      lastTask  = last;
    }
    first = last+1;
  }

}
