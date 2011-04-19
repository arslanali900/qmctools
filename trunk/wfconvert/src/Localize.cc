#include "PlaneWaveClass.h"
#include "ParserClass.h"
#include "ParseCommand.h"
#include <Common/Splines/ComplexMultiTricubicSpline.h>
#include <Common/MatrixOps/MatrixOps.h>
//#include <sys/sysinfo.h>
#include <sys/types.h>
#include <unistd.h>
#include <sstream>

// unsigned long 
// FreeMemory2()
// {
//   struct sysinfo sysInfo;

//   sysinfo (&sysInfo);
//   return sysInfo.mem_unit*
//     (sysInfo.freeram + sysInfo.bufferram + sysInfo.sharedram);
// }


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
//     cerr << "memFree = " << memFree 
// 	 << "   FreeMemory2 = " << FreeMemory2() << endl;
    memParser.CloseFile();
    return memFree;
  }
  else // Default to 256MB
    return 268435456;
}




////////////////////////////////////////////////////////////
//                Localization functions                  //
////////////////////////////////////////////////////////////

bool 
PlaneWaveSystem::ReadCenters(string fname) 
{
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
    SetBitfield (Centers[i]);
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



// The following code roughly implements the nonorthonal localization
// methods of Reboredo and Williamson, PRB 71, 12115(R) (2005)
void
PlaneWaveSystem::Localize (string centersFilename, double radius)
{
  CellClass &cell = WriteSuper ? SuperCell : PrimCell;
  FFTBox &FFT = cell.FFT;
  if (Comm.MyProc() == 0) {
    cell.SetupFFT();
    if (centersFilename == "random") {
      int numElecs = NumElectrons;
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
	SetBitfield (Centers[i]);
      }
    }
    else
      ReadCenters (centersFilename);
  }
  
  if (Comm.NumProcs() > 1) {
    DistributeBands();
    for (int spin=0; spin < 2; spin++) {
      if (OptimizeCenters)
	do {
	  LocalizeMPI(true, spin);
	} while (!UpdateCenters(spin));
      
      LocalizeMPI(Orthogonalize, spin);
    }
  }
  else {
    // Now do the localization  
    for (int spin=0; spin<2; spin++) {
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
PlaneWaveSystem::UpdateCenters(int spin)
{
  CellClass &cell = WriteSuper ? SuperCell : PrimCell;
  TinyVector<Array<OrbitalClass*,2>,2> &orbs =
    WriteSuper ? SuperOrbitals : Orbitals;
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
      cell.PutOrbital(*orbs[spin](k, band));
      FFT.k2r();
      double sx, sy, sz;
      for (int ix=0; ix<nx; ix++) {
	sx = (double)ix * nxInv;
	for (int iy=0; iy<ny; iy++) {
	  sy = (double)iy * nyInv;
	  for (int iz=0; iz<nz; iz++) {
	    sz = (double) iz * nzInv;
	    Vec3 r = sx*r0 + sy*r1 + sz*r2;
	    Vec3 diff = MinImage(r - center.r);
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
    SetBitfield(Centers[ci]);
  }

  return converged;
}


void 
PlaneWaveSystem::SetBitfield (CenterClass &center)
{
  CellClass &cell = WriteSuper ? SuperCell : PrimCell;
  FFTBox &FFT = cell.FFT;
  cell.SetupFFT();

  int nx, ny, nz;
  FFT.GetDims(nx,ny,nz);
//   fprintf (stderr, "FFT box is %dx%d%d in SetBitField.\n",
// 	   nx, ny, nz);
  double sx, sy, sz;
  double nxInv = 1.0/(double)nx;
  double nyInv = 1.0/(double)ny;
  double nzInv = 1.0/(double)nz;
  Vec3 r;
  Vec3 a0 = cell.Lattice.a(0);
  Vec3 a1 = cell.Lattice.a(1);
  Vec3 a2 = cell.Lattice.a(2);
  FFT.GetDims (nx, ny, nz);
  center.Bitfield.Init (nx, ny, nz);
  for (int ix=0; ix<nx; ix++) {
    sx = (double)ix * nxInv;
    for (int iy=0; iy<ny; iy++) {
      sy = (double)iy * nyInv;
      for (int iz=0; iz<nz; iz++) {
	sz = (double) iz * nzInv;
	r = sx*a0 + sy*a1 + sz*a2;
	Vec3 diff = MinImage(r - center.r);
	center.Bitfield.Set
	  (ix,iy,iz, (dot(diff,diff) < center.Radius*center.Radius));
      }
    }
  }
}

void PlaneWaveSystem::LocalFunc(CenterClass center,
				Array<double,3> &func)
{
  CellClass &cell = WriteSuper ? SuperCell : PrimCell;
  FFTBox &FFT = cell.FFT;
  
  int nx, ny, nz;
  FFT.GetDims(nx,ny,nz);
  double sx, sy, sz;
  double nxInv = 1.0/(double)nx;
  double nyInv = 1.0/(double)ny;
  double nzInv = 1.0/(double)nz;
  Vec3 r;
  Vec3 a0 = cell.Lattice.a(0);
  Vec3 a1 = cell.Lattice.a(1);
  Vec3 a2 = cell.Lattice.a(2);

  FFT.GetDims (nx, ny, nz);
  for (int ix=0; ix<nx; ix++) {
    sx = (double)ix * nxInv;
    for (int iy=0; iy<ny; iy++) {
      sy = (double)iy * nyInv;
      for (int iz=0; iz<nz; iz++) {
	sz = (double) iz * nzInv;
	r = (sx+0.5)*a0 + (sy+0.5)*a1 + (sz+0.5)*a2;
	Vec3 diff = MinImage(r - center.r);
	func(ix,iy,iz) =  
	  (dot(diff,diff) < center.Radius*center.Radius) ? 1.0 : 0.0;
      }
    }
  }
  func *= (nxInv*nyInv*nzInv);
}

inline complex<double> promote (complex<float> val)
{ return complex<double>(val.real(), val.imag()); }

inline complex<double> promote (complex<double> val)
{ return val; }


int
PlaneWaveSystem::NearestIon(Vec3 pos)
{
  CellClass &cell = WriteSuper ? SuperCell : PrimCell;

  if (cell.IonPos.size() < 1) 
    return -1;
  Vec3 disp;
  double minDist;
  int minIon;

  disp    = MinImage (pos-cell.IonPos(0));
  minDist = sqrt (dot(disp, disp));
  minIon  = 0;
  for (int i=0; i<cell.IonPos.size(); i++) {
    disp = MinImage (pos-cell.IonPos(i));
    double dist = sqrt (dot (disp, disp));
    if (dist < minDist) {
      minIon = i;
      minDist = dist;
    }
  }
  return minIon;
}
    



void
PlaneWaveSystem::Localize(bool ortho, int spin)
{
  CellClass &cell = WriteSuper ? SuperCell : PrimCell;
  TinyVector<Array<OrbitalClass*,2>,2> &orbs =
    WriteSuper ? SuperOrbitals : Orbitals;
  FFTBox &FFT = cell.FFT;

  int N = orbs[spin](0,0)->GetCoefs().size();
  int numk = orbs[spin].extent(0);
  // HACK HACK HACK
  //int numOcc = orbs[spin].extent(1);
  int numOcc = (spin==0) ? UseUpBands : UseDownBands;
  cerr << "Using " << numOcc << " bands for "
       << ((spin==0) ? "up" : "down") << "-spin electrons.\n";
  int numOrbs = 0;
  vector<CenterClass>::iterator iter;
  int elecs = 0;
  for (iter = Centers.begin(); iter != Centers.end(); iter++) {
    elecs += iter->NumUp + iter->NumDown;
    numOrbs += iter->NumOrbitals;
  }
  if (elecs != NumElectrons) {
    cerr << "Error:  Number of localized orbitals should equal "
	 << "the number of occupied bands.\n";
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

  for (int ki=0; ki < numk; ki++) {
    while (row < numOcc) {
      cerr << "row = " << row << " of " << numOcc << endl;
      // How many rows to do in this iteration
      int todo = min (numFFTs, numOcc-row);
      for (int ri=0; ri<todo; ri++) {
	int band1 = (row+ri) / numk;
	int k1    = (row+ri) % numk;
	cell.PutOrbital (*orbs[spin](k1,band1));
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
	    for (int ix=0; ix<nx; ix++)
	      for (int iy=0; iy<ny; iy++)
		for (int iz=0; iz<nz; iz++) 
		  theta(ci) += nrm*center.Bitfield.GetDouble (ix,iy,iz)*
		    (phi_rows(rowi, ix,iy,iz) * conj(phi_rows(coli, ix,iy,iz)));
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
	cell.PutOrbital(*orbs[spin](k, band));
	FFT.k2r();
	for (int rowi=0; rowi < todo; rowi++) {
	  theta = complex<double>();
	  int ci=0;
	  for (iter = Centers.begin(); iter != Centers.end(); iter++) {
	    CenterClass &center = (*iter);
	    for (int ix=0; ix<nx; ix++)
	      for (int iy=0; iy<ny; iy++)
		for (int iz=0; iz<nz; iz++) 
		  theta(ci) += nrm*center.Bitfield.GetDouble (ix,iy,iz)*
		    (phi_rows(rowi,ix,iy,iz) * promote(FFT.rBox(ix,iy,iz)));
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
      for (int i=0; i<center.NumOrbitals; i++) {
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
  
  // Now, rotate the orbitals with the new coefficients
//   Array<complex<double>,2> dummy(1,N);
//   tmp.WriteVar("NewCoefs", dummy);
//   IO::IOVarBase &coefsVar = *(tmp.GetVarPtr("NewCoefs"));
//   coefsVar.Resize(numOcc);
//   Array<complex<double>,1> coefs(N);
//   for (int orb=0; orb < numOcc; orb++) {
//     coefs = complex<double>();
//     for (int i=0; i<numOcc; i++) {
//       int band = i / numk;
//       int k    = i % numk;
//       coefs += newOrbCoefs(orb,i)*orbs[spin](k, band)->GetCoefs();
//     }
//     coefsVar.Write (coefs, orb, Range::all());
//   }

//   // Finally, overwrite the present orbital coefficients with the
//   // localized ones
//   for (int orb=0; orb < numOcc; orb++) {
//     int band = orb / numk;
//     int k    = orb % numk;
//     coefsVar.Read (coefs, orb, Range::all());
//     orbs[spin](k, band)->GetCoefs() = coefs;
//     orbs[spin](k, band)->SetNorm2(newOrbNorms(orb));
//   }
  tmp.CloseFile();
}

void
CenterClass::Broadcast(CommunicatorClass &comm, int root)
{
  comm.Broadcast (root, r);
  comm.Broadcast (root, Radius);
  comm.Broadcast (root, Spherical);
  comm.Broadcast (root, NumOrbitals);
  comm.Broadcast (root, NumUp);
  comm.Broadcast (root, NumDown);
}

void
CellClass::Broadcast (CommunicatorClass &comm, int root)
{
  comm.Broadcast (root, IsSuper);
  Mat3 lattice = Lattice.GetDirect();
  comm.Broadcast (root, lattice);
  Lattice.SetDirect(lattice);
  int numIons = IonPos.size();
  comm.Broadcast (root, numIons);
  IonPos.resize(numIons);
  AtomTypes.resize(numIons);
  comm.Broadcast (root, IonPos);
  comm.Broadcast (root, AtomTypes);

  // Now send the GVecs
  GVecs.Broadcast (comm, root);
}

void
kPointClass::Broadcast (CommunicatorClass &comm, int root)
{
  // The PrimCell and SuperCell pointers should be set independently.
  // Set the indices
  int  primSize =  PrimIndices.size();
  int superSize = SuperIndices.size();
  comm.Broadcast (root, primSize);
  comm.Broadcast (root, superSize);
  if (root != comm.MyProc()) {
    PrimIndices.resize(primSize);
    SuperIndices.resize(superSize);
  }
  comm.Broadcast (root, PrimIndices);
  comm.Broadcast (root, SuperIndices);
  comm.Broadcast (root, PrimTwist);
  comm.Broadcast (root, SuperTwistInt);
  comm.Broadcast (root, SuperTwistFrac);
}


void
PlaneWaveSystem::DistributeBands()
{
  CellClass &cell = WriteSuper ? SuperCell : PrimCell;
  TinyVector<Array<OrbitalClass*,2>,2> &orbs =
    WriteSuper ? SuperOrbitals : Orbitals;
  FFTBox &FFT = cell.FFT;

  // First, distribute all the essential data
  // Distribute the cells
  PrimCell.Broadcast  (Comm, 0);
  SuperCell.Broadcast (Comm, 0);
  
  Comm.Broadcast (0, NumElectrons);
  Comm.Broadcast (0, UseUpBands);
  Comm.Broadcast (0, UseDownBands);

  // Distribute the kPoints
  int numPrimk = PrimkPoints.size();
  Comm.Broadcast (0, numPrimk);
  if (PrimkPoints.size() != numPrimk) PrimkPoints.resize(numPrimk);
  for (int i=0; i<numPrimk; i++) {
    // The broadcast will fix up the k.  We just need to set the cell pointers.
    Vec3 k(0.0, 0.0, 0.0);
    if (Comm.MyProc() != 0)
      PrimkPoints(i).SetCells (PrimCell, SuperCell, k);
    PrimkPoints(i).Broadcast (Comm, 0);
  }

  // Send orbitals
  for (int spin=0; spin<2; spin++) {
    int numk     = orbs[spin].extent(0);
    int numBands = orbs[spin].extent(1);
    Comm.Broadcast (0, numk);
    Comm.Broadcast (0, numBands);

    if (numBands > 0) {
      Array<int,2> primTwistIndices(numk, numBands);
      if (Comm.MyProc() == 0) 
	for (int ik=0; ik<numk; ik++)
	  for (int iband=0; iband<numBands; iband++)
	    primTwistIndices(ik,iband) = orbs[spin](ik,iband)->GetPrimTwistIndex();
      Comm.Broadcast (0, primTwistIndices);
      if (Comm.MyProc() == 1)
      // Create orbitals if I'm not proc 0
      if ((orbs[spin].extent(0) != numk) ||
	  (orbs[spin].extent(1) != numBands)) {
	orbs[spin].resize(numk, numBands);
	for (int ik=0; ik<numk; ik++) 
	  for (int iband=0; iband<numBands; iband++) {
	    int twistIndex = primTwistIndices(ik,iband);
	    orbs[spin](ik,iband) = 
	      new OrbitalClass(PrimCell, SuperCell, PrimkPoints(twistIndex));
	    orbs[spin](ik,iband)->GetCoefs().resize(cell.GVecs.size());
	  }
      }
      // Now, broadcast eigenvalues and eigenvectors
      for (int ik=0; ik<numk; ik++)
	for (int iband=0; iband<numBands; iband++) {
	  OrbitalClass &orb = *orbs[spin](ik,iband);
	  Comm.Broadcast (0, orb.GetCoefs());
	  double eigval = orb.GetEigVal();
	  Comm.Broadcast (0, eigval);
	  orb.SetEigVal (eigval);
	  orb.SetSpin (spin);
	}
    }
  }
  // Send the centers
  int numCenters = Centers.size();
  Comm.Broadcast (0, numCenters);
  if (Comm.MyProc() != 0) 
    Centers.resize(numCenters);
  for (int i=0; i<numCenters; i++) {
    Centers[i].Broadcast(Comm, 0);
    SetBitfield (Centers[i]);
  }
  
//   else {
//     int numIons, numGVecs;
//     Comm.Broadcast (0, numIons);
//     cell.IonPos.resize(numIons);
//     cell.AtomTypes.resize(numIons);
//     Comm.Broadcast (0, cell.IonPos);
//     // Receive the G-vectors
//     Comm.Broadcast (0, numGVecs);
//     GVecsArray.resize(numGVecs);
//     Comm.Broadcast (0, GVecsArray);
//     Comm.Broadcast (0, FFTFactor);
//     cell.GVecs.Set (cell.Lattice.GetDirect(), GVecsArray, FFTFactor);
//     cell.SetupFFT();
//     // Receive the k-Point objects

//     // Receive the orbitals
//     for (int spin = 0; spin < 2; spin++) {
//       int numk, numBands;
//       Comm.Broadcast (0, numk);
//       Comm.Broadcast (0, numBands);
//       if (numBands > 0) {
// 	orbs[spin].resize(numk, numBands);
// 	for (int ik=0; ik<numk; ik++)
// 	  for (int iband=0; iband<numBands; iband++) {
// 	    orbs[spin](ik,iband) = 
// 	      new OrbitalClass(PrimCell, SuperCell, PrimkPoints(ik));
// 	    OrbitalClass &orb = *orbs[spin](ik,iband);
// 	    orb.GetCoefs().resize(cell.GVecs.size());
// 	    Comm.Broadcast (0, orb.GetCoefs());
// 	    double eigval;
// 	    Comm.Broadcast (0, eigval);
// 	    orb.SetEigVal (eigval);
// 	    Vec3 k;
// 	    Comm.Broadcast (0, k);
// 	    // HACK HACK HACK
// 	    // orb.Setk (k);
// 	  }
//       }
//     }
//     // Receive the centers
//     int numCenters;
//     Comm.Broadcast (0, numCenters);
//     Centers.resize(numCenters);
//     for (int i=0; i<numCenters; i++) {
//       Centers[i].Broadcast(Comm);
//       SetBitfield (Centers[i]);
//     }
//   }
}



void
PlaneWaveSystem::LocalizeMPI(bool ortho, int spin)
{
  CellClass &cell = WriteSuper ? SuperCell : PrimCell;
  TinyVector<Array<OrbitalClass*,2>,2> &orbs =
    WriteSuper ? SuperOrbitals : Orbitals;
  FFTBox &FFT = cell.FFT;

  cell.SetupFFT();

  int N        = orbs[spin](0,0)->GetCoefs().size();
  int numk     = orbs[spin].extent(0);
  int numBands = orbs[spin].extent(1);

  int numOcc = (spin==0) ? UseUpBands : UseDownBands;
  perr << "Using " << numOcc << " bands in localizing spin "
       << ((spin==0) ? "up" : "down") << " orbitals.\n";

  int numOrbs = 0;
  int elecs = 0;
  vector<CenterClass>::iterator iter;
  for (iter = Centers.begin(); iter != Centers.end(); iter++) {
    elecs += iter->NumUp + iter->NumDown;
    numOrbs += (spin == 0) ? iter->NumUp : iter->NumDown;
  }
  if (elecs != NumElectrons) {
    perr << "Error:  Number of localized orbitals should equal "
	 << "the number of occupied bands.\n";
  }
  if (numOrbs > numOcc) {
    perr << "Cannot create more independent localized orbitals than "
	 << "the number of bands we started with.\n";
    abort();
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
  
  // Out loop is over k-points
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
      cell.PutOrbital (*orbs[spin](ki, band));
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
	for (int ix=0; ix<nx; ix++)
	  for (int iy=0; iy<ny; iy++)
	    for (int iz=0; iz<nz; iz++) {
	      complex<double> z = 
		conj(myBands(b1-myFirstBand, ix, iy, iz)) *
		myBands(b2-myFirstBand, ix, iy, iz);
	      for (int ci=0; ci<Centers.size(); ci++)
		A(b1,b2,ci) += z*Centers[ci].Bitfield.GetDouble (ix,iy,iz);
	    }
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
      //////////////////////////////////////
      // Now compute matrix elements of A //
      //////////////////////////////////////
      for (int b1=myFirstBand; b1<=myLastBand; b1++) 
	for (int b2=otherFirst; b2<=otherLast; b2++) {
	  for (int ix=0; ix<nx; ix++)
	    for (int iy=0; iy<ny; iy++)
	      for (int iz=0; iz<nz; iz++) {
		complex<double> z = 
		  conj(myBands(b1-myFirstBand, ix, iy, iz)) *
		  otherBands(b2-otherFirst, ix, iy, iz);
		if (b1 < b2)
		  for (int ci=0; ci<Centers.size()/2; ci++)
		    A(b1,b2,ci) += z*Centers[ci].Bitfield.GetDouble (ix,iy,iz);
		else
		  for (int ci=Centers.size()/2; ci<Centers.size(); ci++)
		    A(b1,b2,ci) += z*Centers[ci].Bitfield.GetDouble (ix,iy,iz);
	      }
	  for (int ci=0; ci<Centers.size(); ci++)
	    A(b2,b1,ci) = conj(A(b1,b2,ci));
	}
    }
    double nrm = 1.0/(double)(nx*ny*nz);
    A *= (-nrm);
    Comm.AllSum (A, Asum);
    A = Asum;
    
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
	fprintf (stderr, "NumUp = %d   NumDown = %d\n", center.NumUp, center.NumDown);
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

    ////////////////////////////////////////////////////////////////
    // UPDATE:  We will now put off the replacing of the          //
    // coefficients until we write the file or have to update the //
    // centers.                                                   // 
    ////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////
    // Now rotate the occupied orbitals to localize //
    //////////////////////////////////////////////////
//     Array<complex<double>,2> dummy(1,N);
//     IO::IOSectionClass tmp;
//     pid_t pid = getpid();
//     stringstream fname;
//     fname << "NewCoefs_" << pid << ".h5";
//     tmp.NewFile (fname.str());
//     tmp.WriteVar("NewCoefs", dummy);
//     IO::IOVarBase &coefsVar = *(tmp.GetVarPtr("NewCoefs"));
//     coefsVar.Resize(numOcc);
//     Array<complex<double>,1> coefs(N);
//     for (int orb=0; orb < numOcc; orb++) {
//       coefs = complex<double>();
//       for (int i=0; i<numOcc; i++) {
// 	int band = i / numk;
// 	int k    = i % numk;
// 	coefs += newOrbCoefs(orb,i)*orbs[spin](k, band)->GetCoefs();
//       }
//       coefsVar.Write (coefs, orb, Range::all());
//     }
    
    // Finally, overwrite the present orbital coefficients with the
    // localized ones
//     for (int orb=0; orb < numOcc; orb++) {
//       int band = orb / numk;
//       int k    = orb % numk;
//       coefsVar.Read (coefs, orb, Range::all());
//       orbs[spin](k, band)->GetCoefs() = coefs;
//       orbs[spin](k, band)->SetNorm2(NewOrbNorms(orb));
//     }
//     tmp.CloseFile();
//     unlink (fname.str().c_str());
//   }
  
//   // Now broadcast the new coefficients to all processors
//   for (int ki=0; ki<numk; ki++)
//     for (int band=0; band<numBands; band++)
//       Comm.Broadcast (0, orbs[spin](ki, band)->GetCoefs());

//   // Associate the centers with each orbital
//   int orb=0;
//   for (int ci=0; ci<Centers.size(); ci++) {
//     int nOrbs = (spin == 0) ? Centers[ci].NumUp: Centers[ci].NumDown;
//     for (int iorb=0; iorb<nOrbs; iorb++) {
//       int band = orb / numk;
//       int k    = orb % numk;
//       orbs[spin](k, band)->SetCenter(Centers[ci]);
//       orb++;
//     }
//   }

//   if (Comm.MyProc() == 0) {
//     FILE *fout = fopen ("A.dat", "w");
//     for (int i=0; i<A.extent(0); i++) {
//       for (int j=0; j<A.extent(1); j++)
// 	fprintf (fout, "%23.16e ", A(i,j,0).real());
//       fprintf (fout, "\n");
//     }
//     fclose (fout);
//   }
  }
  Comm.BarrierSync();
}


void
PlaneWaveSystem::ProcBands(int proc, int numOcc, int &firstBand, int &lastBand)
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
PlaneWaveSystem::ProcTasks (int proc, int numTasks, int &firstTask, int &lastTask)
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
