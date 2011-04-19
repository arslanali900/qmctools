#include "OrbitalSetClass.h"
#include "ReadABINIT.h"

void
ABINIT_header::SkipRecHeader(FILE *fin)
{
  char dummy[100];
  fread (dummy, sizeof(char), recHeaderLength, fin);
}

void
ABINIT_header::Read (FILE* fin, bool print)
{
  recHeaderLength = 0;
  char ch='\0';
  while (ch != '.') {
    fread(&ch, sizeof(char), 1, fin);
    recHeaderLength++;
  }
  fseek (fin, (long)0, SEEK_SET);
  recHeaderLength -= 2;
  char dummy[100];
  fread (dummy, sizeof(char), recHeaderLength, fin);

  char vsn[7];
  fread (vsn, sizeof(char), 6, fin);
  vsn[6] = '\0';
  codvsn = vsn;
  if (print)
    cerr << "codsvn = \"" << codvsn << "\"" << endl;

  fread (&headform, sizeof(int), 1, fin);
  if (print)
    cerr << "headform = " << headform << endl;
  fread (&fform, sizeof(int), 1, fin);
  
  // Skip FORTRAN record end and start:
  fread (dummy, sizeof(char), 2*recHeaderLength, fin);
  
  fread (&bantot,    sizeof(int),    1, fin);
  fread (&date,      sizeof(int),    1, fin);
  fread (&intxc,     sizeof(int),    1, fin);
  fread (&ixc,       sizeof(int),    1, fin);
  fread (&natom,     sizeof(int),    1, fin);
  fread (&ngfft,     sizeof(int),    3, fin);
  fread (&nkpt,      sizeof(int),    1, fin);
  fread (&nspden,    sizeof(int),    1, fin);
  fread (&nspinor,   sizeof(int),    1, fin);
  fread (&nsppol,    sizeof(int),    1, fin);
  fread (&nsym,      sizeof(int),    1, fin);
  fread (&npsp,      sizeof(int),    1, fin);
  fread (&ntypat,    sizeof(int),    1, fin);
  fread (&occopt,    sizeof(int),    1, fin);
  fread (&pertcase,  sizeof(int),    1, fin);
  fread (&usepaw,    sizeof(int),    1, fin);
  fread (&ecut,      sizeof(double), 1, fin);
  fread (&ecutdg,    sizeof(double), 1, fin);
  fread (&ecutsm,    sizeof(double), 1, fin);
  fread (&ecut_eff,  sizeof(double), 1, fin);
  fread (&(qptn[0]), sizeof(double), 3, fin);
  fread (&(rprimd[0][0]), sizeof(double), 9, fin);
  // for (int i=0; i<3; i++)
  //   for (int j=i+1; j<3; j++)
  //     swap (rprimd[i][j], rprimd[j][i]);
  if (print)
    cerr << "rprimd = " << endl << "  " << rprimd[0] << "\n  "
	 << rprimd[1] << "\n  " << rprimd[2] << endl;
  fread (&stmbias,   sizeof(double), 1, fin);
  fread (&tphysel,   sizeof(double), 1, fin);
  fread (&tsmear,    sizeof(double), 1, fin);
  if (headform >= 57)
    fread(&usewvl, sizeof(int), 1, fin);
  //cerr << "tsmear = " << tsmear << endl;
  if (print) {
    cerr << "ecut = " << ecut << endl;
    cerr << "ngfft = " << ngfft[0] << " x " << ngfft[1] 
	 << " x " << ngfft[2] << endl;
  }

  // Skip FORTRAN record end and start:
  fread (dummy, sizeof(char), 2*recHeaderLength, fin);

  istwfk.resize(nkpt);
  fread (&(istwfk[0]), sizeof(int), nkpt, fin);
  for (int i=0; i<nkpt; i++)
    if (istwfk[i] != 1) {
      cerr << "ABINIT has used time-reversale symmetry in this calculation for k-point "
	   << i << ".\nPlease rerun without time-reversal symmetry.\n"
	   << "Use \"istwfk n*1\", where n is the number of k-points in your calculation.\n";
      abort();
    }
  if (print)
    cerr << endl;

  nband.resize(nkpt*nsppol);
  fread (&(nband[0]), sizeof(int), nkpt*nsppol, fin);

  npwarr.resize(nkpt);
  fread (&(npwarr[0]), sizeof(int), nkpt, fin);
  
  so_psp.resize(npsp);
  fread (&(so_psp[0]), sizeof(int), npsp, fin);
  
  symafm.resize(nsym);
  fread (&(symafm[0]), sizeof(int), nsym, fin);
  
  symrel.resize(nsym);
  fread (&(symrel[0](0,0)), sizeof (int), 9*nsym, fin);

  typat.resize(natom);
  fread (&(typat[0]), sizeof(int), natom, fin);
  
  kpt.resize(nkpt);
  fread (&(kpt[0]), sizeof(double), 3*nkpt, fin);
  if (print) {
    for (int i=0; i<nkpt; i++)
      fprintf (stderr, "kpt[%2d] = [%8.4f %8.4f %8.4f]\n",
	       i, kpt[i][0], kpt[i][1], kpt[i][2]);
    
    cerr << "bantot = " << bantot << endl;
  }

  occ.resize(bantot);
  fread (&(occ[0]), sizeof(double), bantot, fin);
  
  tnons.resize(nsym);
  fread (&(tnons[0]), sizeof(double), 3*nsym, fin);
  
  znucltypat.resize(ntypat);
  fread (&(znucltypat[0]), sizeof(double), ntypat, fin);

  if (print) {
    cerr << "Atom   Type #    Znucl\n";
    for (int n=0; n<natom; n++) {
      int i = typat[n] - 1;
      fprintf (stderr, "  %d      %2d       %3.0f\n", n, typat[n], znucltypat[i]);
    }
  }

  wtk.resize(nkpt);
  fread (&(wtk[0]), sizeof(double), nkpt, fin);

  PSPinfo.resize(npsp);
  for (int ipsp=0; ipsp<npsp; ipsp++) {
    // Skip FORTRAN record end and start:
    fread (dummy, sizeof(char), 2*recHeaderLength, fin);
    
    char title[133];
    title[132] = '\0';
    fread (title, sizeof(char), 132, fin);
    PSPinfo[ipsp].title = title;
    fread (&PSPinfo[ipsp].znuclpsp, sizeof(double), 1, fin);
    fread (&PSPinfo[ipsp].zionpsp,  sizeof(double), 1, fin);
    fread (&PSPinfo[ipsp].pspso,    sizeof(int),    1, fin);
    fread (&PSPinfo[ipsp].pspdat,   sizeof(int),    1, fin);
    fread (&PSPinfo[ipsp].pspcod,   sizeof(int),    1, fin);
    fread (&PSPinfo[ipsp].pspxc,    sizeof(int),    1, fin);
    fread (&PSPinfo[ipsp].lmn_size, sizeof(int),    1, fin);
  }
  // Skip FORTRAN record end and start:
  fread (dummy, sizeof(char), 2*recHeaderLength, fin);
  fread (&residm, sizeof(double), 1, fin);

  xred.resize(natom);
  fread (&(xred[0][0]), sizeof(double), 3*natom, fin);

  fread (&etotal, sizeof(double), 1, fin);
  fread (&fermie, sizeof(double), 1, fin);
  if (print) 
    fprintf (stderr, "Total energy = %1.12f, fermi energy = %1.12f\n\n",
	     etotal, fermie);
}

struct Int3Less
{
  inline bool operator()(Int3 a, Int3 b) 
  {
    for (int i=0; i<3; i++)
      if (a[i] < b[i])	        return true;
      else if (a[i] > b[i])	return false;
    return false;
  }
};


bool
OrbitalSetClass::Read_ABINIT_WFK (string fname)
{
  if (Comm.MyProc() != 0)
    return true;

  FILE *fin = fopen (fname.c_str(), "r");
  assert (fin != NULL);

  ABINIT_header header;
  header.Read(fin, true);

  ///////////////////
  // Setup lattice //
  ///////////////////
  Mat3 Aprim, Asuper;
  for (int i=0; i<3; i++)
    for (int j=0; j<3; j++)
      Aprim(i,j) = header.rprimd[i][j];
  Asuper = TileMatrix*Aprim;
  PrimCell.SetLattice (Aprim);
  SuperCell.SetLattice (Asuper);
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

  /////////////////////////
  // Setup ion positions //
  /////////////////////////
  PrimCell.IonPos.resize(header.natom);
  PrimCell.AtomTypes.resize(header.natom);
  PrimCell.Zion.resize(header.natom);
  for (int i=0; i<header.natom; i++) {
    PrimCell.IonPos(i) = PrimCell.Lattice.u2r(header.xred[i]);
    PrimCell.AtomTypes(i) = (int)round(header.znucltypat[header.typat[i]-1]);
    PrimCell.Zion(i) = header.PSPinfo[header.typat[i]-1].zionpsp;
  }

  fprintf (stderr,   "Atom Type           Position:\n");
  for (int atom=0; atom<header.natom; atom++) {
    fprintf (stderr, "   %2d      %9.5f %9.5f %9.5f\n", 
	     PrimCell.AtomTypes(atom),
	     PrimCell.IonPos(atom)[0], PrimCell.IonPos(atom)[1], 
	     PrimCell.IonPos(atom)[2]);
  }
  fprintf (stderr, "\n");

  // First pass:  Read G-vectors for all different k-points, and then
  // take the union of the G-vector sets
  int bandIndex = 0;
  char dummy[1000];
  std::map<Int3,int,Int3Less> gMap;
  long wfStart = ftell(fin);
  for (int isppol=0; isppol<header.nsppol; isppol++) {
    for (int ikpt=0; ikpt < header.nkpt; ikpt++) {
      // Skip FORTRAN record end and start:
      fread (dummy, sizeof(char), 2*header.recHeaderLength, fin);
      int npw, nspinor, nband;
      fread (&npw,     sizeof(int), 1, fin);
      fread (&nspinor, sizeof(int), 1, fin);
      fread (&nband,   sizeof(int), 1, fin);

      // Skip FORTRAN record end and start:
      fread (dummy, sizeof(char), 2*header.recHeaderLength, fin);

      vector<Int3> kg(npw);      
      fread (&(kg[0][0]), sizeof(int), 3*npw, fin);
      
      long skip = (2*header.recHeaderLength + 2*sizeof(double)*nband +
		   2*nband*header.recHeaderLength + 
		   sizeof(double)*2*npw*nspinor*nband);

      // Skip the rest of the data
      fseek (fin, skip, SEEK_CUR);
      
      for (int i=0; i<npw; i++) {
	map<Int3,int,Int3Less>::iterator iter = gMap.find(kg[i]);
	if (iter == gMap.end()) {
	  int n = gMap.size();
	  gMap[kg[i]] = n;
	}
      }
    }
  }

  // Setup the G-vectors
  int numG = gMap.size();
  GVecsArray.resize(numG);
  GInts.resize(numG);
  int ig=0;
  std::map<Int3,int,Int3Less>::iterator iter;
  for (iter=gMap.begin(); iter!=gMap.end(); iter++) {
    Int3 gInt = -1 * iter->first;

    int ig    = iter->second;
    GInts(ig) = gInt;
    GVecsArray(ig) = ((double)gInt[0] * PrimCell.Lattice.b(0) +
		      (double)gInt[1] * PrimCell.Lattice.b(1) +
		      (double)gInt[2] * PrimCell.Lattice.b(2));
    ig++;
  }

  // Now, rewind to beginning of WF file.
  fseek (fin, wfStart, SEEK_SET);
  
  // Resize the orbitals
  int numk = header.nkpt;
  int numBands = 0;
  for (int ik=0; ik<numk; ik++)
    numBands = max(numBands, header.nband[ik]);
  PrimOrbitals[0].resize(numk, numBands);
  for (int ik=0; ik<numk; ik++) {
    for (int band=0; band<numBands; band++) {
      PrimOrbitals[0](ik,band) = new OrbitalClass;
      PrimOrbitals[0](ik,band)->SetCell (PrimCell);
    }
  }

  SpinPolarized = false;
  ECut = header.ecut;
  // If we are spin-polarized
  cerr << "header.nsppol = " << header.nsppol << endl;
  if (header.nsppol > 1) {
    SpinPolarized = true;
    cerr << "Spin-polarized calculation.\n";
    PrimOrbitals[1].resize(numk, numBands);
    for (int ik=0; ik<numk; ik++) {
      for (int band=0; band<numBands; band++) {
	PrimOrbitals[1](ik,band) = new OrbitalClass;
	PrimOrbitals[1](ik,band)->SetCell (PrimCell);
      }
    }
  }
    

  // Normalize the k-point weights
  double wtk_norm = 0.0;
  for (int ikpt=0; ikpt < header.nkpt; ikpt++)
    wtk_norm += header.wtk[ikpt];
  for (int ikpt=0; ikpt < header.nkpt; ikpt++)
    header.wtk[ikpt] /= wtk_norm;

  // Now, actually read the orbitals
  zVec tmpCoefs(numG);
  vector<Int3> kg(numG);
  TinyVector<double,2> numElectrons(0.0, 0.0);

  for (int isppol=0; isppol<header.nsppol; isppol++) {
    for (int ikpt=0; ikpt < header.nkpt; ikpt++) {
      Vec3 twist = -1.0*header.kpt[ikpt];
      //Vec3 twist = header.kpt[ikpt];

      // Skip FORTRAN record end and start:
      fread (dummy, sizeof(char), 2*header.recHeaderLength, fin);
      int npw, nspinor, nband;
      fread (&npw,     sizeof(int), 1, fin);
      assert (npw <= numG);
      fread (&nspinor, sizeof(int), 1, fin);
      fread (&nband,   sizeof(int), 1, fin);
      if (nspinor > 1) 
	cerr << "nspinor > 1!!!!!!!!!!!!!!!!!!!\n";

      // Skip FORTRAN record end and start:
      fread (dummy, sizeof(char), 2*header.recHeaderLength, fin);

      fread (&(kg[0][0]), sizeof(int), 3*npw, fin);

      // Skip FORTRAN record end and start:
      fread (dummy, sizeof(char), 2*header.recHeaderLength, fin);

      // Read eigenvalues
      vector<double> eigen(nband), occ(nband);
      fread (&(eigen[0]), sizeof(double), nband, fin);
      fread (&(  occ[0]), sizeof(double), nband, fin);
      for (int band=0; band<nband; band++) {
	numElectrons[isppol] += header.wtk[ikpt] * occ[band];
	PrimOrbitals[isppol](ikpt,band)->SetEigVal (eigen[band]);
	PrimOrbitals[isppol](ikpt,band)->SetOccupancy (occ[band]);
	PrimOrbitals[isppol](ikpt,band)->SetLabels (isppol, ikpt, band);
	PrimOrbitals[isppol](ikpt,band)->SetTwist (twist);
	PrimOrbitals[isppol](ikpt,band)->SetCell (PrimCell);
	zVec &coefs(*(new zVec(numG)));
	coefs = complex<double>();
	// Skip FORTRAN record end and start:
	fread (dummy, sizeof(char), 2*header.recHeaderLength, fin);
	fread (tmpCoefs.data(), sizeof(double), 2*npw*nspinor, fin);
	for (int i=0; i<npw; i++) {
	  int iG = gMap[kg[i]];
	  coefs(iG) = tmpCoefs(i);
	}
	PrimOrbitals[isppol](ikpt,band)->SetCoefs(coefs);
      }
    }
  }
  if (header.nsppol == 1)
    numElectrons[0] = numElectrons[1] = 0.5*numElectrons[0];

  // Skip FORTRAN record end and start:
  fread (dummy, sizeof(char), 2*header.recHeaderLength, fin);
  // Read forces
  int nxfh;
  fread(&nxfh, sizeof(int), 1, fin);
  fread (dummy, sizeof(char), 2*header.recHeaderLength, fin);
  cerr << "nxfh = " << nxfh << endl;
  if (nxfh == 1) {
    double dummy[12];
    Array<Vec3,1> xred(header.natom);
    cerr << "Reading forces:\n";
    PrimCell.IonForces.resize(header.natom);
    fread(&(xred(0)[0]), sizeof(double), 3*header.natom, fin);
    fread(dummy, sizeof(double), 12, fin);
    fread(&(PrimCell.IonForces(0)[0]), sizeof(double), 3*header.natom, fin);
    for (int i=0; i<header.natom; i++)
      PrimCell.IonForces(i) = PrimCell.Lattice.r2u(PrimCell.IonForces(i));
  }

  fclose(fin);
  NumElectrons[0] = (int)round(numElectrons[0]);
  NumElectrons[1] = (int)round(numElectrons[1]);
  cerr << "numElectrons = " << numElectrons << endl;

  SetupFFTIndices();
  CreateSuperOrbitals();
  if (UseMultiRep)
    CreateFineOrbitals();
  TileIonPos();

  return true;
}



bool
OrbitalSetClass::Read_ABINIT_WF1 (string fname)
{
  if (Comm.MyProc() != 0)
    return true;

  FILE *fin = fopen (fname.c_str(), "r");
  assert (fin != NULL);

  ABINIT_header header;
  header.Read(fin);

  ///////////////////
  // Setup lattice //
  ///////////////////
  Mat3 Aprim, Asuper;
  for (int i=0; i<3; i++)
    for (int j=0; j<3; j++)
      Aprim(i,j) = header.rprimd[i][j];
  Asuper = TileMatrix*Aprim;
  PrimCell.SetLattice (Aprim);
  SuperCell.SetLattice (Asuper);
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

  /////////////////////////
  // Setup ion positions //
  /////////////////////////
  PrimCell.IonPos.resize(header.natom);
  PrimCell.AtomTypes.resize(header.natom);
  for (int i=0; i<header.natom; i++) {
    PrimCell.IonPos(i) = PrimCell.Lattice.u2r(header.xred[i]);
    PrimCell.AtomTypes(i) = (int)round(header.znucltypat[header.typat[i]-1]);
  }

  fprintf (stderr,   "Atom Type           Position:\n");
  for (int atom=0; atom<header.natom; atom++) {
    fprintf (stderr, "   %2d      %9.5f %9.5f %9.5f\n", 
	     PrimCell.AtomTypes(atom),
	     PrimCell.IonPos(atom)[0], PrimCell.IonPos(atom)[1], 
	     PrimCell.IonPos(atom)[2]);
  }
  fprintf (stderr, "\n");

  // First pass:  Read G-vectors for all different k-points, and then
  // take the union of the G-vector sets
  int bandIndex = 0;
  char dummy[1000];
  std::map<Int3,int,Int3Less> gMap;
  long wfStart = ftell(fin);
  for (int isppol=0; isppol<header.nsppol; isppol++) {
    for (int ikpt=0; ikpt < header.nkpt; ikpt++) {
      // Skip FORTRAN record end and start:
      fread (dummy, sizeof(char), 2*header.recHeaderLength, fin);
      int npw, nspinor, nband;
      fread (&npw,     sizeof(int), 1, fin);
      fread (&nspinor, sizeof(int), 1, fin);
      fread (&nband,   sizeof(int), 1, fin);

      // Skip FORTRAN record end and start:
      fread (dummy, sizeof(char), 2*header.recHeaderLength, fin);

      vector<Int3> kg(npw);      
      fread (&(kg[0][0]), sizeof(int), 3*npw, fin);
      
      long skip = (4*nband*header.recHeaderLength + 
		   2*nband*nband*sizeof(double) + 
		   sizeof(double)*2*npw*nspinor*nband);

      // Skip the rest of the data
      fseek (fin, skip, SEEK_CUR);
      
      for (int i=0; i<npw; i++) {
	map<Int3,int,Int3Less>::iterator iter = gMap.find(kg[i]);
	if (iter == gMap.end()) {
	  int n = gMap.size();
	  gMap[kg[i]] = n;
	}
      }
    }
  }

  // Setup the G-vectors
  int numG = gMap.size();
  GVecsArray.resize(numG);
  GInts.resize(numG);
  int ig=0;
  std::map<Int3,int,Int3Less>::iterator iter;
  for (iter=gMap.begin(); iter!=gMap.end(); iter++) {
    Int3 gInt = -1 * iter->first;

    int ig    = iter->second;
    GInts(ig) = gInt;
    GVecsArray(ig) = ((double)gInt[0] * PrimCell.Lattice.b(0) +
		      (double)gInt[1] * PrimCell.Lattice.b(1) +
		      (double)gInt[2] * PrimCell.Lattice.b(2));
    ig++;
  }

  // Now, rewind to beginning of WF file.
  fseek (fin, wfStart, SEEK_SET);
  
  // Resize the orbitals
  int numk = header.nkpt;
  int numBands = 0;
  for (int ik=0; ik<numk; ik++)
    numBands = max(numBands, header.nband[ik]);
  PrimOrbitals[0].resize(numk, numBands);
  for (int ik=0; ik<numk; ik++) {
    for (int band=0; band<numBands; band++) {
      PrimOrbitals[0](ik,band) = new OrbitalClass;
      PrimOrbitals[0](ik,band)->SetCell (PrimCell);
    }
  }

  SpinPolarized = false;
  ECut = header.ecut;
  // If we are spin-polarized
  if (header.nsppol > 1) {
    SpinPolarized = true;
    cerr << "Spin-polarized calculation.\n";
    PrimOrbitals[1].resize(numk, numBands);
    for (int ik=0; ik<numk; ik++) {
      for (int band=0; band<numBands; band++) {
	PrimOrbitals[1](ik,band) = new OrbitalClass;
	PrimOrbitals[1](ik,band)->SetCell (PrimCell);
      }
    }
  }
    

  // Normalize the k-point weights
  double wtk_norm = 0.0;
  for (int ikpt=0; ikpt < header.nkpt; ikpt++)
    wtk_norm += header.wtk[ikpt];
  for (int ikpt=0; ikpt < header.nkpt; ikpt++)
    header.wtk[ikpt] /= wtk_norm;

  // Now, actually read the orbitals
  zVec tmpCoefs(numG);
  vector<Int3> kg(numG);
  TinyVector<double,2> numElectrons(0.0, 0.0);

  for (int isppol=0; isppol<header.nsppol; isppol++) {
    for (int ikpt=0; ikpt < header.nkpt; ikpt++) {
      Vec3 twist = -1.0*header.kpt[ikpt];
      //Vec3 twist = header.kpt[ikpt];

      // Skip FORTRAN record end and start:
      fread (dummy, sizeof(char), 2*header.recHeaderLength, fin);
      int npw, nspinor, nband;
      fread (&npw,     sizeof(int), 1, fin);
      assert (npw <= numG);
      fread (&nspinor, sizeof(int), 1, fin);
      fread (&nband,   sizeof(int), 1, fin);
      if (nspinor > 1) 
	cerr << "nspinor > 1!!!!!!!!!!!!!!!!!!!\n";

      // Skip FORTRAN record end and start:
      fread (dummy, sizeof(char), 2*header.recHeaderLength, fin);

      fread (&(kg[0][0]), sizeof(int), 3*npw, fin);

      // Read eigenvalues
      Array<complex<double>,2> eigen(nband,nband);
      for (int band=0; band<nband; band++) {
	// Skip FORTRAN record end and start:
	fread (dummy, sizeof(char), 2*header.recHeaderLength, fin);
	fread (&(eigen(band,0)), sizeof(complex<double>), nband, fin);
	PrimOrbitals[isppol](ikpt,band)->SetEigVal (eigen(band,0).real());
	PrimOrbitals[isppol](ikpt,band)->SetOccupancy (0.0);
	PrimOrbitals[isppol](ikpt,band)->SetLabels (isppol, ikpt, band);
	PrimOrbitals[isppol](ikpt,band)->SetTwist (twist);
	PrimOrbitals[isppol](ikpt,band)->SetCell (PrimCell);
	zVec &coefs(*(new zVec(numG)));
	coefs = complex<double>();
	// Skip FORTRAN record end and start:
	fread (dummy, sizeof(char), 2*header.recHeaderLength, fin);
	fread (tmpCoefs.data(), sizeof(double), 2*npw*nspinor, fin);
	for (int i=0; i<npw; i++) {
	  int iG = gMap[kg[i]];
	  coefs(iG) = tmpCoefs(i);
	}
	PrimOrbitals[isppol](ikpt,band)->SetCoefs(coefs);
      }
    }
  }
  if (header.nsppol == 1)
    numElectrons[0] = numElectrons[1] = 0.5*numElectrons[0];

  // Skip FORTRAN record end and start:
  fread (dummy, sizeof(char), 2*header.recHeaderLength, fin);
  // Read forces
  int nxfh;
  fread(&nxfh, sizeof(int), 1, fin);
  fread (dummy, sizeof(char), 2*header.recHeaderLength, fin);
  cerr << "nxfh = " << nxfh << endl;
  if (nxfh == 1) {
    double dummy[12];
    Array<Vec3,1> xred(header.natom);
    cerr << "Reading forces:\n";
    PrimCell.IonForces.resize(header.natom);
    fread(&(xred(0)[0]), sizeof(double), 3*header.natom, fin);
    fread(dummy, sizeof(double), 12, fin);
    fread(&(PrimCell.IonForces(0)[0]), sizeof(double), 3*header.natom, fin);
    for (int i=0; i<header.natom; i++)
      PrimCell.IonForces(i) = PrimCell.Lattice.r2u(PrimCell.IonForces(i));
  }

  fclose(fin);
  NumElectrons[0] = (int)round(numElectrons[0]);
  NumElectrons[1] = (int)round(numElectrons[1]);
  cerr << "numElectrons = " << numElectrons << endl;

  SetupFFTIndices();
  CreateSuperOrbitals();
  TileIonPos();

  return true;
}





bool
OrbitalSetClass::Read_ABINIT_First_Order (string prefix)
{
  // First, see how many files ther are
  int n = 1;
  bool done = false;
  while (!done) {
    ostringstream fname;
    fname << prefix << n;
    FILE *fin = fopen (fname.str().c_str(), "r");
    if (fin != NULL) {
      fclose (fin);
      n++;
    }
    else 
      done = true;
  }
  

  int numFiles = n - 1;
  PrimFirstOrder.resize(numFiles);

  if (Comm.MyProc() != 0)
    return true;

  cerr << "Found " << numFiles << " first-order wave function files.\n";

  for (n=0; n< numFiles; n++) {
    ostringstream fname;
    fname << prefix << n+1;
    FILE *fin = fopen (fname.str().c_str(), "r");
    assert (fin != NULL);

    ABINIT_header header;
    header.Read(fin, false);
    TinyVector<Array<OrbitalClass*,2>,2>& orbitals = 
      PrimFirstOrder[n];
        
    // First pass:  Read G-vectors for all different k-points, and then
    // take the union of the G-vector sets
    int bandIndex = 0;
    char dummy[1000];
    std::map<Int3,int,Int3Less> gMap;
    long wfStart = ftell(fin);
    for (int isppol=0; isppol<header.nsppol; isppol++) {
      for (int ikpt=0; ikpt < header.nkpt; ikpt++) {
	// Skip FORTRAN record end and start:
	fread (dummy, sizeof(char), 2*header.recHeaderLength, fin);
	int npw, nspinor, nband;
	fread (&npw,     sizeof(int), 1, fin);
	fread (&nspinor, sizeof(int), 1, fin);
	fread (&nband,   sizeof(int), 1, fin);
	
	// Skip FORTRAN record end and start:
	fread (dummy, sizeof(char), 2*header.recHeaderLength, fin);
	
	vector<Int3> kg(npw);      
	fread (&(kg[0][0]), sizeof(int), 3*npw, fin);
	
	long skip = (4*nband*header.recHeaderLength + 
		     2*nband*nband*sizeof(double) + 
		     sizeof(double)*2*npw*nspinor*nband);
	
	// Skip the rest of the data
	fseek (fin, skip, SEEK_CUR);
	
	for (int i=0; i<npw; i++) {
	  map<Int3,int,Int3Less>::iterator iter = gMap.find(kg[i]);
	  if (iter == gMap.end()) {
	    int n = gMap.size();
	    gMap[kg[i]] = n;
	  }
	}
      }
    }
    
    // Setup the G-vectors
    int numG = gMap.size();
    GVecsArray.resize(numG);
    GInts.resize(numG);
    int ig=0;
    std::map<Int3,int,Int3Less>::iterator iter;
    for (iter=gMap.begin(); iter!=gMap.end(); iter++) {
      Int3 gInt = -1 * iter->first;
      
      int ig    = iter->second;
      GInts(ig) = gInt;
      GVecsArray(ig) = ((double)gInt[0] * PrimCell.Lattice.b(0) +
			(double)gInt[1] * PrimCell.Lattice.b(1) +
			(double)gInt[2] * PrimCell.Lattice.b(2));
      ig++;
    }
    
    // Now, rewind to beginning of WF file.
    fseek (fin, wfStart, SEEK_SET);
    
    // Resize the orbitals
    int numk = header.nkpt;
    int numBands = 0;
    for (int ik=0; ik<numk; ik++)
      numBands = max(numBands, header.nband[ik]);
    orbitals[0].resize(numk, numBands);
    for (int ik=0; ik<numk; ik++) {
      for (int band=0; band<numBands; band++) {
	orbitals[0](ik,band) = new OrbitalClass;
	orbitals[0](ik,band)->SetCell (PrimCell);
      }
    }
    
    SpinPolarized = false;
    ECut = header.ecut;
    // If we are spin-polarized
    if (header.nsppol > 1) {
      SpinPolarized = true;
      cerr << "Spin-polarized calculation.\n";
      orbitals[1].resize(numk, numBands);
      for (int ik=0; ik<numk; ik++) {
	for (int band=0; band<numBands; band++) {
	  orbitals[1](ik,band) = new OrbitalClass;
	  orbitals[1](ik,band)->SetCell (PrimCell);
	}
      }
    }
    

    // Normalize the k-point weights
    double wtk_norm = 0.0;
    for (int ikpt=0; ikpt < header.nkpt; ikpt++)
      wtk_norm += header.wtk[ikpt];
    for (int ikpt=0; ikpt < header.nkpt; ikpt++)
      header.wtk[ikpt] /= wtk_norm;

    // Now, actually read the orbitals
    zVec tmpCoefs(numG);
    vector<Int3> kg(numG);
    TinyVector<double,2> numElectrons(0.0, 0.0);

    for (int isppol=0; isppol<header.nsppol; isppol++) {
      for (int ikpt=0; ikpt < header.nkpt; ikpt++) {
	Vec3 twist = -1.0*header.kpt[ikpt];
	//Vec3 twist = header.kpt[ikpt];

	// Skip FORTRAN record end and start:
	fread (dummy, sizeof(char), 2*header.recHeaderLength, fin);
	int npw, nspinor, nband;
	fread (&npw,     sizeof(int), 1, fin);
	assert (npw <= numG);
	fread (&nspinor, sizeof(int), 1, fin);
	fread (&nband,   sizeof(int), 1, fin);
	if (nspinor > 1) 
	  cerr << "nspinor > 1!!!!!!!!!!!!!!!!!!!\n";

	// Skip FORTRAN record end and start:
	fread (dummy, sizeof(char), 2*header.recHeaderLength, fin);

	fread (&(kg[0][0]), sizeof(int), 3*npw, fin);

	// Read eigenvalues
	Array<complex<double>,2> eigen(nband,nband);
	for (int band=0; band<nband; band++) {
	  // Skip FORTRAN record end and start:
	  fread (dummy, sizeof(char), 2*header.recHeaderLength, fin);
	  fread (&(eigen(band,0)), sizeof(complex<double>), nband, fin);
	  orbitals[isppol](ikpt,band)->SetEigVal (eigen(band,0).real());
	  orbitals[isppol](ikpt,band)->SetOccupancy (0.0);
	  orbitals[isppol](ikpt,band)->SetLabels (isppol, ikpt, band);
	  orbitals[isppol](ikpt,band)->SetTwist (twist);
	  orbitals[isppol](ikpt,band)->SetCell (PrimCell);
	  zVec &coefs(*(new zVec(numG)));
	  coefs = complex<double>();
	  // Skip FORTRAN record end and start:
	  fread (dummy, sizeof(char), 2*header.recHeaderLength, fin);
	  fread (tmpCoefs.data(), sizeof(double), 2*npw*nspinor, fin);
	  for (int i=0; i<npw; i++) {
	    int iG = gMap[kg[i]];
	    coefs(iG) = tmpCoefs(i);
	  }
	  orbitals[isppol](ikpt,band)->SetCoefs(coefs);
	  int atom, dir;
	  atom = (header.pertcase-1)/3;
	  dir  = (header.pertcase-1)%3;
	  orbitals[isppol](ikpt,band)->SetDeriv(atom,dir);
	}
      }
    }
    if (header.nsppol == 1)
      numElectrons[0] = numElectrons[1] = 0.5*numElectrons[0];

    fclose(fin);
    NumElectrons[0] = (int)round(numElectrons[0]);
    NumElectrons[1] = (int)round(numElectrons[1]);


    // SetupFFTIndices();
    // CreateSuperOrbitals();
    // TileIonPos();
  }
  // Set orbitals to indices
  for (int ideriv=0; ideriv<PrimFirstOrder.size(); ideriv++) 
    for (int spin=0; spin<2; spin++) {
      for (int ik=0; ik<PrimFirstOrder[ideriv][spin].extent(0); ik++)
	for (int iband=0; iband<PrimFirstOrder[ideriv][spin].extent(1); iband++)
	  PrimFirstOrder[ideriv][spin](ik,iband)->SetFFTIndices(PrimFFTIndices);
    }
  
  return true;
}

inline complex<double> zdot(zVec &a, zVec&b)
{
  complex<double> sum = 0.0;
  for (int i=0; i<a.size(); i++)
    sum += a(i)*conj(b(i));
  return sum;
}

inline double dot(zVec &a, zVec&b)
{
  double sum = 0.0;
  for (int i=0; i<a.size(); i++)
    sum += real(a(i)*conj(b(i)));
  return sum;
}


bool
OrbitalSetClass::Read_ABINIT_First_Order_FD (string name)
{
  IO::IOSectionClass in;
  assert (in.OpenFile(name));
  
  Array<string,3> fnames;
  double epsilon;
  assert (in.ReadVar("FileNames", fnames));
  assert (in.ReadVar("epsilon", epsilon));

  int numFiles = fnames.size();
  PrimFirstOrder.resize(numFiles/2);

  if (Comm.MyProc() != 0)
    return true;

  cerr << "Found " << numFiles << " first-order wave function files.\n";

  int n = 0;
  for (int atom=0; atom < fnames.extent(0); atom++)
    for (int dim=0; dim < 3; dim++, n++) {
      OrbitalSetClass OrbPlus, OrbMinus;
      OrbPlus.Localized       = OrbMinus.Localized       = Localized;
      OrbPlus.Spline          = OrbMinus.Spline          = Spline;
      OrbPlus.OptimizeCenters = OrbMinus.OptimizeCenters = OptimizeCenters;
      OrbPlus.OptimizeRadii   = OrbMinus.OptimizeRadii   = OptimizeRadii;
      OrbPlus.ShiftOrbitals   = OrbMinus.ShiftOrbitals   = ShiftOrbitals;
      OrbPlus.Real            = OrbMinus.Real            = Real;
      OrbPlus.CheckKE         = OrbMinus.CheckKE         = CheckKE;
      OrbPlus.Orthogonalize   = OrbMinus.Orthogonalize   = Orthogonalize;
      OrbPlus.Truncate        = OrbMinus.Truncate        = Truncate;
      OrbPlus.UseMultiRep     = OrbMinus.UseMultiRep     = UseMultiRep;
      OrbPlus.SetFFTFactor (FFTFactor);   OrbMinus.SetFFTFactor (FFTFactor);
      OrbPlus.SetTileMatrix (TileMatrix); OrbMinus.SetTileMatrix (TileMatrix);

      OrbPlus.Read(fnames(atom,dim,1));
      OrbMinus.Read(fnames(atom,dim,0));

      TinyVector<Array<OrbitalClass*,2>,2>& orbitals = PrimFirstOrder[n];
      
      // Resize the orbitals
      int numk     = PrimOrbitals[0].extent(0);
      int numBands = PrimOrbitals[0].extent(1);
      int numSpins = SpinPolarized ? 2 : 1;
      int numG = PrimCell.GVecs.size();
      assert (OrbPlus.PrimCell.GVecs.size() == numG);
      assert (OrbMinus.PrimCell.GVecs.size() == numG);
      assert (OrbPlus.PrimOrbitals[0].shape() == PrimOrbitals[0].shape());
      assert (OrbMinus.PrimOrbitals[0].shape() == PrimOrbitals[0].shape());
      for (int spin=0; spin < numSpins; spin++)  {
	orbitals[spin].resize(numk, numBands);
	for (int ik=0; ik<numk; ik++) 
	  for (int band=0; band<numBands; band++) {
	    orbitals[spin](ik,band) = new OrbitalClass;
	    zVec *coefs = new zVec(numG);
	    
	    zVec &c      = PrimOrbitals[spin](ik,band)->GetCoefs();
	    zVec &cPlus  = OrbPlus.PrimOrbitals[spin] (ik,band)->GetCoefs();
	    zVec &cMinus = OrbMinus.PrimOrbitals[spin] (ik,band)->GetCoefs();
	    if (dot(cPlus,c) < 0.0)   cPlus  = -1.0*cPlus;
	    if (dot(cMinus,c) < 0.0)  cMinus = -1.0*cMinus;
	    complex<double> phase = zdot(cPlus, c)/
	      sqrt(dot(cPlus,cPlus)*dot(c,c));
	    //	    cPlus *= conj(phase);
	    complex<double> phase2 = zdot (cPlus, c)/
	      sqrt(dot(cPlus,cPlus) * dot(c,c));

	    *coefs = (1.0/(2.0*epsilon))* (cPlus - cMinus);
	    
	    orbitals[spin](ik,band)->SetCoefs(*coefs);
	    orbitals[spin](ik,band)->SetDeriv(atom, dim);
	    orbitals[spin](ik,band)->SetCell(PrimCell);
	    orbitals[spin](ik,band)->SetFFTIndices(PrimFFTIndices);
	    orbitals[spin](ik,band)->SetOccupancy(0.0);
	    orbitals[spin](ik,band)->SetLabels(spin, ik, band);
	    orbitals[spin](ik,band)->SetTwist(PrimOrbitals[spin](ik,band)->GetTwist());
	  }
      }
      
    }
  return true;
}





bool
OrbitalSetClass::Read_ABINIT_DEN (string fname)
{
  PrimDensity.SetLattice  (PrimCell.Lattice.GetDirect());
  SuperDensity.SetLattice (SuperCell.Lattice.GetDirect());

  FILE *fin = fopen (fname.c_str(), "r");
  if (fin == NULL) {
    cerr << "Could not open ABINIT density file \"" << fname 
	 << "\".  Exitting.\n";
    abort();
  }

  ABINIT_header header;
  header.Read(fin);
  header.SkipRecHeader(fin);
  header.SkipRecHeader(fin);
  int nx = header.ngfft[0];
  int ny = header.ngfft[1];
  int nz = header.ngfft[2];
  // The first record is always the total density
  Array<double,3> rhoTotalTran(nz,ny,nx), rhoTotal(nx,ny,nz);
  size_t num_read = fread(rhoTotalTran.data(), sizeof(double), nx*ny*nz, fin);
  if (num_read != (nx*ny*nz)) {
    cerr << "Error reading ABINIT density data.\n";
    abort();
  }
  for (int ix=0; ix<nx; ix++)
    for (int iy=0; iy<ny; iy++)
      for (int iz=0; iz<nz; iz++) 
	rhoTotal(ix,iy,iz) = rhoTotalTran(iz,iy,ix);

  header.SkipRecHeader(fin);
  header.SkipRecHeader(fin);

  // If there is more than one density record, then the next is the
  // up-electron density.
  if (header.nspden > 1) {
    Array<double,3> rhoUpTran(nz,ny,nz), rhoUp(nx,ny,nz), rhoDown(nx,ny,nz);
    fread(rhoUpTran.data(), sizeof(double), nx*ny*nz, fin);
    for (int ix=0; ix<nx; ix++)
      for (int iy=0; iy<ny; iy++)
	for (int iz=0; iz<nz; iz++) 
	  rhoUp(ix,iy,iz) = rhoUpTran(iz,iy,ix);
    
    rhoDown = rhoTotal - rhoUp;
    PrimDensity.Set (rhoUp, rhoDown, header.ecut);
  }
  else 
    PrimDensity.Set (rhoTotal, header.ecut);

  fclose (fin);

  return true;	  
}


bool
OrbitalSetClass::Read_ABINIT_POT (string fname)
{
  PrimVHXC.SetLattice  (PrimCell.Lattice.GetDirect());
  SuperVHXC.SetLattice (SuperCell.Lattice.GetDirect());

  FILE *fin = fopen (fname.c_str(), "r");
  if (fin == NULL) {
    cerr << "Could not open ABINIT potential file \"" << fname 
	 << "\".  Exitting.\n";
    abort();
  }

  ABINIT_header header;
  header.Read(fin);
  header.SkipRecHeader(fin);
  header.SkipRecHeader(fin);
  int nx = header.ngfft[0];
  int ny = header.ngfft[1];
  int nz = header.ngfft[2];

  // cerr << "Potential FFT grid is " << nx << "x" << ny << "x" << nz << endl;

  // The first record is always the up potential
  Array<double,3> potUpTran(nz,ny,nx), potUp(nx,ny,nz);
  size_t num_read = fread(potUpTran.data(), sizeof(double), nx*ny*nz, fin);
  if (num_read != (nx*ny*nz)) {
    cerr << "Error reading ABINIT potential data.\n";
    abort();
  }
  for (int ix=0; ix<nx; ix++)
    for (int iy=0; iy<ny; iy++)
      for (int iz=0; iz<nz; iz++) 
	potUp(ix,iy,iz) = potUpTran(iz,iy,ix);

  header.SkipRecHeader(fin);
  header.SkipRecHeader(fin);

  // If there is more than one potential record, then the next is the
  // down-electron potential.
  if (header.nspden > 1) {
    Array<double,3> potDownTran(nz,ny,nz), potDown(nx,ny,nz);
    fread(potDownTran.data(), sizeof(double), nx*ny*nz, fin);
    for (int ix=0; ix<nx; ix++)
      for (int iy=0; iy<ny; iy++)
	for (int iz=0; iz<nz; iz++) 
	  potDown(ix,iy,iz) = potDownTran(iz,iy,ix);
    
    PrimVHXC.Set (potUp, potDown, header.ecut);
  }
  else 
    PrimVHXC.Set (potUp, header.ecut);

  fclose (fin);

  return true;	  
}
