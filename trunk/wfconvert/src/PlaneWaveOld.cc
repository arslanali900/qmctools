/// This function creates a larger, single k-point set of orbitals by
/// unfolding the plane-wave coefficients from different k-points.
/// This may be accomplished by transfering the plane wave
/// coefficients to a finer G-vector mesh and then filling in zeros
/// between non-zero coefficents.  
/// The new parameter, foldFactor, allow incomplete unfolding of the
/// k-point mesh.  If one had a 4x4x4 k-point mesh, for example, it
/// could be unfolded into a 2x2x2 supercell on a 2x2x2 k-point mesh.
void
PlaneWaveSystem::Unfold(int spin, Int3 foldFactor)
{
  // Just rename to save typing
  Int3 &ff = foldFactor;

  if (Comm.MyProc() != 0)
    return;
  ///////////////////////////////////////////
  // First, make all twist vectors postive //
  ///////////////////////////////////////////
//   for (int ki=0; ki<Orbitals[spin].extent(0); ki++) {
//     Vec3 twist = Orbitals[spin](ki,0)->GetTwist();
//     for (int i=0; i<3; i++)
//       twist[i] -= floor(twist[i]);
//     Orbitals[spin](ki,0)->SetTwist(twist);
//   }

  /////////////////////////////////////////////
  // Check to see if k-points form a lattice //
  /////////////////////////////////////////////
  Vec3 minTwist(1.0e12, 1.0e12, 1.0e12), maxTwist(-1.0e12, -1.0e12, -1.0e12);
  Vec3 mink    (1.0e12, 1.0e12, 1.0e12);
  for (int ki=0; ki<Orbitals[spin].extent(0); ki++) {
    Vec3 twist = Orbitals[spin](ki,0)->GetPrimTwist();
    Vec3 supertwist((double)ff[0]*twist[0], (double)ff[1]*twist[1], (double)ff[2]*twist[2]);
    cerr << "twist = " << twist << "   supertwist = " << supertwist << endl;
    Vec3 k     = Orbitals[spin](ki,0)->GetPrimk();
    minTwist[0] = min(minTwist[0], twist[0]); 
    maxTwist[0] = max(maxTwist[0], twist[0]);
    minTwist[1] = min(minTwist[1], twist[1]); 
    maxTwist[1] = max(maxTwist[1], twist[1]);
    minTwist[2] = min(minTwist[2], twist[2]); 
    maxTwist[2] = max(maxTwist[2], twist[2]);
    mink[0]     = min(k[0],mink[0]);
    mink[1]     = min(k[1],mink[1]);
    mink[2]     = min(k[2],mink[2]);
  }
  // The difference between maxTwist and minTwist should be of the
  // form (n-1)/n.  Therefore, we determine n by
  Vec3 nf;
  nf[0] = -1.0/((maxTwist[0]-minTwist[0]) -1.0);
  nf[1] = -1.0/((maxTwist[1]-minTwist[1]) -1.0);
  nf[2] = -1.0/((maxTwist[2]-minTwist[2]) -1.0);

  // Make sure they are close to integers
  assert (fabs(nf[0] - round(nf[0]))<1.0e-6);
  assert (fabs(nf[1] - round(nf[1]))<1.0e-6);
  assert (fabs(nf[2] - round(nf[2]))<1.0e-6);
  

  Int3 n ((int)round(nf[0]), (int)round(nf[1]), (int)round(nf[2])), newkmesh;

  // Now, make sure we have all the k-points in the lattice
  Vec3 twist;
  for (int ix=0; ix<n[0]; ix++) 
    for (int iy=0; iy<n[1]; iy++)
      for (int iz=0; iz<n[2]; iz++){
	twist[0] = 
	  minTwist[0] + (double)ix/(double)(n[0]-1)*(maxTwist[0]-minTwist[0]);
	twist[1] = 
	  minTwist[1] + (double)iy/(double)(n[1]-1)*(maxTwist[1]-minTwist[1]);
	twist[2] = 
	  minTwist[2] + (double)iz/(double)(n[2]-1)*(maxTwist[2]-minTwist[2]);
	bool twistFound = false;
	for (int ik=0; ik<Orbitals[spin].extent(0); ik++) {
	  Vec3 diff = Orbitals[spin](ik,0)->GetPrimTwist()-twist;
	  if (dot(diff,diff)<1.0e-8) {
	    twistFound = true;
	    for (int band=0; band<Orbitals[spin].extent(1); band++)
	      Orbitals[spin](ik, band)->SetIndex (Int3(ix,iy,iz));
	  }
	}
	if (!twistFound) {
	  fprintf (stderr, "Missing twist vector (%8.4f, %8.4f, %8.4f) "
		   "in unfolding process.\n", twist[0], twist[1], twist[2]);
	  abort();
	}
      }

  if (ff == Int3(0,0,0))
    ff = n;
  else {
    if ((n[0] % ff[0]) != 0 ||
	(n[1] % ff[1]) != 0 ||
	(n[2] % ff[2]) != 0 ) {
      cerr << "The k-point mesh is not a multiple of the foldFactor.\n"
	   << "The k-point mesh is " << n[0] << "x" << n[1] << "x" << n[2] 
	   << ".  The foldFactor is " << ff[0] << "x" << ff[1] << "x" << ff[2] << ".\n";
      abort();
    }
  }
  newkmesh[0] = n[0]/ff[0];
  newkmesh[1] = n[1]/ff[1];
  newkmesh[2] = n[2]/ff[2];

  // The twist vector of the new orbital will be the minimum twist vector.
  Vec3 newTwist = minTwist;

  
  // We have a correct lattice!  Now unfold!
  perr << "Unfolding k-points into a " 
       << ff[0] << "x" << ff[1] << "x" << ff[2] << " supercell.\n";
  perr << "Each supercell will contain a " 
       << newkmesh[0] << "x" << newkmesh[1] << "x" << newkmesh[2] << " k-point mesh.\n";

  // Replication ion positions
  Array<Vec3,1> newIons(   IonPos.size()*ff[0]*ff[1]*ff[2]);
  Array<int,1> newTypes(AtomTypes.size()*ff[0]*ff[1]*ff[2]);
  int newIndex=0;
  for (int i0=0; i0<ff[0]; i0++) 
    for (int i1=0; i1<ff[1]; i1++)
      for (int i2=0; i2<ff[2]; i2++) 
	for (int k=0; k<IonPos.size(); k++) {
	  newIons(newIndex) = IonPos(k) + 
	    (double)i0 * Vec3(Lattice(0,0), Lattice(0,1), Lattice(0,2)) +
	    (double)i1 * Vec3(Lattice(1,0), Lattice(1,1), Lattice(1,2)) +
	    (double)i2 * Vec3(Lattice(2,0), Lattice(2,1), Lattice(2,2));
	  newTypes(newIndex) = AtomTypes(k);
	  newIndex++;
	}
  IonPos.resize(newIons.size());
  IonPos = newIons;
  AtomTypes.resize(newTypes.size());
  AtomTypes = newTypes;

  // Update number of electrons
  NumElectrons *= (ff[0]*ff[1]*ff[2]);
  // Update the energies
  TotalEnergy    *= (double)(ff[0]*ff[1]*ff[2]);
  KineticEnergy  *= (double)(ff[0]*ff[1]*ff[2]);
  LocalPotEnergy *= (double)(ff[0]*ff[1]*ff[2]);
  NonlocalEnergy *= (double)(ff[0]*ff[1]*ff[2]);
  eeEnergy       *= (double)(ff[0]*ff[1]*ff[2]);
  IonIonEnergy   *= (double)(ff[0]*ff[1]*ff[2]);

  // Adjust lattice and reciprocal lattice
  for (int i=0; i<3; i++)
    for (int j=0; j<3; j++) {
      Lattice(i,j)      *= (double)ff[i];
      RecipLattice(i,j) /= (double)ff[i];
      LatticeInv(i,j)   /= (double)ff[i];
    }
  
  // Create new G-vectors
  Array<Vec3,1> newGVecs(GVecsArray.size()*ff[0]*ff[1]*ff[2]);
  int j=0;
  for (int i=0; i<GVecsArray.size(); i++) {
    for (int i0=0; i0<ff[0]; i0++) {
      Vec3 k0 = (double)i0*
	Vec3(RecipLattice(0,0), RecipLattice(0,1), RecipLattice(0,2));
      for (int i1=0; i1<ff[1]; i1++){
	Vec3 k1 = (double)i1*
	  Vec3(RecipLattice(1,0), RecipLattice(1,1), RecipLattice(1,2));
	for (int i2=0; i2<ff[2]; i2++) {
	  Vec3 k2 = (double)i2*
	    Vec3(RecipLattice(2,0), RecipLattice(2,1), RecipLattice(2,2));
	  Vec3 k = k0 + k1 + k2 /*+ mink */;
	  newGVecs(j) = GVecsArray(i) + k;
	  j++;
	}
      }
    }
  }
  // Reset the G-vectors and FFT box
  PrimGVecs.Set  ( PrimCell.Lattice.GetDirect(), newGVecs, FFTFactor);
  SuperGVecs.Set (SuperCell.Lattice.GetDirect(), newGVecs, FFTFactor);

  // Create new orbitals
  int numk     = Orbitals[spin].extent(0);
  int numbands = Orbitals[spin].extent(1);
  int numfold  = ff[0]*ff[1]*ff[2];
  int newnumk  = numk/numfold;
  int newnumbands = numbands*numfold;

  // Create and initialize new orbital set
  Array<OrbitalClass*,2> newOrbitals(newnumk, newnumbands);
  for (int ki=0; ki<newnumk; ki++)
    for (int bi=0; bi<newnumbands; bi++) {
      newOrbitals(ki, bi) = 
	new OrbitalClass (PrimGVecs, PrimFFT, SuperFFT, PrimkPoints(ki));
      zVec &coefs = newOrbitals(ki, bi)->GetCoefs();
      coefs.resize(newGVecs.size());
      coefs = complex<double>();
    }

  Vec3 gamma(0.0, 0.0, 0.0);
  for (int band=0; band<numbands; band++) {
    for (int ki=0; ki<numk; ki++) {
      Int3 index = Orbitals[spin](ki, band)->GetIndex();      
      Int3 newbIndex (index[0]/newkmesh[0], index[1]/newkmesh[1], index[2]/newkmesh[2]);
      Int3 newkIndex (index[0]%newkmesh[0], index[1]%newkmesh[1], index[2]%newkmesh[2]);
      int newki = newkIndex[2] + newkmesh[2]*(newkIndex[1] + newkmesh[1]*newkIndex[0]);
      int newbi = band*ff[0]*ff[1]*ff[2] +
	          newbIndex[2] +  ff[2]*(newbIndex[1] + ff[1]*newbIndex[0]);
      zVec &coefs = newOrbitals(newki,newbi)->GetCoefs();
      Vec3 newTwist; 

      // This offset reflects the k-vector shift we have included
      // implicitly by offseting the G-vector associated with each
      // plane-wave coefficient.  It must be subtracted off.
      Vec3 twistOffset ((double)newbIndex[0],
			(double)newbIndex[1],
			(double)newbIndex[2]);
      
      newTwist = Orbitals[spin](ki,band)->GetPrimTwist();
      newTwist[0] *= (double)ff[0];
      newTwist[1] *= (double)ff[1];  
      newTwist[2] *= (double)ff[2];
      newTwist -= twistOffset;
      cerr << "Original twist = " 
	   << Orbitals[spin](ki,band)->GetPrimTwist() << endl
	   << "New      twist = " << newTwist << endl;

      int offset = newbIndex[2] + ff[2]*(newbIndex[1] + ff[1]*newbIndex[0]);
//       cerr << "index = " << index << endl;
//       cerr << "newbIndex = " << newbIndex << endl;
//       cerr << "newkIndex = " << newkIndex << endl;
//       cerr << "offset = " << offset << endl;
//       cerr << "numfold = " << numfold << endl;
//       cerr << "GVecsArray.size() = " << GVecsArray.size() << endl;
//       cerr << "newki = " << newki << "  newbi = " << newbi << endl;
      for (int gi=0; gi<GVecsArray.size(); gi++)
	coefs (numfold*gi+offset) = Orbitals[spin](ki, band)->GetCoefs()(gi);
//       newOrbitals(newki,band*newnumk+newki)->Setk (newk);
//       newOrbitals(newki,band*newnumk+newki)->SetEigVal (Orbitals[spin](ki,band)->GetEigVal());
      newOrbitals(newki,newbi)->SetTwist (newTwist);
      newOrbitals(newki,newbi)->SetEigVal ((double)numfold*Orbitals[spin](ki,band)->GetEigVal());
      delete Orbitals[spin](ki, band);
    }
  }
  GVecsArray.resize(newGVecs.size());
  GVecsArray = newGVecs;
  Orbitals[spin].resize(newOrbitals.shape());
  Orbitals[spin] = newOrbitals;
}


void
PlaneWaveSystem::Tile(int spin, Int3 tileFactor) 
{
  // Just rename to save typing
  Int3 &tf = tileFactor;

  if (Comm.MyProc() != 0)
    return;
  ///////////////////////////////////////////
  // First, make all twist vectors postive //
  ///////////////////////////////////////////
  for (int ki=0; ki<Orbitals[spin].extent(0); ki++) {
    Vec3 twist = Orbitals[spin](ki,0)->GetPrimTwist();
    for (int i=0; i<3; i++)
      twist[i] -= floor(twist[i]);
    Orbitals[spin](ki,0)->SetTwist(twist);
  }

  /////////////////////////////////////////////
  // Check to see if k-points form a lattice //
  /////////////////////////////////////////////
  Vec3 minTwist(1.0e12, 1.0e12, 1.0e12), maxTwist(-1.0e12, -1.0e12, -1.0e12);
  Vec3 mink    (1.0e12, 1.0e12, 1.0e12);
  for (int ki=0; ki<Orbitals[spin].extent(0); ki++) {
    Vec3 twist = Orbitals[spin](ki,0)->GetPrimTwist();
    cerr << "twist = " << twist << endl;
    Vec3 k     = Orbitals[spin](ki,0)->GetPrimk();
    minTwist[0] = min(minTwist[0], twist[0]); 
    maxTwist[0] = max(maxTwist[0], twist[0]);
    minTwist[1] = min(minTwist[1], twist[1]); 
    maxTwist[1] = max(maxTwist[1], twist[1]);
    minTwist[2] = min(minTwist[2], twist[2]); 
    maxTwist[2] = max(maxTwist[2], twist[2]);
    mink[0]     = min(k[0],mink[0]);
    mink[1]     = min(k[1],mink[1]);
    mink[2]     = min(k[2],mink[2]);
  }
  // The difference between maxTwist and minTwist should be of the
  // form (n-1)/n.  Therefore, we determine n by
  Vec3 nf;
  nf[0] = -1.0/((maxTwist[0]-minTwist[0]) -1.0);
  nf[1] = -1.0/((maxTwist[1]-minTwist[1]) -1.0);
  nf[2] = -1.0/((maxTwist[2]-minTwist[2]) -1.0);

  // Make sure they are close to integers
  assert (fabs(nf[0] - round(nf[0]))<1.0e-6);
  assert (fabs(nf[1] - round(nf[1]))<1.0e-6);
  assert (fabs(nf[2] - round(nf[2]))<1.0e-6);
  

  Int3 n ((int)round(nf[0]), (int)round(nf[1]), (int)round(nf[2])), newkmesh;

  // Now, make sure we have all the k-points in the lattice
  Vec3 twist;
  for (int ix=0; ix<n[0]; ix++) 
    for (int iy=0; iy<n[1]; iy++)
      for (int iz=0; iz<n[2]; iz++){
	twist[0] = 
	  minTwist[0] + (double)ix/(double)(n[0]-1)*(maxTwist[0]-minTwist[0]);
	twist[1] = 
	  minTwist[1] + (double)iy/(double)(n[1]-1)*(maxTwist[1]-minTwist[1]);
	twist[2] = 
	  minTwist[2] + (double)iz/(double)(n[2]-1)*(maxTwist[2]-minTwist[2]);
	bool twistFound = false;
	for (int ik=0; ik<Orbitals[spin].extent(0); ik++) {
	  Vec3 diff = Orbitals[spin](ik,0)->GetPrimTwist()-twist;
	  if (dot(diff,diff)<1.0e-8) {
	    twistFound = true;
	    for (int band=0; band<Orbitals[spin].extent(1); band++)
	      Orbitals[spin](ik, band)->SetIndex (Int3(ix,iy,iz));
	  }
	}
	if (!twistFound) {
	  fprintf (stderr, "Missing twist vector (%8.4f, %8.4f, %8.4f) "
		   "in tiling process.\n", twist[0], twist[1], twist[2]);
	  abort();
	}
      }
  
  newkmesh[0] = n[0]/tf[0];
  newkmesh[1] = n[1]/tf[1];
  newkmesh[2] = n[2]/tf[2];
  
  // Create new orbitals
  int numk     = Orbitals[spin].extent(0);
  int numbands = Orbitals[spin].extent(1);
  int numtiles = tf[0]*tf[1]*tf[2];
  int newnumk  = numk/numtiles;
  int newnumbands = numbands*numtiles;

  // Create and initialize new orbital set
  Array<OrbitalClass*,2> newOrbitals(newnumk, newnumbands);
  for (int ki=0; ki<newnumk; ki++)
    for (int bi=0; bi<newnumbands; bi++) 
      newOrbitals(ki, bi) = 
	new OrbitalClass (GVecs, PrimFFT, SuperFFT, PrimkPoints(ki));
  
  for (int ki=0; ki<numk; ki++)
    for (int bi=0; bi<numbands; bi++) {
      Int3 index = Orbitals[spin](ki,bi)->GetIndex();
      Int3 newbIndex (index[0]/tf[0], index[1]/tf[1], index[2]/tf[2]);
      Int3 newkIndex (index[0]%tf[0], index[1]%tf[1], index[2]%tf[2]);

      int newki = newkIndex[2] + newkmesh[2]*(newkIndex[1] + newkmesh[1]*newkIndex[0]);
      int newbi = bi*numtiles +
	          newbIndex[2] +  tf[2]*(newbIndex[1] + tf[1]*newbIndex[0]);
      Vec3 newk;
      newOrbitals(ki,bi)->GetCoefs() = Orbitals[spin](ki, bi)->GetCoefs();
      newOrbitals(newki,newbi)->Setk (newk);
      newOrbitals(newki,newbi)->SetEigVal (Orbitals[spin](ki,bi)->GetEigVal());
    }
    


}


inline bool
equal (Vec3 a, Vec3 b)
{
  return ((fabs(a[0]-b[0])<1.0e-10) &&
	  (fabs(a[1]-b[1])<1.0e-10) &&
	  (fabs(a[2]-b[2])<1.0e-10));
}


void
PlaneWaveSystem::Tile(int spin, Mat3 tileMatrix) 
{
  TileMatrix = tileMatrix;
  // Alias for brevity
  Mat3 &S = TileMatrix;
  Mat3 Aprim, Asuper, Bprim, Bsuper;
  Aprim  = Lattice;
  Asuper = S * Lattice;
  Bprim  = Transpose(Inverse(Aprim));
  Bsuper = Transpose(Inverse(Asuper));

  int numk = Orbitals[spin].extent(0);
  Array<Vec3,1> primTwist(numk), superTwist(numk);
  for (int ik=0; ik<numk; ik++) {
    primTwist(ik) = Orbitals[spin](ik,0)->GetTwist();
    superTwist(ik) = S*primTwist(ik);
  }

  // Each super twist can now be decomposed into an integer and
  // fractional part.  The integer part reflects a reciprocal lattice
  // vector of the superlattice and the fractional part is a
  // k-vector.  

  // Keep track of distinct fractional parts.  These will be the twist
  // vectors of the unfolded superlattice.
  vector<Vec3> twists;
  // Loop through all superTwists
  for (int ik=0; ik<numk; ik++) {
    Vec3 intPart, fracPart;
    for (int i=0; i<3; i++) {
      intPart[i] = floor(superTwist(ik)[i]);
      fracPart[i] = superTwist(ik)[i] - intPart[i];
    }
    bool found = false;
    for (int j=0; j<twists.size(); j++)
      found = found || equal(twists[j], fracPart);
    if (!found)
      twists.push_back(fracPart);
  }

  Array<Vec3,2> twistGroups;


}


bool
PlaneWaveSystem::Write (string fname)
{
  if (Comm.MyProc() != 0)
    return true;
  if (Spline && !FFTisSetup) {
    FFT.Setup();
    FFTisSetup = true;
  }

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
      lattice(i,j) = Lattice(i,j);
  out.WriteVar("lattice", lattice);
  Array<double,2> reciprocal_lattice(3,3);
  for (int i=0; i<3; i++)
    for (int j=0; j<3; j++)
      reciprocal_lattice(i,j) = RecipLattice(i,j);
  out.WriteVar("reciprocal_lattice", reciprocal_lattice);
  int complex_coefficients = 1;
  out.WriteVar ("complex_coefficients", complex_coefficients);
  out.WriteVar ("num_bands",      Orbitals[0].extent(1)+Orbitals[1].extent(1));
  out.WriteVar ("num_up_bands",   Orbitals[0].extent(1));
  out.WriteVar ("num_down_bands", Orbitals[1].extent(1));
  out.WriteVar ("num_twists",     Orbitals[0].extent(0));
  int num_spins = SpinPolarized ? 2 : 1;
  out.WriteVar ("num_spins", num_spins);
  out.WriteVar ("maximum_ecut", ECut);
  out.WriteVar ("num_electrons", NumElectrons);
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
  Array<double,2> pos(IonPos.size(),3);
  for (int i=0; i<IonPos.size(); i++) {
    
    pos(i,0) = IonPos(i)[0];
    pos(i,1) = IonPos(i)[1];
    pos(i,2) = IonPos(i)[2];
  }
  out.WriteVar ("pos", pos);
  out.WriteVar ("atom_types", AtomTypes);
  out.CloseSection(); // "ions"
  

  out.NewSection ("eigenstates");
  for (int ik=0; ik<Orbitals[0].extent(0); ik++) {
    out.NewSection("twist");
    Vec3 k = Orbitals[0](ik,0)->GetPrimk();
    Array<double,1> twist(3);
    twist(0) = Lattice(0,0)*k[0] + Lattice(0,1)*k[1] + Lattice(0,2)*k[2];
    twist(1) = Lattice(1,0)*k[0] + Lattice(1,1)*k[1] + Lattice(1,2)*k[2];
    twist(2) = Lattice(2,0)*k[0] + Lattice(2,1)*k[1] + Lattice(2,2)*k[2];
    twist = (0.5/M_PI) * twist;
    out.WriteVar ("twist_angle", twist);
    for (int spin=0; spin<2; spin++) {
      for (int band=0; band<Orbitals[spin].extent(1); band++) {
	out.NewSection("band");      
	if (Spline) {
	  if (CheckKE)
	    Orbitals[spin](ik, band)->CheckKineticEnergy();
	  Orbitals[spin](ik, band)->WriteSpline (out, Real, ShiftOrbitals, Truncate);
	}
	else
	  Orbitals[spin](ik, band)->Write (out);
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
PlaneWaveSystem::CheckkMesh(int spin)
{
  bool meshOK = true;

  ///////////////////////////////////////////
  // First, make all twist vectors postive //
  ///////////////////////////////////////////
//   for (int ki=0; ki<Orbitals[spin].extent(0); ki++) {
//     Vec3 twist = Orbitals[spin](ki,0)->GetTwist();
//     for (int i=0; i<3; i++)
//       twist[i] -= floor(twist[i]);
//     Orbitals[spin](ki,0)->SetTwist(twist);
//   }
  /////////////////////////////////////////////
  // Check to see if k-points form a lattice //
  /////////////////////////////////////////////
  Vec3 minTwist(1.0e12, 1.0e12, 1.0e12), maxTwist(-1.0e12, -1.0e12, -1.0e12);
  Vec3 mink    (1.0e12, 1.0e12, 1.0e12);
  for (int ki=0; ki<Orbitals[spin].extent(0); ki++) {
    Vec3 twist = Orbitals[spin](ki,0)->GetPrimTwist();
    cerr << "twist = " << twist << endl;
    Vec3 k     = Orbitals[spin](ki,0)->GetPrimk();
    minTwist[0] = min(minTwist[0], twist[0]); 
    maxTwist[0] = max(maxTwist[0], twist[0]);
    minTwist[1] = min(minTwist[1], twist[1]); 
    maxTwist[1] = max(maxTwist[1], twist[1]);
    minTwist[2] = min(minTwist[2], twist[2]); 
    maxTwist[2] = max(maxTwist[2], twist[2]);
    mink[0]     = min(k[0],mink[0]);
    mink[1]     = min(k[1],mink[1]);
    mink[2]     = min(k[2],mink[2]);
  }
  // The difference between maxTwist and minTwist should be of the
  // form (n-1)/n.  Therefore, we determine n by
  Vec3 nf;
  nf[0] = -1.0/((maxTwist[0]-minTwist[0]) -1.0);
  nf[1] = -1.0/((maxTwist[1]-minTwist[1]) -1.0);
  nf[2] = -1.0/((maxTwist[2]-minTwist[2]) -1.0);
  
  // Make sure they are close to integers
  meshOK = meshOK && (fabs(nf[0] - round(nf[0]))<1.0e-6);
  meshOK = meshOK && (fabs(nf[1] - round(nf[1]))<1.0e-6);
  meshOK = meshOK && (fabs(nf[2] - round(nf[2]))<1.0e-6);

  Int3 n((int)round(nf[0]), (int)round(nf[1]), (int)round(nf[2]));
  kPointMesh = n;
  
  
  // Now, make sure we have all the k-points in the lattice
  Vec3 twist;
  for (int ix=0; ix<n[0]; ix++) 
    for (int iy=0; iy<n[1]; iy++)
      for (int iz=0; iz<n[2]; iz++){
	twist[0] = 
	  minTwist[0] + (double)ix/(double)(n[0]-1)*(maxTwist[0]-minTwist[0]);
	twist[1] = 
	  minTwist[1] + (double)iy/(double)(n[1]-1)*(maxTwist[1]-minTwist[1]);
	twist[2] = 
	  minTwist[2] + (double)iz/(double)(n[2]-1)*(maxTwist[2]-minTwist[2]);
	bool twistFound = false;
	for (int ik=0; ik<Orbitals[spin].extent(0); ik++) {
	  Vec3 diff = Orbitals[spin](ik,0)->GetPrimTwist()-twist;
	  if (dot(diff,diff)<1.0e-8) {
	    twistFound = true;
	    for (int band=0; band<Orbitals[spin].extent(1); band++)
	      Orbitals[spin](ik, band)->SetIndex (Int3(ix,iy,iz));
	  }
	}
	if (!twistFound) {
	  fprintf (stderr, "Missing twist vector (%8.4f, %8.4f, %8.4f) "
		   "in CheckkPointMesh.\n", twist[0], twist[1], twist[2]);
	  return false;
	}
      }

  return meshOK;
}
