#include "OrbitalSetClass.h"
#include <Common/IO/IO.h>


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
OrbitalSetClass::Read_LAPW(string fname)
{
  IO::IOSectionClass in;
  assert (in.OpenFile (fname));

  // Read system information
  assert (in.OpenSection("parameters"));
  Array<double,2> atom_pos, lattice;
  Array<int,1> atom_types;
  assert (in.ReadVar("atom_pos", atom_pos));
  assert (in.ReadVar("atom_types", atom_types));
  PrimCell.IonPos.resize(atom_pos.extent(0));
  PrimCell.AtomTypes.resize(atom_types.size());
  APW.NumAtoms = atom_pos.extent(0);
  for (int atom=0; atom<APW.NumAtoms; atom++) {
    PrimCell.IonPos(atom) = Vec3 (atom_pos(atom,0),
				  atom_pos(atom,1),
				  atom_pos(atom,2));
    PrimCell.AtomTypes(atom) = atom_types(atom);
  }
			       
  assert (in.ReadVar("atom_types", atom_types));

  ////////////////////
  // Setup lattices //
  ////////////////////
  assert (in.ReadVar("lattice", lattice));
  Mat3 Aprim, Asuper;
  for (int i=0; i<3; i++)
    for (int j=0; j<3; j++)
      Aprim(i,j) = lattice(i,j);
  Asuper = TileMatrix*Aprim;
  PrimCell.SetLattice(Aprim);
  SuperCell.SetLattice(Asuper);
  fprintf (stderr, "Superlattice vectors:\n");
  for (int i=0; i<3; i++)
    fprintf (stderr, "  [ %9.6f %9.6f %9.6f ]\n", 
	     Asuper(i,0), Asuper(i,1), Asuper(i,2));

  double norm = 1.0;// /sqrt(fabs(det(Asuper)));
  
  in.CloseSection(); // "parameters"

  //////////////////////////////////////////////////
  // Read G-vectors for all twists and create map //
  //////////////////////////////////////////////////
  int numk, numBands=0;
  std::map<Int3,int,Int3Less> gMap;
  assert (in.OpenSection("eigenstates"));
  // Read local orbital mappings
  assert (in.ReadVar("local_atom", LocalOrbitals.AtomMap));
  assert (in.ReadVar("local_ilo",  LocalOrbitals.iloMap));
  assert (in.ReadVar("local_l",    LocalOrbitals.lMap));
  assert (in.ReadVar("local_m",    LocalOrbitals.mMap));

  int numSpins = in.CountSections("spin");

  for (int spin=0; spin<numSpins; spin++) {
    assert (in.OpenSection("spin",spin));
    numk = in.CountSections("twist");
    for (int ik=0; ik<numk; ik++) {
      assert (in.OpenSection("twist", ik));
      numBands = in.CountSections("band");
      Array<double,2> gvecs;
      assert (in.ReadVar("gvecs_reduced", gvecs));
      for (int ig=0; ig<gvecs.extent(0); ig++) {
	Int3 gint((int)round(gvecs(ig,0)),
		  (int)round(gvecs(ig,1)),
		  (int)round(gvecs(ig,2)));
	map<Int3,int,Int3Less>::iterator iter = gMap.find(gint);
	int oldSize = gMap.size();
	if (iter == gMap.end()) {
	  int size = gMap.size();
	  gMap[gint] = size;
	  int gi = gMap.find(gint)->second;
	  assert (gMap.size() == gMap.find(gint)->second+1);
	}

      }
      in.CloseSection(); // "twist"
    }
    in.CloseSection(); // "spin"
  }
  in.CloseSection(); // "eigenstates"

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

  ///////////////////////////////////////////
  // Now, read the plane-wave coefficients //
  ///////////////////////////////////////////
  // Create a local copy of the coeficients for building
  // the muffin-tin orbitals.  
  // Indices are (spin,ik,iband)(ig)
  APW.Phi.resize(numSpins, numk, numBands);
  int num_local_coefs = LocalOrbitals.AtomMap.size();
  LocalOrbitals.Phi.resize(numSpins, numk, numBands, num_local_coefs);

  assert (in.OpenSection("eigenstates"));
  for (int spin=0; spin<in.CountSections("spin"); spin++) {
    // Create the orbital objects
    PrimOrbitals[spin].resize(numk, numBands);
    for (int ik=0; ik<numk; ik++) {
      for (int band=0; band<numBands; band++) {
	PrimOrbitals[0](ik,band) = new OrbitalClass;
	PrimOrbitals[0](ik,band)->SetCell (PrimCell);
      }
    }

    assert (in.OpenSection("spin",spin));
    for (int ik=0; ik<numk; ik++) {
      assert (in.OpenSection("twist", ik));
      Array<double,1> twist_angle;
      assert (in.ReadVar("twist_angle", twist_angle));
      Vec3 twist(twist_angle(0), twist_angle(1), twist_angle(2));
      twist = -1.0*twist;
      for (int band=0; band<numBands; band++) {
	assert (in.OpenSection("band", band));
	Array<double,2> eigvec;
	assert (in.ReadVar("phi_apw", eigvec));
	Array<double,2> gvecs;
	assert (in.ReadVar("gvecs_reduced", gvecs));
	zVec &coefs(*new zVec(GVecsArray.size()));
	coefs = complex<double>();
	APW.Phi(spin,ik,band).resize(eigvec.extent(0));
	for (int ig=0; ig<gvecs.extent(0); ig++) {
	  complex<double> z (eigvec(ig,0), eigvec(ig,1));
	  APW.Phi(spin,ik,band)(ig) = z;
	  Int3 gint((int)round(gvecs(ig,0)),
		    (int)round(gvecs(ig,1)),
		    (int)round(gvecs(ig,2)));
	  int index = gMap[gint];
	  coefs(index) = norm*z;
	}
	assert (in.ReadVar("phi_local", eigvec));
	assert (eigvec.extent(0) == num_local_coefs);
	for (int i=0; i<num_local_coefs; i++) 
	  LocalOrbitals.Phi(spin,ik,band,i) = 
	    complex<double>(eigvec(i,0), eigvec(i,1));

	double eigval;
	assert (in.ReadVar("eigenvalue", eigval));
	PrimOrbitals[spin](ik,band)->SetCoefs(coefs);
	PrimOrbitals[spin](ik,band)->SetTwist(twist);
	PrimOrbitals[spin](ik,band)->SetEigVal (eigval);
	
	in.CloseSection(); // "band"
      }
      in.CloseSection(); // "twist"
    }
    in.CloseSection(); // "spin"
  }
  in.CloseSection(); // "eigenstates"

  /////////////////////////////////////////////////////
  // Read radial functions and matching coefficients //
  /////////////////////////////////////////////////////
  assert (in.OpenSection("radialfuncs"));
  assert (in.ReadVar("lmax", APW.lMax));
  int numSpecies = in.CountSections("species");
  APW.u.resize(APW.NumAtoms, APW.lMax+1);
  APW.u_data.resize(APW.NumAtoms, APW.lMax+1);
  APW.du_dr_data.resize(APW.NumAtoms, APW.lMax+1);
  APW.r_data.resize(APW.NumAtoms);
  APW.MatchCoefs.resize(numSpins, APW.NumAtoms, numk);
  LocalOrbitals.v.resize(APW.NumAtoms);
  LocalOrbitals.v_data.resize(APW.NumAtoms);
  int atom = 0;
  for (int species=0; species<numSpecies; species++) {
    assert (in.OpenSection("species",species));
    Array<double,1> r;
    assert (in.ReadVar("r", r));
    vector<double> rvec(r.size());
    for (int ir=0; ir<r.size(); ir++)
      rvec[ir] = r(ir);
    SimpleGrid rGrid;
    rGrid.Init (rvec);
    int natom = in.CountSections("atom");
    for (int ia=0; ia<natom; ia++) {
      APW.r_data(atom).resize(r.size());
      APW.r_data(atom) = r;
      assert (in.OpenSection("atom", ia));
      for (int l=0; l<=APW.lMax; l++) {
	assert (in.OpenSection("APW_l", l));
	// First, read radial splines
	Array<double,2>& u_of_r(APW.u_data(atom,l));
	Array<double,1>& du_dr(APW.du_dr_data(atom,l));
	assert (in.ReadVar("u_of_r", u_of_r));
	assert (in.ReadVar("du_dr_final", du_dr));
	APW.u(atom,l).resize(u_of_r.extent(0));
	vector<double> uvals(u_of_r.extent(1));
	int numOrder = u_of_r.extent(0);
	for (int iord=0; iord<numOrder; iord++) {
	  for (int ir=0; ir<uvals.size(); ir++)
	    uvals[ir] = u_of_r(iord,ir);
	  APW.u(atom,l)[iord].Init (rGrid, uvals);
	}
	// Next, read matching coefficients
	for (int spin=0; spin<numSpins; spin++) {
	  assert (in.OpenSection("spin",spin));
	  for (int ik=0; ik<numk; ik++) {
	    int numG = APW.Phi(spin,ik,0).size();
	    // Resize if it hasn't been done yet
	    if (APW.MatchCoefs(spin,atom,ik).getNumG() == 0) {
// 	      fprintf (stderr, "Resizing (%d,%d,%d) to (%d,%d,%d)\n",
// 		       spin,atom,ik,APW.lMax, numOrder, numG);
	      APW.MatchCoefs(spin,atom,ik).resize(APW.lMax, numOrder, numG);
	    }
	    assert (in.OpenSection("twist", ik));
	    Array<double,3> match_imag, match_real;
	    assert (in.ReadVar("match_real", match_real));
	    assert (in.ReadVar("match_imag", match_imag));
	    assert (match_real.extent(0)==(2*l+1));
	    assert (match_real.extent(1)==numOrder);
	    assert (match_real.extent(2)==numG);
	    for (int m=-l; m<=l; m++)
	      for (int iord=0; iord<numOrder; iord++)
		for (int iG=0; iG<numG; iG++)
		  APW.MatchCoefs(spin,atom,ik)(l,m,iord,iG) =
		    complex<double>(match_real(l+m,iord,iG),
				    match_imag(l+m,iord,iG));
	    in.CloseSection(); // "twist"
	  }
	  in.CloseSection(); // "spin"
	}
	in.CloseSection(); // "APW_l"
      }
      // Now, read local orbital information
      int num_local = in.CountSections("local");
      LocalOrbitals.v(atom).resize(num_local);
      LocalOrbitals.v_data(atom).resize(num_local,r.size());
      for (int ilo=0; ilo<num_local; ilo++) {
	assert (in.OpenSection("local", ilo));
	Array<double,1> v_of_r;
	vector<double> vvec;
	assert (in.ReadVar("v_of_r", v_of_r));
	LocalOrbitals.v_data(atom)(ilo,Range::all()) = v_of_r;
	vvec.resize(v_of_r.size());
	for (int ir=0; ir<v_of_r.size(); ir++)
	  vvec[ir] = v_of_r(ir);
	LocalOrbitals.v(atom)[ilo].Init (rGrid, vvec);
	in.CloseSection(); // "local"
      }

      in.CloseSection(); // "atom"
      atom++;
    }
    in.CloseSection(); // "species"
  }

  in.CloseSection(); // "radialfuncs"

  //////////////////////////
  // Read the core states //
  //////////////////////////
  assert (in.OpenSection("corestates"));
  numSpecies = in.CountSections ("species");
  atom = 0;
  for (int species=0; species<numSpecies; species++) {
    assert (in.OpenSection("species",species));
    int natom = in.CountSections("atom");
    Array<double,1> r;
    assert (in.ReadVar("r", r));
    for (int ia=0; ia<natom; ia++,atom++) {
      assert (in.OpenSection("atom",ia));
      int nstates = in.CountSections("state");
      for (int state=0; state<nstates; state++) {
	assert (in.OpenSection("state", state));
	CoreStateClass core;
	core.r.resize(r.size()); 
	core.r = r;
	core.atom = atom;
	assert (in.ReadVar("g0", core.g0));
	assert (in.ReadVar("f0", core.f0));
	assert (in.ReadVar("n" , core.n));
	assert (in.ReadVar("l" , core.l));
	assert (in.ReadVar("k" , core.k));
	assert (in.ReadVar("eigenvalue", core.eigenvalue));
	CoreStates.push_back(core);
	in.CloseSection(); // "state"
      }
      in.CloseSection(); // "atom"
    }
    in.CloseSection(); // "species"
  }

  in.CloseSection(); // "corestates"
  in.CloseFile();


  SetupFFTIndices();
  CreateSuperOrbitals();
  TileIonPos();
  return true;
}
