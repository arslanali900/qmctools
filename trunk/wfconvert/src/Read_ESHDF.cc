#include "OrbitalSetClass.h"
#include <cassert>


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
OrbitalSetClass::Read_ESHDF(string fname)
{
  IO::IOSectionClass in;
  assert (in.OpenFile(fname));
  string format;
  in.ReadVar("format", format);
  if (format != "ES-HDF")
    return false;

  assert (in.OpenSection("supercell"));
  Array<double,2> lattice;
  assert (in.ReadVar("primitive_vectors", lattice));
  in.CloseSection(); // "supercell"
  Mat3 Aprim, Asuper;
  for (int i=0; i<3; i++)
    for (int j=0; j<3; j++)
      Aprim(i,j) = lattice(i,j);
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
  assert (in.OpenSection("atoms"));
  int num_atoms;
  assert (in.ReadVar("number_of_atoms", num_atoms));
  Array<double,2> pos;
  assert (in.ReadVar("positions", pos));
  PrimCell.IonPos.resize(pos.extent(0));
  for (int i=0; i<pos.extent(0); i++)
    for (int j=0; j<3; j++)
      PrimCell.IonPos(i)[j] = pos(i,j);
  PrimCell.AtomTypes.resize(num_atoms);
  PrimCell.Zion.resize(num_atoms);
  Array<int,1> atom_species;
  assert (in.ReadVar("species_ids", atom_species));

  int num_species;
  assert (in.ReadVar("number_of_species", num_species));

  for (int isp=0; isp<num_species; isp++) {
    assert (in.OpenSection("species", isp));
    int Zion, Z;
    assert (in.ReadVar("atomic_number", Z));
    assert (in.ReadVar("valence_charge", Zion));
    for (int iat=0; iat<atom_species.size(); iat++)
      if (atom_species(iat) == isp) {
	PrimCell.AtomTypes(iat) = Z;
	PrimCell.Zion(iat)      = Zion;
      }
    in.CloseSection(); // "species"
  }
  in.CloseSection(); //"atoms"

  assert (in.OpenSection("electrons"));
  int num_spins;
  assert (in.ReadVar ("number_of_spins", num_spins));
  SpinPolarized = (num_spins > 1);

  //////////////////
  // Read Density //
  //////////////////
  PrimDensity.SetLattice  (PrimCell.Lattice.GetDirect());
  SuperDensity.SetLattice (SuperCell.Lattice.GetDirect());
  // First, read density to determine cutoff
  Array<int,2> density_garray;
  Array<TinyVector<int,3>,1> density_gvecs;

  Array<double,2> density_array[2];
  zVec density_g[2];
  if (in.OpenSection("density")) {
    assert (in.ReadVar("gvectors", density_garray));
    density_gvecs.resize(density_garray.extent(0));
    int numG = density_garray.extent(0);
    for (int iG=0; iG<numG; iG++)
      for (int dim=0; dim<3; dim++)
	density_gvecs(iG)[dim] = density_garray(iG, dim);
    
    assert (in.OpenSection("spin", 0));
    assert (in.ReadVar("density_g", density_array[0]));
    density_g[0].resize(numG);
    for (int iG=0; iG<density_array[0].extent(0); iG++)
      density_g[0](iG) = complex<double> (density_array[0](iG,0), 
					  density_array[0](iG,1));
    in.CloseSection(); // "spin"
    if (in.OpenSection("spin", 1)) {
      assert (in.ReadVar("density_g", density_array[1]));
      density_g[1].resize(numG);
      for (int iG=0; iG<density_array[1].extent(0); iG++)
	density_g[1](iG) = complex<double> (density_array[1](iG,0), 
					    density_array[1](iG,1));
      PrimDensity.Set(density_gvecs, density_g[0], density_g[1]);
      in.CloseSection(); // "spin"
    }
    else
      PrimDensity.Set(density_gvecs, density_g[0]);
    in.CloseSection(); // "density"
  }

  // Now, read orbitals
  int num_k, num_states;
  assert (in.ReadVar("number_of_kpoints", num_k));
  assert (in.OpenSection("kpoint",0));
  assert (in.OpenSection("spin", 0));
  assert (in.ReadVar("number_of_states", num_states));
  in.CloseSection(); // "spin"
  in.CloseSection(); // "kpoint"
  PrimOrbitals[0].resize(num_k, num_states);
  if (num_spins == 2)
    PrimOrbitals[1].resize(num_k, num_states);

  ////////////////////
  // Read G-vectors //
  ////////////////////
  // First pass:  Read G-vectors for all different k-points, and then
  // take the union of the G-vector sets
  int bandIndex = 0;
  char dummy[1000];
  std::map<Int3,int,Int3Less> gMap;
  vector<vector<Int3> > gvecs_k (num_k);
  for (int ik=0; ik < num_k; ik++) {
    assert (in.OpenSection("kpoint", ik));
    Array<int,2> gvecs_array;
    assert (in.ReadVar("gvectors", gvecs_array));
    int num_g = gvecs_array.extent(0);
    gvecs_k[ik].resize (num_g);
    for (int ig=0; ig<num_g; ig++) {
      TinyVector<int,3> g(gvecs_array(ig,0), gvecs_array(ig,1),
			  gvecs_array(ig,2));
      gvecs_k[ik][ig] = g;
      map<Int3,int,Int3Less>::iterator iter = gMap.find(g);
      if (iter == gMap.end()) {
	int n = gMap.size();
	gMap[g] = n;
      }
    }
    in.CloseSection(); // "kpoint" 
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



  zVec tmpcoefs;
  Array<double,1> eigenvalues, occupations;
  for (int ik=0; ik<num_k; ik++) {
    assert (in.OpenSection("kpoint", ik));
    Vec3 twist;
    assert (in.ReadVar("reduced_k", twist));
    twist = -1.0 * twist;
    for (int ispin=0; ispin<num_spins; ispin++) {
      assert (in.OpenSection("spin", ispin));
      int nst;
      assert (in.ReadVar("number_of_states", nst));
      assert (nst == num_states);
      assert (in.ReadVar("eigenvalues", eigenvalues));
      cerr << "eignevalues = " << eigenvalues << endl;
      if (!in.ReadVar("occupations", occupations))
	occupations.resize(eigenvalues.size(),1.0);
      
      Array<double,2> coefs_array;
      for (int ist=0; ist<nst; ist++) {
	PrimOrbitals[ispin](ik,ist) = new OrbitalClass;
	assert (in.OpenSection("state", ist));
	assert (in.ReadVar ("psi_g", coefs_array));
	zVec &coefs((*new zVec(coefs_array.extent(0))));
	vector<Int3> &gvecs = gvecs_k[ik];
	assert (coefs.size() == gvecs.size());
	for (int i=0; i<gvecs.size(); i++) {
	  int iG = gMap[gvecs[i]];
	  coefs(iG) = complex<double>(coefs_array(i,0),
				      coefs_array(i,1));
	}
	PrimOrbitals[ispin](ik,ist)->SetEigVal    (eigenvalues(ist));
	PrimOrbitals[ispin](ik,ist)->SetOccupancy (occupations(ist));
	PrimOrbitals[ispin](ik,ist)->SetLabels(ispin, ik, ist);
	PrimOrbitals[ispin](ik,ist)->SetTwist(twist);
	PrimOrbitals[ispin](ik,ist)->SetCell(PrimCell);
	PrimOrbitals[ispin](ik,ist)->SetCoefs(coefs);
	in.CloseSection(); // "state"
      }
      in.CloseSection(); // "spin"
    }
    in.CloseSection(); // "kpoint"
  }
  in.CloseSection(); // "electrons"

  SetupFFTIndices();
  CreateSuperOrbitals();
  if (UseMultiRep)
    CreateFineOrbitals();
  TileIonPos();

  return true;
}
