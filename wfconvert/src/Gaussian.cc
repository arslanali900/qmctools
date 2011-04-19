#include "Gaussian.h"
#include "ParserClass.h"
#include <Common/IO/IO.h>
#include <map>



bool
Shell::ReadGAMESS (ParserClass &log)
{
  std::map<std::string,int> channelMap;
  channelMap["S"] = 0;  channelMap["P"] = 1;
  channelMap["D"] = 2;  channelMap["F"] = 3;
  channelMap["G"] = 4;  channelMap["H"] = 5;

  blitz::TinyVector<double,2> cont;
  // Read first primitive
  int shell, prim;
  string channel;
  assert(log.ReadInt (shell));
  assert(log.ReadWord(channel));
  l = channelMap[channel];
  assert(log.ReadInt (prim));
  assert(log.ReadDouble (cont[1]));
  assert(log.ReadDouble (cont[0]));
  push_back(cont);

  fprintf (stderr, "%5d  %2s l=%d %3d   %12.8f   %16.12f\n",
	   shell, channel.c_str(), l, prim, cont[1], cont[0]);
  
  bool done = false;
  bool lastShell = false;
  while (!done) {
    log.SavePos();
    int sh;
    if (!log.ReadInt(sh)) {
      done = true;
      lastShell = true;
      log.RestorePos();
    }
    else if (sh != shell) {
      fprintf (stderr, "\n");
      done = true;
      log.RestorePos();
    }
    else {
      string ch;
      assert (log.ReadWord(ch));
      assert (ch == channel);
      int p = prim;
      assert (log.ReadInt(prim));
      assert (prim = p+1);
      assert (log.ReadDouble(cont[1]));
      assert (log.ReadDouble(cont[0]));
      fprintf (stderr, "%5d  %2s l=%d %3d   %12.8f   %16.12f\n",
	       shell, channel.c_str(), l, prim, cont[1], cont[0]);
      push_back(cont);
    }

  }
  //  fprintf (stderr, "norm = %1.8f\n", nrm());
  return lastShell;
}


GaussianOrbitalSet::GaussianOrbitalSet()
{
  ZToSymbolMap[1]   = "H";  ZToSymbolMap[2]   = "He";
  ZToSymbolMap[3]   = "Li"; ZToSymbolMap[4]   = "Be";
  ZToSymbolMap[5]   = "B";  ZToSymbolMap[6]   = "C";
  ZToSymbolMap[7]   = "N";  ZToSymbolMap[8]   = "O";
  ZToSymbolMap[9]   = "F";  ZToSymbolMap[10]  = "Ne";
  ZToSymbolMap[11]  = "Na"; ZToSymbolMap[12]  = "Mg";
  ZToSymbolMap[13]  = "Al"; ZToSymbolMap[14]  = "Si";
  ZToSymbolMap[15]  = "P";  ZToSymbolMap[16]  = "S";
  ZToSymbolMap[17]  = "Cl"; ZToSymbolMap[18]  = "Ar";
  ZToSymbolMap[19]  = "K";  ZToSymbolMap[20]  = "Ca";
  ZToSymbolMap[21]  = "Sc"; ZToSymbolMap[22]  = "Ti";
  ZToSymbolMap[23]  = "V";  ZToSymbolMap[24]  = "Cr";
  ZToSymbolMap[25]  = "Mn"; ZToSymbolMap[26]  = "Fe";
  ZToSymbolMap[27]  = "Co"; ZToSymbolMap[28]  = "Ni";
  ZToSymbolMap[29]  = "Cu"; ZToSymbolMap[30]  = "Zn";
  ZToSymbolMap[31]  = "Ga"; ZToSymbolMap[32]  = "Ge";
  ZToSymbolMap[33]  = "As"; ZToSymbolMap[34]  = "Se";
  ZToSymbolMap[35]  = "Br"; ZToSymbolMap[36]  = "Kr";
  ZToSymbolMap[37]  = "Rb"; ZToSymbolMap[38]  = "Sr";
  ZToSymbolMap[39]  = "Y";  ZToSymbolMap[40]  = "Zr";
  ZToSymbolMap[41]  = "Nb"; ZToSymbolMap[42]  = "Mo";
  ZToSymbolMap[43]  = "Tc"; ZToSymbolMap[44]  = "Ru";
  ZToSymbolMap[45]  = "Rh"; ZToSymbolMap[46]  = "Pd";
  ZToSymbolMap[47]  = "Ag"; ZToSymbolMap[48]  = "Cd";
  ZToSymbolMap[49]  = "In"; ZToSymbolMap[50]  = "Sn";
  ZToSymbolMap[51]  = "Sb"; ZToSymbolMap[52]  = "Te";
  ZToSymbolMap[53]  = "I";  ZToSymbolMap[54]  = "Xe";
  ZToSymbolMap[55]  = "Cs"; ZToSymbolMap[56]  = "Ba";
  ZToSymbolMap[57]  = "La"; ZToSymbolMap[58]  = "Ce";
  ZToSymbolMap[59]  = "Pr"; ZToSymbolMap[60]  = "Nd";
  ZToSymbolMap[61]  = "Pm"; ZToSymbolMap[62]  = "Sm";
  ZToSymbolMap[63]  = "Eu"; ZToSymbolMap[64]  = "Gd";
  ZToSymbolMap[65]  = "Tb"; ZToSymbolMap[66]  = "Dy";
  ZToSymbolMap[67]  = "Ho"; ZToSymbolMap[68]  = "Er";
  ZToSymbolMap[69]  = "Tm"; ZToSymbolMap[70]  = "Yb";
  ZToSymbolMap[71]  = "Lu"; ZToSymbolMap[72]  = "Hf";
  ZToSymbolMap[73]  = "Ta"; ZToSymbolMap[74]  = "W";
  ZToSymbolMap[75]  = "Re"; ZToSymbolMap[76]  = "Os";
  ZToSymbolMap[77]  = "Ir"; ZToSymbolMap[78]  = "Pt";
  ZToSymbolMap[79]  = "Au"; ZToSymbolMap[80]  = "Hg";
  ZToSymbolMap[81]  = "Tl"; ZToSymbolMap[82]  = "Pb";
  ZToSymbolMap[83]  = "Bi"; ZToSymbolMap[84]  = "Po";
  ZToSymbolMap[85]  = "At"; ZToSymbolMap[86]  = "Rn";
  ZToSymbolMap[87]  = "Fr"; ZToSymbolMap[88]  = "Ra";
  ZToSymbolMap[89]  = "Ac"; ZToSymbolMap[90]  = "Th";
  ZToSymbolMap[91]  = "Pa"; ZToSymbolMap[92]  = "U";
  ZToSymbolMap[93]  = "Np"; ZToSymbolMap[94]  = "Pu";
  ZToSymbolMap[95]  = "Am"; ZToSymbolMap[96]  = "Cm";
  ZToSymbolMap[97]  = "Bk"; ZToSymbolMap[98]  = "Cf";
  ZToSymbolMap[99]  = "Es"; ZToSymbolMap[100] = "Fm";
  ZToSymbolMap[101] = "Mc"; ZToSymbolMap[102] = "No";
  ZToSymbolMap[103] = "Lw"; 

  ZToMassMap[  1] = 1.00794;
  ZToMassMap[  2] = 4.002602;
  ZToMassMap[  3] = 6.941;
  ZToMassMap[  4] = 9.012182;
  ZToMassMap[  5] = 10.811;
  ZToMassMap[  6] = 12.0107;
  ZToMassMap[  7] = 14.0067;
  ZToMassMap[  8] = 15.9994;
  ZToMassMap[  9] = 18.9984032;
  ZToMassMap[ 10] = 20.1797;
  ZToMassMap[ 11] = 22.98976928;
  ZToMassMap[ 12] = 24.3050;
  ZToMassMap[ 13] = 26.9815386;
  ZToMassMap[ 14] = 28.0855;
  ZToMassMap[ 15] = 30.973762;
  ZToMassMap[ 16] = 32.065;
  ZToMassMap[ 17] = 35.453;
  ZToMassMap[ 18] = 39.948;
  ZToMassMap[ 19] = 39.0983;
  ZToMassMap[ 20] = 40.078;
  ZToMassMap[ 21] = 44.955912;
  ZToMassMap[ 22] = 47.867;
  ZToMassMap[ 23] = 50.9415;
  ZToMassMap[ 24] = 51.9961;
  ZToMassMap[ 25] = 54.938049;
  ZToMassMap[ 26] = 55.845;
  ZToMassMap[ 27] = 58.933200;
  ZToMassMap[ 28] = 58.6934;
  ZToMassMap[ 29] = 63.546;
  ZToMassMap[ 30] = 65.39;
  ZToMassMap[ 31] = 69.723;
  ZToMassMap[ 32] = 72.61;
  ZToMassMap[ 33] = 74.92160;
  ZToMassMap[ 34] = 78.96 ;
  ZToMassMap[ 35] = 79.904;
  ZToMassMap[ 36] = 83.80;
  ZToMassMap[ 37] = 85.4678;
  ZToMassMap[ 38] = 87.62;
  ZToMassMap[ 39] = 88.90585;
  ZToMassMap[ 40] = 91.224;
  ZToMassMap[ 41] = 92.90638;
  ZToMassMap[ 42] = 95.94;
  ZToMassMap[ 43] = 98;
  ZToMassMap[ 44] = 101.07;
  ZToMassMap[ 45] = 102.90550;
  ZToMassMap[ 46] = 106.42;
  ZToMassMap[ 47] = 107.8682;
  ZToMassMap[ 48] = 112.411;
  ZToMassMap[ 49] = 114.818;
  ZToMassMap[ 50] = 118.710;
  ZToMassMap[ 51] = 121.760;
  ZToMassMap[ 52] = 127.60;
  ZToMassMap[ 53] = 126.90447;
  ZToMassMap[ 54] = 131.29;
  ZToMassMap[ 55] = 132.90545;
  ZToMassMap[ 56] = 137.327;
  ZToMassMap[ 57] = 138.9055;
  ZToMassMap[ 58] = 140.116;
  ZToMassMap[ 59] = 140.90765;
  ZToMassMap[ 60] = 144.24;
  ZToMassMap[ 61] = 145;
  ZToMassMap[ 62] = 150.36;
  ZToMassMap[ 63] = 151.964;
  ZToMassMap[ 64] = 157.25;
  ZToMassMap[ 65] = 158.92534;
  ZToMassMap[ 66] = 162.50;
  ZToMassMap[ 67] = 164.93032;
  ZToMassMap[ 68] = 167.26;
  ZToMassMap[ 69] = 168.93421;
  ZToMassMap[ 70] = 173.04;
  ZToMassMap[ 71] = 174.967;
  ZToMassMap[ 72] = 178.49;
  ZToMassMap[ 73] = 180.9479;
  ZToMassMap[ 74] = 183.84;
  ZToMassMap[ 75] = 186.207;
  ZToMassMap[ 76] = 190.23;
  ZToMassMap[ 77] = 192.217;
  ZToMassMap[ 78] = 195.078;
  ZToMassMap[ 79] = 196.96655;
  ZToMassMap[ 80] = 200.59;
  ZToMassMap[ 81] = 204.3833;
  ZToMassMap[ 82] = 207.2;
  ZToMassMap[ 83] = 208.98038;
  ZToMassMap[ 84] = 209;
  ZToMassMap[ 85] = 210;
  ZToMassMap[ 86] = 222;
  ZToMassMap[ 87] = 223;
  ZToMassMap[ 88] = 226;
  ZToMassMap[ 89] = 227;
  ZToMassMap[ 90] = 232.0381;
  ZToMassMap[ 91] = 231.03588;
  ZToMassMap[ 92] = 238.0289;
  ZToMassMap[ 93] = 237;
  ZToMassMap[ 94] = 244;
  ZToMassMap[ 95] = 243;
  ZToMassMap[ 96] = 247;
  ZToMassMap[ 97] = 247;
  ZToMassMap[ 98] = 251;
  ZToMassMap[ 99] = 252;
  ZToMassMap[100] = 257;
  ZToMassMap[101] = 258;
  ZToMassMap[102] = 259;
  ZToMassMap[103] = 262;
  ZToMassMap[104] = 261;
  ZToMassMap[105] = 262;
  ZToMassMap[106] = 263;
  ZToMassMap[107] = 262;
  ZToMassMap[108] = 265;
  ZToMassMap[109] = 266;
  ZToMassMap[110] = 269;
  ZToMassMap[111] = 272;
  ZToMassMap[112] = 277;
  ZToMassMap[113] = 284;
  ZToMassMap[114] = 289;
  ZToMassMap[115] = 288;
  ZToMassMap[116] = 293;
  ZToMassMap[117] = 293; // UNKNOWN
  ZToMassMap[118] = 294;

  for (int i=1; i<103; i++)
    SymbolToZMap[ZToSymbolMap[i]] = i;


}


bool
GaussianOrbitalSet::ReadGAMESS (string logname)
{
  MemParserClass log;
  assert(log.OpenFile (logname));

  /////////////////////////////
  // Read atomic coordinates //
  /////////////////////////////
  assert (log.FindToken ("COORDINATES (BOHR)"));
  assert (log.FindToken ("Z"));

  string token;
  log.ReadWord(token);
  while (token != "INTERNUCLEAR") {
    Atom atom;
    atom.name = token;
    assert (log.ReadDouble (atom.charge));
    assert (log.ReadDouble (atom.pos[0]));
    assert (log.ReadDouble (atom.pos[1]));
    assert (log.ReadDouble (atom.pos[2]));
    atom.Z = SymbolToZMap[token];
    Atoms.push_back(atom);
    log.ReadWord(token);
  }
  fprintf (stderr, "Atom     Position\n");
  for (int iat=0; iat<Atoms.size(); iat++) {
    Atom &atom = Atoms[iat];
    fprintf (stderr, "%3s    %10.5f %10.5f %10.5f\n", atom.name.c_str(), 
	     atom.pos[0], atom.pos[1], atom.pos[2]);
  }

  assert(log.FindToken("CONTRACTION COEFFICIENT(S)"));

  for (int iat=0; iat<Atoms.size(); iat++) {
    fprintf (stderr, "\nAtom %s\n------\n", Atoms[iat].name.c_str());
    Site site;
    site.set_pos (Atoms[iat].pos);
    string name;
    log.ReadWord (name);
    assert (name == Atoms[iat].name);
    Shell shell;
    while (!shell.ReadGAMESS(log)) {
      site.push_back(shell);
      shell.clear();
    }
    // Push back last shell
    site.push_back(shell);
    Basis.push_back(site);
  }

  int nb = Basis.size();
  int num_orbs = Basis.num_orbs();
  fprintf (stderr, "\nBasis has %d elements and %d orbitals.\n", nb, num_orbs);

  // Read number of electrons
  assert (log.FindToken("NUMBER OF ELECTRONS"));
  assert (log.FindToken("="));
  assert (log.ReadInt(NumElecs));
  assert (log.FindToken("SCFTYP="));
  string scftype;
  assert (log.ReadWord(scftype));
  NumSpins = scftype == "RHF" ? 1 : 2;

  if (log.FindToken("NUMBER OF ELECTRONS KEPT IN THE CALCULATION IS ="))
    assert (log.ReadInt(NumElecs));
  if (NumSpins == 2) {
    cerr << "  Reading spin-polarized calculation.\n";
    assert (log.FindToken("(ALPHA)"));
    assert (log.FindToken ("="));
    assert (log.ReadInt(SpinElecs[0]));
    assert (log.FindToken("(BETA )"));
    assert (log.FindToken ("="));
    assert (log.ReadInt(SpinElecs[1]));
    assert (NumElecs == SpinElecs[0] + SpinElecs[1]);
  }
  else {
    SpinElecs[1] = NumElecs/2;
    SpinElecs[0] = NumElecs - SpinElecs[1];
  }

  fprintf (stderr, "  There are %d up electrons and %d down electrons.\n",
	   SpinElecs[0], SpinElecs[1]);

  // Now, read coeficients
  set_num_orbs(num_orbs);
  for (int ispin=0; ispin < NumSpins; ispin++) {
    blitz::Array<double,2> &C = Coefs[ispin];
    assert (log.FindToken("EIGENVECTORS"));
    assert (log.FindToken("------------"));
    
    int num_sets = (num_orbs+4)/5;
    
    int nread=0;
    
    for (int iset=0; iset<num_sets; iset++) {
      int ncols = min (5, num_orbs - 5*iset);
      int iorb;
      for (int icol=0; icol<ncols; icol++) {
	log.ReadInt(iorb);
	iorb--;
	assert (iorb == icol + 5*iset);
      }
      iorb = 5*iset;
      for (int icol=0; icol<ncols; icol++) 
	assert(log.ReadDouble(Eigenvalues[iorb++]));
      
      string s;
      for (int icol=0; icol<ncols; icol++) {
	assert (log.ReadWord(s));
	//assert (s == "A");
      }
      for (int ib=0; ib<nb; ib++) {
	int i;
	assert(log.ReadInt(i));
	assert(i == ib+1);
	assert(log.ReadWord(s));
	assert(log.ReadInt (i));
	assert(log.ReadWord(s));
	iorb = 5*iset;
	for (int icol=0; icol<ncols; icol++) {
	  assert(log.ReadDouble(C(iorb++,ib)));
	  nread++;
	}
      }
    }
    fprintf (stderr, "Read %d coefficients.\n", nread);
    assert (nread == nb*num_orbs);
  }


  return true;
}


void
GaussianOrbitalSet::Write_ESHDF(string fname)
{
  IO::IOSectionClass out;
  assert (out.NewFile(fname));

  TinyVector<double,3> r0(-7.5, -7.5, -7.5);
  TinyVector<double,3> r1( 7.5,  7.5,  7.5);
  TinyVector<int,3> mesh(150,150,150);  
  int num_orb = Coefs[0].extent(0);

  out.WriteVar("format", "ES-HDF");
  out.WriteVar("version", TinyVector<int,3>(2,1,0));
  out.NewSection("creator");
  out.WriteVar("program_name", "wfconv");
  out.WriteVar("version", TinyVector<int,3>(1,0,0));
  out.CloseSection();
  out.NewSection("application");
  out.CloseSection(); // "application"
  out.NewSection("supercell");
  Array<double,2> lattice(3,3);
  for (int i=0; i<3; i++)
    for (int j=0; j<3; j++)
      lattice(i,j) = (i==j) ? (r1[i] - r0[i]) : 0.0;
  out.WriteVar("primitive_vectors", lattice);
  out.CloseSection(); // "supercell"
  out.NewSection("atoms");
  out.WriteVar("number_of_atoms", (int)Atoms.size());

  // Compute number of distinct species
  vector<int> species, zion;
  Array<int,1> atom_species(Atoms.size());
  Array<TinyVector<double,3>,1> atom_pos(Atoms.size());
  for (int iat=0; iat<Atoms.size(); iat++) {
    int Z = Atoms[iat].Z;
    atom_pos(iat) = Atoms[iat].pos - r0;
    bool found = false;
    for (int i=0; i<species.size(); i++) 
      if (Z == species[i]) {
	atom_species(iat) = i;
	found = true;
      }
    if (!found) {
      atom_species(iat) = species.size();
      species.push_back(Z);
      zion.push_back(Atoms[iat].charge);
    }
  }

  out.WriteVar("number_of_species", (int)species.size());
  out.WriteVar("species_ids", atom_species);
  out.SetUnderscores(true);
  for (int isp=0; isp<species.size(); isp++) {
    out.NewSection("species");
    int Z = species[isp];
    out.WriteVar ("atomic_number", Z);
    out.WriteVar ("name", ZToSymbolMap[Z]);
    out.WriteVar ("mass", ZToMassMap[Z]);
    out.WriteVar ("psuedopotential", "unknown");
    if (zion.size() > isp)
      out.WriteVar ("valence_charge", zion[isp]);
    out.CloseSection(); // "species"
  }
  out.SetUnderscores(false);


  out.WriteVar("positions", atom_pos);
//   Array<Vec3,1> reduced_positions(cell.IonPos.size());
//   for (int i=0; i<cell.IonPos.size(); i++)
//     reduced_positions(i) = cell.Lattice.r2u(cell.IonPos(i));
//   out.WriteVar("reduced_positions", reduced_positions);

  out.CloseSection(); // "atoms"

  out.NewSection("electrons");
  out.WriteVar("number_of_electrons", SpinElecs);
  out.WriteVar("psi_r_mesh", mesh);
  out.WriteVar ("number_of_kpoints", 1);
  out.WriteVar ("functional", Functional);
  out.WriteVar ("total_energy", TotalEnergy);
  out.WriteVar ("number_of_spins", NumSpins);
  out.WriteVar ("psi_r_is_complex", (int)0);
  out.WriteVar ("number_of_atomic_orbitals",0);

  out.SetUnderscores(true);
  out.NewSection("kpoint");
  out.WriteVar("reduced_k", TinyVector<double,3>(0.0, 0.0, 0.0));
  out.WriteVar("weight", 1.0);
  for (int ispin=0; ispin<NumSpins; ispin++) {
    out.NewSection("spin");
    out.WriteVar("number_of_states", num_orb);
    
    
    Array<double,4> phi(mesh[0], mesh[1], mesh[2], num_orb);
    
    TinyVector<double,3> dr;
    dr[0] = (r1[0] - r0[0])/mesh[0];
    dr[1] = (r1[1] - r0[1])/mesh[1];
    dr[2] = (r1[2] - r0[2])/mesh[2];
    TinyVector<double,3> r;
    vector<double> phir(num_orb);
    //   for (int ix=0; ix<mesh[0]; ix++) {
    //     r[0] = r0[0] + dr[0]*ix;
    //     for (int iy=0; iy<mesh[1]; iy++) {
    //       r[1] = r0[1] + dr[1]*iy;
    //       for (int iz=0; iz<mesh[2]; iz++) {
    // 	r[2] = r0[2] + dr[2]*iz;
    // 	(*this)(ispin, r, phir);
    // 	for (int iorb=0; iorb<num_orb; iorb++)
    // 	  phi(ix,iy,iz,iorb) = phir[iorb];
    //       }
    //     }
    //   }
    
    vector<TinyVector<double,3> > rlist(mesh[1]*mesh[2]);
    Array<double,2> phivals(mesh[1]*mesh[2],num_orb);
    for (int ix=0; ix<mesh[0]; ix++) {
      r[0] = r0[0] + dr[0]*ix;
      int ir=0;
      for (int iy=0; iy<mesh[1]; iy++) {
	r[1] = r0[1] + dr[1]*iy;
	for (int iz=0; iz<mesh[2]; iz++) {
	  r[2] = r0[2] + dr[2]*iz;
	  rlist[ir++] = r;
	}
      }
      (*this)(ispin, rlist, phivals);
      ir = 0;
      for (int iy=0; iy<mesh[1]; iy++)
	for (int iz=0; iz<mesh[2]; iz++) {
	  for (int iorb=0; iorb<num_orb; iorb++)
	    phi(ix,iy,iz,iorb) = phivals(ir,iorb);
	  ir++;
	}
    }
    
    Array<double,3> orb(mesh);
    Array<double,1> eigvals(num_orb);
    for (int iorb=0; iorb<num_orb; iorb++) {
      out.NewSection("state");
      orb = phi(Range::all(),Range::all(),Range::all(),iorb);
      out.WriteVar("psi_r", orb);
      out.CloseSection(); // "state"
      eigvals(iorb) = Eigenvalues[iorb];
    }
    out.WriteVar("eigenvalues", eigvals);
    out.CloseSection(); // "spin"
  }
  out.CloseSection(); // "kpoint"
  out.CloseSection(); // "electrons"
  out.CloseFile();
}
