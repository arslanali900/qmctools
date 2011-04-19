#include "LAPWClass.h"
#include <Common/IO/IO.h>

using namespace IO;

void
LAPWClass::Read(string fname)
{
  IOSectionClass in;
  assert(in.OpenFile (fname));
  assert(in.OpenSection("parameters"));
  Array<double,2> avec, atom_pos;
  Mat3 latVecs;
  assert(in.ReadVar("lattice",avec));
  for (int i=0; i<3; i++)
    for (int j=0; j<3; j++)
      latVecs(i,j)=avec(i,j);
  Lattice.SetDirect(latVecs);
  Volume = fabs(det(latVecs));
  PWNorm = 1.0/sqrt(Volume);


  assert(in.ReadVar("atom_pos", atom_pos));
  int num_atoms = atom_pos.extent(0);
  AtomPos.resize(num_atoms);
  AtomRadii.resize(num_atoms);
  for (int ia=0; ia<num_atoms; ia++)
    for (int j=0; j<3; j++)
      AtomPos(ia)[j] = atom_pos(ia,j);
  assert(in.ReadVar("atom_types", AtomTypes));
  APWorder.resize(num_atoms);
  in.CloseSection();// "parameters"
  /////////////////////////// 
  // Read radial functions //
  ///////////////////////////
  assert (in.OpenSection("radialfuncs"));
  int iatom = 0;
  int numSpecies = in.CountSections("species");
  int maxAPWorder = 0;
  int maxChannels = 0;
  int maxLocal = 0;
  int maxk=0, maxG=0;
  for (int is=0; is<numSpecies; is++) {
    assert(in.OpenSection("species",is));
    int natom = in.CountSections("atom");
    assert (in.OpenSection("atom"));
    maxChannels = max(maxChannels, in.CountSections("APW_l"));
    maxLocal = max(maxLocal, in.CountSections("local"));
    assert(in.OpenSection("APW_l"));
    Array<double,2> uofr;
    in.ReadVar("u_of_r", uofr);
    maxAPWorder = max(maxAPWorder, uofr.extent(0));
    assert (in.OpenSection("spin"));
    int numk = in.CountSections("twist");
    maxk = max(maxk, numk);
    for (int ik=0; ik<numk; ik++) {
      assert (in.OpenSection("twist", ik));
      Array<double,3> match;
      assert (in.ReadVar("match_real", match));
      maxG = max(maxG, match.extent(2));
      in.CloseSection(); // "twist"
    }
    in.CloseSection(); // "spin"
    in.CloseSection(); // "APW_l"
    in.CloseSection(); // "atom"
    in.CloseSection(); // "species"
  }
  lMax = maxChannels - 1;
  // resize u and actually read
  u.resize(num_atoms, maxChannels, maxAPWorder);
  v.resize(num_atoms, maxLocal);
  // Resize the matching coefficients
  int num_lm = maxChannels*maxChannels;
  //  Ylm.resize(num_lm);
  MatchCoefs.resize(num_atoms, num_lm, maxAPWorder, maxk, maxG);
  int atom=0;
  for (int is=0; is<numSpecies; is++) {
    assert (in.OpenSection("species", is));
    Array<double,1> r;
    assert (in.ReadVar("r", r));
    vector<double> rvec(r.size());
    for (int ir=0; ir<r.size(); ir++)
      rvec[ir] = r(ir);
    SimpleGrid rgrid;
    rgrid.Init (rvec);
    int natom = in.CountSections("atom");

    // Read muffin-tin radius
    double radius;
    assert (in.ReadVar("radius", radius));
    for (int iat=0; iat<natom; iat++, atom++) {
      AtomRadii(atom) = radius;
      assert (in.OpenSection("atom", iat));
      int num_l = in.CountSections("APW_l");
      for (int l=0; l<num_l; l++) {
	assert (in.OpenSection("APW_l",l));
	Array<double,2> u_array;
	assert (in.ReadVar("u_of_r", u_array));
	Array<double,1> du_final;
	assert (in.ReadVar("du_dr_final", du_final));
	vector<double> u1(u_array.extent(1));
	APWorder(atom) = u_array.extent(0);
	for (int ord=0; ord<u_array.extent(0); ord++) {
	  for (int ir=0; ir<u_array.extent(1); ir++)
	    u1[ir] = u_array(ord,ir);
	  u(atom, l, ord).Init(rgrid,u1,1.0e200,du_final(ord));
	  fprintf (stderr, "rgrid[0] = %1.10e\n", rgrid[0]);
	  fprintf (stderr, "ul[0] = %1.10e\n", u1[0]);
	  //u(atom, l, ord).Init(rgrid,u1);
	}
	// Now read matching coefficients
	assert (in.OpenSection("spin"));
	int numk = in.CountSections();
	for (int ik=0; ik<numk; ik++) {
	  assert (in.OpenSection("twist", ik));
	  Array<double,3> match_real, match_imag;
	  assert (in.ReadVar("match_real", match_real));
	  assert (in.ReadVar("match_imag", match_imag));
	  for (int m=-l; m<=l; m++) {
	    int lm = lm_index(l,m);
	    for (int ord=0; ord<match_real.extent(1); ord++)
	      for (int ig=0; ig<match_real.extent(2); ig++)
		MatchCoefs(atom, lm, ord, ik, ig) =
		  complex<double>(match_real(m+l,ord,ig),
				  match_imag(m+l,ord,ig));
	  }
	  in.CloseSection(); // "twist"
	}
	in.CloseSection(); // "spin"
	in.CloseSection(); // "APW_l"
      }
      // Read local radial functions, v(r)
      int numLocal = in.CountSections("local");
      for (int iloc=0; iloc<numLocal; iloc++) {
	assert (in.OpenSection("local",iloc));
	Array<double,1> vArray;
	assert (in.ReadVar("v_of_r", vArray));
	vector<double> vVec(vArray.size());
	for (int i=0; i<vVec.size(); i++)
	  vVec[i] = vArray(i);
	v(atom, iloc).Init (rgrid, vVec);
	in.CloseSection(); // "local"
      }
      in.CloseSection(); // "atom"
    }
    in.CloseSection(); // "species"
  }    
  in.CloseSection(); // "radialfuncs"
  
  // Now read eigenvector coefficients
  assert (in.OpenSection("eigenstates"));
  assert (in.ReadVar("local_atom", LocalAtom));
  assert (in.ReadVar("local_ilo", Local_ilo));
  assert (in.ReadVar("local_l", Local_l));
  assert (in.ReadVar("local_m", Local_m));
  assert (in.ReadVar("local_species", LocalSpecies));
  

  assert (in.OpenSection("spin"));
  int numk = in.CountSections("twist");
  NumG.resize(maxk);
  
  GVecs.resize(maxk, maxG);
  int maxBands = 0;
  kVecs.resize(numk);
  NumBands.resize(numk);
  for (int ik=0; ik<numk; ik++) {
    assert (in.OpenSection("twist", ik));
    Array<double,1> twist_array;
    assert (in.ReadVar("twist_angle", twist_array));
    Vec3 twist(twist_array(0), twist_array(1), twist_array(2));
    kVecs(ik) = Lattice.Twist2k (twist);
    NumBands(ik) = in.CountSections("band");
    maxBands = max(maxBands, in.CountSections("band"));
    Array<double,2> gvecs;
    assert (in.ReadVar("gvecs_cart", gvecs));
    NumG(ik) = gvecs.extent(0);
    for (int ig=0; ig<gvecs.extent(0); ig++)
      for (int dim=0; dim<3; dim++)
	GVecs(ik,ig)[dim] = gvecs(ig,dim);
    in.CloseSection(); // "twist"
  }
  // Resize the coefficients and read them in
  APWcoefs.resize(maxk, maxBands, maxG);
  int numLocal = LocalAtom.size();
  LocalCoefs.resize(maxk, maxBands, numLocal);
  for (int ik=0; ik<numk; ik++) {
    assert (in.OpenSection("twist", ik));
    Array<double,2> phi_apw, phi_local;
    int num_bands = in.CountSections("band");
    for (int iband=0; iband<num_bands; iband++) {
      assert (in.OpenSection("band", iband));
      assert (in.ReadVar("phi_apw", phi_apw));
      double nrm = 0.0;
      for (int ig=0; ig<phi_apw.extent(0); ig++)
	for (int j=0; j<2; j++)
	  nrm += phi_apw(ig,j)*phi_apw(ig,j);
      // fprintf (stderr, "Norm(%d,%d)=%1.8e\n", ik, iband, nrm);
      assert (in.ReadVar("phi_local", phi_local));
      for (int ig=0; ig<phi_apw.extent(0); ig++)
	APWcoefs(ik, iband, ig) = 
	  complex<double>(phi_apw(ig,0), phi_apw(ig,1));
      for (int iloc=0; iloc<numLocal; iloc++) {
	LocalCoefs(ik, iband, iloc) =
	  complex<double>(phi_local(iloc,0), phi_local(iloc,1));
	//fprintf (stderr, "%1.12e + %1.12e\n", phi_local(iloc,0), phi_local(iloc,1));
      }
      in.CloseSection(); // "band"
    }
    in.CloseSection(); // "twist"
  }
  in.CloseSection(); // "spin"
  in.CloseSection(); // "eigenstates"
  in.CloseFile();
}

complex<double>
LAPWClass::evalMuffinTin(int ik, int iband, int atom, Vec3 r, Vec3 L)
{
  complex<double> val;
  val = complex<double>();

  double rmag = sqrt(dot(r,r));
  Vec3 rhat = (1.0/rmag)*r;

  // First, do APW part
  for (int l=0; l<=lMax; l++) 
    for (int ord=0; ord<APWorder(atom); ord++) {
      double uval = u(atom,l,ord)(rmag);
      for (int m=-l; m<=l; m++) {
	complex<double> ylm(Ylm(l,m,rhat));
	int lm = lm_index(l,m);
	for (int ig=0; ig<NumG(ik); ig++) 
	  val += (APWcoefs(ik,iband,ig) * MatchCoefs(atom,lm,ord,ik,ig)
		  * ylm) * uval;
      }
    }

  // Next, do local part
  int numLocal = LocalAtom.size(); 
  for (int iloc=0; iloc<numLocal; iloc++) {
    if (LocalAtom(iloc) == atom) {
      int l = Local_l(iloc);
      int m = Local_m(iloc);
      int ilo = Local_ilo(iloc);
      complex<double> ylm(Ylm(l,m,rhat));
      double vval = v(atom,ilo)(rmag);
      val += (LocalCoefs(ik, iband, iloc) * ylm)*vval;
    }
  }
  double s,c;
  double phase = dot(kVecs(ik),L);
  sincos(phase, &s, &c);
  complex<double> z(c,s);

  return z*val;
}

complex<double>
LAPWClass::evalInterstitial (int ik, int iband, Vec3 r)
{
  Vec3 gpk;
  double phase, s, c;
  complex<double> val = complex<double>();
  for (int ig=0; ig<NumG(ik); ig++) {
    gpk = GVecs(ik,ig) + kVecs(ik);
    phase = dot(r,gpk);
    sincos(phase, &s, &c);
    val += complex<double>(c,s)*APWcoefs(ik,iband,ig);
  }
  return PWNorm*val;
}


void
LAPWClass::SetupMuffinTins(vector<MuffinTinClass>& tins)
{
  tins.resize(AtomPos.size());
  int totalBands = 0;
  for (int ik=0; ik<NumBands.size(); ik++) 
    totalBands += NumBands(ik);

  for (int itin=0; itin<AtomPos.size(); itin++) {
    tins[itin].Center = AtomPos(itin);
    tins[itin].APWRadius = AtomRadii(itin);
    tins[itin].lMax   = lMax;
    tins[itin].NumOrbitals = totalBands;

    ///////////////////
    // Setup splines //
    ///////////////////
    int numYlm = (lMax+1)*(lMax+1);
    tins[itin].YlmVec.resize(numYlm);
    tins[itin].dYlmVec.resize(numYlm);
    int numSplines = totalBands*numYlm;
    tins[itin].RadialVec.resize(numSplines);
    tins[itin].dRadialVec.resize(numSplines);
    tins[itin].d2RadialVec.resize(numSplines);

    tins[itin].kPoints.resize(totalBands);
    

    Ugrid rgrid;
    rgrid.start = 0.0;
    rgrid.end = AtomRadii(itin);
    rgrid.num = 100;
    BCtype_z rBC;
    rBC.lCode = NATURAL;
    // This can be improved later with the exact first derivative
    rBC.rCode = NATURAL;
    cerr << "numSplines = " << numSplines << endl;
    // Create spline object
    multi_UBspline_1d_z *spline = 
      create_multi_UBspline_1d_z (rgrid, rBC, numSplines);
    tins[itin].RadialSplines = spline;


    int numLocal = LocalAtom.size();
    // Now, loop over orbitals and Ylms, to create splines
    Array<complex<double>,1> phi_r(rgrid.num);
    int spline_num = 0;
    for (int ik=0; ik<kVecs.size(); ik++) {
      tins[itin].kPoints[ik] = kVecs(ik);
      for (int band=0; band<NumBands(ik); band++) {
	int lm = 0;
	for (int l=0; l<=lMax; l++) {
	  for (int m=-l; m<=l; m++,lm++) {
	    int lm = lm_index(l,m);
	    for (int ir=0; ir<rgrid.num; ir++) {
	      double x = (double)ir/(double)(rgrid.num-1);
	      double r = (1.0-x)*rgrid.start + x*rgrid.end;
	      phi_r(ir) = 0.0;
	      // Compute APW contribution
	      for (int ord=0; ord<APWorder(itin); ord++) {
		double uval = u(itin,l,ord)(r);
		for (int ig=0; ig<NumG(ik); ig++) {
		  phi_r(ir) += APWcoefs(ik,band,ig)*
		    MatchCoefs(itin,lm, ord, ik, ig)*uval;
		}
	      }
	      // Compute local orbital contribution
	      for (int iloc=0; iloc<numLocal; iloc++) {
		if ((LocalAtom(iloc) == itin) &&
		    (Local_l(iloc)   == l   ) &&
		    (Local_m(iloc)   == m   )) {
		  int ilo = Local_ilo(iloc);
		  phi_r(ir) += LocalCoefs(ik,band,iloc) * v(itin,ilo)(r);
		}
	      }
	    }
	    // Set the spline data
	    set_multi_UBspline_1d_z (spline, spline_num, phi_r.data());
	    spline_num++;
	  }
	}
      }
    }
  }

}

#include <time.h>

main(int argc, char **argv)
{
  LAPWClass lapw;
  vector<MuffinTinClass> tins;
  lapw.Read (argv[1]);
  lapw.SetupMuffinTins(tins);

  Vec3 r(0.1, 0.2, 0.3);
  int N = tins[0].getNumOrbitals();
  vector<complex<double> > phi(N), lapl(N), laplFD(N);
  vector<TinyVector<complex<double>,3> > grad(N), gradFD(N);

  tins[0].evaluate(r, phi);
  fprintf (stderr, "phi(0.1,0.2,0.3) = %16.12e + %16.12e\n",
	   phi[0].real(), phi[0].imag());
  tins[0].evaluate(r, phi, grad, lapl);
  tins[0].evaluateFD(r, phi, gradFD, laplFD);

  fprintf (stderr, "Re(grad[0])   = (%16.12f, %16.12f, %16.12f)\n",
	   grad[0][0].real(), grad[0][1].real(), grad[0][2].real());
  fprintf (stderr, "Re(gradFD[0]) = (%16.12f, %16.12f, %16.12f)\n",
	   gradFD[0][0].real(), gradFD[0][1].real(), gradFD[0][2].real());

  fprintf (stderr, "lapl   = %16.12f + %16.12ei\n", lapl[0].real(), lapl[0].imag());
  fprintf (stderr, "laplFD = %16.12f + %16.12ei\n", laplFD[0].real(), laplFD[0].imag());


  clock_t start, end;
  start = clock();
  for (int i=0; i<100000; i++) {
    r = Vec3(drand48(), drand48(), drand48());
    tins[0].evaluate(r, phi);
  }
  end = clock();
  cerr << "Orbital evaluations per second = " 
       << ((double)CLOCKS_PER_SEC*100000 * 80)/(double)(end-start) << endl;
  
  bool inTin;
  complex<double> val1, val2;
  val1 = lapw(0,3,r, inTin);
  val2 = phi[3];
  fprintf (stderr, "val1 = %18.14e + %18.14e\n", val1.real(), val1.imag());
  fprintf (stderr, "val2 = %18.14e + %18.14e\n", val2.real(), val2.imag());
  cerr << "val1 = " << val1 << endl;
  cerr << "val2 = " << val2 << endl;


  Vec3 p1(0.0, 0.0, 0.0);
  Vec3 p2(6.8408,6.8408,6.8408);
  
//   FILE *fout = fopen ("plot_phi.dat", "w");
//   for (double a=0.0; a<=1.0; a+=0.0001000003) {
//     Vec3 r = (1.0-a)*p1 + a*p2;
//     fprintf (fout, "%1.8f ", a);
//     for (int band=0; band<8; band++) {
//       complex<double> phi = lapw(0,band,r, inTin);
//       fprintf (fout, "%1.14e %1.14e ", real(phi), imag(phi));
//     }
//     fprintf (fout, "%1.1f", inTin ? 1.0 : 0.0);
//     fprintf (fout, "\n");
//   }
//   fclose(fout);


// #pragma omp parallel
//   for (int band=0; band<8; band++) {
//     char fname[20];
//     snprintf (fname, 20,"band%d.dat", band);
//     FILE* fout2 = fopen (fname, "w");
//     int ik = 0;
//     for (double a=0.0; a<=1.0; a+=0.003000003) {
//       for (double b=0.0; b<=1.0; b+=0.003000003) {
//   	// Vec3 r = a*lapw.Lattice.a(0) + b*lapw.Lattice.a(1) 
//   	// 	+ 0.5*lapw.Lattice.a(2);
//   	Vec3 r(6.8408*a, 6.8408*b, 0.25*6.8408);
//   	complex<double> phi = lapw(ik,band,r, inTin);
//   	fprintf (fout2, "%1.14e ", real(phi));
//       }
//       fprintf (fout2, "\n");
//     }
//     fclose(fout2);
//   }
}
