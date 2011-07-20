#include "OrbitalSetClass.h"
#include <libxml/xmlreader.h>

inline Mat3 ToMat3 (string a, string b, string c)
{
  Mat3 cell;
  istringstream as(a), bs(b), cs(c);
  as >> cell(0,0) >> cell(0,1) >> cell(0,2);
  bs >> cell(1,0) >> cell(1,1) >> cell(1,2);
  cs >> cell(2,0) >> cell(2,1) >> cell(2,2);
  return cell;
}

inline Vec3 ToVec3 (string a)
{
  Vec3 v;
  istringstream as(a);
  
  as >> v[0] >> v[1] >> v[2];
  return v;
}


inline Int3 ToInt3 (string a)
{
  Int3 v;
  istringstream as(a);
  
  as >> v[0] >> v[1] >> v[2];
  return v;
}


class Base64
{
private:
  char decode[256];
  char encode[64];
public:
  Base64() {
    for (int j=0; j<256; j++)
      decode[j] = 65;
    unsigned char i = 0;
    for (unsigned char ch='A'; ch <= 'Z'; ch++) 
      decode[ch] = i++;
    for (unsigned char ch='a'; ch <= 'z'; ch++)
      decode[ch] = i++;
    for (unsigned char ch='0'; ch <= '9'; ch++)
      decode[ch] = i++;
    decode['+'] = i++;
    decode['/'] = i++;
  }
  void Decode (string in, vector<double> &out) 
  {
    union {
      char cout[8];
      double x;
    } buff;
    int cindex=0;
    int off = 0;
    unsigned int val;
    cerr << "in.size() = " << in.size() << endl;
    for (int i=0; i<in.size(); i++) {
      unsigned char ch = decode[in[i]];
      if (ch < 65) {
	val |= ch;
	if (off < 3) {
	  val <<= 6;
	  off++;
	}
	else {
	  off = 0;
	  buff.cout[cindex++] = (val>>16) & 0xff;
	  if (cindex == 8) {
	    out.push_back (buff.x);
	    cindex=0;
	  }
	  buff.cout[cindex++] = (val>>8)  & 0xff;
	  if (cindex == 8) {
	    out.push_back (buff.x);
	    cindex=0;
	  }
	  buff.cout[cindex++] = (val>>0)  & 0xff;
	  if (cindex == 8) {
	    out.push_back (buff.x);
	    cindex=0;
	  }
	  val = 0;
	}
      }
    }
    for (int k=0; k<off; k++) {
      buff.cout[cindex++] = (val>>(8*(2-k))) & 0xff;
      if (cindex == 8) {
	out.push_back (buff.x);
	cindex=0;
      }
    }
    
  }
};
      


bool
OrbitalSetClass::Read_FPMD (string fname)
{
  xmlDocPtr doc;
  Mat3 Asuper;
  Base64 decoder;
  vector<vector<xmlNodePtr> > grid_functions;
  vector<Vec3> kPoints;
  CellClass tempCell;


  doc = xmlParseFile (fname.c_str());
  xmlNodePtr cur = xmlDocGetRootElement(doc);
  assert(!xmlStrcmp(cur->name, (const xmlChar *)"sample"));
  cur = cur->xmlChildrenNode;
  Int3 FFTgrid;
  // First pass, read parameters
  while (cur != NULL) {
    if (!xmlStrcmp(cur->name, (const xmlChar *)"atomset")) {
      xmlNodePtr ascur = cur->xmlChildrenNode;
      while (ascur != NULL) {
	//cerr << "ascur->name = " << ascur->name << endl;
	if (!xmlStrcmp(ascur->name, (const xmlChar *)"unit_cell")) {
	  string a = (const char *) xmlGetProp(ascur, (const xmlChar*)"a");
	  string b = (const char *) xmlGetProp(ascur, (const xmlChar*)"b");
	  string c = (const char *) xmlGetProp(ascur, (const xmlChar*)"c");
	  PrimCell.SetLattice(ToMat3(a,b,c));
	  tempCell.SetLattice(ToMat3(a,b,c));
	  Asuper = TileMatrix*(ToMat3(a,b,c));
	  SuperCell.SetLattice(Asuper);
	}
	if (!xmlStrcmp(ascur->name, (const xmlChar *)"atom")) {
	  xmlNodePtr atom = ascur->xmlChildrenNode;
	  while (atom != NULL) {
	    if (!xmlStrcmp(atom->name, (const xmlChar *)"position")) {
	      string pos = (const char*)xmlNodeListGetString 
		(doc, atom->xmlChildrenNode, 1);
	      int n = PrimCell.IonPos.size();
	      PrimCell.IonPos.resizeAndPreserve(n+1);
	      PrimCell.AtomTypes.resizeAndPreserve(n+1);
	      PrimCell.IonPos(n) = ToVec3(pos);
	      PrimCell.AtomTypes(n) = 1;
	      //cerr << "IonPos = " << PrimCell.IonPos(n) << endl;
	    }
	    atom = atom->next;
	  }
	}
	ascur = ascur->next;
      }
    }
    
    if (!xmlStrcmp(cur->name, (const xmlChar *)"wavefunction")) {
      ECut = atof((char *)xmlGetProp(cur, (const xmlChar*)"ecut"));
      xmlNodePtr wfcur = cur->xmlChildrenNode;
      while (wfcur != NULL) {
	//cerr << "  -- name = " << wfcur->name << endl;
	if (!xmlStrcmp(wfcur->name, (const xmlChar *)"grid")) {
	  FFTgrid[0] = 
	    atoi((const char*) xmlGetProp(wfcur, (const xmlChar*)"nx"));
	  FFTgrid[1] = 
	    atoi((const char*) xmlGetProp(wfcur, (const xmlChar*)"ny"));
	  FFTgrid[2] = 
	    atoi((const char*) xmlGetProp(wfcur, (const xmlChar*)"nz"));
	  //cerr << "FFTgrid = " << FFTgrid << endl;
	}
	if (!xmlStrcmp(wfcur->name, (const xmlChar *)"slater_determinant")) {
	  string kstring = (char*)xmlGetProp(wfcur, (const xmlChar*)"kpoint");
	  kPoints.push_back(ToVec3(kstring));
	  vector<xmlNodePtr> grid_funcs;
	  xmlNodePtr sdcur = wfcur->xmlChildrenNode;
	  while (sdcur != NULL) {
	    if (!xmlStrcmp(sdcur->name, 
			   (const xmlChar *)"grid_function")) {
	      grid_funcs.push_back(sdcur);
	      // 	      string wfstring = (const char *) xmlNodeListGetString
	      // 		(doc, sdcur->xmlChildrenNode, 1);
	      // 	      vector<double> data;
	      // 	      decoder.Decode(wfstring, data);
	      // 	      cerr << "data.size() = " << data.size() << endl;
	      // 	      FILE *fout = fopen ("test.dat", "w");
	      // 	      for (int iy=0; iy<FFTgrid[1]; iy++) {
	      // 		for (int iz=0; iz<FFTgrid[0]; iz++)
	      // 		  fprintf (fout, "%14.8f ", data[iy*FFTgrid[0]+iz]);
	      // 		fprintf (fout, "\n");
	      // 	      }
	      // 	      fclose(fout);
	    }
	    //cerr << "      -- name = " << sdcur->name << endl;
	    sdcur = sdcur->next;
	  }
	  grid_functions.push_back(grid_funcs);
	}
	
	wfcur = wfcur->next;
      }
      
    }
    cur = cur->next;
  }
  
  // Setup the FFT box with the appropriate ECut and 
  double kCut = sqrt(2.0*ECut);
  tempCell.GVecs.Set(PrimCell.Lattice.GetDirect(), Vec3(0.0,0.0,0.0), kCut, FFTgrid);
  tempCell.SetupFFT();
  
  GVecsArray.resize(tempCell.GVecs.size());
  for (int iG=0; iG<GVecsArray.size(); iG++)
    GVecsArray(iG) = tempCell.GVecs(iG);
  
  int numk     = grid_functions.size();
  int numBands =0;
  int numG = tempCell.GVecs.size();
  for (int ik=0; ik<numk; ik++)
    numBands = max(numBands, (int)grid_functions[ik].size());
  PrimOrbitals[0].resize(numk, numBands);
  
  // Loop over the grid_function sections, reading them one by one
  for (int ik=0; ik<grid_functions.size(); ik++) {
    for (int iband=0; iband<grid_functions[ik].size(); iband++) {
      xmlNodePtr func = grid_functions[ik][iband];
      string wfstring = (const char *) xmlNodeListGetString
	(doc, func->xmlChildrenNode, 1);
      vector<double> data;
      decoder.Decode(wfstring, data);
      // int N = data.size();
      // fprintf (stderr, "data[N-1] = %1.8e data[N-2] = %1.8e\n");
      int index=0;
      double nrm = 0.0;
      for (int iz=0; iz<FFTgrid[2]; iz++)
	for (int iy=0; iy<FFTgrid[1]; iy++)
	  for (int ix=0; ix<FFTgrid[0]; ix++) {
	    tempCell.FFT.rBox(ix,iy,iz) = data[index++];
	    nrm += norm (tempCell.FFT.rBox(ix,iy,iz));
	  }
      fprintf (stderr, "Real-space norm = %1.14f\n", 
	       nrm/(double)(FFTgrid[0]*FFTgrid[1]*FFTgrid[2]));
      tempCell.FFT.r2k();
      zVec &coefs = *(new zVec(numG));
      tempCell.FFT.GetkVec(coefs);
      double maxcoef = 0.0;
      for (int i=0; i<coefs.size(); i++) {
	maxcoef = max (maxcoef, norm(coefs(i)));
	if (isnan(coefs(i).real()) || isnan(coefs(i).imag()))
	  cerr << "NAN at i = " << i << endl;
      }
      //coefs *= (1.0/sqrt(norm(coefs)));
      fprintf (stderr, "G-space norm    = %1.14f\n", norm(coefs));
      //fprintf (stderr, "maxcoef = %1.8e\n", maxcoef);

      PrimOrbitals[0](ik,iband) = new OrbitalClass;
      PrimOrbitals[0](ik,iband)->SetCell(PrimCell);
      PrimOrbitals[0](ik,iband)->SetCoefs(coefs);
      //Vec3 twist = -1.0*PrimCell.Lattice.k2Twist(kPoints[ik]);
      Vec3 twist = -1.0*kPoints[ik];
      PrimOrbitals[0](ik,iband)->SetTwist(twist);
      PrimOrbitals[0](ik,iband)->SetLabels (0, ik, iband);
    }
  }

  // Now


  fprintf (stderr, "Superlattice vectors:\n");
  for (int i=0; i<3; i++)
    fprintf (stderr, "  [ %9.6f %9.6f %9.6f ]\n", 
	     Asuper(i,0), Asuper(i,1), Asuper(i,2));
  
  NumElectrons[0] = 8;//8.0;
  
  SpinPolarized = false;

  SetupFFTIndices();
  CreateSuperOrbitals();
  TileIonPos();

  return true;

  return doc != NULL;
}
