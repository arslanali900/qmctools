#include "OrbitalSetClass.h"
#include "ParseCommand.h"
#include "Gaussian.h"


Mat3
string2Mat3 (string str)
{
  Mat3 mat;
  stringstream strin;
  strin << str;
  for (int i=0; i<3; i++)
    for (int j=0; j<3; j++) {
      int val;
      strin >> val;
      mat(i,j) = (double) val;
      if (strin.fail()) {
	cerr << "Error reading tilemat string.  Aborting.\n";
	abort();
      }
    } 
  return mat;
}

Int3 
string2Int3 (string str)
{
  Int3 vec;
  stringstream strin;
  strin << str;
  for (int i=0; i<3; i++) {
    int val;
    strin >> val;
    vec[i] = val;
    if (strin.fail()) {
      cerr << "Error reading tilemat string.  Aborting.\n";
      abort();
    }
  }
  return vec;
}

void HelpMultirep()
{
  perr << "\nMultirep usage:\n"
       << "  wfconv --eshdf myoutfile.h5 --multirep myrepfile.in myABINIT.WFK --factor x\n"
       << "    or\n"
       << "  wfconv --eshdf myoutfile.h5 --multirep myrepfile.in pwfn.data    --factor x\n\n"
       << "Example myrepfile.in:\n\n"
       << "Section(Species)          \n"
       << "{                         \n"
       << "  string name = \"Fe\";   \n"
       << "  double radius = 1.35;   \n"
       << "  int spline_points = 100;\n"
       << "  int lmax = 6;           \n"
       << "}                         \n"
       << "                          \n"
       << "Section(Species)          \n"
       << "{                         \n"
       << "  string name = \"O\";    \n"
       << "  double radius = 1.1;    \n"
       << "  int spline_points = 100;\n"
       << "  int lmax = 6;           \n"
       << "}                         \n\n";
}




main(int argc, char **argv)
{
  COMM::Init (argc, argv);

  list<ParamClass> argList;

  argList.push_back(ParamClass("nospline", false));
  argList.push_back(ParamClass("real",   false));
  argList.push_back(ParamClass("check",  false));
  argList.push_back(ParamClass("factor", true));
  argList.push_back(ParamClass("localize", true));
  argList.push_back(ParamClass("radius", true));
  argList.push_back(ParamClass("shift", false));
  argList.push_back(ParamClass("optimize-centers", false));
  argList.push_back(ParamClass("optimize-radii", false));
  argList.push_back(ParamClass("ortho", false));
  argList.push_back(ParamClass("truncate", false));
  argList.push_back(ParamClass("unfold",false));
  argList.push_back(ParamClass("bwfn", true));
  argList.push_back(ParamClass("split", false));
  argList.push_back(ParamClass("pwfn", true));
  argList.push_back(ParamClass("qmcPACK", true));
  argList.push_back(ParamClass("eshdf", true));
  argList.push_back(ParamClass("skin-thickness", true));
  argList.push_back(ParamClass("fold-factor", 3));
  argList.push_back(ParamClass("tile", 3));
  argList.push_back(ParamClass("tilemat", 1));
  argList.push_back(ParamClass("density", 1));
  argList.push_back(ParamClass("vhxc", 1));
  argList.push_back(ParamClass("first_order", 1));
  argList.push_back(ParamClass("first_order_FD", 1));
  argList.push_back(ParamClass("multirep", 1));
  argList.push_back(ParamClass("compare", false));
  argList.push_back(ParamClass("gamess",false));
  argList.push_back(ParamClass("help", 1));
  argList.push_back(ParamClass("interp", false));
  argList.push_back(ParamClass("laplvl", true));


  CommandLineParserClass parser (argList);
  bool success = parser.Parse (argc, argv);

  if (parser.Found("help")) {
    string topic = parser.GetArg("help",0);
    if (topic == "multirep") {
      HelpMultirep();
      exit(-1);
    }
  }
  
  if (!success || parser.NumFiles() != 1) {
    perr << "Usage:  wfconv [options...] infile\n"
	 << "Options:\n"
	 << "  --localize fname        construct localized linear combinations of\n"
	 << "        or random         the occupied bands\n"
	 << "  --radius x              Radius of localization\n"
	 << "  --skin-thickness x      Skin thickness for smooth truncation\n"
	 << "  --nospline              Do not Fourier transform orbitals to real space\n"
	 << "  --real                  take the real part of spline data\n"
         << "  --factor x              Multiply the number of real-space mesh points\n" 
	 << "                          by x in each direction.  x can be fractional\n"
	 << "  --shift                 shift localized and spline orbitals to the \n"
	 << "                          center of the simulation cell\n"
	 << "  --optimize-centers      optimize the localization centers\n" 
	 << "  --optimize-radii x      optimize the localization radii to contain x\n"
	 << "                          fraction of the norm\n"
	 << "  --ortho                 orthogonalize the local orbitals\n"
	 << "  --truncate              truncate localized orbitals beyond localization\n"
	 << "                          radius\n"
	 << "  --unfold                unfold k-point grid into a single k-point  supercell\n"
	 << "  --fold-factor nx ny nz  Partially unfold by factor nx x ny x nz.  Otherwise \n"
         << "                          unfold completely to a single k-point.\n"
	 << "  --tile nx ny nz         Tile the primitive cell by nx x ny x nz.\n"
         << "                          Do not really unfold, just update labeling and twists.\n"
         << "  --tilemat \"S11 S12 S13 S21 S22 S23 S31 S32 S33\"\n"
	 << "                          Tile the primitive cell by the integer matrix S.\n"
	 << "   --density fname        Read the density from file fname.\n"
	 << "   --vhxc fname           Read the VHXC potential from file fname.\n"
	 << "   --first_order fbase    Read ABINIT first-order wavefunctions.\n"
	 << "                          fbase is the file name without the trailing numbers.\n"
         << "   --multirep fname       Use atomic/3D B-spline representation.  The atomic\n"
	 << "                          orbitals should be described in fname"
	 << "   --compare              compare hybrid representation to standard B-spline.\n"
	 << "  --bwfn fname            write CASINO blip wave function file\n"
	 << "  --pwfn fname            write CASINO plane wave function file\n"
	 << "  --interp                with --bwfn, use interpolating coefficients\n"
	 << "  --qmcPACK fname         write qmcPACK wave function file\n"     
    << "  --eshdf fname           write ES-HDF wave function file\n"     
    << "  --laplvl x              write ES-HDF wave function laplacians to level x\n";     

    exit(1); 
  }
    


  string inName = parser.GetFile(0);
  double factor(1.0), radius(3.0);
  int lappwr(0);
  bool real     = parser.Found("real");
  bool checkKE  = parser.Found("check");
  bool shift    = parser.Found("shift");
  if (parser.Found("factor")) {
    factor = atof (parser.GetArg("factor").c_str());
    perr << "factor = " << factor << endl;
  }
  if (parser.Found("radius")) {
    radius = atof (parser.GetArg("radius").c_str());
    perr << "radius = " << radius << endl;
  }
  if (parser.Found("laplvl")) {
    lappwr = atoi (parser.GetArg("laplvl").c_str());
//     perr << "laplvl = " << lappwr << endl;
  }
  
  
  OrbitalSetClass system;
  system.Localized       = parser.Found("localize");
  system.Spline          = !parser.Found("nospline");
  system.OptimizeCenters = parser.Found("optimize-centers");
  system.OptimizeRadii   = parser.Found("optimize-radii");
  system.ShiftOrbitals   = parser.Found("shift");
  system.Real            = parser.Found("real");
  system.CheckKE         = parser.Found("check");
  system.Orthogonalize   = parser.Found("ortho");
  system.Truncate        = parser.Found("truncate");
  system.UseMultiRep     = parser.Found("multirep");
  system.CompareHybrid   = parser.Found("compare");
  system.WriteInterp     = parser.Found("interp");

  bool unfold            = parser.Found("unfold");
  system.SetFFTFactor (factor);
  
  Int3 foldFactor (0,0,0), tileFactor(0,0,0);
  if (parser.Found("tilemat")) 
    system.SetTileMatrix (string2Mat3 (parser.GetArg("tilemat",0)));
  else if (parser.Found("tile")) {
    unfold = true;

    Int3 tile = string2Int3 (parser.GetArg("tile",0));
    Mat3 tileMat;
    tileMat = 0.0;
    tileMat(0,0) = (double)tile[0];
    tileMat(1,1) = (double)tile[1];
    tileMat(2,2) = (double)tile[2];
    system.SetTileMatrix(tileMat);
  }
     
  if (system.ShiftOrbitals && !system.Localized) {
    perr << "Orbitals must be localized to be shifted.\n";
    abort();
  }
  if (system.OptimizeCenters && !system.Localized) {
    perr << "Cannot optimize localization centers without localizing.\n";
    abort();
  }
  if (system.Truncate && !system.ShiftOrbitals) {
    perr << "Orbitals must be shifted to be truncated.\n";
    abort();
  }
  

  system.SetFFTFactor (factor);
  if (inName.find("WFK") != std::string::npos) {
    system.Read_ABINIT_WFK(inName);
    cerr << "Finished reading ABINIT WFK file.\n";
  }
  else if (inName.find("1WF") != std::string::npos) {
    system.Read_ABINIT_WF1(inName);
    cerr << "Finished reading ABINIT WF1 file.\n";
  }
  else if (inName.find("xml") != std::string::npos) {
    system.Read_FPMD(inName);
    cerr << "Finished reading FPMD file.\n";
  }  
  else if (inName.find(".h5") != std::string::npos) {
    if (system.Read_ESHDF(inName)) 
      perr << "Successfully read ESHDF file \"" << inName << "\".\n";
    else if (!system.Read_LAPW(inName)) {
      perr << "Error in reading LAPW file " << inName << ".  Aborting.\n";
      abort();
    }
  }
  else if (parser.Found("gamess")) {
    GaussianOrbitalSet gauss;
    gauss.ReadGAMESS (inName);
    if (!parser.Found("eshdf")) 
      cerr << "GAMESS files can only be converted to ESHDF format.\n";
    else
      gauss.Write_ESHDF (parser.GetArg("eshdf"));
    exit(0);
  }
  else if (!system.Read(inName)) {
    perr << "Could not read file " << inName << ".  Aborting.\n";
    abort();
  }
      
  if (parser.Found("density")) {
    system.Read_ABINIT_DEN(parser.GetArg("density"));
    perr << "Read ABINIT density file.\n";
  }
  if (parser.Found("vhxc")) {
    system.Read_ABINIT_POT(parser.GetArg("vhxc"));
    perr << "Read ABINIT VHXC potential file.\n";
  }


  if (parser.Found("first_order")) {
    system.Read_ABINIT_First_Order (parser.GetArg("first_order"));
  }

  if (parser.Found("first_order_FD")) {
    system.Read_ABINIT_First_Order_FD (parser.GetArg("first_order_FD"));
  }


  if (parser.Found("localize")) {
    string name = parser.GetArg("localize");
    clock_t start = clock();
    system.Localize(name, radius);
    clock_t end = clock();
    perr << "Total time = " << (end-start)/CLOCKS_PER_SEC << " seconds.\n";
  }

  if (parser.Found("skin-thickness")) {
    double thickness = atof (parser.GetArg("skin-thickness").c_str());
    perr << "skin thickness = " << thickness << endl;
    system.SetSkinThickness (thickness);
  }

  if (system.UseMultiRep)
    system.ReadMultiRep (parser.GetArg("multirep"));


  if (parser.Found ("bwfn"))
    system.WriteBWFN (parser.GetArg ("bwfn"), parser.Found("split"));
  if (parser.Found ("pwfn"))
    system.WritePWFN (parser.GetArg ("pwfn"));
  if (parser.Found("qmcPACK"))
    system.Write_qmcPACK (parser.GetArg ("qmcPACK"));
  if (parser.Found("eshdf"))
    system.Write_ESHDF (parser.GetArg ("eshdf"),lappwr);

  COMM::Finalize();
}

