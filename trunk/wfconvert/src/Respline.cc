#include "Respline.h"
#include "ParseCommand.h"

void
OrbSplineClass::Read(IO::IOSectionClass &in)
{
  Array<double,4> rData;
  in.ReadVar ("eigenvector", rData);
  Nx = rData.extent(0);
  Ny = rData.extent(1);
  Nz = rData.extent(2);
  Data.resize (Nx, Ny, Nz);
  for (int ix=0; ix<Nx; ix++)
    for (int iy=0; iy<Ny; iy++)
      for (int iz=0; iz<Nz; iz++)
	Data(ix,iy,iz) = complex<double>(rData(ix,iy,iz,0),
					 rData(ix,iy,iz,1));
  Array<double,1> dummy;
  in.ReadVar("umin", dummy);
  assert (dummy.size() == 3);
  uMin = Vec3(dummy(0), dummy(1), dummy(2));

  in.ReadVar("umax", dummy);
  assert (dummy.size() == 3);
  uMax = Vec3(dummy(0), dummy(1), dummy(2));

  // Create a spline
  if (OldSpline != NULL) 
    destroy_Bspline(OldSpline);
  
  
  Ugrid xgrid, ygrid, zgrid;
  BCtype_z xBC, yBC, zBC;
  xBC.lCode = xBC.rCode = NATURAL;
  yBC.lCode = yBC.rCode = NATURAL;
  zBC.lCode = zBC.rCode = NATURAL;

  xgrid.num = Nx; xgrid.start = 0.0; xgrid.end = 1.0000000001;
  ygrid.num = Ny; ygrid.start = 0.0; ygrid.end = 1.0000000001;
  zgrid.num = Nz; zgrid.start = 0.0; zgrid.end = 1.0000000001;
  OldSpline = create_UBspline_3d_z (xgrid, ygrid, zgrid,
				    xBC, yBC, zBC,
				    Data.data());
  in.ReadVar ("eigenvalue", Eigenvalue);
  in.ReadVar ("radius", Radius);
  in.ReadVar ("spin", Spin);
  Array<double,1> tmp;
  in.ReadVar ("umin", tmp);
  uMin = Vec3 (tmp(0), tmp(1), tmp(2));
  in.ReadVar ("umax", tmp);
  uMax = Vec3 (tmp(0), tmp(1), tmp(2));
  in.ReadVar ("center", tmp);
  Center = Vec3(tmp(0), tmp(1), tmp(2));
  in.ReadVar ("truncated", Truncated, false);

  Array<double,2> centers;
  in.ReadVar("centers"   , centers);
  if (centers.extent(0) != 0) {
    Centers.resize(centers.extent(0));
    for (int i=0; i<centers.extent(0); i++)
      for (int j=0; j<3; j++)
	Centers(i)[j] = centers(i,j);
  }
}

void
OrbSplineClass::GetDims (int &nx, int &ny, int &nz)
{
  nx=Nx; ny=Ny; nz=Nz;
}


void
OrbSplineClass::Write (IO::IOSectionClass &out)
{
  int nx = Nonuniform->x_grid->num_points;
  int ny = Nonuniform->y_grid->num_points;
  int nz = Nonuniform->z_grid->num_points;
  fprintf (stderr, "nx=%d ny=%d nz=%d\n", nx, ny, nz);

  Array<double,1> points;
  out.NewSection ("xgrid");
  out.WriteVar ("clusterfactor", ClusterFactor);
  out.WriteVar ("numpoints", nx);
  points.resize(nx);
  for (int ix=0; ix<nx; ix++)
    points(ix) = Nonuniform->x_grid->points[ix];
  out.WriteVar("points", points);
  out.CloseSection();

  out.NewSection ("ygrid");
  out.WriteVar ("clusterfactor", ClusterFactor);
  out.WriteVar ("numpoints", ny);
  points.resize(ny);
  for (int iy=0; iy<ny; iy++)
    points(iy) = Nonuniform->y_grid->points[iy];
  out.WriteVar("points", points);
  out.CloseSection();

  out.NewSection ("zgrid");
  out.WriteVar ("clusterfactor", ClusterFactor);
  out.WriteVar ("numpoints", nz);
  points.resize(nz);
  for (int iz=0; iz<nz; iz++)
    points(iz) = Nonuniform->z_grid->points[iz];
  out.WriteVar("points", points);
  out.CloseSection();

  Array<double,4> rData (nx, ny, nz, 2);
  for (int ix=0; ix<nx; ix++)
    for (int iy=0; iy<ny; iy++)
      for (int iz=0; iz<nz; iz++) {
	rData(ix, iy, iz, 0) = NonuniformData(ix,iy,iz).real();
	rData(ix, iy, iz, 1) = NonuniformData(ix,iy,iz).imag();
      }
  out.WriteVar("eigenvector", rData);
  
  // Write other data
  out.WriteVar("eigenvalue", Eigenvalue);
  out.WriteVar("radius"    , Radius);
  out.WriteVar("spin"      , Spin);
  out.WriteVar("umin"      , uMin);
  out.WriteVar("umax"      , uMax);
  out.WriteVar("center"    , Center);
  Array<double,2> centers(Centers.size(), 3);
  for (int i=0; i<centers.extent(0); i++) 
    for (int j=0; j<3; j++)
      centers(i,j) = Centers(i)[j];
  out.WriteVar("centers"   , centers);
  out.WriteVar("truncated" , Truncated);

}

void
OrbSplineClass::Respline (int nx, int ny, int nz,
			  double clusterFactor)
{
  ClusterFactor = clusterFactor;
  // Create grids
  if (xNUgrid != NULL)  destroy_grid (xNUgrid);
  if (yNUgrid != NULL)  destroy_grid (yNUgrid);
  if (zNUgrid != NULL)  destroy_grid (zNUgrid);

  xNUgrid = create_center_grid (0.0, 1.0, clusterFactor, nx);
  yNUgrid = create_center_grid (0.0, 1.0, clusterFactor, ny);
  zNUgrid = create_center_grid (0.0, 1.0, clusterFactor, nz);
  
  BCtype_z xBC, yBC, zBC;
  xBC.lCode = xBC.rCode = NATURAL;
  yBC.lCode = yBC.rCode = NATURAL;
  zBC.lCode = zBC.rCode = NATURAL;

//   for (int ix=0; ix<nx; ix++)
//     fprintf (stderr, "xNUgrid[%d] = %1.12f\n", ix, xNUgrid->points[ix]);

  NonuniformData.resize(nx, ny, nz);
  for (int ix=0; ix<nx; ix++) {
    double x = xNUgrid->points[ix];
    for (int iy=0; iy<ny; iy++) {
      double y = yNUgrid->points[iy];
      for (int iz=0; iz<nz; iz++) {
	double z = zNUgrid->points[iz];
	eval_UBspline_3d_z (OldSpline, x, y, z, &NonuniformData(ix,iy,iz));
	if (isnan(NonuniformData(ix,iy,iz).real()) || 
	    isnan(NonuniformData(ix,iy,iz).imag())) 
	  fprintf (stderr, "NAN found at (x,y,z) = (%8.5f, %8.5f, %8.5f in Respline.\n",
		   x, y, z);
      }
    }
  }

  if (Nonuniform != NULL)
    destroy_Bspline (Nonuniform);

  Nonuniform = create_NUBspline_3d_z (xNUgrid, yNUgrid, zNUgrid,
				      xBC, yBC, zBC, NonuniformData.data());
  for (int ix=0; ix<nx+2; ix++)
    for (int iy=0; iy<ny+2; iy++)
      for (int iz=0; iz<nz+2; iz++) {
	complex<double> c = Nonuniform->coefs[ix*Nonuniform->x_stride +
					      iy*Nonuniform->y_stride +
					      iz];
	if (isnan(c.real()) || isnan(c.imag()))
	  fprintf (stderr, "NAN in c(%d,%d,%d)\n", ix,iy,iz);
      }
	
}


void
OrbSplineClass::Respline (int nx, int ny, int nz)
{
  Ugrid xgrid, ygrid, zgrid;
  xgrid.num = nx; xgrid.start = 0.0; xgrid.end = 1.000000000;
  ygrid.num = ny; ygrid.start = 0.0; ygrid.end = 1.000000000;
  zgrid.num = nz; zgrid.start = 0.0; zgrid.end = 1.000000000;

  BCtype_z xBC, yBC, zBC;
  xBC.lCode = xBC.rCode = NATURAL;
  yBC.lCode = yBC.rCode = NATURAL;
  zBC.lCode = zBC.rCode = NATURAL;

  Array<complex<double>,3> newData (nx, ny, nz);
  for (int ix=0; ix<nx; ix++) {
    double x = (double)ix/(double)(nx-1);
    for (int iy=0; iy<ny; iy++) {
      double y = (double)iy/(double)(ny-1);
      for (int iz=0; iz<nz; iz++) {
	double z = (double)iz/(double)(nz-1);
	eval_UBspline_3d_z (OldSpline, x, y, z, &newData(ix,iy,iz));
	if (isnan(newData(ix,iy,iz).real()) || 
	    isnan(newData(ix,iy,iz).imag())) 
	  fprintf (stderr, "NAN found at (x,y,z) = (%8.5f, %8.5f, %8.5f in Respline.\n",
		   x, y, z);
      }
    }
  }

  if (Uniform != NULL)
    destroy_Bspline (Uniform);
  
  Uniform = create_UBspline_3d_z (xgrid, ygrid, zgrid,
				  xBC, yBC, zBC, newData.data());
}


double
OrbSplineClass::NonuniformError (int numPoints)
{
  double weight = 0.0;
  double error  = 0.0;
  for (int i=0; i<numPoints; i++) {
    double x = drand48();
    double y = drand48();
    double z = drand48();
    complex<double> zOld, zNew;
    eval_UBspline_3d_z  (OldSpline, x, y, z, &zOld);
    eval_NUBspline_3d_z (Nonuniform, x, y, z, &zNew);
//     if (isnan (zOld.real()) || isnan(zOld.imag()))
//       fprintf (stderr, "NAN in Old at (x,y,z) = (%1.8f, %1.8f, %1.8f)\n",
//  	       x,y,z);
//     else if (isnan (zNew.real()) || isnan(zNew.imag()))
//       fprintf (stderr, "NAN in New at (x,y,z) = (%1.8f, %1.8f, %1.8f)\n",
//  	       x,y,z);
      
//     else {
      double w = norm (zOld);
      error += w * norm (zOld-zNew);
      weight += w;
      //    }
  }
  return sqrt(error/weight);
}


double
OrbSplineClass::UniformError (int numPoints)
{
  double weight = 0.0;
  double error  = 0.0;
  for (int i=0; i<numPoints; i++) {
    double x = drand48();
    double y = drand48();
    double z = drand48();
    complex<double> zOld, zNew;
    eval_UBspline_3d_z (OldSpline, x, y, z, &zOld);
    eval_UBspline_3d_z (Uniform  , x, y, z, &zNew);
      double w = norm (zOld);
      error += w * norm (zOld-zNew);
      weight += w;
  }
  return sqrt(error/weight);
}


void
OrbSplineClass::Plot (string filename)
{
  FILE *fout = fopen (filename.c_str(), "w");
  double y=0.5, z=0.5;
  for (double x=0.0; x<=1.0; x+=1.0e-4) {
    complex<double> exVal, uniVal, nonuniVal;
    eval_UBspline_3d_z  (OldSpline,  x, y, z, &exVal);
    eval_UBspline_3d_z  (Uniform,    x, y, z, &uniVal);
    eval_NUBspline_3d_z (Nonuniform, x, y, z, &nonuniVal);
    fprintf (fout, "%1.8e %1.16e %1.16e %1.16e %1.16e %1.16e %1.16e\n", x,
	     exVal.real(),     exVal.imag(),
	     uniVal.real(),    uniVal.imag(),
	     nonuniVal.real(), nonuniVal.imag());
  }
  fclose (fout);
  
  fout = fopen ("OptBasis.dat", "w");
  NUBasis* basis = Nonuniform->x_basis;
  NUgrid *grid = Nonuniform->x_grid;
  for (double x=0.0; x<=1.0; x+=1.0e-4) {
    double bfuncs[4];
    get_NUBasis_funcs_d (basis, x, bfuncs);
    int i = (*grid->reverse_map)(grid, x);
    fprintf (fout, "%1.8f %1.12e %1.12e %1.12e %1.12e\n",
	     x, bfuncs[(-i+100)%4], bfuncs[(-i+101)%4], 
	     bfuncs[(-i+102)%4], bfuncs[(-i+103)%4]);
  }
  fclose (fout);
}



void
ResplineClass::CopyHeader()
{
  if (Comm.MyProc() != 0) 
    return;
  // "ions" section
  assert (In.OpenSection ("ions"));
  Out.NewSection         ("ions");
  Array<int,1> atom_types;
  assert(In.ReadVar   ("atom_types", atom_types));
  Out.WriteVar        ("atom_types", atom_types);
  Array<double,2> pos;
  assert(In.ReadVar ("pos", pos));
  Out.WriteVar      ("pos", pos);
  In.CloseSection();  // "ions"
  Out.CloseSection(); // "ions"

  // "parameters" section
  assert (In.OpenSection("parameters"));
  Out.NewSection("parameters");
  int complex_coefficients;
  assert (In.ReadVar("complex_coefficients", complex_coefficients));
  Out.WriteVar      ("complex_coefficients", complex_coefficients);
  Array<double,2> lattice;
  assert (In.ReadVar("lattice", lattice));
  Out.WriteVar      ("lattice", lattice);
  assert (In.ReadVar("reciprocal_lattice", lattice));
  Out.WriteVar      ("reciprocal_lattice", lattice);
  double maximum_ecut;
  assert (In.ReadVar("maximum_ecut", maximum_ecut));
  Out.WriteVar      ("maximum_ecut", maximum_ecut);
  int num_electrons, num_bands, num_up_bands, num_down_bands;
  assert (In.ReadVar("num_electrons",  num_electrons));
  Out.WriteVar      ("num_electrons",  num_electrons);
  assert (In.ReadVar("num_bands",      num_bands));
  Out.WriteVar      ("num_bands",      num_bands);	
  assert (In.ReadVar("num_up_bands",   num_up_bands));
  Out.WriteVar      ("num_up_bands",   num_up_bands);	
  assert (In.ReadVar("num_down_bands", num_down_bands));
  Out.WriteVar      ("num_down_bands", num_down_bands);
  int num_spins, num_twists;
  assert (In.ReadVar("num_spins",  num_spins));
  Out.WriteVar      ("num_spins",  num_spins);
  assert (In.ReadVar("num_twists", num_twists));
  Out.WriteVar      ("num_twists", num_twists);

  In.CloseSection();  // "parameters"
  Out.CloseSection(); // "parameters"

  Array<int,1> version;
  assert (In.ReadVar("version", version));
  Out.WriteVar ("version", version);

  // "basis" section
  assert (In.OpenSection("basis"));
  Out.NewSection        ("basis");
  string type; 
  assert (In.ReadVar("type", type));
  Out.WriteVar      ("type", type);
  In.CloseSection();  // "basis"
  Out.CloseSection(); // "basis"
  Out.FlushFile();
}


void 
ResplineClass::SetFiles (string inName, string outName)
{
  assert (In.OpenFile(inName));
  if (Comm.MyProc() == 0)
    assert (Out.NewFile(outName));
}


void 
ResplineClass::SetSizeFactor (double factor)
{
  SizeFactor = factor;
}

void
ResplineClass::FindOptimalClusterFactor (int twist, int band)
{
  perr << "Finding optimal factor for twist " << twist << ", band " << band << endl;
  const int numPoints = 1000000;
  // const int numPoints = 10000;
  assert (In.OpenSection ("twist", twist));
  assert (In.OpenSection ("band",  band));
  Orb.Read (In);
  In.CloseSection();  // "band"
  In.CloseSection();  // "twist"

  int nx, ny, nz, Nx, Ny, Nz;
  Orb.GetDims (Nx, Ny, Nz);
  nx = (int)round(SizeFactor*(double)Nx);
  ny = (int)round(SizeFactor*(double)Nz);
  nz = (int)round(SizeFactor*(double)Ny);

  ClusterFactor(twist,band) = 1.0;
  Orb.Respline (nx, ny, nz);
  RMSError(twist,band) = Orb.UniformError(numPoints);
		
  for (double factor=1.25; factor <= 25.0; factor+=0.25) {
    Orb.Respline (nx, ny, nz, factor);
    double error = Orb.NonuniformError(numPoints);
    if (Comm.MyProc() == 0)
      fprintf (stderr, "  factor = %5.2f   error = %11.5e\n",
	       factor, error);
    if (error < RMSError(twist, band)) {
      ClusterFactor(twist,band) = factor;
      RMSError(twist,band) = error;
    }
  }
  
}


void
ResplineClass::ProcTasks (int procNum, int numTasks, 
			  int &firstTask, int &lastTask)
{
  int numProcs = Comm.NumProcs();
  firstTask = 0;
  for (int proc=0; proc<procNum; proc++) {
    int tasks = numTasks / numProcs + (proc < (numTasks % numProcs) ? 1 : 0);
    firstTask += tasks;
  }
  int tasks =  numTasks / numProcs + (procNum < (numTasks % numProcs) ? 1 : 0);
  lastTask = firstTask + tasks - 1;
}

void
ResplineClass::DoRespline()
{
  CopyHeader();
  assert(In.OpenSection("eigenstates"));

  int numTwists = In.CountSections("twist");
  In.OpenSection("twist", 0);
  int numBands = In.CountSections("band");
  In.CloseSection(); // "twist"

  perr << "numTwists = " << numTwists << endl;
  perr << "numBands  = " << numBands  << endl;

  Comm.Broadcast (0, numTwists);
  Comm.Broadcast (0, numBands);

  ClusterFactor.resize(numTwists, numBands);
  ClusterFactor = 0.0;
  RMSError.resize(numTwists, numBands);
  RMSError = 0.0;

  int numTasks = numTwists * numBands;
  int firstTask, lastTask;
  ProcTasks (Comm.MyProc(), numTasks, firstTask, lastTask);
  for (int task=firstTask; task <= lastTask; task++) {
    int twist = task / numBands;
    int band  = task % numBands;
    FindOptimalClusterFactor (twist, band);
  }
  Array<double,2> temp(numTwists, numBands);
  Comm.AllSum (ClusterFactor, temp);
  ClusterFactor = temp;
  Comm.AllSum (RMSError, temp);
  RMSError = temp;

  if (Comm.MyProc() == 0) {
    fprintf (stderr, "Twist      Band      Opt. factor     RMS error\n");
    for (int ti=0; ti<numTwists; ti++)
      for (int bi=0; bi<numBands; bi++)
	fprintf (stderr, " %3d        %3d        %5.3f          %8.4e\n",
		 ti, bi, ClusterFactor(ti,bi), RMSError(ti,bi));
  
    Out.NewSection("eigenstates");
    Out.FlushFile();
    // Now write out new orbitals
    for (int ti=0; ti<numTwists; ti++) {
      cerr << "Writing twist " << ti << endl;
      In.OpenSection ("twist", ti);
      Array<double,1> twist_angle;
      In.ReadVar ("twist_angle", twist_angle);
      Out.NewSection ("twist");
      Out.WriteVar ("twist_angle", twist_angle);
      for (int bi=0; bi<numBands; bi++) {
      cerr << "Writing band " << bi << endl;
	In.OpenSection ("band", bi);
	Out.NewSection ("band");
	Orb.Read(In);
	int nx, ny, nz, Nx, Ny, Nz;
	Orb.GetDims (Nx, Ny, Nz);
	nx = (int)round(SizeFactor*(double)Nx);
	ny = (int)round(SizeFactor*(double)Nz);
	nz = (int)round(SizeFactor*(double)Ny);
	
	Orb.Respline (nx, ny, nz, ClusterFactor(ti, bi));
	Orb.Write (Out);
	In.CloseSection (); // band
	Out.CloseSection(); // band
	Out.FlushFile();
      }
      In.CloseSection (); // twist
      Out.CloseSection(); // twist
    }
    Out.CloseFile();
  }
}

void
ResplineClass::Run(int orbNum)
{
  const int numPoints = 1000000;

  CopyHeader();
  assert(In.OpenSection("eigenstates"));

  int numTwists = In.CountSections("twist");
  In.OpenSection("twist", 0);
  int numBands = In.CountSections("bands");
  In.CloseSection(); // "twist"

  Comm.Broadcast (0, numTwists);
  Comm.Broadcast (0, numBands);

  ClusterFactor.resize(numTwists, numBands);
  ClusterFactor = 0.0;
  RMSError.resize(numTwists, numBands);
  RMSError = 0.0;

  for (int ti=0; ti<numTwists; ti++) {
    In.OpenSection("twist", ti);
    int bi = orbNum;
    //for (int bi=0; bi<numBands; bi++) {
    In.OpenSection("band", bi);
    Orb.Read (In);
    int nx, ny, nz;
    Orb.GetDims (nx, ny, nz);
    nx = (int)round(SizeFactor*(double)nx);
    ny = (int)round(SizeFactor*(double)nz);
    nz = (int)round(SizeFactor*(double)ny);
    
    Orb.Respline (nx, ny, nz);
    double uniformError = Orb.UniformError(numPoints);
    fprintf (stdout, "# Cluster factor     Error\n");
    
    for (double factor=1.0001; factor<=15.0001; factor+=0.1) {
      Orb.Respline (nx, ny, nz, factor);
      double error = Orb.NonuniformError(numPoints);
      fprintf (stdout, " %8.3f        %10.12f\n", factor, error);
    }
    In.CloseSection(); // "band";
    In.CloseSection(); // "twist"
  }
  In.CloseSection(); // "eigenstates"
}

void
ResplineClass::RunTest(int orbNum)
{
  const int numPoints = 1000000;
  CopyHeader();
  assert(In.OpenSection("eigenstates"));  
  assert(In.OpenSection("twist", 0));
  assert(In.OpenSection("band", orbNum));
  Orb.Read (In);  
  In.CloseSection(); // "band";
  In.CloseSection(); // "twist"
  In.CloseSection(); // "eigenstates"
  
  fprintf (stdout, "# Nx     uniform error     nonuniform error\n");
  for (int n=40; n<200; n+=5) {
    Orb.Respline(n,n,n);
    double uniformError = Orb.UniformError(numPoints);
    double nonuniformError = uniformError;
    fprintf (stderr, "n = %3d:\n  Factor        error\n", n);
    for (double factor=1.25; factor<=20.0001; factor+=0.25) {
      Orb.Respline (n, n, n, factor);
      double error = Orb.NonuniformError(numPoints);
      nonuniformError = min(error, nonuniformError);
      fprintf (stderr, "   %4.3f     %12.8e\n", 
	       factor, error);
    }
    fprintf (stdout, "  %3d         %12.8e         %12.8e\n", 
	     n, uniformError, nonuniformError);
    fflush(stdout);
  }
}


void
ResplineClass::PlotOrb (int orbNum)
{
  const int numPoints = 1000000;
  assert(In.OpenSection("eigenstates"));  
  assert(In.OpenSection("twist", 0));
  assert(In.OpenSection("band", orbNum));
  Orb.Read (In);  
  In.CloseSection(); // "band";
  In.CloseSection(); // "twist"
  In.CloseSection(); // "eigenstates"

  int n = 45;
  Orb.Respline(n,n,n);
  double minFactor = 1.25;
  double minError = 1.0;
  for (double factor=1.25; factor<=20.0001; factor+=0.25) {
    Orb.Respline (n, n, n, factor);
    double error = Orb.NonuniformError(numPoints);
    fprintf (stderr, " %1.5f  %1.6e\n", factor, error);
    if (error < minError) {
      minFactor = factor;
      minError  = error;
    }
  }
  Orb.Respline (n,n,n,minFactor);
  Orb.Plot ("OrbPlot.dat");
}

int
main(int argc, char **argv)
{
  COMM::Init (argc, argv);

  list<ParamClass> argList;

  argList.push_back(ParamClass("in",     1));
  argList.push_back(ParamClass("out",    1));
  argList.push_back(ParamClass("factor", 1));
  argList.push_back(ParamClass("orbnum", 1));

  CommandLineParserClass parser (argList);
  bool success = parser.Parse (argc, argv);
  success = success && parser.Found("in");
  success = success && parser.Found("out");
  success = success && parser.Found("factor");
  // success = success && parser.Found("orbnum");

  if (!success) {
    cerr << "Usage:\n"
	 << "  respline --in infile.h5 --out outfile.h5 --factor x --orbnum n\n"
	 << "  x <= 1.0 \n";
    exit(-1);
  }
  
  ResplineClass respline;
  
  respline.SetFiles (parser.GetArg("in"),
		     parser.GetArg("out"));
  respline.SetSizeFactor (atof(parser.GetArg("factor").c_str()));
  //respline.RunTest (atol(parser.GetArg("orbnum").c_str()));
  //respline.PlotOrb(atol(parser.GetArg("orbnum").c_str()));
  respline.DoRespline();
  COMM::Finalize();
}
