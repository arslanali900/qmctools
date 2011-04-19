#include "OrbitalSetClass.h"
#include <Common/IO/IO.h>
#include "MuffinTin.h"
#include <einspline/bspline.h>
#include <omp.h>
#include "config.h"
#include "PlaneWaveSum.h"
#include <stdio.h>
#include <stdlib.h>

#ifdef HAVE_CUDA
  #include <cuda_runtime_api.h>
#endif


// class PWOrb_z : public OrbFunctor
// {
// public:
//   zVec &Coefs;
//   GVecsClass &GVecs;
//   complex<double> Prefactor;
//   Vec3 k;
// #ifdef HAVE_CUDA
//   double *GVecs_cuda, *r_cuda;
//   complex<double> *Coefs_cuda, *val_cuda, *grad_cuda,
//     *lapl_cuda;
//   int num_r;
// #endif
// public:
//   complex<double> operator()(Vec3 r)
//   {
//     complex<double> sum;
//     for (int iG=0; iG<GVecs.size(); iG++) {
//       double phase = -dot(GVecs(iG)+k,r);
//       double s,c;
//       sincos(phase, &s, &c);
//       sum += Coefs(iG)*complex<double>(c,s);
//     }
    
//     return Prefactor*sum;
//   }

//   void operator()(Vec3 r, complex<double> &val,
// 		  TinyVector<complex<double>,3> &grad,
// 		  complex<double> &lapl)
//   {
//     val = lapl = complex<double>();
//     grad = TinyVector<complex<double>,3>();
//     for (int iG=0; iG<GVecs.size(); iG++) {
//       double phase = -dot(GVecs(iG)+k,r);
//       double s,c;
//       sincos(phase, &s, &c);
//       complex<double> z(c,s);
//       val +=  Coefs(iG)*z;
//       grad -= Coefs(iG)*complex<double>(-s,c)*(GVecs(iG)+k);
//       lapl -= Coefs(iG)*dot(GVecs(iG)+k,GVecs(iG)+k)*z;
//     }
//     val  *= Prefactor;
//     grad *= Prefactor;
//     lapl *= Prefactor;
//   }
  

// #ifdef HAVE_CUDA
//   void alloc_cuda(int num)
//   {
//     if (val_cuda) {
//       cudaFree (r_cuda);
//       cudaFree (val_cuda);
//       cudaFree (grad_cuda);
//       cudaFree (lapl_cuda);
//     }
//     num_r = num;
//     cudaMalloc((void**)&r_cuda,    3*num*sizeof(double));
//     cudaMalloc((void**)&val_cuda,    num*sizeof(complex<double>));
//     cudaMalloc((void**)&grad_cuda, 3*num*sizeof(complex<double>));
//     cudaMalloc((void**)&lapl_cuda,   num*sizeof(complex<double>));
//   }


//   void operator()(vector<Vec3> &r,
// 		  vector<complex<double> > &val)
//   {
//     if (num_r < r.size())  alloc_cuda(r.size());
//     cudaMemcpy (r_cuda, &r[0][0], 3*r.size()*sizeof(double),
//     		cudaMemcpyHostToDevice);
//     cudaThreadSynchronize();
//     cudaError_t err = cudaGetLastError();
//     if (err != cudaSuccess) {
//       fprintf (stderr, "CUDA error in memcopy:\n  %s\n",
// 	       cudaGetErrorString(err));
//       abort();
//     }
//     plane_wave_sum (r_cuda, GVecs_cuda, k[0], k[1], k[2],
// 		    Coefs_cuda,
//     		    val_cuda, r.size(), GVecs.size());
//     cudaMemcpy (&val[0], val_cuda, val.size()*sizeof(complex<double>),
//     		cudaMemcpyDeviceToHost);
//     for (int i=0; i<val.size(); i++) 
//       val[i] *= Prefactor;

//     // vector<complex<double> > gpu_val = val;
//     // OrbFunctor::operator()(r,val);
//     // for (int i=0; i<val.size(); i++)
//     //   fprintf (stderr, "%12.5e %12.5ei  %12.5e %12.5ei\n",
//     // 	       val[i].real(), val[i].imag(), 
//     // 	       gpu_val[i].real(), gpu_val[i].imag());
  
//   }

//   void operator()(vector<Vec3> &r,
// 		  vector<complex<double> > &val,
// 		  vector<complex<double> > &lapl)
//   {
//     if (num_r != r.size())  alloc_cuda(r.size());
//     cudaMemcpy (r_cuda, &r[0][0], 3*r.size()*sizeof(double),
//     		cudaMemcpyHostToDevice);
//     cudaThreadSynchronize();
//     cudaError_t err = cudaGetLastError();
//     if (err != cudaSuccess) {
//       fprintf (stderr, "CUDA error in memcopy:\n  %s\n",
// 	       cudaGetErrorString(err));
//       abort();
//     }
//     plane_wave_sum (r_cuda, GVecs_cuda, k[0], k[1], k[2],
// 		    Coefs_cuda,
//     		    val_cuda, lapl_cuda, r.size(), GVecs.size());
//     cudaMemcpy (&val[0], val_cuda, val.size()*sizeof(complex<double>),
//     		cudaMemcpyDeviceToHost);
//     cudaMemcpy (&lapl[0], lapl_cuda, val.size()*sizeof(complex<double>),
//     		cudaMemcpyDeviceToHost);
//     for (int i=0; i<val.size(); i++) {
//       val[i]  *= Prefactor;
//       lapl[i] *= Prefactor;
//     }
//   }

// #else
//   void operator()(vector<Vec3> &r,
// 		  vector<complex<double> > &val,
// 		  vector<complex<double> > &lapl)
//   {
//     TinyVector<complex<double>,3> grad;
//     for (int ir=0; ir<r.size(); ir++)
//       (*this)(r[ir], val[ir], grad, lapl[ir]);
//   }
// #endif



//   PWOrb_z (zVec &coefs, GVecsClass &gvecs, 
// 	   Vec3 kpoint, complex<double> pref) :
//     Coefs(coefs), GVecs(gvecs), Prefactor(pref), k(kpoint)
// #ifdef HAVE_CUDA
//     ,val_cuda(NULL), grad_cuda(NULL), lapl_cuda(NULL), r_cuda(NULL),
//     num_r(0)
// #endif
//   {
// #ifdef HAVE_CUDA
//     cudaMalloc((void**)&GVecs_cuda, gvecs.size()*sizeof(Vec3));
//     cudaMalloc((void**)&Coefs_cuda, gvecs.size()*sizeof(complex<double>));
//     cudaMemcpy(GVecs_cuda, &(GVecs(0)[0]), 3*sizeof(double)*gvecs.size(),
// 	       cudaMemcpyHostToDevice);
//     cudaMemcpy(Coefs_cuda, &Coefs(0), 2*sizeof(double)*gvecs.size(),
// 	       cudaMemcpyHostToDevice);
// #endif
//   }

// #ifdef HAVE_CUDA
//   ~PWOrb_z() 
//   {
//     cudaFree (GVecs_cuda);
//     cudaFree (Coefs_cuda);
//     if (r_cuda)       cudaFree (r_cuda);
//     if (val_cuda)     cudaFree (val_cuda);
//     if (grad_cuda)    cudaFree (grad_cuda);
//     if (lapl_cuda)    cudaFree (lapl_cuda);
//   }
// #endif

// //   PWOrb_z (const PWOrb_z &orb) :
// //     Coefs(orb.Coefs), GVecs(orb.GVecs), Prefactor(orb.Prefactor)
// //   {
// // #ifdef HAVE_CUDA
// //     num_r = orb.num_r;
// //     GVecs_cuda = 
// // #endif
// //   }


// };

// class SplineOrb_z : public OrbFunctor
// {
// private:
//   UBspline_3d_z *Spline;
//   LatticeClass &Lattice;
//   Vec3 k;
// public:
//   complex<double> operator()(Vec3 r)
//   {
//     Vec3 u = Lattice.r2u(r);
//     for (int i=0; i<3; i++)
//       u[i] -= floor(u[i]);
//     complex<double> val;
//     eval_UBspline_3d_z (Spline, u[0], u[1], u[2], &val);
//     double s,c, phase;
//     phase = -dot (r, k);
//     sincos(phase, &s,&c);
//     val *= complex<double>(c,s);
//     return val;
//   }

//   void operator()(Vec3 r, complex<double> &val,
// 		  TinyVector<complex<double>,3> &grad,
// 		  complex<double> &lapl)
//   {
//     Vec3 u = Lattice.r2u(r);
//     for (int i=0; i<3; i++)
//       u[i] -= floor(u[i]);
//     TinyVector<complex<double>,3> g;
//     TinyMatrix<complex<double>,3,3> hess;
//     eval_UBspline_3d_z_vgh (Spline, u[0], u[1], u[2], &val,
// 			    &g[0], &hess(0,0));
//     grad = complex<double>();
//     lapl = complex<double>();
//     for (int i=0; i<3; i++)
//       for (int j=0; j<3; j++) {
// 	grad[i] += Lattice.GetB()(i,j) * g[j];
// 	lapl += hess (i,j) * Lattice.GetBBt()(i,j);
//       }
//     double s,c, phase;
//     phase = -dot (r, k);
//     TinyVector<complex<double>,3> ck;
//     ck[0] = k[0];
//     ck[1] = k[1];
//     ck[2] = k[2];
//     sincos(phase, &s,&c);
//     complex<double> emikr(c,s);
//     complex<double> eye(0.0,1.0);
//     complex<double> k_dot_grad = k[0]*grad[0] + k[1]*grad[1] + k[2]*grad[2];
//     lapl = emikr*(lapl - dot(k,k)*val - 2.0*eye*k_dot_grad);
//     grad = emikr*(grad - eye*val*k);
//     val *= emikr;
//   }

//   SplineOrb_z(LatticeClass &lattice, UBspline_3d_z *spline,
// 	      Vec3 kpoint)
//     : Spline(spline), Lattice(lattice), k(kpoint)
//   {
//     // nothing else for now.
//   }
// };
  



void
OrbitalSetClass::WriteESHDFMultiRep (IO::IOSectionClass &out, 
				     int spin, int ki, int band,
				     bool writeComplex)
{
  int nx, ny, nz;
  complex<double> phase_rotation =
    MakeRealSpaceOrbital (spin, ki, band, ShiftOrbitals || Truncate);
  TinyVector<Array<OrbitalClass*,2>,2> &orbs = 
    WriteSuper ? SuperOrbitals : PrimOrbitals;
  vector<TinyVector<Array<OrbitalClass*,2>,2> > &firstOrder = 
    WriteSuper && false ? SuperFirstOrder : PrimFirstOrder;
  CellClass &cell = WriteSuper ? SuperCell : PrimCell;
  FFTBox &FFT = cell.FFT;
  FFTBox &fineFFT = FineCell.FFT;
  FFT.GetDims (nx, ny, nz);


  static FILE *fout = NULL;
  if (!fout) {
    ostringstream fname;
    fname << "multirep_" << nx << ".dat";
    fout = fopen (fname.str().c_str(), "w");
    fprintf (fout, "#  Exact KE      PW KE            3D spline KE     Mixed KE         3D spline RMS    Mixed RMS\n");
  }


  double wfnorm = 0.0;
  double vol = fabs(det(cell.Lattice.GetDirect()));
  double elemvol = vol/(nx*ny*nz);
  double prefact = 1.0; //  double prefact = sqrt(1.0/(vol));
  OrbitalClass &orb = *orbs[spin](ki,band);

  FFT.rBox *= prefact;  
  fineFFT.rBox *= prefact;

  // First, write out the 3D mesh representation
  if (Spline) {
    cout << "Got into section on writing orbitals" << endl;
    if (writeComplex) {
      Array<double,4> eigvec(nx, ny, nz,2);
      for (int ix=0; ix<nx; ix++)
	for (int iy=0; iy<ny; iy++)
	  for (int iz=0; iz<nz; iz++) {
	    eigvec(ix, iy, iz, 0) = FFT.rBox(ix,iy,iz).real();
	    eigvec(ix, iy, iz, 1) = FFT.rBox(ix,iy,iz).imag();
	  }
      // Finally write to file
      out.WriteVar ("psi_r", eigvec);
    }
    else {
      Array<double,3> eigvec(nx, ny, nz);
      for (int ix=0; ix<nx; ix++)
	for (int iy=0; iy<ny; iy++)
	  for (int iz=0; iz<nz; iz++) 
	    eigvec(ix, iy, iz) = FFT.rBox(ix,iy,iz).imag();
      // Finally write to file
      out.WriteVar ("psi_r", eigvec);
    }
  }

  // Now write out the atomic representation
  // Create a bspline object
  
  UBspline_3d_z *spline;
  Ugrid x_grid, y_grid, z_grid;
  BCtype_z xBC, yBC, zBC;
  xBC.lCode = PERIODIC;  xBC.rCode = PERIODIC; 
  yBC.lCode = PERIODIC;  yBC.rCode = PERIODIC; 
  zBC.lCode = PERIODIC;  zBC.rCode = PERIODIC; 
  x_grid.start = 0.0;   x_grid.end = 1.0;  x_grid.num = nx;
  y_grid.start = 0.0;   y_grid.end = 1.0;  y_grid.num = ny;
  z_grid.start = 0.0;   z_grid.end = 1.0;  z_grid.num = nz;
  

  // UBspline_3d_z *finespline, *laplspline;
  // xBC.lCode = PERIODIC;  xBC.rCode = PERIODIC; 
  // yBC.lCode = PERIODIC;  yBC.rCode = PERIODIC; 
  // zBC.lCode = PERIODIC;  zBC.rCode = PERIODIC; 
  // int fnx, fny, fnz;
  // FineCell.FFT.GetDims(fnx, fny, fnz);
  // //  cerr << "Fine cell dimensions are " << fnx << "x" << fny << "x" << fnz << endl;
  // x_grid.start = 0.0;   x_grid.end = 1.0;  x_grid.num = fnx;
  // y_grid.start = 0.0;   y_grid.end = 1.0;  y_grid.num = fny;
  // z_grid.start = 0.0;   z_grid.end = 1.0;  z_grid.num = fnz;
  
  // finespline = create_UBspline_3d_z (x_grid, y_grid, z_grid, 
  // 				     xBC, yBC, zBC, FineCell.FFT.rBox.data());
  // MakeFineLaplacian (spin, ki, band, ShiftOrbitals || Truncate, phase_rotation);
  // fineFFT.rBox *= prefact;
  // laplspline = create_UBspline_3d_z (x_grid, y_grid, z_grid, 
  // 				     xBC, yBC, zBC, FineCell.FFT.rBox.data());

  // SplineOrb_z fineorb(cell.Lattice, finespline);
  // SplineOrb_z laplorb(cell.Lattice, laplspline);

  //PWOrb_z pworb(orb.GetCoefs(), cell.GVecs, orb.Getk(),
  //		phase_rotation*sqrt(1.0/cell.Lattice.GetVolume()));
  PWOrb_z pworb(orb.GetCoefs(), cell.GVecs, orb.Getk(),
		phase_rotation);

  // Now, do angular projections for each atom
  int nat = AtomicOrbitals.size();
  vector<Array<complex<double>,2> > radialFuncs(nat);
  {
    for (int iat=0; iat<nat; iat++) {
     AtomicOrbital &atorb = AtomicOrbitals[iat];
     //atorb.Project (pworb, sporb, radialFuncs[iat]);
     //atorb.ProjectAnalytic (pworb, sporb, radialFuncs[iat]);
   }
  }
  //cerr << "Done projection.\n";
  //complex<double> prefactor = 
  //  phase_rotation*sqrt(1.0/cell.Lattice.GetVolume());
  complex<double> prefactor = phase_rotation;

  for (int iat=0; iat<AtomicOrbitals.size(); iat++) {
    AtomicOrbital &atorb = AtomicOrbitals[iat];
 
    AtomicDataClass &data = orb.AtomicData[iat];
    Array<double,3> splineDataReal(data.SplineData.extent(0),
				   data.SplineData.extent(1), 2);
    for (int i=0; i<data.SplineData.extent(0); i++)
      for (int j=0; j<data.SplineData.extent(1); j++) {
     	splineDataReal(i,j,0) = (prefactor*data.SplineData(i,j)).real();
     	splineDataReal(i,j,1) = (prefactor*data.SplineData(i,j)).imag();
      }

    // for (int i=0; i<radialFuncs[iat].extent(0); i++)
    //   for (int j=0; j<radialFuncs[iat].extent(1); j++) {
    // 	splineDataReal(i,j,0) = radialFuncs[iat](i,j).real();
    // 	splineDataReal(i,j,1) = radialFuncs[iat](i,j).imag();
    //   }

    ostringstream atom_name;
    atom_name << "radial_spline_" << iat;
    out.WriteVar(atom_name.str().c_str(), splineDataReal);

    Array<double,3> polyCoefsReal(data.PolyCoefs.extent(0),
				  data.PolyCoefs.extent(1), 2);
    for (int i=0; i<data.PolyCoefs.extent(0); i++)
      for (int j=0; j<data.PolyCoefs.extent(1); j++) {
	polyCoefsReal(i,j,0) = (prefactor*data.PolyCoefs(i,j)).real();
	polyCoefsReal(i,j,1) = (prefactor*data.PolyCoefs(i,j)).imag();
      }
    ostringstream poly_name;
    poly_name << "poly_coefs_" << iat;
    out.WriteVar(poly_name.str().c_str(), polyCoefsReal); 
    out.FlushFile();
  }
  
  if (CompareHybrid) {
    spline = create_UBspline_3d_z (x_grid, y_grid, z_grid, 
				 xBC, yBC, zBC, FFT.rBox.data());

    SplineOrb_z sporb(cell.Lattice, spline, orb.Getk());

    // Compute kinetic energy
    double fine_norm=0.0, fine_ke=0.0;
    double mixed_norm=0.0,  mixed_ke=0.0;
    double ex_norm = 0.0, ex_ke=0.0;
    double pw_norm = 0.0, pw_ke=0.0;
    double mixed_dev2 = 0.0, fine_dev2 = 0.0;
    
    for (int iG=0; iG<cell.GVecs.size(); iG++) {
      double nrm = norm (orb.GetCoefs()(iG));
      ex_ke += 0.5*dot (cell.GVecs(iG),cell.GVecs(iG)) *nrm;
      ex_norm += nrm;
    }
    

    Vec3 u;
    int numInside=0;
    double delta=0.0411931;
    int numdev = 0;
    
    //#pragma omp reduction(+,fine_norm, fine_ke, pw_nrom, pw_ke, mixed_norm, mixed_ke, mixed_dev2, fine_dev2, num_dev)
    vector<Vec3> rvecs;
    for (u[0]=0.5*delta; u[0]<1.0; u[0]+=delta)
      for (u[1]=0.5*delta; u[1]<1.0; u[1]+=delta)
	for (u[2]=0.5*delta; u[2]<1.0; u[2]+=delta) 
	  rvecs.push_back(cell.Lattice.u2r(u));
    
    int numr = rvecs.size();
    vector<complex<double> > pw_vals(numr), pw_lapl(numr);
    pworb (rvecs, pw_vals, pw_lapl);
    
    
    for (int ir=0; ir<numr; ir++) {
      Vec3 r = rvecs[ir];
      complex<double> fine_val, fine_lapl, mixed_val, mixed_lapl;
      TinyVector<complex<double>,3> fine_grad, mixed_grad;
      sporb(r, fine_val, fine_grad, fine_lapl);
      //pworb(r,   pw_val,   pw_grad,   pw_lapl);	
      // fineorb(r,   pw_val,   pw_grad,   pw_lapl);
      // pw_lapl = laplorb(r);
      fine_norm += norm(fine_val);
      fine_ke   += real(-0.5*conj(fine_val)*fine_lapl);
      pw_norm   += norm (pw_vals[ir]);
      pw_ke     += real(-0.5*conj(pw_vals[ir])*pw_lapl[ir]);
      
      bool inside=false;
      for (int iat=0; iat<AtomicOrbitals.size(); iat++) {
	AtomicOrbital &atorb = AtomicOrbitals[iat];
	
	Vec3 dr = r - atorb.get_pos();
	Vec3 du = cell.Lattice.r2u (dr);
	for (int i=0; i<3; i++) du[i] -= round(du[i]);
	dr = cell.Lattice.u2r(du);
	Vec3 rnew = atorb.get_pos() + dr;
	double dist = sqrt(dot(dr,dr));
	if ((dist < atorb.get_radius())) {
	  AtomicOrbitals[iat].eval(r, cell.Lattice,
				   mixed_val, mixed_grad, mixed_lapl);
	  inside=true;
	  numInside++;
	  break;
	}
	//if (std::isnan(mixed_lapl.real()))
	if (isnan(mixed_lapl.real()))
	  cerr << "NAN at dr =  " << sqrt(dot(dr,dr)) << endl;
      }
      double nl = imag(-0.5*conj(mixed_val)*mixed_lapl);
      if (!inside) {
	mixed_val = fine_val;
	mixed_grad = fine_grad;
	mixed_lapl = fine_lapl;
      }
      nl = real(-0.5*conj(mixed_val)*mixed_lapl);
      
      mixed_norm += norm(mixed_val);
      mixed_ke   += real(-0.5*conj(mixed_val)*mixed_lapl);
      
      double dev = real(-0.5*conj(pw_vals[ir])*pw_lapl[ir]) - real(-0.5*conj(mixed_val)*mixed_lapl);
      mixed_dev2 += dev*dev;
      dev = real(-0.5*conj(pw_vals[ir])*pw_lapl[ir]) - real(-0.5*conj(fine_val)*fine_lapl);
      fine_dev2  += dev*dev;
      numdev++;
    }
    
    //  cerr << "numInside = " << numInside << endl;
    //if (std::isnan(mixed_ke))
    if (isnan(mixed_ke))
      cerr << "mixed_ke = " << mixed_ke << endl;
    fprintf (stderr, "PW KE    = %1.8f\n", pw_ke/pw_norm);
    fprintf (stderr, "Exact KE = %1.8f\n", ex_ke/ex_norm);
    fprintf (stderr, "Fine  KE = %1.8f  numer=%1.8f  denom=%1.8f\n", 
	     fine_ke/fine_norm, fine_ke, fine_norm);
    fprintf (stderr, "Mixed KE = %1.8f  numer=%1.8f  denom=%1.8f\n", 
	     mixed_ke/mixed_norm, mixed_ke, mixed_norm);
    fprintf (stderr, "Fine  RMS KE deviation = %1.8f\n", sqrt(fine_dev2/numdev));
    fprintf (stderr, "Mixed RMS KE deviation = %1.8f\n", sqrt(mixed_dev2/numdev));
    
    fprintf (fout, "%16.10e %16.10e %16.10e %16.10e %16.10e %16.10e\n",
	     ex_ke/ex_norm, pw_ke/pw_norm, fine_ke/fine_norm, mixed_ke/mixed_norm,
	     sqrt(fine_dev2/numdev), sqrt(mixed_dev2/numdev));
    fflush (fout);
    
    destroy_Bspline (spline);
  }
  // destroy_Bspline (finespline);
  // destroy_Bspline (laplspline);
}

bool
OrbitalSetClass::ReadMultiRep (string fname)
{
  CellClass &cell = WriteSuper ? SuperCell : PrimCell;

  IO::IOSectionClass in;
  assert (in.OpenFile (fname));

  int numSpecies = in.CountSections("Species");
  for (int atom=0; atom < cell.IonPos.size(); atom++) {
    for (int isp=0; isp<numSpecies; isp++) {
      assert(in.OpenSection("Species", isp));
      string symbol;
      double radius;
      assert (in.ReadVar ("radius", radius));
      int lMax, spline_points;
      assert (in.ReadVar("spline_points", spline_points));
      assert (in.ReadVar("lmax", lMax));
      assert (in.ReadVar("name", symbol));
      if (SymbolToZMap.find(symbol) == SymbolToZMap.end()) {
	cerr << "Unrecognized element name " << symbol << ".\n";
	abort();
      }
      int Z = SymbolToZMap[symbol];
      if (cell.AtomTypes(atom) == Z) {
	AtomicOrbital orb;
	orb.Set (lMax, radius, spline_points, cell.IonPos(atom));
	AtomicOrbitals.push_back(orb);
	fprintf (stderr, "Adding AtomicOrbital for type %3s at position  %7.3f %7.3f %7.3f  with lmax=%d\n",
		 symbol.c_str(), cell.IonPos(atom)[0], 
		 cell.IonPos(atom)[1], cell.IonPos(atom)[2], lMax);
	// cerr << "Adding AtomicOrbital for type \"" << symbol 
	//      << "\" at position " << cell.IonPos(atom) << endl;
      }
      in.CloseSection(); // "Species"
    }
  }
  
  AtomicOrbitals[0].TimeYlm();
  
  return true;
}


void
OrbitalSetClass::WriteESHDFSpline (IO::IOSectionClass &out, 
				  int spin, int ki, int band,
				  bool writeComplex)
{
  int nx, ny, nz;
  complex<double> phase_rotation =
    MakeRealSpaceOrbital (spin, ki, band, ShiftOrbitals || Truncate);
  TinyVector<Array<OrbitalClass*,2>,2> &orbs = 
    WriteSuper ? SuperOrbitals : PrimOrbitals;
  vector<TinyVector<Array<OrbitalClass*,2>,2> > &firstOrder = 
    WriteSuper && false ? SuperFirstOrder : PrimFirstOrder;
  CellClass &cell = WriteSuper ? SuperCell : PrimCell;
  FFTBox &FFT = cell.FFT;
  FFT.GetDims (nx, ny, nz);

  double wfnorm = 0.0;
  double vol = fabs(det(cell.Lattice.GetDirect()));
  double elemvol = vol/(nx*ny*nz);
  double prefact = sqrt(1.0/(vol));
  OrbitalClass &orb = *orbs[spin](ki,band);

  double eigval = 0.0;
  if (band < NewOrbCoefs[spin].extent(1))
    for (int bi=0; bi<NewOrbCoefs[spin].extent(1); bi++)
      eigval += orbs[spin](ki,bi)->GetEigVal() * 
	norm(NewOrbCoefs[spin](ki,band,bi));
  else
    eigval = orbs[spin](ki,band)->GetEigVal();

  if (Truncate && (&(orb.GetCenter()) != NULL)) {
    // First, find the grid points which must be included
    Vec3 u, du;
    du[0]=1.0/(double)nx; du[1]=1.0/(double)ny; du[2]=1.0/(double)nz;
    double rad  = orb.GetCenter().Radius + orb.GetCenter().SkinThickness;
    double rad2 = rad*rad;
    int minx=nx, maxx=0, miny=ny, maxy=0, minz=nz, maxz=0;
    for (int ix=0; ix<nx; ix++) {
      u[0] = (double)ix*du[0] - 0.5;
      for (int iy=0; iy<ny; iy++) {
	u[1] = (double)iy*du[1] -0.5; 
	for (int iz=0; iz<nz; iz++) {
	  u[2] = (double)iz*du[2] -0.5;
	  Vec3 r = cell.Lattice.u2r(u);
	  if (dot(r,r) <= rad2) {
	    minx = min(minx, ix-4); maxx = max (maxx, ix+4);
	    miny = min(miny, iy-4); maxy = max (maxy, iy+4);
	    minz = min(minz, iz-4); maxz = max (maxz, iz+4);
	  }
	}
      }
    }
    if (writeComplex) {
      Array<double,4> eigvec(maxx-minx+1, maxy-miny+1, maxz-minz+1, 2);
      for (int ix=minx; ix<=maxx; ix++)
	for (int iy=miny; iy<=maxy; iy++)
	  for (int iz=minz; iz<=maxz; iz++) {
	    int jx = (ix+nx)%nx;
	    int jy = (iy+ny)%ny;
	    int jz = (iz+nz)%nz;
	    eigvec(ix-minx, iy-miny, iz-minz,0) = 
	      prefact*FFT.rBox(jx,jy,jz).real();
	    eigvec(ix-minx, iy-miny, iz-minz,1) = 
	      prefact*FFT.rBox(jx,jy,jz).imag();
	  }
      out.WriteVar ("psi_r", eigvec);
    }
    else {
      Array<double,3> eigvec(maxx-minx+1, maxy-miny+1, maxz-minz+1);
      for (int ix=minx; ix<=maxx; ix++)
	for (int iy=miny; iy<=maxy; iy++)
	  for (int iz=minz; iz<=maxz; iz++) {
	    int jx = (ix+nx)%nx;
	    int jy = (iy+ny)%ny;
	    int jz = (iz+nz)%nz;
	    eigvec(ix-minx, iy-miny, iz-minz) = 
	      prefact*FFT.rBox(jx,jy,jz).real();
	  }
      out.WriteVar ("psi_r", eigvec);
    }
    Array<double,1> uMin(3), uMax(3);
    uMin(0) = (double)minx*du[0];    uMax(0) = (double)maxx*du[0];
    uMin(1) = (double)miny*du[1];    uMax(1) = (double)maxy*du[1];
    uMin(2) = (double)minz*du[2];    uMax(2) = (double)maxz*du[2];
    out.WriteVar ("umin", uMin);
    out.WriteVar ("umax", uMax);
    out.WriteVar ("truncated", true);
  }
  else {  // This is not a localized orbital
    if (writeComplex) {
      Array<double,4> eigvec(nx, ny, nz,2);
      for (int ix=0; ix<nx; ix++)
	for (int iy=0; iy<ny; iy++)
	  for (int iz=0; iz<nz; iz++) {
	    eigvec(ix, iy, iz, 0) = prefact*FFT.rBox(ix,iy,iz).real();
	    eigvec(ix, iy, iz, 1) = prefact*FFT.rBox(ix,iy,iz).imag();
	  }
      // Finally write to file
      out.WriteVar ("psi_r", eigvec);
      
      for (int ideriv=0; ideriv<firstOrder.size(); ideriv++) {
	FFTBox &fft = PrimCell.FFT;
	MakeFirstOrderOrbital (ideriv, spin, ki, band, 
			       ShiftOrbitals || Truncate, 
			       phase_rotation);
	ostringstream name;
	int atom, dir;
	firstOrder[ideriv][spin](ki,band)->GetDeriv(atom,dir);
	
	name << "dpsi_" << atom << "_" << dir << "_r";
	for (int ix=0; ix<nx; ix++)
	  for (int iy=0; iy<ny; iy++)
	    for (int iz=0; iz<nz; iz++) {
	      eigvec(ix, iy, iz, 0) = prefact*fft.rBox(ix,iy,iz).real();
	      eigvec(ix, iy, iz, 1) = prefact*fft.rBox(ix,iy,iz).imag();
	    }
	out.WriteVar (name.str(), eigvec);
      }
    }
    else {
      Array<double,3> eigvec(nx, ny, nz);
      for (int ix=0; ix<nx; ix++)
	for (int iy=0; iy<ny; iy++)
	  for (int iz=0; iz<nz; iz++) 
	    eigvec(ix, iy, iz) = prefact*FFT.rBox(ix,iy,iz).imag();
      // Finally write to file
      out.WriteVar ("psi_r", eigvec);
      
      for (int ideriv=0; ideriv<firstOrder.size(); ideriv++) {
	MakeFirstOrderOrbital (ideriv, spin, ki, band, 
			       ShiftOrbitals || Truncate, 
			       phase_rotation);
	FFTBox &fft = PrimCell.FFT;
	ostringstream name;
	int atom, dir;
	firstOrder[ideriv][spin](ki,band)->GetDeriv(atom,dir);
	
	name << "dpsi_" << atom << "_" << dir << "_r";
	for (int ix=0; ix<nx; ix++)
	  for (int iy=0; iy<ny; iy++)
	    for (int iz=0; iz<nz; iz++) 
	      eigvec(ix, iy, iz) = prefact*fft.rBox(ix,iy,iz).real();
	out.WriteVar (name.str(), eigvec);
      }
    }
    
    // out.WriteVar ("eigenvalue", eigval);
    if (&(orb.GetCenter()) != NULL) {
      Vec3 rv = orb.GetCenter().r;
      Array<double,1> r(3);
      r(0) = rv[0];  r(1)=rv[1]; r(2)=rv[2];
      out.WriteVar ("center", r);
      int N = orb.GetCenter().IdenticalSites.size();
      Array<double,2> centers(N,3);
      for (int i=0; i<N; i++)
	for (int j=0; j<3; j++)
	  centers(i,j) = orb.GetCenter().IdenticalSites(i)[j];
      out.WriteVar("centers", centers);
      if (orb.GetCenter().Reflections.size() > 0) {
	int N = orb.GetCenter().Reflections.size();
	Array<double,2> reflections(N,3);
	for (int i=0; i<N; i++)
	  for (int j=0; j<3; j++)
	    reflections(i,j) = orb.GetCenter().Reflections(i)[j];
	out.WriteVar ("reflections", reflections);
      }
      out.WriteVar ("radius", orb.GetCenter().Radius
		    + orb.GetCenter().SkinThickness);
    }
    
    /////////////////////
    // Muffin-tin part //
    /////////////////////
    if (APW.lMax > 0) {
      for (int atom=0; atom<APW.NumAtoms; atom++) {
	out.NewSection("muffin_tin");
	Array<double,1> r(APWRadialPoints);
	double radius=APW.u(atom,0)[0].grid.End();
	double dr = radius/(double)(APWRadialPoints-1);
	for (int ir=0; ir<APWRadialPoints; ir++)
	  r(ir) = (double)ir*dr;
	int num_lm = (APW.lMax+1)*(APW.lMax+1);
	Array<complex<double>,2> u_lm_r(num_lm,APWRadialPoints);
	Array<complex<double>,1> du_lm_dr(num_lm);
	APW.evaluate      (spin, atom, ki, band, u_lm_r, du_lm_dr);
	LocalOrbitals.add (spin, atom, ki, band, u_lm_r);
	
	// Makes sure we use the same phase rotation as we did for the
	// B-spline part
	for (int lm=0; lm<num_lm; lm++) {
	  du_lm_dr (lm) *= phase_rotation;
	  for (int ir=0; ir<u_lm_r.extent(1); ir++)
	    u_lm_r(lm,ir) *= phase_rotation;
	}
	
	// Now make a 3D real array out of the 2D complex one
	Array<double,3> u_lm_real  (u_lm_r.extent(0), u_lm_r.extent(1), 2);
	Array<double,2> du_lm_real (u_lm_r.extent(0), 2);
	for (int lm=0; lm<u_lm_r.extent(0); lm++) {
	  for (int ir=0; ir<u_lm_r.extent(1); ir++) {
	    u_lm_real(lm, ir, 0) = u_lm_r(lm,ir).real();
	    u_lm_real(lm, ir, 1) = u_lm_r(lm,ir).imag();
	  }
	  du_lm_real(lm, 0) = du_lm_dr(lm).real();
	  du_lm_real(lm, 1) = du_lm_dr(lm).imag();
	}
	
	//out.WriteVar("r", r);
	out.WriteVar("r", APW.r_data(atom));
	out.WriteVar("u_lm_r", u_lm_real);
	out.WriteVar("du_lm_dr", du_lm_real);
	
	out.CloseSection(); // "muffin_tin"
      }
    }
    out.FlushFile();
  }
}



bool
OrbitalSetClass::Write_ESHDF (string fname, int lappwr)
{
  if (Comm.MyProc() != 0)
    return true;
  cout<<" Writing laplacians to level "<<lappwr<<endl;
  
  CellClass &cell = WriteSuper ? SuperCell : PrimCell;
  DensityClass &density = WriteSuper ? SuperDensity : PrimDensity;
  DensityClass &vhxc    = WriteSuper ? SuperVHXC : PrimVHXC;
  GVecsClass &GVecs = cell.GVecs;
  TinyVector<Array<OrbitalClass*,2>,2> &orbs =
    WriteSuper ? SuperOrbitals : PrimOrbitals;

  int numSpins = 1 + (orbs[1].size() > 0);

  if (UseMultiRep) {
    int num_k = orbs[0].extent(0);
    int num_bands = orbs[0].extent(1);
    // Resize the atomic data for all orbitals
    for (int ispin=0; ispin<numSpins; ispin++)
      for (int ik=0; ik<num_k; ik++)
	for (int ib=0; ib<num_bands; ib++)
	  orbs[ispin](ik,ib)->AtomicData.resize(AtomicOrbitals.size());
//#pragma omp parallel 
    {
//#pragma omp for schedule(dynamic,1)
#ifdef HAVE_CUDA
      double positions[3*AtomicOrbitals.size()];
      for (int iat = 0; iat<AtomicOrbitals.size(); iat++) {
	//cerr << "  Projecting all orbitals for atom " << iat << endl;
	//AtomicOrbitals[iat].ProjectAnalyticCudaDebug(iat, orbs);
	
	positions[3*iat] = (AtomicOrbitals[iat].get_pos())[0];
	positions[3*iat+1] = (AtomicOrbitals[iat].get_pos())[1];
	positions[3*iat+2] = (AtomicOrbitals[iat].get_pos())[2];
      }
      ProjectAnalyticCuda(AtomicOrbitals.size(), AtomicOrbitals[0].get_lmax(),
      			  AtomicOrbitals[0].get_splinePoints(), AtomicOrbitals[0].get_polyOrder(),
      			  AtomicOrbitals[0].get_polyRadius(), AtomicOrbitals[0].get_outerRadius(),
      			  positions, orbs);
			  
#else
      for (int iat=0; iat<AtomicOrbitals.size(); iat++) {
	cerr << "  Projecting all orbitals for atom " << iat << endl;
	AtomicOrbitals[iat].ProjectAnalytic(iat, orbs);
      }
#endif
    }
  }


  IO::IOSectionClass out;
  if (!out.NewFile (fname)) return false;
  Array<int,1> version(3);
  version(0) = 2;
  version(1) = 1;
  version(2) = 0;
  out.WriteVar("version", version);
  out.WriteVar("format", "ES-HDF");
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
      lattice(i,j) = cell.Lattice.GetDirect()(i,j);
  out.WriteVar("primitive_vectors", lattice);
  out.CloseSection(); // "supercell"

  out.NewSection("atoms");
  out.WriteVar("number_of_atoms", cell.IonPos.size());

  // Compute number of distinct species
  vector<int> species, zion;
  Array<int,1> atom_species(cell.AtomTypes.size());
  for (int iat=0; iat<cell.AtomTypes.size(); iat++) {
    int atType = cell.AtomTypes(iat);
    bool found = false;
    for (int i=0; i<species.size(); i++)
      if (species[i] == atType) {
	atom_species(iat) = i;
	found = true;
      }
    if (!found) {
      atom_species(iat) = species.size();
      species.push_back(atType);
      if (cell.Zion.size())
	zion.push_back (cell.Zion(iat));
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

  out.WriteVar("positions", cell.IonPos);
  Array<Vec3,1> reduced_positions(cell.IonPos.size());
  for (int i=0; i<cell.IonPos.size(); i++)
    reduced_positions(i) = cell.Lattice.r2u(cell.IonPos(i));
  out.WriteVar("reduced_positions", reduced_positions);

  if (cell.IonForces.size())
    out.WriteVar("forces", cell.IonForces);

  out.CloseSection(); // "atoms"
  
  out.NewSection("electrons");
  out.WriteVar ("have_dpsi", (int)(PrimFirstOrder.size() > 0));
  
  PrimDensity.WriteESHDF(out);
  PrimVHXC.WriteESHDF (out, "VHXC");
  TinyVector<int,2> numElecs = NumElectrons;
  if (WriteSuper)
    numElecs *= (int)round(fabs(det(TileMatrix)));
  out.WriteVar ("number_of_electrons", numElecs);
  int numk = orbs[0].extent(0);
  out.WriteVar ("number_of_kpoints", numk);
  if (Functional == "")
    Functional = "unknown";
  out.WriteVar ("functional", Functional);
  out.WriteVar ("total_energy", TotalEnergy);
  out.WriteVar ("number_of_spins", numSpins);
  out.WriteVar ("psi_r_is_complex", (int)(!Real));
  out.WriteVar ("number_of_atomic_orbitals",
		(int)AtomicOrbitals.size());

  out.SetUnderscores(true);
  for (int iat=0; iat<AtomicOrbitals.size(); iat++) {
    out.NewSection("atomic_orbital");
    AtomicOrbitals[iat].Write_ESHDF (out);
    out.CloseSection();
  }


  for (int ik=0; ik<orbs[0].extent(0); ik++) {
    out.NewSection("kpoint");
    Vec3 reduced_k = -1.0*orbs[0](ik,0)->GetTwist();
    out.WriteVar("reduced_k", reduced_k);
    out.WriteVar("weight", 1.0);
    int numG = orbs[0](ik,0)->GetCoefs().size();
    out.WriteVar("number_of_gvectors", numG);
    Array<double,2> cG(numG, 2);
    Array<Int3,1> gints(cell.GVecs.GetMultipliers().size());
    // We inverted the signs of the g-vectors when reading.
    // We need to do it again when writing to be correct.
    for (int i=0; i<gints.size(); i++)
      for (int j=0; j<3; j++)
      gints(i)[j] = -cell.GVecs.GetMultipliers()(i)[j];
    out.WriteVar("gvectors", gints);
    for (int spin=0; spin<numSpins; spin++) {
      out.NewSection("spin");
      int numStates = orbs[spin].extent(1);
      int numOrbitals = numStates*(1+lappwr);
      out.WriteVar("number_of_states", numOrbitals);
      out.WriteVar("number_of_lapl", lappwr);
      Array<double,1> eigvals(numOrbitals), occ(numOrbitals);
      for (int ist=0; ist<numStates; ist++) {
      eigvals(ist*(1+lappwr)) = orbs[spin](ik,ist)->GetEigVal();
      occ(ist*(1+lappwr)) = (numSpins > 1) ? 1.0 : 2.0;
      out.NewSection("state");
      for (int iG=0; iG<numG; iG++) {
        cG(iG,0) = orbs[spin](ik,ist)->GetCoefs()(iG).real();
        cG(iG,1) = orbs[spin](ik,ist)->GetCoefs()(iG).imag();
      }
	//out.WriteVar ("psi_g", orbs[spin](ik,ist)->GetCoefs());
	out.WriteVar ("psi_g", cG);
	if (UseMultiRep)
	  WriteESHDFMultiRep (out, spin, ik, ist, !Real);
	else if (Spline)
	  WriteESHDFSpline (out, spin, ik, ist, !Real);
   out.CloseSection();// "state"
   
   std::vector<double> powers(lappwr);
   powers[0]=-1.0;
   powers[1]=-0.5;
   for(int iL=2; iL<lappwr; iL++) powers[iL]=0.5*(iL-1);
   
   for (int iL=0; iL<lappwr; iL++) 
   {
     out.NewSection("state");
     eigvals(ist*(1+lappwr)+iL+1) = eigvals(ist*(1+lappwr))+ (iL+1)*1000;
     occ(ist*(1+lappwr)+iL+1) = occ(ist*(1+lappwr));
     for (int iG=0; iG<numG; iG++) 
     {
        double pf(0);
        for (int j=0; j<3; j++) pf += (reduced_k[j]+gints(iG)[j])*(reduced_k[j]+gints(iG)[j]);
        if (powers[iL]>0) 
          pf=pow(pf,powers[iL]);
        else
          pf=pow(pf+0.1,powers[iL]);
        cG(iG,0) = pf*orbs[spin](ik,ist)->GetCoefs()(iG).real();
        cG(iG,1) = pf*orbs[spin](ik,ist)->GetCoefs()(iG).imag();
      }
      out.WriteVar ("psi_g", cG);
      if (UseMultiRep)
        WriteESHDFMultiRep (out, spin, ik, ist, !Real);
      else if (Spline)
        WriteESHDFSpline (out, spin, ik, ist, !Real);
      out.CloseSection();// "state"
   }
   
      }
      out.WriteVar("eigenvalues", eigvals);
      out.WriteVar("occupations", occ);

      out.CloseSection(); // "spin"
    }
    out.CloseSection(); // "kpoint";
  }
  out.SetUnderscores(false);
  FFTBox &FFT = cell.FFT;
  Array<int,1> mesh(3);
  FFT.GetDims (mesh(0), mesh(1), mesh(2));
  if (Spline)
    out.WriteVar("psi_r_mesh", mesh);
  out.CloseSection(); // "electrons"


  // //  out.WriteVar ("num_bands",      orbs[0].extent(1)+orbs[1].extent(1));
  // out.WriteVar ("num_up_bands",   orbs[0].extent(1));
  // out.WriteVar ("num_down_bands", orbs[1].extent(1));
  // // Write out the number of core states
  // int numCore = CoreStates.size();
  // out.WriteVar ("num_core_states", numCore);
  // out.WriteVar ("num_twists",     orbs[0].extent(0));
  // int num_spins = SpinPolarized ? 2 : 1;
  // out.WriteVar ("num_spins", num_spins);
  // out.WriteVar ("maximum_ecut", ECut);

  // // Write muffin-tin information
  // if (APW.lMax > 0) {
  //   out.NewSection("muffin_tins");
  //   out.WriteVar("num_tins", APW.NumAtoms);
  //   for (int atom=0; atom<APW.NumAtoms; atom++) {
  //     out.NewSection("muffin_tin");
  //     out.WriteVar("radius", APW.u(atom,0)[0].grid.End());
  //     out.WriteVar("r", APW.r_data(atom));
  //     Array<double,1> center(3);
  //     for (int i=0; i<3; i++)
  // 	center(i) = cell.IonPos(atom)[i];
  //     out.WriteVar("center", center);
  //     out.WriteVar("lmax", APW.lMax);
  //     out.WriteVar("num_radial_points", APWRadialPoints);
  //     out.CloseSection(); // "muffin_tin"
  //   }

  //   out.CloseSection(); // "muffin_tins"
  // }

  // out.CloseSection(); // "parameters"

  // out.NewSection("basis");
  // if (Spline)
  //   out.WriteVar ("type", "spline");
  // else {
  //   out.WriteVar ("type", "planewaves");
  //   out.WriteVar ("num_planewaves", GVecs.size());
  //   Array<double,2> gvecs (GVecs.size(), 3);
  //   Array<int,2> multipliers (GVecs.size(), 3);
  //   for (int i=0; i<GVecs.size(); i++)
  //     for (int j=0; j<3; j++) {
  // 	gvecs(i,j) = GVecs(i)[j];
  // 	multipliers(i,j) = GVecs.Multiplier(i)[j];
  //     }
  //   out.WriteVar ("planewaves", gvecs);
  //   out.WriteVar ("multipliers", multipliers);
  // }
  // out.CloseSection(); // "basis"


  // out.NewSection("ions");
  // Array<double,2> pos(cell.IonPos.size(),3);
  // for (int i=0; i<cell.IonPos.size(); i++) {
  //   pos(i,0) = cell.IonPos(i)[0];
  //   pos(i,1) = cell.IonPos(i)[1];
  //   pos(i,2) = cell.IonPos(i)[2];
  // }
  // out.WriteVar ("pos", pos);
  // out.WriteVar ("atom_types", cell.AtomTypes);
  // out.CloseSection(); // "ions"
  

  // out.NewSection ("eigenstates");
  // int numBandsWritten = 0;
  // for (int icore=0; icore<CoreStates.size(); icore++)
  //   CoreStates[icore].SetupSplines();

  // for (int ik=0; ik<orbs[0].extent(0); ik++) {
  //   out.NewSection("twist");
  //   Vec3 twist = orbs[0](ik,0)->GetTwist();
    
  //   Array<double,1> twist_angle(3);
  //   twist_angle(0)=twist[0]; 
  //   twist_angle(1)=twist[1]; 
  //   twist_angle(2)=twist[2];
  //   out.WriteVar ("twist_angle", twist_angle);
  //   int numLoc=0;
  //   // Write core states, if we have them
  //   for (int icore=0; icore<CoreStates.size(); icore++) {
  //     out.NewSection("core_state");
  //     double rmax = CoreStates[icore].r(CoreStates[icore].r.size()-1);
  //     out.WriteVar("rmax", rmax);
  //     Array<double,1> fArray(CoreRadialPoints), gArray(CoreRadialPoints);
  //     double dr = rmax / (double)(CoreRadialPoints-1);
  //     // Find ir for around r=0.1
  //     int icut = (int)ceil(0.05/dr);
  //     double rcut = dr * icut;
  //     double gcut = CoreStates[icore].g(rcut) / rcut;
  //     for (int ir=1; ir<CoreRadialPoints; ir++) {
  // 	double rval = dr * (double)ir;
  // 	double gval = CoreStates[icore].g(rval) / rval;
  // 	double fval = CoreStates[icore].f(rval) / rval;
  // 	gArray(ir) = CoreStates[icore].g(rval) / rval;
  // 	fArray(ir) = CoreStates[icore].f(rval) / rval;
  //     }
  //     gArray(0) = 2.0 * gArray(1) - gArray(2);
  //     fArray(0) = 2.0 * fArray(1) - fArray(2);
  //     cerr << "Cusp = " 
  // 	   << (gArray(1) - gArray(0))/(gArray(0)*dr) << endl;
  //     double r0 = 1.0e-4;
  //     double g00  = CoreStates[icore].g(r0) / r0;
  //     double dg00 = CoreStates[icore].dgdr(r0) / r0 - 
  // 	CoreStates[icore].g(r0) / (r0*r0);
  //     cerr << "Cusp2 = " << (dg00/g00) << endl;
  //     cerr << "Cusp2 g0  = " << g00  << endl;;
  //     cerr << "Cusp2 dg0 = " << dg00 << endl;
      
  //     int N = CoreStates[icore].r.size();
  //     for (int ir=0; ir<N; ir++) 
  // 	CoreStates[icore].g0(ir) /= CoreStates[icore].r(ir);
  //     out.WriteVar("g",          CoreStates[icore].g0);
  //     out.WriteVar("r",          CoreStates[icore].r);
  //     out.WriteVar("g0",         gArray);
  //     out.WriteVar("f0",         fArray);
  //     out.WriteVar("atom",       CoreStates[icore].atom);
  //     out.WriteVar("n",          CoreStates[icore].n);
  //     out.WriteVar("l",          CoreStates[icore].l);
  //     out.WriteVar("k",          CoreStates[icore].k);
  //     out.WriteVar("eigenvalue", CoreStates[icore].eigenvalue);
      
  //     out.CloseSection(); // "core_state"
  //   }

  //   for (int spin=0; spin<2; spin++) {

  //     // Write localized orbitals first
  //     for (int band=0; band<orbs[spin].extent(1); band++) {
  // 	if (&(orbs[spin](ik,band)->GetCenter()) != NULL) {
  // 	  numLoc += orbs[spin](ik,band)->GetCenter().IdenticalSites.size();
  // 	  out.NewSection("band"); 
  // 	  if (Spline)
  // 	    WriteSpline(out, spin, ik, band);
  // 	  else
  // 	    orbs[spin](ik,band)->Write (out);
  // 	  out.CloseSection(); // "band"
  // 	  numBandsWritten++;
  // 	}
  //     }
  //     if (numLoc == 0)
  // 	NumExtendedBands = orbs[spin].extent(1);
  //     // Now write extended states
  //     for (int band=numLoc; band<numLoc+NumExtendedBands; band++) {
  // 	out.NewSection("band");   
  // 	if (Spline) 
  // 	  WriteSpline (out, spin, ik, band);
  // 	else
  // 	  orbs[spin](ik, band)->Write (out);
  // 	numBandsWritten++;
  // 	out.CloseSection(); // "band"
  //     }
  //   } // spin loop
  //   out.CloseSection(); // "twist"
  // }

  // out.CloseSection(); // "eigenstates"
  // out.OpenSection("parameters");
  // out.WriteVar ("num_bands", numBandsWritten/orbs[0].extent(0));
  // out.CloseSection();

  // PrimDensity.Write (out);

  // out.CloseFile();
  

  return true;
}
