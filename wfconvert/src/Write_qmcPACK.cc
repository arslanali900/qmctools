#include "OrbitalSetClass.h"
#include <Common/IO/IO.h>
#include "MuffinTin.h"


// This function creates the coefficients necessary to shift orbitals
// by rshift by multiplying their FFT box in k-space by these
// coefficients.  They are returned in shiftBox.
void
OrbitalSetClass::MakeBoxShifter (Vec3 rshift, Array<complex<double>,3> &shiftBox)
{
  CellClass &cell = WriteSuper ? SuperCell : PrimCell;
  FFTBox &FFT = cell.FFT;

  Vec3 ushift = cell.Lattice.r2u(rshift);
  Int3 dims;
  FFT.GetDims(dims[0], dims[1], dims[2]);
  if (shiftBox.shape() != dims)
    shiftBox.resize(dims);

  Array<complex<double>,1> shiftVecs[3];
  for (int i=0; i<3; i++) {
    int N = dims[i];
    shiftVecs[i].resize(N);
    int jstart = -(N-1)/2;
    int jend   =      N/2;
    for (int j=jstart; j<=jend; j++) {
      double phase = (double)j*2.0*M_PI*ushift[i];
      shiftVecs[i]((j+N)%N) = complex<double>(cos(phase), sin(phase));
    }
  }
  for (int ix=0; ix<dims[0]; ix++)
    for (int iy=0; iy<dims[1]; iy++)
      for (int iz=0; iz<dims[2]; iz++)
	shiftBox(ix,iy,iz) = shiftVecs[0](ix)*shiftVecs[1](iy)*shiftVecs[2](iz);
  
//   for (int ix=-dims[0]/2; ix<dims[0]/2; ix++)
//     for (int iy=-dims[1]/2; iy<dims[1]/2; iy++)
//       for (int iz=-dims[2]/2; iz<dims[2]/2; iz++) {
// 	double phase = 2.0*M_PI*(ushift[0]*(double)ix+
// 				  ushift[1]*(double)iy+
// 				  ushift[2]*(double)iz);
// 	shiftBox((ix+dims[0])%dims[0],
// 		 (iy+dims[1])%dims[1],
// 		 (iz+dims[2])%dims[2]) = complex<double>(cos(phase), sin(phase));
//      }
}

void
OrbitalSetClass::MakeFineOrbital (int spin, int ki, int band, bool shiftOrbs,
				  complex<double> rotation)
{
  CellClass &cell = FineCell;
  FFTBox &FFT = cell.FFT;
  TinyVector<Array<OrbitalClass*,2>,2> &orbs = FineOrbitals;
  cell.SetupFFT();

  zVec coefs(orbs[spin](0,0)->GetCoefs().size());
  FFT.kBox = complex<double>();

  OrbitalClass &myOrb = *orbs[spin](ki,band);

  if (band < NewOrbCoefs[spin].extent(1)) {
    for (int bi=0; bi<NewOrbCoefs[spin].extent(2); bi++) {
      OrbitalClass &orb = *orbs[spin](ki,bi);
      coefs = orb.GetCoefs();
      coefs *= NewOrbCoefs[spin](ki,band,bi);
      orb.AddToFFTBox(coefs);
    }
  }
  else {
    OrbitalClass &orb = *orbs[spin](ki,band);
    coefs = orb.GetCoefs();
    orb.PutInFFTBox (coefs);
  }

  if (shiftOrbs && (&myOrb.GetCenter()) != NULL) {
    Vec3 rc = myOrb.GetCenter().r;
    Vec3 rhalf = cell.Lattice.u2r(Vec3(0.5, 0.5, 0.5));
    Array<complex<double>,3> shiftBox;
    MakeBoxShifter (rhalf-rc, shiftBox);
    int nx, ny, nz;
    FFT.GetDims(nx,ny,nz);
    for (int ix=0; ix<nx; ix++)
      for (int iy=0; iy<ny; iy++)
	for (int iz=0; iz<nz; iz++)
	  FFT.kBox(ix,iy,iz) = shiftBox(ix,iy,iz) * FFT.kBox(ix,iy,iz);
  }
  FFT.k2r();

  double realNorm = 0.0;
  double imagNorm = 0.0;
  int nx, ny, nz;
  FFT.GetDims(nx,ny,nz);
  double nxInv = 1.0/(double)nx;
  double nyInv = 1.0/(double)ny;
  double nzInv = 1.0/(double)nz;
  for (int ix=0; ix<nx; ix++)
    for (int iy=0; iy<ny; iy++)
      for (int iz=0; iz<nz; iz++)
	FFT.rBox(ix,iy,iz) *= rotation;

  
  realNorm = imagNorm = 0.0;
  double Norm = 0.0;
  for (int ix=0; ix<nx; ix++) {
    double phi_x = 2.0*M_PI*myOrb.GetTwist()[0]*(double)ix*nxInv;
    for (int iy=0; iy<ny; iy++) {
      double phi_y = 2.0*M_PI*myOrb.GetTwist()[1]*(double)iy*nyInv;
      for (int iz=0; iz<nz; iz++) {
	double phi_z = 2.0*M_PI*myOrb.GetTwist()[2]*(double)iz*nzInv;
	double phase = phi_x + phi_y + phi_z;
	double s,c;
	sincos(-phase, &s, &c);
	complex<double> e_mikr(c,s);
	complex<double> z = FFT.rBox(ix,iy,iz)*e_mikr;
	realNorm += z.real()*z.real();
	imagNorm += z.imag()*z.imag();
	Norm += norm(z);
      }
    }
  }

  double vol = cell.Lattice.GetVolume();
  double prefact = 1.0/(double)(nx*ny*nz);

  fprintf (stderr, "ki=%d band=%d real norm = %8.4e  imag norm = %8.4e  ratio=%1.8f   Norm=%1.8f\n",
	   ki, band, realNorm, imagNorm, realNorm/imagNorm, prefact*Norm);
}



void
OrbitalSetClass::MakeFineLaplacian (int spin, int ki, int band, bool shiftOrbs,
				    complex<double> rotation)
{
  CellClass &cell = FineCell;
  FFTBox &FFT = cell.FFT;
  TinyVector<Array<OrbitalClass*,2>,2> &orbs = FineOrbitals;
  cell.SetupFFT();

  zVec coefs(orbs[spin](ki,band)->GetCoefs().size());
  FFT.kBox = complex<double>();

  OrbitalClass &myOrb = *orbs[spin](ki,band);

  if (band < NewOrbCoefs[spin].extent(1)) {
    for (int bi=0; bi<NewOrbCoefs[spin].extent(2); bi++) {
      OrbitalClass &orb = *orbs[spin](ki,bi);
      coefs = orb.GetCoefs();
      coefs *= NewOrbCoefs[spin](ki,band,bi);
      for (int iG=0; iG<coefs.size(); iG++)
	coefs(iG) *= -dot(FFT.GVecs(iG), FFT.GVecs(iG));

      orb.AddToFFTBox(coefs);
    }
  }
  else {
    OrbitalClass &orb = *orbs[spin](ki,band);
    coefs = orb.GetCoefs();
    for (int iG=0; iG<coefs.size(); iG++)
      coefs(iG) *= -dot(FFT.GVecs(iG), FFT.GVecs(iG));
    orb.PutInFFTBox (coefs);
  }

  if (shiftOrbs && (&myOrb.GetCenter()) != NULL) {
    Vec3 rc = myOrb.GetCenter().r;
    Vec3 rhalf = cell.Lattice.u2r(Vec3(0.5, 0.5, 0.5));
    Array<complex<double>,3> shiftBox;
    MakeBoxShifter (rhalf-rc, shiftBox);
    int nx, ny, nz;
    FFT.GetDims(nx,ny,nz);
    for (int ix=0; ix<nx; ix++)
      for (int iy=0; iy<ny; iy++)
	for (int iz=0; iz<nz; iz++)
	  FFT.kBox(ix,iy,iz) = shiftBox(ix,iy,iz) * FFT.kBox(ix,iy,iz);
  }

  FFT.k2r();

  double realNorm = 0.0;
  double imagNorm = 0.0;
  int nx, ny, nz;
  FFT.GetDims(nx,ny,nz);
  double nxInv = 1.0/(double)nx;
  double nyInv = 1.0/(double)ny;
  double nzInv = 1.0/(double)nz;
  for (int ix=0; ix<nx; ix++)
    for (int iy=0; iy<ny; iy++)
      for (int iz=0; iz<nz; iz++)
	FFT.rBox(ix,iy,iz) *= rotation;
}

				  

// Returns the rotation factor, which will will need if we to the APW part
complex<double>
OrbitalSetClass::MakeRealSpaceOrbital (int spin, int ki, int band, 
				       bool shiftOrbs)
{
  CellClass &cell = WriteSuper ? SuperCell : PrimCell;
  FFTBox &FFT = cell.FFT;
  TinyVector<Array<OrbitalClass*,2>,2> &orbs = 
    WriteSuper ? SuperOrbitals : PrimOrbitals;

  cell.SetupFFT();

  zVec coefs(orbs[spin](0,0)->GetCoefs().size());
  FFT.kBox = complex<double>();

  OrbitalClass &myOrb = *orbs[spin](ki,band);

  if (band < NewOrbCoefs[spin].extent(1)) {
    for (int bi=0; bi<NewOrbCoefs[spin].extent(2); bi++) {
      OrbitalClass &orb = *orbs[spin](ki,bi);
      coefs = orb.GetCoefs();
      coefs *= NewOrbCoefs[spin](ki,band,bi);
//       if (shiftOrbs && (&myOrb.GetCenter()) != NULL) {
// 	Vec3 rc = myOrb.GetCenter().r;
// 	Vec3 rhalf = cell.Lattice.u2r(Vec3(0.5, 0.5, 0.5));
// 	for (int i=0; i<coefs.size(); i++) {
// 	  Vec3 G = FFT.GVecs(i)+orb.GetGshift();
// 	  double phase = -dot (G, rc-rhalf);
// 	  double s,c;
// 	  sincos (phase, &s, &c);
// 	  coefs(i) *= complex<double>(c,s);
// 	}
//       }
      orb.AddToFFTBox(coefs);
    }
  }
  else {
    OrbitalClass &orb = *orbs[spin](ki,band);
    coefs = orb.GetCoefs();
//     if (shiftOrbs && (&orb.GetCenter()) != NULL) {
//       Vec3 rc = orb.GetCenter().r;
//       Vec3 rhalf = cell.Lattice.u2r(Vec3(0.5, 0.5, 0.5));
//       for (int i=0; i<coefs.size(); i++) {
// 	Vec3 G = FFT.GVecs(i)+orb.GetGshift();
// 	double phase = -dot (G, rc-rhalf);
// 	double s,c;
// 	sincos (phase, &s, &c);
// 	coefs(i) *= complex<double>(c,s);
//       }
//     }
    orb.PutInFFTBox (coefs);
  }

  if (shiftOrbs && (&myOrb.GetCenter()) != NULL) {
    Vec3 rc = myOrb.GetCenter().r;
    Vec3 rhalf = cell.Lattice.u2r(Vec3(0.5, 0.5, 0.5));
    Array<complex<double>,3> shiftBox;
    MakeBoxShifter (rhalf-rc, shiftBox);
    int nx, ny, nz;
    FFT.GetDims(nx,ny,nz);
    // MakeBoxShifter (rc-rhalf, shiftBox);
    for (int ix=0; ix<nx; ix++)
      for (int iy=0; iy<ny; iy++)
	for (int iz=0; iz<nz; iz++)
	  FFT.kBox(ix,iy,iz) = shiftBox(ix,iy,iz) * FFT.kBox(ix,iy,iz);
  }

//   else {
//     if (band < NewOrbCoefs[spin].extent(1)) {
//       for (int bi=0; bi<NewOrbCoefs[spin].extent(2); bi++)
// 	orbs[spin](ki,bi)->AddToFFTBox(NewOrbCoefs[spin](ki,band,bi));
//     }
//     else 
//      orbs[spin](ki,band)->PutInFFTBox();
//   }

  FFT.k2r();

  double realNorm = 0.0;
  double imagNorm = 0.0;
  int nx, ny, nz;
  FFT.GetDims(nx,ny,nz);
  double nxInv = 1.0/(double)nx;
  double nyInv = 1.0/(double)ny;
  double nzInv = 1.0/(double)nz;
  // Rotate the orbitals by a constant phase to make sure that neither
  // the real nor the imaginary part is strictly zero.  
  // Check real and imaginary norms
  for (int ix=0; ix<nx; ix++) {
    double phi_x = 2.0*M_PI*myOrb.GetTwist()[0]*(double)ix*nxInv;
    for (int iy=0; iy<ny; iy++) {
      double phi_y = 2.0*M_PI*myOrb.GetTwist()[1]*(double)iy*nyInv;
      for (int iz=0; iz<nz; iz++) {
	double phi_z = 2.0*M_PI*myOrb.GetTwist()[2]*(double)iz*nzInv;
	double phase = phi_x + phi_y + phi_z;
	double s,c;
	sincos(-phase, &s, &c);
	complex<double> e_mikr(c,s);
	complex<double> z = FFT.rBox(ix,iy,iz)*e_mikr;
	realNorm += z.real()*z.real();
	imagNorm += z.imag()*z.imag();
      }
    }
  }

  double theta = atan2(imagNorm, realNorm);
  theta = 0.5*(M_PI/4.0 - theta);

  complex<double> rotation (cos(theta), sin(theta));

  for (int ix=0; ix<nx; ix++)
    for (int iy=0; iy<ny; iy++)
      for (int iz=0; iz<nz; iz++)
	FFT.rBox(ix,iy,iz) *= rotation;

  
  realNorm = imagNorm = 0.0;
  double Norm = 0.0;
  for (int ix=0; ix<nx; ix++) {
    double phi_x = 2.0*M_PI*myOrb.GetTwist()[0]*(double)ix*nxInv;
    for (int iy=0; iy<ny; iy++) {
      double phi_y = 2.0*M_PI*myOrb.GetTwist()[1]*(double)iy*nyInv;
      for (int iz=0; iz<nz; iz++) {
	double phi_z = 2.0*M_PI*myOrb.GetTwist()[2]*(double)iz*nzInv;
	double phase = phi_x + phi_y + phi_z;
	double s,c;
	sincos(-phase, &s, &c);
	complex<double> e_mikr(c,s);
	complex<double> z = FFT.rBox(ix,iy,iz)*e_mikr;
	realNorm += z.real()*z.real();
	imagNorm += z.imag()*z.imag();
	Norm += norm(z);
      }
    }
  }

  double vol = cell.Lattice.GetVolume();
  double prefact = 1.0/(double)(nx*ny*nz);

  fprintf (stderr, "ki=%d band=%d real norm = %8.4e  imag norm = %8.4e  ratio=%1.8f   Norm=%1.8f\n",
	   ki, band, realNorm, imagNorm, realNorm/imagNorm, prefact*Norm);

  // if (UseMultiRep)
  //   MakeFineOrbital (spin, ki, band, shiftOrbs, rotation);

  return rotation;
}

// Make a first order wave function in real space
void
OrbitalSetClass::MakeFirstOrderOrbital (int ideriv, int spin, int ki, int band, 
					bool shiftOrbs, complex<double> rotation)
{
  CellClass &cell = WriteSuper && false ? SuperCell : PrimCell;
  FFTBox &FFT = cell.FFT;
  TinyVector<Array<OrbitalClass*,2>,2> &orbs = 
    WriteSuper && false ? SuperFirstOrder[ideriv] : PrimFirstOrder[ideriv];

  cell.SetupFFT();

  zVec coefs(orbs[spin](0,0)->GetCoefs().size());
  FFT.kBox = complex<double>();

  OrbitalClass &myOrb = *orbs[spin](ki,band);

  if (band < NewOrbCoefs[spin].extent(1)) {
    for (int bi=0; bi<NewOrbCoefs[spin].extent(2); bi++) {
      OrbitalClass &orb = *orbs[spin](ki,bi);
      coefs = orb.GetCoefs();
      coefs *= NewOrbCoefs[spin](ki,band,bi);
      orb.AddToFFTBox(coefs);
    }
  }
  else {
    OrbitalClass &orb = *orbs[spin](ki,band);
    coefs = orb.GetCoefs();
    orb.PutInFFTBox (coefs);
  }

  if (shiftOrbs && (&myOrb.GetCenter()) != NULL) {
    Vec3 rc = myOrb.GetCenter().r;
    Vec3 rhalf = cell.Lattice.u2r(Vec3(0.5, 0.5, 0.5));
    Array<complex<double>,3> shiftBox;
    MakeBoxShifter (rhalf-rc, shiftBox);
    int nx, ny, nz;
    FFT.GetDims(nx,ny,nz);
    // MakeBoxShifter (rc-rhalf, shiftBox);
    for (int ix=0; ix<nx; ix++)
      for (int iy=0; iy<ny; iy++)
	for (int iz=0; iz<nz; iz++)
	  FFT.kBox(ix,iy,iz) = shiftBox(ix,iy,iz) * FFT.kBox(ix,iy,iz);
  }

  FFT.k2r();

  double realNorm = 0.0;
  double imagNorm = 0.0;
  int nx, ny, nz;
  FFT.GetDims(nx,ny,nz);
  double nxInv = 1.0/(double)nx;
  double nyInv = 1.0/(double)ny;
  double nzInv = 1.0/(double)nz;

  for (int ix=0; ix<nx; ix++)
    for (int iy=0; iy<ny; iy++)
      for (int iz=0; iz<nz; iz++)
	FFT.rBox(ix,iy,iz) *= rotation;

  
  realNorm = imagNorm = 0.0;
  double Norm = 0.0;
  for (int ix=0; ix<nx; ix++) {
    double phi_x = 2.0*M_PI*myOrb.GetTwist()[0]*(double)ix*nxInv;
    for (int iy=0; iy<ny; iy++) {
      double phi_y = 2.0*M_PI*myOrb.GetTwist()[1]*(double)iy*nyInv;
      for (int iz=0; iz<nz; iz++) {
	double phi_z = 2.0*M_PI*myOrb.GetTwist()[2]*(double)iz*nzInv;
	double phase = phi_x + phi_y + phi_z;
	double s,c;
	sincos(-phase, &s, &c);
	complex<double> e_mikr(c,s);
	complex<double> z = FFT.rBox(ix,iy,iz)*e_mikr;
	realNorm += z.real()*z.real();
	imagNorm += z.imag()*z.imag();
	Norm += norm(z);
      }
    }
  }

  double vol = cell.Lattice.GetVolume();
  double prefact = 1.0/(double)(nx*ny*nz);

  fprintf (stderr, "ki=%d band=%d real norm = %8.4e  imag norm = %8.4e  ratio=%1.8f   Norm=%1.8f\n",
	   ki, band, realNorm, imagNorm, realNorm/imagNorm, prefact*Norm);
}



void
OrbitalSetClass::WriteSpline (IO::IOSectionClass &out, 
			      int spin, int ki, int band)
{
  cerr << "In OrbitalSetClass::WriteSpline\n";
  int nx, ny, nz;
  complex<double> phase_rotation =
    MakeRealSpaceOrbital (spin, ki, band, ShiftOrbitals || Truncate);
  TinyVector<Array<OrbitalClass*,2>,2> &orbs = 
    WriteSuper ? SuperOrbitals : PrimOrbitals;
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
    Array<double,1> uMin(3), uMax(3);
    uMin(0) = (double)minx*du[0];    uMax(0) = (double)maxx*du[0];
    uMin(1) = (double)miny*du[1];    uMax(1) = (double)maxy*du[1];
    uMin(2) = (double)minz*du[2];    uMax(2) = (double)maxz*du[2];
    out.WriteVar ("umin", uMin);
    out.WriteVar ("umax", uMax);
    out.WriteVar ("eigenvector", eigvec);
    out.WriteVar ("truncated", true);
  }
  else {
    Array<double,4> eigvec(nx+1, ny+1, nz+1, 2);
    for (int ix=0; ix<nx+1; ix++)
      for (int iy=0; iy<ny+1; iy++)
	for (int iz=0; iz<nz+1; iz++) {
	  eigvec(ix, iy, iz, 0) = prefact*FFT.rBox(ix%nx,iy%ny,iz%nz).real();
	  eigvec(ix, iy, iz, 1) = prefact*FFT.rBox(ix%nx,iy%ny,iz%nz).imag();
	}
    // Finally write to file
    out.WriteVar ("eigenvector", eigvec);
    out.WriteVar ("umin", 0.0);
    out.WriteVar ("umax", 1.0);
    out.WriteVar ("truncated", true);
  }
  out.WriteVar ("eigenvalue", eigval);
  out.WriteVar ("spin", spin);
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


      // Create a MuffinTin object to test
      // MuffinTinClass tin;
      // Vec3 k = PrimCell.Lattice.Twist2k(orbs[spin](ki,0)->GetTwist());

      // tin.initAPW (radius, APWRadialPoints, APW.lMax, 1);
      // tin.setAPW (0, k, u_lm_r);
      // tin.setCenter (PrimCell.IonPos(atom));
      // tin.setLattice (PrimCell.Lattice.GetDirect());

      // vector<complex<double> > phi(1);
      // TinyVector<double,3> rtry(0.1, 0.2, 0.3);
      // tin.evaluate(rtry, phi);
//       fprintf (stderr, "phi(0.1,0.2,0.3) = %16.12e + %16.12e\n",
// 	       phi[0].real(), phi[0].imag());
    }
  }
  out.FlushFile();
  // out.WriteVar ("prim_twist_index", PrimTwistIndex);
  // out.WriteVar ("prim_band_index",  PrimBandIndex);
}


bool
OrbitalSetClass::Write_qmcPACK (string fname)
{

  if (Comm.MyProc() != 0)
    return true;
  CellClass &cell = WriteSuper ? SuperCell : PrimCell;
  FFTBox &FFT = cell.FFT;
  GVecsClass &GVecs = cell.GVecs;
  TinyVector<Array<OrbitalClass*,2>,2> &orbs =
    WriteSuper ? SuperOrbitals : PrimOrbitals;

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
      lattice(i,j) = cell.Lattice.GetDirect()(i,j);
  out.WriteVar("lattice", lattice);
  Array<double,2> reciprocal_lattice(3,3);
  for (int i=0; i<3; i++)
    for (int j=0; j<3; j++)
      reciprocal_lattice(i,j) = cell.Lattice.GetRecip()(i,j);
  out.WriteVar("reciprocal_lattice", reciprocal_lattice);
  int complex_coefficients = 1;
  out.WriteVar ("complex_coefficients", complex_coefficients);
  //  out.WriteVar ("num_bands",      orbs[0].extent(1)+orbs[1].extent(1));
  out.WriteVar ("num_up_bands",   orbs[0].extent(1));
  out.WriteVar ("num_down_bands", orbs[1].extent(1));
  // Write out the number of core states
  int numCore = CoreStates.size();
  out.WriteVar ("num_core_states", numCore);
  out.WriteVar ("num_twists",     orbs[0].extent(0));
  int num_spins = SpinPolarized ? 2 : 1;
  out.WriteVar ("num_spins", num_spins);
  out.WriteVar ("maximum_ecut", ECut);
  int numElecs = NumElectrons[0] + NumElectrons[1];
  if (WriteSuper)
    numElecs *= (int)round(fabs(det(TileMatrix)));
  out.WriteVar ("num_electrons", numElecs);

  // Write muffin-tin information
  if (APW.lMax > 0) {
    out.NewSection("muffin_tins");
    out.WriteVar("num_tins", APW.NumAtoms);
    for (int atom=0; atom<APW.NumAtoms; atom++) {
      out.NewSection("muffin_tin");
      out.WriteVar("radius", APW.u(atom,0)[0].grid.End());
      out.WriteVar("r", APW.r_data(atom));
      Array<double,1> center(3);
      for (int i=0; i<3; i++)
	center(i) = cell.IonPos(atom)[i];
      out.WriteVar("center", center);
      out.WriteVar("lmax", APW.lMax);
      out.WriteVar("num_radial_points", APWRadialPoints);
      out.CloseSection(); // "muffin_tin"
    }

    out.CloseSection(); // "muffin_tins"
  }

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
  Array<double,2> pos(cell.IonPos.size(),3);
  for (int i=0; i<cell.IonPos.size(); i++) {
    pos(i,0) = cell.IonPos(i)[0];
    pos(i,1) = cell.IonPos(i)[1];
    pos(i,2) = cell.IonPos(i)[2];
  }
  out.WriteVar ("pos", pos);
  out.WriteVar ("atom_types", cell.AtomTypes);
  out.CloseSection(); // "ions"
  

  out.NewSection ("eigenstates");
  int numBandsWritten = 0;
  for (int icore=0; icore<CoreStates.size(); icore++)
    CoreStates[icore].SetupSplines();

  for (int ik=0; ik<orbs[0].extent(0); ik++) {
    out.NewSection("twist");
    Vec3 twist = orbs[0](ik,0)->GetTwist();
    
    Array<double,1> twist_angle(3);
    twist_angle(0)=twist[0]; 
    twist_angle(1)=twist[1]; 
    twist_angle(2)=twist[2];
    out.WriteVar ("twist_angle", twist_angle);
    int numLoc=0;
    // Write core states, if we have them
    for (int icore=0; icore<CoreStates.size(); icore++) {
      out.NewSection("core_state");
      double rmax = CoreStates[icore].r(CoreStates[icore].r.size()-1);
      out.WriteVar("rmax", rmax);
      Array<double,1> fArray(CoreRadialPoints), gArray(CoreRadialPoints);
      double dr = rmax / (double)(CoreRadialPoints-1);
      // Find ir for around r=0.1
      int icut = (int)ceil(0.05/dr);
      double rcut = dr * icut;
      double gcut = CoreStates[icore].g(rcut) / rcut;
      for (int ir=1; ir<CoreRadialPoints; ir++) {
	double rval = dr * (double)ir;
	double gval = CoreStates[icore].g(rval) / rval;
	double fval = CoreStates[icore].f(rval) / rval;
	gArray(ir) = CoreStates[icore].g(rval) / rval;
	fArray(ir) = CoreStates[icore].f(rval) / rval;
      }
      gArray(0) = 2.0 * gArray(1) - gArray(2);
      fArray(0) = 2.0 * fArray(1) - fArray(2);
      cerr << "Cusp = " 
	   << (gArray(1) - gArray(0))/(gArray(0)*dr) << endl;
      double r0 = 1.0e-4;
      double g00  = CoreStates[icore].g(r0) / r0;
      double dg00 = CoreStates[icore].dgdr(r0) / r0 - 
	CoreStates[icore].g(r0) / (r0*r0);
      cerr << "Cusp2 = " << (dg00/g00) << endl;
      cerr << "Cusp2 g0  = " << g00  << endl;;
      cerr << "Cusp2 dg0 = " << dg00 << endl;
      
      int N = CoreStates[icore].r.size();
      for (int ir=0; ir<N; ir++) 
	CoreStates[icore].g0(ir) /= CoreStates[icore].r(ir);
      out.WriteVar("g",          CoreStates[icore].g0);
      out.WriteVar("r",          CoreStates[icore].r);
      out.WriteVar("g0",         gArray);
      out.WriteVar("f0",         fArray);
      out.WriteVar("atom",       CoreStates[icore].atom);
      out.WriteVar("n",          CoreStates[icore].n);
      out.WriteVar("l",          CoreStates[icore].l);
      out.WriteVar("k",          CoreStates[icore].k);
      out.WriteVar("eigenvalue", CoreStates[icore].eigenvalue);
      
      out.CloseSection(); // "core_state"
    }

    for (int spin=0; spin<2; spin++) {

      // Write localized orbitals first
      for (int band=0; band<orbs[spin].extent(1); band++) {
	if (&(orbs[spin](ik,band)->GetCenter()) != NULL) {
	  numLoc += orbs[spin](ik,band)->GetCenter().IdenticalSites.size();
	  out.NewSection("band"); 
	  if (Spline)
	    WriteSpline(out, spin, ik, band);
	  else
	    orbs[spin](ik,band)->Write (out);
	  out.CloseSection(); // "band"
	  numBandsWritten++;
	}
      }
      if (numLoc == 0)
	NumExtendedBands = orbs[spin].extent(1);
      // Now write extended states
      for (int band=numLoc; band<numLoc+NumExtendedBands; band++) {
	out.NewSection("band");   
 	if (Spline) 
	  WriteSpline (out, spin, ik, band);
 	else
 	  orbs[spin](ik, band)->Write (out);
	numBandsWritten++;
	out.CloseSection(); // "band"
      }
    } // spin loop
    out.CloseSection(); // "twist"
  }

  out.CloseSection(); // "eigenstates"
  out.OpenSection("parameters");
  out.WriteVar ("num_bands", numBandsWritten/orbs[0].extent(0));
  out.CloseSection();

  PrimDensity.Write (out);

  out.CloseFile();
  

  return true;
}
