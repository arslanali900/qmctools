#ifndef ATOMIC_ORBITAL_H
#define ATOMIC_ORBITAL_H

#include <Common/Blitz.h>
#include <vector>
#include <einspline/bspline.h>
#include "LatticeClass.h"
#include "ShiftedChebyshev.h"
#include <Common/IO/IO.h>
#include "OrbFunctor.h"


class QuadratureRule {
private:
  void addknot (Vec3 r, double w) {
    rhat.push_back(r);
    weight.push_back(w);
  }
  void MakeLebedev(int rule);
  void Check();
  int lexact;
public:
  std::vector<Vec3> rhat;
  std::vector<double> weight;
  
  QuadratureRule() {}

  QuadratureRule(int nrule) 
  { 
    SetRule (nrule);
  }
  void SetRule (int rule);

};


// class OrbFunctor
// {
// public:
//   virtual complex<double> operator()(Vec3 r) = 0;
//   virtual void operator()(Vec3 r, complex<double> &val,
// 			  TinyVector<complex<double>,3> &grad,
// 			  complex<double> &lapl) = 0;
//   virtual void operator()(vector<Vec3> &r,
// 			  vector<complex<double> > &val)
//   {
//     for (int ir=0; ir<r.size(); ir++)
//       val[ir] = (*this)(r[ir]);
//   }

//   virtual void operator()(vector<Vec3> &r,
// 			  vector<complex<double> > &val,
// 			  vector<complex<double> > &lapl)
//   {
//     TinyVector<complex<double>,3> grad;
//     for (int ir=0; ir<r.size(); ir++)
//       (*this)(r[ir], val[ir], grad, lapl[ir]);
//   }

//   virtual void operator()(vector<Vec3> &r,
// 			  vector<complex<double> > &val,
// 			  vector<TinyVector<complex<double>,3> > &grad,
// 			  vector<complex<double> > &lapl)
//   {
//     for (int ir=0; ir<r.size(); ir++)
//       (*this)(r[ir], val[ir], grad[ir], lapl[ir]);
//   }

// };

class OrbitalClass;

class AtomicOrbital
{
private:
  Vec3 Pos;
  double Radius, OuterRadius, PolyRadius;
  int SplinePoints;
  int PolyOrder;
  int lMax;
  // Store in order 
  // Index = l*(l+1) + m.  There are (lMax+1)^2 Ylm's
  vector<complex<double> > YlmVec, dYlmVec;
  Array<complex<double>,2> YlmArray;
  void CalcYlm(Vec3 rhat);
  QuadratureRule QuadRule;
  multi_UBspline_1d_z *Spline;
  // The first index is n in T_2n(x), the second is lm = l*(l+1)+m
  Array<complex<double>,2> ChebyshevCoefs;
  ShiftedChebyshevPoly<complex<double> > Chebyshev;
  vector<double> Tn, dTn, d2Tn;
public:
  // The first index is n in r^n, the second is lm = l*(l+1)+m
  Array<complex<double>,2> PolyCoefs;

  void TimeYlm();

  void Set (int lmax, double radius, int spline_points, Vec3 pos)
  {
    lMax = lmax;
    YlmVec.resize((lmax+1) * (lmax+1));
    dYlmVec.resize((lmax+1) * (lmax+1));
    YlmArray.resize(YlmVec.size(), QuadRule.rhat.size());
    PolyCoefs.resize(PolyOrder+1, YlmVec.size());
    Radius = radius;
    OuterRadius = 1.2 * radius;
    SplinePoints = spline_points;
    Pos = pos;
    BCtype_z bc;
    bc.lCode = FLAT;
    bc.rCode = NATURAL;
    Ugrid grid;
    grid.start = 0.0;
    grid.end   = OuterRadius;
    grid.num   = SplinePoints;
    if (!Spline)
      Spline = create_multi_UBspline_1d_z (grid, bc, YlmVec.size());
  }

  // This function sets the radial functions array variable to the
  // Ylm projections of the atom centered at Pos.
  // The first index of radialFuncs is the radial distance index, and the second
  // index is (l*(l+1)+m) and the second.
  void Project (OrbFunctor &phi, OrbFunctor &spline,
		Array<complex<double>,2> &radialFuncs);
  void
  ProjectAnalytic (PWOrb_z &phi, OrbFunctor &spline,
		   Array<complex<double>,2> &radialFuncs);
  
  void ProjectAnalytic (int iat, TinyVector<Array<OrbitalClass*,2>,2> &orbitals);
#ifdef HAVE_CUDA
  void ProjectAnalyticCudaDebug (int iat, TinyVector<Array<OrbitalClass*,2>,2> &orbitals);
#endif

  complex<double> eval (Vec3 r);
  complex<double> eval_poly (Vec3 r);
  void eval (Vec3 r,
	     complex<double> &val,
	     TinyVector<complex<double>,3> &grad,
	     complex<double> &lapl);
  void eval_poly (Vec3 r, complex<double> &val,
		  TinyVector<complex<double>,3> &grad,
		  complex<double> &lapl);

  void eval (Vec3 r, LatticeClass &lattice,
	     complex<double> &val,
	     TinyVector<complex<double>,3> &grad,
	     complex<double> &lapl);

  complex<double> eval_cheb(Vec3 r);
  void eval_cheb (Vec3 r, complex<double> &val,
		  TinyVector<complex<double>,3> &grad,
		  complex<double> &lapl);

  Vec3 get_pos () const { return Pos; }
  double get_radius() const { return Radius; }
  int get_lmax() const { return lMax; }
  int get_splinePoints() const { return SplinePoints; }
  int get_polyOrder() const { return PolyOrder; }
  double get_polyRadius() const { return PolyRadius; }
  double get_outerRadius() const { return OuterRadius; }
  void Write_ESHDF (IO::IOSectionClass &out);

  AtomicOrbital (const AtomicOrbital &orb)  : 
    QuadRule(8), Spline((multi_UBspline_1d_z*)1),
    PolyOrder(7), PolyRadius(0.1)
  {
    Set (orb.lMax, orb.Radius, orb.SplinePoints, orb.Pos);
    Spline = orb.Spline;
  }

  AtomicOrbital() : QuadRule(8), Spline(NULL), PolyOrder(7), PolyRadius(0.1)
  {
    // Nothing else for now
  }

};

#ifdef HAVE_CUDA
template <typename T> struct vec3         { typedef float   Type; }; // dummy 
template <>           struct vec3<float>  { typedef float3  Type; };
template <>           struct vec3<double> { typedef double3 Type; };

void ProjectAnalyticCuda(int numats, int lMax, int SplinePoints, int PolyOrder, double PolyRadius, 
			 double OuterRadius, double* positions, 
			 TinyVector<Array<OrbitalClass*,2>,2>& orbitals);

void projectAnalyticCudaDriver (const vec3<double>::Type k, double* pos, const double* gvecs, 
				const double* orbitalCoefs, double* splineData, double* polyData,
				int num_ks, int knum, int numG, int lmax, int numats, int SplinePoints, 
				int PolyPoints, int numSpins, int num_bands, double dr, double poly_dr);

void setupCPUComputeKernel(const typename vec3<double>::Type k, const typename vec3<double>::Type pos,
			   double* Ghats, double* phase_shifts, double* jlArray_spline,
			   double* jlArray_poly, double* ylmptr,  int numG, int SplinePoints, int lmax, 
			   int PolyPoints, double dr, double poly_dr);

void projectAnalyticCPUKernel(const double* const orbitalCoefs, const double* const phase_shifts, 
			      const double* const jlArray_spline, const double* const jlArray_poly,
			      const double* const YlmPtr, double* splineData, double* polyData, 
			      int numSpins, int num_bands, int block_size, int numG, int SplinePoints, 
			      int PolyPoints, int lmax);

//void gsl_sf_bessel_jl_steed_array_redir(int lmax, double val, double* jlvec) {
//  gsl_sf_bessel_jl_steed_array(lmax, val, jlvec);
//}

#endif

#endif
