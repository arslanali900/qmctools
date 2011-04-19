/////////////////////////////////////////////////////////////
// Copyright (C) 2003-2006 Bryan Clark and Kenneth Esler   //
//                                                         //
// This program is free software; you can redistribute it  //
// and/or modify it under the terms of the GNU General     //
// Public License as published by the Free Software        //
// Foundation; either version 2 of the License, or         //
// (at your option) any later version.  This program is    //
// distributed in the hope that it will be useful, but     //
// WITHOUT ANY WARRANTY; without even the implied warranty //
// of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. //  
// See the GNU General Public License for more details.    //
// For more information, please see the PIMC++ Home Page:  //
//           http://pathintegrals.info                     //
/////////////////////////////////////////////////////////////

#ifndef PA_TRICUBIC_FIT2_H
#define PA_TRICUBIC_FIT2_H
#ifdef MAKE_FIT
#include "../MPI/Communication.h"
#endif
#include "PAFitBase.h"
#include "../Splines/MyTricubicSpline.h"

class PAtricubicFit2Class : public PairActionFitClass
{
private:
  bool GridsAreMine;
  CommunicatorClass Comm;

  /** These arrays hold the coefficients of the expansion of the tail
      of the diagonal part of the action/potential in inverse powers
      of q or r.  The first index is for the inverse power, and the
      second gives the level.  That is, outside rcut, 
      \f[ U_\text{diag}(q, \text{level}) \approx Ucoef(0,level)/q +
      Ucoef(1,level)/r^2 + Ucoef(2,leve)/r^3 \f] */
  Array<double,2> Ucoefs, dUcoefs;
  Array<double,1> Vcoefs;
  /// The cutoff radius for the long range breakup.
  double rcut;
  inline void zs2yt(double q, double z, double s, int level,
		    double &y, double &t);
  inline void yt2zs(double q, double y, double t, int level,
		    double &z, double &s);
  /// This computes the partial derivative matrix.  Returs the matrix 
  /// [ dq/dq  dy/dq  dt/dq;
  ///   dq/dz  dy/dz  dt/dz;
  ///   dq/ds  dy/ds  dt/ds ]
  void PartialDerivs (double q, double z, double s, int level,
		      TinyMatrix<double,3,3> &derivs);
  Array<double,1> sMax, sMaxInv;
public:
  Grid *qgrid, *ygrid, *tgrid;
  Array<MyTricubicSpline,1> Usplines, dUsplines;

  void ReadParams  (IOSectionClass &inSection);
  void WriteBetaIndependentInfo (IOSectionClass &outSection);
  /// Returns weighter RMS error
  void Error (Rho &rho, double &Uerror, double &dUerror);
  void DoFit (Rho &rho);
  void WriteFit(IOSectionClass &outSection);

  void Write (IOSectionClass &outSection);
  bool Read  (IOSectionClass &inSection, double lowestBeta,
	      int NumBetas);
  double U(double q, double z, double s2, int level);
  double dU(double q, double z, double s2, int level);
  void Derivs (double q, double z, double s2, int level,
	       double &d_dq, double &d_dz);
  void Derivs (double q, double z, double s2, int level,
	       double &d_dq, double &d_dz, double &d_ds);
  /// The diagonal action only -- used for long-range breakup
  double Udiag(double q, int level);
  /// The q-derivative of the above
  double Udiag_p(double q, int level);
  /// The q-derivative of the above
  double Udiag_pp(double q, int level);
  /// The beta-derivative of the diagonal action
  double dUdiag    (double q, int level);
  /// The q-derivative of the above
  double dUdiag_p  (double q, int level);
  /// The q-derivative of the above
  double dUdiag_pp (double q, int level);

  bool IsLongRange();


  double Xk_U     (double k, int level);
  double dXk_U_dk (double k, int level);
  double Xk_dU    (double k, int level);
  double Xk_V     (double k);
  double Vk       (double k);
  void Setrc (double rc);

  PAtricubicFit2Class()
  { 
#ifdef MAKE_FIT
    Comm.SetWorld();
#endif
    GridsAreMine = false; 
    NumBetas=0;
  }
  ~PAtricubicFit2Class()
  { if (GridsAreMine){ delete qgrid; delete ygrid; } }
};

inline void
PAtricubicFit2Class::zs2yt(double q,  double z, double s, int level,
			   double &y, double &t)
{
  z = fabs(z);
  double smax = min(2.0*q, sMax(level));
  /////////////////////////////////////////
  // CHECK THIS!!!!!!!!!!!!!!!!!!!!!!!!! //
  /////////////////////////////////////////
  //  y = z*sMaxInv(level);
  if (smax == 0.0) {
    y = 0.0;
    t = 0.0;
  }
  else {
    y = z/smax;
    t = (y < 1.0) ? (s-z)/(smax-z) : (s/smax);
  }

//   if (y < 1.0) 
//     t = (s-z)/(smax-z);
//   else
//     t = s/smax;
}

inline void
PAtricubicFit2Class::yt2zs(double q,  double y, double t, int level,
			   double &z, double &s)
{
  double smax = min(2.0*q, sMax(level));
  z = smax*y;
  // s = smax*y + t*(1.0-y)*smax;
  // s = smax*(y+t(1.0-y));
  s = smax*(y+t-t*y);
}

#endif
