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

#ifndef PA_DIPOLE_FIT_H
#define PA_DIPOLE_FIT_H
#include "../MPI/Communication.h" //include not needed
#include "PAFitBase.h"

class PADipoleFitClass : public PairActionFitClass
{
private:
  // Real space cutoff parameter;
  double alpha;
  // The cutoff radius
  double rcut;  
public:

  void ReadParams  (IOSectionClass &inSection);
  void WriteBetaIndependentInfo (IOSectionClass &outSection);
  /// Returns weighter RMS error
  void Error (Rho &rho, double &Uerror, double &dUerror);
  void DoFit (Rho &rho);
  void WriteFit(IOSectionClass &outSection);
  void ReadInput (IOSectionClass &inSection);

  bool Read  (IOSectionClass &inSection, double lowestBeta,
	      int NumBetas);
  double U(double q, double z, double s2, int level);
  double dU(double q, double z, double s2, int level);

  /////////////////////////
  /// Long-ranged stuff ///
  /////////////////////////
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
  /// The potential to which this action corresponds.
  double V  (double r);
  /// The q-derivative of the above
  double Vp (double r);
  /// The q-derivative of the above
  double Vpp(double r);
  /// Set the cutoff radius for the long-range breakup
  void Setrc (double rc);
  /// Modified fourier transform.
  double Xk_U     (double k, int level);
  double dXk_U_dk (double k, int level);
  double Xk_dU    (double k, int level);
  double Xk_V     (double k);

  bool IsLongRange();
  //  void DoBreakup(const Vec3 &box, const Array<Vec3,1> &kVecs);

  PADipoleFitClass()
  { 

  }
};

#endif
