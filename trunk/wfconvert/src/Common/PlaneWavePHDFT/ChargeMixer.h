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

#ifndef CHARGE_MIXER_H
#define CHARGE_MIXER_H
#include "FFTBox.h"
#include <vector>

class ChargeMixerClass
{
protected:
  FFTBox &FFT;
public:
  virtual void Reset() = 0;
  virtual void Mix (const Array<double,3> &newCharge, 
		    Array<double,3>       &mixedCharge_r,
		    zVec                  &mixedCharge_G) = 0;
  ChargeMixerClass (FFTBox &fft) :
    FFT(fft)
  {
    // nothing else for now
  }
};

class KerkerMixerClass : public ChargeMixerClass
{
protected:
  bool HaveLastCharge;
  /// This parameter controls the shift point of the Pade
  /// approximant.  
  double Lambda;
  double NewFraction;
  zVec LastCharge, NewCharge;
public:
  void SetLambda(double lambda);
  double GetLambda();
  void SetNewFraction (double frac);
  void Reset();
  void Mix (const Array<double,3> &newCharge,
   	    Array<double,3>       &mixedCharge_r,
	    zVec                  &mixedCharge_G);
  KerkerMixerClass(FFTBox &fft);
};


class BroydenMixerClass : public ChargeMixerClass
{
protected:
  Array<double,3> LastCharge;
  bool HaveLastCharge;
  std::vector<zVec> DeltaCharge;
public:
  void Reset();
  void Mix (const Array<double,3> &newCharge,
   	    Array<double,3>       &mixedCharge_r,
	    zVec                  &mixedCharge_G);
  BroydenMixerClass(FFTBox &fft);
};

#endif
