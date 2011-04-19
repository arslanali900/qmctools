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

#ifndef MINIMIZE_H
#define MINIMIZE_H

#include "../Blitz.h"

class MinimizeFunction
{
public:
  scalar dummy;  // for satisfying the compiler.
  virtual int NumParams()
  {
    cerr << "MinimizeFunction base class.  Should never get here.\n";
    exit(1);
    return (0);
  }
  virtual scalar &Params(int i)
  {
    cerr << "MinimizeFunction base class.  Should never get here.\n";
    exit(1);
    return (dummy);
  }
  virtual scalar Params(int i) const
  {
    cerr << "MinimizeFunction base class.  Should never get here.\n";
    exit(1);
    return (0.0);
  }
  virtual scalar Cost()
  {
    cerr << "MinimizeFunction base class.  Should never get here.\n";
    exit(1);
    return (0.0);
  }
  virtual void WriteStuff()
  {
    cerr << "MinimizeFunction base class.  Should never get here.\n";
    exit(1);
  }

  //void Minimize();
  void ReadParameters(char *FileName);
  void WriteParameters(char *FileName);
};



class Minimizer
{
public:
  virtual void Minimize(MinimizeFunction &MinFunc)
  {
    cerr << "Minimizer base class: should never get here.\n";
    exit(1);
  }
};


class ConjugateGradient : public Minimizer
{
public:
  scalar epsilon;
  scalar Tolerance;
  scalar StepSize;
  MinimizeFunction *MinFunc;
  void Minimize (MinimizeFunction &MinimFunc);
};


class AnnealingSchedule
{
public:
  scalar StartTemp, EndTemp;
  int NumTemps;
  int StepsPerTemp;
  virtual scalar Temp(int TempNum)
  {
    cerr << "AnnealSchedule base class:  Should never get here.\n";
    exit(1);
    return (0.0);
  }
};

class ExponentialSchedule : public AnnealingSchedule
{
public:
  scalar Chi;
  scalar Temp(int TempNum);
};


class LinearSchedule : public AnnealingSchedule
{
public:
  scalar Temp(int TempNum);
};


class VanderbiltAnnealer : public Minimizer
{
public:
  AnnealingSchedule *Schedule;
  // The Q matrix gives the step for a random uniform vector, u.
  // \Delta x = Q*u
  Array<scalar,2> Q;
  // xMean holds the mean position for the last M steps;
  Array<scalar,1> xMean;
  // Covariance holds the covariance matrix for the last block of steps.
  Array<scalar,2> Covariance;
  // Holds the Minimum vector found
  Array<scalar,1> MinParams;
  // Holds the minimum cost found
  scalar MinCost;

  scalar GrowthFactor;
  scalar AcceptRatio;
  scalar MeanCost;
  Array<scalar,2> s;
  scalar kT;
  scalar HeatCapacity;
  MinimizeFunction *MinFunc;


  void Metropolis(int NumSteps);

  void Minimize (MinimizeFunction &MinimFunc);
  VanderbiltAnnealer(AnnealingSchedule &ASchedule)
  {
    
    Schedule = &ASchedule;
  }
};



#endif
