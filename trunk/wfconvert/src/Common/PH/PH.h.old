#ifndef PH_H
#define PH_H

#include "../Integration/Integrate.h"
#include "../Splines/Grid.h"
#include "Chebyshev.h"
#include "../Splines/CubicSpline.h"

enum PHType {PH_NONE, PH_CHEBYSHEV, PH_CUBIC, PH_CUBICXC, PH_NUCLEAR};

class FullCorePotential
{
private:
  int GridInitialized;
public:
  char FileName[1000];

  CubicSpline V;
  scalar Z;

  inline scalar operator()(scalar r)
  {
    if (r < V.grid->End)
      return (-Z/r + V(r));
    else
      return (-Z/r);
  }
  inline scalar Deriv(scalar r)
  {
    if (r < V.grid->End)
      return (Z/(r*r) + V.Deriv(r));
    else      
    return (Z/(r*r));
  }
  void Read(char *FName);

  FullCorePotential()
  {
    GridInitialized = 0;
  }
  ~FullCorePotential()
  {
    if (GridInitialized)
      delete V.grid;
  }
};




class PseudoHamiltonian
{
public:
  scalar CoreRadius;
  scalar Z, Zion;  // Zion is the unscreened pseudohamiltonian charge
  double Amin, Amax, Bmin, Bmax, Vmin, Vmax;

  FullCorePotential *FullCoreV;
  CubicSpline VHXC;  // The Hartree, exchange, and correlation potentials
  scalar NumElecs;
  int UseVHXC;

  inline void UpdateVHXC (CubicSpline NewVHXC, scalar Nelecs)
  {
    NumElecs = Nelecs;
    VHXC = NewVHXC;
    UseVHXC = 1;
  }

  virtual PHType Type()
  {
    cerr << "PH base class has not type.\n";
    exit(1);
    return (PH_NONE);
  }

  virtual int NumParams()
  {
    cerr << "Should never get to PseudoHamiltonian base class.\n";
    exit(1);
    return(0);
  }

  virtual scalar &Params(int i)
  {
    cerr << "Should never get to PseudoHamiltonian base class.\n";
    exit(1);
    static scalar r = 0.0;
    return(r);
  }
  
  virtual scalar Params(int i) const
  {
    cerr << "Should never get to PseudoHamiltonian base class.\n";
    exit(1);
    return(0.0);
  }
  virtual void ABV(scalar r, scalar &A, scalar &B, scalar &V,
		   scalar &dAdr)
  {
    cerr << "Should never get to PseudoHamiltonian base class.\n";
  }

  virtual scalar d2Adr2(scalar r)
  {
    cerr << "Should never get to PseudoHamiltonian base class.\n";
    return (0.0);
  }

  virtual scalar dVdr(scalar r)
  {
    cerr << "Should never get to PseudoHamiltonian base class.\n";
    exit(1);
    return(0.0);
  }
  virtual scalar    V(scalar r)
  {
    cerr << "Should never get to PseudoHamiltonian base class.\n";
    exit(1);
    return(0.0);
  }
  virtual void Write(char *FileName)
  {
    FILE *fout;
    if ((fout = fopen (FileName, "w")) == NULL)
      {
	cerr << "Can't open " << FileName << " for writing.  Exitting.\n";
	exit(1);
      }
    Write (fout);
  }
  virtual void Write(FILE *fout)
  {
    cerr << "Should never get to PseudoHamiltonian base class.\n";
  }
  virtual void Read (FILE *fin)
  {
    cerr << "Should never get to PseudoHamiltonian base class.\n";
  }

  PseudoHamiltonian()
  {
    UseVHXC = 0;
  }
};


///////////////////////////////////////////////////////////////////
//                        Chebyshev PH                           //
///////////////////////////////////////////////////////////////////


// This class will be used to store the A potential for the
// pseudoHamiltonian.  We have two boundary conditions, which puts
// restrictions on the coefficients of the expansion.  A(r_CA) = 1.0
// and A'(r_CA) = 0.
class AFunction
{
private:
  scalar rA;
  Array<scalar,1> ChebyCoefs;
  void Update();
public:
  Array<scalar,1> Parameters;
  int UpToDate;
  int NumParams;
  scalar AddConst;
  scalar &Params(int i); 
  scalar Params(int i) const;
  scalar operator()(scalar r);
  scalar Deriv(scalar r);
  scalar Deriv2(scalar r);
  void Init(Array<scalar,1> params, scalar rmax);
  AFunction(Array<scalar,1> params, scalar rmax);
  AFunction()
  {
    AddConst = 0.25;
    // Do nothing
  }
};



// This class will be used to store the B potential for the
// pseudoHamiltonian.  We have one boundary condition (for now), which
// puts restrictions on the coefficients of the expansion,  
// B(r_CB) = 1.0.

class BFunction
{
private:
  scalar rB;
  Array<scalar,1> ChebyCoefs;
  void Update();
  int UpToDate;

  AFunction *Afunc;
public:
  Array<scalar,1> Parameters;
  scalar AddConst;
  int NumParams;
  scalar &Params(int i); 
  scalar Params(int i) const;
  scalar operator()(scalar r);
  void Init(Array<scalar,1> params, scalar rmax, AFunction *A);
  BFunction(Array<scalar,1> params, scalar rmax, AFunction *A);
  BFunction()
  {
    AddConst = 0.25;
    // Do nothing
  }
};




class VFunction
{
private:
  scalar V_at_rV, dV_at_rV;
  Array<scalar,1> ChebyCoefs;
  void Update();
  int UpToDate;
public:
  Array<scalar,1> Parameters;
  scalar rV;
  int NumParams;
  scalar &Params(int i);
  scalar Params(int i) const;
  scalar operator()(scalar r);
  scalar Deriv(scalar r);
  void Init(Array<scalar,1> params, 
	    scalar rmax, scalar Vrmax, scalar dVrmax);
  VFunction(Array<scalar,1> params, scalar rmax, 
	    scalar Vrmax, scalar dVrmax);
  VFunction()
  {
    // Do nothing
  }
};


class PH_Chebyshev : public PseudoHamiltonian
{
private:
  int LocalAlloc;

public:
  AFunction *Afunc;
  BFunction *Bfunc;
  VFunction *PseudoCoreV;

  PHType Type()
  {
    return PH_CHEBYSHEV;
  }
  
  int NumParams()
  {
    return (Afunc->NumParams+Bfunc->NumParams+PseudoCoreV->NumParams);
  }

  scalar & Params(int i)
  {
    if (i < Afunc->NumParams)
      return (Afunc->Params(i));
    else if (i < (Afunc->NumParams+Bfunc->NumParams))
      return (Bfunc->Params(i-Afunc->NumParams));
    else
      return (PseudoCoreV->Params(i-Afunc->NumParams -Bfunc->NumParams));
  }

  scalar Params(int i) const
  {
    if (i < Afunc->NumParams)
      return (Afunc->Params(i));
    else if (i < (Afunc->NumParams+Bfunc->NumParams))
      return (Bfunc->Params(i-Afunc->NumParams));
    else
      return (PseudoCoreV->Params(i-Afunc->NumParams -Bfunc->NumParams));
  }

  scalar V(scalar r)
  {
    if (r<CoreRadius)
      return ((*PseudoCoreV)(r));
    else
      return (FullCoreV->Deriv(r));
  }



  scalar dVdr(scalar r)
  {
    if (r<CoreRadius)
      return (PseudoCoreV->Deriv(r));
    else
      return (FullCoreV->Deriv(r));
  }

  
  scalar d2Adr2(scalar r)
  {
    if (r<CoreRadius)
      return (Afunc->Deriv2(r));
    else
      return (0.0);
  }
 
  inline void Init(AFunction *A,
		   BFunction *B,
		   VFunction *V,
		   FullCorePotential *FullCorePot)
  {
    Afunc = A; Bfunc = B; PseudoCoreV = V; FullCoreV = FullCorePot;
    LocalAlloc = 0;
  }

  PH_Chebyshev()
  {
    LocalAlloc = 0;
  }

  ~PH_Chebyshev()
  {
    if (LocalAlloc)
      {
	delete (Afunc);
	delete (Bfunc);
	delete (PseudoCoreV);
	delete (FullCoreV);
      }
  }

  void Write(FILE *fout);
  void Write(char *FileName);
  void Read (FILE *fin);
};



///////////////////////////////////////////////////////////////////
//                       Cubic Spline PH                         //
///////////////////////////////////////////////////////////////////

class PH_CubicSpline : public PseudoHamiltonian
{
private:
  int UpDateB;
public:
  CubicSpline PA, PB, Vfunc;
  LinearGrid Agrid, Bgrid, Vgrid;
  
  PHType Type()
  {
    return PH_CUBIC;
  }

  int NumParams()
  {
    return (PA.NumParams-1 + PB.NumParams-2 + Vfunc.NumParams-1);
  }
  
  scalar &Params(int i)
  {
    if (i==0)
      UpDateB = 1;
    if (i < (PA.NumParams-1))
      return (PA.Params(i));
    else if (i < (PA.NumParams-1 + PB.NumParams-2))
      return (PB.Params(i-PA.NumParams + 2));
    else 
      return (Vfunc.Params(i-PA.NumParams - PB.NumParams + 3));
  }

  scalar Params(int i) const
  {
    if (i < (PA.NumParams-1))
      return (PA.Params(i));
    else if (i < (PA.NumParams-1 + PB.NumParams-2))
      return (PB.Params(i-PA.NumParams + 2));
    else 
      return (Vfunc.Params(i-PA.NumParams - PB.NumParams + 3));
  }

	
  void ABV(scalar r, scalar &A, scalar &B, scalar &V,
		  scalar &dAdr)
  {
    if (r<CoreRadius)
      {
	if (UpDateB)
	  {
	    PB.Params(0) = sqrt(Amin - Bmin + PA.Params(0)*PA.Params(0));
	    UpDateB = 0;
	  }
	scalar pa = PA(r);
	scalar pb = PB(r);
	A = Amin + pa * pa;
	B = Bmin + pb * pb;
	V = Vfunc(r);
	dAdr = 2.0*pa * PA.Deriv(r);
      }
    else
      {
	A = 1.0;
	B = 1.0;
	V = (*FullCoreV)(r);
	dAdr = 0.0;
      }
  }

  scalar V(scalar r)
  {
    if (r<CoreRadius)
      return (Vfunc(r));
    else
      return ((*FullCoreV)(r));
  }

  scalar dVdr(scalar r)
  {
    if (r<CoreRadius)
      return (Vfunc.Deriv(r));
    else
      return (FullCoreV->Deriv(r));
  }

  scalar d2Adr2(scalar r)
  {
    if (r<CoreRadius)
      {
	scalar p = PA(r);
	scalar pp = PA.Deriv(r);
	scalar ppp = PA.Deriv2(r);
	return (2.0*(p*ppp+pp*pp));
      }
    else
      return (0.0);
  }

  void SetCoreRadius (scalar radius)
    {
      CoreRadius = radius;
      Agrid.Init(0.0, CoreRadius, PA.NumParams);
      Bgrid.Init(0.0, CoreRadius, PB.NumParams);
      int NumV = Vfunc.NumParams;
      Vgrid.Init(0.0, CoreRadius, NumV);
      Array<scalar,1> Vparams(NumV);
      for (int i=0; i<NumV; i++)
	Vparams(i) = Vfunc.Params(i);
      Vparams(NumV-1) = (*FullCoreV)(CoreRadius);
      Vfunc.Init (&Vgrid, Vparams, 5.0e30, FullCoreV->Deriv(CoreRadius));
    }
      

  inline void Init (scalar amin, scalar bmin, scalar vmin,
		    scalar amax, scalar bmax, scalar vmax,
		    Array<scalar,1> Aparams,
		    Array<scalar,1> Bparams,
		    Array<scalar,1> Vparams,
		    FullCorePotential *FullCorePot,
		    scalar coreradius)
  {
    CoreRadius = coreradius;
    Amin = amin; Amax = amax; 
    Bmin = bmin; Bmax = bmax;
    Vmin = vmin; Vmax = vmax;
    FullCoreV = FullCorePot;

    int NumA = Aparams.rows();
    Agrid.Init(0.0, CoreRadius, NumA+1);
    Array<scalar,1> Ainit(NumA+1);
    Ainit(Range(0,NumA-1)) = Aparams;
    Ainit(NumA) = sqrt(1.0-amin);
    PA.Init(&Agrid, Ainit, 0.0, 0.0);

    int NumB = Bparams.rows();
    Bgrid.Init(0.0, CoreRadius, NumB+2);
    Array<scalar,1> Binit(NumB+2);
    Binit(Range(1,NumB)) = Bparams;
    Binit(0) = sqrt(Aparams(0)*Aparams(0) + Amin - Bmin);    
    Binit(NumB+1) = sqrt(1.0 - bmin);
    PB.Init(&Bgrid, Binit, 0.0, 0.0);

    int NumV = Vparams.rows();
    Vgrid.Init(0.0, CoreRadius, NumV+1);
    Array<scalar,1> Vinit(NumV+1);
    Vinit(Range(0,NumV-1)) = Vparams;
    Vinit(NumV) = (*FullCoreV)(CoreRadius);
    Vfunc.Init(&Vgrid, Vinit, 5.0e30, FullCoreV->Deriv(CoreRadius));
  }
    
  PH_CubicSpline()
  {
    // Do nothing for now
  }

  void Write(FILE *fout);
  void Write(char *FileName);
  void Read (FILE *fin);
  
};




///////////////////////////////////////////////////////////////////
// Cubic Spline PH with Exchange-Correlation subtracted from V   //
// This is used to store a PH which has had the valance exchange //
// correlation subtracted off.  There is a dense grid for V,     //
// while the A anb B grids remain the same.                      //
///////////////////////////////////////////////////////////////////
class PH_CubicSplineXC : public PseudoHamiltonian
{
private:
  int UpDateB;
public:
  CubicSpline PA, PB, Vfunc;
  Grid *Agrid, *Bgrid, *Vgrid;

  // Nelecs is the number of electrons in the pseudo-atom;
  scalar NumElecs;
  scalar Amin, Bmin;

  PHType Type()
  {
    return PH_CUBICXC;
  }
  
  int NumParams()
  {
    return (PA.NumParams-1 + PB.NumParams-2 + Vfunc.NumParams-1);
  }
  
  scalar &Params(int i)
  {
    if (i==0)
      UpDateB = 1;
    if (i < (PA.NumParams-1))
      return (PA.Params(i));
    else if (i < (PA.NumParams-1 + PB.NumParams-2))
      return (PB.Params(i-PA.NumParams + 2));
    else 
      return (Vfunc.Params(i-PA.NumParams - PB.NumParams + 3));
  }

  scalar Params(int i) const
  {
    if (i < (PA.NumParams-1))
      return (PA.Params(i));
    else if (i < (PA.NumParams-1 + PB.NumParams-2))
      return (PB.Params(i-PA.NumParams + 2));
    else 
      return (Vfunc.Params(i-PA.NumParams - PB.NumParams + 3));
  }

	
  void ABV(scalar r, scalar &A, scalar &B, scalar &V,
		  scalar &dAdr)
  {
    if (r<CoreRadius)
      {
	if (UpDateB)
	  {
	    PB.Params(0) = sqrt(Amin - Bmin + PA.Params(0)*PA.Params(0));
	    UpDateB = 0;
	  }
	scalar pa = PA(r);
	scalar pb = PB(r);
	A = Amin + pa * pa;
	B = Bmin + pb * pb;

	dAdr = 2.0*pa * PA.Deriv(r);
      }
    else
      {
	A = 1.0;
	B = 1.0;
	dAdr = 0.0;
      }
    if (r < Vgrid->End)
      {
	V = Vfunc(r);
	if (UseVHXC)
	  if (r < VHXC.grid->End)
	    V += VHXC(r);
      }
    else
      V = -(Zion - NumElecs) / r;

  }

  scalar V(scalar r)
  {
    scalar v;
    if (r < Vgrid->End)
      {
	v = Vfunc(r);
	if (UseVHXC)
	  if (r < VHXC.grid->End)
	    v += VHXC(r);
      }
    else
      v = -(Zion-NumElecs) / r;
    return (v);
  }

  scalar dVdr(scalar r)
  {
    scalar dV;
    if (r < Vgrid->End)
      {
	dV = Vfunc.Deriv(r);
	if (UseVHXC)
	  if (r < VHXC.grid->End)
	    dV += VHXC.Deriv(r);
      }
    else
      dV = (Zion-NumElecs) / (r*r);

    return (dV);
  }

  
   scalar d2Adr2(scalar r)
  {
    if (r<CoreRadius)
      {
	scalar p = PA(r);
	scalar pp = PA.Deriv(r);
	scalar ppp = PA.Deriv2(r);
	return (2.0*(p*ppp+pp*pp));
      }
    else
      return (0.0);
  }   

  inline void Init (Array<scalar,1> Aparams,
		    Array<scalar,1> Bparams,
		    Array<scalar,1> Vparams,
		    Grid *agrid, Grid *bgrid,
		    Grid *vgrid, scalar zion,
		    scalar coreradius)
  {
    UseVHXC = 0;
    NumElecs = 0;
    CoreRadius = coreradius;
    Agrid = agrid;
    Bgrid = bgrid;
    Vgrid = vgrid;
    Zion = zion;
    Amin = 0.0; 
    Bmin = 0.0;

    int NumA = Aparams.rows();
    Array<scalar,1> Ainit(NumA+1);
    if (Ainit.rows() != Agrid->NumPoints)
      {
	cerr << "Size mismatch in PH_CubicSplineXC.\n";
	exit(1);
      }
    Ainit(Range(0,NumA-1)) = Aparams;
    Ainit(NumA) = sqrt(1.0-Amin);
    PA.Init(Agrid, Ainit, 0.0, 0.0);

    int NumB = Bparams.rows();
    Array<scalar,1> Binit(NumB+2);
    if (Binit.rows() != Bgrid->NumPoints)
      {
	cerr << "Size mismatch in PH_CubicSplineXC.\n";
	exit(1);
      }
    Binit(Range(1,NumB)) = Bparams;
    Binit(0) = sqrt(Aparams(0)*Aparams(0) + Amin - Bmin);    
    Binit(NumB+1) = sqrt(1.0 - Bmin);
    PB.Init(Bgrid, Binit, 0.0, 0.0);

    int NumV = Vparams.rows();
    Array<scalar,1> Vinit(NumV+1);
    if (Vinit.rows() != Vgrid->NumPoints)
      {
	cerr << "Size mismatch in PH_CubicSplineXC.\n";
	cerr << "Vinit.rows() = " << Vinit.rows() << "\n";
	cerr << "Vgrid.NumPoints = " << Vgrid->NumPoints << "\n";
	exit(1);
      }
    Vinit(Range(0,NumV-1)) = Vparams;
    scalar rmax = Vgrid->End; 

    // Force the last point to agree with asymptotic form
    Vinit(NumV) = -Zion/rmax;
    // Make the slope match assymptotic form as well
    Vfunc.Init(Vgrid, Vinit, 5.0e30, Zion/(rmax*rmax));
  }
    
  PH_CubicSplineXC()
  {
    // Do nothing for now
  }

  void Write(FILE *fout);
  void Write(char *FileName);
  void Read (FILE *fin);
  
};




///////////////////////////////////////////////////////////////////
// Nuclear PH:  has no parameters.  The potential is dependent   //
// only on the nuclear charge Z, and possibly on the exchange-   //
// correlation potential.  This is used for full-core atom       //
// calculations.                                                 //
///////////////////////////////////////////////////////////////////
class PH_Nuclear : public PseudoHamiltonian
{
public:
  scalar dummy;  // To satisfy compiler warnings about no return value

  PHType Type()
  {
    return PH_NUCLEAR;
  }
  
  int NumParams()
  {
    return (0);
  }
  
  scalar &Params(int i)
  {
    cerr << "Nuclear PH has no parameters.\n";
    return (dummy);
  }

  scalar Params(int i) const
  {
    cerr << "Nuclear PH has no parameters.\n";
    return (0.0);
  }

	
  void ABV(scalar r, scalar &A, scalar &B, scalar &V,
		  scalar &dAdr)
  {
    A = B = 1.0;
    dAdr = 0.0;

    V = -Zion/r;
    if (UseVHXC)
      if (r < VHXC.grid->End)
	V += VHXC(r);
      else
	V = -(Zion-NumElecs)/r;
  }

  scalar V(scalar r)
  {
    scalar v = -Zion/r;
    if (UseVHXC)
      if (r < VHXC.grid->End)
	v += VHXC(r);
      else
	v = -(Zion-NumElecs) / r;
    return (v);
  }

  scalar dVdr(scalar r)
  {
    scalar dV;
    dV = Zion / (r*r);
    if (UseVHXC)
      if (r < VHXC.grid->End)
	dV += VHXC.Deriv(r);
      else
	dV = (Zion-NumElecs) / (r*r);

    return (dV);
  }

  scalar d2Adr2(scalar r)
  {
    return (0.0);
  }

  PH_Nuclear(scalar z)
  {

    Zion = z;
  }

  PH_Nuclear()
  {
    // Do nothing for now
  }

  //void Write(char *FileName);
  //void Read (FILE *fin);
  
};



///////////////////////////////////////////////////////////////////
// FullCore PH:  has no parameters.  It is use solely to         //
// calculate the full core radial WFs for comparison with the    //
// pseudo radial WFs.                                            //
///////////////////////////////////////////////////////////////////
class PH_FullCore : public PseudoHamiltonian
{
public:
  // Zion is the unscreened pseudohamiltonian charge and Nelecs
  // is the number of electrons in the pseudo-atom;
  scalar dummy;  // To satisfy compiler warning.


  PHType Type()
  {
    return PH_NUCLEAR;
  }
  
  int NumParams()
  {
    return (0);
  }
  
  scalar &Params(int i)
  {
    cerr << "Nuclear PH has no parameters.\n";
    return (dummy);
  }

  scalar Params(int i) const
  {
    cerr << "Nuclear PH has no parameters.\n";
    return (0.0);
  }

	
  void ABV(scalar r, scalar &A, scalar &B, scalar &V,
		  scalar &dAdr)
  {
    A = B = 1.0;
    dAdr = 0.0;

    V = (*FullCoreV)(r);
  }

  scalar V(scalar r)
  {
    return ((*FullCoreV)(r));
  }

  scalar dVdr(scalar r)
  {
    return (FullCoreV->Deriv(r));
  }

  scalar d2Adr2(scalar r)
  {
    return (0.0);
  }
  
  
  PH_FullCore(FullCorePotential *FullCorePot)
  {
    Z = FullCorePot->Z;
    FullCoreV = FullCorePot;
  }

  void Write(FILE *fout);
  void Write(char *FileName);
  void Read (FILE *fin);
  
};




///////////////////////////////////////////////////////////////////
// Zero PH:  has no parameters.  Has unit masses and V(r) = 0.   //
// It is used for testing in matrix squaring.                    //
///////////////////////////////////////////////////////////////////
class PH_Zero : public PseudoHamiltonian
{
public:
  scalar dummy;  // To satisfy compiler warnings about no return value

  PHType Type()
  {
    return PH_NUCLEAR;
  }
  
  int NumParams()
  {
    return (0);
  }
  
  scalar &Params(int i)
  {
    cerr << "Nuclear PH has no parameters.\n";
    return (dummy);
  }

  scalar Params(int i) const
  {
    cerr << "Nuclear PH has no parameters.\n";
    return (0.0);
  }

	
  void ABV(scalar r, scalar &A, scalar &B, scalar &V,
		  scalar &dAdr)
  {
    A = B = 1.0;
    dAdr = 0.0;

    V = 0.0;
  }

  scalar V(scalar r)
  {
    return (0.0);
  }

  scalar dVdr(scalar r)
  {
    return (0.0);
  }

  scalar d2Adr2(scalar r)
  {
    return (0.0);
  }

  PH_Zero()
  {
    // Do nothing for now
  }

  //void Write(char *FileName);
  //void Read (FILE *fin);
  
};



///////////////////////////////////////////////////////////////////
// Gaussian PH:  has two parameters.  Has unit masses and        //
// V(r) = Amp * exp(-r^2/(2.0*sigma^2)).                         //
// It is used for testing in matrix squaring.                    //
///////////////////////////////////////////////////////////////////
class PH_Gaussian : public PseudoHamiltonian
{
public:
  scalar Amp, sigma;  // To satisfy compiler warnings about no return value

  PHType Type()
  {
    return PH_NUCLEAR;
  }
  
  int NumParams()
  {
    return (0);
  }
  
  scalar &Params(int i)
  {
    if (i==0)
      return (Amp);
    else
      return (sigma);
  }

  scalar Params(int i) const
  {
    if (i==0)
      return (Amp);
    else
      return (sigma);
  }

	
  void ABV(scalar r, scalar &A, scalar &B, scalar &V,
		  scalar &dAdr)
  {
    A = B = 1.0;
    dAdr = 0.0;

    V = Amp*exp(-r*r/(2.0*sigma*sigma));
  }

  scalar V(scalar r)
  {
    return (Amp*exp(-r*r/(2.0*sigma*sigma)));
  }

  scalar dVdr(scalar r)
  {
    return (-r/(sigma*sigma) * Amp*exp(-r*r/(2.0*sigma*sigma)));
  }

  scalar d2Adr2(scalar r)
  {
    return (0.0);
  }

  PH_Gaussian(double Amp_, double sigma_)
  {
    Amp = Amp_;
    sigma = sigma_;
    // Do nothing for now
  }

  //void Write(char *FileName);
  //void Read (FILE *fin);
  
};



PseudoHamiltonian *ReadPH (InputBuffer &SectionBuf);
PseudoHamiltonian *Read_PH(char *FileName);

#endif
