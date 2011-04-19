#ifndef PH_H
#define PH_H

#include "../Integration/Integrate.h"
#include "../Splines/CubicSpline.h"
#include "../IO/InputOutput.h"

enum PHType {PH_NONE, PH_CHEBYSHEV, PH_CUBIC, PH_CUBICXC, PH_NUCLEAR,
PH_COULOMB3D, PH_COULOMB};

class FullCorePotential
{
private:
  int GridInitialized;
public:
  char FileName[1000];

  CubicSpline V;
  double Z;

  inline double operator()(double r)
  {
    if (r < V.grid->End)
      return (-Z/r + V(r));
    else
      return (-Z/r);
  }
  inline double Deriv(double r)
  {
    if (r < V.grid->End)
      return (Z/(r*r) + V.Deriv(r));
    else      
      return (Z/(r*r));
  }

  void Write (IOSectionClass &outSection)
  {
    outSection.WriteVar ("Z", Z);
    outSection.NewSection("V");
    V.Write(outSection);
    outSection.CloseSection();
  }
  void Read  (IOSectionClass &inSection)
  {
    assert(inSection.ReadVar("Z", Z));
    assert(inSection.OpenSection("V"));
    V.Read (inSection);
    inSection.CloseSection();
  }

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
  double CoreRadius;
  double Z, Zion;  // Zion is the unscreened pseudohamiltonian charge
  double Amin, Amax, Bmin, Bmax, Vmin, Vmax;

  FullCorePotential *FullCoreV;
  CubicSpline VHXC;  // The Hartree, exchange, and correlation potentials
  double NumElecs;
  bool UseVHXC;

  inline void UpdateVHXC (CubicSpline NewVHXC, double Nelecs)
  {
    NumElecs = Nelecs;
    VHXC = NewVHXC;
    UseVHXC = true;
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

  virtual double &Params(int i)
  {
    cerr << "Should never get to PseudoHamiltonian base class.\n";
    exit(1);
    static double r = 0.0;
    return(r);
  }
  
  virtual double Params(int i) const
  {
    cerr << "Should never get to PseudoHamiltonian base class.\n";
    exit(1);
    return(0.0);
  }
  virtual double    V(double r)
  {
    cerr << "Should never get to PseudoHamiltonian base class.\n";
    exit(1);
    return(0.0);
  }
  virtual void ABV(double r, double &A, double &B, double &Vval,
		   double &dAdr)
  {
    A = 1.0; B = 1.0; dAdr = 0.0;
    Vval = V(r);
  }

  virtual double d2Adr2(double r)
  {
    cerr << "Should never get to PseudoHamiltonian base class.\n";
    return (0.0);
  }

  virtual double dVdr(double r)
  {
    cerr << "Should never get to PseudoHamiltonian base class.\n";
    exit(1);
    return(0.0);
  }

  virtual void Write (IOSectionClass &outSection)
  {
    outSection.WriteVar ("CoreRadius", CoreRadius);
    outSection.WriteVar ("Z", Z);
    outSection.WriteVar ("Zion", Zion);
    outSection.WriteVar ("NumElecs", NumElecs);
    outSection.WriteVar ("Amin", Amin);
    outSection.WriteVar ("Amax", Amax);
    outSection.WriteVar ("Bmin", Bmin);
    outSection.WriteVar ("Bmax", Bmax);
    outSection.WriteVar ("Vmin", Vmin);
    outSection.WriteVar ("Vmax", Vmax);
    outSection.WriteVar ("UseVHXC", UseVHXC);
    if (UseVHXC) {
      outSection.NewSection("VHXC");
      VHXC.Write (outSection);
      outSection.CloseSection();
    }
    outSection.NewSection ("FullCoreV");
    FullCoreV->Write(outSection);
    outSection.CloseSection();
  }

  virtual bool Read (IOSectionClass &inSection)
  {
    assert (inSection.ReadVar ("CoreRadius", CoreRadius));
    assert (inSection.ReadVar ("Z", Z));
    assert (inSection.ReadVar ("Zion", Zion));
    assert (inSection.ReadVar ("NumElecs", NumElecs));
    assert (inSection.ReadVar ("Amin", Amin));
    assert (inSection.ReadVar ("Amax", Amax));
    assert (inSection.ReadVar ("Bmin", Bmin));
    assert (inSection.ReadVar ("Bmax", Bmax));
    assert (inSection.ReadVar ("Vmin", Vmin));
    assert (inSection.ReadVar ("Vmax", Vmax));
    assert (inSection.ReadVar ("UseVHXC", UseVHXC));
    if (UseVHXC){
      assert(inSection.OpenSection("VHXC"));
      VHXC.Read (inSection);
      inSection.CloseSection();
    }
    assert(inSection.OpenSection ("FullCoreV"));
    FullCoreV = new FullCorePotential;
    FullCoreV->Read (inSection);
    inSection.CloseSection();
    return (true);
  } 



  PseudoHamiltonian()
  {
    UseVHXC = false;
  }
};


class Potential : public PseudoHamiltonian
{
public:
  // These five functions must be specialized
  virtual PHType Type() = 0;
  virtual double V(double r) = 0;
  virtual double dVdr(double r) = 0;
  virtual void Write(IOSectionClass &outSection) = 0;
  virtual bool Read (IOSectionClass &inSection)  = 0;

  int NumParams()
  { return 0; }
  double Params(int i) const 
  { return 0.0; }
  double &Params(int i) 
  { return CoreRadius; }
  
  void ABV(double r, double &A, double &B, double &Vval, double &dAdr)
  {
    A = 1.0; B = 1.0; dAdr = 0.0;
    Vval = V(r);
  }
  double d2Adr2 (double r)
  { return 0.0; }
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
  
  double &Params(int i)
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

  double Params(int i) const
  {
    if (i < (PA.NumParams-1))
      return (PA.Params(i));
    else if (i < (PA.NumParams-1 + PB.NumParams-2))
      return (PB.Params(i-PA.NumParams + 2));
    else 
      return (Vfunc.Params(i-PA.NumParams - PB.NumParams + 3));
  }

	
  void ABV(double r, double &A, double &B, double &V,
		  double &dAdr)
  {
    if (r<CoreRadius)
      {
	if (UpDateB)
	  {
	    PB.Params(0) = sqrt(Amin - Bmin + PA.Params(0)*PA.Params(0));
	    UpDateB = 0;
	  }
	double pa = PA(r);
	double pb = PB(r);
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

  double V(double r)
  {
    if (r<CoreRadius)
      return (Vfunc(r));
    else
      return ((*FullCoreV)(r));
  }

  double dVdr(double r)
  {
    if (r<CoreRadius)
      return (Vfunc.Deriv(r));
    else
      return (FullCoreV->Deriv(r));
  }

  double d2Adr2(double r)
  {
    if (r<CoreRadius)
      {
	double p = PA(r);
	double pp = PA.Deriv(r);
	double ppp = PA.Deriv2(r);
	return (2.0*(p*ppp+pp*pp));
      }
    else
      return (0.0);
  }

  void SetCoreRadius (double radius)
    {
      CoreRadius = radius;
      Agrid.Init(0.0, CoreRadius, PA.NumParams);
      Bgrid.Init(0.0, CoreRadius, PB.NumParams);
      int NumV = Vfunc.NumParams;
      Vgrid.Init(0.0, CoreRadius, NumV);
      Array<double,1> Vparams(NumV);
      for (int i=0; i<NumV; i++)
	Vparams(i) = Vfunc.Params(i);
      Vparams(NumV-1) = (*FullCoreV)(CoreRadius);
      Vfunc.Init (&Vgrid, Vparams, 5.0e30, FullCoreV->Deriv(CoreRadius));
    }
      

  inline void Init (double amin, double bmin, double vmin,
		    double amax, double bmax, double vmax,
		    Array<double,1> Aparams,
		    Array<double,1> Bparams,
		    Array<double,1> Vparams,
		    FullCorePotential *FullCorePot,
		    double coreradius)
  {
    CoreRadius = coreradius;
    Amin = amin; Amax = amax; 
    Bmin = bmin; Bmax = bmax;
    Vmin = vmin; Vmax = vmax;
    FullCoreV = FullCorePot;

    int NumA = Aparams.rows();
    Agrid.Init(0.0, CoreRadius, NumA+1);
    Array<double,1> Ainit(NumA+1);
    Ainit(Range(0,NumA-1)) = Aparams;
    Ainit(NumA) = sqrt(1.0-amin);
    PA.Init(&Agrid, Ainit, 0.0, 0.0);

    int NumB = Bparams.rows();
    Bgrid.Init(0.0, CoreRadius, NumB+2);
    Array<double,1> Binit(NumB+2);
    Binit(Range(1,NumB)) = Bparams;
    Binit(0) = sqrt(Aparams(0)*Aparams(0) + Amin - Bmin);    
    Binit(NumB+1) = sqrt(1.0 - bmin);
    PB.Init(&Bgrid, Binit, 0.0, 0.0);

    int NumV = Vparams.rows();
    Vgrid.Init(0.0, CoreRadius, NumV+1);
    Array<double,1> Vinit(NumV+1);
    Vinit(Range(0,NumV-1)) = Vparams;
    Vinit(NumV) = (*FullCoreV)(CoreRadius);
    Vfunc.Init(&Vgrid, Vinit, 5.0e30, FullCoreV->Deriv(CoreRadius));
  }
    
  PH_CubicSpline()
  {
    // Do nothing for now
  }

  void Write (IOSectionClass &outSection)
  {
    outSection.WriteVar ("Type", "PH_CubicSpline");
    // Write common variables
    PseudoHamiltonian::Write(outSection);
    outSection.NewSection ("Agrid");
    Agrid.Write(outSection);
    outSection.CloseSection();
    outSection.NewSection ("Bgrid");
    Bgrid.Write(outSection);
    outSection.CloseSection();
    outSection.NewSection ("Vgrid");
    Vgrid.Write(outSection);
    outSection.CloseSection();
    Array<double,1> Params;
    Params.resize(Agrid.NumPoints);
    for (int i=0; i<Agrid.NumPoints; i++)
      Params(i) = PA.Params(i);
    outSection.WriteVar ("PAparams", Params);

    Params.resize(Bgrid.NumPoints);
    for (int i=0; i<Bgrid.NumPoints; i++)
      Params(i) = PB.Params(i);
    outSection.WriteVar ("PBparams", Params);

    Params.resize(Vgrid.NumPoints);
    for (int i=0; i<Vgrid.NumPoints; i++)
      Params(i) = Vfunc.Params(i);
    outSection.WriteVar ("Vfuncparams", Params);
  }

  bool Read (IOSectionClass &inSection)
  {
    PseudoHamiltonian::Read(inSection);
    assert(inSection.OpenSection("Agrid"));
    Agrid.Read(inSection);
    inSection.CloseSection();
    
    assert(inSection.OpenSection("Bgrid"));
    Bgrid.Read(inSection);
    inSection.CloseSection();

    assert(inSection.OpenSection("Vgrid"));
    Vgrid.Read(inSection);
    inSection.CloseSection();
   
    Array<double,1> Params;
    assert(inSection.ReadVar ("PAparams", Params));
    PA.Init (&Agrid, Params, 0.0, 0.0);
    assert(inSection.ReadVar ("PBparams", Params));
    PB.Init (&Bgrid, Params, 0.0, 0.0);
    assert(inSection.ReadVar ("Vfuncparams", Params));
    Vfunc.Init(&Vgrid, Params, 5.0e30, FullCoreV->Deriv(CoreRadius));
    return (true);
   }
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

  double Amin, Bmin;

  PHType Type()
  {
    return PH_CUBICXC;
  }
  
  int NumParams()
  {
    return (PA.NumParams-1 + PB.NumParams-2 + Vfunc.NumParams-1);
  }
  
  double &Params(int i)
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

  double Params(int i) const
  {
    if (i < (PA.NumParams-1))
      return (PA.Params(i));
    else if (i < (PA.NumParams-1 + PB.NumParams-2))
      return (PB.Params(i-PA.NumParams + 2));
    else 
      return (Vfunc.Params(i-PA.NumParams - PB.NumParams + 3));
  }

	
  void ABV(double r, double &A, double &B, double &V,
		  double &dAdr)
  {
    if (r<CoreRadius)
      {
	if (UpDateB)
	  {
	    PB.Params(0) = sqrt(Amin - Bmin + PA.Params(0)*PA.Params(0));
	    UpDateB = 0;
	  }
	double pa = PA(r);
	double pb = PB(r);
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

  double V(double r)
  {
    double v;
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

  double dVdr(double r)
  {
    double dV;
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

  
   double d2Adr2(double r)
  {
    if (r<CoreRadius)
      {
	double p = PA(r);
	double pp = PA.Deriv(r);
	double ppp = PA.Deriv2(r);
	return (2.0*(p*ppp+pp*pp));
      }
    else
      return (0.0);
  }   

  inline void Init (Array<double,1> Aparams,
		    Array<double,1> Bparams,
		    Array<double,1> Vparams,
		    Grid *agrid, Grid *bgrid,
		    Grid *vgrid, double zion,
		    double coreradius)
  {
    UseVHXC = false;
    NumElecs = 0;
    CoreRadius = coreradius;
    Agrid = agrid;
    Bgrid = bgrid;
    Vgrid = vgrid;
    Zion = zion;
    Amin = 0.0; 
    Bmin = 0.0;

    int NumA = Aparams.rows();
    Array<double,1> Ainit(NumA+1);
    if (Ainit.rows() != Agrid->NumPoints)
      {
	cerr << "Size mismatch in PH_CubicSplineXC.\n";
	exit(1);
      }
    Ainit(Range(0,NumA-1)) = Aparams;
    Ainit(NumA) = sqrt(1.0-Amin);
    PA.Init(Agrid, Ainit, 0.0, 0.0);

    int NumB = Bparams.rows();
    Array<double,1> Binit(NumB+2);
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
    Array<double,1> Vinit(NumV+1);
    if (Vinit.rows() != Vgrid->NumPoints)
      {
	cerr << "Size mismatch in PH_CubicSplineXC.\n";
	cerr << "Vinit.rows() = " << Vinit.rows() << "\n";
	cerr << "Vgrid.NumPoints = " << Vgrid->NumPoints << "\n";
	exit(1);
      }
    Vinit(Range(0,NumV-1)) = Vparams;
    double rmax = Vgrid->End; 

    // Force the last point to agree with asymptotic form
    Vinit(NumV) = -Zion/rmax;
    // Make the slope match assymptotic form as well
    Vfunc.Init(Vgrid, Vinit, 5.0e30, Zion/(rmax*rmax));
  }

  void Write (IOSectionClass &outSection)
  {
    outSection.WriteVar ("Type", "PH_CubicSplineXC");
    // Write common variables
    outSection.WriteVar ("CoreRadius", CoreRadius);
    outSection.WriteVar ("Zion", Zion);
    outSection.NewSection ("Agrid");
    Agrid->Write(outSection);
    outSection.CloseSection();
    outSection.NewSection ("Bgrid");
    Bgrid->Write(outSection);
    outSection.CloseSection();
    outSection.NewSection ("Vgrid");
    Vgrid->Write(outSection);
    outSection.CloseSection();
    Array<double,1> Params;
    Params.resize(Agrid->NumPoints);
    for (int i=0; i<Agrid->NumPoints; i++)
      Params(i) = PA.Params(i);
    outSection.WriteVar ("PAparams", Params);

    Params.resize(Bgrid->NumPoints);
    for (int i=0; i<Bgrid->NumPoints; i++)
      Params(i) = PB.Params(i);
    outSection.WriteVar ("PBparams", Params);

    Params.resize(Vgrid->NumPoints);
    for (int i=0; i<Vgrid->NumPoints; i++)
      Params(i) = Vfunc.Params(i);
    outSection.WriteVar ("Vfuncparams", Params);
  }

  bool Read (IOSectionClass &inSection)
  {
    Amin = 0.0; Bmin = 0.0;
    NumElecs = 0.0;
    assert (inSection.ReadVar("CoreRadius", CoreRadius));
    assert (inSection.ReadVar("Zion", Zion));
    assert(inSection.OpenSection("Agrid"));
    Agrid = ReadGrid(inSection);
    inSection.CloseSection();
    
    assert(inSection.OpenSection("Bgrid"));
    Bgrid = ReadGrid(inSection);
    inSection.CloseSection();

    assert(inSection.OpenSection("Vgrid"));
    Vgrid = ReadGrid(inSection);
    inSection.CloseSection();

    Array<double,1> Params;
    assert(inSection.ReadVar ("PAparams", Params));
    PA.Init (Agrid, Params, 0.0, 0.0);
    assert(inSection.ReadVar ("PBparams", Params));
    PB.Init (Bgrid, Params, 0.0, 0.0);
    assert(inSection.ReadVar ("Vfuncparams", Params));
    Vfunc.Init(Vgrid, Params, 5.0e30, 5.0e30);
    UseVHXC = false;
    return (true);
   }

    
  PH_CubicSplineXC()
  {
    // Do nothing for now
  }
  
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
  double dummy;  // To satisfy compiler warnings about no return value

  PHType Type()
  {
    return PH_NUCLEAR;
  }
  
  int NumParams()
  {
    return (0);
  }
  
  double &Params(int i)
  {
    cerr << "Nuclear PH has no parameters.\n";
    return (dummy);
  }

  double Params(int i) const
  {
    cerr << "Nuclear PH has no parameters.\n";
    return (0.0);
  }

	
  void ABV(double r, double &A, double &B, double &V,
		  double &dAdr)
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

  double V(double r)
  {
    double v = -Zion/r;
    if (UseVHXC)
      if (r < VHXC.grid->End)
	v += VHXC(r);
      else
	v = -(Zion-NumElecs) / r;
    return (v);
  }

  double dVdr(double r)
  {
    double dV;
    dV = Zion / (r*r);
    if (UseVHXC)
      if (r < VHXC.grid->End)
	dV += VHXC.Deriv(r);
      else
	dV = (Zion-NumElecs) / (r*r);

    return (dV);
  }

  double d2Adr2(double r)
  {
    return (0.0);
  }

  bool Read (IOSectionClass &inSection)
  {
    assert(inSection.ReadVar ("Zion", Zion));
    assert(inSection.ReadVar ("UseVHXC", UseVHXC));
    if (UseVHXC)
      {
	assert(inSection.ReadVar ("NumElecs", NumElecs));
	assert(inSection.OpenSection("VHXC"));
	VHXC.Read(inSection);
	inSection.CloseSection();
      }
    return (true);
  }

  void Write (IOSectionClass &outSection)
  {
    outSection.WriteVar ("Type", "Nuclear");
    outSection.WriteVar ("Zion", Zion);
    outSection.WriteVar ("UseVHXC", UseVHXC);
    if (UseVHXC) {
      outSection.WriteVar ("NumElecs", NumElecs);
      outSection.NewSection("VHXC");
      VHXC.Write(outSection);
      outSection.CloseSection();
    } 
  }

  PH_Nuclear(double z)
  {

    Zion = z;
  }

  PH_Nuclear()
  {
    // Do nothing for now
  }
  
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
  double dummy;  // To satisfy compiler warning.


  PHType Type()
  {
    return PH_NUCLEAR;
  }
  
  int NumParams()
  {
    return (0);
  }
  
  double &Params(int i)
  {
    cerr << "Nuclear PH has no parameters.\n";
    return (dummy);
  }

  double Params(int i) const
  {
    cerr << "Nuclear PH has no parameters.\n";
    return (0.0);
  }

	
  void ABV(double r, double &A, double &B, double &V,
		  double &dAdr)
  {
    A = B = 1.0;
    dAdr = 0.0;

    V = (*FullCoreV)(r);
  }

  double V(double r)
  {
    return ((*FullCoreV)(r));
  }

  double dVdr(double r)
  {
    return (FullCoreV->Deriv(r));
  }

  double d2Adr2(double r)
  {
    return (0.0);
  }
  
  
  PH_FullCore(FullCorePotential *FullCorePot)
  {
    Z = FullCorePot->Z;
    FullCoreV = FullCorePot;
  }

  PH_FullCore()
  {

  }

};




///////////////////////////////////////////////////////////////////
// Zero PH:  has no parameters.  Has unit masses and V(r) = 0.   //
// It is used for testing in matrix squaring.                    //
///////////////////////////////////////////////////////////////////
class PH_Zero : public PseudoHamiltonian
{
public:
  double dummy;  // To satisfy compiler warnings about no return value

  PHType Type()
  {
    return PH_NUCLEAR;
  }
  
  int NumParams()
  {
    return (0);
  }
  
  double &Params(int i)
  {
    cerr << "Nuclear PH has no parameters.\n";
    return (dummy);
  }

  double Params(int i) const
  {
    cerr << "Nuclear PH has no parameters.\n";
    return (0.0);
  }

	
  void ABV(double r, double &A, double &B, double &V,
		  double &dAdr)
  {
    A = B = 1.0;
    dAdr = 0.0;

    V = 0.0;
  }

  double V(double r)
  {
    return (0.0);
  }

  double dVdr(double r)
  {
    return (0.0);
  }

  double d2Adr2(double r)
  {
    return (0.0);
  }

  bool Read (IOSectionClass &in) 
  {
    return true;
  }
  void Write (IOSectionClass &outSection)
  {
    outSection.WriteVar ("Type", "Zero");
  }

  PH_Zero()
  {
    // Do nothing for now
  }
};



///////////////////////////////////////////////////////////////////
// Gaussian PH:  has two parameters.  Has unit masses and        //
// V(r) = Amp * exp(-r^2/(2.0*sigma^2)).                         //
// It is used for testing in matrix squaring.                    //
///////////////////////////////////////////////////////////////////
class PH_Gaussian : public PseudoHamiltonian
{
public:
  double Amp, sigma;  // To satisfy compiler warnings about no return value

  PHType Type()
  {
    return PH_NUCLEAR;
  }
  
  int NumParams()
  {
    return (0);
  }
  
  double &Params(int i)
  {
    if (i==0)
      return (Amp);
    else
      return (sigma);
  }

  double Params(int i) const
  {
    if (i==0)
      return (Amp);
    else
      return (sigma);
  }

	
  void ABV(double r, double &A, double &B, double &V,
		  double &dAdr)
  {
    A = B = 1.0;
    dAdr = 0.0;

    V = Amp*exp(-r*r/(2.0*sigma*sigma));
  }

  double V(double r)
  {
    return (Amp*exp(-r*r/(2.0*sigma*sigma)));
  }

  double dVdr(double r)
  {
    return (-r/(sigma*sigma) * Amp*exp(-r*r/(2.0*sigma*sigma)));
  }

  double d2Adr2(double r)
  {
    return (0.0);
  }

  void Write (IOSectionClass &outSection)
  {
    outSection.WriteVar ("Type", "Gaussian");
    outSection.WriteVar ("Amp", Amp);
    outSection.WriteVar ("sigma", sigma);
  }

  bool Read (IOSectionClass &inSection)
  {
    assert (inSection.ReadVar ("Amp", Amp));
    assert (inSection.ReadVar ("sigma", sigma));
    return (true);
  }

  PH_Gaussian(double Amp_, double sigma_)
  {
    Amp = Amp_;
    sigma = sigma_;
    // Do nothing for now
  }

  PH_Gaussian ()
  {
    Amp = 0.0;
    sigma = 1.0;
  }

};


// class Aziz92 : public PseudoHamiltonian
// {
// private:
//   static const double Astar = 1.8443101e5;
//   static const double alphastar = 10.43329537;
//   static const double c6 = 1.36745214;
//   static const double c8 = 0.42123807;
//   static const double c10 = 0.17473318;
// public:


// };


// This potential provides a good primitive action for the coulomb
// potential by making the potential linear at the origin instead of
// divergent. 
class Coulomb3D : public PseudoHamiltonian
{
public:
  double lambda, beta, Z1Z2;
  double U0, dU0;

  // Use Pade form for short-range approximation
  inline double V1(double r)
  {
    double sigma = sqrt(2.0*beta*lambda);
    double alpha = 1.0;
    double a,b,c,d;
    a = (1.0-alpha*Z1Z2*beta/(sigma*U0) + alpha*Z1Z2*beta*dU0/(U0*U0))/
      (alpha*Z1Z2*beta/U0 + sigma*(alpha-1.0));
    b = beta/U0;
    c = a*b - dU0/beta*b*b;
    d = a/Z1Z2;
    
    return (1.0+a*r)/(b+c*r+d*r*r);
  }
  // The above pade has converges slowly for large r, so we use
  // a primitive approximation for large r.
  inline double V2(double r)
  {
    if (r > 1.0e-4)
      return (Z1Z2/r);
    else
      return (1e4*Z1Z2);
  }
  
  // This is a mixing function to transition from V1 to V2;
  inline double f(double r)
  {
    double sigma = sqrt(2.0*beta*lambda);
    if (r < sigma)
      return (1.0);
    else
      return (exp(-3.0*(r/sigma - 1.0)*(r/sigma-1.0)));
  }

  double V(double r)
  {
    double F = f(r);
    return (F*V1(r) + (1.0-F)*V2(r));
  }

  PHType Type()
  {
    return PH_COULOMB3D;
  }
  
  int NumParams()
  {
    return (0);
  }
  
  double &Params(int i)
  {
    return Z1Z2;
  }

  double Params(int i) const
  {
    return (0.0);
  }
	
  void ABV(double r, double &A, double &B, double &Vr,
	   double &dAdr)
  {
    A = B = 1.0;
    dAdr = 0.0;

    Vr=V(r);
  }
  // We don't need this
  double dVdr(double r)
  {    return (0.0); }

  double d2Adr2(double r)
  {
    return (0.0);
  }

  void Write (IOSectionClass &outSection)
  {
    outSection.WriteVar ("Type", "Coulomb3D");
    outSection.WriteVar ("Z1Z2", Z1Z2);
    outSection.WriteVar ("lambda", lambda);
    outSection.WriteVar ("beta", beta);
  }

  bool Read (IOSectionClass &inSection)
  {
    assert (inSection.ReadVar ("lambda", lambda));
    assert (inSection.ReadVar ("Z1Z2", Z1Z2));
    assert (inSection.ReadVar ("beta", beta));
    Init(lambda, beta, Z1Z2);
    return (true);
  }

  void Init (double lambda_, double beta_, double Z1Z2_)
  {
    lambda = lambda_;
    beta = beta_;
    Z1Z2 = Z1Z2_;
    
    // Calculate potential action at r=rp=0 -- U0;
    TinyVector<double,7> Pj;
    Pj = 1.772453851, -0.074137740, 0.005834805, -0.000382686,
      0.000008738, 0.000002138, 0.000000356;
    double Z = Z1Z2/lambda;
    U0 = 0.0;
    double gamma = lambda*beta*Z*Z;
    double sqrtgamma = sqrt(gamma);
    double gamma2jby2 = sqrtgamma;
    double s = Z / fabs(Z);
    double sign = s;
    for (int j=0; j<7; j++) {
      U0 += sign * Pj(j) * gamma2jby2;
      gamma2jby2 *= sqrtgamma;
      sign*= s;
    }
    // Calculate slope
    dU0 = -Z;
  }
  Coulomb3D(double lambda_, double beta_, double Z1Z2_)
  {
    Init (lambda_, beta_, Z1Z2_);
  }
  Coulomb3D() { }
};





class CoulombPot : public PseudoHamiltonian
{
public:
  double Z1Z2;

  double V(double r)
  {
    // HACK 
    if (r > 0.001)
      return (Z1Z2/r);
    else
      return (1e3*Z1Z2);
  }

  PHType Type()
  {
    return PH_COULOMB;
  }
  
  int NumParams()
  {
    return (0);
  }
  
  double &Params(int i)
  {
    return Z1Z2;
  }

  double Params(int i) const
  {
    return (0.0);
  }
	
  void ABV(double r, double &A, double &B, double &Vr,
	   double &dAdr)
  {
    A = B = 1.0;
    dAdr = 0.0;

    Vr=V(r);
  }
  // We don't need this
  double dVdr(double r)
  {    return (0.0); }

  double d2Adr2(double r)
  {
    return (0.0);
  }

  void Write (IOSectionClass &outSection)
  {
    outSection.WriteVar ("Type", "Coulomb");
    outSection.WriteVar ("Z1Z2", Z1Z2);
  }

  bool Read (IOSectionClass &inSection)
  {
    assert (inSection.ReadVar ("Z1Z2", Z1Z2));
    return (true);
  }

  CoulombPot() { }
};





PseudoHamiltonian *ReadPH (IOSectionClass &inSection);

#endif
