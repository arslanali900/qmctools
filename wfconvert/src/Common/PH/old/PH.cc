#include "PH.h"
#include "../IO/InputFile.h"


PseudoHamiltonian* ReadPH (IOSectionClass &inSection)
{
  string Type;
  PseudoHamiltonian *PH;
  assert (inSection.ReadVar ("Type", Type));
  if (Type == "PH_CubicSpline")
    PH = new PH_CubicSpline;
  else if (Type == "Nuclear")
    PH = new PH_Nuclear;
  else if (Type == "Zero")
    PH = new PH_Zero;
  else if (Type == "Gaussian")
    PH = new PH_Gaussian;
  else if (Type == "FullCore")
    PH = new PH_FullCore;
  else if (Type == "PH_CubicSplineXC")
    PH = new PH_CubicSplineXC;
  else if (Type == "Coulomb3D")
    PH = new Coulomb3D;
  else if (Type == "Coulomb")
    PH = new CoulombPot;
  else {
    cerr << "Unrecognized PH type \"" << Type << "\".  Aborting.\n";
    abort();
  }

  PH->Read(inSection);
  return (PH);
}







  
