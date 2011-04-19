#ifndef COST_H
#define COST_H

#include "Atom.h"
#include "../Optimize/Minimize.h"

class PHCost : public MinimizeFunction
{
public:
  Atom *FCatom;
  Atom *Patom;
  scalar NegativityCoef, PartialNormCoef, NodeCoef, LogDerivCoef,
    OverMaxCoef, CurveCoef, MatchCoef; //, MaxCurve;
  bool UseLog;
  char TempPHName[300];

  scalar WFMatchCost();

  int NumParams()
  {
    return (Patom->PH->NumParams());
  }
  scalar &Params(int i)
  {
    return (Patom->PH->Params(i));
  }
  scalar Params(int i) const
  {
    return (Patom->PH->Params(i));
  }

  scalar Cost();
};

#endif
