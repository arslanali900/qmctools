#ifndef EWALD_BASE_H
#define EWALD_BASE_H

#include "../PH/Potential.h"

class EwaldClass
{
protected:
  Potential *Pot;
  TinyMatrix<double,3,3> LatVecs;
  double CutOff;
  double Z;
public:
  Array<double,4> kVecs;
  void SetPot(Potential &pot);
  void SetBox(TinyVector<double, 3> box);
  void SetBox(TinyMatrix<double,3,3> latVecs);
  void SetZ (double z);
  void MakekVecs (double kcutoff);

  virtual double Vshort(double r) = 0;
  virtual void BreakUp (double rcutoff) = 0;
  virtual double Vlong_k (double k) = 0;
  EwaldClass();
};


#endif
