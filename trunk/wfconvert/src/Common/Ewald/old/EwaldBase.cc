#include "EwaldBase.h"

void EwaldClass::SetPot (Potential &pot)
{
  Pot = &pot;
}

void EwaldClass::SetBox(TinyMatrix<double,3,3> latvecs)
{
  LatVecs = latvecs;
}

void EwaldClass::SetBox(TinyVector<double,3> box)
{
  LatVecs = 0.0;
  LatVecs(0,0) = box(0);
  LatVecs(1,1) = box(1);
  LatVecs(2,2) = box(2);
}

void EwaldClass::SetZ(double Z)
{
  Z = 0.0;
}

EwaldClass::EwaldClass()
{
  LatVecs = 0.0;
  CutOff = 0.0;
  Z = 0.0;
  Pot = NULL;
}
