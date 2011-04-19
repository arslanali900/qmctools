#ifndef CBN_MODEL_H
#define CBN_MODEL_H

#include "VinetEOS2.h"
#include "NonlinearFit.h"
#include <vector>
#include "PhononFreeEnergy.h"
#include "DebyeModel.h"
#include "RamanFreq.h"
#include "Raman.h"



class cBNModel
{
public:
  VinetEOSClass Static;
  SplineFreeEnergy Phonon;
  DebyeModel Debye;
  
  RamanModel Raman;
  inline double P(double V, double T)
  { return Static.Pressure(V) + Phonon.P(V,T); }
  double FindV(double P, double T);

  void SetStatic (string fname);
  void SetPhonon (string fname);
  void SetRaman  (string fname);
  void SetAll (string staticName, 
	       string phononName, 
	       string ramanName);

  double MeanFrequency_PT (double P, double T);
  void Write_nuPT_Table();
};

#endif
