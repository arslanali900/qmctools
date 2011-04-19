#include "Ewald.h"

double NaCl_Madelung()
{
  CoulombPot NaCl, NaNa, ClCl;
  NaCl.Z1Z2 = -1.0;
  ClCl.Z1Z2 = 1.0;
  NaNa.Z1Z2 = 1.0;

  TinyVector<double,3> box;
  box[0] = 1.0; box[1] = 1.0; box[2] = 1.0;
  SimpleEwald ewald_NaCl, ewald_NaNa, ewald_ClCl;
  ewald_NaCl.SetPot (NaCl);
  ewald_NaNa.SetPot (NaNa);
  ewald_ClCl.SetPot (ClCl);

  ewald_NaCl.SetBox(box);
  ewald_NaNa.SetBox(box);
  ewald_ClCl.SetBox(box);

  ewald_NaCl.SetZ(-1.0);
  ewald_NaCl.SetZ(-1.0);
  

}


main()
{
  double Madelung = NaCl_Madelung();
  cerr << "NaCl Madelung constant = " << Madelung << endl;
}
