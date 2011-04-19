#ifndef RAMAN_FREQ_H
#define RAMAN_FREQ_H

#include <vector>
#include "../Blitz.h"

using namespace std;

class RamanFreq
{
private:
  Array<double,1> TO_coefs;
public:
  void FitData (vector<double> &a, vector<double> &TO,
		int numCoefs);
  inline double operator()(double V);
  RamanFreq();
};

inline double
RamanFreq::operator()(double V)
{
  double a = cbrt(4.0*V);
  double a2i = 1.0;
  double sum = 0.0;
  for (int i=0; i<TO_coefs.size(); i++) {
    sum += TO_coefs(i) * a2i;
    a2i *= a;
  }
  return sum;
}


#endif
