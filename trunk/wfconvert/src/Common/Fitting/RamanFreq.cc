#include "RamanFreq.h"
#include "Fitting.h"

void
RamanFreq::FitData (vector<double> &a, vector<double> &TO,
		    int numCoefs)
{
  TO_coefs.resize(numCoefs);
    int N = a.size();
  Array<double,2> basis(N, numCoefs);
  for (int i=0; i<N; i++) {
    basis(i,0) = 1.0;
    for (int j=1; j<numCoefs; j++)
      basis(i,j) = a[i]*basis(i,j-1);
  }
  TO_coefs.resize(numCoefs);
  Array<double,1> sigma(N), freq(N), errors(numCoefs);
  for (int i=0; i<N; i++) {
    sigma(i) = 1.0e-7;
    freq(i) = TO[i];
  }
  LinFitLU(freq, sigma, basis, TO_coefs, errors);
 
}


RamanFreq::RamanFreq()
{
  vector<double> a, TO, LO;
  a.push_back(7.1138   );  TO.push_back( 885.8568);  LO.push_back(1128.7219);
  a.push_back(6.8326828);  TO.push_back(1032.7117);  LO.push_back(1266.0320);
  a.push_back(6.5989   );  TO.push_back(1171.7002);  LO.push_back(1395.6597);
  a.push_back(6.3684   );  TO.push_back(1326.3245);  LO.push_back(1539.8697);
  a.push_back(6.1459   );  TO.push_back(1493.4728);  LO.push_back(1696.1205);
  a.push_back(6.0000   );  TO.push_back(1613.6059);  LO.push_back(1808.7275);
  a.push_back(5.4288835);  TO.push_back(2185.4867);  LO.push_back(2348.5609);
  a.push_back(6.69433  );  TO.push_back(1112.8076);  LO.push_back(1340.7611);
  a.push_back(6.97953  );  TO.push_back( 953.4158);  LO.push_back(1191.9357);
  a.push_back(7.243156 );  TO.push_back( 825.3917);  LO.push_back(1071.7322);
  
  FitData(a, TO, 3);
}
