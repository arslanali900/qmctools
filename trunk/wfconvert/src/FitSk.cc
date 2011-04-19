#include "FitSk.h"
#include "ParserClass.h"
#include "ParseCommand.h"
#include <Common/Fitting/Fitting.h>
#include <Common/Fitting/NonlinearFit.h>


template<typename T> void
RecipSet<T>::AddVal (double Gmag, T sk)
{
  typename std::vector<RecipVal<T> >::iterator iter = Vals.begin();
  while ((iter!=Vals.end()) && (Gmag > (iter->Gmag + 1.0e-10)))
    iter++;
  if (iter == Vals.end()) {
    RecipVal<T> newVal;
    newVal.Gmag = Gmag;
    newVal.Degeneracy = 1.0;
    newVal.ValSum = sk;
    Vals.push_back (newVal);
  }
  else if (fabs (iter->Gmag - Gmag) < 1.0e-9) {
    iter->Degeneracy += 1.0;
    iter->ValSum += sk;
  }
  else {
    RecipVal<T> newVal;
    newVal.Gmag = Gmag;
    newVal.Degeneracy = 1.0;
    newVal.ValSum = sk;
    Vals.insert (iter, newVal);
  }
}


bool
SkClass::Read (string fname)
{
  MemParserClass memParser;
  FileParserClass2 fileParser;
  streamsize fsize = memParser.FileSize(fname);
  if (fsize == -1)
    return false;
  
  ParserClass *parserPtr;
  if (fsize < (streamsize)(1<<28))
    parserPtr = &memParser;
  else
    parserPtr = &fileParser;
  
  ParserClass &parser = *parserPtr;
  if (!parser.OpenFile (fname))                                  return false;
  if (!parser.FindToken ("File version"))                        return false;
  if (!parser.FindToken ("\n"))                                  return false;
  int version;
  if (!parser.ReadInt (version))                                 return false;
  if (version != 1)                                              return false;
  int numTypes, numPtcls;
  if (!parser.FindToken ("Number of particle types"))            return false;
  if (!parser.FindToken ("\n"))                                  return false;
  if (!parser.ReadInt (numTypes))                                return false;
  if (!parser.FindToken ("Number of each type of particle"))     return false;
  NumParticles = 0;
  for (int i=0; i<numTypes; i++) {
    if (!parser.ReadInt (numPtcls))                              return false;
    NumParticles += numPtcls;
  }
  


  if (!parser.FindToken ("Primitive translation vectors (au)"))  return false;
  for (int i=0; i<3; i++)
    for (int j=0; j<3; j++)
      if (!parser.ReadDouble (PrimLattice(i,j)))                 return false;
  if (!parser.FindToken ("Multiples of primitive translation vectors"))
    return false;
  for (int i=0; i<3; i++)
    if (!parser.ReadInt (Supercell[i]))                          return false;
  if (!parser.FindToken ("Volume of simulation cell"))           return false;
  if (!parser.ReadDouble (BoxVol))                               return false;
  
  int numG;
  if (!parser.FindToken ("START GVECTOR SET 1"))                 return false;
  if (!parser.FindToken ("Number of G-vectors in set"))          return false;
  if (!parser.ReadInt (numG))                                    return false;
  if (!parser.FindToken ("G-vector components Gx, Gy, Gz (au)"))
    return false;
  Array<Vec3,1> gvecs(numG);
  for (int i=0; i<numG; i++)
    for (int j=0; j<3; j++)
      if (!parser.ReadDouble (gvecs(i)[j]))                      return false;

  if (!parser.FindToken ("START STRUCTURE FACTOR"))              return false;
  if (!parser.FindToken ("Number of sets"))                      return false;
  if (!parser.FindToken ("\n"))                                  return false;
  int numSets;
  if (!parser.ReadInt (numSets))                                 return false;
  Sk.resize (numSets, numG);
  SkSets.resize(numSets);
  for (int set=0; set<numSets; set++) {
    if (!parser.FindToken ("START SET"))                         return false;
    if (!parser.FindToken ("Types of particles in pair"))        return false;
    if (!parser.ReadInt (SkSets[set].ParticleType1))             return false;
    if (!parser.ReadInt (SkSets[set].ParticleType2))             return false;
    SkSets[set].ParticleType1--;
    SkSets[set].ParticleType2--;
    if (!parser.FindToken ("rho_a(G)*rho_b(-G)"))                return false;
    for (int i=0; i<numG; i++) {
      if (!parser.ReadDouble (Sk(set,i)))                        return false;
      double Gmag = sqrt(dot (gvecs(i), gvecs(i)));
    }
  }

  if (!parser.FindToken ("Number of sets for spin density part")) return false;
  int numRho;
  if (!parser.ReadInt (numRho))                                   return false;
  Rhok.resize (numRho, numG);
  RhoSets.resize(numRho);
  for (int set=0; set<numRho; set++) {
    if (!parser.FindToken ("START SET"))                          return false;
    if (!parser.FindToken("Complex charge-density coefficients"))return false;
    if (!parser.FindToken ("\n"))                                 return false;
    for (int i=0; i<numG; i++) {
      double rpart, ipart;
      if (!parser.ReadDouble (rpart))                             return false;
      if (!parser.ReadDouble (ipart))                             return false;
      Rhok (set, i) = complex<double>(rpart, ipart);
    }
  }
  
  for (int set=0; set<numSets; set++) {
    int ptype1 = SkSets[set].ParticleType1;
    int ptype2 = SkSets[set].ParticleType2;
    double prefactor = (ptype1 == ptype2) ? 1.0 : 2.0;
    for (int iG=0; iG<numG; iG++) {
      double gmag = sqrt (dot (gvecs(iG), gvecs(iG)));
      complex<double> val = prefactor * ((2.0*Sk(set,iG)) - 
	2.0*Rhok(ptype1,iG)*conj(Rhok(ptype2,iG)));
      int numCells = Supercell[0]*Supercell[1]*Supercell[2];
      val *= (double) numCells*numCells;
      val /= (double) NumParticles;
      
      SkSets[set].AddVal (gmag, 0.5*val.real());
    }
  }
        
  FILE *fout = fopen ("Sk.dat", "w");
  for (int i=0; i<numG; i++) {
    double Gmag = sqrt(dot (gvecs(i), gvecs(i)));
    fprintf (fout, "%20.14f ", Gmag);
    for (int set=0; set<numSets; set++)
      fprintf (fout, "%20.14f ", Sk(set, i));
    fprintf (fout, "\n");
  }
  fclose (fout);
  
  fout = fopen ("SkAvg.dat", "w");
  for (int i=0; i<SkSets[0].Vals.size(); i++) {
    fprintf (fout, "%21.14e ", SkSets[0].Vals[i].Gmag);
    for (int set=0; set<SkSets.size(); set++)
      fprintf (fout, "%21.14e ", 
	       SkSets[set].Vals[i].ValSum / 
	       SkSets[set].Vals[i].Degeneracy);
    fprintf (fout, "\n");
  }
  fclose (fout);    
    

  return true;
}

double
SkClass::BlendFunc (double k)
{
  //return (1.0/(exp((k-kcut)/(0.5*width))+1.0));
  return exp (-k*k/(0.5*width*width));
}

void
SkClass::Fit (int numKnots, double k_c)
{
  Basis.SetNumKnots (numKnots);
  Basis.SetCutoff (k_c);
  FitCoefs.resize(SkSets.size(), Basis.NumElements());
  for (int set=0; set < SkSets.size(); set++) {
    int N = SkSets[set].Vals.size();
    Array<double,1> y(N), sigma(N), error(N);
    Array<double,2> F(N, Basis.NumElements());
    Array<double,1> coefs (Basis.NumElements());
    Array<bool,1> adjust (Basis.NumElements());
    adjust = true;
    adjust(0) = false;
    adjust(1) = false;
    coefs = 0.0;
    y.resize (N);
    sigma.resize(N);
    error.resize(N);
    for (int i=0; i<N; i++) {
      double val =  SkSets[set].Vals[i].ValSum/SkSets[set].Vals[i].Degeneracy;
      y(i) = val;
      sigma(i) = 0.01;
      double g = SkSets[set].Vals[i].Gmag;
      for (int j=0; j<Basis.NumElements(); j++) 
	F(i,j) = Basis (j, g);
    }
    LinFitSVD (y, sigma, F, adjust, coefs, error, 1.0e-14);
    FitCoefs(set, Range::all()) = coefs;
  }

  FILE *fout = fopen ("SkFit.dat", "w");
  for (double k=0.0; k<5.0; k+=0.001) {
    fprintf (fout, "%10.6f ", k);
    for (int set=0; set < SkSets.size(); set++) {
      double val = 0.0;
      for (int i=0; i<Basis.NumElements(); i++) 
	val += FitCoefs(set, i) * Basis(i, k);
      fprintf (fout, "%24.16e ", val);
    }
    fprintf (fout, "\n");
  }
  fclose (fout);

  double upper_bound = k_c;

  kcut  = 300.0;
  width = 1.0e20;;
  for (int set=0; set < SkSets.size(); set++) {
    double integ = 0.0;
    for (int i=0; i<Basis.NumElements(); i++)
      integ += FitCoefs(set,i) * Basis.Integral (i, upper_bound);
    integ -= upper_bound;
    integ *= 1.0/M_PI;

    fprintf (stderr, "Integral of fit    for set %d = %1.10f\n",
	     set, integ);
    int N = SkSets[set].Vals.size();
    double sum = Sum (set, 0, N-1);
    fprintf (stderr, "Sum of data points for set %d = %1.10f\n",
	     set, sum);
    fprintf (stderr, "Correction for set         %d = %1.10f\n",
	     set, integ - sum);
  }
  kcut  = 2.0;
  width = 0.22;
  double total_corr = 0.0;
  for (int set=0; set < SkSets.size(); set++) {
    double integ = 0.0;
    for (double k=0.0; k<=upper_bound; k+=0.0001) {
      double sk = 0.0;
      for (int i=0; i<Basis.NumElements(); i++) 
	sk += (FitCoefs(set,i) * Basis(i, k));
      integ += 0.0001 * BlendFunc(k)*(sk);
    }
    integ *= 1.0/M_PI;
    int N = SkSets[set].Vals.size();
    double sum = Sum (set, 0, N-1);    
    fprintf (stderr, "Num. Integral of fit for set %d = %1.10f\n",
	     set, integ);
    fprintf (stderr, "Sum of data points   for set %d = %1.10f\n",
	     set, sum);
    fprintf (stderr, "Correction           for set %d = %1.10f\n",
	     set, integ - sum);
    total_corr += (integ-sum);
  }

  // Total fit
  kcut  = 300.0;
  width = 0.5;
  int N = SkSets[0].Vals.size();
  Array<double,1> y(N), sigma(N), error(N);
  Array<double,2> F(N, Basis.NumElements());
  Array<double,1> coefs (Basis.NumElements());
  Array<bool,1> adjust (Basis.NumElements());
  adjust = true;
  adjust(0) = false;
  adjust(1) = false;
  adjust(Basis.NumElements()-2) = false;
  coefs = 0.0;
  y.resize (N);
  sigma.resize(N);
  error.resize(N);
  for (int i=0; i<N; i++) {
    y(i) =  0.0;
    y(i) += SkSets[0].Vals[i].ValSum/SkSets[0].Vals[i].Degeneracy;
    y(i) += SkSets[1].Vals[i].ValSum/SkSets[1].Vals[i].Degeneracy;
    y(i) += SkSets[2].Vals[i].ValSum/SkSets[2].Vals[i].Degeneracy;
    sigma(i) = 0.01;
    double g = SkSets[0].Vals[i].Gmag;
    for (int j=0; j<Basis.NumElements(); j++) 
      F(i,j) = Basis (j, g);
  }
  LinFitSVD (y, sigma, F, adjust, coefs, error, 1.0e-14);
  coefs;
  double integ = 0.0;
  for (int i=0; i<Basis.NumElements(); i++)
    integ += coefs(i) * Basis.Integral (i, upper_bound);
  //integ -= upper_bound;
  integ *= 1.0/M_PI;
  
  fprintf (stderr, "Integral of fit  for total = %1.10f\n",
	   integ);
  double sum = TotalSum (0, N-1);
  fprintf (stderr, "Sum of data points for set = %1.10f\n",
	   sum);
  fprintf (stderr, "Correction for total       = %1.10f\n",
	   integ - sum);
  FILE *tout = fopen ("Sktotal.dat", "w");
  for (double k=0.0; k<=k_c; k+=0.001) {
    double fit = 0.0;
    for (int i=0; i<Basis.NumElements(); i++)
      fit += coefs(i) * Basis (i, k);
    fprintf (tout, "%24.16e %24.16e\n", k, fit);
  }
  fclose (tout);

  double k = 1.0e-5;
  double val = 0.0;
  for (int i=0; i<Basis.NumElements(); i++)
    val += coefs(i)*Basis (i,k);
  double correction = 2.0*M_PI/ BoxVol * val/(k*k);
  fprintf (stderr, "Analytic correction  = %1.8f\n", correction);
  fprintf (stderr, "Numerical correction = %1.8f\n", total_corr);

}


double
SkClass::Sum (int set, int firstk, int lastk)
{
  double sum = 0.0;
  double prefactor = 2.0*M_PI/BoxVol;
  for (int ki= firstk; ki<=lastk; ki++) {
    double gmag = SkSets[set].Vals[ki].Gmag;
    double degen = SkSets[set].Vals[ki].Degeneracy;
    double sk   = SkSets[set].Vals[ki].ValSum;
    if (fabs(gmag) > 1.0e-14)
      sum += prefactor * BlendFunc(gmag)*(sk)/(gmag*gmag);
  }
  return sum;
}


double
SkClass::TotalSum (int firstk, int lastk)
{
  double sum = 0.0;
  double prefactor = 2.0*M_PI/BoxVol;
  for (int ki= firstk; ki<=lastk; ki++) {
    double degen = SkSets[0].Vals[ki].Degeneracy;
    double gmag  = SkSets[0].Vals[ki].Gmag;
    for (int set=0; set<SkSets.size(); set++) {
      double sk   = SkSets[set].Vals[ki].ValSum;
      if (fabs(gmag) > 1.0e-10)
	sum += prefactor * (sk)/(gmag*gmag);
    }
//     if (fabs(gmag) > 1.0e-10)
//       sum -= prefactor * degen/(gmag*gmag);
  }
  return sum;
}

template<int M>
class ExpModel
{
private:
  TinyVector<double,M> Coefs;
public:
  inline void SetParams (TinyVector<double,M> params)
  { Coefs = params; }
  inline double GetParm (int i) { return Coefs[i]; }
  inline double operator()(double k);
  inline TinyVector<double,M> Grad (double k);
};

template<int M>
inline double ExpModel<M>::operator()(double k)
{
  double sum = 0.0;
  double prod = k*k;
  for (int i=0; i<M; i++) {
    sum += Coefs[i] * prod;
    prod *= k;
  }
  return 1.0-exp(-sum);
}

// template<int M>
// inline TinyVector<double,M> ExpModel<M>::Grad (double k)
// {
//   const double eps = 1.0e-7;
//   TinyVector<double,M> grad;
//   for (int i=0; i<M; i++) {
//     double save = Coefs[i];
//     Coefs[i] = save + eps;
//     double plus  = (*this)(k);
//     Coefs[i] = save - eps;
//     double minus = (*this)(k);
//     Coefs[i] = save;
//     grad[i] = (plus-minus)/(2.0*eps);
//   }
//   return grad;
// }

template<int M>
inline TinyVector<double,M> ExpModel<M>::Grad (double k)
{
  TinyVector<double,M> grad;
  for (int i=0; i<M; i++)
    grad[i] = 0.0;
  double sum = 0.0;
  double prod = k*k;
  for (int i=0; i<M; i++) {
    sum += Coefs[i] * prod;
    prod *= k;
  }
  prod = k*k;
  for (int i=0; i<M; i++) {
    grad[i] = -prod;
    prod *= k;
  }
  grad = -exp(-sum) * grad;
  return grad;
}





// Fit to form 1-exp(-\sum_{n=2}^{2+numTerms-1}a_n k^n)
template<int M> void
SkClass::FitNonlinear ()
{
  typedef ExpModel<M> ModelType;
  ModelType model;
  NonlinearFitClass<M,ModelType> NLfit(model);
  
  int N = SkSets[0].Vals.size();
  Array<double,1> x(N-1), y(N-1), sigma(N-1), error(N-1);
  double alpha = 1.0;
  double salpha = 1.0;
  for (int i=1; i<N; i++) {
    x(i-1) = SkSets[0].Vals[i].Gmag;
    y(i-1) =  0.0;
    y(i-1) += SkSets[0].Vals[i].ValSum/SkSets[0].Vals[i].Degeneracy;
    y(i-1) += SkSets[1].Vals[i].ValSum/SkSets[1].Vals[i].Degeneracy;
    y(i-1) += SkSets[2].Vals[i].ValSum/SkSets[2].Vals[i].Degeneracy;
    //    sigma(i) = SkSets[0].Vals[i].Gmag * SkSets[0].Vals[i].Gmag +
    //    0.001;
    //sigma(i) = exp (x(i)*x(i)/(4.0*salpha*salpha));
    sigma(i-1) = 0.001+ x(i-1)*x(i-1);
    //    sigma(i) = 0.1;
  }
  TinyVector<double,M> params;
  params = 0.001;
  params(0) = 0.1;
  NLfit.Fit (x, y, sigma, params);
  model.SetParams (params);
  char name[100];
  snprintf (name, 100, "nonlin%d.dat", M);
  FILE *fout = fopen (name, "w");
  for (double k=0.0; k<SkSets[0].Vals[N-1].Gmag; k+=0.001) {
    double val = model (k);
    fprintf (fout, "%24.16e %24.16e\n", k, val);
  }
  fclose (fout);
  double corr = params[0] * 2.0*M_PI/BoxVol;
  fprintf (stderr, 
	   "Nonlinear fit with %2d terms, leading-order correction = %1.8f\n",
	   M, corr);
  fprintf (stderr, "  Correction per prim. cell = %1.8f\n",
	   corr*NumParticles/(double)(Supercell[0]*Supercell[1]*Supercell[2]));

  // Now try a numerical integral - sum difference
  double sum = 0.0, integral = 0.0;
  for (int i=1; i<N; i++) {
    double val;
    val  = 1.0 * SkSets[0].Vals[i].ValSum;
    val += 1.0 * SkSets[1].Vals[i].ValSum;
    val += 1.0 * SkSets[2].Vals[i].ValSum;
    double g = SkSets[0].Vals[i].Gmag;
    val *= exp (-g*g/(2.0* alpha*alpha));
    val /= g*g;
    sum += val;
  }
  sum *= 2.0*M_PI/BoxVol;
  double gmax = SkSets[0].Vals[N-1].Gmag;
  for (double k=0.0; k<=gmax; k+=0.00001) {
    double val = model(k)*exp(-k*k/(2.0*alpha*alpha));
    integral += 0.00001*val;
  }
  integral *= 1.0/M_PI;
//   fprintf (stderr, "integral = %1.8f  sum = %1.8f  diff = %1.8f\n",
// 	   integral, sum, integral-sum);
}



void
SkClass::FitRon (int numStars)
{
  for (int set=0; set < SkSets.size(); set++) {
    int N = numStars;
    Array<double,1> y(N), sigma(N), error(N);
    Array<double,2> F(N, 3);
    Array<double,1> coefs (3);
    Array<bool,1> adjust (2);
    adjust = true;
    coefs = 0.0;
    y.resize (N);
    sigma.resize(N);
    error.resize(N);
    for (int i=0; i<N; i++) {
      double val =  SkSets[set].Vals[i].ValSum/SkSets[set].Vals[i].Degeneracy;
      y(i) = val;
      sigma(i) = 0.01;
      double g = SkSets[set].Vals[i].Gmag;
      F(i,0) = g*g;
      F(i,1) = g*g*g;
      F(i,2) = g*g*g*g;
    }
    LinFitSVD (y, sigma, F, adjust, coefs, error, 1.0e-14);
    FitCoefs(set, Range::all()) = coefs;
    double gmax = 0.5*(SkSets[set].Vals[N-1].Gmag +
		       SkSets[set].Vals[N].Gmag);
    double integral = 1.0/(M_PI) * 
      (1.0/3.0*gmax*gmax*gmax      * FitCoefs(set, 0) +
       1.0/4.0*gmax*gmax*gmax*gmax * FitCoefs(set, 1));
    double sum = Sum (set, 0, numStars);
    // fprintf (stderr, "Set = %d, num G points = %d:\n", set, numStars);
    // fprintf (stderr, "  Ron's integral part = %12.8f\n", integral);
    // fprintf (stderr, "  Ron's sum part      = %12.8f\n", sum);
    // fprintf (stderr, "  Ron's correction          = %12.8f\n", integral - sum);
  }
  fprintf (stderr, "  Ron's analytic correction for %2d stars = %12.8f\n", 
	   numStars,
	   2.0*M_PI/ BoxVol * (FitCoefs (0,0)+FitCoefs(1,0)+FitCoefs(2,0))); 

}


void
SkClass::FitRon2 (double kcut)
{
  int numStars = 0;
  for (int i=0; i<SkSets[0].Vals.size(); i++)
    if (SkSets[0].Vals[i].Gmag <= kcut)
      numStars++;
  
  for (int set=0; set < SkSets.size(); set++) {
    int N = numStars;
    Array<double,1> y(N), sigma(N), error(N);
    Array<double,2> F(N, 3);
    Array<double,1> coefs (3);
    Array<bool,1> adjust (2);
    adjust = true;
    coefs = 0.0;
    y.resize (N);
    sigma.resize(N);
    error.resize(N);
    for (int i=0; i<N; i++) {
      double val =  SkSets[set].Vals[i].ValSum/SkSets[set].Vals[i].Degeneracy;
      y(i) = val;
      sigma(i) = 0.01;
      double g = SkSets[set].Vals[i].Gmag;
      F(i,0) = g*g;
      F(i,1) = g*g*g;
      F(i,2) = g*g*g*g;
    }
    LinFitSVD (y, sigma, F, adjust, coefs, error, 1.0e-14);
    FitCoefs(set, Range::all()) = coefs;
    //    double gmax = 0.5*(SkSets[set].Vals[N-1].Gmag +
    //    SkSets[set].Vals[N].Gmag);
    double gmax = kcut;
    double integral = 1.0/(M_PI) * 
      (1.0/3.0*gmax*gmax*gmax      * FitCoefs(set, 0) +
       1.0/4.0*gmax*gmax*gmax*gmax * FitCoefs(set, 1));
    double sum = Sum (set, 0, numStars);
    // fprintf (stderr, "Set = %d, num G points = %d:\n", set, numStars);
    // fprintf (stderr, "  Ron's integral part = %12.8f\n", integral);
    // fprintf (stderr, "  Ron's sum part      = %12.8f\n", sum);
    // fprintf (stderr, "  Ron's correction          = %12.8f\n", integral - sum);
  }
  fprintf (stderr, "  Ron's analytic correction for %2d stars = %12.8f\n", 
	   numStars,
	   2.0*M_PI/ BoxVol * (FitCoefs (0,0)+FitCoefs(1,0)+FitCoefs(2,0))); 

}



main(int argc, char **argv) 
{
  list<ParamClass> argList;
  argList.push_back(ParamClass("order", false));

  CommandLineParserClass parser (argList);
  bool success = parser.Parse (argc, argv);

  if (!success || parser.NumFiles() != 1) {
    cerr << "Usage:  fitsk [options...] infile\n";
    exit (1);
  }

  string inName = parser.GetFile(0);

  SkClass Sk;
  bool goodRead = Sk.Read (inName);
  if (goodRead) {
//     cerr << "Read succeeded.\n";
//     for (int set=0; set<Sk.SkSets.size(); set++) {
//       double sum = Sk.Sum (set, 0, Sk.SkSets[set].Vals.size()-1);
//       fprintf (stderr, "Set = %d, sum = %1.8f\n", set, sum);
//     }
    int N = Sk.SkSets[0].Vals.size();
    double gmax = Sk.SkSets[0].Vals[N-1].Gmag;
    Sk.Fit (2, gmax);
    Sk.FitNonlinear<1>();
    Sk.FitNonlinear<2>();
    Sk.FitNonlinear<3>();
    Sk.FitNonlinear<4>();
    Sk.FitNonlinear<5>();
    Sk.FitNonlinear<6>();
    Sk.FitNonlinear<7>();
    Sk.FitNonlinear<8>();
    Sk.FitNonlinear<9>();
    Sk.FitNonlinear<10>();
//     Sk.FitRon (3);
//     Sk.FitRon (4);
//     for (int i=2; i<100; i++)
//       Sk.FitRon (i);
//     Sk.FitRon (7);
//     Sk.FitRon (8);
//     Sk.FitRon (9);
//     Sk.FitRon (10);
  }
  else
    cerr << "Read failed.\n";
  
}
