#include "NonlinearFit.h"
#include "../MatrixOps/MatrixOps.h"

template<int M, typename ModelType> double
NonlinearFitClass<M,ModelType>::Chi2(const Array<ArgType,1> &x, 
				     const Array<double,1> &y,
				     const Array<double,1> &sigma,
				     TinyVector<double,M> params)
{
  int N = x.size();
  assert (y.size()     == N);
  assert (sigma.size() == N);
  Model.SetParams(params);

  double chi2 = 0;
  for (int i=0; i<N; i++) {
    double val = Model (x(i), Params);
    chi2 += (val-y(i))*(val-y(i))/(sigma(i)*sigma(i));
  }
  return chi2;
}

template<int M, typename ModelType> void
NonlinearFitClass<M,ModelType>::CalcAlphaBeta (const Array<ArgType,1> &x, 
					       const Array<double,1> &y,
					       const Array<double,1> &sigma,
					       TinyVector<double,M> params)
{
  Model.SetParams(params);
  int N = x.size();
  assert (y.size()     == N);
  assert (sigma.size() == N);
  
  for (int k=0; k<M; k++) {
    Beta[k] = 0.0;
    fori (int l=0; l<M; l++)
      Alpha(k,l) = 0.0;
  }
  
  for (int i=0; i<N; i++) {
    double val = Model(x(i), params);
    TinyVector<double,M> grad = Model.Grad(x(i),params);
    for (int k=0; k<M; k++)
      Beta[k] += grad[k]*(y(i)-val)/(sigma(i)*sigma(i));
  }
  for (int i=0; i<N; i++) {
    double val = Model(x(i), params);
    TinyVector<double,M> grad = Model.Grad(x(i),params);
    for (int k=0; k<M; k++)
      for (int l=0; l<M; l++)
	Alpha(k,l) += grad[k]*grad[l]/(sigma(i)*sigma(i));
  }
}


template<int M, typename ModelType> void
NonlinearFitClass<M,ModelType>::Solve()
{
  for (int i=0; i<M; i++) {
    dParams[i] = 0.0;
    for (int j=0; j<M; j++)
      AlphaInv(i,j) = Alpha(i,j);
  }

  GJInverse(AlphaInv);

  for (int i=0; i<M; i++)
    for (int j=0; j<M; j++)
      dParams[i] += AlphaInv(i,j)*Beta(j);

  // Do simple Gaussian elimination
//   double dInv = 1.0/Alpha(0,0);
//   for (int col=0; col<M; col++)
//     Alpha(0,col) *= dInv;
//   Beta(0) *= dInv;
//   for (int row=1; row<M; row++) {
//     for (int col=0; col<row; col++) {
//       Beta(row) -= Alpha(row,col)*Beta(col);
//       for (int k=0;
//     double diag = Alpha(row,row);
//     double diagInv = 1.0/diag;
//     for (int col=0; col<M; col++)
//     Beta(row) *= diagInv;
//     for (int col=row+1; col<


}


template<int M, typename ModelType> void
NonlinearFitClass<M,ModelType>::Fit (const Array<ArgType,1> &x,
				     const Array<double,1> &y,
				     const Array<double,1> &sigma,
				     TinyVector<double,M> &params)
{
  double chiNow = Chi2 (x, y, sigma, params);
  double lambda = 0.001;
  
  bool done = false;
  while (!done) {
    CalcAlphaBeta (x, y, sigma, params);
    for (int i=0; i<M; i++)
      Alpha(i,i) *= (1.0+lambda);
    cerr << "params = " << params << endl;
    Solve();
    cerr << "dParams = " << dParams << endl;
    done = true;
    for (int i=0; i<M; i++)
      if (fabs(dParams[i]/params[i]) > 1.0e-5)
	done = false;
    if (!done) {
      double chiNew = Chi2 (x, y, sigma, params + dParams);
      if (chiNew < chiNow) {
	Lambda *= 5.0;
	params = params + dParams;
      }
      else {
	Lambda *= 0.2;
      }
    }
  }
}
