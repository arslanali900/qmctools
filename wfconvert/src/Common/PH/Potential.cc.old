#include "PH.h"

Array<scalar, 2> Inverse (Array<scalar,2> A)
{
  Array<scalar,2> AInv(A.rows(), A.cols());
  AInv = 0.0;
  Array<scalar,1> v(A.cols());
  Array<scalar,1> v2(A.cols());

  for (int i=0; i<A.rows(); i++)
    AInv(i,i) = 1.0;

  for (int row=0; row<A.rows(); row++)
    {
      scalar temp = A(row,row);
      for (int col=0; col<A.cols(); col++)
	{
	  A(row,col) /= temp;
	  AInv(row,col) /= temp;
	}
      for (int k=0; k<A.rows(); k++)
	if (k!=row)
	  {
	    scalar temp = A(k,row);
	    v = temp  * A(row,Range::all());
	    v2 = temp * AInv(row,Range::all());
	    A(k,Range::all()) -=  v;
	    AInv(k,Range::all()) -= v2;
	  }
    }
  return (AInv);
}
	

Array<scalar,1> Prod(Array<scalar,2> A, Array<scalar,1>b)
{
  if (A.cols() != b.rows())
    {
      cerr << "Dimension mismatch if matrix-vector product.\n";
      exit(1);
    }
  Array<scalar,1> Ab(A.rows());

  Ab = 0.0;
  for (int row=0; row<Ab.rows(); row++)
    for (int k=0; k<A.cols(); k++)
      Ab(row) += A(row,k) * b(k);
  return (Ab);
}




//////////////////////////////////////////////////////////////////
// The three sets of routines that follow are for the storage   //
// classes for the A, B, and V potentials.  Strictly speaking   //
// only V is a potential, while A and B are inverse masses.     //
// Each of these potentials are represented in a Chebyshev      //
// polynomial expansion and are subject to different            //
// constraints.  These constraints produce linear relations     //
// between the polynomial coefficients, reducing the total      //
// number of required parameters.  The following function sets  //
// enforce the constraints on these potentials in a black-box   //
// manner, so that the calling routines need not worry about    //
// their enforcement.                                           //
//////////////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////////////
//                     AFunction Routines                       //
//////////////////////////////////////////////////////////////////

//  void 
//  AFunction::Update()
//  {
//    for (int i=2; i<(NumParams+2); i++)
//      ChebyCoefs(i) = Parameters(i-2);
//    ChebyCoefs(1) = 0.0;
//    for (int i=2; i<(NumParams+2); i++)
//      ChebyCoefs(1) -= (scalar) (i*i) * ChebyCoefs(i);
//    ChebyCoefs(0) = sqrt(1.0-Addconst);
//    for (int i=1; i<(NumParams+2); i++)
//      ChebyCoefs(0) -= ChebyCoefs(i);
//    ChebyCoefs(0) *= 2.0;
//    UpToDate = 1;
//  }


void 
AFunction::Update()
{
  const int NumConstraints = 3;
  // Solve constraint equations
  // Put them in the form Ax = b;
  Array<scalar,2> A(NumConstraints,NumConstraints);
  Array<scalar,1> b(NumConstraints);
  // indices holds which coeficients will be adjusted to enforce the
  // constraints 
  Array<int,1> indices(NumConstraints);
  
  int m = NumParams + NumConstraints-1;
  //indices(0) = m-2;
  //indices(1) = m-1;
  //indices(2) = m;

  indices(0) = 0;
  indices(1) = 1;
  indices(2) = 2;
 

  int index=0;
  for (int i=0; i<=m; i++)
    {
      int IsParam = 1;
      for (int k=0; k<NumConstraints; k++)
	if (indices(k) == i)
	  IsParam = 0;
      if (IsParam)
	{
	  ChebyCoefs(i) = Params(index);
	  index++;
	}
    }

 
  // Sum0: P(1) = sqrt(1-AddConst);
  b(0) = sqrt(1-AddConst);

  for (int j=0; j<=m; j++)
    {
      scalar coef;
      if (j==0)
	coef = 0.5;
      else
	coef = 1.0;
      
      int InSum=1;
      for (int k=0; k<NumConstraints; k++)
	{
	  if (j == indices(k))
	    {
	      A(0,k) = coef;
	      InSum=0;
	    }
	}
      if (InSum)
	b(0) -= coef * ChebyCoefs(j);
    }
// Sum1: P'(x=1) = 0;
  b(1) = 0.0;

  for (int j=0; j<=m; j++)
    {
      scalar coef;
      coef = (scalar)(j*j);
      
      int InSum=1;
      for (int k=0; k<NumConstraints; k++)
	{
	  if (j == indices(k))
	    {
	      A(1,k) = coef;
	      InSum=0;
	    }
	}
      if (InSum)
	b(1) -= coef * ChebyCoefs(j);
    }


// Sum2: P'(x=-1) = 0;
  b(2) = 0.0;
  for (int j=0; j<=m; j++)
    {
      scalar coef;
      coef = (scalar)(j*j);
      if (!(j%2))
	coef *= -1.0;
      
      int InSum=1;
      for (int k=0; k<NumConstraints; k++)
	{
	  if (j == indices(k))
	    {
	      A(2,k) = coef;
	      InSum=0;
	    }
	}
      if (InSum)
	b(2) -= coef * ChebyCoefs(j);
    }

  // Now invert A;
  Array<scalar,2> AInverse(NumConstraints, NumConstraints);

  AInverse = Inverse(A);
  Array<scalar,1> x(NumConstraints);

  x = Prod(AInverse,b);

  for (int k=0; k<NumConstraints; k++)
    ChebyCoefs(indices(k)) = x(k);

  UpToDate = 1;

//    for (int i=0; i<NumParams; i++)
//      ChebyCoefs(i) = Parameters(i);
  
//    scalar alpha = 0.0;
//    scalar beta = 0.0;
//    int m = NumParams+2;
//    scalar gamma = (scalar)((m-2)*(m-2));
//    scalar delta = (scalar)((m-1)*(m-1));
//    scalar sign = 1.0;
//    scalar minus1tom = 1.0;
//    if (m%2)
//      minus1tom = -1.0;
//    for (int j=1; j<NumParams; j++)
//      {
//        scalar j2 = (scalar) j*j;
//        alpha -= sign * j2 * ChebyCoefs(j);
//        beta -= j2 * ChebyCoefs(j);
//        sign *= -1.0;
//      } 
//    alpha *= minus1tom;
//    ChebyCoefs(m-2) = 1.0/(2*gamma) * (alpha+beta);
//    ChebyCoefs(m-1) = 1.0/(2.0*delta)* (beta-alpha);

//    ChebyCoefs(m) = sqrt(1.0-AddConst);
//    ChebyCoefs(m) -= 2.0*ChebyCoefs(0);
//    for (int i=1; i<m; i++)
//      ChebyCoefs(m) -= ChebyCoefs(i);
//    UpToDate = 1;
}


void
AFunction::Init(Array<scalar,1> params, scalar rmax)
{
  rA = rmax;
  NumParams = params.rows();
  Parameters.resize(NumParams);
  Parameters = params;
  ChebyCoefs.resize(NumParams+3);

  Update();
  //cerr << "A Coefs:\n";
  //for (int i=0; i<NumParams+2; i++)
  // fprintf (stderr, "%1.12f\n", 
}


AFunction::AFunction(Array<scalar,1> params, scalar rmax)
{
  Init(params, rmax);
}


scalar &
AFunction::Params(int i)
{
  UpToDate=0;
  return (Parameters(i));
}

scalar 
AFunction::Params(int i) const
{
  return (Parameters(i));
}

scalar 
AFunction::operator()(scalar r)
{
  scalar x = 2.0*r/rA - 1.0;

  if (!UpToDate)
    Update();

  scalar P = Chebyshev(ChebyCoefs, x);
  scalar val = P*P + AddConst;
  return (P*P + AddConst);
  //return (1.0/val);
}

scalar 
AFunction::Deriv(scalar r)
{
  scalar x = 2.0*r/rA - 1.0;
  
  if (!UpToDate)
    Update();
  
  Array<scalar,1> deriv(ChebyCoefs.rows()-1);

  deriv = ChebyshevDeriv(ChebyCoefs);

  scalar P = Chebyshev(ChebyCoefs, x);
  scalar dPdx = Chebyshev(deriv, x);
  scalar dxdr = 2.0/rA;

  scalar dAdr = 2.0*P*dPdx*dxdr;
  scalar val = P*P + AddConst;
  return (2.0*P*dPdx*dxdr);
  //return (-dAdr/(val*val));
}


scalar 
AFunction::Deriv2(scalar r)
{
  cerr << "Deriv2 not implemented for Chebychev AFunction.\n";
  exit (1);
  return (0.0);
  //return (-dAdr/(val*val));
}



//////////////////////////////////////////////////////////////////
//                     BFunction Routines                       //
//////////////////////////////////////////////////////////////////


void 
BFunction::Update()
{
  const int NumConstraints = 4;
  // Solve constraint equations
  // Put them in the form Ax = b;
  Array<scalar,2> A(NumConstraints,NumConstraints);
  Array<scalar,1> b(NumConstraints);
  // indices holds which coeficients will be adjusted to enforce the
  // constraints 
  Array<int,1> indices(NumConstraints);
  
  int m = NumParams + NumConstraints-1;
  //indices(0) = m-2;
  //indices(1) = m-1;
  //indices(2) = m;

  indices(0) = 0;
  indices(1) = 1;
  indices(2) = 2;
  indices(3) = 3;
 

  int index=0;
  for (int i=0; i<=m; i++)
    {
      int IsParam = 1;
      for (int k=0; k<NumConstraints; k++)
	if (indices(k) == i)
	  IsParam = 0;
      if (IsParam)
	{
	  ChebyCoefs(i) = Params(index);
	  index++;
	}
    }

 
  // Sum0: P(1) = sqrt(1-AddConst);
  b(0) = sqrt(1-AddConst);

  for (int j=0; j<=m; j++)
    {
      scalar coef;
      if (j==0)
	coef = 0.5;
      else
	coef = 1.0;
      
      int InSum=1;
      for (int k=0; k<NumConstraints; k++)
	{
	  if (j == indices(k))
	    {
	      A(0,k) = coef;
	      InSum=0;
	    }
	}
      if (InSum)
	b(0) -= coef * ChebyCoefs(j);
    }

  // Sum3: P(-1) = sqrt(A(-1)-AddConst);
  b(3) = sqrt((*Afunc)(0.0)-AddConst);

  for (int j=0; j<=m; j++)
    {
      scalar coef;
      if (j==0)
	coef = 0.5;
      else
	coef = 1.0;
      if (j%2)
	coef *= -1.0;
      
      int InSum=1;
      for (int k=0; k<NumConstraints; k++)
	{
	  if (j == indices(k))
	    {
	      A(3,k) = coef;
	      InSum=0;
	    }
	}
      if (InSum)
	b(3) -= coef * ChebyCoefs(j);
    }


// Sum1: P'(x=1) = 0;
  b(1) = 0.0;

  for (int j=0; j<=m; j++)
    {
      scalar coef;
      coef = (scalar)(j*j);
      
      int InSum=1;
      for (int k=0; k<NumConstraints; k++)
	{
	  if (j == indices(k))
	    {
	      A(1,k) = coef;
	      InSum=0;
	    }
	}
      if (InSum)
	b(1) -= coef * ChebyCoefs(j);
    }


// Sum2: P'(x=-1) = 0;
  b(2) = 0.0;
  for (int j=0; j<=m; j++)
    {
      scalar coef;
      coef = (scalar)(j*j);
      if (!(j%2))
	coef *= -1.0;
      
      int InSum=1;
      for (int k=0; k<NumConstraints; k++)
	{
	  if (j == indices(k))
	    {
	      A(2,k) = coef;
	      InSum=0;
	    }
	}
      if (InSum)
	b(2) -= coef * ChebyCoefs(j);
    }

  // Now invert A;
  Array<scalar,2> AInverse(NumConstraints, NumConstraints);

  AInverse = Inverse(A);
  Array<scalar,1> x(NumConstraints);

  x = Prod(AInverse,b);

  for (int k=0; k<NumConstraints; k++)
    ChebyCoefs(indices(k)) = x(k);

  UpToDate = 1;
}



//  void 
//  BFunction::Update()
//  {
//    for (int i=1; i<(NumParams+2); i++)
//      ChebyCoefs(i) = Parameters(i-2);
//    ChebyCoefs(1) = 0.0;
//    for (int i=2; i<(NumParams+2); i++)
//      ChebyCoefs(1) -= (scalar) (i*i) * ChebyCoefs(i);
//    ChebyCoefs(0) = sqrt(1.0-AddConst);
//    for (int i=1; i<(NumParams+2); i++)
//      ChebyCoefs(0) -= ChebyCoefs(i);
//    ChebyCoefs(0) *= 2.0;
//    UpToDate = 1;  
//  }


void
BFunction::Init(Array<scalar,1> params, scalar rmax, AFunction *A)
{
  rB = rmax;
  NumParams = params.rows();
  Parameters.resize(NumParams);
  Parameters = params;
  ChebyCoefs.resize(NumParams+4);
  Afunc = A;
  
  Update();
}


BFunction::BFunction(Array<scalar,1> params, scalar rmax, AFunction *A)
{
  Init(params, rmax, A);
}


scalar &
BFunction::Params(int i)
{
  UpToDate=0;
  return (Parameters(i));
}

scalar 
BFunction::Params(int i) const
{
  return (Parameters(i));
}

scalar 
BFunction::operator()(scalar r)
{
  scalar x = 2.0*r/rB - 1.0;

  if((!Afunc->UpToDate) || (!UpToDate))
    Update();

  scalar P = Chebyshev(ChebyCoefs, x);
  return (P*P + AddConst);
}





//////////////////////////////////////////////////////////////////
//                     VFunction Routines                       //
//////////////////////////////////////////////////////////////////

void 
VFunction::Update()
{
  UpToDate = 1;

  for (int i=2; i<(NumParams+2); i++)
    ChebyCoefs(i) = Parameters(i-2);
  ChebyCoefs(1) = dV_at_rV;
  //cerr << "dV_at_rV = " << dV_at_rV << "\n";
  //cerr << "V_at_rV = " << V_at_rV << "\n";
  for (int i=2; i<(NumParams+2); i++)
    ChebyCoefs(1) -= (scalar)(i*i) * ChebyCoefs(i);
  ChebyCoefs(0) = V_at_rV;
  for (int i=1; i<(NumParams+2); i++)
    ChebyCoefs(0) -= ChebyCoefs(i);
  ChebyCoefs(0) *= 2.0;
}


void
VFunction::Init(Array<scalar,1> params, scalar rmax, 
		scalar Vmax, scalar dVmax)
{
  rV = rmax;
  V_at_rV = Vmax;
  dV_at_rV = dVmax;
  NumParams = params.rows();
  Parameters.resize(NumParams);
  Parameters = params;
  ChebyCoefs.resize(NumParams+2);
  Update();
}


VFunction::VFunction(Array<scalar,1> params, scalar rmax,
		     scalar Vmax, scalar dVmax)
{
  Init(params, rmax, Vmax, dVmax);
}


scalar &
VFunction::Params(int i)
{
  UpToDate=0;
  return (Parameters(i));
}

scalar 
VFunction::Params(int i) const
{
  return (Parameters(i));
}

scalar 
VFunction::operator()(scalar r)
{
  scalar x = 2.0*r/rV - 1.0;

  if (!UpToDate)
    Update();

  scalar P = Chebyshev(ChebyCoefs, x);
  return (P);
}


scalar 
VFunction::Deriv(scalar r)
{
  scalar x = 2.0*r/rV - 1.0;
  
  if (!UpToDate)
    Update();
  
  Array<scalar,1> deriv(ChebyCoefs.rows()-1);

  deriv = ChebyshevDeriv(ChebyCoefs);

  scalar P = Chebyshev(ChebyCoefs, x);
  scalar dPdx = Chebyshev(deriv, x);
  scalar dxdr = 2.0/rV;

  return (dPdx*dxdr);
}



//////////////////////////////////////////////////////////////////
//                FullCorePotential Routines                    //
//////////////////////////////////////////////////////////////////

void
FullCorePotential::Read(char *FName)
{
  strncpy (FileName, FName, 1000);
  FILE *fin;
  if ((fin = fopen(FileName, "r")) == NULL)
    {
      cerr << "Cannot open " << FileName << " for reading. Exitting.";
      exit(1);
    }

  fscanf (fin, " %lf ", &Z);
  Grid *TempGrid;
  TempGrid = ReadGrid(fin);
  if (GridInitialized)
    delete(V.grid);
  GridInitialized = 1;
  scalar dummy;
  Array<scalar, 1> PotVals(TempGrid->NumPoints);
  for (int i=0; i<TempGrid->NumPoints; i++)
    fscanf(fin, " %lf %lf ", &PotVals(i), &dummy);
  V.Init(TempGrid, PotVals, 5.0e30, 5.0e30);

  fclose(fin);
}
