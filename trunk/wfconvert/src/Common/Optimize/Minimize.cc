/////////////////////////////////////////////////////////////
// Copyright (C) 2003-2006 Bryan Clark and Kenneth Esler   //
//                                                         //
// This program is free software; you can redistribute it  //
// and/or modify it under the terms of the GNU General     //
// Public License as published by the Free Software        //
// Foundation; either version 2 of the License, or         //
// (at your option) any later version.  This program is    //
// distributed in the hope that it will be useful, but     //
// WITHOUT ANY WARRANTY; without even the implied warranty //
// of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. //  
// See the GNU General Public License for more details.    //
// For more information, please see the PIMC++ Home Page:  //
//           http://pathintegrals.info                     //
/////////////////////////////////////////////////////////////

#include "Minimize.h"
#include <gsl/gsl_multimin.h>

double MCost (const gsl_vector *v, void *params)
{
  ConjugateGradient &CG = *(ConjugateGradient *)params;
  MinimizeFunction &MF = *CG.MinFunc;
  for (int i=0; i<MF.NumParams(); i++)
    MF.Params(i) = gsl_vector_get(v, i);
  return (MF.Cost());
}


void MGradCost(const gsl_vector *v, void *params, gsl_vector *df)
{
  ConjugateGradient &CG = *(ConjugateGradient *)params;
  MinimizeFunction &MF = *CG.MinFunc;
  double epsilon = CG.epsilon;  
  Array<scalar,1> gradient(MF.NumParams());
  for (int i=0; i<MF.NumParams(); i++)
    {
      for (int j=0; j<MF.NumParams(); j++)
	MF.Params(j) = gsl_vector_get(v,j);
      scalar CostPlus, CostMinus;
      MF.Params(i) = gsl_vector_get(v,i) + epsilon;
      CostPlus = MF.Cost();
      MF.Params(i) = gsl_vector_get(v,i) - epsilon;
      CostMinus = MF.Cost();
      gradient(i) = (CostPlus-CostMinus)/(2.0*epsilon);
      gsl_vector_set(df, i, (CostPlus-CostMinus)/(2.0*epsilon));
    }
  fprintf (stderr, "Gradient = \n");
  for (int i=0;i<MF.NumParams(); i++)
    fprintf (stderr, "%1.14e\n", gradient(i));
}

void MBoth (const gsl_vector *v, void *params, double *f,
		gsl_vector *df)
{
  *f = MCost(v, params);
  MGradCost(v, params, df);
}


void
ConjugateGradient::Minimize(MinimizeFunction &MinimFunc)
{
  size_t iter = 0;
  int status;

  MinFunc = &MinimFunc;
  const gsl_multimin_fdfminimizer_type *T;
  gsl_multimin_fdfminimizer *Minimizer;

  gsl_multimin_function_fdf gslMinFunc;
 
  gslMinFunc.f = &MCost;
  gslMinFunc.df = &MGradCost;
  gslMinFunc.fdf = &MBoth;
  gslMinFunc.n = MinimFunc.NumParams();
  gslMinFunc.params = this;

  gsl_vector *StartPoint;
  StartPoint = gsl_vector_alloc(MinimFunc.NumParams());
  for (int i=0; i<MinimFunc.NumParams(); i++)
    gsl_vector_set(StartPoint, i, MinimFunc.Params(i));

  T = gsl_multimin_fdfminimizer_conjugate_fr;
  Minimizer = gsl_multimin_fdfminimizer_alloc(T,MinimFunc.NumParams());

  gsl_multimin_fdfminimizer_set(Minimizer, &gslMinFunc, StartPoint,
				StepSize, Tolerance);
  
  do
    {
      iter++;
      status = gsl_multimin_fdfminimizer_iterate(Minimizer);
      
      if (status)
	break;
      
      status = gsl_multimin_test_gradient (Minimizer->gradient,
					   Tolerance);

      if (status == GSL_SUCCESS)
	cerr << "Minimum found at:\n";

      fprintf (stderr, "%5d ", iter);
      for (int i=0; i<MinimFunc.NumParams(); i++)
	fprintf (stderr, "%15.12e ", gsl_vector_get(Minimizer->x, i));
      fprintf (stderr, "\n");
      for (int i=0; i<MinimFunc.NumParams(); i++)
	MinimFunc.Params(i) = gsl_vector_get(Minimizer->x,i);
      MinimFunc.WriteStuff();
    }
  while (status == GSL_CONTINUE && iter < 100);

  for (int i=0; i<MinimFunc.NumParams(); i++)
    MinimFunc.Params(i) = gsl_vector_get(Minimizer->x,i);

  gsl_multimin_fdfminimizer_free (Minimizer);
  gsl_vector_free (StartPoint);

}






void
MinimizeFunction::ReadParameters(char *FileName)
{
  FILE *fin;
  if ((fin = fopen (FileName, "r")) == NULL)
    {
      cerr << "Can't open parmeters file.  Exitting.\n";
      exit(1);
    }

  for (int i=0; i<NumParams(); i++)
    {
      scalar temp;
      fscanf (fin, " %lf ", &temp);
      Params(i) = temp;
    }
}


class TestMinim : public MinimizeFunction
{
public:
  Array<scalar, 1> x;
  int NumParams()
  {
    return (2);
  }
  scalar &Params(int i)
  {
    return(x(i));
  }
  scalar Params(int i) const
  {
    return (x(i));
  }
  scalar Cost()
  {
    return ((x(0)-1.0)*(x(0)-1.0) + (x(1)-2.0)*(x(1)-2.0));
  }
};



//  main()
//  {
//    TestMinim Minim;

//    Minim.Tolerance = 1.0e-10;
//    Minim.StepSize = 0.01;
//    Minim.epsilon = 1.0e-6;
//    Minim.x.resize(2);
//    Minim.x(0) = 3.19;
//    Minim.x(1) = -2.56;

//    Minim.Minimize();
//  }
