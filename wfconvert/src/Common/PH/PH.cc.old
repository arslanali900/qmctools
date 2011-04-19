#include "PH.h"
#include "../IO/InputFile.h"


PseudoHamiltonian *ReadPH (InputBuffer &SectionBuf)
{
  PseudoHamiltonian *PH;
  char Type[50];
  char FullCoreName[500];
  scalar CoreRadius;
  int Success;


  if (SectionBuf.ReadVar("Type", Type, 50))
    {
      if (!strcmp (Type, "CUBIC"))
	{
	  int NumAParams, NumBParams, NumVParams;
	  scalar Amin, Bmin, Vmin, Amax, Bmax, Vmax;
	  Array<scalar,1> Ainit, Binit, Vinit;
	  SectionBuf.Rewind();
	  if (!SectionBuf.ReadVar("FullCoreFileName", FullCoreName, 500))
	    Abort ("Cannot find variable FullCoreFileName in section PH.");
	  SectionBuf.Rewind();
	  if (!SectionBuf.ReadVar("CoreRadius", CoreRadius))
	    Abort ("Cannot find variable CoreRadius in section PH.");
	  SectionBuf.Rewind();
	  if (!SectionBuf.ReadVar ("Aparams", Ainit))
	    Abort ("Cannot find variable Aparams in section PH.");
	  SectionBuf.Rewind();
	  if (!SectionBuf.ReadVar ("Bparams", Binit))
	    Abort ("Cannot find variable Bparams in section PH.");
	  SectionBuf.Rewind();
	  if (!SectionBuf.ReadVar ("Vparams", Vinit))
	    Abort ("Cannot find variable Vparams in section PH.");
	  SectionBuf.Rewind();
	  if (!SectionBuf.ReadVar ("Amin", Amin))
	    Abort ("Cannot find variable Amin in section PH.");
	  SectionBuf.Rewind();
	  if (!SectionBuf.ReadVar ("Bmin", Bmin))
	    Abort ("Cannot find variable Bmin in section PH.");
	  SectionBuf.Rewind();
	  if (!SectionBuf.ReadVar ("Amax", Amax))
	    Abort ("Cannot find variable Amax in section PH.");
	  SectionBuf.Rewind();
	  if (!SectionBuf.ReadVar ("Bmax", Bmax))
	    Abort ("Cannot find variable Bmax in section PH.");
	  SectionBuf.Rewind();
	  if (!SectionBuf.ReadVar ("Vmin", Vmin))
	    Abort ("Cannot find variable Vmin in section PH.");
	  SectionBuf.Rewind();
	  if (!SectionBuf.ReadVar ("Vmax", Vmax))
	    Abort ("Cannot find variable Vmax in section PH.");
	  SectionBuf.Rewind();

	  FullCorePotential *FullCorePot = new FullCorePotential;
	  FullCorePot->Read(FullCoreName);
	  PH_CubicSpline *PHC;
	  PHC = new PH_CubicSpline;
	  PHC->Init (Amin, Bmin, Vmin, Amax, Bmax, Vmax,
		     Ainit, Binit, Vinit, FullCorePot,
		     CoreRadius);
	  PHC->Z = FullCorePot->Z;
	  PH = PHC;
	}
      else if (!strcmp (Type, "CUBICXC"))
	{
	  Array<scalar,1> Aparams, Bparams, Vparams;
	  Grid *Agrid, *Bgrid, *Vgrid;
	  scalar Zion;
	  SectionBuf.Rewind();
	  if (!SectionBuf.ReadVar("CoreRadius", CoreRadius))
	    Abort ("Cannot find variable CoreRadius in section PH.");
	  SectionBuf.Rewind();
	  if (!SectionBuf.ReadVar ("Aparams", Aparams))
	    Abort ("Cannot find variable Ainit in section PH.");
	  SectionBuf.Rewind();
	  if (!SectionBuf.ReadVar ("Bparams", Bparams))
	    Abort ("Cannot find variable Binit in section PH.");
	  SectionBuf.Rewind();
	  if (!SectionBuf.ReadVar ("Vparams", Vparams))
	    Abort ("Cannot find variable Vinit in section PH.");
	  SectionBuf.Rewind();
	  if (!SectionBuf.ReadVar ("Zion", Zion))
	    Abort ("Cannot find variable Amin in section PH.");
	  SectionBuf.Rewind();
	  InputBuffer GridBuf;
	  if (!SectionBuf.FindSection("Agrid", GridBuf))
	    Abort ("Cannot find section Agrid in section PH.");
	  Agrid = ReadGrid (GridBuf);
	  if (Agrid == NULL)
	    Abort ("Error in Agrid section of section PH.");

	  SectionBuf.Rewind();
	  if (!SectionBuf.FindSection("Bgrid", GridBuf))
	    Abort ("Cannot find section Bgrid in section PH.");
	  Bgrid = ReadGrid (GridBuf);
	  if (Bgrid == NULL)
	    Abort ("Error in Bgrid section of section PH.");
	  
	  SectionBuf.Rewind();
	  if (!SectionBuf.FindSection("Vgrid", GridBuf))
	    Abort ("Cannot find section Vgrid in section PH.");
	  Vgrid = ReadGrid (GridBuf);
	  if (Vgrid == NULL)
	    Abort ("Error in Vgrid section of section PH.");
	  
	  // DEBUG
	  cerr << "Vgrid.End = " << Vgrid->End << ".\n";

	  	  
	  PH_CubicSplineXC *PHXC;
	  PHXC = new PH_CubicSplineXC;
	  PHXC->Init (Aparams, Bparams, Vparams, 
		      Agrid, Bgrid, Vgrid, Zion,
		      CoreRadius);
	  PH = PHXC;
	}
      else if (!strcmp (Type, "NUCLEAR"))
	{
	  scalar Z;
	  SectionBuf.Rewind();
	  if (!SectionBuf.ReadVar ("Z", Z))
	    Abort ("Cannot find Z in Nuclear_PH section.\n");
	  PH_Nuclear *PHNuc;
	  PHNuc = new PH_Nuclear(Z);
	  PH = PHNuc;
	}
      else if (!strcmp (Type, "FULLCORE"))
	{
	  SectionBuf.Rewind();
	  if (!SectionBuf.ReadVar("FullCoreFileName", FullCoreName, 300))
	    Abort ("Cannot find variable FullCoreFileName in section PH.");
	  FullCorePotential *FullCorePot = new FullCorePotential;
	  FullCorePot->Read(FullCoreName);
	  PH_FullCore *PHFC;
	  PHFC = new PH_FullCore(FullCorePot);
	  PHFC->Zion = FullCorePot->Z;
	  PH = PHFC;
	}
      else if (!strcmp (Type, "CHEBYSHEV"))
	{
	  Abort ("CHEBYSHEV type not yet implemented.");
	}
      else
	Abort ("Unrecognized PH Type in ReadPH.");
    }
  else
    Abort ("Couldn't find Type variable in ReadPH.");

  return PH;
}


void
PH_Chebyshev::Write (char *FileName)
{
  FILE *fout;

  if ((fout = fopen(FileName, "w"))==NULL)
    {
      cerr << "Cannot open " << FileName << " for writing.\n";
      exit(1);
    }
  Write(fout);
  fclose(fout);
}


void
PH_Chebyshev::Write(FILE *fout)
{
  // First, write the core radius:
  fprintf (fout, "Chebyshev\n");
  fprintf (fout, "%1.16e\n", CoreRadius);
  // Next, write the FullCore filename
  fprintf (fout, "%s\n", FullCoreV->FileName);

  // Now write the parameters for the A potential
  fprintf (fout, "%d\n", Afunc->NumParams);
  fprintf (fout, "%1.16e\n", Afunc->AddConst);
  for (int i=0; i<Afunc->NumParams; i++)
    fprintf (fout, "%1.16e ", Afunc->Params(i));
  fprintf (fout, "\n");

  // Now write the parameters for the B potential
  fprintf (fout, "%d\n", Bfunc->NumParams);
  fprintf (fout, "%1.16e\n", Bfunc->AddConst);
  for (int i=0; i<Bfunc->NumParams; i++)
    fprintf (fout, "%1.16e ", Bfunc->Params(i));
  fprintf (fout, "\n");

  // Now write the parameters for the V potential
  fprintf (fout, "%d\n", PseudoCoreV->NumParams);
  for (int i=0; i<PseudoCoreV->NumParams; i++)
    fprintf (fout, "%1.16e ", PseudoCoreV->Params(i));
  fprintf (fout, "\n");

}



void
PH_CubicSpline::Write (char *FileName)
{
  FILE *fout;

  if ((fout = fopen(FileName, "w"))==NULL)
    {
      cerr << "Cannot open " << FileName << " for writing.\n";
      exit(1);
    }
  Write(fout);
  fclose (fout);
}

void
PH_CubicSpline::Write(FILE *fout)
{
  fprintf (fout, "PH\n{\n");
  fprintf (fout, "  Type = \"CUBIC\";\n");
  fprintf (fout, "  CoreRadius = %1.5f;\n",  CoreRadius);
  fprintf (fout, "  FullCoreFileName = \"%s\";\n", FullCoreV->FileName);
  fprintf (fout, "  Amin = %1.16e;\n", Amin);
  fprintf (fout, "  Bmin = %1.16e;\n", Bmin);
  fprintf (fout, "  Vmin = %1.16e;\n", Vmin);
  fprintf (fout, "  Amax = %1.16e;\n", Amax);
  fprintf (fout, "  Bmax = %1.16e;\n", Bmax);
  fprintf (fout, "  Vmax = %1.16e;\n", Vmax);
  
  fprintf (fout, "  Aparams = [ ");
  for (int i=0; i < PA.NumParams-1; i++)
    fprintf (fout, "%1.16e ", PA.Params(i));
  fprintf (fout, " ];\n  Bparams = [ ");
  for (int i=1; i < PB.NumParams-1; i++)
    fprintf (fout, "%1.16e ", PB.Params(i));
  fprintf (fout, " ];\n  Vparams = [ ");
  for (int i=0; i < Vfunc.NumParams-1; i++)
    fprintf (fout, "%1.16e ", Vfunc.Params(i));
  fprintf (fout, " ];\n}\n");
}

/*
void
PH_CubicSpline::Write(char *FileName)
{
  FILE *fout;

  if ((fout = fopen(FileName, "w"))==NULL)
    {
      cerr << "Cannot open " << FileName << " for writing.\n";
      exit(1);
    }

  // First, write the core radius:
  fprintf (fout, "CubicSpline\n");
  fprintf (fout, "%1.16e\n", CoreRadius);
  // Next, write the FullCore filename
  fprintf (fout, "%s\n", FullCoreV->FileName);

  // Now write the parameters for the A potential
  fprintf (fout, "%d\n", PA.NumParams-1);
  fprintf (fout, "%1.16e\n", Amin);
  for (int i=0; i<(PA.NumParams-1); i++)
    fprintf (fout, "%1.16e ", PA.Params(i));
  fprintf (fout, "\n");

  // Now write the parameters for the B potential
  fprintf (fout, "%d\n", PB.NumParams-2);
  fprintf (fout, "%1.16e\n", Bmin);
  for (int i=1; i<(PB.NumParams-1); i++)
    fprintf (fout, "%1.16e ", PB.Params(i));
  fprintf (fout, "\n");

  // Now write the parameters for the V potential
  fprintf (fout, "%d\n", Vfunc.NumParams-1);
  for (int i=0; i<(Vfunc.NumParams-1); i++)
    fprintf (fout, "%1.16e ", Vfunc.Params(i));
  fprintf (fout, "\n");

  fclose (fout);
}
*/


void
PH_Chebyshev::Read(FILE *fin)
{
  LocalAlloc = 1;
  Afunc = new AFunction;
  Bfunc = new BFunction;
  PseudoCoreV = new VFunction;
  FullCoreV = new FullCorePotential;

  Array<scalar,1> temp;

  // First, write the core radius:
  fscanf (fin, " %lf ", &CoreRadius);
  fscanf (fin, " %s ", FullCoreV->FileName);
  // Remove trailing newline from filename;
  FullCoreV->FileName[strlen(FullCoreV->FileName)-1] = '\0';
  FullCoreV->Read(FullCoreV->FileName);
  scalar Vmax  = (*FullCoreV)(CoreRadius);
  scalar dVmax = FullCoreV->Deriv(CoreRadius);

  fscanf (fin, " %d ", &(Afunc->NumParams));
  fscanf (fin, " %lf ", &(Afunc->AddConst));
  temp.resize(Afunc->NumParams);
  for (int i=0; i<Afunc->NumParams; i++)
    fscanf (fin, " %lf ", &(temp(i)));
  Afunc->Init(temp, CoreRadius);

  fscanf (fin, " %d ", &(Bfunc->NumParams));
  fscanf (fin, " %lf ", &(Bfunc->AddConst));
  temp.resize(Bfunc->NumParams);
  for (int i=0; i<Bfunc->NumParams; i++)
    fscanf (fin, " %lf ", &(temp(i)));
  Bfunc->Init(temp, CoreRadius, Afunc);
	
  fscanf (fin, " %d ", &(PseudoCoreV->NumParams));
  temp.resize(PseudoCoreV->NumParams);
  for (int i=0; i<PseudoCoreV->NumParams; i++)
    fscanf (fin, " %lf ", &(temp(i)));
  PseudoCoreV->Init(temp, CoreRadius, Vmax, dVmax);  

  fclose (fin);
} 




void
PH_CubicSpline::Read(FILE *fin)
{
  cerr << "Reading PH_CubicSpline.\n";
  FullCoreV = new FullCorePotential;
  
  Array<scalar,1> tempA, tempB, tempV;
  int NumA, NumB, NumV;
  
  // First, write the core radius:
  fscanf (fin, " %lf ", &CoreRadius);
  fscanf (fin, " %s ", FullCoreV->FileName);
  FullCoreV->Read(FullCoreV->FileName);
  scalar Vmax  = (*FullCoreV)(CoreRadius);
  scalar dVmax = FullCoreV->Deriv(CoreRadius);
  
  // Read A parameters
  fscanf (fin, " %d ", &NumA);
  fscanf (fin, " %lf ", &Amin);
  tempA.resize(NumA);
  for (int i=0; i<NumA; i++)
    fscanf (fin, " %lf ", &(tempA(i)));
  
  // Read B parameters
  fscanf (fin, " %d ", &NumB);
  fscanf (fin, " %lf ", &Bmin);
  tempB.resize(NumB);
  for (int i=0; i<NumB; i++)
    fscanf (fin, " %lf ", &(tempB(i)));
  
  // Read V parameters
  fscanf (fin, " %d ", &NumV);
  tempV.resize(NumV);
  for (int i=0; i<NumV; i++)
    fscanf (fin, " %lf ", &(tempV(i)));
  
  //Init (Amin, Bmin, tempA, tempB, tempV, FullCoreV, CoreRadius);

  Write("Test.PH");
  
  fclose (fin);
} 


void
PH_CubicSplineXC::Write (char *FileName)
{
  FILE *fout;

  if ((fout = fopen(FileName, "w"))==NULL)
    {
      cerr << "Cannot open " << FileName << " for writing.\n";
      exit(1);
    }
  Write(fout);
  fclose(fout);
}



void PH_CubicSplineXC::Write(FILE *fout)
{
  fprintf (fout, "PH\n{\n");
  fprintf (fout, "  Type = \"CUBICXC\";\n");
  fprintf (fout, "  CoreRadius = %1.16e;\n", CoreRadius);
  fprintf (fout, "  Zion = %1.16f;\n", Zion);

  fprintf (fout, "  Agrid\n");
  Agrid->WriteInput(fout);
  fprintf (fout, "  Bgrid\n");
  Bgrid->WriteInput(fout);
  fprintf (fout, "  Vgrid\n");
  Vgrid->WriteInput(fout);

  fprintf (fout, "  Aparams = [ %23.16e\n", PA.Params(0));
  for (int i=1; i<PA.NumParams-2; i++)
    fprintf (fout, "              %23.16e\n", PA.Params(i));
  fprintf (fout, "              %23.16e ];\n\n", PA.Params(PA.NumParams-2));
  
  fprintf (fout, "  Bparams = [ %23.16e\n", PB.Params(1));
  for (int i=2; i<PB.NumParams-2; i++)
    fprintf (fout, "              %23.16e\n", PB.Params(i));
  fprintf (fout, "              %23.16e ];\n\n", PB.Params(PB.NumParams-2));
  
  fprintf (fout, "  Vparams = [ %23.16e\n", Vfunc.Params(0));
  for (int i=1; i<Vfunc.NumParams-2; i++)
    fprintf (fout, "              %23.16e\n", Vfunc.Params(i));
  fprintf (fout, "              %23.16e ];\n\n", 
	   Vfunc.Params(Vfunc.NumParams-2));

  fprintf (fout, "}\n");

}


void
PH_CubicSplineXC::Read(FILE *fin)
{

}


void
PH_FullCore::Write (char *FileName)
{
  FILE *fout;

  if ((fout = fopen(FileName, "w"))==NULL)
    {
      cerr << "Cannot open " << FileName << " for writing.\n";
      exit(1);
    }
  Write(fout);
  fclose(fout);
}



void PH_FullCore::Write(FILE *fout)
{

}

void
PH_FullCore::Read(FILE *fin)
{

}


/* PseudoHamiltonian *Read_PH(char *FileName)
{
  cerr << "Got to Read_PH.\n";

  FILE *fin;

  if ((fin = fopen (FileName, "r")) == NULL)
    {
      cerr << "Can't open " << FileName << " for reading.  Exitting.\n";
      exit(1);
    }
  char PHtype[500];
  
  fscanf(fin, " %s ", PHtype);
  
  PseudoHamiltonian *PH;
  if (strcmp(PHtype, "Chebyshev") == 0)
    {
      PH = new PH_Chebyshev;
      PH->Read(fin);
    }
  else if (strcmp(PHtype, "CubicSpline") == 0)
    {
      PH = new PH_CubicSpline;
      PH->Read(fin);
    }
  else
    {
      fclose(fin);
      cerr << "Unrecognize PseudoHamiltonian type: " << PHtype <<
	"\n";
      cerr << "Exitting.\n";
      exit(1);
    }
  return (PH);
}*/









  
