#ifndef RESPLINE_H
#define RESPLINE_H

#include "config.h"
#include <einspline/bspline.h>
#include <einspline/nubspline.h>
#include <Common/IO/IO.h>
#include <Common/Blitz.h>
#include <Common/MPI/Communication.h>

class OrbSplineClass
{
private:
  UBspline_3d_z  *OldSpline, *Uniform;
  NUBspline_3d_z *Nonuniform;
  NUgrid *xNUgrid, *yNUgrid, *zNUgrid;
  int Nx, Ny, Nz;
  Array<complex<double>,3> Data, NonuniformData;
  double ClusterFactor;
  // Other data we need to copy
  double Eigenvalue, Radius;
  int Spin;
  Vec3 uMin, uMax, Center;
  bool Truncated;
  Array<Vec3,1> Centers;
public:
  void Read  (IO::IOSectionClass &in);
  void Write (IO::IOSectionClass &out);
  void Respline (int nx, int ny, int nz, double clusterFactor);
  void Respline (int nx, int ny, int nz);
  void GetDims (int &nx, int &ny, int &nz);
  void Plot (string filename);
  double UniformError(int numPoints);
  double NonuniformError(int numPoints);
  OrbSplineClass() : OldSpline(NULL), Uniform(NULL), Nonuniform(NULL),
		     xNUgrid(NULL), yNUgrid(NULL), zNUgrid(NULL)
  {

  }
};


class ResplineClass
{
private:
  IO::IOSectionClass In, Out;
  double SizeFactor;
  void CopyHeader();
  OrbSplineClass Orb;
  CommunicatorClass Comm;
  Array<double,2> ClusterFactor, RMSError;
  void FindOptimalClusterFactor (int twist, int band);
  void ProcTasks (int procNum, int numTasks,
		  int &firstTask, int &lastTask);
public:
  void SetFiles (string inName, string outName);
  void SetSizeFactor (double factor);
  void DoRespline();
  void Run(int orbNum);
  void RunTest(int orbNum);
  void PlotOrb (int orbNum);
  ResplineClass() 
  {
    Comm.SetWorld();
  }
};


#endif
