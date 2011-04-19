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

#ifndef GRACE_PLOT_H
#define GRACE_PLOT_H
#include "../Blitz.h"
#include <grace_np.h>
#include <list>


class GraceCurve
{
private:

public:
  int PlotNum, CurveNum;
  void SetLineWidth(double width);
  void SetColor(string color);
  void SetData(const Array<double,1> &x, const Array<double,1> &y);
  void Redraw();
};

class GracePlot
{
private:
  static list<int> Plots;
  int MyPlot, MyNumCurves;
  void SendCommand (stringstream &cmd);
public:
  void MakeCurrent();
  Array<GraceCurve, 1> Curves;
  inline int NumCurves() { return MyNumCurves; }
  void AddCurve();
  void SetTitle(string title);
  void Redraw();
  void SetBounds (double xmin, double ymin, double xmax, double ymax);

  GracePlot();
  ~GracePlot();
};

#endif
