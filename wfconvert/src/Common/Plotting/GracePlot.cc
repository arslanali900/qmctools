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

#include "GracePlot.h"

list<int> GracePlot::Plots;

void GraceCurve::SetData(const Array<double,1> &x, 
			 const Array<double,1> &y)
{
  GracePrintf ("kill g%d.s%d", PlotNum, CurveNum);
  GracePrintf ("g%d.s%d on", PlotNum, CurveNum);
  GracePrintf ("focus g%d", PlotNum);
  for (int i=0; i<x.size(); i++) 
    GracePrintf ("g%d.s%d point %1.8e, %1.8e", PlotNum, CurveNum,
		 x(i), y(i));
  GracePrintf ("with g%d; autoscale", PlotNum);
  GraceCommand("redraw");
}

void GraceCurve::SetColor (string color)
{
  GracePrintf ("g%d.s%d color \"%s\"", PlotNum, CurveNum, color.c_str());
}

void GraceCurve::SetLineWidth(double width)
{
  GracePrintf ("g%d.s%d linewidth %1.2f", PlotNum, CurveNum, width);
}



GracePlot::GracePlot() : MyNumCurves(0)
{
  if (!GraceIsOpen()) {
    if (GraceOpen (4096) == -1)
      cerr << "Can't run Grace.\n";
  }
  if (Plots.empty()) 
    MyPlot = 0;
  else
    MyPlot = Plots.back()+1;
  Plots.push_back(MyPlot);
  MakeCurrent();
}

GracePlot::~GracePlot()
{
  GracePrintf ("kill g%d", MyPlot);
  Plots.remove(MyPlot);
  if (Plots.empty())
    GraceClose();
}

void GracePlot::SetBounds(double xmin, double ymin, double xmax, double ymax)
{
  MakeCurrent();
  GracePrintf ("world xmin %1.8e", xmin);
  GracePrintf ("world xmax %1.8e", xmax);
  GracePrintf ("world ymin %1.8e", ymin);
  GracePrintf ("world ymax %1.8e", ymax);
//   GracePrintf ("view xmin %1.8e", xmin);
//   GracePrintf ("view xmax %1.8e", xmax);
//   GracePrintf ("view ymin %1.8e", ymin);
//   GracePrintf ("view ymax %1.8e", ymax);
}

void GracePlot::Redraw()
{
  GracePrintf ("redraw");
}

void GracePlot::SendCommand (stringstream &cmd)
{
  string s;
  cmd >> s;
  GracePrintf (s.c_str());
}

void GracePlot::MakeCurrent()
{
  GracePrintf ("with g%d", MyPlot);
  list<int>::iterator iter;
  for (iter=Plots.begin(); iter != Plots.end(); iter++)
    GracePrintf("g%d hidden true", *iter);
  GracePrintf("g%d on", MyPlot);
  GracePrintf("focus g%d", MyPlot);
}

void GracePlot::SetTitle (string title)
{
  MakeCurrent();
  stringstream cmd;
  cmd << "TITLE " << title;
  SendCommand(cmd);
}


void GracePlot::AddCurve ()
{
  MyNumCurves++;
  Curves.resizeAndPreserve(MyNumCurves);
  Curves(MyNumCurves-1).CurveNum = MyNumCurves-1;
  Curves(MyNumCurves-1).PlotNum = MyPlot;
}
