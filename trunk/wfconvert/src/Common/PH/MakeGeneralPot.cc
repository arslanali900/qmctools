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

#include "../IO/IO.h"
#include <vector>
using namespace IO;
using namespace std;

main(int argc, char **argv)
{
  if (argc < 3) {
    cerr << "Usage:\n   MakeGeneralPot pot.dat outfile.h5\n";
    exit(1);
  }
  IOSectionClass out;
  out.NewFile ((string)argv[2]);
  std::vector<double> rvec, Vvec;
  Array<double,1> grid, V;
  
  FILE *fin;
  assert ((fin = fopen (argv[1], "r"))!= NULL);
  double r, c, v;
  bool done = false;
  while (!done) {
    done = (fscanf (fin, "%lf %lf %lf\n", &r, &c, &v) == EOF);
    if (!done) {
      rvec.push_back(r);
      Vvec.push_back(v);
    }
  }
  fclose (fin);
  int N = rvec.size();
  grid.resize(N);
  V.resize(N);
  out.NewSection("Grid");
  for (int i=0; i<N; i++) {
    grid(i) = rvec[i];
    V(i) = Vvec[i];
  }
  double Z = -round(V(N-1)*grid(N-1));
  out.WriteVar ("Points", grid);
  out.WriteVar ("Type", string("General"));
  out.CloseSection();
  out.WriteVar ("Type", string("General"));
  out.WriteVar("V", V);
  out.WriteVar("Z", Z);
  out.CloseFile();
}
