#include "ColorMap.h"
#include <Common/IO/IO.h>
#include <gtkmm.h>

Array<double,3> MapData;

string FindFullPath(string filename)
{
  string fullName;

  fullName = filename;
  if (Glib::file_test(fullName, Glib::FILE_TEST_EXISTS))
    return fullName;
  else {
    fullName = PKG_DATA_DIR + filename;
    if (Glib::file_test(fullName, Glib::FILE_TEST_EXISTS))
      return fullName;
    else {
      cerr << "Cannot find \"" << filename << "\" anywhere.\n";
      return filename;
    }
  }
}

bool ColorMap::MapsRead = false;

void
ColorMap::ReadMaps()
{
  IO::IOSectionClass in;
  if (in.OpenFile (FindFullPath ("colormaps.in"))) {
    MapData.resize(WINTER+1, 64, 3);
    Array<double,2> map;
    in.ReadVar ("autumn", map); 
    MapData(AUTUMN,Range::all(),Range::all()) = map;
    in.ReadVar ("bone", map); 
    MapData(BONE,  Range::all(),Range::all()) = map;
    in.ReadVar ("colorcube", map);
    MapData(COLORCUBE,  Range::all(),Range::all()) = map;
    in.ReadVar ("cool", map);
    MapData(COOL,  Range::all(),Range::all()) = map;
    in.ReadVar ("copper", map);
    MapData(COPPER,  Range::all(),Range::all()) = map;
    in.ReadVar ("flag", map);
    MapData(FLAG,  Range::all(),Range::all()) = map;
    in.ReadVar ("gray", map);
    MapData(GRAY,  Range::all(),Range::all()) = map;
    in.ReadVar ("hot", map);
    MapData(HOT,  Range::all(),Range::all()) = map;
    in.ReadVar ("hsv", map);
    MapData(HSV,  Range::all(),Range::all()) = map;
    in.ReadVar ("jet", map);
    MapData(JET,  Range::all(),Range::all()) = map;
    in.ReadVar ("lines", map);
    MapData(LINES,  Range::all(),Range::all()) = map;
    in.ReadVar ("pink", map);
    MapData(PINK,  Range::all(),Range::all()) = map;
    in.ReadVar ("prism", map);
    MapData(PRISM,  Range::all(),Range::all()) = map;
    in.ReadVar ("spring", map);
    MapData(SPRING,  Range::all(),Range::all()) = map;
    in.ReadVar ("summer", map);
    MapData(SUMMER,  Range::all(),Range::all()) = map;
    in.ReadVar ("white", map);
    MapData(WHITE,  Range::all(),Range::all()) = map;
    in.ReadVar ("winter", map);
    MapData(WINTER,  Range::all(),Range::all()) = map;
    
    in.CloseFile();
    MapsRead = true;
  }
}

void
ColorMap::Init (double min, double max, ColorMapType map)
{
  Min = 0.9999*min;
  Max = 0.9999*max;

  Array<double,1> r, g, b, a;
  if (map <= WINTER) {
    r.resize(64); g.resize(64); b.resize(64); a.resize(64);
    r = 0.99*MapData(map,Range::all(), 0);
    g = 0.99*MapData(map,Range::all(), 1);
    b = 0.99*MapData(map,Range::all(), 2);
    a = 0.99;
  }
  else if (map == BLUE_WHITE_RED) {
    r.resize(10); g.resize(10); b.resize(10); a.resize(10);
    r = 0.00, 0.00, 0.00, 0.10, 0.50, 0.80, 0.85, 0.85, 0.70, 0.50;
    g = 0.00, 0.05, 0.20, 0.50, 0.85, 0.85, 0.50, 0.20, 0.10, 0.05;
    b = 0.50, 0.80, 0.85, 0.85, 0.80, 0.50, 0.20, 0.05, 0.00, 0.00;
    a = 0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99;
  }
  else if (map == GRAY) {
    r.resize(5); g.resize(5); b.resize(5); a.resize(5);
    r = 0.00, 0.25, 0.50, 0.75, 1.00;
    g = 0.00, 0.25, 0.50, 0.75, 1.00;
    b = 0.00, 0.25, 0.50, 0.75, 1.00;
    a = 0.99, 0.99, 0.99, 0.99, 0.99;
  }
  else {
    cerr << "Unknown map type in ColorMap::Init.\n";
    abort();
  }
  BoundaryCondition<double> fBC(FLAT), nBC(NATURAL);
  Splines[0].Init (min, max, r, true, fBC, nBC);
  Splines[1].Init (min, max, g, true, fBC, fBC);
  Splines[2].Init (min, max, b, true, nBC, fBC);
  Splines[3].Init (min, max, a, true, fBC, fBC);
  Initialized = true;
//   FILE *fout = fopen ("colormap.dat", "w");
//   for (double x=min; x<=max; x+=0.001)
//     fprintf (fout, "%10.6e %10.6e %10.6e %10.6e %10.6e\n",
// 	     x, Splines[0](x), Splines[1](x), Splines[2](x), Splines[3](x));
//   fclose (fout);
}
