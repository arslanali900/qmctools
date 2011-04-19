#include "PlaneObject.h"
#include <sstream>
#include <gtkmm.h>

void
PlaneObject::SetCenter (Vec3 ucenter, Vec3 umin, Vec3 umax)
{
  uCenter = ucenter;
  uMin    = umin;
  uMax    = umax;
}

void
PlaneObject::SetPosition(int dir, double pos)
{
  Direction = dir;
  Position = pos;
  Set();
}

void
PlaneObject::SetIsocontours (bool active)
{
  UseContours = active;
}


void
PlaneObject::SetColorMap(ColorMapType map)
{
  MapType = map;
  CMap.Init (MinVal, MaxVal, map);
}

void
PlaneObject::Init()
{
  MinVal = Spline(0,0,0);
  MaxVal = Spline(0,0,0);
  for (int ix=0; ix<Spline.Nx; ix++) 
    for (int iy=0; iy<Spline.Ny; iy++)
      for (int iz=0; iz<Spline.Nz; iz++) {
	double val = Spline (ix,iy,iz);
	MinVal = (val < MinVal) ? val : MinVal;
	MaxVal = (val > MaxVal) ? val : MaxVal;
      }
  CMap. Init (MinVal, MaxVal, MapType);
  IsInitialized = true;
}

PlaneObject::~PlaneObject()
{
  if (HaveTexture) 
    glDeleteTextures(1, &TextureNum);
}

Vec3
PlaneObject::FindEdge(int is, int it, int edgeNum,
		      Vec3 u0, Vec3 s, Vec3 t,
		      double isoVal)
{
  double Ns = ValData.extent(0);
  double Nt = ValData.extent(1);
  int dim, is1, is2, it1, it2;
  dim = EdgeTable[edgeNum][0];
  is1 = is+EdgeTable[edgeNum][1];
  is2 = is+EdgeTable[edgeNum][2];
  it1 = it+EdgeTable[edgeNum][3];
  it2 = it+EdgeTable[edgeNum][4];

//   cerr << "is1=" << is1 << "  is2=" << is2  << "  "
//        << "it1=" << it1 << "  it2=" << it2 << endl;
//   cerr << "isoval=" << isoVal << endl;

//   cerr << "edgeNum = " << edgeNum << endl;
//   cerr << "is1 = " << is1 << endl;
//   cerr << "it1 = " << it1 << endl;
//   cerr << "is2 = " << is2 << endl;
//   cerr << "it2 = " << it2 << endl;

  Vec3 u1 = u0 + (double)(is1)/(Ns-1.0) * s + 
    (double)(it1)/(Nt-1.0) *t;
  Vec3 u2 = u0 + (double)(is2)/(Ns-1.0) * s + 
    (double)(it2)/(Nt-1.0) *t;
//   cerr << "is = " << is << endl;
//   cerr << "it = " << it << endl;
  double v1 = ValData (is1, it1) - isoVal;
  double v2 = ValData (is2, it2) - isoVal;
  //   cerr << "v1 = " << v1 << "   v2 = " << v2 << endl;
  

  if (v1*v2 >0.0)
    return (Vec3 (0.0, 0.0, 0.0));

  double frac = fabs(v1/(v2-v1));
  //  cerr << "frac = " << frac << endl;

  Vec3 u = (1.0-frac)*u1 + frac*u2;
  return u;
}


void
PlaneObject::Set()
{
  if (!IsInitialized)
    Init();
  if (!HaveTexture) {
    // Allocate a texture
    glGenTextures (1, &TextureNum);
    glBindTexture (GL_TEXTURE_2D, TextureNum);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    HaveTexture = true;
  }


  const int N = 256;
  if (ValData.size() != N*N) 
    ValData.resize(N,N);
  Array<TinyVector<GLubyte,4>,2> texData(N,N);
  
  Vec3 u0, sVec, tVec;
  Vec3 dr[3];
  dr[0] = (Spline.Xgrid->End - Spline.Xgrid->Start) * Vec3(1.0, 0.0, 0.0);
  dr[1] = (Spline.Ygrid->End - Spline.Ygrid->Start) * Vec3(0.0, 1.0, 0.0);
  dr[2] = (Spline.Zgrid->End - Spline.Zgrid->Start) * Vec3(0.0, 0.0, 1.0);
  u0 = Vec3 (Spline.Xgrid->Start, Spline.Ygrid->Start, Spline.Zgrid->Start);
  u0 += Position*dr[Direction];
  if (Direction == 0) {
    sVec = dr[1];
    tVec = dr[2];
  }
  else if (Direction == 1) {
    sVec = dr[2];
    tVec = dr[0];
  }
  else if (Direction == 2) {
    sVec = dr[0];
    tVec = dr[1];
  }

  // Create the data
  double scale = 1.0/(MaxVal - MinVal);
  double nInv = 1.0/(double)(N-1);
  for (int is=0; is<N; is++) {
    double s = 0.999999*nInv * (double)is;
    for (int it=0; it<N; it++) {
      double t = 0.999999*nInv * (double)it;
      Vec3 r = u0 + s*sVec + t*tVec;
      double val = Spline(r[0], r[1], r[2]);
      val = min(val, MaxVal);
      val = max(val, MinVal);
      
      ValData(is, it) = (val - MinVal)*scale;
      // cerr << "val=" << ValData(is,it) << "   r =" << r << endl;
      
      TinyVector<double,4> color;
      CMap (val, color);
      texData(is, it)[0] = (GLubyte) floor (256.0 * color[0]);
      texData(is, it)[1] = (GLubyte) floor (256.0 * color[1]);
      texData(is, it)[2] = (GLubyte) floor (256.0 * color[2]);
      texData(is, it)[3] = (GLubyte) floor (256.0 * color[3]);
    }
  }
  
  glBindTexture(GL_TEXTURE_2D, TextureNum);
  if (!BuiltTexture) {
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, N, N, 0, GL_RGBA,
		 GL_UNSIGNED_BYTE, texData.data());
    BuiltTexture = true;
  }
  else
    glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, N, N, 
		    GL_RGBA, GL_UNSIGNED_BYTE, texData.data());
  

  // Create the polygons
  Start();
  glEnable(GL_TEXTURE_2D);
  glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
  glBindTexture(GL_TEXTURE_2D, TextureNum);
  glBegin (GL_QUADS);
  Vec3 u1 = u0 + sVec;
  Vec3 u2 = u0 + sVec + tVec;
  Vec3 u3 = u0 + tVec;
  Vec3 u0save = u0;
  for (int i=0; i<3; i++) {
    u0[i] = uCenter[i] + (uMax[i]-uMin[i])*u0[i];
    u1[i] = uCenter[i] + (uMax[i]-uMin[i])*u1[i];
    u2[i] = uCenter[i] + (uMax[i]-uMin[i])*u2[i];
    u3[i] = uCenter[i] + (uMax[i]-uMin[i])*u3[i];
  }

  Vec3 r0 = u0*Lattice;
  Vec3 r1 = u1*Lattice;
  Vec3 r2 = u2*Lattice;
  Vec3 r3 = u3*Lattice;
  // When commensurate with the face of the cell, place slightly
  // outside the cell.
  if ((fabs(Position)<1.0e-6) || (fabs(1.0-Position)<1.0e-6)) {
    r0 = 1.0001*r0;  r1 = 1.0001*r1;  r2 = 1.0001*r2;  r3 = 1.0001*r3;
  }
  glTexCoord2f(0.0, 0.0); glVertex3f(r0[0], r0[1], r0[2]);
  glTexCoord2f(0.0, 1.0); glVertex3f(r1[0], r1[1], r1[2]);
  glTexCoord2f(1.0, 1.0); glVertex3f(r2[0], r2[1], r2[2]);
  glTexCoord2f(1.0, 0.0); glVertex3f(r3[0], r3[1], r3[2]);
  glEnd();
  glDisable(GL_TEXTURE_2D);

  // Now add contours
  if (UseContours) {
    int numContours = 20;
    glColor4d (1.0, 1.0, 1.0, 1.0);
    glBegin(GL_LINES);
    for (int cont=0; cont<numContours; cont++) {
      double isoVal = ((double)cont+0.5)/(double)(numContours+1);
      for (int is=0; is<(N-1); is++) {
	for (int it=0; it<(N-1); it++) {
	  int index = 0;
	  index |= ((ValData(is+0,it+0)> isoVal) << 3);
	  index |= ((ValData(is+1,it+0)> isoVal) << 2);
	  index |= ((ValData(is+1,it+1)> isoVal) << 1);
	  index |= ((ValData(is+0,it+1)> isoVal) << 0);
	  int ei=0;
	  int edge;
	  while ((edge=EdgeData[index][ei]) != -1) {
	    Vec3 reduced = 
	      FindEdge (is, it, edge, u0save, sVec, tVec, isoVal);
	    for (int i=0; i<3; i++) 
	      reduced[i] = uCenter[i] + (uMax[i]-uMin[i])*reduced[i];
	    Vec3 vertex = reduced * Lattice;
	    //  cerr << "vertex = " << vertex << endl;
	    glVertex3dv (&(vertex[0]));
	    ei++;
	  }
	}
      }
    }
  }
    glEnd();

  End();
}
	

void
PlaneObject::DrawPOV(FILE *fout, string rotMatrix)
{
  // First, we have to create an image file for the texture map
  const int N = 1024;
  Array<TinyVector<guint8,4>,2> texData(N,N);
  if (ValData.size() != N*N) 
    ValData.resize(N,N);

  
  Vec3 u0, sVec, tVec;
  Vec3 dr[3];
  dr[0] = (Spline.Xgrid->End - Spline.Xgrid->Start) * Vec3(1.0, 0.0, 0.0);
  dr[1] = (Spline.Ygrid->End - Spline.Ygrid->Start) * Vec3(0.0, 1.0, 0.0);
  dr[2] = (Spline.Zgrid->End - Spline.Zgrid->Start) * Vec3(0.0, 0.0, 1.0);
  u0 = Vec3 (Spline.Xgrid->Start, Spline.Ygrid->Start, Spline.Zgrid->Start);
  u0 += Position*dr[Direction];
  if (Direction == 0) {
    sVec = dr[1];
    tVec = dr[2];
  }
  else if (Direction == 1) {
    sVec = dr[2];
    tVec = dr[0];
  }
  else if (Direction == 2) {
    sVec = dr[0];
    tVec = dr[1];
  }

  // Create the data
  double scale = 1.0/(MaxVal - MinVal);
  double nInv = 1.0/(double)(N-1);
  for (int is=0; is<N; is++) {
    double s = nInv * (double)is;
    for (int it=0; it<N; it++) {
      double t = nInv * (double)it;
      Vec3 r = u0 + s*sVec + t*tVec;
      double val = Spline(r[0], r[1], r[2]);
      ValData(is, it) = (val - MinVal)*scale;
      TinyVector<double,4> color;
      CMap (val, color);
      texData(is, it)[0] = (guint8) floor (256.0 * color[0]);
      texData(is, it)[1] = (guint8) floor (256.0 * color[1]);
      texData(is, it)[2] = (guint8) floor (256.0 * color[2]);
      texData(is, it)[3] = (guint8) floor (256.0 * color[3]);
    }
  }

  // Old way: create pixbuf directly
  //   Glib::RefPtr<Gdk::Pixbuf> pixbuf = 
  //     Gdk::Pixbuf::create_from_data((guint8*)texData.data(), 
  // 				  Gdk::COLORSPACE_RGB,
  // 				  true, 8, N, N, 4*N);

  // New way:  create pixmap, draw to it, and then render to the pixbuf.
  Glib::RefPtr<Gdk::Drawable> drawable = (Glib::RefPtr<Gdk::Drawable>)NULL;
  Glib::RefPtr<Gdk::Pixmap> pixmap = 
    Gdk::Pixmap::create(drawable, N, N, 24);
  Glib::RefPtr<Gdk::Colormap> cmap = Gdk::Colormap::get_system();

  Glib::RefPtr<Gdk::GC> gc = 
    Gdk::GC::create((Glib::RefPtr<Gdk::Drawable>)pixmap);

  pixmap->draw_rgb_32_image (gc, 0, 0, N, N, Gdk::RGB_DITHER_NONE, 
			     (guchar*)texData.data(), 4*N);

  if (UseContours) {
    // Create cairo context
    Cairo::RefPtr<Cairo::Context> context = pixmap->create_cairo_context();
    context->set_line_width(1.5);
    context->set_source_rgb(1.0, 1.0, 1.0);
    // Now add contours
    int numContours = 20;
    bool close = false;
    int numseg = 0;
    for (int cont=0; cont<numContours; cont++) {
      double isoVal = ((double)cont+0.5)/(double)(numContours+1);
      for (int is=0; is<(N-1); is++) {
	for (int it=0; it<(N-1); it++) {
	  int index = 0;
	  index |= ((ValData(is+0,it+0)> isoVal) << 3);
	  index |= ((ValData(is+1,it+0)> isoVal) << 2);
	  index |= ((ValData(is+1,it+1)> isoVal) << 1);
	  index |= ((ValData(is+0,it+1)> isoVal) << 0);
	  int ei=0;
	  int edge;
	  while ((edge=EdgeData[index][ei]) != -1) {
	    Vec3 reduced = 
	      FindEdge (is, it, edge, u0, sVec, tVec, isoVal);
	    reduced -= u0;
	    double x = (double)N*dot (tVec, reduced);
	    double y = (double)N*dot (sVec, reduced);
	    // x = max(min(x, 1022.0), 1.0);
	    // y = max(min(y, 1022.0), 1.0);
	    //	    cerr << "Before drawing (" << x << ", " << y << ").\n";
	    if (close) {
	      context->line_to (x, y);
	      numseg++;
	      if (numseg > 1000) {
		context->stroke();
		numseg = 0;
	      }
	    }
	    else
	      context->move_to (x, y);
	    //	  cerr << "After drawing.\n";
	    close = !close;
	    ei++;
	  }
	}
      }
    }
    if (numseg)
      context->stroke();
  }
		      
  //  cerr << "Before copying into pixbuf.\n";
  // Copy into a pixbuf
  Glib::RefPtr<Gdk::Pixbuf> pixbuf = 
  Gdk::Pixbuf::create ((Glib::RefPtr<Gdk::Drawable>)pixmap, cmap, 
		       0, 0, 0, 0, N, N);

  // Write to a file
  stringstream fname;
  fname << "ColorPlane" << Direction << ".png";
  //  cerr << "Before writing .png file.\n";
  pixbuf->save (fname.str(), (Glib::ustring)"png"); 

  // Now, we actually create a box geometry onto which to apply this
  // texture.  
//   fprintf (fout, "intersection {\n");
//   fprintf (fout, "  box {\n");
//   fprintf (fout, "    <%10.8f, %10.8f, %10.8f>,\n",
// 	   Spline.Xgrid->Start, Spline.Ygrid->Start, Spline.Zgrid->Start);
//   fprintf (fout, "    <%10.8f, %10.8f, %10.8f>\n",
//  	   Spline.Xgrid->End, Spline.Ygrid->End, Spline.Zgrid->End);
//   fprintf (fout, "%s", rotMatrix.c_str());
//   fprintf (fout, "  }\n");
  fprintf (fout, "box {\n");
  fprintf (fout, "  <0.0000,0.0000,0.0000>, \n");
  fprintf (fout, "  <1.0000,1.0000,0.0001>  \n");
  fprintf (fout, "  pigment {\n");
  fprintf (fout, "    image_map {\"%s\"}\n", fname.str().c_str());
  fprintf (fout, "  }\n"); // pigment
  fprintf (fout, "  translate <-0.5, -0.5, 0.0>\n");
  if (Direction == 0)
    fprintf (fout, "  scale <-1, -1, 1>\n");
  else if (Direction == 1)
    fprintf (fout, "  scale <1, -1, 1>\n");
  else
    fprintf (fout, "  scale <1, 1, 1>\n");
  // Now rotate
  if (Direction == 0)
    fprintf (fout, "  rotate <0, 90, 0>\n");
  else if (Direction == 1)
    fprintf (fout, "  rotate <90, 0, 0>\n");
  else if (Direction == 2)
    fprintf (fout, "  rotate <0, 0, 90>\n");

  // First scale
  fprintf (fout, "  scale <%1.6f, %1.6f, %1.6f>\n",
	   (Spline.Xgrid->End-Spline.Xgrid->Start),
	   (Spline.Ygrid->End-Spline.Ygrid->Start),
	   (Spline.Zgrid->End-Spline.Zgrid->Start));
  fprintf (fout, "  matrix <%12.8f, %12.8f, %12.8f, \n",
	   Lattice(0,0), Lattice(0,1), Lattice(0,2));
  fprintf (fout, "          %12.8f, %12.8f, %12.8f, \n",
	   Lattice(1,0), Lattice(1,1), Lattice(1,2));
  fprintf (fout, "          %12.8f, %12.8f, %12.8f, \n",
	   Lattice(2,0), Lattice(2,1), Lattice(2,2));
  fprintf (fout, "          %12.8f, %12.8f, %12.8f> \n",
	   0.0, 0.0, 0.0);

  fprintf (fout, "  translate <%1.6f, %1.6f, %1.6f>\n",
	   (Position-0.5) * Lattice(Direction,0), 
	   (Position-0.5) * Lattice(Direction,1), 
	   (Position-0.5) * Lattice(Direction,2));
  

  
//   if (Direction == 0) 
//     fprintf (fout, "  translate <%1.6f, %1.6f, %1.6f>\n",
//  	     (1.0-Position)*Spline.Xgrid->Start
//  	     +Position*Spline.Xgrid->End, 0.0, 0.0);
//   else if (Direction == 1)
//      fprintf (fout, "  translate <%1.6f, %1.6f, %1.6f>\n",
//  	     0.0, 
//  	     (1.0-Position)*Spline.Ygrid->Start
//  	     +Position*Spline.Ygrid->End, 0.0);
//   else if (Direction == 2)
//     fprintf (fout, "  translate <%1.6f, %1.6f, %1.6f>\n",
//  	     0.0, 0.0, (1.0-Position)*Spline.Zgrid->Start
//  	     + Position*Spline.Zgrid->End);
  


  // Now translate
  //   fprintf (fout, "  translate <%1.6f, %1.6f, %1.6f>\n",
  // 	   Spline.Xgrid->Start, Spline.Ygrid->Start, Spline.Zgrid->Start);

  fprintf (fout, "%s", rotMatrix.c_str());
  fprintf (fout, "  }\n"); // box
  //  fprintf (fout, "}\n");   // intersection
}


void
PlaneObject::SetLattice (Mat3 lattice)
{
  Lattice = lattice;
}


int PlaneObject::EdgeTable[4][5] = {
//  dim x1 x2 y1 y2
  {  0,  0, 1, 0, 0 },
  {  1,  1, 1, 0, 1 },
  {  0,  0, 1, 1, 1 },
  {  1,  0, 0, 0, 1 }
};

int PlaneObject::EdgeData[16][5]=
{
  {-1,-1,-1,-1,-1},
  { 2, 3,-1,-1,-1},
  { 1, 2,-1,-1,-1},
  { 1, 3,-1,-1,-1},
  { 0, 1,-1,-1,-1},
  { 3, 0, 1, 3,-1},
  { 0, 2,-1,-1,-1},
  { 3, 0,-1,-1,-1},
  { 3, 0,-1,-1,-1},
  { 0, 2,-1,-1,-1},
  { 2, 3, 0, 1,-1},
  { 0, 1,-1,-1,-1},
  { 1, 3,-1,-1,-1},
  { 1, 2,-1,-1,-1},
  { 2, 3,-1,-1,-1},
  {-1,-1,-1,-1,-1}
};
