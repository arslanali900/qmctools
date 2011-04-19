#include "SphereObject.h"
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>

int  SphereObject::SphereListNum(0);
bool SphereObject::SphereListCreated(false);
int  SphereObject::OffScreenListNum(0);
bool SphereObject::OffScreenListCreated(false);

SphereObject::SphereObject(bool offScreen) : 
  Radius(1.0), Color(1.0, 0.0, 0.0, 1.0), Pos(0.0, 0.0, 0.0),
  OffScreen(offScreen), NumFacets (20)
{
  if (offScreen) {
    if (!OffScreenListCreated) {
      glEnable(GL_NORMALIZE);
//       GLUquadricObj* qobj = gluNewQuadric();
//       gluQuadricDrawStyle(qobj, GLU_FILL);

//       OffScreenListNum = glGenLists(1);
//       glNewList (OffScreenListNum, GL_COMPILE);
//       gluSphere(qobj, 1.0, 30, 30);
//       gluDeleteQuadric(qobj);
//       glEndList();
//       OffScreenListCreated = true;
    }
  }
  else if (!SphereListCreated) {
    glEnable(GL_NORMALIZE);
    GLUquadricObj* qobj = gluNewQuadric();
    gluQuadricDrawStyle(qobj, GLU_FILL);
    SphereListNum = glGenLists(1);

    glNewList (SphereListNum, GL_COMPILE);
    gluSphere(qobj, 1.0, NumFacets, NumFacets);
    gluDeleteQuadric(qobj);
    glEndList();
    SphereListCreated = true;
  }
}

void 
SphereObject::Set()
{
//   GLUquadricObj* qobj = gluNewQuadric();
//   gluQuadricDrawStyle(qobj, GLU_FILL);
  Start();
  // Set color
  float fcolor[4];
  fcolor[0] = Color[0]; fcolor[1] = Color[1]; 
  fcolor[2] = Color[2]; fcolor[3] = Color[3];
  glColor4f (Color[0], Color[1], Color[2], Color[3]);
  glMaterialfv (GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, fcolor);
  float spec[4] = { 1.0, 1.0, 1.0, 1.0 };
  glMaterialfv (GL_FRONT_AND_BACK, GL_SPECULAR, spec);
  glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 30.0);
  if (fabs(Color[3]-1.0) > 1.0e-3)
    glDepthMask(GL_FALSE);
  // Create sphere
  glPushMatrix();
  glTranslated (Pos[0], Pos[1], Pos[2]);
 
  if (OffScreen) {
    //    glCallList(OffScreenListNum);
    GLUquadricObj* qobj = gluNewQuadric();
    gluSphere(qobj, Radius, 20, 20);
    gluDeleteQuadric(qobj);
  }
  else {
    glScaled (Radius, Radius, Radius);
    glCallList(SphereListNum);
  }
  glPopMatrix();
  //  gluDeleteQuadric (qobj);
  if (fabs(Color[3]-1.0) > 1.0e-3)
  glDepthMask(GL_TRUE);
  End();
}

void
SphereObject::SetNumFacets(int num)
{
  NumFacets = num;
}

void 
SphereObject::SetPos (Vec3 pos)
{
  Pos = pos;
  Set();
}

void 
SphereObject::SetRadius (double radius)
{
  Radius = radius;
  Set();
}

void 
SphereObject::SetColor (Vec3 color)
{
  Color[0] = color[0];
  Color[1] = color[1];
  Color[2] = color[2];
  Color[3] = 1.0;
  Set();
}

void 
SphereObject::SetColor (Vec4 color)
{
  Color = color;
  Set();
}

void
SphereObject::SetBox (Vec3 box)
{
  Mat3 lattice;
  lattice = box[0], 0.0, 0.0, 0.0, box[1], 0.0, 0.0, 0.0, box[2];
  SetBox (lattice);
}

void
SphereObject::SetBox (Mat3 lattice)
{
  Lattice = lattice;
}


void
SphereObject::DrawPOV (FILE *fout, string rotString)
{
  fprintf (fout, "intersection {\n");
  fprintf (fout, "  box {\n");
//   fprintf (fout, "    <%10.8f, %10.8f, %10.8f>,\n",
// 	   -0.5*Box[0], -0.5*Box[1], -0.5*Box[2]);
//   fprintf (fout, "    <%10.8f, %10.8f, %10.8f>\n",
// 	   0.5*Box[0],   0.5*Box[1],  0.5*Box[2]);
  fprintf (fout, "    <-0.5, -0.5, -0.5>,\n");
  fprintf (fout, "    < 0.5,  0.5,  0.5>\n");
  fprintf (fout, "    matrix < %10.8f, %10.8f, %10.8f,\n",
	   Lattice(0,0), Lattice(0,1), Lattice(0,2));
  fprintf (fout, "             %10.8f, %10.8f, %10.8f,\n",
	   Lattice(1,0), Lattice(1,1), Lattice(1,2));
  fprintf (fout, "             %10.8f, %10.8f, %10.8f,\n",
	   Lattice(2,0), Lattice(2,1), Lattice(2,2));
  fprintf (fout, "             %10.8f, %10.8f, %10.8f>\n",
	   0.0, 0.0, 0.0);
  fprintf (fout, "%s", rotString.c_str());
  fprintf (fout, "  }\n");
  fprintf (fout, "  sphere {\n");
  fprintf (fout, "    <%10.8f, %10.8f, %10.8f>, %10.8f\n",
	   Pos[0], Pos[1], Pos[2], Radius);
  fprintf (fout, "%s", rotString.c_str());
  fprintf (fout, "  }\n");    // Sphere
  if (Color[3] < 0.99999)
    fprintf (fout, "    pigment { color rgbt <%1.5f %1.5f %1.5f %1.5f> }\n", 
	     0.6*Color[0], 0.6*Color[1], 0.6*Color[2], 1.0-0.7*Color[3]);
  else
    fprintf (fout, "    pigment { color rgb <%1.5f %1.5f %1.5f> }\n", 
	     Color[0], Color[1], Color[2]);
  fprintf (fout, "    finish {\n");
  fprintf (fout, "      ambient 0.1\n  diffuse 0.8\n");
  //  fprintf (fout, "    reflection 0.5\n");
  fprintf (fout, "      specular 0.6\n");
  fprintf (fout, "      roughness 0.025 \n");
  fprintf (fout, "    }\n"); //
  fprintf (fout, "}\n");    // Intersection

}
