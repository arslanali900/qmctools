#include "CylinderObject.h"
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>

int  CylinderObject::CylinderListNum(0);
bool CylinderObject::CylinderListCreated(false);
int  CylinderObject::OffScreenListNum(0);
bool CylinderObject::OffScreenListCreated(false);

CylinderObject::CylinderObject(bool offScreen) : 
  Radius(1.0), Length(1.0), Color(1.0, 0.0, 0.0, 1.0), 
  Pos1(0.0, 0.0, 0.0),  Pos2(0.0, 0.0, 0.0),
  OffScreen(offScreen), NumFacets (20)
{
  if (offScreen) {
    if (!OffScreenListCreated) {
      glEnable(GL_NORMALIZE);
//       GLUquadricObj* qobj = gluNewQuadric();
//       gluQuadricDrawStyle(qobj, GLU_FILL);

//       OffScreenListNum = glGenLists(1);
//       glNewList (OffScreenListNum, GL_COMPILE);
//       gluCylinder(qobj, 1.0, 30, 30);
//       gluDeleteQuadric(qobj);
//       glEndList();
//       OffScreenListCreated = true;
    }
  }
  else if (!CylinderListCreated) {
    glEnable(GL_NORMALIZE);
    GLUquadricObj* qobj = gluNewQuadric();
    gluQuadricDrawStyle(qobj, GLU_FILL);
    CylinderListNum = glGenLists(1);

    glNewList (CylinderListNum, GL_COMPILE);
    gluCylinder(qobj, 1.0, 1.0, 1.0, NumFacets, 1);
    gluDeleteQuadric(qobj);
    glEndList();
    CylinderListCreated = true;
  }
}

void 
CylinderObject::Set()
{
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
  
  // Create cylinder
  glPushMatrix();
  Vec3 pos = 0.5*(Pos1 + Pos2);
  glTranslated (pos[0], pos[1], pos[2]);
 
  Vec3 dr =  Pos2 - Pos1;
  dr /= sqrt(dot(dr,dr));
  double theta = 180.0/M_PI * acos(dr[2]);
  double phi   = 180.0/M_PI * atan2 (dr[1], dr[0]);

  if (OffScreen) {
    //    glCallList(OffScreenListNum);
    GLUquadricObj* qobj = gluNewQuadric();
    gluCylinder(qobj, Radius, Radius, Length, 20, 1);
    gluDeleteQuadric(qobj);
  }
  else {
    glRotated (phi,   0.0, 0.0, 1.0);
    glRotated (theta, 0.0, 1.0, 0.0);
    glScaled (Radius, Radius, Length);
    
    glTranslated (0.0, 0.0, -0.5);
    glCallList(CylinderListNum);
  }
  glPopMatrix();
  if (fabs(Color[3]-1.0) > 1.0e-3)
  glDepthMask(GL_TRUE);
  End();
}

void
CylinderObject::SetNumFacets(int num)
{
  NumFacets = num;
}

void 
CylinderObject::SetPos (Vec3 pos1, Vec3 pos2)
{
  Pos1 = pos1;
  Pos2 = pos2;
  Length = sqrt(dot(pos1-pos2, pos1-pos2));
  Set();
}

void 
CylinderObject::SetRadius (double radius)
{
  Radius = radius;
  Set();
}

void 
CylinderObject::SetColor (Vec3 color)
{
  Color[0] = color[0];
  Color[1] = color[1];
  Color[2] = color[2];
  Color[3] = 1.0;
  Set();
}

void 
CylinderObject::SetColor (Vec4 color)
{
  Color = color;
  Set();
}

void
CylinderObject::SetBox (Vec3 box)
{
  Mat3 lattice;
  lattice = box[0], 0.0, 0.0, 0.0, box[1], 0.0, 0.0, 0.0, box[2];
  SetBox (lattice);
}

void
CylinderObject::SetBox (Mat3 lattice)
{
  Lattice = lattice;
}


void
CylinderObject::DrawPOV (FILE *fout, string rotString)
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
  fprintf (fout, "  cylinder {\n");
  fprintf (fout, "    <%10.8f, %10.8f, %10.8f>,\n",
	   Pos1[0], Pos1[1], Pos1[2]);
  fprintf (fout, "    <%10.8f, %10.8f, %10.8f>, %10.8f\n",
	   Pos2[0], Pos2[1], Pos2[2], Radius);
  fprintf (fout, "%s", rotString.c_str());
  fprintf (fout, "  }\n");    // Cylinder
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
