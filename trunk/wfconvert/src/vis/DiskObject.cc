#include "DiskObject.h"
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>

int  DiskObject::DiskListNum(0);
bool DiskObject::DiskListCreated(false);
int  DiskObject::OffScreenListNum(0);
bool DiskObject::OffScreenListCreated(false);

DiskObject::DiskObject(bool offScreen) : 
  Radius(1.0), Color(1.0, 0.0, 0.0), Pos(0.0, 0.0, 0.0),
  OffScreen(offScreen), Axis(0)
{
  if (offScreen) {
    if (!OffScreenListCreated) {
      glEnable(GL_NORMALIZE);
      GLUquadricObj* qobj = gluNewQuadric();
      gluQuadricDrawStyle(qobj, GLU_FILL);
      OffScreenListNum = glGenLists(1);
      
      glNewList (OffScreenListNum, GL_COMPILE);
      gluDisk(qobj, 0.0, 1.0, 30, 1);
      gluDeleteQuadric(qobj);
      glEndList();
      OffScreenListCreated = true;
    }
  }
  else if (!DiskListCreated) {
    glEnable(GL_NORMALIZE);
    GLUquadricObj* qobj = gluNewQuadric();
    gluQuadricDrawStyle(qobj, GLU_FILL);
    DiskListNum = glGenLists(1);

    glNewList (DiskListNum, GL_COMPILE);
    gluDisk(qobj, 0.0, 1.0, 30, 1);
    gluDeleteQuadric(qobj);
    glEndList();
    DiskListCreated = true;
  }
}

// void 
// DiskObject::Set()
// {
//   Start();
//   // Set color
//   float fcolor[4];
//   fcolor[0] = Color[0]; fcolor[1] = Color[1]; fcolor[2] = Color[2];
//   fcolor[3] = 1.0;
//   glColor3f (Color[0], Color[1], Color[2]);
//   glMaterialfv (GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, fcolor);
//   float spec[4] = { 1.0, 1.0, 1.0, 1.0 };
//   glMaterialfv (GL_FRONT_AND_BACK, GL_SPECULAR, spec);
//   glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 30.0);

//   // Create disk
//   glPushMatrix();
//   glTranslated (Pos[0], Pos[1], Pos[2]);
//   glScaled (Radius, Radius, Radius);
//   if (Axis == 1)
//     glRotated (270.0, 0.0, 1.0, 0.0);
//   else if (Axis == 0)
//     glRotated (90.0,  0.0, 1.0, 0.0);
//   else if (Axis == 2)
//     glRotated (270.0, 1.0, 0.0, 0.0);
//   else if (Axis == 3)
//     glRotated (90.0,  1.0, 0.0, 0.0);
//   else if (Axis == 5)
//     glRotated (180.0, 1.0, 0.0, 0.0);
//   //  gluDisk(qobj, Radius, 20, 20);
//   if (OffScreen) {
//     // glCallList(OffScreenListNum);
//     GLUquadricObj* qobj = gluNewQuadric();
//     gluDisk(qobj, 0.0, 1.0, 30, 1);
//     gluDeleteQuadric(qobj);
//   }
//   else
//     glCallList(DiskListNum);
//   glPopMatrix();
//   //  gluDeleteQuadric (qobj);
//   End();
// }


void 
DiskObject::Set()
{
  Start();
  // Set color
  float fcolor[4];
  fcolor[0] = Color[0]; fcolor[1] = Color[1]; fcolor[2] = Color[2];
  fcolor[3] = 1.0;
  glColor3f (Color[0], Color[1], Color[2]);
  glMaterialfv (GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, fcolor);
  float spec[4] = { 1.0, 1.0, 1.0, 1.0 };
  glMaterialfv (GL_FRONT_AND_BACK, GL_SPECULAR, spec);
  glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 30.0);

  // Create disk
  glPushMatrix();
  glTranslated (Pos[0], Pos[1], Pos[2]);
  glScaled (Radius, Radius, Radius);
  Vec3 rotAxis = cross (NormVec, Vec3 (0.0, 0.0, 1.0));
  
  double angle = -180.0/M_PI * acos(dot(NormVec, Vec3(0.0, 0.0, 1.0)));
  glRotated (angle, rotAxis[0], rotAxis[1], rotAxis[2]);
  //  gluDisk(qobj, Radius, 20, 20);
  if (OffScreen) {
    // glCallList(OffScreenListNum);
    GLUquadricObj* qobj = gluNewQuadric();
    gluDisk(qobj, 0.0, 1.0, 30, 1);
    gluDeleteQuadric(qobj);
  }
  else
    glCallList(DiskListNum);
  glPopMatrix();
  //  gluDeleteQuadric (qobj);
  End();
}

void 
DiskObject::SetPos (Vec3 pos)
{
  Pos = pos;
  Set();
}

void
DiskObject::SetAxis(int axis)
{
  Axis = axis;
}

void 
DiskObject::SetRadius (double radius)
{
  Radius = radius;
}

void
DiskObject::SetNormVec (Vec3 norm)
{
  norm = norm/sqrt(dot(norm,norm));
  NormVec = norm;
}

void 
DiskObject::SetColor (Vec3 color)
{
  Color = color;
}


void
DiskObject::DrawPOV (FILE *fout, string rotString)
{
//   fprintf (fout, "disk {\n");
//   fprintf (fout, "  <%10.8f, %10.8f, %10.8f>, %10.8f\n",
// 	   Pos[0], Pos[1], Pos[2], Radius);
//   fprintf (fout, "%s", rotString.c_str());
//   fprintf (fout, "  pigment { color rgb <%1.5f %1.5f %1.5f> }\n", 
// 	   Color[0], Color[1], Color[2]);
//   fprintf (fout, "  finish {\n");
//   fprintf (fout, "    ambient 0.1\n  diffuse 0.8\n");
//   //  fprintf (fout, "    reflection 0.5\n");
//   fprintf (fout, "    specular 0.6\n");
//   fprintf (fout, "    roughness 0.025 \n");
//   fprintf (fout, "  }\n");
//   fprintf (fout, "}\n");

}
