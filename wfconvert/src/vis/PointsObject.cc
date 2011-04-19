#include "PointsObject.h"
#include <GL/gl.h>

PointsObject::PointsObject(bool offscreen) :
  Color(0.0, 0.0, 1.0, 1.0)
{
}
  
void
PointsObject::DrawPOV (FILE *fout, string rotString)
{

}

void
PointsObject::SetColor (Vec3 color)
{
  Color = Vec4(color[0], color[1], color[2], 1.0);
}

void
PointsObject::SetColor (Vec4 color)
{
  Color = Vec4(color[0], color[1], color[2], color[3]);
}


void
PointsObject::SetPoints (Array<Vec3,2> &points)
{
  Start();
  glEnable(GL_POINT_SMOOTH);
  glPointSize(1.0);
  // Set color
  float fcolor[4];
  fcolor[0] = Color[0]; fcolor[1] = Color[1]; 
  fcolor[2] = Color[2]; fcolor[3] = Color[3];
  glColor4f (Color[0], Color[1], Color[2], Color[3]);
  glMaterialfv (GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, fcolor);
  // float spec[4] = { 1.0, 1.0, 1.0, 1.0 };
  // glMaterialfv (GL_FRONT_AND_BACK, GL_SPECULAR, spec);
  // glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 30.0);
  if (fabs(Color[3]-1.0) > 1.0e-3)
    glDepthMask(GL_FALSE);

  glBegin(GL_POINTS);
  for (int i=0; i<points.extent(0); i++)
    for (int j=0; j<points.extent(1); j++) {
      float pos[3];
      glVertex3dv (&(points(i,j)[0]));
    }
  glEnd();
  if (fabs(Color[3]-1.0) > 1.0e-3)
    glDepthMask(GL_TRUE);
  End();
}
