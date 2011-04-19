#include "BoxObject.h"
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/gle.h>
#include <GL/glut.h>

void 
BoxObject::Set(Vec3 box, bool useClip)
{
  Set (box[0], box[1], box[2], useClip);
}

void 
BoxObject::Set (double lx, double ly, double lz, bool useClip)
{
  Mat3 lattice;
  lattice = lx, 0.0, 0.0, 0.0, ly, 0.0, 0.0, 0.0, lz;
  Set (lattice, useClip);
}

void
BoxObject::Set (Mat3 lattice, bool useClip)
{
  Set (lattice, true, useClip);
}

void
BoxObject::Set (Mat3 lattice, bool visible, bool useClip)
{
  Visible = visible;
  Lx = lattice(0,0); Ly=lattice(1,1); Lz=lattice(2,2);
  Mat3 &l = lattice;
  LatticeVecs[0] = Vec3 (lattice(0,0), lattice(0,1), lattice(0,2));
  LatticeVecs[1] = Vec3 (lattice(1,0), lattice(1,1), lattice(1,2));
  LatticeVecs[2] = Vec3 (lattice(2,0), lattice(2,1), lattice(2,2));
  Start();
  glLineWidth (2.0);
  float fcolor[4];
  fcolor[0] = Color[0]; fcolor[1] = Color[1]; fcolor[2] = Color[2];
  fcolor[3] = 1.0;
  glMaterialfv (GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, fcolor);
  glColor3d (Color[0], Color[1], Color[2]);

  // Create the clippling planes.  These planes are constructed in
  // pairs by computing the normals to the planes by cross products.
  Vec3 p01 = cross (LatticeVecs[0], LatticeVecs[1]);
  Vec3 p12 = cross (LatticeVecs[1], LatticeVecs[2]);
  Vec3 p20 = cross (LatticeVecs[2], LatticeVecs[0]);
  p01 = 1.0/sqrt(dot(p01, p01)) * p01;
  p12 = 1.0/sqrt(dot(p12, p12)) * p12;
  p20 = 1.0/sqrt(dot(p20, p20)) * p20;
  double d01 = fabs(0.5001*dot(LatticeVecs[2], p01));
  double d12 = fabs(0.5001*dot(LatticeVecs[0], p12));
  double d20 = fabs(0.5001*dot(LatticeVecs[1], p20));

  if (useClip) {
    GLdouble eqn0[4] = { p01[0], p01[1], p01[2], d01};
    GLdouble eqn1[4] = {-p01[0],-p01[1],-p01[2], d01};
    GLdouble eqn2[4] = { p12[0], p12[1], p12[2], d12};
    GLdouble eqn3[4] = {-p12[0],-p12[1],-p12[2], d12};
    GLdouble eqn4[4] = { p20[0], p20[1], p20[2], d20};
    GLdouble eqn5[4] = {-p20[0],-p20[1],-p20[2], d20};
    glClipPlane(GL_CLIP_PLANE0, eqn0);
    glClipPlane(GL_CLIP_PLANE1, eqn1);
    glClipPlane(GL_CLIP_PLANE2, eqn2);
    glClipPlane(GL_CLIP_PLANE3, eqn3);
    glClipPlane(GL_CLIP_PLANE4, eqn4);
    glClipPlane(GL_CLIP_PLANE5, eqn5);
    glEnable(GL_CLIP_PLANE0);
    glEnable(GL_CLIP_PLANE1);
    glEnable(GL_CLIP_PLANE2);
    glEnable(GL_CLIP_PLANE3);
    glEnable(GL_CLIP_PLANE4);
    glEnable(GL_CLIP_PLANE5);
  }
  if (visible) {
    // Now draw the cell
    Vec3 a[3], ma[3];
    a[0] = Vec3 (lattice(0,0), lattice(0,1), lattice(0,2));
    a[1] = Vec3 (lattice(1,0), lattice(1,1), lattice(1,2));
    a[2] = Vec3 (lattice(2,0), lattice(2,1), lattice(2,2));
    ma[0] = -1.0*a[0];  ma[1] = -1.0*a[1];  ma[2] = -1.0*a[2];
    Vec3 r[8];
    r[0] = 0.5*(ma[0] + ma[1] + ma[2]);
    r[1] = 0.5*(ma[0] + ma[1] +  a[2]);
    r[2] = 0.5*(ma[0] +  a[1] +  a[2]);
    r[3] = 0.5*(ma[0] +  a[1] + ma[2]);
    r[4] = 0.5*( a[0] + ma[1] + ma[2]);
    r[5] = 0.5*( a[0] + ma[1] +  a[2]);
    r[6] = 0.5*( a[0] +  a[1] +  a[2]);
    r[7] = 0.5*( a[0] +  a[1] + ma[2]);
    glBegin(GL_LINES);
    glVertex3dv(&(r[0][0])); glVertex3dv (&(r[1][0]));
    glVertex3dv(&(r[1][0])); glVertex3dv (&(r[2][0]));
    glVertex3dv(&(r[2][0])); glVertex3dv (&(r[3][0]));
    glVertex3dv(&(r[3][0])); glVertex3dv (&(r[0][0]));
    glVertex3dv(&(r[4][0])); glVertex3dv (&(r[5][0]));
    glVertex3dv(&(r[5][0])); glVertex3dv (&(r[6][0]));
    glVertex3dv(&(r[6][0])); glVertex3dv (&(r[7][0]));
    glVertex3dv(&(r[7][0])); glVertex3dv (&(r[4][0]));
    glVertex3dv(&(r[0][0])); glVertex3dv (&(r[4][0]));
    glVertex3dv(&(r[1][0])); glVertex3dv (&(r[5][0]));
    glVertex3dv(&(r[2][0])); glVertex3dv (&(r[6][0]));
    glVertex3dv(&(r[3][0])); glVertex3dv (&(r[7][0]));
    glEnd();
  }
  End();
}

// void 
// BoxObject::Set (double lx, double ly, double lz, bool useClip)
// {

//   Lx=lx; Ly=ly; Lz=lz;
//   Start();
//   if (useClip) {
//     GLdouble eqn0[4] = {  1.0,  0.0,  0.0,  0.5001*lx };
//     GLdouble eqn1[4] = { -1.0,  0.0,  0.0,  0.5001*lx };
//     GLdouble eqn2[4] = {  0.0,  1.0,  0.0,  0.5001*ly };
//     GLdouble eqn3[4] = {  0.0, -1.0,  0.0,  0.5001*ly };
//     GLdouble eqn4[4] = {  0.0,  0.0,  1.0,  0.5001*lz };
//     GLdouble eqn5[4] = {  0.0,  0.0, -1.0,  0.5001*lz };
//     glClipPlane(GL_CLIP_PLANE0, eqn0);
//     glClipPlane(GL_CLIP_PLANE1, eqn1);
//     glClipPlane(GL_CLIP_PLANE2, eqn2);
//     glClipPlane(GL_CLIP_PLANE3, eqn3);
//     glClipPlane(GL_CLIP_PLANE4, eqn4);
//     glClipPlane(GL_CLIP_PLANE5, eqn5);
//     glEnable(GL_CLIP_PLANE0);
//     glEnable(GL_CLIP_PLANE1);
//     glEnable(GL_CLIP_PLANE2);
//     glEnable(GL_CLIP_PLANE3);
//     glEnable(GL_CLIP_PLANE4);
//     glEnable(GL_CLIP_PLANE5);
//   }
//   else {
//     glDisable(GL_CLIP_PLANE0);
//     glDisable(GL_CLIP_PLANE1);
//     glDisable(GL_CLIP_PLANE2);
//     glDisable(GL_CLIP_PLANE3);
//     glDisable(GL_CLIP_PLANE4);
//     glDisable(GL_CLIP_PLANE5);
//   }
//   glPushMatrix();
//   glLineWidth (2.0);
//   float fcolor[4];
//   fcolor[0] = Color[0]; fcolor[1] = Color[1]; fcolor[2] = Color[2];
//   fcolor[3] = 1.0;
//   glMaterialfv (GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, fcolor);
//   glColor3d (Color[0], Color[1], Color[2]);
//   glScaled (lx, ly, lz);
//   glutWireCube (1.0);
//   glTranslated (-0.5, -0.5, -0.5);
//   glPopMatrix();
//   End();
// }

void BoxObject::SetColor (double red, double green, double blue)
{
  Color[0] = red;
  Color[1] = green;
  Color[2] = blue;
}


void 
BoxObject::POVLine (FILE *fout, 
		    double x1, double y1, double z1,
		    double x2, double y2, double z2,
		    double radius, string rotString)
{
  fprintf (fout, "cylinder {\n");
  fprintf (fout, "  <%12.8f, %12.8f, %12.8f>,\n",
	   x1, y1, z1);
  fprintf (fout, "  <%12.8f, %12.8f, %12.8f>,\n",
	   x2, y2, z2);
  fprintf (fout, "  %10.8f\n", radius);
  fprintf (fout, "  pigment { color rgb <%1.5f %1.5f %1.5f> }\n", 
	   0.0, 0.0, 0.0);
  fprintf (fout, "%s", rotString.c_str());
  fprintf (fout, "}\n\n");
}

void
BoxObject::DrawPOV (FILE *fout, string rotString)
{
  if (!Visible)
    return;

  double l[3];
  l[0] = sqrt(dot(LatticeVecs[0], LatticeVecs[0]));
  l[1] = sqrt(dot(LatticeVecs[1], LatticeVecs[1]));
  l[2] = sqrt(dot(LatticeVecs[2], LatticeVecs[2]));
  double maxDim = max(max(l[0],l[1]),l[2]);
  double minDim = min(min(l[0],l[1]),l[2]);
  double radius = (minDim + maxDim)/900.0;
  Vec3 a[3], ma[3];
  a[0] = LatticeVecs[0];
  a[1] = LatticeVecs[1];
  a[2] = LatticeVecs[2];
  ma[0] = -1.0*a[0];  ma[1] = -1.0*a[1];  ma[2] = -1.0*a[2];
  Vec3 r[8];
  r[0] = 0.5*(ma[0] + ma[1] + ma[2]);
  r[1] = 0.5*(ma[0] + ma[1] +  a[2]);
  r[2] = 0.5*(ma[0] +  a[1] +  a[2]);
  r[3] = 0.5*(ma[0] +  a[1] + ma[2]);
  r[4] = 0.5*( a[0] + ma[1] + ma[2]);
  r[5] = 0.5*( a[0] + ma[1] +  a[2]);
  r[6] = 0.5*( a[0] +  a[1] +  a[2]);
  r[7] = 0.5*( a[0] +  a[1] + ma[2]);
  
  POVLine (fout, r[0][0], r[0][1], r[0][2],
 	         r[1][0], r[1][1], r[1][2], radius, rotString);
  POVLine (fout, r[1][0], r[1][1], r[1][2],
 	         r[2][0], r[2][1], r[2][2], radius, rotString);
  POVLine (fout, r[2][0], r[2][1], r[2][2],
 	         r[3][0], r[3][1], r[3][2], radius, rotString);
  POVLine (fout, r[3][0], r[3][1], r[3][2],
 	         r[0][0], r[0][1], r[0][2], radius, rotString);

  POVLine (fout, r[4][0], r[4][1], r[4][2],
 	         r[5][0], r[5][1], r[5][2], radius, rotString);
  POVLine (fout, r[5][0], r[5][1], r[5][2],
 	         r[6][0], r[6][1], r[6][2], radius, rotString);
  POVLine (fout, r[6][0], r[6][1], r[6][2],
 	         r[7][0], r[7][1], r[7][2], radius, rotString);
  POVLine (fout, r[7][0], r[7][1], r[7][2],
 	         r[4][0], r[4][1], r[4][2], radius, rotString);

  POVLine (fout, r[0][0], r[0][1], r[0][2],
 	         r[4][0], r[4][1], r[4][2], radius, rotString);
  POVLine (fout, r[1][0], r[1][1], r[1][2],
 	         r[5][0], r[5][1], r[5][2], radius, rotString);
  POVLine (fout, r[2][0], r[2][1], r[2][2],
 	         r[6][0], r[6][1], r[6][2], radius, rotString);
  POVLine (fout, r[3][0], r[3][1], r[3][2],
 	         r[7][0], r[7][1], r[7][2], radius, rotString);

//   POVLine (fout, -0.5*Lx, -0.5*Ly, -0.5*Lz,
// 	          0.5*Lx ,-0.5*Ly, -0.5*Lz, radius, rotString);
//   POVLine (fout, -0.5*Lx,  0.5*Ly, -0.5*Lz,
// 	          0.5*Lx , 0.5*Ly, -0.5*Lz, radius, rotString);
//   POVLine (fout, -0.5*Lx, -0.5*Ly, -0.5*Lz,
// 	         -0.5*Lx , 0.5*Ly, -0.5*Lz, radius, rotString);
//   POVLine (fout,  0.5*Lx, -0.5*Ly, -0.5*Lz,
// 	          0.5*Lx , 0.5*Ly, -0.5*Lz, radius, rotString);

//   POVLine (fout, -0.5*Lx, -0.5*Ly,  0.5*Lz,
// 	          0.5*Lx ,-0.5*Ly,  0.5*Lz, radius, rotString);
//   POVLine (fout, -0.5*Lx,  0.5*Ly,  0.5*Lz,
// 	          0.5*Lx , 0.5*Ly,  0.5*Lz, radius, rotString);
//   POVLine (fout, -0.5*Lx, -0.5*Ly,  0.5*Lz,
// 	         -0.5*Lx , 0.5*Ly,  0.5*Lz, radius, rotString);
//   POVLine (fout,  0.5*Lx, -0.5*Ly,  0.5*Lz,
// 	          0.5*Lx , 0.5*Ly,  0.5*Lz, radius, rotString);

//   POVLine (fout, -0.5*Lx, -0.5*Ly, -0.5*Lz,
// 	         -0.5*Lx ,-0.5*Ly,  0.5*Lz, radius, rotString);
//   POVLine (fout, -0.5*Lx,  0.5*Ly, -0.5*Lz,
// 	         -0.5*Lx , 0.5*Ly,  0.5*Lz, radius, rotString);
//   POVLine (fout,  0.5*Lx, -0.5*Ly, -0.5*Lz,
// 	          0.5*Lx ,-0.5*Ly,  0.5*Lz, radius, rotString);
//   POVLine (fout,  0.5*Lx,  0.5*Ly, -0.5*Lz,
// 	          0.5*Lx , 0.5*Ly,  0.5*Lz, radius, rotString);

}
