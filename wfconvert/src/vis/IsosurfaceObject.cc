#include "IsosurfaceObject.h"
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/gle.h>


extern "C" void
GTS_IsoCartesianWrapper(gdouble **a, GtsCartesianGrid g,
			guint i, gpointer data)
{
  ((Isosurface *)data)->IsoCartesian(a, g, i);
}


void 
Isosurface::Construct()
{
  
  Surface = gts_surface_new (gts_surface_class(),
			     gts_face_class(),
			     gts_edge_class(),
			     gts_vertex_class());
  GtsCartesianGrid g;
  g.nx = InterpFactor*Spline.Nx;
  g.ny = InterpFactor*Spline.Ny;
  g.nz = InterpFactor*Spline.Nz;
  g.x  = Spline.Xgrid->Start;
  g.y  = Spline.Ygrid->Start;
  g.z  = Spline.Zgrid->Start;
  g.dx = (Spline.Xgrid->End - Spline.Xgrid->Start)/(g.nx-1);
  g.dy = (Spline.Xgrid->End - Spline.Ygrid->Start)/(g.nx-1);
  g.dz = (Spline.Xgrid->End - Spline.Zgrid->Start)/(g.nx-1);

  GtsIsoCartesianFunc f = GTS_IsoCartesianWrapper;
  gts_isosurface_cartesian (Surface, g, f, this, IsoVal);

}


void
Isosurface::IsoCartesian (gdouble **a, GtsCartesianGrid g, guint i)
{
  double z = g.z+g.dz*i;

  for (int ix=0; ix<g.nx; ix++) {
    double x = g.x+g.dx*ix; 
    for (int iy=0; iy<g.ny; iy++) {
      double y = g.y+g.dy*iy; 
      a[ix][iy] = Spline(x,y,z);
    }
  }
}

void
Isosurface::DrawGL()
{
  glBegin (GL_TRIANGLES);
  GSList *saveList = gts_surface_strip (Surface);
  GSList *stripList = saveList;
  while (stripList != NULL) {
    GSList *triList = (GSList*) stripList->data;
    while (triList != NULL) {
      GtsTriangle *triangle = (GtsTriangle*)triList->data;
      GtsPoint &p1 = triangle->e1->segment.v1->p;
      GtsPoint &p2 = triangle->e2->segment.v1->p;
      GtsPoint &p3 = triangle->e3->segment.v1->p;
      glVertex3d(p1.x, p1.y, p1.z);
      glVertex3d(p2.x, p2.y, p2.z);
      glVertex3d(p3.x, p3.y, p3.z);
//       v1[0] = triangle->e1->segment.v1.p->x;
//       v1[1] = triangle->e1->segment.v1.p->y;
//       v1[2] = triangle->e1->segment.v1.p->z;

      triList = triList->next;
    }
    stripList = stripList->next;
  }
  glEnd();
}
