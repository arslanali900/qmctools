#include <Common/IO/IO.h>
#include <iostream>
#include <cstdlib>

#include <gtkmm.h>

#include <gtkglmm.h>
#include <GL/gl.h>
#include <GL/glu.h>
#include <Common/Blitz.h>
#include "GLObject.h"
#include "BoxObject.h"
#include "SphereObject.h"
#include "DiskObject.h"
#include "Isosurface.h"


class PathVisClass;

class ViewClass : public sigc::trackable
{
private:
  friend class PathVisClass;
  double StartX, StartY;
  bool Button1Pressed, Button2Pressed, Button3Pressed; 
  double MinScale, MaxScale;
  PathVisClass &PathVis;
  double Distance;
  bool UsePerspective;
  double xTrans, yTrans;
public:
  double Scale, OldScale;
  double Quaternion[4];
  double RotMat[4][4];
  Vec3 LightPos;

  bool OnButtonPress   (GdkEventButton* event);
  bool OnButtonRelease (GdkEventButton* event);
  bool OnMotion        (GdkEventMotion* event);
  bool OnScroll        (GdkEventScroll* event);
  void SetDistance (double dist);

  inline void SetPerspective (bool usePersp) 
  { UsePerspective = usePersp; }
  void GLtransform();
  void GLtransform_axes();
  void POVtransform (FILE *fout);
  string RotationString();
  void Write (IOSectionClass &out);
  void Read  (IOSectionClass &in);
  void Reset();

  ViewClass (PathVisClass &pathVis);
};

class PathVisClass : public Gtk::DrawingArea,
		     public Gtk::GL::Widget<PathVisClass>
{
  friend class ViewClass;
protected:
  Glib::RefPtr<Gdk::GL::Window> GLwindow;

  virtual void on_realize();
  virtual bool on_configure_event(GdkEventConfigure* event);
  virtual bool on_expose_event(GdkEventExpose* event);
  int NumLists;
public:
  ViewClass View;
  vector<GLObject *> Objects;
  void GLRender();
  void POVRender(string filename);
  void Invalidate();

  // Constructor
  PathVisClass();
  // Destructor
  virtual ~PathVisClass();
};

