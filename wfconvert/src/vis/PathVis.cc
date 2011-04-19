
///////////////////////////////////////////////////////////////////////////////
//
// Simple OpenGL scene.
//
///////////////////////////////////////////////////////////////////////////////

#include "PathVis.h"
#include "GLObject.h"
#include "PathObject.h"
#include "BoxObject.h"

namespace Trackball {
  extern "C" {
    #include "trackball.h"
  }
}

PathVisClass::PathVisClass() : 
  View(*this)
{
  NumLists = 0;
  //
  // Configure OpenGL-capable visual.
  //

  Glib::RefPtr<Gdk::GL::Config> glconfig;

  // Try double-buffered visual
  glconfig = Gdk::GL::Config::create(Gdk::GL::MODE_RGB    |
                                     Gdk::GL::MODE_DEPTH  |
                                     Gdk::GL::MODE_DOUBLE);
  if (!glconfig)  {
    std::cerr << "*** Cannot find the double-buffered visual.\n"
	      << "*** Trying single-buffered visual.\n";
    
    // Try single-buffered visual
    glconfig = Gdk::GL::Config::create(Gdk::GL::MODE_RGB   |
				       Gdk::GL::MODE_DEPTH);
    if (!glconfig) {
      std::cerr << "*** Cannot find any OpenGL-capable visual.\n";
      std::exit(1);
    }
  }

  // Set OpenGL-capability to the widget.
    
  set_gl_capability(glconfig);
}


PathVisClass::~PathVisClass()
{
}


void PathVisClass::Invalidate()
{
  get_window()->invalidate_rect(get_allocation(), false);
}



void PathVisClass::on_realize()
{
  // We need to call the base on_realize()
  Gtk::DrawingArea::on_realize();

  glShadeModel(GL_SMOOTH);
  glEnable (GL_LIGHTING);
  glEnable (GL_LINE_SMOOTH);
  glEnable (GL_POLYGON_SMOOTH);
  glDisable (GL_POLYGON_SMOOTH);
  glEnable (GL_MULTISAMPLE);
  glEnable (GL_COLOR_MATERIAL);

  // Get GL::Window.
  GLwindow = get_gl_window();

  // *** OpenGL BEGIN ***
  if (!GLwindow->gl_begin(get_gl_context()))
    return;

  static GLfloat light_diffuse[] = {1.0, 1.0, 1.0, 1.0};
  static GLfloat light_ambient[] = {0.2, 0.2, 0.2, 1.0};
  static GLfloat light_specular[]= {1.0, 1.0, 1.0, 1.0};
  static GLfloat light_position[] = {1.0, 1.0, 1.0, 0.0};
  glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse);
  glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient);
  glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular);
  glLightfv(GL_LIGHT0, GL_POSITION, light_position);
  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);
  glEnable(GL_DEPTH_TEST);
  glLightModeli(GL_LIGHT_MODEL_COLOR_CONTROL, GL_SEPARATE_SPECULAR_COLOR);

  glClearColor(1.0, 1.0, 1.0, 1.0);
  glClearDepth(1.0);

  glViewport(0, 0, get_width(), get_height());

  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(40.0, 1.0, 1.0, 10.0);

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  gluLookAt(0.0, 0.0, 3.0,
            0.0, 0.0, 0.0,
            0.0, 1.0, 0.0);
  
  glTranslatef(0.0, 0.0, -3.0);

  GLwindow->gl_end();
  // *** OpenGL END ***
}

bool PathVisClass::on_configure_event(GdkEventConfigure* event)
{
  // *** OpenGL BEGIN ***

  GLwindow = get_gl_window();

  if (!GLwindow->gl_begin(get_gl_context()))
    return false;

  glViewport(0, 0, get_width(), get_height());

  GLwindow->gl_end();
  // *** OpenGL END ***

  return true;
}

void PathVisClass::GLRender()
{  
  glClearColor(1.0, 1.0, 1.0, 1.0);
  glClearDepth(1.0);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  View.GLtransform_axes();
  vector<GLObject*>::iterator iter = Objects.begin();
  while (iter != Objects.end()) {
    if ((*iter)->IsAxes)
      (*iter)->Draw();
    iter++;
  }

  View.GLtransform();
  glEnable (GL_LINE_SMOOTH);
  glEnable (GL_BLEND);
  glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glHint (GL_LINE_SMOOTH_HINT, GL_NICEST);
  iter = Objects.begin();
  while (iter != Objects.end()) {
    if (!(*iter)->IsAxes)
      (*iter)->Draw();
    iter++;
  }
}


void PathVisClass::POVRender(string filename)
{
  FILE *fout;
  fout = fopen (filename.c_str(), "w");
  if (fout == NULL) {
    cerr << "Cannot open " << filename << " for writing.\n";
    return;
  }
  
  // do global settings
  fprintf (fout, "background { color rgbt <1.0, 1.0, 1.0, 1.0> }\n");
  fprintf (fout, "global_settings { ambient_light rgb <4.0, 4.0, 4.0> }\n");
  double light_position[] = {50.0, 50.0, 60.0, 0.0};
  fprintf (fout, "light_source { <%8.5f, %8.5f, %8.5f> \n",
	   light_position[0], light_position[1], light_position[2]);  
  fprintf (fout, "    color rgb <1.0, 1.0, 1.0> \n");
  fprintf (fout, "    area_light <-8,8,-8>, <8,-8,-8>, 5, 5\n");
  fprintf (fout, "    jitter\n");
  fprintf (fout, "    adaptive 5\n");
  fprintf (fout, "}\n\n");

  fprintf (fout, "light_source { <%8.5f, %8.5f, %8.5f> \n",
	   light_position[0], light_position[1], light_position[2]);  
  fprintf (fout, "    color rgb <0.5, 0.5, 0.5> \n");
  fprintf (fout, "    shadowless\n");
  fprintf (fout, "}\n\n");

  double light_position2[] = {-150.0, 100.0, -50.0, 0.0};
  fprintf (fout, "light_source { <%8.5f, %8.5f, %8.5f> \n",
	   light_position2[0], light_position2[1], light_position2[2]);  
  fprintf (fout, "    color rgb <0.75, 0.75, 0.75> \n");
  fprintf (fout, "    area_light <-18,18,-18>, <18,-18,-18>, 5, 5\n");
  fprintf (fout, "    jitter\n");
  fprintf (fout, "    adaptive 5\n");
  fprintf (fout, "}\n\n");

//   double light_position2[] = {-80.0, 30.0, 50.0, 0.0};
//   fprintf (fout, "light_source { <%8.5f, %8.5f, %8.5f> \n",
// 	   light_position2[0], light_position2[1], light_position2[2]);
//   fprintf (fout, "    color rgb <1.333, 1.333, 1.333> \n}\n\n");
  View.POVtransform (fout);

  string rotString = View.RotationString();

  vector<GLObject*>::iterator iter = Objects.begin();
  while (iter != Objects.end()) {
    (*iter)->DrawPOV(fout, rotString);
    iter++;
  }


  fclose (fout);
}


bool PathVisClass::on_expose_event(GdkEventExpose* event)
{
  // *** OpenGL BEGIN ***
  if (!GLwindow->gl_begin(get_gl_context()))
    return false;

  GLRender();

  // Swap buffers.
  if (GLwindow->is_double_buffered())
    GLwindow->swap_buffers();
  else
    glFlush();

  GLwindow->gl_end();
  // *** OpenGL END ***

  return true;
}


// void PathVisClass::AddPath (Array<Vec3,1> &path, bool closed)
// {
//   NumLists++;
//   glNewList(NumLists, GL_COMPILE);
//   glLineWidth (5.0);
//   glColor3f (1.0, 0.1, 0.1);
//   glBegin(GL_LINE_STRIP);
//     for (int i=0; i<path.size(); i++)
//       //      glVertex3dv(path(i)[0], path(i)[1], path(i)[2]);
//       glVertex3dv((double *)(&path(i)));
//     glEnd();
//   glEndList();
// }


///////////////////////////////////////////////////////////////////////////////
//
// The application class.
//
///////////////////////////////////////////////////////////////////////////////

class Simple : public Gtk::Window
{

protected:
  // signal handlers:
  void on_button_quit_clicked();

protected:
  // member widgets:
  Gtk::VBox m_VBox;
  Gtk::Button m_ButtonQuit;

public:
  PathVisClass PathVis;
  
  Simple();
  virtual ~Simple();

};

Simple::Simple()
  : m_VBox(false, 0), m_ButtonQuit("Quit")
{
  //
  // Top-level window.
  //

  set_title("Simple");

  // Get automatically redrawn if any of their children changed allocation.
  set_reallocate_redraws(true);

  add(m_VBox);

  //
  // Simple OpenGL scene.
  //

  PathVis.set_size_request(400, 400);

  m_VBox.pack_start(PathVis);

  //
  // Simple quit button.
  //

  m_ButtonQuit.signal_clicked().connect(
    sigc::mem_fun(*this, &Simple::on_button_quit_clicked));

  m_VBox.pack_start(m_ButtonQuit, Gtk::PACK_SHRINK, 0);

  //
  // Show window.
  //

  show_all();
}

Simple::~Simple()
{}

void Simple::on_button_quit_clicked()
{
  Gtk::Main::quit();
}


///////////////////////////////////////////////////////////////////////////////
//
// Main.
//
///////////////////////////////////////////////////////////////////////////////

// int main(int argc, char** argv)
// {
//   Gtk::Main kit(argc, argv);

//   //
//   // Init gtkglextmm.
//   //

//   Gtk::GL::init(argc, argv);

//   //
//   // Query OpenGL extension version.
//   //

//   int major, minor;
//   Gdk::GL::query_version(major, minor);
//   std::cout << "OpenGL extension version - "
//             << major << "." << minor << std::endl;

//   //
//   // Instantiate and run the application.
//   //

//   Simple simple;

//   Array<Vec3, 1> path(5);
//   path(0) = Vec3(-0.5,  0.5, -0.5);
//   path(1) = Vec3( 0.5,  0.5, -0.5);
//   path(2) = Vec3( 0.5, -0.5, -0.5);
//   path(3) = Vec3(-0.5, -0.5, -0.5);
//   path(4) = Vec3(-0.5,  0.5, -0.5);

//   PathObject *p1 = new PathObject();
//   p1->SetColor (0.0, 0.0, 1.0);
//   p1->Set (path);

//   simple.PathVis.Objects.push_back(p1);

//   // simple.PathVis.AddPath (path);
//   for (int i=0; i<5; i++)
//     path(i) += Vec3(0.0, 0.0, 1.0);
//   //  simple.PathVis.AddPath (path);
  
//   PathObject *p2 = new PathObject();
//   p2->SetColor (1.0, 0.0, 0.0);
//   p2->Set (path);
//   simple.PathVis.Objects.push_back(p2);
  
//   BoxObject *box = new BoxObject;
//   box->Set (2.0, 1.0, 0.5);
//   simple.PathVis.Objects.push_back(box);


//   kit.run(simple);

//   return 0;
// }
