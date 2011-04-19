// -*- C++ -*-
/*
 * simple.cc:
 * Simple gtkglextmm example.
 *
 * written by Naofumi Yasufuku  <naofumi@users.sourceforge.net>
 */

#include <iostream>
#include <cstdlib>

#include <gtkmm.h>

#include <gtkglmm.h>

#ifdef G_OS_WIN32
#define WIN32_LEAN_AND_MEAN 1
#include <windows.h>
#endif

#include <GL/gl.h>
#include <GL/glu.h>


///////////////////////////////////////////////////////////////////////////////
//
// OpenGL frame buffer configuration utilities.
//
///////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////
//
// Simple OpenGL scene.
//
///////////////////////////////////////////////////////////////////////////////

class SimpleGLScene : public Gtk::DrawingArea,
                      public Gtk::GL::Widget<SimpleGLScene>
{
public:
  SimpleGLScene();
  virtual ~SimpleGLScene();

protected:
  virtual void on_realize();
  virtual bool on_configure_event(GdkEventConfigure* event);
  virtual bool on_expose_event(GdkEventExpose* event);
};

SimpleGLScene::SimpleGLScene()
{
  //
  // Configure OpenGL-capable visual.
  //

  Glib::RefPtr<Gdk::GL::Config> glconfig;

  // Try double-buffered visual
  glconfig = Gdk::GL::Config::create(Gdk::GL::MODE_RGB    |
                                     Gdk::GL::MODE_DEPTH  |
                                     Gdk::GL::MODE_DOUBLE);
  if (!glconfig)
    {
      std::cerr << "*** Cannot find the double-buffered visual.\n"
                << "*** Trying single-buffered visual.\n";

      // Try single-buffered visual
      glconfig = Gdk::GL::Config::create(Gdk::GL::MODE_RGB   |
                                         Gdk::GL::MODE_DEPTH);
      if (!glconfig)
        {
          std::cerr << "*** Cannot find any OpenGL-capable visual.\n";
          std::exit(1);
        }
    }

  //
  // Set OpenGL-capability to the widget.
  //

  set_gl_capability(glconfig);
}

SimpleGLScene::~SimpleGLScene()
{
}

void SimpleGLScene::on_realize()
{
  // We need to call the base on_realize()
  Gtk::DrawingArea::on_realize();

  glShadeModel(GL_SMOOTH);
  glEnable (GL_LIGHTING);

  //
  // Get GL::Window.
  //

  Glib::RefPtr<Gdk::GL::Window> glwindow = get_gl_window();

  //
  // GL calls.
  //

  // *** OpenGL BEGIN ***
  if (!glwindow->gl_begin(get_gl_context()))
    return;

  GLUquadricObj* qobj = gluNewQuadric();
  gluQuadricDrawStyle(qobj, GLU_FILL);
  glNewList(1, GL_COMPILE);

  glColor3f (0.0, 0.0, 1.0);


  GLfloat mat_specular[] = { 1.0, 1.0, 1.0, 1.0 };
  GLfloat mat_diffuse[]  = { 0.1, 0.1, 1.0, 1.0 };
  GLfloat mat_shininess[] = { 50.0 };
  glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
  glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
  glMaterialfv(GL_FRONT, GL_SHININESS, mat_shininess);

  glMaterialfv(GL_BACK, GL_SPECULAR, mat_specular);
  glMaterialfv(GL_BACK, GL_DIFFUSE, mat_diffuse);
  glMaterialfv(GL_BACK, GL_SHININESS, mat_shininess);
  for (double theta=0.0; theta<=360.0; theta+=6.0) {
    glPushMatrix();
    glRotated(theta, 0.0, 1.0, 0.0);
    glTranslated(0.0,0.0, 1.5);
    glBegin(GL_QUAD_STRIP);
    const int N = 100;
    double delta = 2.0*M_PI/(double)N;
    for (double phi=0; phi <= 2.0*M_PI; phi+= delta) {
      glVertex3f(-0.1, 0.2*cos(phi), 0.2*sin(phi));
      glVertex3f( 0.1, 0.2*cos(phi), 0.2*sin(phi));
    }
    glEnd();
    glPopMatrix();
  }


  gluSphere(qobj, 1.0, 20, 20);
  // draw quad strip




  glEndList();

  static GLfloat light_diffuse[] = {1.0, 1.0, 1.0, 1.0};
  static GLfloat light_ambient[] = {1.0, 1.0, 1.0, 1.0};
  static GLfloat light_position[] = {1.0, 1.0, 1.0, 0.0};
  glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse);
  glLightfv(GL_LIGHT0, GL_POSITION, light_position);
  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);
  glEnable(GL_DEPTH_TEST);

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
  glRotated (-45.0, 1.0, 0.0, 0.0);


  glwindow->gl_end();
  // *** OpenGL END ***
}

bool SimpleGLScene::on_configure_event(GdkEventConfigure* event)
{
  //
  // Get GL::Window.
  //

  Glib::RefPtr<Gdk::GL::Window> glwindow = get_gl_window();

  //
  // GL calls.
  //

  // *** OpenGL BEGIN ***
  if (!glwindow->gl_begin(get_gl_context()))
    return false;

  glViewport(0, 0, get_width(), get_height());

  glwindow->gl_end();
  // *** OpenGL END ***

  return true;
}

bool SimpleGLScene::on_expose_event(GdkEventExpose* event)
{
  //
  // Get GL::Window.
  //

  Glib::RefPtr<Gdk::GL::Window> glwindow = get_gl_window();

  //
  // GL calls.
  //

  // *** OpenGL BEGIN ***
  if (!glwindow->gl_begin(get_gl_context()))
    return false;

  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  glCallList(1);

  // Swap buffers.
  if (glwindow->is_double_buffered())
    glwindow->swap_buffers();
  else
    glFlush();

  glwindow->gl_end();
  // *** OpenGL END ***

  return true;
}


///////////////////////////////////////////////////////////////////////////////
//
// The application class.
//
///////////////////////////////////////////////////////////////////////////////

class Simple : public Gtk::Window
{
public:
  Simple();
  virtual ~Simple();

protected:
  // signal handlers:
  void on_button_quit_clicked();

protected:
  // member widgets:
  Gtk::VBox m_VBox;
  SimpleGLScene m_SimpleGLScene;
  Gtk::Button m_ButtonQuit;
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

  m_SimpleGLScene.set_size_request(200, 200);

  m_VBox.pack_start(m_SimpleGLScene);

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

int main(int argc, char** argv)
{
  Gtk::Main kit(argc, argv);

  //
  // Init gtkglextmm.
  //

  Gtk::GL::init(argc, argv);

  //
  // Query OpenGL extension version.
  //

  int major, minor;
  Gdk::GL::query_version(major, minor);
  std::cout << "OpenGL extension version - "
            << major << "." << minor << std::endl;

  //
  // Instantiate and run the application.
  //

  Simple simple;

  kit.run(simple);

  return 0;
}
