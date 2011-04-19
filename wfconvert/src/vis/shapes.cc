// -*- C++ -*-
/*
 * shapes.cc:
 * shapes demo.
 *
 * written by Naofumi Yasufuku  <naofumi@users.sourceforge.net>
 */

#include <iostream>
#include <cstdlib>
#include <cstring>
#include <cmath>

#ifdef G_OS_WIN32
#define WIN32_LEAN_AND_MEAN 1
#include <windows.h>
#endif

#include <GL/gl.h>
#include <GL/glu.h>

#include "shapes.h"

//
// Trackball utilities.
//
namespace Trackball {
  extern "C" {
    #include "trackball.h"
  }
}

#define DIG_2_RAD (G_PI / 180.0)
#define RAD_2_DIG (180.0 / G_PI)


///////////////////////////////////////////////////////////////////////////////
//
// OpenGL frame buffer configuration utilities.
//
///////////////////////////////////////////////////////////////////////////////

struct GLConfigUtil
{
  static void print_gl_attrib(const Glib::RefPtr<const Gdk::GL::Config>& glconfig,
                              const char* attrib_str,
                              int attrib,
                              bool is_boolean);

  static void examine_gl_attrib(const Glib::RefPtr<const Gdk::GL::Config>& glconfig);
};

//
// Print a configuration attribute.
//
void GLConfigUtil::print_gl_attrib(const Glib::RefPtr<const Gdk::GL::Config>& glconfig,
                                   const char* attrib_str,
                                   int attrib,
                                   bool is_boolean)
{
  int value;

  if (glconfig->get_attrib(attrib, value))
    {
      std::cout << attrib_str << " = ";
      if (is_boolean)
        std::cout << (value == true ? "true" : "false") << std::endl;
      else
        std::cout << value << std::endl;
    }
  else
    {
      std::cout << "*** Cannot get "
                << attrib_str
                << " attribute value\n";
    }
}

//
// Print configuration attributes.
//
void GLConfigUtil::examine_gl_attrib(const Glib::RefPtr<const Gdk::GL::Config>& glconfig)
{
  std::cout << "\nOpenGL visual configurations :\n\n";

  std::cout << "glconfig->is_rgba() = "
            << (glconfig->is_rgba() ? "true" : "false")
            << std::endl;
  std::cout << "glconfig->is_double_buffered() = "
            << (glconfig->is_double_buffered() ? "true" : "false")
            << std::endl;
  std::cout << "glconfig->is_stereo() = "
            << (glconfig->is_stereo() ? "true" : "false")
            << std::endl;
  std::cout << "glconfig->has_alpha() = "
            << (glconfig->has_alpha() ? "true" : "false")
            << std::endl;
  std::cout << "glconfig->has_depth_buffer() = "
            << (glconfig->has_depth_buffer() ? "true" : "false")
            << std::endl;
  std::cout << "glconfig->has_stencil_buffer() = "
            << (glconfig->has_stencil_buffer() ? "true" : "false")
            << std::endl;
  std::cout << "glconfig->has_accum_buffer() = "
            << (glconfig->has_accum_buffer() ? "true" : "false")
            << std::endl;

  std::cout << std::endl;

  print_gl_attrib(glconfig, "Gdk::GL::USE_GL",           Gdk::GL::USE_GL,           true);
  print_gl_attrib(glconfig, "Gdk::GL::BUFFER_SIZE",      Gdk::GL::BUFFER_SIZE,      false);
  print_gl_attrib(glconfig, "Gdk::GL::LEVEL",            Gdk::GL::LEVEL,            false);
  print_gl_attrib(glconfig, "Gdk::GL::RGBA",             Gdk::GL::RGBA,             true);
  print_gl_attrib(glconfig, "Gdk::GL::DOUBLEBUFFER",     Gdk::GL::DOUBLEBUFFER,     true);
  print_gl_attrib(glconfig, "Gdk::GL::STEREO",           Gdk::GL::STEREO,           true);
  print_gl_attrib(glconfig, "Gdk::GL::AUX_BUFFERS",      Gdk::GL::AUX_BUFFERS,      false);
  print_gl_attrib(glconfig, "Gdk::GL::RED_SIZE",         Gdk::GL::RED_SIZE,         false);
  print_gl_attrib(glconfig, "Gdk::GL::GREEN_SIZE",       Gdk::GL::GREEN_SIZE,       false);
  print_gl_attrib(glconfig, "Gdk::GL::BLUE_SIZE",        Gdk::GL::BLUE_SIZE,        false);
  print_gl_attrib(glconfig, "Gdk::GL::ALPHA_SIZE",       Gdk::GL::ALPHA_SIZE,       false);
  print_gl_attrib(glconfig, "Gdk::GL::DEPTH_SIZE",       Gdk::GL::DEPTH_SIZE,       false);
  print_gl_attrib(glconfig, "Gdk::GL::STENCIL_SIZE",     Gdk::GL::STENCIL_SIZE,     false);
  print_gl_attrib(glconfig, "Gdk::GL::ACCUM_RED_SIZE",   Gdk::GL::ACCUM_RED_SIZE,   false);
  print_gl_attrib(glconfig, "Gdk::GL::ACCUM_GREEN_SIZE", Gdk::GL::ACCUM_GREEN_SIZE, false);
  print_gl_attrib(glconfig, "Gdk::GL::ACCUM_BLUE_SIZE",  Gdk::GL::ACCUM_BLUE_SIZE,  false);
  print_gl_attrib(glconfig, "Gdk::GL::ACCUM_ALPHA_SIZE", Gdk::GL::ACCUM_ALPHA_SIZE, false);

  std::cout << std::endl;
}


///////////////////////////////////////////////////////////////////////////////
//
// Shapes classes.
//
///////////////////////////////////////////////////////////////////////////////

namespace Shapes
{

  //
  // View class implementation.
  //

  const float View::NEAR_CLIP   = 5.0;
  const float View::FAR_CLIP    = 60.0;

  const float View::INIT_POS_X  = 0.0;
  const float View::INIT_POS_Y  = 0.0;
  const float View::INIT_POS_Z  = -10.0;

  const float View::INIT_AXIS_X = 1.0;
  const float View::INIT_AXIS_Y = 0.0;
  const float View::INIT_AXIS_Z = 0.0;
  const float View::INIT_ANGLE  = 0.0;

  const float View::INIT_SCALE  = 1.0;

  const float View::SCALE_MAX   = 2.0;
  const float View::SCALE_MIN   = 0.5;

  const float View::ANIMATE_THRESHOLD = 25.0;

  View::View()
    : m_Scale(INIT_SCALE),
      m_BeginX(0.0), m_BeginY(0.0),
      m_DX(0.0), m_DY(0.0),
      m_Animate(false)
  {
    reset();
  }

  View::~View()
  {
  }

  void View::frustum(int w, int h)
  {
    glViewport(0, 0, w, h);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();

    if (w > h) {
      float aspect = static_cast<float>(w) / static_cast<float>(h);
      glFrustum(-aspect, aspect, -1.0, 1.0, NEAR_CLIP, FAR_CLIP);
    } else {
      float aspect = static_cast<float>(h) / static_cast<float>(w);
      glFrustum(-1.0, 1.0, -aspect, aspect, NEAR_CLIP, FAR_CLIP);
    }

    glMatrixMode(GL_MODELVIEW);
  }

  void View::xform()
  {
    glTranslatef(m_Pos[0], m_Pos[1], m_Pos[2]);

    glScalef(m_Scale, m_Scale, m_Scale);

    float m[4][4];
    Trackball::add_quats(m_QuatDiff, m_Quat, m_Quat);
    Trackball::build_rotmatrix(m, m_Quat);
    glMultMatrixf(&m[0][0]);
  }

  void View::reset()
  {
    m_Pos[0] = INIT_POS_X;
    m_Pos[1] = INIT_POS_Y;
    m_Pos[2] = INIT_POS_Z;

    float sine = sin(0.5 * INIT_ANGLE * DIG_2_RAD);
    m_Quat[0] = INIT_AXIS_X * sine;
    m_Quat[1] = INIT_AXIS_Y * sine;
    m_Quat[2] = INIT_AXIS_Z * sine;
    m_Quat[3] = cos(0.5 * INIT_ANGLE * DIG_2_RAD);

    m_Scale = INIT_SCALE;

    m_QuatDiff[0] = 0.0;
    m_QuatDiff[1] = 0.0;
    m_QuatDiff[2] = 0.0;
    m_QuatDiff[3] = 1.0;
  }

  void View::enable_animation()
  {
    m_Animate = true;
  }

  void View::disable_animation()
  {
    m_Animate = false;

    m_QuatDiff[0] = 0.0;
    m_QuatDiff[1] = 0.0;
    m_QuatDiff[2] = 0.0;
    m_QuatDiff[3] = 1.0;
  }

  bool View::on_button_press_event(GdkEventButton* event,
                                   Scene* scene)
  {
    if (is_animate()) {
      if (event->button == 1) {
        disable_animation();
        scene->idle_remove();
        scene->invalidate();
      }
    } else {
      m_QuatDiff[0] = 0.0;
      m_QuatDiff[1] = 0.0;
      m_QuatDiff[2] = 0.0;
      m_QuatDiff[3] = 1.0;
    }

    m_BeginX = event->x;
    m_BeginY = event->y;

    // don't block
    return false;
  }

  bool View::on_button_release_event(GdkEventButton* event,
                                     Scene* scene)
  {
    if (!is_animate()) {
      if (event->button == 1 &&
          ((m_DX*m_DX + m_DY*m_DY) > ANIMATE_THRESHOLD)) {
        enable_animation();
        scene->idle_add();
      }
    }

    m_DX = 0.0;
    m_DY = 0.0;

    // don't block
    return false;
  }

  bool View::on_motion_notify_event(GdkEventMotion* event,
                                    Scene* scene)
  {
    float w = scene->get_width();
    float h = scene->get_height();
    float x = event->x;
    float y = event->y;
    bool redraw = false;

    // Rotation.
    if (event->state & GDK_BUTTON1_MASK) {
      Trackball::trackball(m_QuatDiff,
                           (2.0 * m_BeginX - w) / w,
                           (h - 2.0 * m_BeginY) / h,
                           (2.0 * x - w) / w,
                           (h - 2.0 * y) / h);

      m_DX = x - m_BeginX;
      m_DY = y - m_BeginY;

      redraw = true;
    }

    // Scaling.
    if (event->state & GDK_BUTTON2_MASK) {
      m_Scale = m_Scale * (1.0 + (y - m_BeginY) / h);
      if (m_Scale > SCALE_MAX)
        m_Scale = SCALE_MAX;
      else if (m_Scale < SCALE_MIN)
        m_Scale = SCALE_MIN;

      redraw = true;
    }

    m_BeginX = x;
    m_BeginY = y;

    if (redraw)
      scene->invalidate();

    // don't block
    return false;
  }


  //
  // Model class implementation.
  //

  const unsigned int Model::NUM_SHAPES = 9;

  const Model::ShapeType Model::SHAPE_CUBE         = CUBE;
  const Model::ShapeType Model::SHAPE_SPHERE       = SPHERE;
  const Model::ShapeType Model::SHAPE_CONE         = CONE;
  const Model::ShapeType Model::SHAPE_TORUS        = TORUS;
  const Model::ShapeType Model::SHAPE_TETRAHEDRON  = TETRAHEDRON;
  const Model::ShapeType Model::SHAPE_OCTAHEDRON   = OCTAHEDRON;
  const Model::ShapeType Model::SHAPE_DODECAHEDRON = DODECAHEDRON;
  const Model::ShapeType Model::SHAPE_ICOSAHEDRON  = ICOSAHEDRON;
  const Model::ShapeType Model::SHAPE_TEAPOT       = TEAPOT;

  const Model::MaterialProp Model::MAT_EMERALD = {
    {0.0215, 0.1745, 0.0215, 1.0},
    {0.07568, 0.61424, 0.07568, 1.0},
    {0.633, 0.727811, 0.633, 1.0},
    0.6
  };

  const Model::MaterialProp Model::MAT_JADE = {
    {0.135, 0.2225, 0.1575, 1.0},
    {0.54, 0.89, 0.63, 1.0},
    {0.316228, 0.316228, 0.316228, 1.0},
    0.1
  };

  const Model::MaterialProp Model::MAT_OBSIDIAN = {
    {0.05375, 0.05, 0.06625, 1.0},
    {0.18275, 0.17, 0.22525, 1.0},
    {0.332741, 0.328634, 0.346435, 1.0},
    0.3
  };

  const Model::MaterialProp Model::MAT_PEARL = {
    {0.25, 0.20725, 0.20725, 1.0},
    {1.0, 0.829, 0.829, 1.0},
    {0.296648, 0.296648, 0.296648, 1.0},
    0.088
  };

  const Model::MaterialProp Model::MAT_RUBY = {
    {0.1745, 0.01175, 0.01175, 1.0},
    {0.61424, 0.04136, 0.04136, 1.0},
    {0.727811, 0.626959, 0.626959, 1.0},
    0.6
  };

  const Model::MaterialProp Model::MAT_TURQUOISE = {
    {0.1, 0.18725, 0.1745, 1.0},
    {0.396, 0.74151, 0.69102, 1.0},
    {0.297254, 0.30829, 0.306678, 1.0},
    0.1
  };

  const Model::MaterialProp Model::MAT_BRASS = {
    {0.329412, 0.223529, 0.027451, 1.0},
    {0.780392, 0.568627, 0.113725, 1.0},
    {0.992157, 0.941176, 0.807843, 1.0},
    0.21794872
  };

  const Model::MaterialProp Model::MAT_BRONZE = {
    {0.2125, 0.1275, 0.054, 1.0},
    {0.714, 0.4284, 0.18144, 1.0},
    {0.393548, 0.271906, 0.166721, 1.0},
    0.2
  };

  const Model::MaterialProp Model::MAT_CHROME = {
    {0.25, 0.25, 0.25, 1.0},
    {0.4, 0.4, 0.4, 1.0},
    {0.774597, 0.774597, 0.774597, 1.0},
    0.6
  };

  const Model::MaterialProp Model::MAT_COPPER = {
    {0.19125, 0.0735, 0.0225, 1.0},
    {0.7038, 0.27048, 0.0828, 1.0},
    {0.256777, 0.137622, 0.086014, 1.0},
    0.1
  };

  const Model::MaterialProp Model::MAT_GOLD = {
    {0.24725, 0.1995, 0.0745, 1.0},
    {0.75164, 0.60648, 0.22648, 1.0},
    {0.628281, 0.555802, 0.366065, 1.0},
    0.4
  };

  const Model::MaterialProp Model::MAT_SILVER = {
    {0.19225, 0.19225, 0.19225, 1.0},
    {0.50754, 0.50754, 0.50754, 1.0},
    {0.508273, 0.508273, 0.508273, 1.0},
    0.4
  };

  Model::Model()
    : m_ListBase(0),
      m_CurrentShape(TEAPOT),
      m_CurrentMat(&MAT_SILVER)
  {
  }

  Model::~Model()
  {
  }

  void Model::init_gl(Glib::RefPtr<Gdk::GL::Drawable>& gldrawable)
  {
    /* Shape display lists */
    m_ListBase = glGenLists(NUM_SHAPES);

    /* Cube */
    glNewList(m_ListBase + CUBE, GL_COMPILE);
      gldrawable->draw_cube(true, 1.5);
    glEndList();

    /* Sphere */
    glNewList(m_ListBase + SPHERE, GL_COMPILE);
      gldrawable->draw_sphere(true, 1.0, 30, 30);
    glEndList();

    /* Cone */
    glNewList(m_ListBase + CONE, GL_COMPILE);
      glPushMatrix();
        glTranslatef(0.0, 0.0, -1.0);
        gldrawable->draw_cone(true, 1.0, 2.0, 30, 30);
      glPopMatrix();
    glEndList();

    /* Torus */
    glNewList(m_ListBase + TORUS, GL_COMPILE);
      gldrawable->draw_torus(true, 0.4, 0.8, 30, 30);
    glEndList();

    /* Tetrahedron */
    glNewList(m_ListBase + TETRAHEDRON, GL_COMPILE);
      glPushMatrix();
        glScalef(1.2, 1.2, 1.2);
        gldrawable->draw_tetrahedron(true);
      glPopMatrix();
    glEndList();

    /* Octahedron */
    glNewList(m_ListBase + OCTAHEDRON, GL_COMPILE);
      glPushMatrix();
        glScalef(1.2, 1.2, 1.2);
        gldrawable->draw_octahedron(true);
      glPopMatrix();
    glEndList();

    /* Dodecahedron */
    glNewList(m_ListBase + DODECAHEDRON, GL_COMPILE);
      glPushMatrix();
        glScalef(0.7, 0.7, 0.7);
        gldrawable->draw_dodecahedron(true);
      glPopMatrix();
    glEndList();

    /* Icosahedron */
    glNewList(m_ListBase + ICOSAHEDRON, GL_COMPILE);
      glPushMatrix();
        glScalef(1.2, 1.2, 1.2);
        gldrawable->draw_icosahedron(true);
      glPopMatrix();
    glEndList();

    /* Teapot */
    glNewList(m_ListBase + TEAPOT, GL_COMPILE);
      gldrawable->draw_teapot(true, 1.0);
    glEndList();
  }

  void Model::draw(Glib::RefPtr<Gdk::GL::Drawable>& gldrawable)
  {
    // Init GL context.
    static bool initialized = false;
    if (!initialized) {
      init_gl(gldrawable);
      initialized = true;
    }

    // Render shape
    glMaterialfv(GL_FRONT, GL_AMBIENT, m_CurrentMat->ambient);
    glMaterialfv(GL_FRONT, GL_DIFFUSE, m_CurrentMat->diffuse);
    glMaterialfv(GL_FRONT, GL_SPECULAR, m_CurrentMat->specular);
    glMaterialf(GL_FRONT, GL_SHININESS, m_CurrentMat->shininess * 128.0);
    glCallList(m_ListBase + m_CurrentShape);
  }


  //
  // Scene class implementation.
  //

  const float Scene::CLEAR_COLOR[4] = { 0.5, 0.5, 0.8, 1.0 };
  const float Scene::CLEAR_DEPTH    = 1.0;

  const float Scene::LIGHT0_POSITION[4] = { 0.0, 3.0, 3.0, 0.0 };
  const float Scene::LIGHT0_AMBIENT[4]  = { 0.0, 0.0, 0.0, 1.0 };
  const float Scene::LIGHT0_DIFFUSE[4]  = { 1.0, 1.0, 1.0, 1.0 };

  const float Scene::LIGHT_MODEL_AMBIENT[4]       = { 0.2, 0.2, 0.2, 1.0 };
  const float Scene::LIGHT_MODEL_LOCAL_VIEWER[1]  = { 0.0 };

  Scene::Scene()
    : m_Menu(0)
  {
    //
    // Configure OpenGL-capable visual.
    //

    Glib::RefPtr<Gdk::GL::Config> glconfig;

    // Try double-buffered visual
    glconfig = Gdk::GL::Config::create(Gdk::GL::MODE_RGB    |
                                       Gdk::GL::MODE_DEPTH  |
                                       Gdk::GL::MODE_DOUBLE);
    if (!glconfig) {
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

    // print frame buffer attributes.
    GLConfigUtil::examine_gl_attrib(glconfig);

    //
    // Set OpenGL-capability to the widget.
    //

    set_gl_capability(glconfig);

    //
    // Add events.
    //
    add_events(Gdk::BUTTON1_MOTION_MASK    |
               Gdk::BUTTON2_MOTION_MASK    |
               Gdk::BUTTON_PRESS_MASK      |
               Gdk::BUTTON_RELEASE_MASK    |
               Gdk::VISIBILITY_NOTIFY_MASK);

    // View transformation signals.
    signal_button_press_event().connect(
      sigc::bind(sigc::mem_fun(m_View, &View::on_button_press_event), this));
    signal_button_release_event().connect(
      sigc::bind(sigc::mem_fun(m_View, &View::on_button_release_event), this));
    signal_motion_notify_event().connect(
      sigc::bind(sigc::mem_fun(m_View, &View::on_motion_notify_event), this));

    //
    // Popup menu.
    //

    m_Menu = create_popup_menu();
  }

  Scene::~Scene()
  {
  }

  void Scene::on_realize()
  {
    // We need to call the base on_realize()
    Gtk::DrawingArea::on_realize();

    //
    // Get GL::Drawable.
    //

    Glib::RefPtr<Gdk::GL::Drawable> gldrawable = get_gl_drawable();

    //
    // GL calls.
    //

    // *** OpenGL BEGIN ***
    if (!gldrawable->gl_begin(get_gl_context()))
      return;

    glClearColor(CLEAR_COLOR[0], CLEAR_COLOR[1], CLEAR_COLOR[2], CLEAR_COLOR[3]);
    glClearDepth(CLEAR_DEPTH);

    glLightfv(GL_LIGHT0, GL_POSITION, LIGHT0_POSITION);
    glLightfv(GL_LIGHT0, GL_AMBIENT,  LIGHT0_AMBIENT);
    glLightfv(GL_LIGHT0, GL_DIFFUSE,  LIGHT0_DIFFUSE);
    glLightModelfv(GL_LIGHT_MODEL_AMBIENT, LIGHT_MODEL_AMBIENT);
    glLightModelfv(GL_LIGHT_MODEL_LOCAL_VIEWER, LIGHT_MODEL_LOCAL_VIEWER);

    glFrontFace(GL_CW);
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    glEnable(GL_AUTO_NORMAL);
    glEnable(GL_NORMALIZE);
    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LESS);

    gldrawable->gl_end();
    // *** OpenGL END ***
  }

  bool Scene::on_configure_event(GdkEventConfigure* event)
  {
    //
    // Get GL::Drawable.
    //

    Glib::RefPtr<Gdk::GL::Drawable> gldrawable = get_gl_drawable();

    //
    // GL calls.
    //

    // *** OpenGL BEGIN ***
    if (!gldrawable->gl_begin(get_gl_context()))
      return false;

    m_View.frustum(get_width(), get_height());

    gldrawable->gl_end();
    // *** OpenGL END ***

    return true;
  }

  bool Scene::on_expose_event(GdkEventExpose* event)
  {
    //
    // Get GL::Drawable.
    //

    Glib::RefPtr<Gdk::GL::Drawable> gldrawable = get_gl_drawable();

    //
    // GL calls.
    //

    // *** OpenGL BEGIN ***
    if (!gldrawable->gl_begin(get_gl_context()))
      return false;

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glLoadIdentity();

    // View transformation.
    m_View.xform();

    // Logo model.
    m_Model.draw(gldrawable);

    // Swap buffers.
    if (gldrawable->is_double_buffered())
      gldrawable->swap_buffers();
    else
      glFlush();

    gldrawable->gl_end();
    // *** OpenGL END ***

    return true;
  }

  bool Scene::on_button_press_event(GdkEventButton* event)
  {
    if (event->button == 3) {
      m_Menu->popup(event->button, event->time);
      return true;
    }

    // don't block
    return false;
  }

  bool Scene::on_unmap_event(GdkEventAny* event)
  {
    idle_remove();

    return true;
  }

  bool Scene::on_visibility_notify_event(GdkEventVisibility* event)
  {
    if (m_View.is_animate()) {
      if (event->state == GDK_VISIBILITY_FULLY_OBSCURED)
        idle_remove();
      else
        idle_add();
    }

    return true;
  }

  bool Scene::on_idle()
  {
    // Invalidate whole window.
    invalidate();
    // Update window synchronously (fast).
    update();

    return true;
  }

  void Scene::idle_add()
  {
    if (!m_ConnectionIdle.connected())
      m_ConnectionIdle = Glib::signal_idle().connect(
        sigc::mem_fun(*this, &Scene::on_idle), GDK_PRIORITY_REDRAW);
  }

  void Scene::idle_remove()
  {
    if (m_ConnectionIdle.connected())
      m_ConnectionIdle.disconnect();
  }

  void Scene::change_shape(Model::ShapeType shape)
  {
    m_Model.set_shape(shape);
    m_View.reset();
  }

  void Scene::change_material(const Model::MaterialProp* material)
  {
    m_Model.set_material(material);
  }

  Gtk::Menu* Scene::create_popup_menu()
  {

    // Shapes submenu
    Gtk::Menu* shapes_menu = Gtk::manage(new Gtk::Menu());
    {
      Gtk::Menu::MenuList& menu_list = shapes_menu->items();

      // Cube
      menu_list.push_back(Gtk::Menu_Helpers::MenuElem("Cube",
        sigc::bind(sigc::mem_fun(*this, &Scene::change_shape),
                   Model::SHAPE_CUBE)));

      // Sphere
      menu_list.push_back(Gtk::Menu_Helpers::MenuElem("Sphere",
        sigc::bind(sigc::mem_fun(*this, &Scene::change_shape),
                   Model::SHAPE_SPHERE)));

      // Cone
      menu_list.push_back(Gtk::Menu_Helpers::MenuElem("Cone",
        sigc::bind(sigc::mem_fun(*this, &Scene::change_shape),
                   Model::SHAPE_CONE)));

      // Torus
      menu_list.push_back(Gtk::Menu_Helpers::MenuElem("Torus",
        sigc::bind(sigc::mem_fun(*this, &Scene::change_shape),
                   Model::SHAPE_TORUS)));

      // Tetrahedron
      menu_list.push_back(Gtk::Menu_Helpers::MenuElem("Tetrahedron",
        sigc::bind(sigc::mem_fun(*this, &Scene::change_shape),
                   Model::SHAPE_TETRAHEDRON)));

      // Octahedron
      menu_list.push_back(Gtk::Menu_Helpers::MenuElem("Octahedron",
        sigc::bind(sigc::mem_fun(*this, &Scene::change_shape),
                   Model::SHAPE_OCTAHEDRON)));

      // Dodecahedron
      menu_list.push_back(Gtk::Menu_Helpers::MenuElem("Dodecahedron",
        sigc::bind(sigc::mem_fun(*this, &Scene::change_shape),
                   Model::SHAPE_DODECAHEDRON)));

      // Icosahedron
      menu_list.push_back(Gtk::Menu_Helpers::MenuElem("Icosahedron",
        sigc::bind(sigc::mem_fun(*this, &Scene::change_shape),
                   Model::SHAPE_ICOSAHEDRON)));

      // Teapot
      menu_list.push_back(Gtk::Menu_Helpers::MenuElem("Teapot",
        sigc::bind(sigc::mem_fun(*this, &Scene::change_shape),
                   Model::SHAPE_TEAPOT)));

    }

    // Materials submenu
    Gtk::Menu* materials_menu = Gtk::manage(new Gtk::Menu());
    {
      Gtk::Menu::MenuList& menu_list = materials_menu->items();

      // Emerald
      menu_list.push_back(Gtk::Menu_Helpers::MenuElem("Emerald",
        sigc::bind(sigc::mem_fun(*this, &Scene::change_material),
                   &Model::MAT_EMERALD)));

      // Jade
      menu_list.push_back(Gtk::Menu_Helpers::MenuElem("Jade",
        sigc::bind(sigc::mem_fun(*this, &Scene::change_material),
                   &Model::MAT_JADE)));

      // Obsidian
      menu_list.push_back(Gtk::Menu_Helpers::MenuElem("Obsidian",
        sigc::bind(sigc::mem_fun(*this, &Scene::change_material),
                   &Model::MAT_OBSIDIAN)));

      // Pearl
      menu_list.push_back(Gtk::Menu_Helpers::MenuElem("Pearl",
        sigc::bind(sigc::mem_fun(*this, &Scene::change_material),
                   &Model::MAT_PEARL)));

      // Ruby
      menu_list.push_back(Gtk::Menu_Helpers::MenuElem("Ruby",
        sigc::bind(sigc::mem_fun(*this, &Scene::change_material),
                   &Model::MAT_RUBY)));

      // Turquoise
      menu_list.push_back(Gtk::Menu_Helpers::MenuElem("Turquoise",
        sigc::bind(sigc::mem_fun(*this, &Scene::change_material),
                   &Model::MAT_TURQUOISE)));

      // Brass
      menu_list.push_back(Gtk::Menu_Helpers::MenuElem("Brass",
        sigc::bind(sigc::mem_fun(*this, &Scene::change_material),
                   &Model::MAT_BRASS)));

      // Bronze
      menu_list.push_back(Gtk::Menu_Helpers::MenuElem("Bronze",
        sigc::bind(sigc::mem_fun(*this, &Scene::change_material),
                   &Model::MAT_BRONZE)));

      // Chrome
      menu_list.push_back(Gtk::Menu_Helpers::MenuElem("Chrome",
        sigc::bind(sigc::mem_fun(*this, &Scene::change_material),
                   &Model::MAT_CHROME)));

      // Copper
      menu_list.push_back(Gtk::Menu_Helpers::MenuElem("Copper",
        sigc::bind(sigc::mem_fun(*this, &Scene::change_material),
                   &Model::MAT_COPPER)));

      // Gold
      menu_list.push_back(Gtk::Menu_Helpers::MenuElem("Gold",
        sigc::bind(sigc::mem_fun(*this, &Scene::change_material),
                   &Model::MAT_GOLD)));

      // Silver
      menu_list.push_back(Gtk::Menu_Helpers::MenuElem("Silver",
        sigc::bind(sigc::mem_fun(*this, &Scene::change_material),
                   &Model::MAT_SILVER)));

    }

    // Root popup menu
    Gtk::Menu* menu = Gtk::manage(new Gtk::Menu());
    {
      Gtk::Menu::MenuList& menu_list = menu->items();

      // Shapes submenu
      menu_list.push_back(Gtk::Menu_Helpers::MenuElem("Shapes",
                                                      *shapes_menu));

      // Materials submenu
      menu_list.push_back(Gtk::Menu_Helpers::MenuElem("Materials",
                                                      *materials_menu));

      // Quit
      menu_list.push_back(Gtk::Menu_Helpers::MenuElem("Quit",
        sigc::ptr_fun(&Gtk::Main::quit)));
    }

    return menu;
  }


  //
  // Application class implementation.
  //

  const Glib::ustring Application::APP_NAME = "Shapes";

  Application::Application()
    : m_VBox(false, 0), m_ButtonQuit("Quit")
  {
    //
    // Top-level window.
    //

    set_title(APP_NAME);

    // Get automatically redrawn if any of their children changed allocation.
    set_reallocate_redraws(true);

    add(m_VBox);

    //
    // Scene.
    //

    m_Scene.set_size_request(300, 300);

    m_VBox.pack_start(m_Scene);

    //
    // Simple quit button.
    //

    m_ButtonQuit.signal_clicked().connect(
      sigc::mem_fun(*this, &Application::on_button_quit_clicked));

    m_VBox.pack_start(m_ButtonQuit, Gtk::PACK_SHRINK, 0);

    //
    // Show window.
    //

    show_all();
  }

  Application::~Application()
  {
  }

  void Application::on_button_quit_clicked()
  {
    Gtk::Main::quit();
  }

  bool Application::on_key_press_event(GdkEventKey* event)
  {
    switch (event->keyval) {
    case GDK_Escape:
      Gtk::Main::quit();
      break;
    default:
      return true;
    }

    m_Scene.invalidate();

    return true;
  }


} // namespace Shapes


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

  Shapes::Application application;

  kit.run(application);

  return 0;
}
