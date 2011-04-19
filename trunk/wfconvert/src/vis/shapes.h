// -*- C++ -*-
/*
 * shapes.h:
 * shapes demo.
 *
 * written by Naofumi Yasufuku  <naofumi@users.sourceforge.net>
 */

#ifndef _SHAPES_H
#define _SHAPES_H

#include <gtkmm.h>

#include <gtkglmm.h>


///////////////////////////////////////////////////////////////////////////////
//
// Shapes classes.
//
///////////////////////////////////////////////////////////////////////////////

namespace Shapes
{

  class Scene;

  //
  // View class.
  //

  class View : public sigc::trackable
  {
    friend class Scene;

  public:
    static const float NEAR_CLIP;
    static const float FAR_CLIP;

    static const float INIT_POS_X;
    static const float INIT_POS_Y;
    static const float INIT_POS_Z;

    static const float INIT_AXIS_X;
    static const float INIT_AXIS_Y;
    static const float INIT_AXIS_Z;
    static const float INIT_ANGLE;

    static const float INIT_SCALE;

    static const float SCALE_MAX;
    static const float SCALE_MIN;

    static const float ANIMATE_THRESHOLD;

  public:
    View();
    virtual ~View();

  public:
    void frustum(int w, int h);

    void xform();

    void reset();

    void set_pos(float x, float y, float z)
    { m_Pos[0] = x; m_Pos[1] = y; m_Pos[2] = z; }

    void set_quat(float q0, float q1, float q2, float q3)
    { m_Quat[0] = q0; m_Quat[1] = q1; m_Quat[2] = q2; m_Quat[3] = q3; }

    void set_scale(float scale)
    { m_Scale = scale; }

    void enable_animation();

    void disable_animation();

    bool is_animate() const
    { return m_Animate; }

  protected:
    // Signal handlers:
    virtual bool on_button_press_event(GdkEventButton* event, Scene* scene);
    virtual bool on_button_release_event(GdkEventButton* event, Scene* scene);
    virtual bool on_motion_notify_event(GdkEventMotion* event, Scene* scene);

  private:
    float m_Pos[3];
    float m_Quat[4];
    float m_Scale;

    float m_QuatDiff[4];
    float m_BeginX;
    float m_BeginY;
    float m_DX;
    float m_DY;

    bool m_Animate;

  };


  //
  // Model class.
  //

  class Model
  {
    friend class Scene;

  public:
    static const unsigned int NUM_SHAPES;

    enum ShapeType
      {
        CUBE,
        SPHERE,
        CONE,
        TORUS,
        TETRAHEDRON,
        OCTAHEDRON,
        DODECAHEDRON,
        ICOSAHEDRON,
        TEAPOT,
      };

    static const ShapeType SHAPE_CUBE;
    static const ShapeType SHAPE_SPHERE;
    static const ShapeType SHAPE_CONE;
    static const ShapeType SHAPE_TORUS;
    static const ShapeType SHAPE_TETRAHEDRON;
    static const ShapeType SHAPE_OCTAHEDRON;
    static const ShapeType SHAPE_DODECAHEDRON;
    static const ShapeType SHAPE_ICOSAHEDRON;
    static const ShapeType SHAPE_TEAPOT;

  public:
    
    struct MaterialProp
    {
      GLfloat ambient[4];
      GLfloat diffuse[4];
      GLfloat specular[4];
      GLfloat shininess;
    };

    static const MaterialProp MAT_EMERALD;
    static const MaterialProp MAT_JADE;
    static const MaterialProp MAT_OBSIDIAN;
    static const MaterialProp MAT_PEARL;
    static const MaterialProp MAT_RUBY;
    static const MaterialProp MAT_TURQUOISE;
    static const MaterialProp MAT_BRASS;
    static const MaterialProp MAT_BRONZE;
    static const MaterialProp MAT_CHROME;
    static const MaterialProp MAT_COPPER;
    static const MaterialProp MAT_GOLD;
    static const MaterialProp MAT_SILVER;

  public:
    Model();
    virtual ~Model();

  private:
    void init_gl(Glib::RefPtr<Gdk::GL::Drawable>& gldrawable);

  public:
    void draw(Glib::RefPtr<Gdk::GL::Drawable>& gldrawable);

    void set_shape(ShapeType shape)
    { m_CurrentShape = shape; }

    void set_material(const MaterialProp* material)
    { m_CurrentMat = material; }

  private:
    unsigned int m_ListBase;
    ShapeType m_CurrentShape;
    const MaterialProp* m_CurrentMat;

  };


  //
  // Scene class.
  //

  class Scene : public Gtk::GL::DrawingArea
  {
    friend class View;
    friend class Model;

  public:
    // OpenGL scene related constants:
    static const float CLEAR_COLOR[4];
    static const float CLEAR_DEPTH;

    static const float LIGHT0_POSITION[4];
    static const float LIGHT0_AMBIENT[4];
    static const float LIGHT0_DIFFUSE[4];

    static const float LIGHT_MODEL_AMBIENT[4];
    static const float LIGHT_MODEL_LOCAL_VIEWER[1];

  public:
    explicit Scene();
    virtual ~Scene();

  protected:
    // signal handlers:
    virtual void on_realize();
    virtual bool on_configure_event(GdkEventConfigure* event);
    virtual bool on_expose_event(GdkEventExpose* event);
    virtual bool on_button_press_event(GdkEventButton* event);
    virtual bool on_unmap_event(GdkEventAny* event);
    virtual bool on_visibility_notify_event(GdkEventVisibility* event);
    virtual bool on_idle();

  public:
    // Invalidate whole window.
    void invalidate() {
      get_window()->invalidate_rect(get_allocation(), false);
    }

    // Update window synchronously (fast).
    void update()
    { get_window()->process_updates(false); }

  protected:
    // idle signal connection:
    sigc::connection m_ConnectionIdle;

    void idle_add();
    void idle_remove();

  protected:
    void change_shape(Model::ShapeType shape);
    void change_material(const Model::MaterialProp* material);

  protected:
    Gtk::Menu* create_popup_menu();

  protected:
    // Popup menu:
    Gtk::Menu* m_Menu;

  protected:
    // OpenGL scene related objects:
    View m_View;
    Model m_Model;

  };


  //
  // Application class.
  //

  class Application : public Gtk::Window
  {
  public:
    static const Glib::ustring APP_NAME;

  public:
    Application();
    virtual ~Application();

  protected:
    // signal handlers:
    virtual void on_button_quit_clicked();
    virtual bool on_key_press_event(GdkEventKey* event);

  protected:
    // member widgets:
    Gtk::VBox m_VBox;
    Scene m_Scene;
    Gtk::Button m_ButtonQuit;
  };


} // namespace Shapes


#endif // _SHAPES_H
