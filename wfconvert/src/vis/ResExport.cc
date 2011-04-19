#include "ResExport.h"
#include "ResVis.h"
#include <GL/gl.h>
#include <GL/glu.h>
#include <Common/Splines/CubicSpline.h>
#include <sstream>

using namespace IO;

ResExportClass::ResExportClass(ResVisualClass &visual) :
    Visual(visual), Width(2000), Height(2000),
    BaseNameChooser("Export filename", Gtk::FILE_CHOOSER_ACTION_SAVE),
    WidthLabel("Width:"), HeightLabel("Height:"),
    POVTolAdjust(0.0, 0.0, 1.0, 0.05, 0.2)
{
  // Set up base name entry and chooser
  BaseNameFrame.set_label("Base filename");
  BaseNameFrame.add (BaseNameHBox);
  BaseNameHBox.pack_start(BaseNameEntry, Gtk::PACK_EXPAND_WIDGET, 5);
  BaseNameHBox.pack_start(BaseNameBrowseButton, Gtk::PACK_SHRINK);
  BaseNameBrowseButton.set_label("Browse");
  MainVBox.pack_start (BaseNameFrame, Gtk::PACK_SHRINK, 5);
  BaseNameChooser.add_button (Gtk::Stock::CANCEL, Gtk::RESPONSE_CANCEL);
  BaseNameChooser.add_button (Gtk::Stock::SAVE,   Gtk::RESPONSE_OK);
  
  SizeFrame.set_label ("Size");
  SizeFrame.add (SizeBox);
  WidthButton.set_range(1.0, 3000.0);
  HeightButton.set_range(1.0, 3000.0);
  WidthButton.set_increments (10.0, 100.0);
  HeightButton.set_increments (10.0, 100.0);
  WidthButton.set_value (700.0);
  HeightButton.set_value (700.0);
  WidthButton.set_digits(0);
  HeightButton.set_digits(0);
  WidthBox.pack_start (WidthLabel, Gtk::PACK_SHRINK,5);
  WidthBox.pack_start (WidthButton, Gtk::PACK_SHRINK);
  HeightBox.pack_start (HeightLabel, Gtk::PACK_SHRINK,5);
  HeightBox.pack_start (HeightButton, Gtk::PACK_SHRINK);
  WidthHeightBox.pack_start (WidthBox, Gtk::PACK_SHRINK, 5);
  WidthHeightBox.pack_start (HeightBox, Gtk::PACK_SHRINK, 5);
  SizeBox.pack_start (WidthHeightBox, Gtk::PACK_SHRINK);
  SizeBox.pack_start (RatioButton,    Gtk::PACK_SHRINK);
  RatioButton.set_label ("Fixed ratio");
  
  // Setup Rendering Type combo
  TypeFrame.set_label("Rendering Type");
  TypeFrame.add(TypeBox);
  TypeBox.pack_start(TypeCombo, Gtk::PACK_SHRINK);
  TypeCombo.insert_text(0, "Offscreen GL");
  TypeCombo.insert_text(1, "POVray");
  TypeCombo.set_active (0);
  ExportButton.set_label ("Export");
  CancelButton.set_label ("Cancel");

  // Setup POV widgets
  POVFrame.set_label ("POVray options");
  POVRenderButton.set_label("Render");
  POVAntiAliasButton.set_label("Anti-alias");
  POVOnScreenButton.set_label("On-screen");
  POVTransparentBGButton.set_label ("Transparent BG");
  POVTolLabel.set_text("AA tolerance");
  POVVBox1.pack_start(POVRenderButton, Gtk::PACK_SHRINK, 5);
  POVVBox1.pack_start(POVAntiAliasButton, Gtk::PACK_SHRINK, 5);
  POVVBox1.pack_start(POVOnScreenButton, Gtk::PACK_SHRINK, 5);
  POVVBox1.pack_start(POVTransparentBGButton, Gtk::PACK_SHRINK, 5);
  POVVBox2.pack_start(POVTolLabel,  Gtk::PACK_SHRINK, 1);
  POVVBox2.pack_start(POVTolerance, Gtk::PACK_SHRINK, 1);
  POVHBox.pack_start(POVVBox1, Gtk::PACK_SHRINK, 5);
  POVHBox.pack_start(POVVBox2, Gtk::PACK_SHRINK, 5);
  POVFrame.add (POVHBox);
  POVFrame.set_sensitive(false);
  POVTolerance.set_adjustment (POVTolAdjust);
  POVTolerance.set_digits(2);
  
  StillMovieFrame.set_label("Still/Movie");
  StillMovieFrame.add (StillMovieBox);
  StillMovieBox.pack_start (StillButton, Gtk::PACK_SHRINK, 5);
  StillMovieBox.pack_start (MovieButton, Gtk::PACK_SHRINK, 5);
  Gtk::RadioButtonGroup group = StillButton.get_group();
  MovieButton.set_group (group);
  StillButton.set_label("Current Frame (png)");
  MovieButton.set_label("Multi-frame movie (MPEG4)");
  
  MovieParamFrame.set_label("Movie parameters");
  MovieParamFrame.add (MovieParamBox);
  MovieParamBox.pack_start (FirstLastBox, Gtk::PACK_SHRINK, 5);
  MovieParamBox.pack_start (InterpBox, Gtk::PACK_SHRINK, 5);
  FirstFrameLabel.set_text ("First frame:");
  FirstFrameBox.pack_start (FirstFrameLabel,  Gtk::PACK_SHRINK, 5);
  FirstFrameBox.pack_start (FirstFrameButton, Gtk::PACK_SHRINK, 5);
  LastFrameLabel.set_text  ("Last frame:");
  LastFrameBox.pack_start  (LastFrameLabel,   Gtk::PACK_SHRINK, 5);
  LastFrameBox.pack_start  (LastFrameButton,  Gtk::PACK_SHRINK, 5);
  FirstLastBox.pack_start (FirstFrameBox, Gtk::PACK_EXPAND_PADDING);
  FirstLastBox.pack_start (LastFrameBox, Gtk::PACK_EXPAND_PADDING);
  InterpFactorLabel.set_text ("Interpolation factor:");
  InterpFactorBox.pack_start (InterpFactorLabel,  Gtk::PACK_SHRINK, 5);
  InterpFactorBox.pack_start (InterpFactorButton, Gtk::PACK_SHRINK, 5);
  InterpBox.pack_start (InterpFactorBox, Gtk::PACK_EXPAND_PADDING);
  
  MovieParamFrame.set_sensitive(false);
  
  ButtonBox.pack_start(ExportButton,    Gtk::PACK_SHRINK, 5);
  ButtonBox.pack_start(CancelButton,    Gtk::PACK_SHRINK, 5);
  MainVBox.pack_start (SizeFrame,       Gtk::PACK_SHRINK, 5);
  MainVBox.pack_start (TypeFrame,       Gtk::PACK_SHRINK, 5);
  MainVBox.pack_start (POVFrame,        Gtk::PACK_SHRINK, 5);
  MainVBox.pack_start (StillMovieFrame, Gtk::PACK_SHRINK, 5);
  MainVBox.pack_start (MovieParamFrame, Gtk::PACK_SHRINK, 5);
  MainVBox.pack_start (ButtonBox,       Gtk::PACK_SHRINK, 10);
  add (MainVBox);
  set_title ("Export");
  // Set signal handlers
  ExportButton.signal_clicked().connect
    (sigc::mem_fun(*this, &ResExportClass::OnExportButton));
  CancelButton.signal_clicked().connect
    (sigc::mem_fun(*this, &ResExportClass::OnCancelButton));
  BaseNameBrowseButton.signal_clicked().connect
    (sigc::mem_fun(*this, &ResExportClass::OnBrowseButton));
  RatioButton.signal_toggled().connect
    (sigc::mem_fun(*this, &ResExportClass::OnRatioToggle));
  WidthButton.signal_value_changed().connect
    (sigc::mem_fun(*this, &ResExportClass::OnWidthAdjust));
  HeightButton.signal_value_changed().connect
    (sigc::mem_fun(*this, &ResExportClass::OnHeightAdjust));
  StillButton.signal_toggled().connect
    (sigc::mem_fun(*this, &ResExportClass::OnStillMovie));
  BaseNameChooser.signal_selection_changed().connect
    (sigc::mem_fun(*this, &ResExportClass::OnChooserChange));
  BaseNameEntry.signal_activate().connect
    (sigc::mem_fun(*this, &ResExportClass::OnEntryChange));
  TypeCombo.signal_changed().connect
    (sigc::mem_fun(*this, &ResExportClass::OnTypeChange));
}


void ResExportClass::SetupWidgets()
{
  FirstFrameButton.set_digits(0);
  FirstFrameButton.set_range(1.0, Visual.GetNumSteps());
  FirstFrameButton.set_increments (1.0, 10.0);
  LastFrameButton.set_digits(0);
  LastFrameButton.set_range (1.0, Visual.GetNumSteps());
  LastFrameButton.set_increments (1.0, 10.0);
  InterpFactorButton.set_range(1.0, 500.0);
  InterpFactorButton.set_digits(0);
  InterpFactorButton.set_increments (1.0, 10.0);
}  

void ResExportClass::InitGLStuff()
{
//   glShadeModel(GL_SMOOTH);
//   glEnable (GL_LIGHTING);
//   glEnable (GL_LINE_SMOOTH);
//   glEnable (GL_POLYGON_SMOOTH);

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
  glEnable (GL_MULTISAMPLE);
  glEnable (GL_COLOR_MATERIAL);
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  

  glViewport(0, 0, Width, Height);
  
}


void ResExportClass::MakePixmap()
{
  GdkPixmap.clear();
  GLConfig.clear();

  GLConfig = Gdk::GL::Config::create (Gdk::GL::MODE_RGB   |
				      Gdk::GL::MODE_DEPTH |
				      Gdk::GL::MODE_SINGLE);
  if (!GLConfig) {
    cerr << "Cannot find a visual capable of OpenGL.\n";
    return;
  }

  GdkPixmap = Gdk::Pixmap::create (Visual.get_window(),
				   Width, Height, GLConfig->get_depth());
  
  Glib::RefPtr<Gdk::GL::Pixmap> GLPixmap = 
    Gdk::GL::ext(GdkPixmap).set_gl_capability (GLConfig);

  GLContext = Gdk::GL::Context::create (GLPixmap, false);

  GLPixmap->make_current(GLContext);
  GLPixmap->gl_begin(GLContext);

  // Render here
  InitGLStuff();
  Visual.DrawFrame (true);
  Visual.PathVis.GLRender();

  glFlush();
  GLPixmap->gl_end();
  GLPixmap->wait_gl();
}

void ResExportClass::Export (string filename)
{
//   GdkPixmap.clear();
//   GLConfig.clear();

//   GLConfig = Gdk::GL::Config::create (Gdk::GL::MODE_RGB   |
// 				      Gdk::GL::MODE_DEPTH |
// 				      Gdk::GL::MODE_SINGLE);
//   if (!GLConfig) {
//     cerr << "Cannot find a visual capable of OpenGL.\n";
//     return;
//   }

//   GdkPixmap = Gdk::Pixmap::create (Visual.get_window(),
// 				   Width, Height, GLConfig->get_depth());
  
//   Glib::RefPtr<Gdk::GL::Pixmap> GLPixmap = 
//     Gdk::GL::ext(GdkPixmap).set_gl_capability (GLConfig);

//   GLContext = Gdk::GL::Context::create (GLPixmap, false);

//   GLPixmap->make_current(GLContext);
//   GLPixmap->gl_begin(GLContext);

//   // Render here
//   InitGLStuff();
//   Visual.MakeFrame ((int)floor(Visual.FrameAdjust.get_value()));
//   Visual.PathVis.GLRender();

//   glFlush();
//   GLPixmap->gl_end();
//   GLPixmap->wait_gl();

  MakePixmap();
  Glib::RefPtr<Gdk::Drawable> drawable =  GdkPixmap;
  Glib::RefPtr<Gdk::Drawable> drawable2 =  GdkPixmap;

  cerr << "drawable = " << drawable << endl;

  Glib::RefPtr<Gdk::Pixbuf> Pbuf = Gdk::Pixbuf::create
    (drawable, GdkPixmap->get_colormap(), 0, 0, 0, 0, Width, Height);

  if (drawable == drawable2)
    cerr << "no change\n";
			
  Pbuf->save(filename, Glib::ustring("png"));
}



void
ResExportClass::ExportPOV(string basename)
{
  string filename = basename;
  filename.append(".pov");
  Visual.PathVis.POVRender (filename);
  if (POVRenderButton.get_active()) {
    POVFile = filename;
    Glib::Thread::create
      (sigc::mem_fun (*this, &ResExportClass::RenderPOV), true);
  }
}


void
ResExportClass::RenderPOV()
{
  stringstream povcmd;
  povcmd << "povray ";
  if (!POVOnScreenButton.get_active())
    povcmd << "-D ";
  if (POVAntiAliasButton.get_active()) {
    povcmd.setf(ios_base::floatfield, ios_base::fixed);
    povcmd.precision(3);
    povcmd << "+A" << (double)POVTolerance.get_value() << " ";
  }
  if (POVTransparentBGButton.get_active())
    povcmd << "+ua ";
  povcmd << "+W" << (int)WidthButton.get_value() 
	 << " +H" << (int)HeightButton.get_value() << " ";
  povcmd << POVFile;

  cerr << "povcmd = " << povcmd.str() << endl;

  system (povcmd.str().c_str());
}

void
ResExportClass::ExportMovie (string basename, 
			     int firstFrame, int lastFrame, 
			     int interpFactor)
{
  int i=0;
  int numFrames = lastFrame - firstFrame + 1;
  int numDigits = (int) ceil (log(numFrames+1)/log(10.0));
  for (int frame=firstFrame; frame <= lastFrame; frame++) {
    Visual.SetStep(frame);
    Visual.DrawFrame();
    string name = basename + "_";
    stringstream numstream;
    numstream << i;
    string numstr = numstream.str();
    int numzeros = numDigits - numstr.size();
    for (int j=0; j<numzeros; j++)
      name = name + "0";
    name = name + numstr;

    ExportPOV(name);

    i++;
  }
}



void
ResExportClass::OnExportButton()
{
  Width  = (int)round(WidthButton.get_value());
  Height = (int)round(HeightButton.get_value());
  //  string basename = BaseNameChooser.get_filename();
  string basename = BaseNameEntry.get_text();

  if (MovieButton.get_active()) {
    int firstFrame   = (int)round(FirstFrameButton.get_value())-1;
    int lastFrame    = (int)round(LastFrameButton.get_value())-1;
    int interpFactor = (int)round(InterpFactorButton.get_value());
    ExportMovie (basename, firstFrame, lastFrame, interpFactor);
  }
  else {
    if (TypeCombo.get_active_text() == "POVray") 
      ExportPOV(basename);
    else {
      basename.append(".png");
      Export(basename);
    }
  }
}

void
ResExportClass::OnCancelButton()
{
  hide();
}

void
ResExportClass::OnBrowseButton()
{
  int result = BaseNameChooser.run();
  if (result ==  Gtk::RESPONSE_OK) {
    string fname = BaseNameChooser.get_filename();
    BaseNameEntry.set_text(fname);
  }

  else if (result == Gtk::RESPONSE_CANCEL)
    cerr << "Response was OK.\n";
  BaseNameChooser.hide();
}

void 
ResExportClass::OnHeightAdjust()
{
  if (RatioButton.get_active()) 
    WidthButton.set_value(HeightButton.get_value()*Ratio);
}

void 
ResExportClass::OnWidthAdjust()
{
  if (RatioButton.get_active()) 
    HeightButton.set_value(WidthButton.get_value()/Ratio);
}

void 
ResExportClass::OnRatioToggle()
{
  if (RatioButton.get_active())
    Ratio = WidthButton.get_value()/HeightButton.get_value();
}

void
ResExportClass::OnStillMovie()
{
  bool do_movie = !StillButton.get_active();
  MovieParamFrame.set_sensitive(do_movie);
}

void
ResExportClass::OnChooserChange()
{
  BaseNameEntry.set_text(BaseNameChooser.get_filename());
}

void
ResExportClass::OnEntryChange()
{
  BaseNameChooser.select_filename(BaseNameEntry.get_text());
}


void
ResExportClass::OnTypeChange()
{
  if (TypeCombo.get_active_text() == "POVray") 
    POVFrame.set_sensitive(true);
  else 
    POVFrame.set_sensitive(false);
}
