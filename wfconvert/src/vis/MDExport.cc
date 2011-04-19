#include "MDExport.h"
#include "MDVis.h"
#include <GL/gl.h>
#include <GL/glu.h>
#ifdef HAVE_LIBREVEL
  #include <revel.h>
#endif
#include <Common/Splines/CubicSpline.h>

using namespace IO;


MDExportClass::MDExportClass(MDVisualClass &mdvisual) :
    MDVisual(mdvisual), Width(2000), Height(2000),
    BaseNameChooser("Export filename", Gtk::FILE_CHOOSER_ACTION_SAVE),
    WidthLabel("Width:"), HeightLabel("Height:"),
    GLContextSetup(false)
{
  BaseNameFrame.set_label("Base filename");
  BaseNameFrame.add (BaseNameHBox);
  BaseNameHBox.pack_start(BaseNameEntry, Gtk::PACK_SHRINK, 5);
  BaseNameHBox.pack_start(BaseNameBrowseButton, Gtk::PACK_SHRINK);
  BaseNameBrowseButton.set_label("Browse");
  MainVBox.pack_start (BaseNameFrame, Gtk::PACK_SHRINK, 5);
  
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
  
  TypeFrame.set_label("Rendering Type");
  TypeFrame.add(TypeBox);
  TypeBox.pack_start(TypeCombo, Gtk::PACK_SHRINK);
  TypeCombo.insert_text(0, "Offscreen GL");
  TypeCombo.insert_text(1, "POVray");
  TypeCombo.set_active (0);
  ExportButton.set_label ("Export");
  CancelButton.set_label ("Cancel");
  
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
  MainVBox.pack_start (StillMovieFrame, Gtk::PACK_SHRINK, 5);
  MainVBox.pack_start (MovieParamFrame, Gtk::PACK_SHRINK, 5);
  MainVBox.pack_start (ButtonBox,       Gtk::PACK_SHRINK, 10);
  add (MainVBox);
  set_title ("Export");
  // Set signal handlers
  ExportButton.signal_clicked().connect
    (sigc::mem_fun(*this, &MDExportClass::OnExportButton));
  CancelButton.signal_clicked().connect
    (sigc::mem_fun(*this, &MDExportClass::OnCancelButton));
  BaseNameBrowseButton.signal_clicked().connect
    (sigc::mem_fun(*this, &MDExportClass::OnBrowseButton));
  RatioButton.signal_toggled().connect
    (sigc::mem_fun(*this, &MDExportClass::OnRatioToggle));
  WidthButton.signal_value_changed().connect
    (sigc::mem_fun(*this, &MDExportClass::OnWidthAdjust));
  HeightButton.signal_value_changed().connect
    (sigc::mem_fun(*this, &MDExportClass::OnHeightAdjust));
  StillButton.signal_toggled().connect
    (sigc::mem_fun(*this, &MDExportClass::OnStillMovie));
  BaseNameChooser.signal_selection_changed().connect
    (sigc::mem_fun(*this, &MDExportClass::OnChooserChange));
  BaseNameEntry.signal_activate().connect
    (sigc::mem_fun(*this, &MDExportClass::OnEntryChange));

#ifndef HAVE_LIBREVEL
  MovieButton.set_sensitive(false);
#endif
}


void MDExportClass::SetupWidgets()
{
  FirstFrameButton.set_digits(0);
  FirstFrameButton.set_range(1.0, MDVisual.Trajectory.extent(0)+1);
  FirstFrameButton.set_increments (1.0, 10.0);
  LastFrameButton.set_digits(0);
  LastFrameButton.set_range (1.0, MDVisual.Trajectory.extent(0)+1);
  LastFrameButton.set_increments (1.0, 10.0);
  InterpFactorButton.set_range(1.0, 500.0);
  InterpFactorButton.set_digits(0);
  InterpFactorButton.set_increments (1.0, 10.0);
}  

void MDExportClass::InitGLStuff()
{
  glShadeModel(GL_SMOOTH);
  glEnable (GL_LIGHTING);
  glEnable (GL_LINE_SMOOTH);
  glEnable (GL_POLYGON_SMOOTH);

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


void 
MDExportClass::MakePixmap (int frame)
{
  if (!GLContextSetup) {
    cerr << "Setting up GL context.\n";
    GLContextSetup = true;
    GdkPixmap.clear();
    GLConfig.clear();
    
    GLConfig = Gdk::GL::Config::create (Gdk::GL::MODE_RGB   |
					Gdk::GL::MODE_DEPTH |
					Gdk::GL::MODE_SINGLE);
    if (!GLConfig) {
      cerr << "Cannot find a visual capable of OpenGL.\n";
      return;
    }
    
    cerr << "Width = " << Width << "   Height = " << Height << endl;
    GdkPixmap = Gdk::Pixmap::create (MDVisual.get_window(),
				     Width, Height, GLConfig->get_depth());
    
    GLPixmap = 
      Gdk::GL::ext(GdkPixmap).set_gl_capability (GLConfig);

    GLContext = Gdk::GL::Context::create (GLPixmap, false);
  }
  assert (GdkPixmap);
  assert (GLPixmap);
  assert (GLContext);

  GLPixmap->make_current(GLContext);
  GLPixmap->gl_begin(GLContext);

  // Render here
  InitGLStuff();
  int saveFrame = MDVisual.CurrentFrame;
  MDVisual.CurrentFrame = frame;
  MDVisual.DrawFrame (true);
  MDVisual.PathVis.GLRender();
  MDVisual.CurrentFrame = saveFrame;

  glFlush();
  GLPixmap->gl_end();
  GLPixmap->wait_gl();
}

#include <stdlib.h>
//Glib::RefPtr<Gdk::Pixbuf>
void
MDExportClass::MakePOVPixbuf(int frame)
{
  MDVisual.CurrentFrame = frame;
  MDVisual.DrawFrame();
  MDVisual.PathVis.POVRender("POVframe.pov");
  char povString[1000];
  snprintf (povString, 1000, "povray -D +A0.01 +W%d +H%d POVframe.pov", 
	    Width, Height);
  system(povString);
  Pixbuf.clear();
  Pixbuf = Gdk::Pixbuf::create_from_file ("POVframe.png", Width, Height);
  guint8 *data = Pixbuf->get_pixels();
  /// Swap from RGB to BGR format
  for (int i=0; i<Width*Height; i++) {
    guint8 tmp    = *(data+3*i);
    *(data+3*i)   = *(data+3*i+2);
    *(data+3*i+2) = tmp;
  }
}

void 
MDExportClass::Export (string filename)
{
  MakePixmap ((int) floor (MDVisual.FrameAdjust.get_value()));
  Glib::RefPtr<Gdk::Drawable> drawable =  GdkPixmap;
  Glib::RefPtr<Gdk::Drawable> drawable2 =  GdkPixmap;

  cerr << "drawable = " << drawable << endl;

  Glib::RefPtr<Gdk::Pixbuf> Pbuf = Gdk::Pixbuf::create
    (drawable, GdkPixmap->get_colormap(), 0, 0, 0, 0, Width, Height);

  if (drawable == drawable2)
    cerr << "no change\n";
						  
  Pbuf->save(filename, "png");
}



void
MDExportClass::ExportPOV(string basename)
{
  string filename = basename;
  filename.append(".pov");
  MDVisual.PathVis.POVRender (filename);
}

void
MDExportClass::ExportMovie (string basename, 
			    int firstFrame, int lastFrame, 
			    int interpFactor)
{
#ifdef HAVE_LIBREVEL
  ///////////////////
  // Encoder setup //
  ///////////////////
  int encoderHandle;
  Revel_CreateEncoder(&encoderHandle);
  Revel_Params params;
  Revel_InitializeParams(&params);
  params.width  = Width;
  params.height = Height;
  params.frameRate = 60.0;
  params.quality = 1.0;
  params.codec = REVEL_CD_XVID;
  params.hasAudio = 0;
  basename.append(".avi");
  cerr << "filename = " << basename << endl;
  Revel_EncodeStart(encoderHandle, basename.c_str(), &params);
  Revel_VideoFrame videoFrame;
  videoFrame.width  = Width;
  videoFrame.height = Height;

  ////////////////
  // Data Setup //
  ////////////////
  for (int frame=0; frame < MDVisual.Trajectory.extent(0); frame++) {
  
    if (TypeCombo.get_active_text() == "POVray") {
      if (MDVisual.IsoButton.get_active()) {
	MDVisual.RhoVar->Read(MDVisual.RhoData, frame, 
			      Range::all(), Range::all(), Range::all()); 
	MDVisual.MaxRho = MDVisual.FindMaxRho();
      }
      videoFrame.bytesPerPixel = 3;  
      videoFrame.pixelFormat = REVEL_PF_BGR;
      MakePOVPixbuf(frame);
      videoFrame.pixels = Pixbuf->get_pixels();
    }
    else {
      videoFrame.bytesPerPixel = 4;  
      videoFrame.pixelFormat = REVEL_PF_BGRA;
      MakePixmap(frame);
      Glib::RefPtr<Gdk::Image> image = GdkPixmap->get_image(0,0,Width,Height);
      videoFrame.pixels = image->get_mem();
    }
    int frameSize;
    Revel_EncodeFrame(encoderHandle, &videoFrame, &frameSize);
    cerr << "After encoding frame.\n";
//     char fname[100];
//     snprintf (fname, 100, "frame%d.png", frame);
//     Glib::RefPtr<Gdk::Drawable> drawable =  GdkPixmap;
//     Glib::RefPtr<Gdk::Pixbuf> Pbuf = Gdk::Pixbuf::create
//       (drawable, GdkPixmap->get_colormap(), 0, 0, 0, 0, Width, Height);
//     Pbuf->save(fname, "png");
//     Pbuf.clear();
    

  }
  int totalSize;
  Revel_EncodeEnd(encoderHandle, &totalSize);
  Revel_DestroyEncoder(encoderHandle);
#endif
}


void
MDExportClass::OnExportButton()
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
MDExportClass::OnCancelButton()
{
  hide();
}

void
MDExportClass::OnBrowseButton()
{
  BaseNameChooser.run();
}

void 
MDExportClass::OnHeightAdjust()
{
  if (RatioButton.get_active()) 
    WidthButton.set_value(HeightButton.get_value()*Ratio);
}

void 
MDExportClass::OnWidthAdjust()
{
  if (RatioButton.get_active()) 
    HeightButton.set_value(WidthButton.get_value()/Ratio);
}

void 
MDExportClass::OnRatioToggle()
{
  if (RatioButton.get_active())
    Ratio = WidthButton.get_value()/HeightButton.get_value();
}

void
MDExportClass::OnStillMovie()
{
  bool do_movie = !StillButton.get_active();
  MovieParamFrame.set_sensitive(do_movie);
}

void
MDExportClass::OnChooserChange()
{
  BaseNameEntry.set_text(BaseNameChooser.get_filename());
}

void
MDExportClass::OnEntryChange()
{
  BaseNameChooser.select_filename(BaseNameEntry.get_text());
}
