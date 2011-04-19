#include "Export.h"
#include "Visual.h"
#include <GL/gl.h>
#include <GL/glu.h>
#include <Common/Splines/CubicSpline.h>

using namespace IO;

ExportClass::ExportClass(VisualClass &visual) :
    Visual(visual), Width(2000), Height(2000),
    BaseNameChooser("Export filename", Gtk::FILE_CHOOSER_ACTION_SAVE),
    WidthLabel("Width:"), HeightLabel("Height:")
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
    (sigc::mem_fun(*this, &ExportClass::OnExportButton));
  CancelButton.signal_clicked().connect
    (sigc::mem_fun(*this, &ExportClass::OnCancelButton));
  BaseNameBrowseButton.signal_clicked().connect
    (sigc::mem_fun(*this, &ExportClass::OnBrowseButton));
  RatioButton.signal_toggled().connect
    (sigc::mem_fun(*this, &ExportClass::OnRatioToggle));
  WidthButton.signal_value_changed().connect
    (sigc::mem_fun(*this, &ExportClass::OnWidthAdjust));
  HeightButton.signal_value_changed().connect
    (sigc::mem_fun(*this, &ExportClass::OnHeightAdjust));
  StillButton.signal_toggled().connect
    (sigc::mem_fun(*this, &ExportClass::OnStillMovie));
  BaseNameChooser.signal_selection_changed().connect
    (sigc::mem_fun(*this, &ExportClass::OnChooserChange));
  BaseNameEntry.signal_activate().connect
    (sigc::mem_fun(*this, &ExportClass::OnEntryChange));
}


void ExportClass::SetupWidgets()
{
  FirstFrameButton.set_digits(0);
  FirstFrameButton.set_range(1.0, Visual.PathArray.extent(0)+1);
  FirstFrameButton.set_increments (1.0, 10.0);
  LastFrameButton.set_digits(0);
  LastFrameButton.set_range (1.0, Visual.PathArray.extent(0)+1);
  LastFrameButton.set_increments (1.0, 10.0);
  InterpFactorButton.set_range(1.0, 500.0);
  InterpFactorButton.set_digits(0);
  InterpFactorButton.set_increments (1.0, 10.0);

#ifndef HAVE_LIBREVEL
  MovieButton.set_sensitive(false);
#endif
}  

void ExportClass::InitGLStuff()
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


void ExportClass::MakePixmap (int frame)
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
  Visual.MakeFrame (frame, true);
  Visual.PathVis.GLRender();

  glFlush();
  GLPixmap->gl_end();
  GLPixmap->wait_gl();
}

void ExportClass::Export (string filename)
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

  MakePixmap ((int) floor (Visual.FrameAdjust.get_value()));
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
ExportClass::ExportPOV(string basename)
{
  string filename = basename;
  filename.append(".pov");
  Visual.PathVis.POVRender (filename);
}

#ifdef HAVE_LIBREVEL
#include <revel.h>
#endif

void
ExportClass::ExportMovie (string basename, 
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
  params.frameRate = 30.0;
  params.quality = 1.0;
  params.codec = REVEL_CD_XVID;
  params.hasAudio = 0;
  basename.append(".avi");
  Revel_EncodeStart(encoderHandle, basename.c_str(), &params);
  Revel_VideoFrame videoFrame;
  videoFrame.width  = Width;
  videoFrame.height = Height;
  videoFrame.bytesPerPixel = 4;
  videoFrame.pixelFormat = REVEL_PF_BGRA;

  ////////////////
  // Data Setup //
  ////////////////
  Array<double,4> &PathArray = Visual.PathArray;
  Array<double,3> &NodeData = Visual.ANodeData;
  int numPaths = Visual.PathArray.extent(0) - 1;
  LinearGrid tGrid (0.0, (double)(numPaths-1), numPaths);
  CubicSpline interpSpline;
  int numFrames = interpFactor*(lastFrame-firstFrame)+1;
  Array<double,1> pathData(numPaths);
  for (int frame=0; frame<numFrames; frame++) {
    /// First interpolate
    for (int ptcl=0; ptcl<PathArray.extent(1); ptcl++){
      for (int slice=0; slice<PathArray.extent(2); slice++)
	for (int dim=0; dim<3; dim++) {
	  for (int path=0; path<numPaths; path++) 
	    pathData(path) = PathArray(path, ptcl, slice, dim);
	  interpSpline.Init (&tGrid, pathData);
	  // Put interpolation in last slice
	  PathArray(numPaths, ptcl, slice, dim) = 
	    interpSpline((double)frame/(double)interpFactor+firstFrame);
	}
    }
    if (Visual.OpenPtcl(0)!=-1)
	for (int dim=0; dim<3; dim++) {
	  for (int path=0; path<numPaths; path++) 
	    pathData(path) = Visual.Tail(path, dim);
	  interpSpline.Init (&tGrid, pathData);
	  // Put interpolation in last slice
	  Visual.Tail(numPaths, dim) = 
	    interpSpline((double)frame/(double)interpFactor+firstFrame);
	}
    
    /// Interpolate node data
    for (int ix=0; ix<NodeData.extent(0); ix++)
      for (int iy=0; iy<NodeData.extent(1); iy++)
	for (int iz=0; iz<NodeData.extent(2); iz++) {
	  Visual.ANodeVar->Read(pathData, Range::all(), ix, iy, iz);
	  interpSpline.Init(&tGrid,pathData);
	  // Put interpolation in last slice
	  Visual.ANodeData(ix, iy, iz) =
	    interpSpline((double)frame/(double)interpFactor+firstFrame);
	  if (Visual.HaveBNodeData) {
	    Visual.BNodeVar->Read(pathData, Range::all(), ix, iy, iz);
	    interpSpline.Init(&tGrid,pathData);
	    // Put interpolation in last slice
	    Visual.BNodeData(ix, iy, iz) =
	      interpSpline((double)frame/(double)interpFactor+firstFrame);
	  }
	}

    /// Now, create image from interpolated frame
    MakePixmap(numPaths);
    Glib::RefPtr<Gdk::Image> image = GdkPixmap->get_image(0,0,Width,Height);
    videoFrame.pixels = image->get_mem();
    cerr << "Bits per pixel = " << image->get_bits_per_pixel() << endl;
    int frameSize;
    Revel_EncodeFrame(encoderHandle, &videoFrame, &frameSize);
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
ExportClass::OnExportButton()
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
ExportClass::OnCancelButton()
{
  hide();
}

void
ExportClass::OnBrowseButton()
{
  BaseNameChooser.run();
}

void 
ExportClass::OnHeightAdjust()
{
  if (RatioButton.get_active()) 
    WidthButton.set_value(HeightButton.get_value()*Ratio);
}

void 
ExportClass::OnWidthAdjust()
{
  if (RatioButton.get_active()) 
    HeightButton.set_value(WidthButton.get_value()/Ratio);
}

void 
ExportClass::OnRatioToggle()
{
  if (RatioButton.get_active())
    Ratio = WidthButton.get_value()/HeightButton.get_value();
}

void
ExportClass::OnStillMovie()
{
  bool do_movie = !StillButton.get_active();
  MovieParamFrame.set_sensitive(do_movie);
}

void
ExportClass::OnChooserChange()
{
  BaseNameEntry.set_text(BaseNameChooser.get_filename());
}

void
ExportClass::OnEntryChange()
{
  BaseNameChooser.select_filename(BaseNameEntry.get_text());
}
