#ifndef EXPORT_H
#define EXPORT_H

#include <gtkmm.h>
#include <gtkglmm.h>

using namespace std;

class VisualClass;

class ExportClass : public Gtk::Window
{
private:
  VisualClass &Visual;
  /// These are used for GL-based exporting
  Glib::RefPtr<Gdk::GL::Pixmap> GLPixmap;
  Glib::RefPtr<Gdk::Pixmap>     GdkPixmap;
  Glib::RefPtr<Gdk::Drawable>   Draw;
  Glib::RefPtr<Gdk::GL::Config> GLConfig;
  Glib::RefPtr<Gdk::GL::Context> GLContext;
  //////////////////////////
  /// Widgets for the UI ///
  //////////////////////////
  Gtk::VBox MainVBox;

  /// Size adjustment
  Gtk::SpinButton WidthButton, HeightButton;
  Gtk::HBox WidthBox, HeightBox, WidthHeightBox;
  Gtk::VBox SizeBox;
  Gtk::Frame SizeFrame;
  Gtk::Label WidthLabel, HeightLabel;
  Gtk::CheckButton RatioButton;

  /// Export type
  Gtk::Frame TypeFrame;
  Gtk::HBox TypeBox;
  Gtk::ComboBoxText TypeCombo;
  Gtk::CheckButton AntialiasButton;

  /// Export and cancel buttons
  Gtk::HButtonBox ButtonBox;
  Gtk::Button ExportButton, CancelButton;

  /// File base
  Gtk::Button BrowseButton;
  Gtk::Frame BaseNameFrame;
  Gtk::HBox BaseNameHBox;
  Gtk::Entry BaseNameEntry;
  Gtk::Button BaseNameBrowseButton;
  Gtk::FileChooserDialog BaseNameChooser;

  /// Still/Video selection
  Gtk::RadioButton StillButton;
  Gtk::RadioButton MovieButton;
  Gtk::Frame StillMovieFrame;
  Gtk::VBox StillMovieBox;

  /// Movie parameters
  Gtk::Frame MovieParamFrame;
  Gtk::SpinButton FirstFrameButton, LastFrameButton;
  Gtk::HBox FirstFrameBox, LastFrameBox, InterpFactorBox;
  Gtk::Label FirstFrameLabel, LastFrameLabel, InterpFactorLabel;
  Gtk::HBox FirstLastBox, InterpBox;
  Gtk::VBox MovieParamBox;
  Gtk::SpinButton InterpFactorButton;
  

  /// Signal handlers
  void OnExportButton(), OnCancelButton(), OnBrowseButton();
  void OnWidthAdjust(), OnHeightAdjust(), OnRatioToggle();
  void OnStillMovie(), OnChooserChange(), OnEntryChange();

  int Width, Height;
  double Ratio;
  int AttribList;
  void InitGLStuff();
  void MakePixmap(int frame);
public:
  void Export (string filename);
  void ExportPOV (string filename);
  void ExportMovie (string filename, int firstFrame, int lastFrame,
		    int interpFactor);
  void SetupWidgets();
  ExportClass (VisualClass &visual);
};


#endif
