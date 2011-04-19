#ifndef WF_EXPORT_H
#define WF_EXPORT_H

#include <gtkmm.h>
#include <gtkglmm.h>

using namespace std;

class WFVisualClass;

class WFExportClass : public Gtk::Window
{
private:
  WFVisualClass &Visual;
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
  
  /// POV parameters
  Gtk::Frame POVFrame;
  Gtk::VBox  POVVBox1, POVVBox2;
  Gtk::HBox  POVHBox;
  Gtk::CheckButton POVRenderButton;
  Gtk::CheckButton POVAntiAliasButton;
  Gtk::CheckButton POVOnScreenButton;
  Gtk::CheckButton POVTransparentBGButton;
  Gtk::SpinButton  POVTolerance;
  Gtk::Adjustment  POVTolAdjust;
  Gtk::Label       POVTolLabel;
  string POVFile;
  void RenderPOV();

  /// Signal handlers
  void OnExportButton(), OnCancelButton(), OnBrowseButton();
  void OnWidthAdjust(), OnHeightAdjust(), OnRatioToggle();
  void OnStillMovie(), OnChooserChange(), OnEntryChange();
  void OnTypeChange();

  int Width, Height;
  double Ratio;
  int AttribList;
  void InitGLStuff();
  void MakePixmap();
public:
  void Export (string filename);
  void ExportPOV (string filename);
  void SetupWidgets();
  WFExportClass (WFVisualClass &visual);
};


#endif
