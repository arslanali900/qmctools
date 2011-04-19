#ifndef WF_VIS_H
#define WF_VIS_H

#include "PathVis.h"
#include "BoxClass.h"
#include "CoordObject.h"
#include "Isosurface.h"
#include "PlaneObject.h"
#include "CylinderObject.h"
#include "WFExport.h"
#include <Common/IO/IO.h>
#include <gtkmm/adjustment.h>


using namespace IO;

typedef enum {REAL_PART, IMAG_PART, MAG2} WFDisplayType;

class WFVisualClass;
class BandRow
{
public:
  Gtk::CheckButton Check;
  Gtk::Label kLabel, BandLabel;
  Isosurface *Iso;
  BandRow(WFVisualClass &wfvis) : Iso(NULL)
  {
    Check.set_active(false);
  }
};

class WFVisualClass : public Gtk::Window
{
protected:
  //////////
  // Data //
  //////////
  Array<double,3> WFData;
  int CurrBand, Currk, NumBands, Numk, NumElectrons;
  Vec3 SuperTwistInt;
  Array<Vec3,1> AtomPos;
  Array<int,1> AtomTypes;
  // First index in the walker number, the second is the
  // electron number
  Array<Vec3,2> Walkers;
  BoxClass Box;
  int Currentk, CurrentBand;
  WFDisplayType WFDisplay;
  // This stores the amount to shift the coordinates by
  Vec3 Shift;
  bool DoShift;

  /////////////
  // Widgets //
  /////////////
  Gtk::VBox MainVBox;
  Gtk::Button QuitButton;
  Gtk::HScale kScale, BandScale;
  Gtk::Adjustment kAdjust, BandAdjust;
  Gtk::Frame kFrame, BandFrame;
  Glib::RefPtr<Gtk::ToggleAction> CoordToggle, SphereToggle, BoxToggle,
    TruncRadiiToggle, IsocontourToggle, FullscreenToggle, BondsToggle;
;
  Gtk::RadioButtonGroup DisplayGroup, ColorMapGroup;
  Glib::RefPtr<Gtk::RadioAction> RealRadio, ImagRadio, Mag2Radio;
  
  Gtk::Toolbar Tools;
  Gtk::HBox ToolBox;
  Gtk::RadioToolButton OrthoButton, PerspectButton;
  Gtk::Image OrthoImage, PerspectImage;
  Gtk::ToggleToolButton ClipButton;
  Gtk::Image ClipImage;
  Gtk::ToggleToolButton IsoButton;
  Gtk::Image IsoImage;

  Glib::RefPtr<Gtk::ActionGroup> Actions;
  Glib::RefPtr<Gtk::UIManager> Manager;
  Gtk::HBox MiddleBox;
  Gtk::VBox OptionsBox;
  Gtk::Frame RadiusFrame;
  Gtk::HScale RadiusScale;
  Gtk::Adjustment RadiusAdjust;
  Gtk::HBox RadiusBox;

  //////////////////////
  // Isosurface stuff //
  //////////////////////
  IOSectionClass Infile;
  bool FileIsOpen;
  Array<double,3> RhoData;
  // Localize orbitals
  bool Localized, Truncated;
  Vec3 Center, uMin, uMax, uCenter;
  double TruncRadius;
  bool Nonuniform;
  
  Isosurface WFIso;
  LinearGrid Xgrid, Ygrid, Zgrid;
  GeneralGrid NUXgrid, NUYgrid, NUZgrid;
  Gtk::VBox IsoBox, DensityBox;
  Gtk::Label rsLabel;
  Gtk::HScale IsoScale;
  Gtk::Adjustment IsoAdjust;
  Gtk::Frame IsoFrame;
  double FindMaxVal();
  double MaxVal;
  double FindMaxBandVal();
  double MaxBandVal;

  ///////////////////////
  // Color plane stuff //
  ///////////////////////
  PlaneObject xPlane, yPlane, zPlane;
  Gtk::Frame PlaneFrame;
  Gtk::VBox PlaneBox;
  Gtk::Adjustment xPlaneAdjust, yPlaneAdjust, zPlaneAdjust;
  Gtk::HScale xPlaneScale, yPlaneScale, zPlaneScale;
  Gtk::CheckButton xPlaneButton, yPlaneButton, zPlaneButton;
  Gtk::HBox xPlaneBox, yPlaneBox, zPlaneBox;
  ColorMapType CMap;
  vector<Glib::RefPtr<Gtk::RadioAction> > CMapActions;

  //////////////////////////
  // Multi-band selection //
  //////////////////////////
  Gtk::Frame VisibleBandFrame;
  Gtk::VBox  VisibleBandBox;
  Gtk::ScrolledWindow VisibleBandWindow;
  Gtk::Table VisibleBandTable;
  vector<BandRow*> VisibleBandRows;
  Gtk::Label kLabel, BandLabel;
  Gtk::CheckButton MultiBandButton;
  void SetupBandTable();
  void SetupBandTable_ESHDF();
  void OnMultiBandToggle(), OnBandToggle(int row);
  void UpdateMultiIsos();

  /////////////////
  // State flags //
  /////////////////
  bool UpdateIso, UpdateIsoVal, UpdateIsoType, ResetIso;
  bool UpdatePlane[3];
  
  /////////////////////////////////////
  // Saving and opening viewer state //
  /////////////////////////////////////
  bool WriteState(string fname);
  void OnOpenState();
  void OnSaveState();
  Gtk::FileChooserDialog OpenStateChooser, SaveStateChooser;


  ////////////
  // Export //
  ////////////
  WFExportClass Export;

  //////////////////////
  // Callback methods //
  //////////////////////
  void OnExport();
  void Quit();
  void OnSpeedChange();
  void OnIsoChange();
  void OnBandChange();
  void OnkChange();
  void OnPlaneChange(int dir);
  void OnPerspectiveToggle();
  void OnPlayToggle();
  void OnClipToggle();
  void OnIsoToggle();
  void OnViewReset();
  void OnCoordToggle();
  void OnSphereToggle();
  void OnBondsToggle();
  void OnBoxToggle();
  void OnTruncRadiiToggle();
  void OnIsocontourToggle();
  void OnRadiusChange();
  void OnOpen();
  void OnDisplayRadio(WFDisplayType type);
  void OnColorMapRadio(ColorMapType type);
  void OnFullscreenToggle();
  bool UpToDate;

  ///////////////////
  // Other methods //
  ///////////////////
  string FindFullPath(string filename);
  bool ReadWF (int kPoint, int band);
  bool ReadWF_ESHDF (int kPoint, int band);
  bool IsESHDF;
public:
  PathVisClass PathVis;

  ////////////////////
  // Public methods //
  ////////////////////
  void Read(string filename);  
  void ReadWalkers(string filename);  
  void Read_ESHDF();
  bool ReadState (string fname);
  bool DrawFrame(bool offScreen=false);
  void SetShift (Vec3 shift);
  void SetViewportSize (int size);
  
  
  WFVisualClass ();
  virtual ~WFVisualClass();
};

#endif
