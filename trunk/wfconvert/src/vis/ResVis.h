#ifndef RES_VIS_H
#define RES_VIS_H

#include "PathVis.h"
#include "BoxClass.h"
#include "CoordObject.h"
#include "Isosurface.h"
#include "PlaneObject.h"
#include "CylinderObject.h"
#include "ResExport.h"
#include <Common/IO/IO.h>
#include <gtkmm/adjustment.h>


using namespace IO;

typedef enum {PRESSURE, SATURATION} PropertyType;
typedef enum {GAS, OIL, WATER} FluidType;


class ResVisualClass;

class ResVisualClass : public Gtk::Window
{
protected:
  //////////
  // Data //
  //////////
  Array<double,3> GasData, OilData, WaterData;
  int CurrentStep, CurrentFluid, NumSteps;
  Array<Vec3,1> AtomPos;
  Array<int,1> AtomTypes;
  BoxClass Box;
  PropertyType DisplayProperty;
  // This stores the amount to shift the coordinates by
  Vec3 Shift;
  bool DoShift;

  /////////////
  // Widgets //
  /////////////
  Gtk::VBox MainVBox;
  Gtk::Button QuitButton;
  Gtk::HScale StepScale;
  Gtk::Adjustment StepAdjust;
  Gtk::Frame StepFrame;
  Glib::RefPtr<Gtk::ToggleAction> CoordToggle, SphereToggle, BoxToggle,
    TruncRadiiToggle, IsocontourToggle, FullscreenToggle, BondsToggle;
  Gtk::RadioButtonGroup DisplayGroup, ColorMapGroup;
  Glib::RefPtr<Gtk::RadioAction> SaturationRadio, PressureRadio;
  
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
  Vec3 Center, uMin, uMax, uCenter;
  
  Isosurface ResIso;
  LinearGrid Xgrid, Ygrid, Zgrid;
  GeneralGrid NUXgrid, NUYgrid, NUZgrid;
  Gtk::VBox IsoBox, DensityBox;
  Gtk::HScale IsoScale;
  Gtk::Adjustment IsoAdjust;
  Gtk::Frame IsoFrame;
  double FindMaxVal();
  double MinVal, MaxVal;
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
  ///////////////////
  // Fluid control //
  ///////////////////
  Gtk::Frame FluidFrame;
  Gtk::VBox FluidBox;
  Gtk::Adjustment GasAdjust, OilAdjust, WaterAdjust;
  Gtk::HScale GasScale, OilScale, WaterScale;
  Gtk::CheckButton GasButton, OilButton, WaterButton;
  Gtk::HBox GasBox, OilBox, WaterBox;

  /////////////////
  // State flags //
  /////////////////
  bool UpdateIso, UpdateIsoVal, UpdateIsoType, ResetIso;
  bool UpdatePlane[3], UpdateGas, UpdateOil, UpdateWater;
  
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
  ResExportClass Export;

  //////////////////////
  // Callback methods //
  //////////////////////
  void OnExport();
  void Quit();
  void OnSpeedChange();
  void OnIsoChange();
  void OnStepChange();
  void OnPlaneChange(int dir);
  void OnFluidChange(FluidType fl);
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
  void OnDisplayRadio(PropertyType type);
  void OnColorMapRadio(ColorMapType type);
  void OnFullscreenToggle();
  bool UpToDate;

  ///////////////////
  // Other methods //
  ///////////////////
  string FindFullPath(string filename);
  bool ReadRes (int step, int fluid, PropertyType prop);
public:
  PathVisClass PathVis;

  ////////////////////
  // Public methods //
  ////////////////////
  void Read(string filename);  
  bool ReadState (string fname);
  bool DrawFrame(bool offScreen=false);
  void SetShift (Vec3 shift);
  void SetViewportSize (int size);
  void SetStep (int step);
  int  GetNumSteps() { return NumSteps; }
  
  ResVisualClass ();
  virtual ~ResVisualClass();
};

#endif
