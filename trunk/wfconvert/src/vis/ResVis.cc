#include "ResVis.h"
#include <GL/glut.h>
#include "ElementData.h"
#include "ParseCommand.h"
#include "PointsObject.h"



ResVisualClass::ResVisualClass() :
  MainVBox(false, 0), 
  QuitButton("Quit"),
  IsoAdjust  (0.01, 0.0, 1.0, 0.01, 0.1),
  StepAdjust (0.0, 0.0, 8.0, 1.0, 1.0),
  xPlaneAdjust(0.0, 0.0, 1.0, 0.01, 0.1),
  yPlaneAdjust(0.0, 0.0, 1.0, 0.01, 0.1),
  zPlaneAdjust(0.0, 0.0, 1.0, 0.01, 0.1),
  GasAdjust(0.0, 0.0, 1.0, 0.01, 0.1),
  OilAdjust(0.0, 0.0, 1.0, 0.01, 0.1),
  WaterAdjust(0.0, 0.0, 1.0, 0.01, 0.1),
  RadiusAdjust(0.4, 0.0, 1.0, 0.01, 0.1),
  UpToDate(true),
  FileIsOpen(false),
  xPlane(ResIso),
  yPlane(ResIso),
  zPlane(ResIso),
  Export(*this),
  DisplayProperty(SATURATION),
  ResetIso(false),
  SaveStateChooser ("State filename", Gtk::FILE_CHOOSER_ACTION_SAVE),
  OpenStateChooser ("State filename", Gtk::FILE_CHOOSER_ACTION_OPEN),
  UpdateIsoType(false), 
  UpdateIsoVal(false),
  DoShift (false),
  Shift (0.0, 0.0, 0.0),
  CMap(BLUE_WHITE_RED),
  uCenter(0.0, 0.0, 0.0),
  uMin(0.0, 0.0, 0.0),
  uMax(1.0, 1.0, 1.0)
{
  if (!Glib::thread_supported())
    Glib::thread_init();
  ResIso.Dynamic = false;
  xPlane.Dynamic = false;
  yPlane.Dynamic = false;
  zPlane.Dynamic = false;
  IsoScale.set_adjustment(IsoAdjust);
  IsoScale.set_digits(2);
  IsoScale.signal_value_changed().connect
    (sigc::mem_fun(*this, &ResVisualClass::OnIsoChange));
  IsoFrame.set_label("Density");
  IsoFrame.add(DensityBox);
  DensityBox.pack_start(IsoScale);
  IsoScale.property_width_request().set_value(75);

  StepScale.set_adjustment (StepAdjust);
  StepScale.set_digits(0);
  StepScale.signal_value_changed().connect
    (sigc::mem_fun(*this, &ResVisualClass::OnStepChange));
  StepScale.property_width_request().set_value(75);
  StepFrame.set_label("Step");
  StepFrame.add(StepScale);

  RadiusScale.set_adjustment(RadiusAdjust);
  RadiusScale.set_digits(2);
  RadiusScale.signal_value_changed().connect
    (sigc::mem_fun(*this, &ResVisualClass::OnRadiusChange));
  RadiusScale.property_width_request().set_value(75);
  RadiusFrame.set_label("Ion radii");
  RadiusFrame.add(RadiusScale);
  RadiusBox.pack_start(RadiusFrame, Gtk::PACK_SHRINK, 5);
  OptionsBox.pack_start(RadiusBox,  Gtk::PACK_SHRINK, 5);

  OrthoImage.set(FindFullPath("orthographic.png"));
  OrthoButton.set_icon_widget(OrthoImage);
  OrthoButton.set_label("Ortho");
  PerspectImage.set(FindFullPath("perspective.png"));
  PerspectButton.set_icon_widget(PerspectImage);
  PerspectButton.set_label("Perspect");

  ClipImage.set(FindFullPath("clipping.png"));
  ClipButton.set_icon_widget(ClipImage);
  ClipButton.set_label("Clip");

  IsoImage.set(FindFullPath("isoButton.png"));
  IsoButton.set_icon_widget(IsoImage);
  IsoButton.set_label("Isosurf");


  //////////////////////////////
  // Color plane widget setup //
  //////////////////////////////
  xPlaneScale.set_value_pos(Gtk::POS_RIGHT);
  yPlaneScale.set_value_pos(Gtk::POS_RIGHT);
  zPlaneScale.set_value_pos(Gtk::POS_RIGHT);
  xPlaneScale.set_adjustment(xPlaneAdjust);
  yPlaneScale.set_adjustment(yPlaneAdjust);
  zPlaneScale.set_adjustment(zPlaneAdjust);
  xPlaneScale.property_width_request().set_value(75);
  yPlaneScale.property_width_request().set_value(75);
  zPlaneScale.property_width_request().set_value(75);
  xPlaneScale.set_digits(2);
  yPlaneScale.set_digits(2);
  zPlaneScale.set_digits(2);
  xPlaneButton.set_label("x Plane"); 
  ((std::vector<Gtk::Widget*>)xPlaneButton.get_children())[0]->
    modify_fg(Gtk::STATE_NORMAL, Gdk::Color("red"));
  yPlaneButton.set_label("y Plane");
  ((std::vector<Gtk::Widget*>)yPlaneButton.get_children())[0]->
    modify_fg(Gtk::STATE_NORMAL, Gdk::Color("green"));
  zPlaneButton.set_label("z Plane");
  ((std::vector<Gtk::Widget*>)zPlaneButton.get_children())[0]->
    modify_fg(Gtk::STATE_NORMAL, Gdk::Color("blue"));
  xPlaneBox.pack_start (xPlaneButton, Gtk::PACK_SHRINK, 5);
  xPlaneBox.pack_start (xPlaneScale,  Gtk::PACK_SHRINK, 5);
  yPlaneBox.pack_start (yPlaneButton, Gtk::PACK_SHRINK, 5);
  yPlaneBox.pack_start (yPlaneScale,  Gtk::PACK_SHRINK, 5);
  zPlaneBox.pack_start (zPlaneButton, Gtk::PACK_SHRINK, 5);
  zPlaneBox.pack_start (zPlaneScale,  Gtk::PACK_SHRINK, 5);
  PlaneBox.pack_start(xPlaneBox, Gtk::PACK_SHRINK);
  PlaneBox.pack_start(yPlaneBox, Gtk::PACK_SHRINK);
  PlaneBox.pack_start(zPlaneBox, Gtk::PACK_SHRINK);
  PlaneFrame.add (PlaneBox);
  xPlaneScale.signal_value_changed().connect
    (sigc::bind<int>(sigc::mem_fun(*this, &ResVisualClass::OnPlaneChange), 0));
  yPlaneScale.signal_value_changed().connect
    (sigc::bind<int>(sigc::mem_fun(*this, &ResVisualClass::OnPlaneChange), 1));
  zPlaneScale.signal_value_changed().connect
    (sigc::bind<int>(sigc::mem_fun(*this, &ResVisualClass::OnPlaneChange), 2));
  xPlaneButton.signal_toggled().connect
    (sigc::bind<int>(sigc::mem_fun(*this, &ResVisualClass::OnPlaneChange), 0));
  yPlaneButton.signal_toggled().connect
    (sigc::bind<int>(sigc::mem_fun(*this, &ResVisualClass::OnPlaneChange), 1));
  zPlaneButton.signal_toggled().connect
    (sigc::bind<int>(sigc::mem_fun(*this, &ResVisualClass::OnPlaneChange), 2));

  ///////////////////////////////////
  // Fluid isosuface control setup //
  ///////////////////////////////////
  GasScale.set_value_pos(Gtk::POS_RIGHT);
  OilScale.set_value_pos(Gtk::POS_RIGHT);
  WaterScale.set_value_pos(Gtk::POS_RIGHT);
  GasScale.set_adjustment(GasAdjust);
  OilScale.set_adjustment(OilAdjust);
  WaterScale.set_adjustment(WaterAdjust);
  GasScale.property_width_request().set_value(75);
  OilScale.property_width_request().set_value(75);
  WaterScale.property_width_request().set_value(75);
  GasScale.set_digits(2);
  OilScale.set_digits(2);
  WaterScale.set_digits(2);
  GasButton.set_label("Gas"); 
  ((std::vector<Gtk::Widget*>)GasButton.get_children())[0]->
    modify_fg(Gtk::STATE_NORMAL, Gdk::Color("red"));
  OilButton.set_label("Oil");
  ((std::vector<Gtk::Widget*>)OilButton.get_children())[0]->
    modify_fg(Gtk::STATE_NORMAL, Gdk::Color("green"));
  WaterButton.set_label("Water");
  ((std::vector<Gtk::Widget*>)WaterButton.get_children())[0]->
    modify_fg(Gtk::STATE_NORMAL, Gdk::Color("blue"));
  GasBox.pack_start (GasButton, Gtk::PACK_SHRINK, 5);
  GasBox.pack_start (GasScale,  Gtk::PACK_SHRINK, 5);
  OilBox.pack_start (OilButton, Gtk::PACK_SHRINK, 5);
  OilBox.pack_start (OilScale,  Gtk::PACK_SHRINK, 5);
  WaterBox.pack_start (WaterButton, Gtk::PACK_SHRINK, 5);
  WaterBox.pack_start (WaterScale,  Gtk::PACK_SHRINK, 5);
  FluidBox.pack_start(GasBox, Gtk::PACK_SHRINK);
  FluidBox.pack_start(OilBox, Gtk::PACK_SHRINK);
  FluidBox.pack_start(WaterBox, Gtk::PACK_SHRINK);
  FluidFrame.add (FluidBox);
  GasScale.signal_value_changed().connect
    (sigc::bind<FluidType>(sigc::mem_fun(*this, &ResVisualClass::OnFluidChange), GAS));
  OilScale.signal_value_changed().connect
    (sigc::bind<FluidType>(sigc::mem_fun(*this, &ResVisualClass::OnFluidChange), OIL));
  WaterScale.signal_value_changed().connect
    (sigc::bind<FluidType>(sigc::mem_fun(*this, &ResVisualClass::OnFluidChange), WATER));
  GasButton.signal_toggled().connect
    (sigc::bind<FluidType>(sigc::mem_fun(*this, &ResVisualClass::OnFluidChange), GAS));
  OilButton.signal_toggled().connect
    (sigc::bind<FluidType>(sigc::mem_fun(*this, &ResVisualClass::OnFluidChange), OIL));
  WaterButton.signal_toggled().connect
    (sigc::bind<FluidType>(sigc::mem_fun(*this, &ResVisualClass::OnFluidChange), WATER));


  set_reallocate_redraws(true);
  PathVis.set_size_request(800, 800);
  ////////////////////
  // Setup tool bar //
  ////////////////////
  Gtk::RadioButtonGroup group = OrthoButton.get_group();
  PerspectButton.set_group(group);
  
  Tools.append(OrthoButton);
  Tools.append(PerspectButton);
  Tools.append(ClipButton);
  Tools.append(IsoButton);
  
  /////////////////
  // Setup menus //
  /////////////////
  Actions = Gtk::ActionGroup::create();
  Actions->add (Gtk::Action::create("MenuFile", "_File"));
  Actions->add (Gtk::Action::create("Open", "_Open"),
		sigc::mem_fun(*this, &ResVisualClass::OnOpen));
  Actions->add (Gtk::Action::create("Export", "_Export Image"),
		sigc::mem_fun(*this, &ResVisualClass::OnExport));
  Actions->add (Gtk::Action::create("OpenState", "Open State"),
		sigc::mem_fun(*this, &ResVisualClass::OnOpenState));
  Actions->add (Gtk::Action::create("SaveState", "Save State"),
		sigc::mem_fun(*this, &ResVisualClass::OnSaveState));
  Actions->add (Gtk::Action::create("Quit", "_Quit"),
		sigc::mem_fun(*this, &ResVisualClass::Quit));
  Actions->add (Gtk::Action::create("MenuView", "View"));
  Actions->add (Gtk::Action::create("Reset", "Reset"),
		sigc::mem_fun(*this, &ResVisualClass::OnViewReset));
  CoordToggle = Gtk::ToggleAction::create("Axes", "Coordinate axes",
					  "Show coordinate axes", true);
  Actions->add (CoordToggle,
		sigc::mem_fun(*this, &ResVisualClass::OnCoordToggle));
  SphereToggle = Gtk::ToggleAction::create("Nuclei", "Nuclei",
					  "Show nuclei", true);
  Actions->add (SphereToggle,
		sigc::mem_fun(*this, &ResVisualClass::OnSphereToggle));

  BondsToggle = Gtk::ToggleAction::create("Bonds", "Bonds",
					  "Show bonds", false);
  Actions->add (BondsToggle,
		sigc::mem_fun(*this, &ResVisualClass::OnBondsToggle));

  BoxToggle = Gtk::ToggleAction::create("Box", "Box", "Show box", true);
  Actions->add (BoxToggle,
		sigc::mem_fun(*this, &ResVisualClass::OnBoxToggle));

  TruncRadiiToggle = 
    Gtk::ToggleAction::create("Trunc", "Truncation radii",
			      "Show truncation radii", true);
  Actions->add (TruncRadiiToggle,
		sigc::mem_fun(*this, &ResVisualClass::OnTruncRadiiToggle));

  IsocontourToggle = 
    Gtk::ToggleAction::create("Isocontours", "Plane isocontours",
			      "Show isocontours on color planes", true);
  Actions->add (IsocontourToggle,
		sigc::mem_fun(*this, &ResVisualClass::OnIsocontourToggle));
  FullscreenToggle = Gtk::ToggleAction::create("Fullscreen", "Fullscreen mode",
					       "Display this window in fullscreen mode", false);
  Actions->add (FullscreenToggle,
		sigc::mem_fun(*this, &ResVisualClass::OnFullscreenToggle));

  Actions->add (Gtk::Action::create("MenuDisplay", "Display"));
  SaturationRadio = Gtk::RadioAction::create
    (DisplayGroup, "Saturation", "Saturation",
     "Display fluid saturation");
  PressureRadio = Gtk::RadioAction::create
    (DisplayGroup, "Pressure", "Pressure", "Display fluid pressure");
  Actions->add
    (SaturationRadio, sigc::bind<PropertyType> 
     (sigc::mem_fun(*this, &ResVisualClass::OnDisplayRadio),SATURATION));
  Actions->add
    (PressureRadio, sigc::bind<PropertyType> 
     (sigc::mem_fun(*this, &ResVisualClass::OnDisplayRadio),PRESSURE));

  ///////////////////
  // Colormap menu //
  ///////////////////
  vector<string> mapNames;
  mapNames.push_back ("Autumn");
  mapNames.push_back ("Bone");
  mapNames.push_back ("Colorcube");
  mapNames.push_back ("Cool");
  mapNames.push_back ("Copper");
  mapNames.push_back ("Flag");
  mapNames.push_back ("Gray");
  mapNames.push_back ("Hot");
  mapNames.push_back ("HSV");
  mapNames.push_back ("Jet");
  mapNames.push_back ("Lines");
  mapNames.push_back ("Pink");
  mapNames.push_back ("Prism");
  mapNames.push_back ("Spring");
  mapNames.push_back ("Summer");
  mapNames.push_back ("White");
  mapNames.push_back ("Winter");
  mapNames.push_back ("BlueWhiteRed");
  Actions->add (Gtk::Action::create("MenuColormap", "Colormap"));
  CMapActions.resize(BLUE_WHITE_RED+1);
  for (int i=0; i <= BLUE_WHITE_RED; i++) {
    CMapActions[i] = Gtk::RadioAction::create
      (ColorMapGroup, mapNames[i], mapNames[i]);
    Actions->add
      (CMapActions[i], sigc::bind<ColorMapType>
       (sigc::mem_fun (*this, &ResVisualClass::OnColorMapRadio),
	(ColorMapType)i));
  }

  Glib::ustring ui_info =
    "<ui>"
    "  <menubar name='MenuBar'>"
    "    <menu action='MenuFile'>"
    "      <menuitem action='Open'/>"
    "      <menuitem action='Export'/>"
    "      <menuitem action='SaveState'/>"
    "      <menuitem action='OpenState'/>"
    "      <separator/>"
    "      <menuitem action='Quit'/>"
    "    </menu>"
    "    <menu action='MenuView'>"
    "      <menuitem action='Reset'/>"
    "      <menuitem action='Nuclei'/>"
    "      <menuitem action='Bonds'/>"
    "      <menuitem action='Axes'/>"
    "      <menuitem action='Box'/>"
    "      <menuitem action='Trunc'/>"
    "      <menuitem action='Isocontours'/>"
    "      <menuitem action='Fullscreen'/>"
    "    </menu>"
    "    <menu action='MenuDisplay'>"
    "      <menuitem action='Saturation'/>"
    "      <menuitem action='Pressure'/>"
    "    </menu>"
    "    <menu action='MenuColormap'>"
    "       <menuitem action='Autumn'/>"
    "       <menuitem action='Bone'/>"
    "       <menuitem action='Colorcube'/>"
    "       <menuitem action='Cool'/>"
    "       <menuitem action='Copper'/>"
    "       <menuitem action='Flag'/>"
    "       <menuitem action='Gray'/>"
    "       <menuitem action='Hot'/>"
    "       <menuitem action='HSV'/>"
    "       <menuitem action='Jet'/>"
    "       <menuitem action='Lines'/>"
    "       <menuitem action='Pink'/>"
    "       <menuitem action='Prism'/>"
    "       <menuitem action='Spring'/>"
    "       <menuitem action='Summer'/>"
    "       <menuitem action='White'/>"
    "       <menuitem action='Winter'/>"
    "       <menuitem action='BlueWhiteRed'/>"
    "    </menu>"
    "  </menubar>"
    "  <toolbar  name='ToolBar'>"
    "    <toolitem action='Open'/>"
    "    <toolitem action='Quit'/>"
    "  </toolbar>"
    "</ui>";
  Manager = Gtk::UIManager::create();
  Manager->insert_action_group(Actions);
  add_accel_group (Manager->get_accel_group());
  Manager->add_ui_from_string (ui_info);

  ////////////////////
  // Setup choosers //
  ////////////////////
  SaveStateChooser.add_button (Gtk::Stock::CANCEL, Gtk::RESPONSE_CANCEL);
  SaveStateChooser.add_button (Gtk::Stock::SAVE,   Gtk::RESPONSE_OK);
  OpenStateChooser.add_button (Gtk::Stock::CANCEL, Gtk::RESPONSE_CANCEL);
  OpenStateChooser.add_button (Gtk::Stock::OPEN,   Gtk::RESPONSE_OK);

  /////////////////////
  // Connect signals //
  /////////////////////
  OrthoButton.signal_toggled().connect
    (sigc::mem_fun(*this, &ResVisualClass::OnPerspectiveToggle));
  ClipButton.signal_toggled().connect
    (sigc::mem_fun(*this, &ResVisualClass::OnClipToggle));
  IsoButton.signal_toggled().connect
    (sigc::mem_fun(*this, &ResVisualClass::OnIsoToggle));
  QuitButton.signal_clicked().connect
    (sigc::mem_fun(*this, &ResVisualClass::Quit));


  ////////////////////
  // Pack the boxes //
  ////////////////////
  ToolBox.pack_start(Tools);
  ToolBox.pack_start(FluidFrame, Gtk::PACK_SHRINK, 5);
  ToolBox.pack_start(PlaneFrame, Gtk::PACK_SHRINK, 5);
  ToolBox.pack_start(IsoFrame,  Gtk::PACK_SHRINK, 5);
  ToolBox.pack_start(StepFrame, Gtk::PACK_SHRINK, 5);
  MainVBox.pack_start(*Manager->get_widget("/MenuBar"), Gtk::PACK_SHRINK,0);
  MainVBox.pack_start(ToolBox, Gtk::PACK_SHRINK, 0);
  MiddleBox.pack_start(PathVis, Gtk::PACK_SHRINK, 5);
  MiddleBox.pack_start(OptionsBox, Gtk::PACK_SHRINK, 5);
  //  MainVBox.pack_start(PathVis);
  MainVBox.pack_start(MiddleBox, Gtk::PACK_SHRINK, 5);
  MainVBox.pack_start(QuitButton, Gtk::PACK_SHRINK, 0);

  add (MainVBox);
  set_title ("resvis++");
  show_all();

  PerspectButton.set_active(true);
  TruncRadiiToggle->set_active(false);
  IsocontourToggle->set_active(true);
  xPlane.SetIsocontours(true);
  yPlane.SetIsocontours(true);
  zPlane.SetIsocontours(true);

  UpdateIso = true;
  UpdatePlane[0] = UpdatePlane[1] = UpdatePlane[2] = true;
  CMapActions[BLUE_WHITE_RED]->set_active(true);
}

void
ResVisualClass::OnViewReset()
{
  PathVis.View.Reset();
  PathVis.Invalidate();
}


void
ResVisualClass::OnOpen()
{

}

void
ResVisualClass::Quit()
{
  Gtk::Main::quit();
}

void
ResVisualClass::OnExport()
{
  Export.SetupWidgets();
  Export.show_all();
}





string 
ResVisualClass::FindFullPath(string filename)
{
  string fullName;

  fullName = filename;
  if (Glib::file_test(fullName, Glib::FILE_TEST_EXISTS))
    return fullName;
  else {
    fullName = PKG_DATA_DIR + filename;
    if (Glib::file_test(fullName, Glib::FILE_TEST_EXISTS))
      return fullName;
    else {
      cerr << "Cannot find \"" << filename << "\" anywhere.\n";
      return filename;
    }
  }
}



void
ResVisualClass::OnClipToggle()
{
  Glib::signal_idle().connect
    (sigc::bind<bool>(mem_fun(*this, &ResVisualClass::DrawFrame), false));
}

void
ResVisualClass::OnIsoToggle()
{
  Glib::signal_idle().connect
    (sigc::bind<bool>(mem_fun(*this, &ResVisualClass::DrawFrame), false));
}

void
ResVisualClass::OnPerspectiveToggle()
{
  bool persp = !OrthoButton.get_active();
  PathVis.View.SetPerspective(persp);
  //  PathVis.Invalidate();
  DrawFrame();
  //  PathVis.Invalidate();
}


// void
// ResVisualClass::ClipSpheres(list<Vec3>& sphereList, double radius)
// {
//    list<Vec3>::iterator iter;
//    /// First, replicate spheres
//    for (iter=sphereList.begin(); iter != sphereList.end(); iter++) {
//      Vec3 &r = (*iter);
//      bool makeDisks = false;
//      if ((r[0]+radius) > 0.5*Box[0]) 
//        sphereList.push_front(Vec3(r[0]-Box[0], r[1], r[2]));
//       if ((r[0]-radius) < -0.5*Box[0]) 
// 	sphereList.push_front(Vec3(r[0]+Box[0], r[1], r[2]));
//     }
    
//     for (iter=sphereList.begin(); iter != sphereList.end(); iter++) {
//       Vec3 &r = (*iter);
//       if ((r[1]+radius) > 0.5*Box[1])
// 	sphereList.push_front(Vec3(r[0], r[1]-Box[1], r[2]));
//       if ((r[1]-radius) < -0.5*Box[1])
// 	sphereList.push_front(Vec3(r[0], r[1]+Box[1], r[2]));
//     }
    
//     for (iter=sphereList.begin(); iter != sphereList.end(); iter++) {
//       Vec3 &r = (*iter);
//       if ((r[2]+radius) > 0.5*Box[2])
// 	sphereList.push_front(Vec3(r[0], r[1], r[2]-Box[2]));
//       if ((r[2]-radius) < -0.5*Box[2])
// 	sphereList.push_front(Vec3(r[0], r[1], r[2]+Box[2]));
//     }
//     // Now make disks
//     for (iter=sphereList.begin(); iter != sphereList.end(); iter++) {
//       Vec3 &r = (*iter);
//       for (int dim=0; dim<3; dim++) {
// 	if ((r[dim]+radius) > Box[dim]) {
// 	  double l = 0.5*Box[dim]-fabs(r[dim]);
// 	  double diskRad = sqrt(radius*radius-l*l);
// 	  DiskObject *disk1 = new DiskObject(offScreen);
// 	  DiskObject *disk2 = new DiskObject(offScreen);
// 	  disk1->SetRadius(diskRad);
// 	  disk2->SetRadius(diskRad);
// 	  disk1->SetAxis(dim);
// 	  disk2->SetAxis(-dim);
// 	  Vec3 r1, r2;
// 	  r1 = r; r2 = r;
// 	  r1[dim] =  0.5*Box[dim];
// 	  r2[dim] = -0.5*Box[dim];
// 	  disk1->SetPos(r1);
// 	  disk2->SetPos(r2);
// 	  PathVis.Objects.push_back(disk1);
// 	  PathVis.Objects.push_back(disk2);
//       }
// }

class AtomClass
{
public:
  Vec3 Pos;
  int Type;
  AtomClass(Vec3 pos, int type) {
    Pos = pos;
    Type = type;
  }
};


void
ResVisualClass::SetStep (int step)
{
  StepScale.set_value(step);
}

bool
ResVisualClass::DrawFrame(bool offScreen)
{
  bool clipping = ClipButton.get_active();
  bool boxVisible = BoxToggle->get_active();
  for (int i=0; i<PathVis.Objects.size(); i++) 
    if (PathVis.Objects[i]->Dynamic) 
      delete PathVis.Objects[i];
  PathVis.Objects.resize(0);

  if (CoordToggle->get_active()) {
    CoordObject *coord = new CoordObject;
    Vec3 box = Box(0) + Box(1) + Box(2);
    coord->Set (box);
    PathVis.Objects.push_back(coord);
  }

  BoxObject *boxObject = new BoxObject;
  boxObject->SetColor (0.5, 0.5, 1.0);
  //boxObject->Set (Box, clipping);
  boxObject->Set (Box.GetLattice(), boxVisible, clipping);
  PathVis.Objects.push_back(boxObject);
  
  Vec3 nLattice[3];
  double length[3];
  length[0] = sqrt(dot(Box(0), Box(0)));
  length[1] = sqrt(dot(Box(1), Box(1)));
  length[2] = sqrt(dot(Box(2), Box(2)));
  nLattice[0] = Box(0)/length[0];
  nLattice[1] = Box(1)/length[1];
  nLattice[2] = Box(2)/length[2];

  Vec3 unitVecs[3], normVecs[3];
  for (int i=0; i<3; i++)
    unitVecs[i] = Box(i)/sqrt(dot(Box(i),Box(i)));
  normVecs[0] = cross (unitVecs[1], unitVecs[2]);
  normVecs[1] = cross (unitVecs[2], unitVecs[0]);
  normVecs[2] = cross (unitVecs[0], unitVecs[1]);
  for (int i=0; i<3; i++) {
    normVecs[i] /= sqrt(dot(normVecs[i], normVecs[i]));
    if (dot(normVecs[i], unitVecs[i]) < 0.0)
      normVecs[i] = -1.0*normVecs[i];
  }

  Array<double,3> &ResData = OilData;

  if (FileIsOpen) {
    int step;
    step = (int)round(StepScale.get_value());
    if ((CurrentStep != step)) {
      CurrentStep = step;
      ReadRes (step, CurrentFluid, DisplayProperty);
      Xgrid.Init(-0.5, 0.5, ResData.extent(0));
      Ygrid.Init(-0.5, 0.5, ResData.extent(1));
      Zgrid.Init(-0.5, 0.5, ResData.extent(2));
      ResIso.Init(&Xgrid, &Ygrid, &Zgrid, ResData, true);
      ResIso.SetLattice (Box.GetLattice());
      xPlane.SetCenter (uCenter, uMin, uMax);
      yPlane.SetCenter (uCenter, uMin, uMax);
      zPlane.SetCenter (uCenter, uMin, uMax);
      xPlane.Init(); yPlane.Init(); zPlane.Init();
    }
    if (ResetIso) {
      Xgrid.Init(-0.5, 0.5, ResData.extent(0));
      Ygrid.Init(-0.5, 0.5, ResData.extent(1));
      Zgrid.Init(-0.5, 0.5, ResData.extent(2));
      ResIso.Init(&Xgrid, &Ygrid, &Zgrid, ResData, true);
      ResIso.SetLattice (Box.GetLattice());
      ResIso.SetCenter (uCenter, uMin, uMax);
      ResetIso = false;
    }
    if (UpdateIso) {
      ResIso.SetColor (0.0, 0.8, 0.0);
      ResIso.SetCenter (uCenter, uMin, uMax);
      ResIso.SetIsoval(MinVal + (MaxVal-MinVal)*IsoAdjust.get_value());
      xPlane.SetCenter (uCenter, uMin, uMax);
      yPlane.SetCenter (uCenter, uMin, uMax);
      zPlane.SetCenter (uCenter, uMin, uMax);
      xPlane.Init(); yPlane.Init(); zPlane.Init();
    }

    if (UpdatePlane[0] && xPlaneButton.get_active()) 
      xPlane.SetPosition (0, xPlaneScale.get_value());
    if (xPlaneButton.get_active())
      PathVis.Objects.push_back(&xPlane);
    
    if (UpdatePlane[1] && yPlaneButton.get_active())
      yPlane.SetPosition (1, yPlaneScale.get_value());
    if (yPlaneButton.get_active())
      PathVis.Objects.push_back(&yPlane);
    
    if (UpdatePlane[2] && zPlaneButton.get_active())
      zPlane.SetPosition (2, zPlaneScale.get_value());
    if (zPlaneButton.get_active())
      PathVis.Objects.push_back(&zPlane);
    if (IsoButton.get_active()) 
      PathVis.Objects.push_back(&ResIso);
  }

  PathVis.Invalidate();
  UpToDate = true;
  UpdateIso = false;
  UpdatePlane[0] = false;
  UpdatePlane[1] = false;
  UpdatePlane[2] = false;
  return false;
}


bool 
IsDiag(Array<double,2> &lattice) 
{
  assert (lattice.extent(0) == 3);
  assert (lattice.extent(1) == 3);
  return ((fabs(lattice(0,1)) < 1.0e-16) &&
	  (fabs(lattice(1,0)) < 1.0e-16) &&
	  (fabs(lattice(0,2)) < 1.0e-16) &&
	  (fabs(lattice(2,0)) < 1.0e-16) &&
	  (fabs(lattice(1,2)) < 1.0e-16) &&
	  (fabs(lattice(2,1)) < 1.0e-16));
}

Mat3 ToMat3 (Array<double,2> &a)
{
  assert (a.rows()==3);
  assert (a.cols()==3);
  Mat3 b;
  b = 
    a(0,0), a(0,1), a(0,2), 
    a(1,0), a(1,1), a(1,2),
    a(2,0), a(2,1), a(2,2);
  return b;
}


void
ResVisualClass::Read(string filename)
{
  if (FileIsOpen)
    Infile.CloseFile();
  assert (Infile.OpenFile(filename));
  FileIsOpen = true;
//   string format;
//   Infile.ReadVar("format", format);

  /// Read box dimensions
  assert (Infile.OpenSection("geometry"));
  Vec3 box_dim;
  assert (Infile.ReadVar("box_dim", box_dim));
  Box.Set (box_dim);
  xPlane.SetLattice (Box.GetLattice());
  yPlane.SetLattice (Box.GetLattice());
  zPlane.SetLattice (Box.GetLattice());
  
  Infile.CloseSection(); // geometry


  /// Read first saturation
  Array<double,3> &ResData = OilData;
  NumSteps = Infile.CountSections("step");
  CurrentStep = 0; 
  CurrentFluid = 1;
  ReadRes(CurrentStep, CurrentFluid, SATURATION);
  Xgrid.Init(-0.5, 0.5, ResData.extent(0));
  Ygrid.Init(-0.5, 0.5, ResData.extent(1));
  Zgrid.Init(-0.5, 0.5, ResData.extent(2));
  ResIso.Init(&Xgrid, &Ygrid, &Zgrid, ResData, true);
  ResIso.SetLattice(Box.GetLattice());
  StepAdjust.set_upper(NumSteps-1.0);

  double maxDim = 1.2*sqrt(dot(box_dim, box_dim));
  PathVis.View.SetDistance (1.2*maxDim);
  
  IsoButton.set_active(false);
  IsoButton.set_sensitive(true);
  IsoFrame.set_sensitive(true);
  DrawFrame();
}


	
void
ResVisualClass::OnIsoChange()
{
  double rho = IsoAdjust.get_value() * MaxVal;
  double rs = pow (3.0/(4.0*M_PI*rho),1.0/3.0);
  char rstext[100];
  snprintf (rstext, 100, "rs = %1.3f", rs);
  
  UpdateIso = true;
  UpdateIsoVal = true;
  DrawFrame();
}

void
ResVisualClass::OnStepChange()
{
  UpdateIso = true;
  UpdatePlane[0] = true;
  UpdatePlane[1] = true;
  UpdatePlane[2] = true;
  Glib::signal_idle().connect
    (sigc::bind<bool>(mem_fun(*this, &ResVisualClass::DrawFrame), false));
}

void
ResVisualClass::OnPlaneChange(int dir)
{
  if (dir==0 && xPlaneButton.get_active())
    UpdatePlane[0] = true;
  if (dir==1 && yPlaneButton.get_active())
    UpdatePlane[1] = true;
  if (dir==2 && zPlaneButton.get_active())
    UpdatePlane[2] = true;
  DrawFrame();
}

void
ResVisualClass::OnFluidChange(FluidType fl)
{
  if (fl==GAS && GasButton.get_active())
    UpdateGas = true;
  if (fl==OIL && OilButton.get_active())
    UpdateOil = true;
  if (fl==WATER && WaterButton.get_active())
    UpdateWater = true;
  DrawFrame();
}


ResVisualClass::~ResVisualClass()
{

}


bool
ResVisualClass::ReadRes (int step, int fluid, PropertyType prop)
{
  Array<double,3> &ResData = OilData;
  assert (Infile.OpenSection("step", step));
  assert (Infile.OpenSection("fluid", fluid));
  string typeString;
  if (prop == PRESSURE)
    typeString = "pressure";
  else if (prop == SATURATION)
    typeString = "saturation";
  else {
    cerr << "Unknown property " << prop << endl;
    abort();
  }
//   string typeString("porosity");
//   assert (Infile.OpenSection("geometry"));
  assert (Infile.ReadVar(typeString, ResData));
//   Infile.CloseSection(); // "geometry"
  Infile.CloseSection(); // "fluid"
  Infile.CloseSection(); // "step"
  uMin = Vec3 (0.0, 0.0, 0.0);
  uMax = Vec3 (1.0, 1.0, 1.0);

  Mat3 l = Box.GetLattice();

  MaxVal = -1.0e50;
  MinVal = 1.0e50;
  int xShift, yShift, zShift;
  Vec3 u(0.0, 0.0, 0.0);
  Vec3 du(1.0/(double)(ResData.extent(0)-1),
	  1.0/(double)(ResData.extent(1)-1),
	  1.0/(double)(ResData.extent(2)-1));
  
  for (int ix=0; ix<ResData.extent(0); ix++) 
    for (int iy=0; iy<ResData.extent(1); iy++) 
      for (int iz=0; iz<ResData.extent(2); iz++) {
	MaxVal = max(MaxVal, ResData(ix,iy,iz));	
	MinVal = min(MinVal, ResData(ix,iy,iz));
      }
  return true;
}



void
ResVisualClass::OnCoordToggle()
{
  DrawFrame();
}

void
ResVisualClass::OnSphereToggle()
{
  DrawFrame();
}


void
ResVisualClass::OnBondsToggle()
{
  DrawFrame();
}


void
ResVisualClass::OnBoxToggle()
{
  DrawFrame();
}

void
ResVisualClass::OnTruncRadiiToggle()
{
  DrawFrame();
}

void
ResVisualClass::OnIsocontourToggle()
{
  bool active = IsocontourToggle->get_active();
  xPlane.SetIsocontours (active);
  yPlane.SetIsocontours (active);
  zPlane.SetIsocontours (active);
  UpdatePlane[0] = xPlaneButton.get_active();
  UpdatePlane[1] = yPlaneButton.get_active();
  UpdatePlane[2] = zPlaneButton.get_active();
  DrawFrame();
}

void
ResVisualClass::OnRadiusChange()
{
  DrawFrame();
}

void
ResVisualClass::OnDisplayRadio(PropertyType type)
{
  PropertyType newtype;
  if (SaturationRadio->get_active() && type==SATURATION)
    newtype = SATURATION;
  if (PressureRadio->get_active() && type==PRESSURE)
    newtype = PRESSURE;

  if (DisplayProperty != newtype) {
    DisplayProperty = newtype;
    ReadRes (CurrentStep, CurrentFluid, newtype);
    UpdatePlane[0] = UpdatePlane[1] = UpdatePlane[2] = true;
    UpdateIso = true;  ResetIso = true;
    DrawFrame ();
  }
  UpdateIsoType = true;
}


bool
ResVisualClass::WriteState(string fname)
{
  IOSectionClass out;
  
  if (out.NewFile(fname) == false)
    return false;

  out.NewSection ("Flags");
  out.WriteVar ("xPlane", xPlaneButton.get_active());
  out.WriteVar ("yPlane", yPlaneButton.get_active());
  out.WriteVar ("zPlane", zPlaneButton.get_active());
  out.WriteVar ("Nuclei", SphereToggle->get_active());
  out.WriteVar ("CoordAxes", CoordToggle->get_active());
  out.WriteVar ("Isosurface", IsoButton.get_active());
  out.WriteVar ("Clip", ClipButton.get_active());
  out.WriteVar ("Perspective", PerspectButton.get_active());

  string restype;
  if (DisplayProperty == SATURATION)
    restype = "Saturation";
  if (DisplayProperty == PRESSURE)
    restype = "Pressure";
  out.WriteVar ("DisplayProperty", restype);
  out.CloseSection();

  out.WriteVar ("Step", (int)round(StepAdjust.get_value()));
  out.WriteVar ("IsoVal", IsoAdjust.get_value());
  out.WriteVar ("xPlanePos", xPlaneAdjust.get_value());
  out.WriteVar ("yPlanePos", yPlaneAdjust.get_value());
  out.WriteVar ("zPlanePos", zPlaneAdjust.get_value());
  out.WriteVar ("Radius", RadiusScale.get_value());
  out.NewSection("View");
  PathVis.View.Write(out);
  out.CloseSection();
  out.CloseFile();
  return true;
}

bool
ResVisualClass::ReadState (string fname)
{
  IOSectionClass in;
  if (in.OpenFile(fname) == false)
    return false;
  bool active;
  
  assert(in.OpenSection ("Flags"));
  in.ReadVar ("xPlane", active);
  xPlaneButton.set_active(active);
  assert(in.ReadVar ("yPlane", active));
  yPlaneButton.set_active(active);
  assert(in.ReadVar ("zPlane", active));
  zPlaneButton.set_active(active);
  assert(in.ReadVar ("Nuclei", active));
  SphereToggle->set_active(active);
  assert(in.ReadVar ("Isosurface", active));
  IsoButton.set_active(active);
  assert(in.ReadVar ("Clip", active));
  ClipButton.set_active(active);
  assert(in.ReadVar ("Perspective", active));
  PerspectButton.set_active(active);
  if (in.ReadVar ("CoordAxes", active))
    CoordToggle->set_active(active);

  string restype;
  assert(in.ReadVar ("DisplayProperty", restype));
  if (restype == "SATURATION") {
    SaturationRadio->set_active(true);
    DisplayProperty = SATURATION;
  }
  else if (restype == "PRESSURE") {
    PressureRadio->set_active(true);
    DisplayProperty = PRESSURE;
  }
  in.CloseSection(); // flags

  int step;
  assert(in.ReadVar ("Step", step));
  double val;
  assert(in.ReadVar ("IsoVal", val));
  IsoAdjust.set_value(val);

  assert(in.ReadVar ("xPlanePos", val));
  xPlaneAdjust.set_value(val);
  assert(in.ReadVar ("yPlanePos", val));
  yPlaneAdjust.set_value(val);
  assert(in.ReadVar ("zPlanePos", val));
  zPlaneAdjust.set_value(val);
  if (in.ReadVar("Radius", val))
    RadiusScale.set_value(val);

  StepAdjust.set_value((double)step);
  assert(in.OpenSection("View"));
  PathVis.View.Read(in);
  in.CloseSection();
  in.CloseFile();
  Glib::signal_idle().connect
    (sigc::bind<bool>(mem_fun(*this, &ResVisualClass::DrawFrame), false));
  return true;
}

void
ResVisualClass::OnSaveState()
{
  int result = SaveStateChooser.run();
  if (result == Gtk::RESPONSE_OK) {
    string fname = SaveStateChooser.get_filename();
    WriteState (fname);
  }
  SaveStateChooser.hide();
}

void
ResVisualClass::OnOpenState()
{
  int result = OpenStateChooser.run();
  if (result == Gtk::RESPONSE_OK) {
    string fname = OpenStateChooser.get_filename();
    ReadState (fname);
  }
  OpenStateChooser.hide();
}  


void
ResVisualClass::SetShift(Vec3 shift)
{
  DoShift = true;
  Shift = shift;
}


void
ResVisualClass::SetViewportSize (int size)
{
  PathVis.set_size_request(size, size);
  resize(10,10);
}


void
ResVisualClass::OnFullscreenToggle()
{
  if (FullscreenToggle->get_active())
    fullscreen();
  else
    unfullscreen();
}


void
ResVisualClass::OnColorMapRadio (ColorMapType type)
{
  if (CMapActions[type]->get_active()) {
    CMap = type;
    xPlane.SetColorMap(type);
    yPlane.SetColorMap(type);
    zPlane.SetColorMap(type);
    if (xPlaneButton.get_active())    xPlane.Set();
    if (yPlaneButton.get_active())    yPlane.Set();
    if (zPlaneButton.get_active())    zPlane.Set();
    if (xPlaneButton.get_active() || yPlaneButton.get_active() || zPlaneButton.get_active())
      DrawFrame();
  }
}


vector<string>
BreakString (string str, char sep)
{
  vector<string> strvec;
  int len = str.size();
  int i = 0;
  while (i < len) {
    char s[2];
    s[1] = '\0';
    string item;
    while (str[i] != sep && i < len) {
      s[0] = str[i];
      item.append(s);
      i++;
    }
    strvec.push_back(item);
    i++;
  }
  return strvec;
}


int main(int argc, char** argv)
{
  Gtk::Main kit(argc, argv);

  // Init gtkglextmm.
  Gtk::GL::init(argc, argv);
  glutInit(&argc, argv);

  list<ParamClass> optionList;
  optionList.push_back(ParamClass("shift", true));
  optionList.push_back(ParamClass("small", false));
  optionList.push_back(ParamClass("remote", false));
  optionList.push_back(ParamClass("walkers", true));
  CommandLineParserClass parser (optionList);
  bool success = parser.Parse (argc, argv);
  if (!success || parser.NumFiles() < 1 || parser.NumFiles() > 2) {
    cerr << "Usage:\n  resvis++ [options...] myfile.h5 [statefile.h5]\n"
	 << "Options:\n"
	 << "  --shift x,y,z          shift by reduced coordinates\n"
	 << "  --small                reduce size for small displays\n"
	 << "  --remote               reduce data transfer for remote operation\n"
	 << "  --walkers f.config.h5  read electron positions from a QMCPACK config dump\n"
	 << "                      over a network connection\n";
    exit (1);
  }
  
  // Query OpenGL extension version.
  int major, minor;
  Gdk::GL::query_version(major, minor);
  std::cout << "OpenGL extension version - "
	    << major << "." << minor << std::endl;

  // Instantiate and run the application.
  ResVisualClass resvisual;

  if (parser.Found("shift")) {
    string shiftStr = parser.GetArg("shift");
    Vec3 shift;
    vector<string> components = BreakString (shiftStr, ',');
    if (components.size() != 3) {
      cerr << "Expected 3 components for shift.\n";
      abort();
    }
    for (int i=0; i<3; i++)
      shift[i] = strtod (components[i].c_str(), NULL);
    resvisual.SetShift (shift);
  }
  if (parser.Found("small"))
    resvisual.SetViewportSize (600);

  resvisual.Read (parser.GetFile(0));

  if (parser.NumFiles() == 2)
    resvisual.ReadState (parser.GetFile(1));
  kit.run(resvisual);

  return 0;
}
