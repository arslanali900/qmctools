#include "WFVis.h"
#include <GL/glut.h>
#include "ElementData.h"
#include "ParseCommand.h"
#include "PointsObject.h"



WFVisualClass::WFVisualClass() :
  MainVBox(false, 0), 
  QuitButton("Quit"),
  IsoAdjust  (0.01, 0.0, 1.0, 0.01, 0.1),
  BandAdjust (0.0, 0.0, 8.0, 1.0, 1.0),
  kAdjust (0.0, 0.0, 16.0, 1.0, 1.0),
  xPlaneAdjust(0.0, 0.0, 1.0, 0.01, 0.1),
  yPlaneAdjust(0.0, 0.0, 1.0, 0.01, 0.1),
  zPlaneAdjust(0.0, 0.0, 1.0, 0.01, 0.1),
  RadiusAdjust(0.4, 0.0, 1.0, 0.01, 0.1),
  UpToDate(true),
  FileIsOpen(false),
  xPlane(WFIso),
  yPlane(WFIso),
  zPlane(WFIso),
  Export(*this),
  WFDisplay(MAG2),
  ResetIso(false),
  SaveStateChooser ("State filename", Gtk::FILE_CHOOSER_ACTION_SAVE),
  OpenStateChooser ("State filename", Gtk::FILE_CHOOSER_ACTION_OPEN),
  UpdateIsoType(false), 
  UpdateIsoVal(false),
  DoShift (false),
  Shift (0.0, 0.0, 0.0),
  CMap(BLUE_WHITE_RED),
  Nonuniform(false)
{
  if (!Glib::thread_supported())
    Glib::thread_init();
  WFIso.Dynamic = false;
  xPlane.Dynamic = false;
  yPlane.Dynamic = false;
  zPlane.Dynamic = false;
  IsoScale.set_adjustment(IsoAdjust);
  IsoScale.set_digits(2);
  IsoScale.signal_value_changed().connect
    (sigc::mem_fun(*this, &WFVisualClass::OnIsoChange));
  IsoFrame.set_label("Density");
  IsoFrame.add(DensityBox);
  DensityBox.pack_start(rsLabel);
  DensityBox.pack_start(IsoScale);
  rsLabel.set_text("rs = ");
  IsoScale.property_width_request().set_value(75);

  BandScale.set_adjustment (BandAdjust);
  BandScale.set_digits(0);
  BandScale.signal_value_changed().connect
    (sigc::mem_fun(*this, &WFVisualClass::OnBandChange));
  BandScale.property_width_request().set_value(75);
  BandFrame.set_label("Band");
  BandFrame.add(BandScale);

  kScale.set_adjustment (kAdjust);
  kScale.set_digits(0);
  kScale.signal_value_changed().connect
    (sigc::mem_fun(*this, &WFVisualClass::OnkChange));
  kScale.property_width_request().set_value(75);
  kFrame.set_label("k point");
  kFrame.add(kScale);

  RadiusScale.set_adjustment(RadiusAdjust);
  RadiusScale.set_digits(2);
  RadiusScale.signal_value_changed().connect
    (sigc::mem_fun(*this, &WFVisualClass::OnRadiusChange));
  RadiusScale.property_width_request().set_value(75);
  RadiusFrame.set_label("Ion radii");
  RadiusFrame.add(RadiusScale);
  RadiusBox.pack_start(RadiusFrame, Gtk::PACK_SHRINK, 5);
  OptionsBox.pack_start(RadiusBox,  Gtk::PACK_SHRINK, 5);

  ///////////////////////////////////////
  // Setup multiband selection widgets //
  ///////////////////////////////////////
  VisibleBandFrame.add (VisibleBandBox);
  VisibleBandBox.pack_start (MultiBandButton, Gtk::PACK_SHRINK, 5);
  VisibleBandBox.pack_start (VisibleBandWindow);
  VisibleBandWindow.add (VisibleBandTable);
  VisibleBandFrame.set_label("Multiband");
  OptionsBox.pack_start(VisibleBandFrame, Gtk::PACK_SHRINK, 5);
  VisibleBandWindow.property_height_request().set_value(680);
  VisibleBandWindow.set_policy (Gtk::POLICY_NEVER, Gtk::POLICY_AUTOMATIC);
  MultiBandButton.set_label ("Enable");
  MultiBandButton.signal_toggled().connect
    (sigc::mem_fun(*this, &WFVisualClass::OnMultiBandToggle));

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
    (sigc::bind<int>(sigc::mem_fun(*this, &WFVisualClass::OnPlaneChange), 0));
  yPlaneScale.signal_value_changed().connect
    (sigc::bind<int>(sigc::mem_fun(*this, &WFVisualClass::OnPlaneChange), 1));
  zPlaneScale.signal_value_changed().connect
    (sigc::bind<int>(sigc::mem_fun(*this, &WFVisualClass::OnPlaneChange), 2));
  xPlaneButton.signal_toggled().connect
    (sigc::bind<int>(sigc::mem_fun(*this, &WFVisualClass::OnPlaneChange), 0));
  yPlaneButton.signal_toggled().connect
    (sigc::bind<int>(sigc::mem_fun(*this, &WFVisualClass::OnPlaneChange), 1));
  zPlaneButton.signal_toggled().connect
    (sigc::bind<int>(sigc::mem_fun(*this, &WFVisualClass::OnPlaneChange), 2));


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
		sigc::mem_fun(*this, &WFVisualClass::OnOpen));
  Actions->add (Gtk::Action::create("Export", "_Export Image"),
		sigc::mem_fun(*this, &WFVisualClass::OnExport));
  Actions->add (Gtk::Action::create("OpenState", "Open State"),
		sigc::mem_fun(*this, &WFVisualClass::OnOpenState));
  Actions->add (Gtk::Action::create("SaveState", "Save State"),
		sigc::mem_fun(*this, &WFVisualClass::OnSaveState));
  Actions->add (Gtk::Action::create("Quit", "_Quit"),
		sigc::mem_fun(*this, &WFVisualClass::Quit));
  Actions->add (Gtk::Action::create("MenuView", "View"));
  Actions->add (Gtk::Action::create("Reset", "Reset"),
		sigc::mem_fun(*this, &WFVisualClass::OnViewReset));
  CoordToggle = Gtk::ToggleAction::create("Axes", "Coordinate axes",
					  "Show coordinate axes", true);
  Actions->add (CoordToggle,
		sigc::mem_fun(*this, &WFVisualClass::OnCoordToggle));
  SphereToggle = Gtk::ToggleAction::create("Nuclei", "Nuclei",
					  "Show nuclei", true);
  Actions->add (SphereToggle,
		sigc::mem_fun(*this, &WFVisualClass::OnSphereToggle));

  BondsToggle = Gtk::ToggleAction::create("Bonds", "Bonds",
					  "Show bonds", false);
  Actions->add (BondsToggle,
		sigc::mem_fun(*this, &WFVisualClass::OnBondsToggle));

  BoxToggle = Gtk::ToggleAction::create("Box", "Box", "Show box", true);
  Actions->add (BoxToggle,
		sigc::mem_fun(*this, &WFVisualClass::OnBoxToggle));

  TruncRadiiToggle = 
    Gtk::ToggleAction::create("Trunc", "Truncation radii",
			      "Show truncation radii", true);
  Actions->add (TruncRadiiToggle,
		sigc::mem_fun(*this, &WFVisualClass::OnTruncRadiiToggle));

  IsocontourToggle = 
    Gtk::ToggleAction::create("Isocontours", "Plane isocontours",
			      "Show isocontours on color planes", true);
  Actions->add (IsocontourToggle,
		sigc::mem_fun(*this, &WFVisualClass::OnIsocontourToggle));
  FullscreenToggle = Gtk::ToggleAction::create("Fullscreen", "Fullscreen mode",
					       "Display this window in fullscreen mode", false);
  Actions->add (FullscreenToggle,
		sigc::mem_fun(*this, &WFVisualClass::OnFullscreenToggle));

  Actions->add (Gtk::Action::create("MenuDisplay", "Display"));
  Mag2Radio = Gtk::RadioAction::create
    (DisplayGroup, "Mag2", "Magnitude squared",
     "Display square magnitude of WF");
  RealRadio = Gtk::RadioAction::create
    (DisplayGroup, "Real", "Real part", "Display real part of WF");
  ImagRadio = Gtk::RadioAction::create
    (DisplayGroup, "Imag", "Imag part", "Display imaginary part of WF");
  Actions->add
    (RealRadio, sigc::bind<WFDisplayType> 
     (sigc::mem_fun(*this, &WFVisualClass::OnDisplayRadio),REAL_PART));
  Actions->add
    (ImagRadio, sigc::bind<WFDisplayType> 
     (sigc::mem_fun(*this, &WFVisualClass::OnDisplayRadio),IMAG_PART));
  Actions->add
    (Mag2Radio, sigc::bind<WFDisplayType> 
     (sigc::mem_fun(*this, &WFVisualClass::OnDisplayRadio),MAG2));  


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
       (sigc::mem_fun (*this, &WFVisualClass::OnColorMapRadio),
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
    "      <menuitem action='Real'/>"
    "      <menuitem action='Imag'/>"
    "      <menuitem action='Mag2'/>"
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
    (sigc::mem_fun(*this, &WFVisualClass::OnPerspectiveToggle));
  ClipButton.signal_toggled().connect
    (sigc::mem_fun(*this, &WFVisualClass::OnClipToggle));
  IsoButton.signal_toggled().connect
    (sigc::mem_fun(*this, &WFVisualClass::OnIsoToggle));
  QuitButton.signal_clicked().connect
    (sigc::mem_fun(*this, &WFVisualClass::Quit));


  ////////////////////
  // Pack the boxes //
  ////////////////////
  ToolBox.pack_start(Tools);
  ToolBox.pack_start(PlaneFrame, Gtk::PACK_SHRINK, 5);
  ToolBox.pack_start(IsoFrame,  Gtk::PACK_SHRINK, 5);
  ToolBox.pack_start(BandFrame, Gtk::PACK_SHRINK, 5);
  ToolBox.pack_start(kFrame,    Gtk::PACK_SHRINK, 5);
  MainVBox.pack_start(*Manager->get_widget("/MenuBar"), Gtk::PACK_SHRINK,0);
  MainVBox.pack_start(ToolBox, Gtk::PACK_SHRINK, 0);
  MiddleBox.pack_start(PathVis, Gtk::PACK_SHRINK, 5);
  MiddleBox.pack_start(OptionsBox, Gtk::PACK_SHRINK, 5);
  //  MainVBox.pack_start(PathVis);
  MainVBox.pack_start(MiddleBox, Gtk::PACK_SHRINK, 5);
  MainVBox.pack_start(QuitButton, Gtk::PACK_SHRINK, 0);

  add (MainVBox);
  set_title ("wfvis++");
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
WFVisualClass::OnViewReset()
{
  PathVis.View.Reset();
  PathVis.Invalidate();
}


void
WFVisualClass::OnOpen()
{

}

void
WFVisualClass::Quit()
{
  Gtk::Main::quit();
}

void
WFVisualClass::OnExport()
{
  Export.SetupWidgets();
  Export.show_all();
}





string 
WFVisualClass::FindFullPath(string filename)
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
WFVisualClass::OnClipToggle()
{
  Glib::signal_idle().connect
    (sigc::bind<bool>(mem_fun(*this, &WFVisualClass::DrawFrame), false));
}

void
WFVisualClass::OnIsoToggle()
{
  Glib::signal_idle().connect
    (sigc::bind<bool>(mem_fun(*this, &WFVisualClass::DrawFrame), false));
}

void
WFVisualClass::OnPerspectiveToggle()
{
  bool persp = !OrthoButton.get_active();
  PathVis.View.SetPerspective(persp);
  //  PathVis.Invalidate();
  DrawFrame();
  //  PathVis.Invalidate();
}


// void
// WFVisualClass::ClipSpheres(list<Vec3>& sphereList, double radius)
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


bool
WFVisualClass::DrawFrame(bool offScreen)
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

  if (Walkers.size()) {
    PointsObject *points = new PointsObject();
    Array<Vec3,2> p(Walkers.shape());
    Mat3 Linv = Box.GetLatticeInv();
    Mat3 L = Box.GetLattice();
    for (int i=0; i<Walkers.extent(0); i++)
      for (int j=0; j<Walkers.extent(1); j++) {
	Vec3 u = Linv * Walkers(i,j);
	u[0] -= floor(u[0]);
	u[1] -= floor(u[1]);
	u[2] -= floor(u[2]);
	u -= Vec3(0.5,0.5,0.5);
	p(i,j) = L * u;
      }
    points->SetPoints (p);
    PathVis.Objects.push_back(points);
  }


  list<AtomClass> sphereList;
  if (SphereToggle->get_active()) {
    list<AtomClass>::iterator iter;
    for (int ptcl=0; ptcl < AtomPos.extent(0); ptcl++) 
      sphereList.push_back(AtomClass(AtomPos(ptcl), AtomTypes(ptcl)));
    if (clipping) {
      for (int dim=0; dim<3; dim++) {
	for (iter=sphereList.begin(); iter != sphereList.end(); iter++) {
	  int type = (*iter).Type;
	  double radius = RadiusScale.get_value() *
	    ElementData::GetRadius(type);
	  
	  Vec3 r = (*iter).Pos;
	  Vec3 n = Box.GetLatticeInv() * r; 
	  r = Box.GetLattice()*n;

	  Vec3 rplus = r + radius*normVecs[dim];
	  Vec3 nplus = Box.GetLatticeInv() * rplus;
	  Vec3 rminus = r - radius*normVecs[dim];
	  Vec3 nminus = Box.GetLatticeInv() * rminus;

	  if (nplus[dim] < nminus[dim])
	    swap(nplus, nminus);

	  if ((nplus[dim] > 0.5) || nminus[dim] > 0.5) {
	    Vec3 nNew = n;
	    nNew[dim] -= 1.0;
	    Vec3 rNew = Box.GetLattice()*nNew;
	    sphereList.push_front(AtomClass(rNew, type));
	  }
	  if (nminus[dim] < -0.5 || nplus[dim] < -0.5) {
	    Vec3 nNew = n;
	    nNew[dim] += 1.0;
	    Vec3 rNew = Box.GetLattice()*nNew;
	    sphereList.push_front(AtomClass(rNew, type));
	  }
	}
      }
//       for (iter=sphereList.begin(); iter != sphereList.end(); iter++) {
// 	Vec3 &r = (*iter).Pos;
// 	Vec3 n = Box.GetLatticeInv() * r; 
// 	int type = (*iter).Type;
// 	double radius = RadiusScale.get_value() *
// 	  ElementData::GetRadius(type)/length[0];
// 	if ((n[1]+radius/length[1]) > 0.5)
// 	  sphereList.push_front(AtomClass(r-Box(1),type));
// 	if ((n[1]-radius/length[1]) < -0.5)
// 	  sphereList.push_front(AtomClass(r+Box(1),type));
//       }
      
//       for (iter=sphereList.begin(); iter != sphereList.end(); iter++) {
// 	Vec3 &r = (*iter).Pos;
// 	Vec3 n = Box.GetLatticeInv() * r; 
// 	int type = (*iter).Type;
// 	double radius = RadiusScale.get_value() *
// 	  ElementData::GetRadius(type);
// 	if ((n[2]+radius/length[2]) > 0.5)
// 	  sphereList.push_front(AtomClass(r-Box(2),type));
// 	if ((n[2]-radius/length[2]) < -0.5)
// 	  sphereList.push_front(AtomClass(r+Box(2),type));
//       }
      // Now make disks to close spheres
      for (iter=sphereList.begin(); iter != sphereList.end(); iter++) {
	Vec3 &r = (*iter).Pos;
	int type = (*iter).Type;
	double radius = RadiusScale.get_value() *
	  ElementData::GetRadius(type);
	
	Mat3 lattice = Box.GetLattice();
	Vec3 unitVecs[3], normVecs[3];
	for (int i=0; i<3; i++)
	  unitVecs[i] = Box(i)/sqrt(dot(Box(i),Box(i)));
	normVecs[0] = cross (unitVecs[1], unitVecs[2]);
	normVecs[1] = cross (unitVecs[2], unitVecs[0]);
	normVecs[2] = cross (unitVecs[0], unitVecs[1]);
	for (int i=0; i<3; i++)
	  normVecs[i] = normVecs[i]/sqrt(dot(normVecs[i], normVecs[i]));

	Vec3 n = Box.GetLatticeInv() * r;
	for (int dim=0; dim<3; dim++) {
	  Vec3 rplus = r + radius*normVecs[dim];
	  Vec3 nplus = Box.GetLatticeInv() * rplus;
	  Vec3 rminus = r - radius*normVecs[dim];
	  Vec3 nminus = Box.GetLatticeInv() * rminus;
	  if (fabs(nplus[dim]) > 0.5) {
	    double nsign = nplus[dim]/fabs(nplus[dim]);
	    double dr = dot(normVecs[dim], (0.5*nsign -n[dim])*Box(dim));
	    Vec3 delta = 0.9999*dr*normVecs[dim];
	    Vec3 c1, c2;
	    c1 = r + delta;
	    c2 = c1 - nsign*Box(dim);
	    double diskRad = sqrt(radius*radius - dr*dr);
	    DiskObject *disk1 = new DiskObject(offScreen);
	    DiskObject *disk2 = new DiskObject(offScreen);
	    disk1->SetRadius(diskRad);
	    disk2->SetRadius(diskRad);
	    disk1->SetNormVec (+1.0*normVecs[dim]);
	    disk2->SetNormVec (-1.0*normVecs[dim]);
	    Vec3 color = ElementData::GetColor (type);
 	    disk1->SetColor(color);
 	    disk2->SetColor(color);
	    disk1->SetPos(c1);
	    disk2->SetPos(c2);
	    PathVis.Objects.push_back(disk1);
	    PathVis.Objects.push_back(disk2);
	  }

// 	  if ((r[dim]+radius) > 0.5*Box[dim]) {
// 	    double l = 0.5*Box[dim]-fabs(r[dim]);
// 	    double diskRad = sqrt(radius*radius-l*l);
// 	    DiskObject *disk1 = new DiskObject(offScreen);
// 	    DiskObject *disk2 = new DiskObject(offScreen);
// 	    Vec3 color = ElementData::GetColor (type);
// 	    disk1->SetRadius(diskRad);
// 	    disk2->SetRadius(diskRad);
// 	    disk1->SetAxis(2*dim);
// 	    disk2->SetAxis(2*dim+1);
// 	    disk1->SetColor(color);
// 	    disk2->SetColor(color);
// 	    Vec3 r1, r2;
// 	    r1 = r; r2 = r;
// 	    r1[dim] =  0.5*Box[dim];
// 	    r2[dim] = -0.5*Box[dim];
// 	    disk1->SetPos(r1);
// 	    disk2->SetPos(r2);
// 	    PathVis.Objects.push_back(disk1);
// 	    PathVis.Objects.push_back(disk2);
// 	  }
	}
      }
    }
    /// Add sphere objects from list
    for (iter=sphereList.begin(); iter!=sphereList.end(); iter++) {
      SphereObject *sphere = new SphereObject (offScreen);
      sphere->SetPos((*iter).Pos);
      Vec3 color = ElementData::GetColor (iter->Type);
      double radius = RadiusScale.get_value() *
	ElementData::GetRadius(iter->Type);
      sphere->SetColor(color);
      sphere->SetBox(Box.GetLattice());
      sphere->SetRadius(radius);
      PathVis.Objects.push_back(sphere);
    }
  }

  if (BondsToggle->get_active()) {
    list<AtomClass>::iterator iter1, iter2;
    //    const double bondLength = 2.0;
    for (iter1=sphereList.begin(); iter1!=sphereList.end(); iter1++) {
      // for (int i=0; i < AtomPos.extent(0); i++) {
      // Vec3 ri = AtomPos(i);
      // double rad1 = ElementData::GetRadius(AtomTypes(i));
      Vec3 ri = iter1->Pos;
      double rad1 = ElementData::GetRadius(iter1->Type);
      // for (int j=i+1; j < AtomPos.extent(0); j++) {
      // 	Vec3 rj = AtomPos(j);
      // 	double rad2 = ElementData::GetRadius(AtomTypes(j));
      iter2 = iter1;
      iter2++;
      for (; iter2!=sphereList.end(); iter2++) {
	Vec3 rj = iter2->Pos;
	double rad2 = ElementData::GetRadius(iter2->Type);
	double bondLength = 0.9 * (rad1+rad2);
	if (dot(ri-rj, ri-rj) < bondLength*bondLength) {
	  CylinderObject* bond = new CylinderObject(offScreen);
	  bond->SetPos (ri, rj);
	  bond->SetRadius(0.2);
	  bond->SetColor (Vec3 (0.5, 0.5, 0.5));
	  bond->SetBox (Box.GetLattice());
	  PathVis.Objects.push_back(bond);
	}
      }
    }
  }

  if (FileIsOpen && !MultiBandButton.get_active()) {
    int band, k;
    band = (int)round(BandScale.get_value());
    k    = (int)round(kScale.get_value());
    if ((CurrBand != band) || (Currk != k)) {
      CurrBand = band;
      Currk    = k;
      if (IsESHDF) 
	ReadWF_ESHDF (k, band);
      else
	ReadWF (k, band);
      if (Nonuniform)
	WFIso.Init(&NUXgrid, &NUYgrid, &NUZgrid, WFData, true);
      else {
	Xgrid.Init(-0.5, 0.5, WFData.extent(0));
	Ygrid.Init(-0.5, 0.5, WFData.extent(1));
	Zgrid.Init(-0.5, 0.5, WFData.extent(2));
	WFIso.Init(&Xgrid, &Ygrid, &Zgrid, WFData, true);
      }
      // WFIso.Init (-0.5, 0.5, -0.5, 0.5, -0.5, 0.5, WFData);
      WFIso.SetLattice (Box.GetLattice());
      xPlane.SetCenter (uCenter, uMin, uMax);
      yPlane.SetCenter (uCenter, uMin, uMax);
      zPlane.SetCenter (uCenter, uMin, uMax);
      xPlane.Init(); yPlane.Init(); zPlane.Init();
    }
    if (ResetIso) {
      if (Nonuniform) 
	WFIso.Init(&NUXgrid, &NUYgrid, &NUZgrid, WFData, true);
      else {
	Xgrid.Init(-0.5, 0.5, WFData.extent(0));
	Ygrid.Init(-0.5, 0.5, WFData.extent(1));
	Zgrid.Init(-0.5, 0.5, WFData.extent(2));
	WFIso.Init(&Xgrid, &Ygrid, &Zgrid, WFData, true);
	// WFIso.Init (-0.5, 0.5, -0.5, 0.5, -0.5, 0.5, WFData);
      }
      WFIso.SetLattice (Box.GetLattice());
      WFIso.SetCenter (uCenter, uMin, uMax);
      ResetIso = false;
    }
    if (UpdateIso) {
      if (WFDisplay == MAG2) {
	WFIso.SetColor (0.0, 0.8, 0.0);
	WFIso.SetCenter (uCenter, uMin, uMax);
	WFIso.SetIsoval(MaxVal*IsoAdjust.get_value());
      }
      else {
	vector<TinyVector<double,3> > colors;
	colors.push_back(TinyVector<double,3> (0.8, 0.0, 0.0));
	colors.push_back(TinyVector<double,3> (0.0, 0.0, 0.8));
	WFIso.SetColor(colors);
	WFIso.SetCenter (uCenter, uMin, uMax);
	vector<double> vals;
	vals.push_back(+MaxVal*sqrt(IsoAdjust.get_value()));
	vals.push_back(-MaxVal*sqrt(IsoAdjust.get_value()));
	WFIso.SetIsoval(vals);
      }
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
      PathVis.Objects.push_back(&WFIso);

    // If localized, add a translucent sphere to signify the
    // truncation radius
    if (Localized && TruncRadiiToggle->get_active()) {
      SphereObject* truncSphere = new SphereObject(offScreen);
      truncSphere->SetRadius (TruncRadius);
      truncSphere->SetPos (Center);
      truncSphere->SetColor (Vec4 (0.8, 0.8, 0.0, 0.5));
      PathVis.Objects.push_back(truncSphere);
    }

  }
  if (MultiBandButton.get_active()) {
    if (UpdateIsoVal || UpdateIsoType)
      UpdateMultiIsos();
    for (int i=0; i<VisibleBandRows.size(); i++) {
      BandRow& band = *(VisibleBandRows[i]);
      if (band.Iso != NULL)
	PathVis.Objects.push_back(band.Iso);
    }
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
WFVisualClass::SetupBandTable()
{
  kLabel.set_text ("k");
  BandLabel.set_text ("Band");
  VisibleBandTable.resize(Numk*NumBands+1,3);
  VisibleBandTable.attach (kLabel,    1, 2, 0, 1, Gtk::EXPAND, Gtk::SHRINK, 3);
  VisibleBandTable.attach (BandLabel, 2, 3, 0, 1, Gtk::EXPAND, Gtk::SHRINK, 3);
  VisibleBandRows.resize (Numk*NumBands);
  int row = 0;
  assert (Infile.OpenSection("eigenstates"));
  for (int ki=0; ki<Numk; ki++) {
    assert (Infile.OpenSection("twist", ki));
    int numSpin[2] = {0, 0};
    for (int bi=0; bi<NumBands; bi++) {
      assert (Infile.OpenSection("band", bi));
      int spin;
      Infile.ReadVar ("spin", spin, 0);
      BandRow &band = *(new BandRow(*this));
      band.Check.signal_toggled().connect 
	(sigc::bind<int>
	 (sigc::mem_fun(*this, &WFVisualClass::OnBandToggle),row));
      VisibleBandTable.attach(band.Check,     0, 1, row+1, row+2, 
			      Gtk::EXPAND, Gtk::SHRINK);
      VisibleBandTable.attach(band.kLabel,    1, 2, row+1, row+2,
			      Gtk::EXPAND, Gtk::SHRINK);
      VisibleBandTable.attach(band.BandLabel, 2, 3, row+1, row+2, 
			      Gtk::EXPAND, Gtk::SHRINK);
      char kstr[50], bstr[250];
      snprintf (kstr, 50, "%d", ki+1);
      snprintf (bstr, 250, "%d<span  font_family=\"Standard Symbols L\">&#%d;</span>", numSpin[spin]+1, 0x2191+2*spin); 
      band.kLabel.set_text(kstr);
      band.BandLabel.set_markup(bstr);
      if (numSpin[spin] < NumElectrons/2) {
	band.kLabel.modify_fg(Gtk::STATE_NORMAL, Gdk::Color("black"));
	band.BandLabel.modify_fg(Gtk::STATE_NORMAL, Gdk::Color("black"));
      }
      else {
	band.kLabel.modify_fg(Gtk::STATE_NORMAL, Gdk::Color("red"));
	band.BandLabel.modify_fg(Gtk::STATE_NORMAL, Gdk::Color("red"));
      }
      VisibleBandRows[row] = &band;
      row++;
      Infile.CloseSection (); //"band"
      numSpin[spin]++;
    }
    Infile.CloseSection(); // "twist"
  }
  Infile.CloseSection(); // "eigenstates"
  VisibleBandTable.show_all();
}


void
WFVisualClass::SetupBandTable_ESHDF()
{
  kLabel.set_text ("k");
  BandLabel.set_text ("Band");
  VisibleBandTable.resize(Numk*NumBands+1,3);
  VisibleBandTable.attach (kLabel,    1, 2, 0, 1, Gtk::EXPAND, Gtk::SHRINK, 3);
  VisibleBandTable.attach (BandLabel, 2, 3, 0, 1, Gtk::EXPAND, Gtk::SHRINK, 3);
  VisibleBandRows.resize (Numk*NumBands);
  int row = 0;
  assert (Infile.OpenSection("electrons"));
  int numSpins;
  assert (Infile.ReadVar("number_of_spins", numSpins));
  for (int ki=0; ki<Numk; ki++) {
    assert (Infile.OpenSection("kpoint", ki));
    for (int spin=0; spin<numSpins; spin++) {
      assert (Infile.OpenSection("spin", spin));
      for (int bi=0; bi<NumBands; bi++) {
	assert (Infile.OpenSection("state", bi));
	BandRow &band = *(new BandRow(*this));
	band.Check.signal_toggled().connect 
	  (sigc::bind<int>
	   (sigc::mem_fun(*this, &WFVisualClass::OnBandToggle),row));
	VisibleBandTable.attach(band.Check,     0, 1, row+1, row+2, 
				Gtk::EXPAND, Gtk::SHRINK);
	VisibleBandTable.attach(band.kLabel,    1, 2, row+1, row+2,
				Gtk::EXPAND, Gtk::SHRINK);
	VisibleBandTable.attach(band.BandLabel, 2, 3, row+1, row+2, 
				Gtk::EXPAND, Gtk::SHRINK);
	char kstr[50], bstr[250];
	snprintf (kstr, 50, "%d", ki+1);
	snprintf (bstr, 250, "%d<span  font_family=\"Standard Symbols L\">&#%d;</span>", bi+1, 0x2191+2*spin); 
	band.kLabel.set_text(kstr);
	band.BandLabel.set_markup(bstr);
	if (bi < NumElectrons/2) {
	  band.kLabel.modify_fg(Gtk::STATE_NORMAL, Gdk::Color("black"));
	  band.BandLabel.modify_fg(Gtk::STATE_NORMAL, Gdk::Color("black"));
	}
	else {
	  band.kLabel.modify_fg(Gtk::STATE_NORMAL, Gdk::Color("red"));
	  band.BandLabel.modify_fg(Gtk::STATE_NORMAL, Gdk::Color("red"));
	}
	VisibleBandRows[row] = &band;
	row++;
	Infile.CloseSection (); // "band"
      }
      Infile.CloseSection();    // "spin"
    }
    Infile.CloseSection();    // "kpoint"
  }
  Infile.CloseSection();      // "eigenstates"
  VisibleBandTable.show_all();
}


void
WFVisualClass::Read_ESHDF ()
{
  
  assert (Infile.OpenSection("supercell"));
  Array<double,2> lattice;
  assert (Infile.ReadVar("primitive_vectors", lattice));
  Box.Set (ToMat3(lattice));
  xPlane.SetLattice (ToMat3(lattice));
  yPlane.SetLattice (ToMat3(lattice));
  zPlane.SetLattice (ToMat3(lattice));
  Infile.CloseSection(); // "supercell"
  
  assert(Infile.OpenSection("atoms"));
  Array<double,2> pos;
  assert (Infile.ReadVar("positions", pos));
  AtomPos.resize(pos.extent(0));
  for (int i=0; i<pos.extent(0); i++) {
    AtomPos(i) = Vec3(pos(i,0), pos(i,1), pos(i,2));
    for (int j=0; j<3; j++)
      AtomPos(i) -= (0.5-Shift[j])*Box(j);
  }
  Array<int,1> species_ids;
  assert (Infile.ReadVar("species_ids", species_ids));
  int num_species;
  assert (Infile.ReadVar("number_of_species", num_species));
  Array<int,1> atomic_numbers(num_species);
  for (int isp=0; isp < num_species; isp++) {
    assert (Infile.OpenSection("species", isp));
    assert (Infile.ReadVar("atomic_number", atomic_numbers(isp)));
    Infile.CloseSection(); // "species"
  }
  AtomTypes.resize(species_ids.size());
  for (int iat=0; iat<species_ids.size(); iat++)
    AtomTypes(iat) = atomic_numbers(species_ids(iat));
  Infile.CloseSection(); // "atoms"

  assert (Infile.OpenSection("electrons"));
  TinyVector<int,2> num_elecs;
  assert (Infile.ReadVar("number_of_electrons", num_elecs));
  NumElectrons = num_elecs[0] + num_elecs[1];
  Infile.ReadVar("number_of_kpoints", Numk);
  assert (Infile.OpenSection("kpoint", 0));
  assert (Infile.OpenSection("spin", 0));
  assert (Infile.ReadVar("number_of_states", NumBands));
  Infile.CloseSection(); // "spin"
  Infile.CloseSection(); // "kpoint"

  Infile.CloseSection(); // "electrons"
  SetupBandTable_ESHDF();


  /// Read first wave function
  ReadWF_ESHDF(0,0);
  WFIso.SetCenter (uCenter, uMin, uMax);
  Xgrid.Init(-0.5, 0.5, WFData.extent(0));
  Ygrid.Init(-0.5, 0.5, WFData.extent(1));
  Zgrid.Init(-0.5, 0.5, WFData.extent(2));
  if (Nonuniform) 
    WFIso.Init(&NUXgrid, &NUYgrid, &NUZgrid, WFData, true);
  else 
    WFIso.Init(&Xgrid, &Ygrid, &Zgrid, WFData, true);
  //  WFIso.Init (-0.5, 0.5, -0.5, 0.5, -0.5, 0.5, WFData);
  WFIso.SetLattice(Box.GetLattice());
  CurrBand = 0; 
  Currk = 0;
  BandAdjust.set_upper(NumBands-1.0);
  kAdjust.set_upper(Numk-1.0);

  Vec3 a[3];
  a[0] = Vec3 (lattice(0,0), lattice(0,1), lattice(0,2));
  a[1] = Vec3 (lattice(1,0), lattice(1,1), lattice(1,2));
  a[2] = Vec3 (lattice(2,0), lattice(2,1), lattice(2,2));
  double maxDim = 1.2*
    max(sqrt(dot(a[0],a[0])), max(sqrt(dot(a[1],a[1])),
				  sqrt(dot(a[2],a[2]))));
  PathVis.View.SetDistance (1.2*maxDim);
  
  IsoButton.set_active(true);
  IsoButton.set_sensitive(true);
  IsoFrame.set_sensitive(true);
  DrawFrame();

}

void
WFVisualClass::Read(string filename)
{
  if (FileIsOpen)
    Infile.CloseFile();
  assert (Infile.OpenFile(filename));
  FileIsOpen = true;
  string format;
  Infile.ReadVar("format", format);
  IsESHDF = format == "ES-HDF";
  if (IsESHDF) 
    return Read_ESHDF();

  /// Read lattice vectors
  assert (Infile.OpenSection("parameters"));
  Array<double,2> lattice;
  assert (Infile.ReadVar("lattice", lattice));
  assert (Infile.ReadVar("num_electrons", NumElectrons));
//   if (!IsDiag(lattice)) { 
//     cerr << "wfvis does not currently work with non-orthorhombic cells.\n";
//     abort();
//     //Box.Set (ToMat3(lattice));
//   }
  //  Box.Set (lattice(0,0), lattice(1,1), lattice(2,2));
  Box.Set (ToMat3(lattice));
  xPlane.SetLattice (ToMat3(lattice));
  yPlane.SetLattice (ToMat3(lattice));
  zPlane.SetLattice (ToMat3(lattice));
  
  Infile.CloseSection(); // parameters

  /// Read ion positions
  assert (Infile.OpenSection("ions"));
  Array<double,2> pos;
  assert (Infile.ReadVar("pos", pos));
  AtomPos.resize(pos.extent(0));
//   for (int i=0; i<pos.extent(0); i++)
//     for (int j=0; j<3; j++)
//       AtomPos(i)[j] = (pos(i,j) - 0.5*Box[j]);
  for (int i=0; i<pos.extent(0); i++) {
    AtomPos(i) = Vec3(pos(i,0), pos(i,1), pos(i,2));
    for (int j=0; j<3; j++)
      AtomPos(i) -= (0.5-Shift[j])*Box(j);
  }
  assert (Infile.ReadVar("atom_types", AtomTypes));
  Infile.CloseSection (); // "ions"

  /// Count k-points and bands
  assert (Infile.OpenSection("eigenstates"));
  Numk = Infile.CountSections ("twist");
  assert (Infile.OpenSection("twist", 0));
  NumBands = Infile.CountSections("band");
  Infile.CloseSection(); // "twist"
  Infile.CloseSection(); // "eigenstates"
  SetupBandTable();

  /// Read first wave function
  ReadWF(0,0);
  WFIso.SetCenter (uCenter, uMin, uMax);
//   Xgrid.Init(-0.5*Box[0], 0.5*Box[0], WFData.extent(0));
//   Ygrid.Init(-0.5*Box[1], 0.5*Box[1], WFData.extent(1));
//   Zgrid.Init(-0.5*Box[2], 0.5*Box[2], WFData.extent(2));
  Xgrid.Init(-0.5, 0.5, WFData.extent(0));
  Ygrid.Init(-0.5, 0.5, WFData.extent(1));
  Zgrid.Init(-0.5, 0.5, WFData.extent(2));
  if (Nonuniform) 
    WFIso.Init(&NUXgrid, &NUYgrid, &NUZgrid, WFData, true);
  else 
    WFIso.Init(&Xgrid, &Ygrid, &Zgrid, WFData, true);
  //  WFIso.Init (-0.5, 0.5, -0.5, 0.5, -0.5, 0.5, WFData);
  WFIso.SetLattice(Box.GetLattice());
  CurrBand = 0; 
  Currk = 0;
  BandAdjust.set_upper(NumBands-1.0);
  kAdjust.set_upper(Numk-1.0);

  Vec3 a[3];
  a[0] = Vec3 (lattice(0,0), lattice(0,1), lattice(0,2));
  a[1] = Vec3 (lattice(1,0), lattice(1,1), lattice(1,2));
  a[2] = Vec3 (lattice(2,0), lattice(2,1), lattice(2,2));
  double maxDim = 1.2*
    max(sqrt(dot(a[0],a[0])), max(sqrt(dot(a[1],a[1])),
				  sqrt(dot(a[2],a[2]))));
  PathVis.View.SetDistance (1.2*maxDim);
  
  IsoButton.set_active(true);
  IsoButton.set_sensitive(true);
  IsoFrame.set_sensitive(true);
  DrawFrame();
}


	
void
WFVisualClass::OnIsoChange()
{
  double rho = IsoAdjust.get_value() * MaxVal;
  double rs = pow (3.0/(4.0*M_PI*rho),1.0/3.0);
  char rstext[100];
  snprintf (rstext, 100, "rs = %1.3f", rs);
  rsLabel.set_text(rstext);
  
  UpdateIso = true;
  UpdateIsoVal = true;
  DrawFrame();
}

void
WFVisualClass::OnBandChange()
{
  UpdateIso = true;
  UpdatePlane[0] = true;
  UpdatePlane[1] = true;
  UpdatePlane[2] = true;
  Glib::signal_idle().connect
    (sigc::bind<bool>(mem_fun(*this, &WFVisualClass::DrawFrame), false));
}

void
WFVisualClass::OnkChange()
{
  UpdateIso = true;
  UpdatePlane[0] = true;
  UpdatePlane[1] = true;
  UpdatePlane[2] = true;
  DrawFrame();
}

void
WFVisualClass::OnPlaneChange(int dir)
{
  if (dir==0 && xPlaneButton.get_active())
    UpdatePlane[0] = true;
  if (dir==1 && yPlaneButton.get_active())
    UpdatePlane[1] = true;
  if (dir==2 && zPlaneButton.get_active())
    UpdatePlane[2] = true;
  DrawFrame();
}


WFVisualClass::~WFVisualClass()
{

}


bool
WFVisualClass::ReadWF_ESHDF (int kpoint, int band)
{
  int is_complex(1);

  assert (Infile.OpenSection("electrons"));
  assert (Infile.ReadVar("psi_r_is_complex", is_complex));
  assert (Infile.OpenSection("kpoint", kpoint));
  Array<double,1> twist_angle;
  assert (Infile.ReadVar("reduced_k", twist_angle));
  bool gammaPoint = 
    (fabs(twist_angle(0)) < 1.0e-12) &&
    (fabs(twist_angle(1)) < 1.0e-12) &&
    (fabs(twist_angle(2)) < 1.0e-12);
  assert (Infile.OpenSection("spin", 0));
  assert (Infile.OpenSection("state", band));
  SuperTwistInt = Vec3(0.0, 0.0, 0.0);
  uMin = Vec3 (0.0, 0.0, 0.0);
  uMax = Vec3 (1.0, 1.0, 1.0);
  uCenter = Vec3 (0.0, 0.0, 0.0);
  Array<double,3> wfdata;
  if (is_complex) {
    Array<double,4> zdata;
    assert (Infile.ReadVar ("psi_r", zdata));
    if (gammaPoint && ((WFDisplay == REAL_PART) || (WFDisplay == IMAG_PART))) {
      double maxRho = 0.0;
      int ixMax, iyMax, izMax;
      for (int ix=0; ix<zdata.extent(0); ix+=5)
	for (int iy=0; iy<zdata.extent(1); iy+=5)
	  for (int iz=0; iz<zdata.extent(2); iz+=5) {
	    double rho = (zdata(ix,iy,iz,0)*zdata(ix,iy,iz,0) +
			  zdata(ix,iy,iz,1)*zdata(ix,iy,iz,1));
	    if (rho > maxRho) {
	      maxRho = rho;
	      ixMax = ix; iyMax = iy; izMax=iz;
	    }
	  }
      double phase = atan2 (zdata(ixMax, iyMax, izMax,1),
			    zdata(ixMax, iyMax, izMax,0));
      complex<double> factor (cos(phase), -sin(phase));
      for (int ix=0; ix<zdata.extent(0); ix++)
	for (int iy=0; iy<zdata.extent(1); iy++)
	  for (int iz=0; iz<zdata.extent(2); iz++) {
	    complex<double> wf(zdata(ix,iy,iz,0), zdata(ix,iy,iz,1));	  
	    wf *= factor;
	    zdata(ix,iy,iz,0) = wf.real();
	    zdata(ix,iy,iz,1) = wf.imag();
	  }
    }
    wfdata.resize(zdata.extent(0), zdata.extent(1), zdata.extent(2));
    if (WFDisplay == MAG2) 
      for (int ix=0; ix<zdata.extent(0); ix++)
	for (int iy=0; iy<zdata.extent(1); iy++)
	  for (int iz=0; iz<zdata.extent(2); iz++)
	    wfdata(ix,iy,iz) = (zdata(ix,iy,iz,0)*zdata(ix,iy,iz,0) + 
				zdata(ix,iy,iz,1)*zdata(ix,iy,iz,1));
    else if (WFDisplay == REAL_PART)
      for (int ix=0; ix<zdata.extent(0); ix++)
	for (int iy=0; iy<zdata.extent(1); iy++)
	  for (int iz=0; iz<zdata.extent(2); iz++)
	    wfdata(ix,iy,iz) = zdata(ix,iy,iz,0);
    else
      for (int ix=0; ix<zdata.extent(0); ix++)
	for (int iy=0; iy<zdata.extent(1); iy++)
	  for (int iz=0; iz<zdata.extent(2); iz++)
	    wfdata(ix,iy,iz) = zdata(ix,iy,iz,1);
  }
  else 
    assert (Infile.ReadVar ("psi_r", wfdata));
    
  int Nx = wfdata.extent(0) + 1;
  int Ny = wfdata.extent(1) + 1;
  int Nz = wfdata.extent(2) + 1;
  WFData.resize(Nx,Ny,Nz);


  MaxVal = 0.0;
  int xShift, yShift, zShift;
  xShift = (int)round(Shift[0]*wfdata.extent(0));
  yShift = (int)round(Shift[1]*wfdata.extent(1));
  zShift = (int)round(Shift[2]*wfdata.extent(2));
  Vec3 u(0.0, 0.0, 0.0);
  Vec3 du(1.0/(double)(wfdata.extent(0)),
	  1.0/(double)(wfdata.extent(1)),
	  1.0/(double)(wfdata.extent(2)));
  
  for (int ix=0; ix<Nx; ix++) {
    u[1] = 0.0;
    for (int iy=0; iy<Ny; iy++) {
      u[2] = 0.0;
      for (int iz=0; iz<Nz; iz++) {
	// WF data is store from 0 to Lx, not -Lx/2 to Lx/2
	int jx = (ix-xShift+Nx-1)%(Nx-1);
	int jy = (iy-yShift+Ny-1)%(Ny-1);
	int jz = (iz-zShift+Nz-1)%(Nz-1);
	WFData(ix,iy,iz) = wfdata(jx,jy,jz);
	MaxVal = max(MaxVal, fabs(WFData(ix,iy,iz)));
	u[2] += du[2];
      }
      u[1] += du[1];
    }
    u[0] += du[0];
  }

  Infile.CloseSection(); // "spin"
  Infile.CloseSection(); // "state"
  Infile.CloseSection(); // "kpoint"
  Infile.CloseSection(); // "electrons"

  Nonuniform = false;

  return true;
}





bool
WFVisualClass::ReadWF (int kpoint, int band)
{
  Array<double,4> wfdata;
  assert (Infile.OpenSection("eigenstates"));
  assert (Infile.OpenSection("twist", kpoint));
  Array<double,1> twist_angle;
  assert (Infile.ReadVar("twist_angle", twist_angle));
  bool gammaPoint = 
    (fabs(twist_angle(0)) < 1.0e-12) &&
    (fabs(twist_angle(1)) < 1.0e-12) &&
    (fabs(twist_angle(2)) < 1.0e-12);
  assert (Infile.OpenSection("band", band));
  Array<double,1> super_twist_int;
  if (Infile.ReadVar("super_twist_int", super_twist_int)) {
    SuperTwistInt[0] = super_twist_int(0);
    SuperTwistInt[1] = super_twist_int(1);
    SuperTwistInt[2] = super_twist_int(2);
  }
  else
    SuperTwistInt = Vec3(0.0, 0.0, 0.0);
  assert (Infile.ReadVar ("eigenvector", wfdata));
  Array<double,1> center;
  Localized = Infile.ReadVar ("center", center);
  uMin = Vec3 (0.0, 0.0, 0.0);
  uMax = Vec3 (1.0, 1.0, 1.0);
  // CHECK THIS!!!!
  uCenter = Vec3 (0.0, 0.0, 0.0);
  Array<double,1> umin, umax;
  Truncated = (Infile.ReadVar ("umin", umin) && 
	       Infile.ReadVar("umax", umax));
  if (Truncated) { 
    uMin = Vec3 (umin(0), umin(1), umin(2));
    uMax = Vec3 (umax(0), umax(1), umax(2));
  }

  Mat3 l = Box.GetLattice();
  if (Localized) {
    Center[0]=center(0); 
    Center[1]=center(1);  
    Center[2]=center(2);
    // Shift by half a box length
    Center -= 0.5*Vec3(l(0,0)+l(1,0)+l(2,0),
		       l(0,1)+l(1,1)+l(2,1),
		       l(0,2)+l(1,2)+l(2,2));
  }
  if (Truncated) {
    uCenter = Box.GetLatticeInv () * Center;
  }
  else
    uCenter = Vec3 (0.0, 0.0, 0.0);
  Localized = Localized && Infile.ReadVar ("radius", TruncRadius);

  Nonuniform = Infile.OpenSection ("xgrid");
  if (Nonuniform) {
    Array<double,1> points;
    Infile.ReadVar ("points", points);
    points = points - 0.5;
    NUXgrid.Init (points);
    Infile.CloseSection();
    Infile.OpenSection ("ygrid");
    Infile.ReadVar ("points", points);
    points = points - 0.5;
    NUYgrid.Init (points);
    Infile.CloseSection();
    Infile.OpenSection ("zgrid");
    Infile.ReadVar("points", points);
    points = points - 0.5;
    NUZgrid.Init (points);
    Infile.CloseSection();
  }


  Infile.CloseSection(); // "eigenstates"
  Infile.CloseSection(); // "twist"
  Infile.CloseSection(); // "band"
  WFData.resize(wfdata.extent(0), 
		wfdata.extent(1),
		wfdata.extent(2));
  int Nx = wfdata.extent(0);
  int Ny = wfdata.extent(1);
  int Nz = wfdata.extent(2);
  // If we're at the gamma point, adjust the phase of the WF so that
  // it is purely real.
  if (gammaPoint && ((WFDisplay == REAL_PART) || (WFDisplay == IMAG_PART))) {
    double maxRho = 0.0;
    int ixMax, iyMax, izMax;
    for (int ix=0; ix<Nx; ix+=5)
      for (int iy=0; iy<Ny; iy+=5)
	for (int iz=0; iz<Nz; iz+=5) {
	  double rho = (wfdata(ix,iy,iz,0)*wfdata(ix,iy,iz,0) +
			wfdata(ix,iy,iz,1)*wfdata(ix,iy,iz,1));
	  if (rho > maxRho) {
	    maxRho = rho;
	    ixMax = ix; iyMax = iy; izMax=iz;
	  }
	}
    double phase = atan2 (wfdata(ixMax, iyMax, izMax,1),
			  wfdata(ixMax, iyMax, izMax,0));
    complex<double> factor (cos(phase), -sin(phase));
    for (int ix=0; ix<wfdata.extent(0); ix++)
      for (int iy=0; iy<wfdata.extent(1); iy++)
	for (int iz=0; iz<wfdata.extent(2); iz++) {
	  complex<double> wf(wfdata(ix,iy,iz,0), wfdata(ix,iy,iz,1));	  
	  wf *= factor;
	  wfdata(ix,iy,iz,0) = wf.real();
	  wfdata(ix,iy,iz,1) = wf.imag();
	}
  }
    
  MaxVal = 0.0;
  int xShift, yShift, zShift;
  xShift = (int)round(Shift[0]*wfdata.extent(0));
  yShift = (int)round(Shift[1]*wfdata.extent(1));
  zShift = (int)round(Shift[2]*wfdata.extent(2));
  Vec3 u(0.0, 0.0, 0.0);
  Vec3 du(1.0/(double)(wfdata.extent(0)-1),
	  1.0/(double)(wfdata.extent(1)-1),
	  1.0/(double)(wfdata.extent(2)-1));
  
  for (int ix=0; ix<wfdata.extent(0); ix++) {
    u[1] = 0.0;
    for (int iy=0; iy<wfdata.extent(1); iy++) {
      u[2] = 0.0;
      for (int iz=0; iz<wfdata.extent(2); iz++) {
	// WF data is store from 0 to Lx, not -Lx/2 to Lx/2
	int jx = (ix-xShift+Nx-1/*+Nx/2*/)%(Nx-1);
	int jy = (iy-yShift+Ny-1/*+Ny/2*/)%(Ny-1);
	int jz = (iz-zShift+Nz-1/*+Nz/2*/)%(Nz-1);
	complex<double> zval(wfdata(jx,jy,jz,0), wfdata(jx,jy,jz,1));

	if (false) {
	  double s,c;
	  sincos(-2.0*M_PI*dot(u,SuperTwistInt), &s, &c);
	  zval *= complex<double>(c,s);
	}

	if (WFDisplay == MAG2)
	  WFData(ix,iy,iz) = norm(zval);
	else if (WFDisplay == REAL_PART)
	  WFData(ix,iy,iz) = zval.real();
	else if (WFDisplay == IMAG_PART)
	  WFData(ix,iy,iz) = zval.imag();
	MaxVal = max(MaxVal, fabs(WFData(ix,iy,iz)));
	u[2] += du[2];
      }
      u[1] += du[1];
    }
    u[0] += du[0];
  }
  return true;
}



void
WFVisualClass::OnCoordToggle()
{
  DrawFrame();
}

void
WFVisualClass::OnSphereToggle()
{
  DrawFrame();
}


void
WFVisualClass::OnBondsToggle()
{
  DrawFrame();
}


void
WFVisualClass::OnBoxToggle()
{
  DrawFrame();
}

void
WFVisualClass::OnTruncRadiiToggle()
{
  DrawFrame();
}

void
WFVisualClass::OnIsocontourToggle()
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
WFVisualClass::OnRadiusChange()
{
  DrawFrame();
}

void
WFVisualClass::OnDisplayRadio(WFDisplayType type)
{
  WFDisplayType newtype;
  if (RealRadio->get_active() && type==REAL_PART)
    newtype = REAL_PART;
  if (ImagRadio->get_active() && type==IMAG_PART)
    newtype = IMAG_PART;
  if (Mag2Radio->get_active() && type==MAG2)
    newtype = MAG2;

  if (WFDisplay != newtype) {
    WFDisplay = newtype;
    if (IsESHDF)
      ReadWF_ESHDF(Currk,CurrBand);
    else
      ReadWF (Currk, CurrBand);
    UpdatePlane[0] = UpdatePlane[1] = UpdatePlane[2] = true;
    UpdateIso = true;  ResetIso = true;
    DrawFrame ();
  }
  UpdateIsoType = true;
  if (MultiBandButton.get_active())
    DrawFrame();
}


bool
WFVisualClass::WriteState(string fname)
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
  out.WriteVar ("MultiBands", MultiBandButton.get_active());
  Array<bool,1> visibleBands (VisibleBandRows.size());
  for (int i=0; i<VisibleBandRows.size(); i++)
    visibleBands(i) = VisibleBandRows[i]->Check.get_active();
  out.WriteVar ("VisibleBands", visibleBands);

  string wftype;
  if (WFDisplay == REAL_PART)
    wftype = "real";
  if (WFDisplay == IMAG_PART)
    wftype = "imag";
  if (WFDisplay == MAG2)
    wftype = "mag2";
  out.WriteVar ("WFDisplay", wftype);
  out.CloseSection();

  out.WriteVar ("kPoint", (int)round(kAdjust.get_value()));
  out.WriteVar ("Band", (int)round(BandAdjust.get_value()));
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
WFVisualClass::ReadState (string fname)
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
  if (in.ReadVar ("MultiBands", active))
    MultiBandButton.set_active(active);
  Array<bool,1> visibleBands;
  if (in.ReadVar ("VisibleBands", visibleBands))
    for (int i=0; i<visibleBands.size(); i++) 
      VisibleBandRows[i]->Check.set_active(visibleBands(i));


  string wftype;
  assert(in.ReadVar ("WFDisplay", wftype));
  if (wftype == "real") {
    RealRadio->set_active(true);
    WFDisplay = REAL_PART;
  }
  else if (wftype == "imag") {
    ImagRadio->set_active(true);
    WFDisplay = IMAG_PART;
  }
  else if (wftype == "mag2") {
    Mag2Radio->set_active(true);
    WFDisplay = MAG2;
  }
  in.CloseSection(); // flags

  int k, band;
  assert(in.ReadVar ("kPoint", k));
  assert(in.ReadVar ("Band", band));
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

  kAdjust.set_value((double)k);
  BandAdjust.set_value((double)band);
  assert(in.OpenSection("View"));
  PathVis.View.Read(in);
  in.CloseSection();
  in.CloseFile();
  Glib::signal_idle().connect
    (sigc::bind<bool>(mem_fun(*this, &WFVisualClass::DrawFrame), false));
  return true;
}

void
WFVisualClass::OnSaveState()
{
  int result = SaveStateChooser.run();
  if (result == Gtk::RESPONSE_OK) {
    string fname = SaveStateChooser.get_filename();
    WriteState (fname);
  }
  SaveStateChooser.hide();
}

void
WFVisualClass::OnOpenState()
{
  int result = OpenStateChooser.run();
  if (result == Gtk::RESPONSE_OK) {
    string fname = OpenStateChooser.get_filename();
    ReadState (fname);
  }
  OpenStateChooser.hide();
}  



void
WFVisualClass::OnMultiBandToggle()
{
  if (MultiBandButton.get_active()) {
    IsoButton.set_active(false);
    xPlaneButton.set_active(false);
    yPlaneButton.set_active(false);
    zPlaneButton.set_active(false);
    IsoButton.set_sensitive (false);
    PlaneFrame.set_sensitive(false);
    BandFrame.set_sensitive(false);
    kFrame.set_sensitive(false);
  }
  else {
    IsoButton.set_sensitive (true);
    PlaneFrame.set_sensitive(true);
    BandFrame.set_sensitive (true);
    kFrame.set_sensitive    (true);
  }
}

void
WFVisualClass::OnBandToggle (int row)
{
  int ki = row/NumBands;
  int bi = row % NumBands;
  BandRow &band = *(VisibleBandRows[row]);
  if (band.Check.get_active()) {
    if (band.Iso == NULL) {
      band.Iso = new Isosurface;
      if (FileIsOpen) {
	if (IsESHDF)
	  ReadWF_ESHDF (ki, bi);
	else
	  ReadWF (ki, bi);
	Xgrid.Init(-0.5, 0.5, WFData.extent(0));
	Ygrid.Init(-0.5, 0.5, WFData.extent(1));
	Zgrid.Init(-0.5, 0.5, WFData.extent(2));
	band.Iso->Init(&Xgrid, &Ygrid, &Zgrid, WFData, true);
	if (Truncated)
	  band.Iso->SetCenter (uCenter, uMin, uMax);
	//	band.Iso.Init (-0.5, 0.5, -0.5, 0.5, -0.5, 0.5, WFData);
	band.Iso->SetLattice (Box.GetLattice());
	if (WFDisplay == MAG2) {
	  band.Iso->SetColor (0.0, 0.8, 0.0);
	  band.Iso->SetIsoval(MaxVal*IsoAdjust.get_value());
	}
	else {
	  vector<TinyVector<double,3> > colors;
	  colors.push_back(TinyVector<double,3> (0.8, 0.0, 0.0));
	  colors.push_back(TinyVector<double,3> (0.0, 0.0, 0.8));
	  band.Iso->SetColor(colors);
	  vector<double> vals;
	  vals.push_back(+MaxVal*sqrt(IsoAdjust.get_value()));
	  vals.push_back(-MaxVal*sqrt(IsoAdjust.get_value()));
	  band.Iso->SetIsoval(vals);
	}
      }
      band.Iso->Dynamic = false;
    } 
  }
  else {
    if (band.Iso != NULL) {
      delete band.Iso;
      band.Iso = NULL;
    }
  }
  if (MultiBandButton.get_active())
    DrawFrame();
}

void
WFVisualClass::UpdateMultiIsos()
{
  for (int i=0; i<VisibleBandRows.size(); i++) {
    int ki = i / NumBands;
    int bi = i % NumBands;
    if (VisibleBandRows[i]->Iso != NULL) {
      Isosurface &iso = *(VisibleBandRows[i]->Iso);
      if (UpdateIsoType) 
	if (IsESHDF)
	  ReadWF_ESHDF(ki, bi);
	else
	  ReadWF(ki, bi);
      if (UpdateIsoType || UpdateIsoVal) {
	if (WFDisplay == MAG2) {
	  iso.SetCenter (uCenter, uMin, uMax);
	  iso.SetColor (0.0, 0.8, 0.0);
	  iso.SetIsoval(MaxVal*IsoAdjust.get_value());
	}
	else {
	  vector<TinyVector<double,3> > colors;
	  colors.push_back(TinyVector<double,3> (0.8, 0.0, 0.0));
	  colors.push_back(TinyVector<double,3> (0.0, 0.0, 0.8));
	  iso.SetColor(colors);
	  vector<double> vals;
	  vals.push_back(+MaxVal*sqrt(IsoAdjust.get_value()));
	  vals.push_back(-MaxVal*sqrt(IsoAdjust.get_value()));
	  iso.SetIsoval(vals);
	}
      }
    }
  }
  UpdateIsoVal = false;
  UpdateIsoType = false;
}

void
WFVisualClass::SetShift(Vec3 shift)
{
  DoShift = true;
  Shift = shift;
}


void
WFVisualClass::SetViewportSize (int size)
{
  VisibleBandWindow.property_height_request().set_value(size-100);
  PathVis.set_size_request(size, size);
  resize(10,10);
}


void
WFVisualClass::OnFullscreenToggle()
{
  if (FullscreenToggle->get_active())
    fullscreen();
  else
    unfullscreen();
}


void
WFVisualClass::OnColorMapRadio (ColorMapType type)
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


void
WFVisualClass::ReadWalkers (string fname)
{
  IOSectionClass in;
  assert (in.OpenFile (fname));
  assert (in.OpenSection("state"));
  Array<double,3> w;
  assert (in.ReadVar("walkers", w));
  assert (w.extent(2) == 3);
  int nw = w.extent(0);
  int N  = w.extent(1);
  Walkers.resize(nw, N);
  for (int iw=0; iw<nw; iw++)
    for (int i=0; i<N; i++)
      for (int j=0; j<3; j++)
	Walkers(iw,i)[j] = w(iw,i,j);
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
    cerr << "Usage:\n  wfvis++ [options...] myfile.h5 [statefile.h5]\n"
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
  // std::cout << "OpenGL extension version - "
  //           << major << "." << minor << std::endl;

  // Instantiate and run the application.
  WFVisualClass wfvisual;

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
    wfvisual.SetShift (shift);
  }
  if (parser.Found("small"))
    wfvisual.SetViewportSize (600);

  wfvisual.Read (parser.GetFile(0));

  if (parser.Found("walkers"))
    wfvisual.ReadWalkers (parser.GetArg("walkers"));

  if (parser.NumFiles() == 2)
    wfvisual.ReadState (parser.GetFile(1));
  kit.run(wfvisual);

  return 0;
}
