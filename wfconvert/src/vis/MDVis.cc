#include "MDVis.h"
#include <GL/glut.h>

MDVisualClass::MDVisualClass() :
  MainVBox(false, 0), 
  QuitButton("Quit"),
  FrameAdjust(0.0, 0.0, 0.0),
  IsoAdjust  (0.01, 0.0, 1.0, 0.01, 0.1),
  BandAdjust (0.0, 0.0, 8.0, 1.0, 1.0),
  SpeedAdjust(5.0, 1.0, 10.0),
  CurrentFrame(0),
  PlayDirection(1),
  TimeoutDelay(10),
  MDExport(*this),
  UpToDate(true),
  FileIsOpen(false)
{
  FrameScale.set_adjustment(FrameAdjust);
  FrameScale.signal_value_changed().connect
    (sigc::mem_fun(*this, &MDVisualClass::OnFrameChange));

  SpeedScale.set_adjustment(SpeedAdjust);
  SpeedScale.signal_value_changed().connect
    (sigc::mem_fun(*this, &MDVisualClass::OnSpeedChange));
  SpeedFrame.set_label("Speed");
  SpeedFrame.add(SpeedScale);
  SpeedScale.property_width_request().set_value(75);
  TimeoutDelay = (int) round (20.0/SpeedAdjust.get_value());

  IsoScale.set_adjustment(IsoAdjust);
  IsoScale.set_digits(2);
  IsoScale.signal_value_changed().connect
    (sigc::mem_fun(*this, &MDVisualClass::OnIsoChange));
  IsoFrame.set_label("Density");
  IsoFrame.add(DensityBox);
  DensityBox.pack_start(rsLabel);
  DensityBox.pack_start(IsoScale);
  rsLabel.set_text("rs = ");
  IsoScale.property_width_request().set_value(75);

  BandScale.set_adjustment (BandAdjust);
  BandScale.set_digits(0);
  BandScale.signal_value_changed().connect
    (sigc::mem_fun(*this, &MDVisualClass::OnBandChange));
  BandScale.property_width_request().set_value(75);
  BandFrame.set_label("Band");
  BandFrame.add(BandScale);

  IsoBox.pack_start(IsoFrame);
  IsoBox.pack_start(BandFrame);

  PlayImage.set(FindFullPath("player_play.png"));
  PlayButton.set_icon_widget(PlayImage);
  PlayButton.set_label("Play");
  RevImage.set (FindFullPath("player_rev.png"));
  RevButton.set_icon_widget(RevImage);
  RevButton.set_label("Reverse");
  PauseImage.set (FindFullPath("player_pause.png"));
  PauseButton.set_icon_widget(PauseImage);
  PauseButton.set_label("Pause");

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

  set_reallocate_redraws(true);
  PathVis.set_size_request(800, 800);
  ////////////////////
  // Setup tool bar //
  ////////////////////
  Gtk::RadioButtonGroup group = PauseButton.get_group();
  RevButton.set_group(group);
  group = RevButton.get_group();
  PlayButton.set_group(group);

  group = OrthoButton.get_group();
  PerspectButton.set_group(group);
  
  Tools.append(RevButton);
  Tools.append(PauseButton);
  Tools.append(PlayButton);
  Tools.append(ToolSep1);
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
		sigc::mem_fun(*this, &MDVisualClass::OnOpen));
  Actions->add (Gtk::Action::create("Export", "_Export Image"),
		sigc::mem_fun(*this, &MDVisualClass::OnExport));
  Actions->add (Gtk::Action::create("Quit", "_Quit"),
		sigc::mem_fun(*this, &MDVisualClass::Quit));
  Actions->add (Gtk::Action::create("MenuView", "View"));
  Actions->add (Gtk::Action::create("Reset", "Reset"),
		sigc::mem_fun(*this, &MDVisualClass::OnViewReset));

  Glib::ustring ui_info =
    "<ui>"
    "  <menubar name='MenuBar'>"
    "    <menu action='MenuFile'>"
    "      <menuitem action='Open'/>"
    "      <menuitem action='Export'/>"
    "      <separator/>"
    "      <menuitem action='Quit'/>"
    "    </menu>"
    "    <menu action='MenuView'>"
    "      <menuitem action='Reset'/>"
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

  /////////////////////
  // Connect signals //
  /////////////////////
  PlayButton.signal_toggled().connect
    (sigc::mem_fun(*this, &MDVisualClass::OnPlayToggle));
  RevButton.signal_toggled().connect
    (sigc::mem_fun(*this, &MDVisualClass::OnPlayToggle));
  PauseButton.signal_toggled().connect
    (sigc::mem_fun(*this, &MDVisualClass::OnPlayToggle));
  OrthoButton.signal_toggled().connect
    (sigc::mem_fun(*this, &MDVisualClass::OnPerspectiveToggle));
  ClipButton.signal_toggled().connect
    (sigc::mem_fun(*this, &MDVisualClass::OnClipToggle));
  IsoButton.signal_toggled().connect
    (sigc::mem_fun(*this, &MDVisualClass::OnIsoToggle));
  QuitButton.signal_clicked().connect
    (sigc::mem_fun(*this, &MDVisualClass::Quit));

  ////////////////////
  // Pack the boxes //
  ////////////////////
  ToolBox.pack_start(Tools);
  ToolBox.pack_start(IsoBox, Gtk::PACK_SHRINK,20);
  ToolBox.pack_start(SpeedFrame, Gtk::PACK_SHRINK,20);
  MainVBox.pack_start(*Manager->get_widget("/MenuBar"), Gtk::PACK_SHRINK,0);
  MainVBox.pack_start(ToolBox, Gtk::PACK_SHRINK, 0);
  MainVBox.pack_start(PathVis);
  MainVBox.pack_start(FrameScale, Gtk::PACK_SHRINK, 0);
  MainVBox.pack_start(QuitButton, Gtk::PACK_SHRINK, 0);

  add (MainVBox);
  set_title ("mdvis++");
  show_all();
}

void
MDVisualClass::OnViewReset()
{
  PathVis.View.Reset();
  PathVis.Invalidate();
}


void
MDVisualClass::OnOpen()
{

}

void
MDVisualClass::Quit()
{
  Gtk::Main::quit();
}

void
MDVisualClass::OnExport()
{
  MDExport.SetupWidgets();
  MDExport.show_all();
}





string 
MDVisualClass::FindFullPath(string filename)
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
MDVisualClass::OnPlayToggle()
{
  TimeoutConnection.disconnect();
  if (PlayButton.get_active()) 
    PlayDirection = 1;
  else if (RevButton.get_active()) {
    PlayDirection = -1;
  }
  else
    return;
  
  /// Connect idle timeout
  sigc::slot<bool> slot = sigc::mem_fun(*this, &MDVisualClass::OnTimeout);
  TimeoutConnection = Glib::signal_timeout().connect (slot, TimeoutDelay);

}

void
MDVisualClass::OnClipToggle()
{
  Glib::signal_idle().connect
    (sigc::bind<bool>(mem_fun(*this, &MDVisualClass::DrawFrame), false));
}

void
MDVisualClass::OnIsoToggle()
{
  Glib::signal_idle().connect
    (sigc::bind<bool>(mem_fun(*this, &MDVisualClass::DrawFrame), false));
}

bool
MDVisualClass::OnTimeout()
{
  if (UpToDate) {
    UpToDate = false;
    Glib::signal_idle().connect
      (sigc::bind<bool>(mem_fun(*this, &MDVisualClass::DrawFrame), false));
  }
  //  DrawFrame();
  // PathVis.Invalidate();
  if (PlayDirection == 1) {
    if (CurrentFrame < (Trajectory.extent(0)-1)) {
      CurrentFrame++;
      //      FrameAdjust.set_value(CurrentFrame);
      return true;
    }
    else {
      PauseButton.set_active();
      return false;
    }
  }
  else if (PlayDirection == -1) {
    if (CurrentFrame > 0) {
      CurrentFrame--;
      //      FrameAdjust.set_value(CurrentFrame);
      return true;
    }
    else {
      PauseButton.set_active();
      return false;
    }
  }
}

void
MDVisualClass::OnPerspectiveToggle()
{
  bool persp = !OrthoButton.get_active();
  cerr << "Now using " << (persp ? "perspective" : "orthographic") 
       << " projection.\n";
  PathVis.View.SetPerspective(persp);
  //  PathVis.Invalidate();
  DrawFrame();
  //  PathVis.Invalidate();
}

void
MDVisualClass::OnFrameChange()
{
  int tempFrame =  (int)round(FrameScale.get_value());
  if (IsoButton.get_active()) {
    if (RhoVar != NULL) {
      RhoVar->Read(RhoData, tempFrame, 
		   Range::all(), Range::all(), Range::all()); 
    }
    if (BandRhoVar != NULL)
      BandRhoVar->Read(BandRhoData, tempFrame, Range::all(), Range::all(), Range::all(), Range::all());
  }  

  if (CurrentFrame != (int)round(FrameScale.get_value())) {
    CurrentFrame = (int)round(FrameScale.get_value());
    DrawFrame();
  }
  //  PathVis.Invalidate();
}


// void
// MDVisualClass::ClipSpheres(list<Vec3>& sphereList, double radius)
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

bool
MDVisualClass::DrawFrame(bool offScreen)
{
  /// Update frame adjust
  FrameAdjust.set_value(CurrentFrame);

  bool clipping = ClipButton.get_active();
  for (int i=0; i<PathVis.Objects.size(); i++) {
    if (dynamic_cast<SphereObject*> (PathVis.Objects[i]) != NULL)
      delete dynamic_cast<SphereObject*> (PathVis.Objects[i]);
    else if (dynamic_cast<DiskObject*> (PathVis.Objects[i]) != NULL)
      delete dynamic_cast<DiskObject*> (PathVis.Objects[i]);
    else if (dynamic_cast<Isosurface*> (PathVis.Objects[i]) != NULL)
      delete dynamic_cast<Isosurface*> (PathVis.Objects[i]);
    else if (dynamic_cast<BoxObject*> (PathVis.Objects[i]) != NULL)
      delete dynamic_cast<BoxObject*> (PathVis.Objects[i]);

    //    delete (PathVis.Objects[i]);
  }
  PathVis.Objects.resize(0);
  BoxObject *boxObject = new BoxObject;
  boxObject->SetColor (0.5, 0.5, 1.0);
  boxObject->Set ((Vec3)Box, clipping);
  PathVis.Objects.push_back(boxObject);
  

  const double radius = 1.0;

  list<Vec3>::iterator iter;
  list<Vec3> sphereList;
  for (int ptcl=0; ptcl < Trajectory.extent(1); ptcl++) 
    sphereList.push_back(Vec3(Trajectory(CurrentFrame,ptcl,0),
			      Trajectory(CurrentFrame,ptcl,1),
			      Trajectory(CurrentFrame,ptcl,2)));
  if (clipping) {
    for (iter=sphereList.begin(); iter != sphereList.end(); iter++) {
      Vec3 &r = (*iter);
      bool makeDisks = false;
      if ((r[0]+radius) > 0.5*Box[0]) {
	sphereList.push_front(Vec3(r[0]-Box[0], r[1], r[2]));
	makeDisks = true;
      }
      if ((r[0]-radius) < -0.5*Box[0]) {
	sphereList.push_front(Vec3(r[0]+Box[0], r[1], r[2]));
	makeDisks = true;
      }
    }
    
    for (iter=sphereList.begin(); iter != sphereList.end(); iter++) {
      Vec3 &r = (*iter);
      if ((r[1]+radius) > 0.5*Box[1])
	sphereList.push_front(Vec3(r[0], r[1]-Box[1], r[2]));
      if ((r[1]-radius) < -0.5*Box[1])
	sphereList.push_front(Vec3(r[0], r[1]+Box[1], r[2]));
    }
    
    for (iter=sphereList.begin(); iter != sphereList.end(); iter++) {
      Vec3 &r = (*iter);
      if ((r[2]+radius) > 0.5*Box[2])
	sphereList.push_front(Vec3(r[0], r[1], r[2]-Box[2]));
      if ((r[2]-radius) < -0.5*Box[2])
	sphereList.push_front(Vec3(r[0], r[1], r[2]+Box[2]));
    }
    // Now make disks to close spheres
    for (iter=sphereList.begin(); iter != sphereList.end(); iter++) {
      Vec3 &r = (*iter);
      for (int dim=0; dim<3; dim++) {
	if ((r[dim]+radius) > 0.5*Box[dim]) {
	  double l = 0.5*Box[dim]-fabs(r[dim]);
	  double diskRad = sqrt(radius*radius-l*l);
	  DiskObject *disk1 = new DiskObject(offScreen);
	  DiskObject *disk2 = new DiskObject(offScreen);
	  disk1->SetRadius(diskRad);
	  disk2->SetRadius(diskRad);
	  disk1->SetAxis(2*dim);
	  disk2->SetAxis(2*dim+1);
	  Vec3 r1, r2;
	  r1 = r; r2 = r;
	  r1[dim] =  0.5*Box[dim];
	  r2[dim] = -0.5*Box[dim];
	  disk1->SetPos(r1);
	  disk2->SetPos(r2);
	  PathVis.Objects.push_back(disk1);
	  PathVis.Objects.push_back(disk2);
	}
      }
    }
  }

  /// Add sphere objects from list
  for (iter=sphereList.begin(); iter!=sphereList.end(); iter++) {
    SphereObject *sphere = new SphereObject (offScreen);
    sphere->SetPos(*iter);
    sphere->SetBox(Box);
    sphere->SetRadius(radius);
    PathVis.Objects.push_back(sphere);
  }

  if (IsoButton.get_active()) {
    int band;
    band = (int)round(BandScale.get_value());
    bool doTotal = (BandRhoVar == NULL) || (band == BandRhoData.extent(3));
    if (doTotal) {
      Isosurface *rhoIso = new Isosurface;
      rhoIso->Init(&Xgrid, &Ygrid, &Zgrid, RhoData, true);
      rhoIso->SetIsoval(MaxRho*IsoAdjust.get_value());
      PathVis.Objects.push_back(rhoIso);
    }
    else if (BandRhoVar != NULL) {
      Array<double,3> bandData;
      bandData.reference (BandRhoData(Range::all(), Range::all(), Range::all(), band));
      Isosurface *bandIso = new Isosurface;
      bandIso->Init (&Xgrid, &Ygrid, &Zgrid, bandData);
      Vec3 color;
      int n = BandRhoData.extent(3);
      double x = (double)band/(double)(n-1);
      color[0] = exp(-10.0*x*x);
      color[1] = exp(-10.0*(x-0.5)*(x-0.5));
      color[2] = exp(-10.0*(1.0-x)*(1.0-x));
      bandIso->SetColor(color);
      bandIso->SetIsoval(MaxBandRho*IsoAdjust.get_value());
      PathVis.Objects.push_back(bandIso);
    }
  }

  PathVis.Invalidate();
  UpToDate = true;
  return false;
}



void
MDVisualClass::Read(string filename)
{
  if (FileIsOpen)
    Infile.CloseFile();
  assert (Infile.OpenFile(filename));
  FileIsOpen = true;
  Array<double,1> box;
  assert (Infile.OpenSection("System"));
  assert (Infile.ReadVar ("Box", box));
  Box.Set (box(0), box(1), box(2));
  Infile.CloseSection();

  double maxDim = max(max(box(0), box(1)), box(2));
  PathVis.View.SetDistance (1.2*maxDim);
  
  assert(Infile.OpenSection("Moves"));
  assert(Infile.OpenSection("Langevin"));
  bool extraSec = Infile.OpenSection("Langevin");
  assert(Infile.ReadVar("R", Trajectory));
  assert(Infile.ReadVar("Time", Time));

  RhoVar = Infile.GetVarPtr("Rho");
  BandRhoVar = Infile.GetVarPtr("BandRho");
  bool haveRho = RhoVar != NULL;
  MaxRho = 0.0;
  if (haveRho) {
    for (int frame=0; frame < Trajectory.extent(0); frame+=20) {
      cerr << "Checking frame " << frame << endl;
      RhoVar->Read(RhoData, frame, Range::all(), Range::all(), Range::all());
      MaxRho = max(FindMaxRho(), MaxRho);
    }
    RhoVar->Read(RhoData, 0, Range::all(), Range::all(), Range::all());
    Xgrid.Init(-0.5*Box[0], 0.5*Box[0], RhoData.extent(0));
    Ygrid.Init(-0.5*Box[1], 0.5*Box[1], RhoData.extent(1));
    Zgrid.Init(-0.5*Box[2], 0.5*Box[2], RhoData.extent(2));
  }
  bool haveBandRho = BandRhoVar != NULL;
  if (haveBandRho) {
    BandRhoVar->Read(BandRhoData, 0, Range::all(), Range::all(), Range::all(), Range::all());
    BandAdjust.set_upper(BandRhoData.extent(3));
    MaxBandRho = FindMaxBandRho();
    Xgrid.Init(-0.5*Box[0], 0.5*Box[0], RhoData.extent(0));
    Ygrid.Init(-0.5*Box[1], 0.5*Box[1], RhoData.extent(1));
    Zgrid.Init(-0.5*Box[2], 0.5*Box[2], RhoData.extent(2));
  }
  else
    BandFrame.hide();
  IsoButton.set_active(haveRho);
  IsoButton.set_sensitive(haveRho);
  IsoFrame.set_sensitive(haveRho);
    

  if (extraSec) 
    Infile.CloseSection(); // "Langevin"
  Infile.CloseSection(); // "Langevin"
  Infile.CloseSection(); // "Moves"

  FrameAdjust.set_upper(Trajectory.extent(0)-1);
  FrameScale.set_digits(0);

  DrawFrame();
}

double
MDVisualClass::FindMaxRho()
{
  double maxRho = 0.0;
  double totalRho = 0.0;
  for (int ix=0; ix<RhoData.extent(0); ix++)
    for (int iy=0; iy<RhoData.extent(1); iy++)
      for (int iz=0; iz<RhoData.extent(2); iz++) 
	maxRho = max(maxRho, RhoData(ix,iy,iz));
  for (int ix=0; ix<RhoData.extent(0)-1; ix++)
    for (int iy=0; iy<RhoData.extent(1)-1; iy++)
      for (int iz=0; iz<RhoData.extent(2)-1; iz++) 
	totalRho += RhoData(ix, iy, iz);
  
  totalRho *= Box[0]*Box[1]*Box[2]/(double)((RhoData.extent(0)-1)*
					    (RhoData.extent(1)-1)*
					    (RhoData.extent(2)-1));
  cerr << "TotalRho = " << totalRho << endl;
  return maxRho;
}

double
MDVisualClass::FindMaxBandRho()
{
  double maxRho = 0.0;
  for (int ix=0; ix<BandRhoData.extent(0); ix++)
    for (int iy=0; iy<BandRhoData.extent(1); iy++)
      for (int iz=0; iz<BandRhoData.extent(2); iz++)
	for (int band=0; band<BandRhoData.extent(3); band++)
	  maxRho = max(maxRho, BandRhoData(ix,iy,iz,band));
  return maxRho;
}
	

void
MDVisualClass::OnSpeedChange()
{
  TimeoutDelay = (int) round (20.0/SpeedAdjust.get_value());
  OnPlayToggle();
}

void
MDVisualClass::OnIsoChange()
{
  double rho = IsoAdjust.get_value() * MaxRho;
  double rs = pow (3.0/(4.0*M_PI*rho),1.0/3.0);
  char rstext[100];
  snprintf (rstext, 100, "rs = %1.3f", rs);
  rsLabel.set_text(rstext);
  
  DrawFrame();
}

void
MDVisualClass::OnBandChange()
{
  DrawFrame();
}


MDVisualClass::~MDVisualClass()
{

}

int main(int argc, char** argv)
{
  Gtk::Main kit(argc, argv);

  // Init gtkglextmm.
  Gtk::GL::init(argc, argv);
  glutInit(&argc, argv);

  if (argc < 2) {
    cerr << "Usage:\n  mdvis++ myfile.h5\n";
    exit (1);
  }
  
  // Query OpenGL extension version.
  int major, minor;
  Gdk::GL::query_version(major, minor);
  std::cout << "OpenGL extension version - "
            << major << "." << minor << std::endl;

  // Instantiate and run the application.
  MDVisualClass mdvisual;
  mdvisual.Read (argv[1]);
  kit.run(mdvisual);

  return 0;
}
