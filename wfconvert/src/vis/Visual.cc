#include "Visual.h"
#include "ParseCommand.h"
#include <GL/glut.h>
#include "CoordObject.h"

void 
VisualClass::ReadFrameData(int frame)
{
  if (HaveANodeData) {
    if (ANodeVar != NULL)
      ANodeVar->Read(ANodeData, frame, Range::all(),Range::all(),Range::all());
    if (BNodeVar != NULL)
      BNodeVar->Read(BNodeData, frame, Range::all(),Range::all(),Range::all());
  }
  if (RhoVar != NULL) {
    RhoVar->Read(RhoData, frame, Range::all(), Range::all(), Range::all());
    /// HACKY inversion of nodes
    //     Array<double,3> rho;
    //     RhoVar->Read(rho, frame, Range::all(), Range::all(), Range::all());
    //     RhoData.resize(rho.shape());
    //     int nx = rho.extent(0);
    //     int ny = rho.extent(1);
    //     int nz = rho.extent(2);
    //     for (int ix=0; ix<nx; ix++)
    //       for (int iy=0; iy<ny; iy++)
    // 	for (int iz=0; iz<nz; iz++)
    // 	  RhoData(ix,iy,iz) = rho(nx-ix-1,ny-iy-1,nz-iz-1);
    //     RhoData = rho;

    RhoData = log (RhoData);

    MaxRho=0.0;
    MinRho = 1.0e100;
    for (int ix=0; ix<RhoData.extent(0); ix++)
      for (int iy=0; iy<RhoData.extent(1); iy++)
	for (int iz=0; iz<RhoData.extent(2); iz++) {
	  MaxRho = max(RhoData(ix,iy,iz), MaxRho);
	  MinRho = min(RhoData(ix,iy,iz), MinRho);
	}
  }
}

void VisualClass::MakeFrame(int frame, bool offScreen)
{
  int numSpecies = Species.size();

  for (int i=0; i<PathVis.Objects.size(); i++)
    if (PathVis.Objects[i]->Dynamic)
      delete PathVis.Objects[i];
  
  PathVis.Objects.resize(0);
  CoordObject *coord = new CoordObject;
  Vec3 box = Box(0) + Box(1) + Box(2);
  coord->Set (box);
  PathVis.Objects.push_back(coord);

  ReadFrameData (frame);

  MakePaths(frame);

  int numPtcls = PathArray.extent(1);
  int numSlices = PathArray.extent(2);


//   Array<Vec3,1> onePath(numSlices);
//   for (int si=0; si<numSpecies; si++) {
//     for (int ptcl=Species(si).FirstParticle; 
// 	 ptcl<=Species(si).LastParticle; ptcl++) {
//       if (Species(si).lambda != 0.0) {
// 	for (int slice=0; slice<numSlices; slice++) {
// 	  onePath(slice)[0] = PathArray(frame,ptcl,slice,0);
// 	  onePath(slice)[1] = PathArray(frame,ptcl,slice,1);
// 	  onePath(slice)[2] = PathArray(frame,ptcl,slice,2);
// 	}
// 	PathObject* pathObj = new PathObject;
// 	if (si == 0)
// 	  pathObj->SetColor (0.3, 0.3, 1.0);
// 	else
// 	  pathObj->SetColor (0.0, 1.0, 0.0);
// 	if (PathType == TUBES)
// 	  pathObj->TubesSet (onePath);
// 	else
// 	  pathObj->LinesSet (onePath);
// 	PathVis.Objects.push_back(pathObj);
//       }
//       else {
// 	Vec3 pos;
// 	pos[0] = PathArray(frame, ptcl, 0, 0);
// 	pos[1] = PathArray(frame, ptcl, 0, 1);
// 	pos[2] = PathArray(frame, ptcl, 0, 2);

// 	SphereObject* sphere = new SphereObject;
// 	sphere->SetPos (pos);
// 	sphere->SetColor (Vec3(1.0, 0.0, 1.0));
// 	PathVis.Objects.push_back(sphere);
//       }
//     }
//   }
  for (int li=0; li<Paths.size(); li++) {
    PathObject* pathObj = new PathObject;
// // <<<<<<< .mine
// //     pathObj->Closed = Paths[li]->Closed;
// //     if (Paths[li]->Closed)
// //       pathObj->SetColor (0.0, 0.0, 1.0);
// //     else
// //       pathObj->SetColor (1.0, 0.0, 0.0);
// =======
    pathObj->Closed = Paths[li]->Closed;
    if (Paths[li]->Winding)
      pathObj->SetColor (0.6, 0.6, 0.0);
    else
      pathObj->SetColor (0.0, 0.0, 1.0);
//     if (li < 8)
//       pathObj->SetColor (0.0, 0.0, 1.0);
//     else
//       pathObj->SetColor (0.7, 0.7, 0.0);
// >>>>>>> .r523
    pathObj->SetRadius (min(min(Box[0], Box[1]), Box[2])*0.003);
    if (PathType == TUBES)
      pathObj->TubesSet (Paths[li]->Path);
    else
      pathObj->LinesSet (Paths[li]->Path);
    PathVis.Objects.push_back(pathObj);
  }
  for (int si=0; si<numSpecies; si++) {
    if(!processH2O || (processH2O && Species(si).Name!="e")){
      if (Species(si).lambda == 0){
        for (int ptcl=Species(si).FirstParticle; ptcl<=Species(si).LastParticle;
	   ptcl++) {
	        SphereObject* sphere = new SphereObject(offScreen);
	        dVec pos;
	        pos[0] = PathArray(frame, ptcl, 0, 0);
 	        pos[1] = PathArray(frame, ptcl, 0, 1);
 	        pos[2] = PathArray(frame, ptcl, 0, 2);
	        /// HACK HACK HACK HACK
//         	pos[0] += 0.5*Box[0];
//         	pos[1] += 0.5*Box[1];
//         	pos[2] += 0.5*Box[2];
	        Box.PutInBox(pos);
	        sphere->SetPos (pos);
	        sphere->SetBox(Box);
//         	if (ptcl==0)
//         	  sphere->SetColor (Vec3(1.0, 0.0, 1.0));
//         	else
            Vec3 colorCode(1.0, 0.0, 0.0);
            if(processH2O && (Species(si).Name == "p" || Species(si).Name == "H")){
              colorCode(0) = 0.8;
              colorCode(1) = 1.0;
              colorCode(2) = 0.8;
            }
	          //sphere->SetColor (Vec3(1.0, 0.0, 0.0));
	          sphere->SetColor (colorCode);
	        PathVis.Objects.push_back(sphere);
        }
      }
    }
  }

  ///Addition of sphere for helium droplets
  if (false) {
    SphereObject* sphere = new SphereObject;
    dVec pos;
    pos[0] = 0;
    pos[1] = 0;
    pos[2] = 0;
    sphere->SetPos (pos);
    sphere->SetColor (Vec3(1.0, 0.0, 0.0));
    sphere->SetRadius(31.0);
    //  PathVis.Objects.push_back(sphere);
    ///end addition of sphere for helium droplets
  }
  
  BoxObject *boxObject = new BoxObject;
  boxObject->SetColor (0.5, 0.5, 1.0);
  boxObject->Set (Box[0], Box[1], Box[2]);
  PathVis.Objects.push_back(boxObject);


  /// HACKED IN ISOSURFACE OBJECT
//   Isosurface *isoPtr = new Isosurface;
//   Isosurface &iso = *(isoPtr);
//   LinearGrid *xgrid = new LinearGrid(-8.0, 8.0, 10);
//   LinearGrid *ygrid = new LinearGrid(-8.0, 8.0, 10);
//   LinearGrid *zgrid = new LinearGrid(-8.0, 8.0, 10);
//   Array<double,3> initData(10,10,10);
//   for (int ix=0; ix<10; ix++) {
//     double x = (*xgrid)(ix);
//     for (int iy=0; iy<10; iy++) {
//       double y = (*ygrid)(iy);
//       for (int iz=0; iz<10; iz++) {
// 	double z = (*zgrid)(iz);
// 	initData(ix,iy,iz) = x*x+y*y+z*z;
//       }
//     }
//   }
//   iso.Init (xgrid, ygrid, zgrid, initData);
//   iso.SetIsoval (36.0);

//   PathVis.Objects.push_back(isoPtr);
  if (RhoVar != NULL) {
    Isosurface *isoPtr = new Isosurface;
    Isosurface &iso = *isoPtr;
    Xgrid.Init(-0.5, 0.5, RhoData.extent(0));
    Ygrid.Init(-0.5, 0.5, RhoData.extent(1));
    Zgrid.Init(-0.5, 0.5, RhoData.extent(2));
    iso.Init (&Xgrid, &Ygrid, &Zgrid, RhoData, true);
    Mat3 lattice;
    lattice =  Box[0], 0.0, 0.0, 0.0, Box[1], 0.0, 0.0, 0.0, Box[2];
    Vec3 uCenter(0.0, 0.0, 0.0), uMin(0.0, 0.0, 0.0), uMax(1.0, 1.0, 1.0);
    iso.SetLattice (lattice);
    iso.SetCenter (uCenter, uMin, uMax);
    iso.SetIsoval(0.01*MaxRho*RhoAdjust.get_value());
    //iso.SetIsoval(7.0e-9);
    PathVis.Objects.push_back(isoPtr);
  }
  else if (HaveANodeData) {
    Isosurface *isoPtr = new Isosurface;
    Isosurface &iso = *isoPtr;
    iso.Init (&Xgrid, &Ygrid, &Zgrid, ANodeData, true);
    iso.SetIsoval(IsoAdjust.get_value());
    PathVis.Objects.push_back(isoPtr);
    // HACK HACK HACK 
    // Draw electron positions
    for (int ptcl=16; ptcl<32; ptcl++) {
      	SphereObject* sphere = new SphereObject;
	dVec pos;
	pos[0] = PathArray(frame, ptcl, NodeSlice(frame), 0);
 	pos[1] = PathArray(frame, ptcl, NodeSlice(frame), 1);
 	pos[2] = PathArray(frame, ptcl, NodeSlice(frame), 2);
	sphere->SetPos (pos);
	if (ptcl == NodePtcl(frame))
	  sphere->SetColor (Vec3(0.9, 0.9, 0.0));
	else
	  sphere->SetColor (Vec3(0.0, 0.0, 0.9));
	sphere->SetRadius(0.25);
	PathVis.Objects.push_back(sphere);
    }
    dVec pos( 0.20003,  0.70002,  1.10501);
    SphereObject* sphere = new SphereObject;
    sphere->SetPos(pos);
    sphere->SetColor(Vec3(0.9, 0.9, 0.0));
    sphere->SetRadius(1.0);
    PathVis.Objects.push_back(sphere);
  }
  if (HaveBNodeData) {
    Isosurface *isoPtr = new Isosurface;
    Isosurface &iso = *isoPtr;
    iso.SetColor (0.8, 0.0, 0.0);
    iso.Init (&Xgrid, &Ygrid, &Zgrid, BNodeData, true);
    iso.SetIsoval(IsoAdjust.get_value());

    PathVis.Objects.push_back(isoPtr);
    // HACK HACK HACK 
    // Draw electron positions
    for (int ptcl=16; ptcl<24; ptcl++) {
      	SphereObject* sphere = new SphereObject;
	dVec pos;
	pos[0] = PathArray(frame, ptcl, 0, 0);
 	pos[1] = PathArray(frame, ptcl, 0, 1);
 	pos[2] = PathArray(frame, ptcl, 0, 2);
	sphere->SetPos (pos);
	if (ptcl == 16)
	  sphere->SetColor (Vec3(0.9, 0.0, 0.9));
	else
	  sphere->SetColor (Vec3(0.0, 0.0, 0.9));
	sphere->SetRadius(0.25);
	PathVis.Objects.push_back(sphere);
    }
  }
}

void VisualClass::Read(string fileName)
{
  assert(Infile.OpenFile (fileName));

  Array<double,1> box;
  assert (Infile.OpenSection("System"));
  assert (Infile.ReadVar ("Box", box));
  Box.Set (box(0), box(1), box(2));

  double maxDim = max(max(box(0), box(1)), box(2));
  PathVis.View.SetDistance (1.2*maxDim);
  //PathVis.View.SetDistance (0.2*maxDim);

  int numSpecies = Infile.CountSections ("Species");  
  Species.resize(numSpecies);
  for (int i=0; i<numSpecies; i++)
    Species(i).FirstParticle = 0;
  for (int i=0; i<numSpecies; i++) {
    Infile.OpenSection("Species",i);
    assert (Infile.ReadVar("lambda", Species(i).lambda));
    assert (Infile.ReadVar("Name", Species(i).Name));
    assert (Infile.ReadVar("NumParticles", Species(i).NumParticles));
    for (int j=i+1; j<numSpecies; j++)
      Species(j).FirstParticle += Species(i).NumParticles;
    Species(i).LastParticle=Species(i).FirstParticle+Species(i).NumParticles-1;

    Infile.CloseSection(); // Species
  }

//   for (int i=0; i<numSpecies; i++) 
//     cerr << "Species:  Name = " << Species(i).Name << "    First Ptcl = " 
// 	 << Species(i).FirstParticle 
// 	 << "   Last Ptcl = " << Species(i).LastParticle << "\n";

  Infile.CloseSection (); // "System"

  assert(Infile.OpenSection("Observables"));
  assert(Infile.OpenSection("PathDump"));
  assert(Infile.ReadVar ("Path", PathArray));
  assert(Infile.ReadVar ("Permutation", PermArray));
  RhoVar = Infile.GetVarPtr("Rho");
  PutInBox();

  if (Infile.ReadVar ("OpenPtcl", OpenPtcl)) {
    assert (OpenPtcl.size() == PathArray.extent(0));
    assert (Infile.ReadVar ("TailLocation", Tail));
    assert (Tail.extent(0) == PathArray.extent(0));
  }
  else {
    OpenPtcl.resize (PathArray.extent(0));
    OpenPtcl = -1;
  }
  HaveANodeData = Infile.OpenSection("Xgrid");
  if (HaveANodeData) {
    HaveWarpPos = Infile.ReadVar("WarpPos", WarpPos);
    Xgrid.Read(Infile);
    Infile.CloseSection();
    assert(Infile.OpenSection("Ygrid"));
    Ygrid.Read(Infile);
    Infile.CloseSection();
    assert(Infile.OpenSection("Zgrid"));
    Zgrid.Read(Infile);
    Infile.CloseSection();
    Infile.ReadVar("NodePtcl", NodePtcl);
    Infile.ReadVar("NodeSlice", NodeSlice);
    ANodeVar = Infile.GetVarPtr ("ANodes");
    if (ANodeVar == NULL)
      assert((ANodeVar = Infile.GetVarPtr("Nodes")) != NULL);
    BNodeVar = Infile.GetVarPtr ("BNodes");
    HaveBNodeData = (BNodeVar != NULL);

    double maxVal = -1.0e100;
    double minVal = 1.0e100;
    //    for (int frame=0; frame<ANodeVar->GetExtent(0); frame++) {
    //    ReadFrameData(frame);
    ReadFrameData(0);
    for (int ix=0; ix<ANodeData.extent(0); ix++)
      for (int iy=0; iy<ANodeData.extent(1); iy++)
	for (int iz=0; iz<ANodeData.extent(2); iz++) {
	  double mx, mn;
	  if (HaveBNodeData) {
	    mx = max(ANodeData(ix,iy,iz), BNodeData(ix,iy,iz));
	    mn = min(ANodeData(ix,iy,iz), BNodeData(ix,iy,iz));
	  }
	  else {
	    mx = ANodeData(ix,iy,iz);
	    mn = ANodeData(ix,iy,iz);
	  }
	  maxVal = max(maxVal, mx);
	  minVal = min(minVal, mn);
	}
    //    }
    IsoAdjust.set_lower(minVal);
    IsoAdjust.set_upper(maxVal);
  }
  else
    HaveBNodeData = false;

  /// Make sure images are continuous from frame to frame
  for (int frame=0; frame<PathArray.extent(0)-1; frame++) 
    for (int ptcl=0; ptcl<PathArray.extent(1); ptcl++) 
      for (int dim=0; dim<3; dim++) {
	while ((PathArray(frame+1,ptcl,0,dim)-PathArray(frame,ptcl,0,dim)) 
	       > 0.5*Box[dim])
	  for (int slice=0; slice<PathArray.extent(2); slice++)
	    PathArray(frame+1,ptcl,slice,dim) -= Box[dim];
	while ((PathArray(frame+1,ptcl,0,dim)-PathArray(frame,ptcl,0,dim)) 
	       < -0.5*Box[dim])
	  for (int slice=0; slice<PathArray.extent(2); slice++)
	    PathArray(frame+1,ptcl,slice,dim) += Box[dim];
      }
  FileIsRead = true;
	       

  FrameAdjust.set_upper(PathArray.extent(0)-1);
  DetailAdjust.set_upper(PathArray.extent(2)/2);
  
  Infile.CloseSection();
  FrameScale.set_value(0.0);
  MakeFrame (0);

}

void VisualClass::SetFlag(string newFlag){
  if(newFlag == "water"){
    processH2O = true;
    cerr << "setting processH2O true " << processH2O << endl;
  }
}

void VisualClass::PutInBox()
{
  for (int frame=0; frame<PathArray.extent(0); frame++) 
    for (int ptcl=0; ptcl<PathArray.extent(1); ptcl++) 
      for (int dim=0; dim<3; dim++) {
	while (PathArray(frame,ptcl,0,dim) > 0.5*Box[dim])
	  for (int slice=0; slice<PathArray.extent(2); slice++)
	    PathArray(frame,ptcl,slice,dim) -= Box[dim];
	while (PathArray(frame,ptcl,0,dim) < -0.5*Box[dim])
	  for (int slice=0; slice<PathArray.extent(2); slice++)
	    PathArray(frame,ptcl,slice,dim) += Box[dim];
      }
	  
}


#include "tubes.xpm"
#include "lines.xpm"

string VisualClass::FindFullPath(string filename)
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


VisualClass::VisualClass()
  : m_VBox(false, 0), m_ButtonQuit("Quit"), 
    FrameAdjust (0.0, 0.0, 0.0),
    PathType (LINES),
//     TubesImage("tubes.png"), LinesImage("lines.png"),
//     StraightImage("straight.png"), SmoothImage("smooth.png"),
//     NoWrapImage("nowrap2.png"), WrapImage("wrap.png"),
//     OrthoImage("orthographic.png"), PerspectImage("perspective.png"),
    FileChooser ("Choose an output file"),
    Wrap(false), Smooth(false), 
    DetailFrame ("Detail"),  DetailAdjust (1.0, 1.0, 2.0), 
    IsoFrame    ("Isosurf"),    IsoAdjust (0.0, 0.0, 5.0, 0.1),
    RhoFrame    (/*"Rho"*/"V (% max)"),        RhoAdjust (0.0, 0.0, 100.0, 1.0, 10.0),
    Export(*this), RhoVar(NULL), FileIsRead(false)
{
  TubesImage.set(FindFullPath("tubes.png"));
  LinesImage.set(FindFullPath("lines.png"));
  StraightImage.set(FindFullPath("straight.png"));
  SmoothImage.set(FindFullPath("smooth.png"));
  NoWrapImage.set(FindFullPath("nowrap2.png"));
  WrapImage.set(FindFullPath("wrap.png"));
  OrthoImage.set(FindFullPath("orthographic.png"));
  PerspectImage.set(FindFullPath("perspective.png"));

  // Top-level window.
  set_title("VisualClass");

  // Get automatically redrawn if any of their children changed allocation.
  set_reallocate_redraws(true);
  add(m_VBox);

  // VisualClass OpenGL scene.
  // PathVis.set_size_request(400, 400);
  PathVis.set_size_request(800, 800);


  // VisualClass quit button.

  m_ButtonQuit.signal_clicked().connect
    (sigc::mem_fun(*this, &VisualClass::Quit));
  FrameScale.set_adjustment (FrameAdjust);
  FrameScale.signal_value_changed().connect
    (sigc::mem_fun(*this, &VisualClass::FrameChanged));
  FrameAdjust.set_step_increment(1.0);
  FrameScale.set_digits(0);

  // Setup tool bar
  Gtk::RadioButtonGroup group = LinesButton.get_group();
  TubesButton.set_group (group);
  group = StraightButton.get_group();
  SmoothButton.set_group(group);
  group = NoWrapButton.get_group();
  WrapButton.set_group(group);

  LinesButton.set_label("Lines");
  TubesButton.set_label("Tubes");
  StraightButton.set_label("Straight");
  SmoothButton.set_label("Smooth");
  NoWrapButton.set_label("No Wrap");
  WrapButton.set_label("Wrap");


  LinesButton.signal_toggled().connect
    (sigc::mem_fun(*this, &VisualClass::LineToggle));
  WrapButton.signal_toggled().connect
    (sigc::mem_fun(*this, &VisualClass::WrapToggle));
  SmoothButton.signal_toggled().connect
    (sigc::mem_fun(*this, &VisualClass::SmoothToggle));
  OrthoButton.signal_toggled().connect
    (sigc::mem_fun(*this, &VisualClass::PerspectiveToggle));

  Gtk::RadioButtonGroup group2 = PerspectButton.get_group();
  OrthoButton.set_group (group2);
  OrthoButton.set_label ("Ortho");
  PerspectButton.set_label ("Persp");

  TubesButton.set_icon_widget (TubesImage);
  LinesButton.set_icon_widget (LinesImage);
  StraightButton.set_icon_widget(StraightImage);
  SmoothButton.set_icon_widget(SmoothImage);
  NoWrapButton.set_icon_widget(NoWrapImage);
  WrapButton.set_icon_widget(WrapImage);
  OrthoButton.set_icon_widget(OrthoImage);
  PerspectButton.set_icon_widget(PerspectImage);
  Tools.append (LinesButton);
  Tools.append (TubesButton);
  Tools.append (ToolSep1);
  Tools.append (StraightButton);
  Tools.append (SmoothButton);
  Tools.append (ToolSep2);
  Tools.append (NoWrapButton);
  Tools.append (WrapButton);
  Tools.append (ToolSep3);
  Tools.append (OrthoButton);
  Tools.append (PerspectButton);

  // Setup detail stuff
  DetailScale.set_adjustment(DetailAdjust);
  DetailScale.set_digits(0);
  DetailAdjust.set_step_increment(1.0);
  DetailAdjust.signal_value_changed().connect
    (sigc::mem_fun(*this, &VisualClass::OnDetailChange));
  DetailFrame.add(DetailScale);
  DetailScale.set_size_request(75,-1);

  // Setup iso stuff
  IsoScale.set_adjustment(IsoAdjust);
  IsoScale.set_digits(1);
  IsoAdjust.set_step_increment(0.1);
  IsoAdjust.signal_value_changed().connect
    (sigc::mem_fun(*this, &VisualClass::OnIsoChange));
  IsoFrame.add(IsoScale);
  IsoScale.set_size_request(75,-1);

  // Setup rho stuff
  RhoScale.set_adjustment(RhoAdjust);
  RhoScale.set_digits(0);
  RhoAdjust.set_step_increment(1.0);
  RhoAdjust.signal_value_changed().connect
    (sigc::mem_fun(*this, &VisualClass::OnRhoChange));
  RhoFrame.add(RhoScale);
  RhoScale.set_size_request(75,-1);

  // Setup the file chooser
  FileChooser.set_select_multiple(false);
  FileChooser.set_action(Gtk::FILE_CHOOSER_ACTION_OPEN);
  FileChooser.add_button(Gtk::Stock::CANCEL, Gtk::RESPONSE_CANCEL);
  FileChooser.add_button("OK", Gtk::RESPONSE_OK);


//   Gtk::Color black ("Black");
//   Glib::RefPtr<Pixmap> linesPM = 
//     Gdk::Pixmap::create_from_xpm (TubesButton.window, NULL, black, lines);


  Actions = Gtk::ActionGroup::create();
  Actions->add (Gtk::Action::create("MenuFile", "_File"));
  Actions->add (Gtk::Action::create("Open", "_Open"),
		sigc::mem_fun(*this, &VisualClass::OnOpen));
  Actions->add (Gtk::Action::create("Export", "_Export Image"),
		sigc::mem_fun(*this, &VisualClass::OnExport));
  Actions->add (Gtk::Action::create("ExportPOV", "_Export POV"),
		sigc::mem_fun(*this, &VisualClass::OnExportPOV));
  Actions->add (Gtk::Action::create("Quit", "_Quit"),
		sigc::mem_fun(*this, &VisualClass::Quit));
  Actions->add (Gtk::Action::create("MenuView", "View"));
  Actions->add (Gtk::Action::create("Reset", "Reset"),
		sigc::mem_fun(*this, &VisualClass::ResetView));
  Actions->add (Gtk::Action::create("Fullscreen", "Fullscreen"),
		sigc::mem_fun(*this, &Window::fullscreen));
  FullscreenToggle = Gtk::ToggleAction::create("Fullscreen", "Fullscreen mode",
					       "Display this window in fullscreen mode", false);
  Actions->add (FullscreenToggle,
		sigc::mem_fun(*this, &VisualClass::OnFullscreenToggle));

  Manager = Gtk::UIManager::create();
  Manager->insert_action_group(Actions);
  add_accel_group (Manager->get_accel_group());

  Glib::ustring ui_info =
    "<ui>"
    "  <menubar name='MenuBar'>"
    "    <menu action='MenuFile'>"
    "      <menuitem action='Open'/>"
    "      <menuitem action='Export'/>"
    "      <menuitem action='ExportPOV'/>"
    "      <separator/>"
    "      <menuitem action='Quit'/>"
    "    </menu>"
    "    <menu action='MenuView'>"
    "      <menuitem action='Reset'/>"
    "      <menuitem action='Fullscreen'/>"
    "    </menu>"
    "  </menubar>"
    "  <toolbar  name='ToolBar'>"
    "    <toolitem action='Open'/>"
    "    <toolitem action='Quit'/>"
    "  </toolbar>"
    "</ui>";
  
  Manager->add_ui_from_string (ui_info);
  m_VBox.pack_start (*Manager->get_widget("/MenuBar"), Gtk::PACK_SHRINK, 0);
  ToolBox.pack_start (Tools);
  ToolBox.pack_start (DetailFrame, Gtk::PACK_SHRINK, 0);
  ToolBox.pack_start (IsoFrame, Gtk::PACK_SHRINK, 0);
  ToolBox.pack_start (RhoFrame, Gtk::PACK_SHRINK, 0);
  m_VBox.pack_start(ToolBox, Gtk::PACK_SHRINK, 0);
  m_VBox.pack_start(PathVis);
  m_VBox.pack_start(FrameScale, Gtk::PACK_SHRINK,0);
  m_VBox.pack_start(m_ButtonQuit, Gtk::PACK_SHRINK, 0);


  // Show window.
  show_all();
  // HACK HACK HACK  
  IsoFrame.hide();
  processH2O = false;
}

VisualClass::~VisualClass()
{}
  
void VisualClass::OnOpen()
{
  int result = FileChooser.run();
  switch (result) {
    case (Gtk::RESPONSE_OK): {
      Read (FileChooser.get_filename());
      FrameChanged();
      break;
    }
    case (Gtk::RESPONSE_CANCEL): {
      cerr << "Cancel.\n";
      break;
    }
  }
  
  FileChooser.hide();
}

void VisualClass::OnExport()
{
  Export.SetupWidgets();
  Export.show_all();
  //  Export.Export ("frame.png");
}

void VisualClass::OnExportPOV()
{
  //  Export.ExportPOV("frame.pov");
}

void VisualClass::ResetView()
{
  PathVis.View.Reset();
  double maxDim = max(max(Box[0], Box[1]), Box[2]);
  PathVis.View.SetDistance (1.2*maxDim);
  PathVis.Invalidate();
}

void
VisualClass::OnFullscreenToggle()
{
  if (FullscreenToggle->get_active())
    fullscreen();
  else
    unfullscreen();
}

// bool VisualClass::on_delete_event()
// {
//   cerr << "delete event called.\n";
//   return true;
// }

void VisualClass::Quit()
{
  Gtk::Main::quit();
}

void VisualClass::FrameChanged()
{
  if (FileIsRead) {
    MakeFrame ((int)floor(FrameAdjust.get_value()));
    PathVis.Invalidate();
  }
}




void VisualClass::LineToggle()
{
  if (LinesButton.get_active())
    PathType = LINES;
  else
    PathType = TUBES;

  FrameChanged();
}


void VisualClass::SmoothToggle()
{
  Smooth = SmoothButton.get_active();
  FrameChanged();
}


void VisualClass::WrapToggle()
{
  Wrap = WrapButton.get_active();
  FrameChanged();
}

void VisualClass::PerspectiveToggle()
{
  bool persp = !OrthoButton.get_active();
  PathVis.View.SetPerspective(persp);
  //  PathVis.Invalidate();
  FrameChanged();
}


void VisualClass::MakePaths(int frame)
{
  int numPtcls  = PathArray.extent(1);
  int numSlices = PathArray.extent(2);
  Array<bool,1> used(numPtcls);
  used = false;

  for (int i=0; i<Paths.size(); i++)
    delete Paths[i];
  Paths.resize(0);

  vector<vector<int> > loopList;
  // Constuct list of permuting loops.  Ignore classical particles.
  for (int ptcl=0; ptcl<numPtcls; ptcl++) {
    if (!used(ptcl) && (PtclSpecies(ptcl).lambda!=0.0)) {
      vector<int> loop;
      int permPtcl = ptcl;
      bool haveOpen = false;
       do {
	 loop.push_back(permPtcl);
	 if (permPtcl == OpenPtcl(frame))
	   haveOpen = true;
	 used(permPtcl) = true;
	 permPtcl = PermArray(frame, permPtcl);
       } while (permPtcl != ptcl);

       if (haveOpen) {  // make sure open ptcl is las
	 if (loop.size()==2)
	   std::swap (loop[0], loop[1]);
// 	 while (loop[loop.size()-1] != OpenPtcl(frame) {
// 	   // cyclic permute
// 	 }
       }
       //>>>>>>> .r523
       loopList.push_back(loop);
    }
  }
  //  Paths.resize(loopList.size());
  for (int li=0; li<loopList.size(); li++) {
    vector<int> &loop = loopList[li];
    OnePath &path = (*new OnePath);
    path.Closed = true;
    path.Path.resize(numSlices*loop.size()+1);
    path.Color.resize(numSlices*loop.size()+1);
    int offset = 0;
    for (int pi=0; pi<loop.size(); pi++) {
      int ptcl = loop[pi];
      for (int slice=0; slice < numSlices; slice++) {
	path.Path[slice+offset][0] = PathArray(frame,ptcl,slice,0);
	path.Path[slice+offset][1] = PathArray(frame,ptcl,slice,1);
	path.Path[slice+offset][2] = PathArray(frame,ptcl,slice,2);
      }
      offset +=  numSlices;
      if (ptcl == OpenPtcl(frame)) 
	path.Closed = false;
    }
    // Close the path!!!
    if (path.Closed) {
      path.Path[path.Path.size()-1][0] = PathArray(frame,loop[0],0,0);
      path.Path[path.Path.size()-1][1] = PathArray(frame,loop[0],0,1);
      path.Path[path.Path.size()-1][2] = PathArray(frame,loop[0],0,2);
    }
    else {
      path.Path[path.Path.size()-1][0] = Tail(frame,0);
      path.Path[path.Path.size()-1][1] = Tail(frame,1);
      path.Path[path.Path.size()-1][2] = Tail(frame,2);
    }

    
    // First, put first slice in box
    Box.PutInBox(path.Path[0]);
    // Now, make sure the path doesn't have any discontinuities 
    // because having different period images.
    for (int slice=0; slice<path.Path.size()-1; slice++) 
      for (int dim=0; dim<3; dim++) {
	while ((path.Path[slice+1][dim] - path.Path[slice][dim]) >0.5*Box[dim])
	  path.Path[slice+1][dim] -= Box[dim];
	while ((path.Path[slice+1][dim] - path.Path[slice][dim])<-0.5*Box[dim])
	  path.Path[slice+1][dim] += Box[dim];
      }
    // Check to see if the path is closed or is a winding path
    int last = path.Path.size()-1;
    path.Closed = true;
    path.Winding = false;
    for (int dim=0; dim<3; dim++) {
      if (fabs(path.Path[last][dim] - path.Path[0][dim]) > 0.5*Box[dim]) {
	cerr << "Winding path.\n";
	path.Closed = false;
	path.Winding = true;
      }
    }
    // If we're not close, pop off the closing point
//   if (!path.Closed) 
//     path.Path.pop_back();

    Paths.push_back (&path);
  }
  // Do Fourier smoothing if desired
  if (Smooth)
    Smoother.SmoothPaths(Paths);

  // Wrap paths into box if desired
  if (Wrap)
    for (int i=0; i<3; i++)
      Box.PutPathsInBox (Paths);

}




void VisualClass::OnDetailChange()
{
  double val = DetailAdjust.get_value();
  Smoother.SetLevel(val);
  if (Smooth)
    FrameChanged();
}

void
VisualClass::OnIsoChange()
{
  FrameChanged();
}


void
VisualClass::OnRhoChange()
{
  FrameChanged();
}

void
VisualClass::SetViewportSize (int xsize, int ysize)
{
  PathVis.set_size_request(xsize, ysize);
  resize(10,10);
}




//////////
// Main.//
//////////

int main(int argc, char** argv)
{
  Gtk::Main kit(argc, argv);

  // Init gtkglextmm.
  Gtk::GL::init(argc, argv);
  glutInit (&argc, argv);

  list<ParamClass> optionList;
  optionList.push_back(ParamClass("small", false));
  optionList.push_back(ParamClass("water", false));
  CommandLineParserClass parser (optionList);
  bool success = parser.Parse (argc, argv);

  if (!success || parser.NumFiles() < 1) {
    cerr << "Usage:\n  pathvis++ [--small] [--water] myfile.h5\n";
    exit (1);
  }

  // Query OpenGL extension version.
  int major, minor;
  Gdk::GL::query_version(major, minor);
  std::cout << "OpenGL extension version - "
            << major << "." << minor << std::endl;

  // Instantiate and run the application.
  VisualClass visual;

  // John's addition to read in special flags
  // int index = 1;
  // if (argc == 3){
  //   cerr << "Read in flag " << argv[1] << endl;
  //   visual.SetFlag(argv[1]);
  //   index = 2;
  // }

  if (parser.Found("water")) {
    // John, add your stuff here:

  }
  if (parser.Found("small"))
    visual.SetViewportSize(800, 600);

  visual.Read (parser.GetFile(0));
  kit.run(visual);

  return 0;
}
