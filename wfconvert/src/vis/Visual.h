#ifndef VISUAL_H
#define VISUAL_H

#include "../Common.h"
#include "PathVis.h"
#include <Common/IO/IO.h>
#include <gtkmm/adjustment.h>
#include "OnePath.h"
#include "BoxClass.h"
#include "SmoothClass.h"
#include "Export.h"
// #include "ExportVideo.h"

/// This species class stores info about the species for the path
/// visualization.  Not the same as the pimc++ version.
class SpeciesClass
{
public:
  double lambda;
  string Name;
  int NumParticles;
  int FirstParticle, LastParticle;
  SpeciesClass () : FirstParticle(0)
  { 

  }
};

typedef enum {LINES, TUBES} PathTypeType;


using namespace IO;

/// VisualClass is the main GUI interaction for the path visualization
/// code.  It allows the user to rotate around and zoom in and out of
/// a box containing paths and classical particles dumped from
/// pimc++.  In includes features such as a frame slider, line or tube
/// rendering of paths, and optional fourier smoothing with detail
/// control.  
class VisualClass : public Gtk::Window
{
protected:
  friend class ExportClass;
  //////////
  // Data //
  //////////
  /// This stores the raw paths from the pimc++ output file.  The 4
  /// dimension are (frame, ptcl, slice, dim)
  Array<double,4> PathArray;

  ///////////////
  // Node data //
  ///////////////
  bool HaveANodeData, HaveBNodeData;
  bool HaveWarpPos;
  /// This stores raw node data from the pimc++ output file.  The 4
  /// dimensions are (x, y, z)
  void ReadFrameData (int frame);
  IOSectionClass Infile;
  IOVarBase *ANodeVar, *BNodeVar, *RhoVar;
  Array<int,1> NodePtcl, NodeSlice;
  Array<double,3> ANodeData, BNodeData, RhoData;
  Array<double,2> WarpPos;
  LinearGrid Xgrid, Ygrid, Zgrid;

  /// These are used only for open loops.
  Array<int,1>    OpenPtcl;
  Array<double,2>   Tail;

  //////////////////////////////////////
  // items to process water molecules //
  //////////////////////////////////////
  bool processH2O;

  /// This stores the global permutations for the paths.
  /// (frame, ptcl).  The permutation acts between the last and first
  /// slice.  
  Array<int,2> PermArray;
  Array<SpeciesClass,1> Species;
  vector<OnePath*> Paths;
  BoxClass Box;
  SmoothClass Smoother;
  void MakePaths(int frame);

  // member widgets:
  Gtk::VBox m_VBox;
  Gtk::Button m_ButtonQuit;
  Gtk::HScale FrameScale;
  Gtk::Adjustment FrameAdjust;
  Gtk::Toolbar Tools;
  Gtk::HBox ToolBox;
  Gtk::RadioToolButton LinesButton, TubesButton, 
    SmoothButton, StraightButton,
    WrapButton, NoWrapButton,
    OrthoButton, PerspectButton;
  void FrameChanged();
  PathTypeType PathType; 
  void LineToggle(), WrapToggle(), SmoothToggle(), PerspectiveToggle();
  Gtk::Image TubesImage, LinesImage, StraightImage, SmoothImage,
    NoWrapImage, WrapImage, OrthoImage, PerspectImage;
  Gtk::SeparatorToolItem ToolSep1, ToolSep2, ToolSep3;

  // Detail control
  Gtk::Frame DetailFrame;
  Gtk::HScale DetailScale;
  Gtk::Adjustment DetailAdjust;

  // Isosurface control
  Gtk::Frame IsoFrame;
  Gtk::HScale IsoScale;
  Gtk::Adjustment IsoAdjust;
  // Density isosurface control
  Gtk::Frame RhoFrame;
  Gtk::HScale RhoScale;
  Gtk::Adjustment RhoAdjust;
  double MaxRho, MinRho;

  Glib::RefPtr<Gtk::ToggleAction> FullscreenToggle;
  void OnFullscreenToggle();
  Glib::RefPtr<Gtk::ActionGroup> Actions;
  Glib::RefPtr<Gtk::UIManager> Manager;

  void OnOpen();
  void OnExport();
  void OnExportPOV();
  //  bool on_delete_event();
  void Quit();
  void ResetView();
  void SetFullscreen ();

  void PutInBox();
  void OnDetailChange();
  void OnIsoChange();
  void OnRhoChange();

  Gtk::FileChooserDialog FileChooser;

  bool Wrap, Smooth;
  ExportClass Export;
  string FindFullPath(string filename);
  //  ExportVideoClass ExportVideo;
  bool FileIsRead;
public:
  PathVisClass PathVis;

  inline SpeciesClass& PtclSpecies (int ptcl);
  void Read(string fileName);
  void MakeFrame (int frame, bool offScreen=false);
  void SetFlag(string newFlag);
  void SetViewportSize(int xsize, int ysize);

  VisualClass();
  virtual ~VisualClass();

};


inline SpeciesClass& VisualClass::PtclSpecies (int ptcl)
{
  int si = 0;
  while (si < Species.size() && (Species(si).LastParticle < ptcl))
    si++;
  return Species(si);
}

#endif
