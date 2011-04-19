AM_LDFLAGS = $(ALL_STATIC)
bin_PROGRAMS = wfvis resvis

LDADD = $(GTKMM_LIBS) $(GTKGLEXTMM_LIBS) $(COMMON_LIBS) $(GTS_LIBS) $(GTHREAD_LIBS) $(CAIROMM_LIBS)
AM_CXXFLAGS = $(GTKMM_CFLAGS) $(GTKGLEXTMM_CFLAGS) $(COMMON_CFLAGS)\
              $(GTS_CFLAGS) $(GTHREAD_CFLAGS) $(CAIROMM_CFLAGS) \
              -DPKG_DATA_DIR=\"$(pkgdatadir)/\"

wfvis_SOURCES =         \
    BoxClass.cc         \
    BoxClass.h          \
    BoxObject.cc        \
    BoxObject.h         \
    ColorMap.h          \
    ColorMap.cc         \
    CoordObject.cc      \
    CoordObject.h       \
    CylinderObject.cc   \
    CylinderObject.h    \
    DiskObject.cc       \
    DiskObject.h        \
    ElementData.cc      \
    ElementData.h       \
    Isosurface.h        \
    Isosurface.cc       \
    GLObject.cc         \
    GLObject.h          \
    WFVis.cc            \
    WFVis.h             \
    ParseCommand.cc     \
    ParseCommand.h      \
    PathVis.cc          \
    PathVis.h           \
    PlaneObject.cc      \
    PlaneObject.h       \
    PointsObject.h      \
    PointsObject.cc     \
    SphereObject.cc     \
    SphereObject.h      \
    ViewClass.cc        \
    WFExport.cc         \
    WFExport.h          \
    trackball.h         \
    trackball.c    

resvis_SOURCES =         \
    BoxClass.cc         \
    BoxClass.h          \
    BoxObject.cc        \
    BoxObject.h         \
    ColorMap.h          \
    ColorMap.cc         \
    CoordObject.cc      \
    CoordObject.h       \
    CylinderObject.cc   \
    CylinderObject.h    \
    DiskObject.cc       \
    DiskObject.h        \
    ElementData.cc      \
    ElementData.h       \
    Isosurface.h        \
    Isosurface.cc       \
    GLObject.cc         \
    GLObject.h          \
    ResVis.cc            \
    ResVis.h             \
    ParseCommand.cc     \
    ParseCommand.h      \
    PathVis.cc          \
    PathVis.h           \
    PlaneObject.cc      \
    PlaneObject.h       \
    PointsObject.h      \
    PointsObject.cc     \
    SphereObject.cc     \
    SphereObject.h      \
    ViewClass.cc        \
    ResExport.cc         \
    ResExport.h          \
    trackball.h         \
    trackball.c    


dist_pkgdata_DATA = tubes.png lines.png straight.png smooth.png nowrap2.png \
                    wrap.png orthographic.png perspective.png player_play.png \
                    player_pause.png player_rev.png clipping.png              \
                    isoButton.png colormaps.in
