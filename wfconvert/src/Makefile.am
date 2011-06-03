AM_LDFLAGS = $(ALL_STATIC)

check_PROGRAMS = TestPBC

if WANT_VIS
   MAYBE_VIS = vis
endif

SUBDIRS = $(MAYBE_VIS) #common

if PARALLEL
  bin_PROGRAMS = wfconvert_mpi ppconvert TestSpline TestNU fitsk wfconv_mpi \
                 respline
else
  bin_PROGRAMS = wfconvert ppconvert TestSpline TestNU fitsk kMesh TestTile \
                 wfconv respline plot_lapw
endif

if HAVE_CUDA
   MAYBE_CUDA = PlaneWaveSum.cu AtomicOrbitalCuda.cu
endif


AM_CXXFLAGS = $(COMMON_CFLAGS) $(EINSPLINE_CFLAGS) $(LIBXML_CFLAGS) $(CUDA_CFLAGS) -g

TestPBC_SOURCES   =    \
    TestPBC.cc         \
    LatticeClass.h     \
    PBC.h 

plot_lapw_SOURCES =    \
    CubicSpline.cc     \
    CubicSpline.h      \
    LAPWClass.cc       \
    LAPWClass.h        \
    MuffinTin.cc       \
    MuffinTin.h        

respline_SOURCES  =    \
    config.h           \
    ParseCommand.cc    \
    ParseCommand.h     \
    Respline.cc        \
    Respline.h

wfconvert_SOURCES =    \
    ExportCasino.cc    \
    Localize.cc        \
    ParserClass.cc     \
    ParserClass.h      \
    ParseCommand.h     \
    ParseCommand.cc    \
    PlaneWaveClass.h   \
    PlaneWaveClass.cc  

wfconv_SOURCES =       \
    AtomicOrbital.h    \
    AtomicOrbital.cc   \
    APWClass.h         \
    APWClass.cc        \
    CellClass.h        \
    CellClass.cc       \
    CenterClass.h      \
    CenterClass.cc     \
    CubicSpline.h      \
    CubicSpline.cc     \
    Gaussian.h         \
    Gaussian.cc        \
    LatticeClass.h     \
    LocalizeOrbs.cc    \
    MuffinTin.h        \
    MuffinTin.cc       \
    OrbitalClass.h     \
    OrbitalClass.cc    \
    OrbitalSetClass.h  \
    OrbitalSetClass.cc \
    ParserClass.cc     \
    ParserClass.h      \
    ParseCommand.h     \
    ParseCommand.cc    \
    ReadABINIT.h       \
    ReadABINIT.cc      \
    Read_ESHDF.cc      \
    Read_FPMD.cc       \
    ReadLAPW.cc        \
    TileMapClass.h     \
    TileMapClass.cc    \
    WFConv.cc          \
    WriteCASINO.cc     \
    Write_qmcPACK.cc   \
    Write_ESHDF.cc     \
    $(MAYBE_CUDA)
#    LocalizeOrbs.cc    

wfconv_mpi_SOURCES =   $(wfconv_SOURCES)



wfconvert_mpi_SOURCES = \
    ExportCasino.cc    \
    Localize.cc        \
    OrbitalClass.h     \
    OrbitalClass.h     \
    ParserClass.cc     \
    ParserClass.h      \
    ParseCommand.h     \
    ParseCommand.cc    \
    PlaneWaveClass.h   \
    PlaneWaveClass.cc  

ppconvert_SOURCES =    \
    CubicSpline.cc     \
    CubicSpline.h      \
    NLPPClass.cc       \
    NLPPClass.h        \
    ParseCommand.h     \
    ParseCommand.cc    \
    ParserClass.cc     \
    ParserClass.h      \
    XMLWriterClass2.cc \
    XMLWriterClass2.h  

fitsk_SOURCES =        \
    FitSk.cc           \
    FitSk.h            \
    ParseCommand.h     \
    ParseCommand.cc    \
    ParserClass.cc     \
    ParserClass.h      

TestTile_SOURCES =     \
    TestTile.cc        \
    ParseCommand.h     \
    ParseCommand.cc    

kMesh_SOURCES =        \
    kMesh.cc

TestSpline_SOURCES =   \
    CubicSpline.cc     \
    CubicSpline.h      \
    TestSpline.cc    

TestNU_SOURCES =       \
    TryNonuniform.cc


#INCLUDES = -I$(QMCPACK_HOME)/src -I$(QMCPACK_ROOT)/src 

respline_LDADD      = $(COMMON_LIBS) $(EINSPLINE_LIBS)
wfconvert_LDADD     = $(COMMON_LIBS) 
wfconvert_mpi_LDADD = $(COMMON_LIBS) 
wfconv_LDADD        = $(COMMON_LIBS) $(EINSPLINE_LIBS) $(LIBXML_LIBS) $(CUDA_LIBS)
wfconv_mpi_LDADD    = $(COMMON_LIBS) $(EINSPLINE_LIBS) $(LIBXML_LIBS) $(CUDA_LIBS)
ppconvert_LDADD     = $(COMMON_LIBS) 
fitsk_LDADD         = $(COMMON_LIBS)
TestNU_LDADD        = $(COMMON_LIBS)
TestTile_LDADD      = $(COMMON_LIBS) $(EINSPLINE_LIBS)
TestPBC_LDADD       = $(COMMON_LIBS)
plot_lapw_LDADD     = $(COMMON_LIBS) $(EINSPLINE_LIBS)

.cu.o:
	$(NVCC) -c $(NVCCFLAGS) $<

.cu.lo:
	$(top_builddir)/cudalt.py $@ $(NVCC) -c $(NVCCFLAGS) $<
