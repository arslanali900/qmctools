check_PROGRAMS = TestDFTAtom TestRadialWF TestNLPPAtom

TestDFTAtom_SOURCES  = TestDFTAtom.cc 
TestNLPPAtom_SOURCES = TestNLPPAtom.cc 
TestRadialWF_SOURCES = TestRadialWF.cc

LDADD   = libatom.la                 ../IO/libIO.la           \
          ../DFT/libdft.la           ../PH/libPH.la             \
          ../Splines/libspline.la    ../Integration/libquad.la  \
          ../Fitting/libfitting.la   ../MatrixOps/libmatops.la  \
          $(BLITZ_LIBS) $(GSL_LIBS) $(FLIBS) $(HDF5_LIBS) $(LAPACK_LIBS)
AM_CXXFLAGS = $(BLITZ_CFLAGS) $(HDF5_CFLAGS)

includedir=$(prefix)/include/Common/Atom

noinst_LTLIBRARIES = libatom.la
libatom_la_SOURCES = AtomBase.h          \
                     DFTAtom.h           \
                     DFTAtom.cc          \
                     NewAtom.h           \
                     RadialWF.h          \
                     NewAtom.h           \
                     RadialWF.cc

libatom_la_CXXFLAGS = $(BLITZ_CFLAGS) $(AM_CXXFLAGS)
#libatom_la_LIBADD = $(GMP_LIBS) 



include_HEADERS = AtomBase.h DFTAtom.h RadialWF.h NewAtom.h


