AM_LDFLAGS = $(ALL_STATIC)

check_PROGRAMS = TestkSpacePH MakeGeneralPot TestNLPP

TestkSpacePH_SOURCES = TestkSpacePH.cc
MakeGeneralPot_SOURCES = MakeGeneralPot.cc
TestNLPP_SOURCES       = TestNLPP.cc

LDADD = libPH.la ../IO/libIO.la   ../MatrixOps/libmatops.la ../Fitting/libfitting.la ../Integration/libquad.la ../Splines/libspline.la $(HDF5_LIBS) $(BLITZ_LIBS) $(G2C_LIBS) $(GSL_LIBS) $(LAPACK_LIBS) $(FLIBS)

AM_CXXFLAGS = $(HDF5_CFLAGS) $(BLITZ_CFLAGS) $(GSL_CFLAGS)

noinst_LTLIBRARIES = libPH.la
libPH_la_SOURCES = CoreTransform.h    \
                   CoreTransform.cc   \
                   CoulombPot.h       \
                   CoulombPot.cc      \
                   GeneralPot.h       \
                   GeneralPot.cc      \
                   HeAzizPot.h        \
                   HeAzizPot.cc       \
                   NLPP.h             \
                   NLPP.cc            \
                   Potential.h        \
                   Potential.cc       \
                   PotentialBase.h    \
                   PotentialBase.cc   \
                   QuinticPH.h        \
                   QuinticPH.cc       \
                   ScreenedPot.h      \
                   ScreenedPot.cc     \
                   SplinePot.h        \
                   SplinePot.cc       \
                   ToppHopfield.h     \
                   ToppHopfield.cc    \
                   kSpacePH.h         \
                   kSpacePH.cc      


includedir=$(prefix)/include/Common/PH


#if WITH_F77
#     libsprng_la_SOURCES += common/fwrap_mpi.c
#endif
#libPH_la_CXXFLAGS = $(BLITZ_CFLAGS) $(HDF5_CFLAGS)
libPH_la_LIBADD   =

#if HAVE_GMP
#     libsprng_la_SOURCES += pmlcg/pmlcg.h pmlcg/pmlcg.c pmlcg/info.h pmlcg/basic.h 
#endif


include_HEADERS = CoreTransform.h CoulombPot.h HeAzizPot.h Potential.h  \
                  PotentialBase.h QuinticPH.h ScreenedPot.h SplinePot.h \
	          kSpacePH.h ToppHopfield.h GeneralPot.h NLPP.h

