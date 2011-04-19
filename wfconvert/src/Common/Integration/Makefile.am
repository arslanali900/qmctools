AM_LDFLAGS = $(ALL_STATIC)

check_PROGRAMS = TestIntegrate TestHermite

TestIntegrate_SOURCES = TestIntegrate.cc
TestHermite_SOURCES = TestHermite.cc 

LDADD = libquad.la ../IO/libIO.la $(BLITZ_LIBS) $(HDF5_LIBS)
AM_CXXFLAGS = $(BLITZ_CFLAGS) $(HDF5_CFLAGS)

noinst_LTLIBRARIES = libquad.la
libquad_la_SOURCES = GKIntegration.h    \
                     GKIntegration.cc   \
                     HermiteQuad.h      \
                     HermiteQuad.cc     \
                     Integrate.h        \
                     Integrate.cc       \
                     RungeKutta.h       \
                     Standard.h         

libquad_la_CXXFLAGS = $(HDF5_CLFAGS) $(AM_CXXFLAGS)

includedir=$(prefix)/include/Common/Integration

include_HEADERS = GKIntegration.h HermiteQuad.h Integrate.h RungeKutta.h Standard.h
