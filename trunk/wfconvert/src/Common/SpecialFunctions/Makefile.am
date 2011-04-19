AM_LDFLAGS = $(ALL_STATIC)

check_PROGRAMS = TestHermite TestLegendre TestPoly TestPolySet

TestHermite_SOURCES = TestHermite.cc

TestLegendre_SOURCES = TestLegendre.cc 

TestPoly_SOURCES = TestPoly.cc

TestPolySet_SOURCES = TestPolySet.cc 

LDADD = ../Integration/libquad.la libspecfun.la $(BLITZ_LIBS) $(GSL_LIBS) $(FLIBS)
AM_CXXFLAGS = $(BLITZ_CFLAGS) $(GSL_CFLAGS) 

noinst_LTLIBRARIES = libspecfun.la
libspecfun_la_SOURCES = HermitePoly.h       \
                        SpecialFunctions.h  \
                        SpecialFunctions.cc \
                        Polynomial.h        \
                        PolynomialSet.h     \
                        PolynomialSet.cc    \
                        LegendrePoly.h      

libspecfun_la_CXXFLAGS = $(GSL_CFLAGS) $(BLITZ_CFLAGS)

includedir=$(prefix)/include/Common/SpecialFunctions

include_HEADERS = HermitePoly.h LegendrePoly.h   SpecialFunctions.h \
                  Polynomial.h  PolynomialSet.h 
