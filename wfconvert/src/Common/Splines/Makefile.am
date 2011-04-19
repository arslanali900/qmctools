AM_LDFLAGS = $(ALL_STATIC)

check_PROGRAMS =  TestFastCubic TestFastMultiCubic TestBicubic TestComplexMultiTricubicSpline \
                  TestGrid TestMultiTricubicSpline TestMyTricubic TestPeriodic TestQuintic    \
                  TestCubicBspline TestTricubicBspline TestBlips TestNonuniform               \
                  TestNUBspline TestTricubicNUBspline TestBspline

TestFastCubic_SOURCES                  = TestFastCubic.cc
TestFastMultiCubic_SOURCES             = TestFastMultiCubic.cc
TestBicubic_SOURCES                    = TestBicubic.cc 
TestComplexMultiTricubicSpline_SOURCES = TestComplexMultiTricubicSpline.cc
TestGrid_SOURCES                       = TestGrid.cc
TestMultiTricubicSpline_SOURCES        = TestMultiTricubicSpline.cc
TestMyTricubic_SOURCES                 = TestMyTricubic.cc
TestPeriodic_SOURCES                   = TestPeriodic.cc
TestQuintic_SOURCES                    = TestQuintic.cc
TestCubicBspline_SOURCES               = TestCubicBspline.cc
TestTricubicBspline_SOURCES            = TestTricubicBspline.cc
TestBlips_SOURCES                      = TestBlips.cc
TestNonuniform_SOURCES                 = TestNonuniform.cc
TestNUBspline_SOURCES                  = TestNUBspline.cc
TestTricubicNUBspline_SOURCES          = TestTricubicNUBspline.cc
TestBspline_SOURCES                    = TestBspline.c

LDADD = ../IO/libIO.la ../PlaneWavePHDFT/libpw.la ../FFT/libfft.la libspline.la  $(HDF5_LIBS) $(FLIBS) $(BLITZ_LIBS) $(FFTW3_LIBS) $(LAPACK_LIBS)
AM_CXXFLAGS = -g $(BLITZ_CFLAGS) $(HDF5_CFLAGS) $(FFTW3_CFLAGS) $(FFT_DEF)

TestBspline_LDADD = libbspline.la 

noinst_LTLIBRARIES = libspline.la libbspline.la
libspline_la_SOURCES = BicubicSpline.h                \
                       BicubicSpline.cc               \
                       BsplineHelper.h                \
                       BoundaryCondition.h            \
                       CubicBspline.cc                \
                       CubicBspline.h                 \
                       CubicNUBspline.h               \
                       CubicSpline.h                  \
                       CubicSpline.cc                 \
                       FastCubicSpline.cc             \
                       FastCubicSpline.h              \
                       FastMultiCubicSpline.cc        \
                       FastMultiCubicSpline.h         \
                       Grid.h                         \
                       Grid.cc                        \
                       ComplexMultiTricubicSpline.h   \
                       ComplexMultiTricubicSpline.cc  \
                       LinearSpline.h                 \
                       MultiTricubicSpline.h          \
                       MultiTricubicSpline.cc         \
                       MyTricubicSpline.h             \
                       MyTricubicSpline.cc            \
                       NUBsplineBasis.h               \
                       PeriodicSpline.h               \
                       PeriodicSpline.cc              \
                       QuinticSpline.h                \
                       QuinticSpline.cc               \
                       QuinticSplines.f               \
                       TricubicSpline.h               \
                       TricubicBspline.cc             \
                       TricubicBspline.h              \
                       TricubicNUBspline.h            \
                       r3spline.f                     \
                       z3spline.f                     \
                       bcspeval.f                     \
                       bcspline.f                     \
                       evtricub.f                     \
                       herm3ev.f                      \
                       ibc_ck.f                       \
                       mktricubw.f                    \
                       splinck.f                      \
                       tcspeval.f                     \
                       tcspline.f                     \
                       v_spline.f                     \
                       zonfind.f                      
libbspline_la_SOURCES =                               \
                       bspline_base.h                 \
                       bspline_create.c               \
                       bspline_create.h               \
                       bspline_structs_sse.h          \
                       bspline_structs_std.h          \
                       bspline_data.c                 \
                       bspline_eval_std_s.h           \
                       bspline_eval_sse_s.h           \
                       bspline_eval_std_c.h           \
                       bspline_eval_sse_c.h           \
                       bspline_eval_std_d.h           \
                       bspline_eval_sse_d.h           \
                       bspline_eval_std_z.h           \
                       bspline_eval_sse_z.h           


includedir=$(prefix)/include/Common/Splines


#if WITH_F77
#     libsprng_la_SOURCES += common/fwrap_mpi.c
#endif
libspline_la_CXXFLAGS = $(BLITZ_CFLAGS) $(AM_CXXFLAGS) $(FFT_DEF)
libspline_la_LIBADD   =

#if HAVE_GMP
#     libsprng_la_SOURCES += pmlcg/pmlcg.h pmlcg/pmlcg.c pmlcg/info.h pmlcg/basic.h 
#endif


include_HEADERS = BicubicSpline.h ComplexMultiTricubicSpline.h CubicSpline.h Grid.h MultiTricubicSpline.h \
                  MyTricubicSpline.h PeriodicSpline.h QuinticSpline.h TricubicSpline.h FastCubicSpline.h  \
                  CubicBspline.h BsplineHelper.h BoundaryCondition.h NUBsplineBasis.h TricubicBspline.h   \
                  TricubicNUBspline.h bspline_base.h bspline_create.h bspline_eval_sse_c.h                \
                  bspline_eval_sse_d.h bspline_eval_sse_s.h bspline_eval_sse_z.h bspline_eval_std_c.h     \
                  bspline_eval_std_d.h bspline_eval_std_s.h bspline_eval_std_z.h bspline.h                \
                  bspline_structs_sse.h bspline_structs_std.h LinearSpline.h CubicNUBspline.h


EXTRA_DIST = CubicSplines.tex
