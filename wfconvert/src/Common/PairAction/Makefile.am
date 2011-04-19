AM_LDFLAGS = $(ALL_STATIC)

check_PROGRAMS = TestGradients TestHDF

TestGradients_SOURCES = TestGradients.cc

TestHDF_SOURCES = testHDF.cc

LDADD = libpairact.la  ../Fitting/libfitting.la ../MatrixOps/libmatops.la ../PH/libPH.la ../Splines/libspline.la  ../IO/libIO.la ../Integration/libquad.la ../Distributed/libdistrib.la ../SpecialFunctions/libspecfun.la  ../MPI/libcomm.la $(HDF5_LIBS) $(BLITZ_LIBS) $(GSL_LIBS) $(FLIBS) $(MPI_LIBS) $(LAPACK_LIBS)
AM_CXXFLAGS = $(MPI_DEF) $(BLITZ_CFLAGS) $(HDF5_CFLAGS) $(GSL_CFLAGS) $(MPI_CFLAGS)
noinst_LTLIBRARIES = libpairact.la
libpairact_la_SOURCES = DavidPAClass.cc    \
                        DavidPAClass.h     \
                        PADipoleFit.cc     \
                        PADipoleFit.h      \
                        PAFit.h            \
                        PAFitBase.h        \
                        PAMonopoleFit.cc   \
                        PAMonopoleFit.h    \
                        PATripoleFit.cc    \
                        PATripoleFit.h     \
                        PAclassicalFit.cc  \
                        PAclassicalFit.h   \
                        PAcoulombBCFit.cc  \
                        PAcoulombBCFit.h   \
                        PAtricubicFit.cc   \
                        PAtricubicFit.h    \
                        PAtricubicFit2.cc  \
                        PAtricubicFit2.h   \
                        PAzeroFit.cc       \
                        PAzeroFit.h        \
                        Particle.h         \
                        FreeParticles.h    \
                        FreeParticles.cc   \
                        U_l.h              \
                        U_l.cc             \
                        Rho.h              \
                        Rho.cc 


includedir=$(prefix)/include/Common/PairAction


#if WITH_F77
#     libsprng_la_SOURCES += common/fwrap_mpi.c
#endif
libpairact_la_CXXFLAGS =  $(MPI_CFLAGS) $(BLITZ_CFLAGS) $(HDF5_CFLAGS) $(GSL_CFLAGS) $(MPI_DEF)
libpairact_la_LIBADD   =

#if HAVE_GMP
#     libsprng_la_SOURCES += pmlcg/pmlcg.h pmlcg/pmlcg.c pmlcg/info.h pmlcg/basic.h 
#endif


include_HEADERS = DavidPAClass.h    PAMonopoleFit.h   PAtricubicFit2.h \
                  PADipoleFit.h     PATripoleFit.h    PAzeroFit.h      \
                  PAFit.h           PAclassicalFit.h  Particle.h       \
                  PAFitBase.h       PAcoulombBCFit.h  PAtricubicFit.h  \
                  U_l.h             Rho.h             FreeParticles.h


