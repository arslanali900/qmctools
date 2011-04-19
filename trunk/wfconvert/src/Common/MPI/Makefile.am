AM_LDFLAGS = $(ALL_STATIC)

check_PROGRAMS = TestComm

TestComm_SOURCES = TestComm.cc

AM_CXXFLAGS = $(MPI_CFLAGS) $(MPI_DEF)  $(BLITZ_CFLAGS) $(SPRNG2_CFLAGS)
LDADD = libcomm.la $(MPI_LIBS) $(SPRNG2_LIBS)

noinst_LTLIBRARIES = libcomm.la
libcomm_la_SOURCES = Communication.cc

includedir=$(prefix)/include/Common/MPI

#if WITH_F77
#     libsprng_la_SOURCES += common/fwrap_mpi.c
#endif
libcomm_la_CXXFLAGS = $(MPI_CFLAGS) $(MPI_DEF) $(BLITZ_CFLAGS)  $(AM_CXXFLAGS)
libcomm_la_LIBADD   =

#if HAVE_GMP
#     libsprng_la_SOURCES += pmlcg/pmlcg.h pmlcg/pmlcg.c pmlcg/info.h pmlcg/basic.h 
#endif


include_HEADERS = Communication.h

