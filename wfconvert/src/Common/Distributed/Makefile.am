AM_LDFLAGS = $(ALL_STATIC)

check_PROGRAMS = DistributedTest

DistributedTest_SOURCES = TestDistributedMat.cc
LDADD = ../MPI/libcomm.la libdistrib.la $(BLITZ_LIBS) $(MPI_LIBS)
AM_CXXFLAGS = $(BLITZ_CFLAGS) $(MPI_CFLAGS) $(MPI_DEF)

noinst_LTLIBRARIES = libdistrib.la
libdistrib_la_SOURCES = DistributedMat.cc
libdistrib_la_CXXFLAGS = $(BLITZ_CFLAGS)  $(MPI_CFLAGS) $(MPI_DEF) $(AM_CXXFLAGS)

includedir = $(prefix)/include/Common/Distributed

include_HEADERS = DistributedMat.h

