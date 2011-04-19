AM_LDFLAGS = $(ALL_STATIC)

check_PROGRAMS = TestRandom

TestRandom_SOURCES = TestRandom.cc

AM_CXXFLAGS = $(MPI_CFLAGS) $(MPI_DEF) $(BLITZ_CFLAGS) $(SPRNG2_CFLAGS)
LDADD    = ../MPI/libcomm.la $(SPRNG2_LIBS) $(MPI_LIBS)

includedir=$(prefix)/include/Common/Random

include_HEADERS = Random.h
