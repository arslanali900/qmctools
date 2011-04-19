AM_LDFLAGS = $(ALL_STATIC)

check_PROGRAMS = TestOptimizedBreakup TestLattice

TestOptimizedBreakup_SOURCES = TestOptimizedBreakup.cc
TestLattice_SOURCES          = TestLattice.cc
LDADD =  libbreakup.la ../Integration/libquad.la ../MatrixOps/libmatops.la ../MPI/libcomm.la $(MPI_LIBS) $(BLITZ_LIBS) $(LAPACK_LIBS) $(FLIBS)
AM_CXXFLAGS = $(BLITZ_CFLAGS) $(MPI_CFLAGS) $(MPI_DEF)

noinst_LTLIBRARIES = libbreakup.la
libbreakup_la_SOURCES = OptimizedBreakup2.cc OptimizedBreakup2.h LatticeClass.h
libbreakup_la_CXXFLAGS = $(BLITZ_CFLAGS)  $(MPI_CFLAGS) $(MPI_DEF) $(AM_CXXFLAGS)

includedir = $(prefix)/include/Common/Ewald/

include_HEADERS = OptimizedBreakup.h

EXTRA_DIST = Ewald.tex Simone_long_range_notes.tex
