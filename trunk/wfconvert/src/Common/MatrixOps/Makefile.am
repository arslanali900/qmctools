AM_LDFLAGS = $(ALL_STATIC)

check_PROGRAMS = Test

Test_SOURCES = Test.cc 

LDADD = libmatops.la $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS)
AM_CXXFLAGS = $(LAPACK_LIBS) $(BLAS_LIBS) $(BLITZ_CFLAGS)

noinst_LTLIBRARIES = libmatops.la
libmatops_la_SOURCES = MatrixOps.h MatrixOps.cc

includedir=$(prefix)/include/Common/MatrixOps

libmatops_la_CXXFLAGS = $(BLITZ_CFLAGS)

include_HEADERS = MatrixOps.h

