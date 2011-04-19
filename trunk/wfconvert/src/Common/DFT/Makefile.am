AM_LDFLAGS = $(ALL_STATIC)

noinst_LTLIBRARIES = libdft.la
libdft_la_SOURCES = ExCorr.f        \
                    Functionals.cc  \
                    Functionals.h   

includedir=$(prefix)/include/Common/DFT

libdft_la_CXXFLAGS = $(AM_CXXFLAGS) $(BLITZ_CFLAGS)
libdft_la_LIBADD   =

include_HEADERS = Functionals.h

