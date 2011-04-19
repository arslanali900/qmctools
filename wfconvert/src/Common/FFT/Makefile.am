AM_LDFLAGS = $(ALL_STATIC)

check_PROGRAMS = TestFFT

TestFFT_SOURCES = TestFFT.cc
LDADD =  libfft.la $(FFTW3_LIBS) $(BLITZ_LIBS) 
AM_CXXFLAGS =  $(FFTW3_CFLAGS) $(BLITZ_CFLAGS) $(FFT_DEF)

noinst_LTLIBRARIES = libfft.la
libfft_la_SOURCES = FFT.cc
libfft_la_CXXFLAGS = $(BLITZ_CFLAGS)  $(FFTW3_CFLAGS) $(FFT_DEF) $(AM_CXXFLAGS) 

includedir = $(prefix)/include/Common/FFT/

include_HEADERS = FFT.h

