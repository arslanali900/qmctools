AM_LDFLAGS = $(ALL_STATIC)

check_PROGRAMS = TestIO
#check_PROGRAMS = TestHDF5 TestIO TestXML GrabRho

LDADD = $(HDF5_LIBS) $(BLITZ_LIBS) $(LIBS)
AM_CXXFLAGS = $(HDF5_CFLAGS) $(BLITZ_CFLAGS)

#TestXML_SOURCES = TestXML.cc
#TestXML_LDADD = libciao.la $(XML2_LIBS) $(HDF5_LIBS) $(BLITZ_LIBS) $(LIBS)

#TestHDF5_SOURCES = TestHDF5.cc
#TestHDF5_LDADD = libciao.la $(XML2_LIBS) $(HDF5_LIBS) $(BLITZ_LIBS) $(LIBS)

TestIO_SOURCES = TestIO.cc 
TestIO_LDADD = libIO.la $(HDF5_LIBS) $(BLITZ_LIBS) $(LIBS)

#GrabRho_SOURCES = GrabRho.cc 
#GrabRho_LDADD = libIO.la $(HDF5_LIBS) $(BLITZ_LIBS) $(LIBS)

includedir=$(prefix)/include/Common/IO

#noinst_LTLIBRARIES = libciao.la libIO.la
noinst_LTLIBRARIES  = libIO.la
#libciao_la_SOURCES = FileExpand.h        \
#                     InputOutput.h       \
#                     InputOutput.cc      \
#                     InputOutputBase.h   \
#                     InputOutputASCII.h  \
#                     InputOutputASCII.cc \
#                     InputOutputHDF5.h   \
#                     InputOutputHDF5.cc  \
#                     InputOutputXML.h    \
#                     InputOutputXML.cc   

libIO_la_SOURCES =   IO.cc               \
                     IO.h                \
                     IOASCII.cc          \
                     IOASCII.h           \
                     IOBase.h            \
                     IOHDF5.cc           \
                     IOHDF5.h            \
                     IOVar.h             \
                     IOVarASCII.h        \
                     IOVarBase.h         \
                     IOVarHDF5.cc        \
                     IOVarHDF5.h          


if WANT_PYTHON_WRAPPERS
  # Python wrappers
  pyIOexecdir = $(pyexecdir)

  pyIOexec_PYTHON = IO.py

  pyIOexec_LTLIBRARIES = IOSection.la
  IOSection_la_LDFLAGS  = -module -avoid-version -export-symbols-regex init_IO
  IOSection_la_LIBADD   = libIO.la $(HDF5_LIBS) $(BLITZ_LIBS) $(LIBS)
  IOSection_la_CPPFLAGS = $(PYTHON_INCLUDES) $(IO_CFLAGS) $(NUMPY_INC)
  IOSection_la_SOURCES  = IONumPyWrapper.cc
endif


#if WITH_F77
#     libsprng_la_SOURCES += common/fwrap_mpi.c
#endif
libciao_la_CXXFLAGS = $(BLITZ_CFLAGS) $(HDF5_CFLAGS)

#if HAVE_GMP
#     libsprng_la_SOURCES += pmlcg/pmlcg.h pmlcg/pmlcg.c pmlcg/info.h pmlcg/basic.h 
#endif


#include_HEADERS = InputOutput.h     InputOutputASCII.h InputOutputBase.h \
#                  InputOutputHDF5.h InputOutputXML.h   FileExpand.h      \
#                  IO.h              IOASCII.h          IOBase.h          \
#                  IOHDF5.h          IOVar.h            IOVarASCII.h      \
#                  IOVarBase.h       IOVarHDF5.h

include_HEADERS = FileExpand.h      IO.h              IOASCII.h          \
                  IOBase.h          IOHDF5.h          IOVar.h            \
                  IOVarASCII.h      IOVarBase.h       IOVarHDF5.h

