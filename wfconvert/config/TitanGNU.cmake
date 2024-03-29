SET(CMAKE_SYSTEM_PROCESSOR "XK6")
#2011-12-06

set(CMAKE_C_COMPILER  /opt/cray/xt-asyncpe/5.05/bin/cc)
set(CMAKE_CXX_COMPILER  /opt/cray/xt-asyncpe/5.05/bin/CC)
set(CMAKE_Fortran_COMPILER /opt/cray/xt-asyncpe/5.05/bin/ftn)

set(GNU_OPTS "-DADD_ -DINLINE_ALL=inline")
set(GNU_FLAGS "-fopenmp -O3 -Drestrict=__restrict__ -finline-limit=1000 -fstrict-aliasing -funroll-all-loops -Wno-deprecated ")
set(XT_FLAGS "-march=bdver1 -msse3 -D_CRAYMPI")
#set(XT_FLAGS "-msse3 -D_CRAYMPI")
set(CMAKE_CXX_FLAGS "${XT_FLAGS} ${GNU_FLAGS} -ftemplate-depth-60 ${GNU_OPTS}")
set(CMAKE_C_FLAGS "${XT_FLAGS} ${GNU_FLAGS} -std=c99")
set(CMAKE_Fortran_FLAGS "-O3 -march=bdver1 -funroll-all-loops -fno-f2c")
set(CMAKE_Fortran_FLAGS_RELEASE ${CMAKE_Fortran_FLAGS})
set(CMAKE_Fortran_FLAGS_DEBUG  "-march=bdver1 -fopenmp  -msse3 -fno-f2c -O0 -g")

SET(QMC_BUILD_STATIC 1)
SET(ENABLE_OPENMP 1)
SET(HAVE_MPI 1)
SET(HAVE_SSE 1)
SET(HAVE_SSE2 1)
SET(HAVE_SSE3 1)
SET(HAVE_SSSE3 1)
SET(USE_PREFETCH 1)
SET(PREFETCH_AHEAD 12)

set(CMAKE_FIND_ROOT_PATH_MODE_PROGRAM NEVER)
set(CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY)
set(CMAKE_FIND_ROOT_PATH_MODE_INCLUDE ONLY)
set(CMAKE_SHARED_LINKER_FLAGS "")

FOREACH(type SHARED_LIBRARY SHARED_MODULE EXE)
  SET(CMAKE_${type}_LINK_STATIC_C_FLAGS "-Wl,-Bstatic")
  SET(CMAKE_${type}_LINK_DYNAMIC_C_FLAGS "-static")
  SET(CMAKE_${type}_LINK_STATIC_CXX_FLAGS "-Wl,-Bstatic")
  SET(CMAKE_${type}_LINK_DYNAMIC_CXX_FLAGS "-static")
ENDFOREACH(type)

set(CMAKE_FIND_ROOT_PATH
  /opt/cray/hdf5/1.8.6/gnu/46
  /opt/fftw/3.3.0.0/interlagos
  /sw/xk6/boost/1.44.0/cle4.0_gnu4.5.3
  /ccs/proj/mat034/jnkim/xk6/gnu45/libxml2
  /sw/xk6/gsl/1.15/cle4.0_gcc4.6.2/
)

set(EINSPLINE_SSE_BUG 1)
#set(HAVE_EINSPLINE 1)
#set(HAVE_EINSPLINE_EXT 0)
