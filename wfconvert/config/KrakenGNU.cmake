# tool chain for Jaguar, Cray XT5  using PrgEnv-gnu
# Fri Feb 20 15:50:10 EST 2009 by Jeongnim Kim
# module list
# Currently Loaded Modulefiles:
# 1) modules/3.1.6                           13) Base-opts/2.1.41HD
# 2) DefApps                                 14) subversion/1.5.0
# 3) torque/2.3.2-snap.200807092141          15) gcc/4.2.0.quadcore
# 4) moab/5.2.4                              16) fftw/3.1.1
# 5) xtpe-quadcore                           17) xt-libsci/10.3.1
# 6) MySQL/5.0.45                            18) xt-mpt/3.1.0
# 7) xt-service/2.1.41HD                     19) xt-pe/2.1.41HD
# 8) xt-libc/2.1.41HD                        20) xt-asyncpe/2.0
# 9) xt-os/2.1.41HD                          21) PrgEnv-gnu/2.1.41HD
# 10) xt-boot/2.1.41HD                        22) cmake/2.6.1
# 11) xt-lustre-ss/2.1.UP00_ORNL.nic52_1.6.5  23) acml/4.1.0
# 12) xtpe-target-cnl                         24) hdf5/1.6.8

SET(CMAKE_SYSTEM_PROCESSOR "XT5")
SET_PROPERTY(GLOBAL PROPERTY TARGET_SUPPORTS_SHARED_LIBS FALSE)

set(CMAKE_C_COMPILER  /opt/cray/xt-asyncpe/4.0/bin/cc)
set(CMAKE_CXX_COMPILER  /opt/cray/xt-asyncpe/4.0/bin/CC)
set(CMAKE_Fortran_COMPILER /opt/cray/xt-asyncpe/4.0/bin/ftn)

set(GNU_OPTS "-DADD_ -DINLINE_ALL=inline")
set(GNU_FLAGS "-fopenmp -O3 -ftemplate-depth-60 -Drestrict=__restrict__  -finline-limit=1000 -fstrict-aliasing -funroll-all-loops -Wno-deprecated ")
set(XT_FLAGS "-march=amdfam10 -msse3 -D_CRAYMPI")
set(CMAKE_CXX_FLAGS "${XT_FLAGS} ${GNU_FLAGS} -ftemplate-depth-60 ${GNU_OPTS}")
set(CMAKE_C_FLAGS "${XT_FLAGS} ${GNU_FLAGS}")
set(CMAKE_Fortran_FLAGS "-O3 -march=amdfam10 -funroll-all-loops -fno-f2c")
set(CMAKE_Fortran_FLAGS_RELEASE ${CMAKE_Fortran_FLAGS})
set(CMAKE_Fortran_FLAGS_DEBUG  "-march=amdfam10 -fopenmp  -msse3 -fno-f2c -O0 -g")

FOREACH(type SHARED_LIBRARY SHARED_MODULE EXE)
  SET(CMAKE_${type}_LINK_STATIC_C_FLAGS "-Wl,-Bstatic")
  SET(CMAKE_${type}_LINK_DYNAMIC_C_FLAGS "-Wl,-Bstatic")
  SET(CMAKE_${type}_LINK_STATIC_CXX_FLAGS "-Wl,-Bstatic")
  SET(CMAKE_${type}_LINK_DYNAMIC_CXX_FLAGS "-Wl,-Bstatic")
ENDFOREACH(type)

set(ACML_HOME /opt/acml/4.2.0/gfortran64)

set(CMAKE_FIND_ROOT_PATH
      /opt/fftw/3.2.2.1
      /sw/xt/hdf5/1.8.5/cnl2.2_gnu4.4.3
      /sw/xt/szip/2.1/sles10.1_gnu4.4.3
      /nics/a/proj/qmc/boost_1_38_0
      /nics/a/proj/qmc/kraken/einspline
      /sw/xt/gsl/1.14/cnl2.2_gnu4.4.3
     )

SET(ENABLE_OPENMP 1)
SET(HAVE_MPI 1)
SET(HAVE_SSE 1)
SET(HAVE_SSE2 1)
SET(HAVE_SSE3 1)
SET(HAVE_SSSE3 1)
SET(USE_PREFETCH 1)
SET(PREFETCH_AHEAD 12)
SET(ACML_LIBRARIES ${ACML_HOME}/lib/libacml.a ${ACML_HOME}/lib/libacml_mv.a)
link_libraries(${ACML_LIBRARIES})

