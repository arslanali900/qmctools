PROJECT(qmctools)

CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
enable_language (Fortran)

#####################################################
SET(CMAKE_FALLBACK FALSE)
if(${CMAKE_MINOR_VERSION} VERSION_LESS "8")
  MESSAGE(WARNING ">>>> Recommend using cmake version 2.8.0 and higher. Using old cmake scripts. <<<< ")
  SET(CMAKE_FALLBACK TRUE)
endif(${CMAKE_MINOR_VERSION} VERSION_LESS "8")

######################################################################
# Version 
######################################################################
SET(QMCTOOLS_VERSION_MAJOR 0)
SET(QMCTOOLS_VERSION_MINOR 0)
SET(QMCTOOLS_VERSION_PATCH 0)
SET(QMCTOOLS_VERSION
  "${QMCTOOLS_VERSION_MAJOR}.${TOOLSLUS_VERSION_MINOR}.${QMCTOOLS_VERSION_PATH}")

######################################################################
#define the paths for library and executable for out-of-soruce compilation
######################################################################
SET (LIBRARY_OUTPUT_PATH ${qmctools_BINARY_DIR}/lib CACHE PATH "Single output directory for building all libraries.")
SET (EXECUTABLE_OUTPUT_PATH ${qmctools_BINARY_DIR}/bin CACHE PATH "Single output directory for building all executables.")

######################################################################
# build include and library paths required for the projects
# For XYZ library, set an environment as
# export XYZ_HOME='root-of-the-library'
# XYZ_HOME/include
# XYZ_HOME/lib
# IF XYZ library is found, XYZ_FOUND is set to true
#
# find_package is used with cmake 2.8 and higher
######################################################################
SET(QMC_INCLUDE_PATHS 
  $ENV{HDF5_HOME}/include 
  $ENV{HDF5_DIR}/include 
  $ENV{FFTW_HOME}/include
  $ENV{LIBXML2_HOME}/include
  $ENV{EINSPLINE_HOME}/include
  )
SET(QMC_LIBRARY_PATHS 
  $ENV{HDF5_HOME}/lib 
  $ENV{HDF5_DIR}/lib 
  $ENV{FFTW_HOME}/lib
  $ENV{LIBXML2_HOME}/lib
  $ENV{EINSPLINE_HOME}/lib
  )

if(APPLE)
  set(QMC_INCLUDE_PATHS ${QMC_INCLUDE_PATHS} /sw/include)
  set(QMC_LIBRARY_PATHS ${QMC_LIBRARY_PATHS} /sw/lib)
endif(APPLE)

######################################################################
# enable MPI and OPNEMP  if possible
######################################################################
SET(QMC_MPI 1 CACHE BOOL "Enable/disable MPI")
SET(QMC_OMP 1 CACHE BOOL "Enable/disable OpenMP")
SET(QMC_BITS 64 CACHE INTEGER "Select OS bit")
SET(HAVE_CUDA 0 CACHE BOOL "Build with GPU support through CUDA")
SET(PRINT_DEBUG 0 CACHE BOOL "Enable/disable debug printing")
######################################################################
# check environment variables 
######################################################################
IF($ENV{QMC_MPI})
  SET(QMC_MPI $ENV{QMC_MPI}) 
  MESSAGE("-- Found QMC_MPI=${QMC_MPI}")
ENDIF($ENV{QMC_MPI})
IF($ENV{QMC_OMP})
  SET(QMC_OMP $ENV{QMC_OMP}) 
  MESSAGE("-- Found QMC_OMP=${QMC_OMP}")
ENDIF($ENV{QMC_OMP})
#IF($ENV{QMC_BITS})
#  SET(QMC_BITS $ENV{QMC_BITS})
#  MESSAGE("-- Found QMC_BITS=${QMC_BITS}")
#ENDIF($ENV{QMC_BITS})
IF($ENV{QMC_COMPLEX})
  SET(QMC_COMPLEX $ENV{QMC_COMPLEX})
  MESSAGE("-- Found QMC_COMPLEX=${QMC_COMPLEX}")
ENDIF($ENV{QMC_COMPLEX})
#IF($ENV{QMC_CUDA})
#  SET(QMC_CUDA $ENV{QMC_CUDA} CACHE BOOL "Build with GPU support through CUDA")
#  MESSAGE("-- Found QMC_CUDA=${QMC_CUDA}")
#ENDIF($ENV{QMC_CUDA})

######################################################################
# Standard test
######################################################################
SET(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/CMake" ${CMAKE_MODULE_PATH})
INCLUDE(${CMAKE_ROOT}/Modules/CheckFunctionExists.cmake)
INCLUDE(${CMAKE_ROOT}/Modules/CheckIncludeFile.cmake)
INCLUDE(${CMAKE_ROOT}/Modules/CheckTypeSize.cmake)
CHECK_FUNCTION_EXISTS(sincos HAVE_SINCOS)
CHECK_FUNCTION_EXISTS(round HAVE_STD_ROUND)

######################################################################
# add macros to compile einspline if necessary
######################################################################
CHECK_FUNCTION_EXISTS(floor HAVE_FLOOR)
CHECK_FUNCTION_EXISTS(posix_memalign HAVE_POSIX_MEMALIGN)
CHECK_FUNCTION_EXISTS(pow HAVE_POW)
CHECK_FUNCTION_EXISTS(sqrt HAVE_SQRT)

CHECK_INCLUDE_FILE(dlfcn.h HAVE_DLFCN_H)
CHECK_INCLUDE_FILE(inttypes.h HAVE_INTTYPES_H)
CHECK_INCLUDE_FILE(memory.h HAVE_MEMORY_H)
CHECK_INCLUDE_FILE(pmmintrin.h HAVE_PMMINTRIN_H)
CHECK_INCLUDE_FILE(emmintrin.h HAVE_EMMINTRIN_H)
CHECK_INCLUDE_FILE(sys/stat.h HAVE_SYS_STAT_H)
CHECK_INCLUDE_FILE(sys/time.h HAVE_SYS_TIME_H)
CHECK_INCLUDE_FILE(sys/types.h HAVE_SYS_TYPES_H)
CHECK_INCLUDE_FILE(unistd.h HAVE_UNISTD_H)

######################################################################
# Build options
# QMC_BUILD_STATIC build static/dynamic  libraries
# BUILD_QMCTOOLS   build utility programs
# BUILD_SANDBOX    build test programs
# MPIP_PROFILE     profile mpi performance
######################################################################
SET(QMC_BUILD_STATIC 0 CACHE BOOL "Build static libraries and binary")


######################################################################
#list of options to determine libraries. First, try to use libraries available
#FindXYZ.cmake will set the values if successful
######################################################################
SET(ENABLE_PHDF5 0 CACHE BOOL "Enable phdf5 for output")
SET(HAVE_LIBXML2 1)
SET(HAVE_LIBFFTW 0)
SET(HAVE_LIBXMLPP 0)
SET(HAVE_LIBHDF5 0)
SET(HAVE_LIBSPRNG 0)
SET(HAVE_LIBBLITZ 0)
SET(HAVE_LIBGSL 0)
SET(MAC_VECLIB 0)
SET(ENABLE_GTK_GUI 0)
SET(HAVE_EINSPLINE 0)
SET(HAVE_EINSPLINE_EXT 0)
set(ENABLE_TAU_PROFILE 0 CACHE BOOL "Enable tau for profiling")
SET(HAVE_COMMON 1)
set(TAU_HOME $ENV{TAU_HOME})

######################################################
# Directory where customize cmake files reside
######################################################
SET (PROJECT_CMAKE ${qmctools_SOURCE_DIR}/CMake)

######################################################
# Get the current SVN version
######################################################
SET(QMCTOOLS_SVN_LAST_CHANGED "2011-04-18")
FIND_PACKAGE(Subversion) 
IF(Subversion_FOUND) 
  Subversion_WC_INFO(${PROJECT_SOURCE_DIR} Project) 
  MESSAGE("-- Current revision : ${Project_WC_REVISION} modified on ${Project_WC_LAST_CHANGED_DATE}") 
  SET(QMCPLUSPLUS_BRANCH ${Project_WC_REVISION})
  SET(QMCPLUSPLUS_LAST_CHANGED_DATE ${Project_WC_LAST_CHANGED_DATE})
ENDIF(Subversion_FOUND)

######################################################
#COMPILER choose one of the cmake files to customize the compiler options
#If nothing is chosen, default settings by cmake will  be used.
#If the automatic detection does not work, comment out everything 
#upto COMPILER MANUAL SELECTION and use one of the customized cmake file.
######################################################
INCLUDE(${CMAKE_ROOT}/Modules/TestCXXAcceptsFlag.cmake)
INCLUDE(${CMAKE_ROOT}/Modules/CheckCCompilerFlag.cmake)

######################################################
# DO NOT MODIFY BELOW
# If you like to customize your build, create a toolchain as instructed
######################################################
######################################################
# now set the compiler flags and the macros in src/config.h
######################################################
IF(CMAKE_TOOLCHAIN_FILE)
  MESSAGE(STATUS "-- Using ${CMAKE_TOOLCHAIN_FILE} toolchain ")
ELSE(CMAKE_TOOLCHAIN_FILE)
  #fun begins
  SET(FOUND_CXXENV 0)
  SET(ENABLE_OPENMP 0)
  SET(HAVE_MPI 0)
  SET(INTEL_COMPILER 0)
  SET(HAVE_MKL 0)
  SET(HAVE_MKL_VML 0)
  SET(HAVE_ACML 0)
  SET(HAVE_SSE 0)
  SET(HAVE_SSE2 0)
  SET(HAVE_SSE3 0)
  SET(HAVE_SSSE3 0)
  SET(USE_PREFETCH 0)
  SET(PREFETCH_AHEAD 12)
  SET(HAVE_C_VARARRAYS 1)

  ######################################################
  # Try to identify CPU identity
  ######################################################
  SET(CPU_IDENTITY "generic")
  INCLUDE(${PROJECT_CMAKE}/CheckProcessorID.cmake)

  #------------------------------------
  # On Jaguar use CrayXT.cmake
  #------------------------------------
  IF($ENV{HOST} MATCHES "jaguar")
    MESSAGE("  Working on jaguar. Using CMake/CrayXT.cmake")	 
    SET(CPU_IDENTITY "barcelona")
    INCLUDE(${PROJECT_CMAKE}/CrayXT.cmake)
    SET(FOUND_CXXENV 1)
    SET(CMAKE_SKIP_RPATH TRUE)
  ENDIF($ENV{HOST} MATCHES "jaguar")

  #------------------------------------
  # Check if using IBM compilers
  #------------------------------------
  IF($ENV{CXX} MATCHES "xlC")
    INCLUDE(${PROJECT_CMAKE}/IBMCompilers.cmake)
    SET(FOUND_CXXENV 1)
  ENDIF($ENV{CXX} MATCHES "xlC")

  #------------------------------------
  # Check if using Intel compilers
  #------------------------------------
  IF($ENV{CXX} MATCHES "icpc")
    INCLUDE(${PROJECT_CMAKE}/IntelCompilers.cmake)
    SET(FOUND_CXXENV 1)
  ENDIF($ENV{CXX} MATCHES "icpc")

  #------------------------------------
  # other compilers, e.g., mpicxx 
  # most likely *unix with gnu or intel compilers
  # using "-restrict" option to find out if intel compilers are backend.
  #------------------------------------
  IF(NOT FOUND_CXXENV)
    IF(CMAKE_COMPILER_IS_GNUCXX)
      INCLUDE(${PROJECT_CMAKE}/GNUCompilers.cmake)
    ELSE(CMAKE_COMPILER_IS_GNUCXX)
      SET(CMAKE_TRY_INTEL_CXX_FLAGS "-restrict")
      CHECK_CXX_ACCEPTS_FLAG(${CMAKE_TRY_INTEL_CXX_FLAGS} INTEL_CXX_FLAGS)
      IF(INTEL_CXX_FLAGS)
        INCLUDE(${PROJECT_CMAKE}/IntelCompilers.cmake)
        SET(FOUND_CXXENV 1)
      ENDIF(INTEL_CXX_FLAGS)
    ENDIF(CMAKE_COMPILER_IS_GNUCXX)
  ENDIF(NOT FOUND_CXXENV)

  #-------------------------------------------------------------------
  #  check MPI installation
  #-------------------------------------------------------------------
  IF(QMC_MPI)
    ## mpi compilers
    if($ENV{CXX} MATCHES "mp")
      SET(HAVE_MPI 1)
      SET(HAVE_OOMPI 1)
      SET(MPI_FOUND TRUE)
    else($ENV{CXX} MATCHES "mp")
      include(FindMPI)
      IF(MPI_FOUND)
        MESSAGE(STATUS "Found mpi.h and libraries. Check for working mpi compiler ")
        FILE(WRITE ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/testmpi.cxx 
          "#include <mpi.h>"
          "int main(int argc, char** argv){MPI_Init(&argc,&argv); MPI_Finalize();return 0;}\n")
        TRY_COMPILE(MPI_WORKS ${CMAKE_BINARY_DIR} 
          ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/testmpi.cxx
          OUTPUT_VARIABLE OUTPUT)
        IF(MPI_WORKS)
          MESSAGE(STATUS "MPI is enabled")
          SET(HAVE_MPI 1)
          SET(HAVE_OOMPI 1)
          LINK_LIBRARIES(${MPI_LIBRARY})
        ELSE(MPI_WORKS)
          MESSAGE(STATUS "MPI is disabled")
          SET(HAVE_MPI 0)
          SET(HAVE_OOMPI 0)
        ENDIF(MPI_WORKS)
      ENDIF(MPI_FOUND)
    endif($ENV{CXX} MATCHES "mp")
  ENDIF(QMC_MPI)

  #-------------------------------------------------------------------
  # check lapack/blas
  #-------------------------------------------------------------------
  INCLUDE(${PROJECT_CMAKE}/FindLapack.cmake)


ENDIF(CMAKE_TOOLCHAIN_FILE)

MESSAGE("-- CXX compiler flags: ${CMAKE_CXX_FLAGS}")
MESSAGE("-- CC  compiler flags: ${CMAKE_C_FLAGS}")

SET(HAVE_OOMPI ${HAVE_MPI})
INCLUDE(${CMAKE_ROOT}/Modules/FindThreads.cmake)

####################################################################
#First check the required libraries. Abort if these are not found.
#Changed the fall-back versions to use similar conventions to 2.8.0 and higher
####################################################################
# use customized cmake file for boost so that we can deal with multiple
# installataions
INCLUDE(${PROJECT_CMAKE}/FindBoost.cmake)
INCLUDE(${PROJECT_CMAKE}/FindLibxml2.cmake)
INCLUDE(${PROJECT_CMAKE}/FindHDF5.cmake)
INCLUDE(${PROJECT_CMAKE}/FindZLIB.cmake)
INCLUDE(${PROJECT_CMAKE}/FindSZLIB.cmake)
INCLUDE(${PROJECT_CMAKE}/FindGSL.cmake)
#INCLUDE(${PROJECT_CMAKE}/FindBlitzpp.cmake)
#IF(BLITZ_FOUND)
#INCLUDE_DIRECTORIES(${BLITZ_INCLUDE_DIR} ${BLITZ_PLATFORM_INCLUDE_DIR})
#ENDIF(BLITZ_FOUND)

IF(LIBXML2_FOUND)
  INCLUDE_DIRECTORIES(${LIBXML2_INCLUDE_DIR})
  LINK_LIBRARIES(${LIBXML2_LIBRARIES})
  SET(HAVE_LIBXML2 1)
ENDIF(LIBXML2_FOUND)

IF(HDF5_FOUND)
  SET(HAVE_LIBHDF5 1)
  INCLUDE_DIRECTORIES(${HDF5_INCLUDE_DIR})
  LINK_LIBRARIES(${HDF5_LIBRARIES})
  IF(ZLIB_FOUND)
    LINK_LIBRARIES(${ZLIB_LIBRARIES})
  ENDIF(ZLIB_FOUND)
  IF(SZLIB_FOUND)
    LINK_LIBRARIES(${SZLIB_LIBRARIES})
  ENDIF(SZLIB_FOUND)
ELSE(HDF5_FOUND)
  MESSAGE(FATAL_ERROR "Require hdf5 1.6.4 or higher. Set HDF5_HOME")
ENDIF(HDF5_FOUND)

IF(GSL_FOUND)
  set(HAVE_LIBGSL 1)
  INCLUDE_DIRECTORIES(${GSL_INCLUDE_DIR})
  LINK_LIBRARIES(${GSL_LIBRARY})
ENDIF(GSL_FOUND)
 
#IF(Boost_FOUND)
#  SET(HAVE_LIBBOOST 1)
#  MESSAGE(STATUS "Setting Boost_INCLUDE_DIR="${Boost_INCLUDE_DIR})
#  INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIR})
#ELSE()
#  MESSAGE(FATAL_ERROR "Require boost 1.33.x or higher. Set BOOST_HOME")
#ENDIF()

IF(HAVE_CUDA)
  SET(CUDA_NVCC_FLAGS "-arch;sm_13")
  FIND_PACKAGE(CUDA REQUIRED)
  SET(CUDA_NVCC_FLAGS "-arch;sm_13")
ENDIF(HAVE_CUDA)

######################################################
# Add optional macros
######################################################
INCLUDE(${PROJECT_CMAKE}/FindFFTW.cmake)
IF(FFTW_FOUND)
  SET(HAVE_LIBFFTW 1)
  INCLUDE_DIRECTORIES(${FFTW_INCLUDE_DIR})
  LINK_LIBRARIES(${FFTW_LIBRARIES})
ENDIF(FFTW_FOUND)

######################################################################
# Manage Einspline library
# - check if libeinspline can be found in $EINSPLINE_HOME/include and
# $EINSPLINE_HOME/lib
# - if not found, build it from the source directory in $EINSPLINE_HOME
# - if everything fails, do not use it
######################################################################
INCLUDE(${PROJECT_CMAKE}/FindEinspline.cmake)
if(EINSPLINE_FOUND)
  set(HAVE_EINSPLINE 1)
  IF(HAVE_EINSPLINE_EXT)
    INCLUDE_DIRECTORIES(${EINSPLINE_INCLUDE_DIR})
    LINK_LIBRARIES(${EINSPLINE_LIBRARIES})  
  ENDIF(HAVE_EINSPLINE_EXT)
  MESSAGE(STATUS "EINSPLINE_INCLUDE_DIR="${EINSPLINE_INCLUDE_DIR})
endif()

#######################################################################
## Manage GTK+ and GTKMM packages for GUI
## search headers and libraries
##   1) -lXYZ-${XYZ_VERSION} 
##   2) -lXYZ
## example: -lgtk-${GTK_VERSION} | -lgtk
#######################################################################
#SET(GTK_VERSION 2.0)
#SET(PANGO_VERSION 1.0)
#
#SET(GTKMM_VERSION 2.4)
#SET(PANGOMM_VERSION 1.4)
#SET(SIGCPP_VERSION 2.0)
#IF(ENABLE_GTK_GUI)
#  INCLUDE(${PROJECT_CMAKE}/FindGTK.cmake)
#  #THIS SHOULD GO TO the directory where GUI is built
#  IF(GTKMM_FOUND)
#   INCLUDE_DIRECTORIES(${GTKMM_INCLUDE_DIR})
#   LINK_LIBRARIES(${GTKMM_LIBRARIES})
#  ENDIF(GTKMM_FOUND)
#ENDIF(ENABLE_GTK_GUI)

###################################################
# Done searching libraries. Start building qmcpack
###################################################

MESSAGE("-- Ready to build qmcpack")

INCLUDE_DIRECTORIES(
  ${qmctools_SOURCE_DIR}/src
  ${qmctools_BINARY_DIR}/src
  )

####################################################################
# add definitions
#  -DH5_USE_16_API for HDF5
#  -DMPICH_SKIP_MPICXX : disable C++ binding of MPICH
####################################################################
ADD_DEFINITIONS(-DHAVE_CONFIG_H -DH5_USE_16_API)
if(HAVE_MPI)
  ADD_DEFINITIONS(-DMPICH_SKIP_MPICXX)
endif()

#ADD_DEFINITIONS(-DDISABLE_WS_CELL)
CONFIGURE_FILE(${PROJECT_SOURCE_DIR}/src/config.h.cmake.in
               ${PROJECT_BINARY_DIR}/src/Common/config.h)

SUBDIRS(src/Common src)
