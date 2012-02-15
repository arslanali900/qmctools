PROJECT(einspline)

CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

SET (LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib CACHE PATH "Single output directory for building all libraries.")
SET (EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin CACHE PATH "Single output directory for building all executables.")

SET(HAVE_CUDA 0)

INCLUDE(${CMAKE_ROOT}/Modules/CheckFunctionExists.cmake)
INCLUDE(${CMAKE_ROOT}/Modules/CheckIncludeFile.cmake)
INCLUDE(${CMAKE_ROOT}/Modules/CheckTypeSize.cmake)

CHECK_FUNCTION_EXISTS(sincos HAVE_SINCOS)
CHECK_FUNCTION_EXISTS(round HAVE_STD_ROUND)
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

INCLUDE(${CMAKE_ROOT}/Modules/TestCXXAcceptsFlag.cmake)
INCLUDE(${CMAKE_ROOT}/Modules/CheckCCompilerFlag.cmake)

SET (PROJECT_CMAKE ${PROJECT_SOURCE_DIR}/CMake)

IF(CMAKE_TOOLCHAIN_FILE)
  MESSAGE(STATUS "Using ${CMAKE_TOOLCHAIN_FILE} toolchain ")
ELSE(CMAKE_TOOLCHAIN_FILE)
  SET(FOUND_CXXENV 0)
  SET(ENABLE_OPENMP 0)
  SET(HAVE_MPI 0)
  SET(INTEL_COMPILER 0)
  SET(HAVE_SSE 0)
  SET(HAVE_SSE2 0)
  SET(HAVE_SSE3 0)
  SET(HAVE_SSSE3 0)
  SET(USE_PREFETCH 0)
  SET(PREFETCH_AHEAD 12)
  SET(HAVE_C_VARARRAYS 1)

  #######################################################
  ## Try to identify CPU identity
  #######################################################
  #SET(CPU_IDENTITY "generic")
  #INCLUDE(${PROJECT_CMAKE}/CheckProcessorID.cmake)

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

ENDIF(CMAKE_TOOLCHAIN_FILE)


IF(ENABLE_CUDA)
  #if(CMAKE_COMPILER_IS_GNUCXX)
  #  SET(CUDA_NVCC_FLAGS "-arch;sm_13;")
  #else(CMAKE_COMPILER_IS_GNUCXX)
  #  set(CUDA_PROPAGATE_HOST_FLAGS "OFF")
  #  SET(CUDA_NVCC_FLAGS "-arch;sm_13;-Drestrict=__restrict__")
  #endif(CMAKE_COMPILER_IS_GNUCXX)

  FIND_PACKAGE(CUDA REQUIRED)
  if(CUDA_FOUND)
    set(HAVE_CUDA 1)
  endif()

  if(CMAKE_COMPILER_IS_GNUCXX)
    SET(CUDA_NVCC_FLAGS "-arch;sm_13;")
  else(CMAKE_COMPILER_IS_GNUCXX)
    set(CUDA_PROPAGATE_HOST_FLAGS "OFF")
    SET(CUDA_NVCC_FLAGS "-arch;sm_13;-Drestrict=__restrict__")
  endif(CMAKE_COMPILER_IS_GNUCXX)

  MESSAGE(STATUS "Setting CUDA FLAGS="${CUDA_NVCC_FLAGS})

ENDIF(ENABLE_CUDA)

CONFIGURE_FILE(${PROJECT_SOURCE_DIR}/src/config.h.cmake.in
  ${PROJECT_BINARY_DIR}/src/config.h)

INCLUDE_DIRECTORIES(
  ${PROJECT_SOURCE_DIR}/src
  ${PROJECT_BINARY_DIR}/src
  )

SUBDIRS(src)
