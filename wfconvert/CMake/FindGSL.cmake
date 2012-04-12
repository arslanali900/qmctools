#
# this module look for gsl (http://www.gnu.org/software/gsl) support
# it will define the following values
#
# GSL_INCLUDE_DIR = where gsl/gsl_version.h can be found
# GSL_LIBRARY     = the library to link against libgsl
# FOUND_GSL       = set to 1 if gsl is found
#
set(Libgsl gsl)
IF(QMC_BUILD_STATIC)
  set(Libgsl libgsl.a)
ENDIF(QMC_BUILD_STATIC)


IF(EXISTS ${PROJECT_CMAKE}/GslConfig.cmake)
  INCLUDE(${PROJECT_CMAKE}/GslConfig.cmake)
ENDIF(EXISTS ${PROJECT_CMAKE}/GslConfig.cmake)

IF(Gsl_INCLUDE_DIRS)

  FIND_PATH(GSL_INCLUDE_DIR gsl/gsl_version.h ${Gsl_INCLUDE_DIRS})
  FIND_LIBRARY(GSL_LIBRARY gsl ${Gsl_LIBRARY_DIRS})

ELSE(Gsl_INCLUDE_DIRS)

  FIND_LIBRARY(GSL_LIBRARY gsl ${GSL_HOME}/lib $ENV{GSL_HOME}/lib)
  FIND_LIBRARY(GSLBLAS_LIBRARY gslcblas ${GSL_HOME}/lib $ENV{GSL_HOME}/lib)
  FIND_PATH(GSL_INCLUDE_DIR gsl/gsl_version.h ${GSL_HOME}/include $ENV{GSL_HOME}/include)

ENDIF(Gsl_INCLUDE_DIRS)

IF(GSL_INCLUDE_DIR AND GSL_LIBRARY)
  SET(GSL_FOUND 1 CACHE BOOL "Found gsl library")
  SET(GSL_LIBRARY ${GSL_LIBRARY} ${GSLBLAS_LIBRARY})
ELSE(GSL_INCLUDE_DIR AND GSL_LIBRARY)
  SET(GSL_FOUND 0 CACHE BOOL "Not fount gsl library")
ENDIF(GSL_INCLUDE_DIR AND GSL_LIBRARY)

MARK_AS_ADVANCED(
  GSL_INCLUDE_DIR 
  GSL_LIBRARY 
  GSL_FOUND
  )
