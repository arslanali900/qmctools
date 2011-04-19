#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
AC_PREREQ(2.59)
AC_INIT(wfconv, 0.5, kesler@ciw.edu, wfconvert)
AM_INIT_AUTOMAKE(wfconv, 0.5)
AC_CONFIG_SRCDIR([src/PlaneWaveClass.h])
AC_CONFIG_HEADER([src/config.h])

# Check which compilers we're using
AX_CXX_MAXOPT
AC_MSG_CHECKING([optimal compiler flags])
AC_MSG_RESULT([${CXXFLAGS}])

AC_CHECK_LIB([m], [sqrt])
AC_MSG_CHECKING([LIBS = $LIBS])
# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([stdlib.h sys/time.h unistd.h])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_PROG_F77
AC_F77_LIBRARY_LDFLAGS
AC_LANG_PUSH(C++)
AC_C_RESTRICT
AC_LANG_POP(C++)
AC_PROG_LIBTOOL
PKG_PROG_PKG_CONFIG

AC_CANONICAL_HOST
AC_MSG_CHECKING([host cpu = ])
AC_MSG_RESULT([${host_cpu}])
case "${host_cpu}" in
   "powerpc")
      CPPFLAGS="-DMAC $CPPFLAGS"
      CXXFLAGS="$CXX_FLAGS -O3 -mcpu=970 -mtune=970 -mpowerpc64 -mpowerpc-gpopt -ffast-math -g"
      FFLAGS="-O3 -mcpu=970 -mtune=970 -mpowerpc64 -mpowerpc-gpopt -ffast-math -g";;
   "i686")
      AC_MSG_CHECKING([ i686 cpu type])
      if grep -q XP /proc/cpuinfo; then
        AC_MSG_RESULT([Athlon XP])
        AX_CC_OPTION(TUNEOPT, [-mtune=athlon-xp], [-mtune=athlon-xp], [-mcpu=athlon-xp])
        AX_CC_OPTION(ARCHOPT, [-march=athlon-xp], [-march=athlon-xp], [-march=athlon])
	FFLAG="-O3 -march=athlon -ffast-math -msse -m3dnow -g"
      fi
      if grep -q "Pentium(R) M" /proc/cpuinfo; then
	AC_MSG_RESULT([Pentium M])
#        AX_CC_OPTION(TUNEOPT, [-mtune=pentium-m], [-mtune=pentium-m], [-mcpu=pentium-m])
#        AX_CC_OPTION(ARCHOPT, [-march=pentium-m], [-march=pentium-m], [-march=pentium-m])
#	FFLAGS="-O3 $TUNEOPT $ARCHOPT -ffast-math -g"
        AX_CC_OPTION(TUNEOPT, [-mtune=pentium-m], [-mtune=pentium-m], [-mcpu=pentium4])
	AX_CC_OPTION(ARCHOPT, [-march=pentium-m], [-march=pentium-m], [-mcpu=pentium4])
        AX_F77_OPTION(FARCHOPT, [-march=pentium4], [-march=pentium4], [ ])
	FFLAGS="-O3 $FARCHOPT -ffast-math -g"
      fi
      if grep -q "Intel(R) Xeon" /proc/cpuinfo; then
	AC_MSG_RESULT([Intel Xeon])
        AX_CC_OPTION(TUNEOPT, [-tpp7], [-tpp7], [-tpp7])
        AX_CC_OPTION(ARCHOPT, [-xW], [-xW], [-xW])
	FFLAGS="-O3 $TUNEOPT $ARCHOPT -g"
      fi
      CXXFLAGS="$CXXFLAGS -O3 $TUNEOPT $ARCHOPT -ffast-math -g"
      if grep -q sse /proc/cpuinfo ; then
        CXXFLAGS="$CXXFLAGS -msse"
      fi
      if grep -q sse2 /proc/cpuinfo ; then
        CXXFLAGS="$CXXFLAGS -msse2"
      fi
      if grep -q sse3 /proc/cpuinfo ; then
        CXXFLAGS="$CXXFLAGS -msse3"
      fi
      if grep -q 3dnow /proc/cpuinfo ; then
        CXXFLAGS="$CXXFLAGS -m3dnow"
      fi;;
      
    "ia64")
       AC_MSG_CHECKING([ Itanium cpu type])
       if grep -q "Itanium 2" /proc/cpuinfo; then
         AC_MSG_RESULT([Intel Itanium 2])
         AX_CC_OPTION(TUNEOPT, [-mcpu=itanium2], [-mcpu=itanium2], [-mcpu=itanium2])
       fi
       CXXFLAGS="$CXXFLAGS -O3 $TUNEOPT -g"
       FFLAGS="-O3 $TUNEOPT -g";;
esac


AC_ARG_ENABLE(openmp, AC_HELP_STRING([--enable-openmp], [enable OpenMP for thread-level parallelism]))
if test "$enable_openmp" = "yes"; then
   AX_OPENMP
fi

CXXFLAGS="$CXXFLAGS $OPENMP_CXXFLAGS $OPENMP_CFLAGS"
LDFLAGS="$LDFLAGS $OPENMP_CXXFLAGS"
#CXXFLAGS="$CXXFLAGS -pthread"

# Check if we want static-only binaries
AC_ARG_ENABLE(all-static, AC_HELP_STRING([--enable-all-static], [build static binaries (default=no)]))
if test "$enable_all_static" = "yes"; then
  ALL_STATIC="-all-static"
else
  ALL_STATIC=""
fi
AC_SUBST(ALL_STATIC)


# Compile CUDA versions
AC_ARG_ENABLE(cuda, AC_HELP_STRING([--enable-cuda], [compile CUDA routines]))
AM_CONDITIONAL([HAVE_CUDA], [test "$enable_cuda" = "yes"])
if test "x$enable_cuda" = "xyes"; then
   NVCC="nvcc"
   NVCCFLAGS="-Xcomplier -fopenmp -DNO_CUDA_MAIN -g -O3 -arch sm_13 -Drestrict=__restrict__"
   CUDA_LIBS="-L/usr/local/cuda/lib64 -lcudart"
   LIBS="$LIBS -L/usr/local/cuda/lib64 -lcudart"
   AC_DEFINE(HAVE_CUDA, [1], [Compile with CUDA extensions])
fi
########################
# Check for CUDA
########################
AC_ARG_WITH([cuda],
    [  --with-cuda=PATH        prefix where cuda is installed [default=auto]])
if test "x$with_cuda" != "x"; then
  CUDA_CFLAGS="-I$with_cuda/include"
  CUDA_LIBS="-L$with_cuda/lib64 -lcudart"
  NVCC="$with_cuda/bin/nvcc"
else 
   if test "x$enable_cuda" = "xyes"; then
     CUDA_CFLAGS="-I/usr/local/cuda/include"
     CUDA_LIBS="-L/usr/local/cuda/lib64 -lcudart"
     NVCC="nvcc"
   fi
fi 

echo $CUDA_LIBS

NVCCFLAGS="-Xcompiler -fopenmp -DNO_CUDA_MAIN -O3 -use_fast_math -g -arch sm_13 -Drestrict=__restrict__" 
AC_SUBST(CUDA_CFLAGS)
AC_SUBST(CUDA_LIBS)
AC_SUBST(NVCC)
AC_SUBST(NVCCFLAGS)



AC_SUBST(NVCC)
AC_SUBST(NVCCFLAGS)
AC_SUBST(CUDA_LIBS)


AC_ARG_ENABLE(profile,[  --enable-profile        instrument code with profiling information])
if test "x$enable_profile" = "xyes"; then
   CXXFLAGS="-pg $CXXFLAGS"
   FFLAGS="-pg $FFLAGS"
   LDFLAGS="-pg $LDFLAGS"
fi

AC_ARG_ENABLE(debug,[  --enable-debug          enable code for debugging checks])
if test "x$enable_debug" = "xyes"; then
   CXXFLAGS="-DDEBUG -DBZ_DEBUG $CXXFLAGS"
   FFLAGS="-DDEBUG $FFLAGS"
fi

# Find einspline library
PKG_CHECK_MODULES(EINSPLINE, einspline, einspline_ok=yes, einspline_ok=no)
if test "x$einspline_ok" != "xyes"; then
   AC_MSG_ERROR([Could not file einspline library.  Please download from einspline.sf.net.])
fi

# Find LIBXML
PKG_CHECK_MODULES(LIBXML, libxml-2.0, libxml_ok=yes, libxml_ok=no)
if test "x$libxml_ok" != "xyes"; then
   AC_MSG_ERROR([Could not file libxml library.  Please download from libxml2.])
fi

# Checks if we want MPI support
AC_ARG_ENABLE(mpi,    [  --enable-mpi            enable MPI support (default=no)])

# If so, use MPI version of libcommon.  Otherwise, use serial version
if test "x$enable_mpi" = "xyes"; then
  PKG_CHECK_MODULES(COMMON, mpicommon, common_ok=yes, common_ok=no)
  AC_MSG_CHECKING([Linking against parallel version of libmpicommon.])
else
  PKG_CHECK_MODULES(COMMON, common, common_ok=yes, common_ok=no)
  AC_MSG_CHECKING([Linking against serial version of libcommon.])
fi

if test "x$common_ok" != "xyes"; then
   AC_MSG_ERROR([Could not find common library.  Please install first.])
else
   AC_DEFINE([HAVE_COMMON],[1],[libcommon is available])
fi


AC_ARG_WITH(boost,    [  --with-boost=/path/to/my/boost	include boost libraries; needed for qmcpack])
if (test "x$with_boost" != "x"); then
	 	BOOST_LIB="-I$with_boost/include -L$with_boost/lib"
else # don't link with boost
		BOOST_LIB=" "
fi

AC_SUBST(BOOST_LIB)
AM_CONDITIONAL([PARALLEL], [test "$enable_mpi" = "yes"])
AC_F77_WRAPPERS



# Optionally compile visual
AC_ARG_ENABLE(vis, AC_HELP_STRING([--enable-vis], [compile visualization tool]))
AM_CONDITIONAL([WANT_VIS], [test "$enable_vis" = "yes"])
if test "x$enable_vis" = "xyes"; then
  PKG_CHECK_MODULES(GTKMM, gtkmm-2.4, gtkmm_ok=yes, gtkmm_ok=no)
  if test "x$gtkmm_ok" != "xyes"; then
     AC_MSG_ERROR([gtkmm libraries missing.  Please install.])
  fi
  PKG_CHECK_MODULES(GTKGLEXTMM, gtkglextmm-1.2, gtkglextmm_ok=yes, gtglextkmm_ok=no)
  if test "x$gtkglextmm_ok" != "xyes"; then
     AC_MSG_ERROR([gtkglextmm libraries missing.  Please install.])
  fi
  PKG_CHECK_MODULES(GTHREAD, gthread-2.0, gthread_ok=yes, gthread_ok=no)
  if test "x$gthread_ok" != "xyes"; then
     AC_MSG_ERROR([gthread libraries missing.  Please install.])
  fi
#  AC_CHECK_LIB([xvidcore],[xvid_decore])
  AC_LANG_SAVE
  AC_LANG_CPLUSPLUS
#  AC_CHECK_LIB([revel],[Revel_EncodeStart])
  AC_LANG_RESTORE
  AC_CHECK_LIB([glut],[glutSolidTeapot])
  AC_CHECK_LIB([gle],[gleSetNumSides])
fi



# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST
AC_C_INLINE
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_ERROR_AT_LINE
AC_FUNC_STRTOD
AC_CHECK_FUNCS([clock_gettime floor pow sqrt strtol])

AC_SUBST(PKGDATADEF)

AC_CONFIG_FILES([ Makefile                  \
		  src/Makefile              \
		  src/vis/Makefile])
AC_MSG_RESULT([ COMMON_LIBS  = $COMMON_LIBS ])
AC_OUTPUT
