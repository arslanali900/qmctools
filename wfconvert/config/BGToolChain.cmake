# the name of the target operating system
#SET(CMAKE_SYSTEM_NAME BlueGeneP)
SET(BGP 1 CACHE BOOL "On BlueGeneP")
SET(Linux 0)

SET(ENABLE_OPENMP 1)
SET(HAVE_MPI 1)
SET(HAVE_LIBESSL 1)

# set the compiler
#set(CMAKE_C_COMPILER  /opt/ibmcmp/vacpp/bg/9.0/bin/bgxlc_r)
#set(CMAKE_CXX_COMPILER  /opt/ibmcmp/vacpp/bg/9.0/bin/bgxlC_r)
set(CMAKE_C_COMPILER  /soft/apps/darshan/bin/default/mpixlc_r)
set(CMAKE_CXX_COMPILER  /soft/apps/darshan/bin/default/mpixlcxx_r)
set(CMAKE_Fortran_COMPILER /soft/apps/darshan/bin/default/mpixlf95_r)


# set the search path for the environment coming with the compiler
# and a directory where you can install your own compiled software
set(CMAKE_FIND_ROOT_PATH
/home/jnkim/share/boost_1_45_0
/intrepid-fs0/users/jnkim/persistent/qmcpack/extlibs
/soft/apps/hdf5-1.6.6
/soft/apps/fftw-3.1.2-double
)

# adjust the default behaviour of the FIND_XXX() commands:
# search headers and libraries in the target environment, search
# programs in the host environment
set(CMAKE_FIND_ROOT_PATH_MODE_PROGRAM NEVER)
set(CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY)
set(CMAKE_FIND_ROOT_PATH_MODE_INCLUDE ONLY)
set(CMAKE_SHARED_LINKER_FLAGS " ")

set(AIX_ARCH "450d")
SET(AIX_ARCH_FLAGS "-qarch=${AIX_ARCH}  -qsmp=omp -qthreaded")
SET(AIX_CXX_COMMON_FLAGS " -qkeyword=restrict -qstrict -qhot -qnoeh -qsuppress=1540-1090:1540-1088 ")
SET(AIX_OPT_FLAGS "-O3 -qmaxmem=-1  -qlargepage -qprefetch ")

#SET(AIX_OPT_FLAGS "-O3 -Q -qmaxmem=-1 -qipa=inline -qinline -qlargepage -qprefetch ")
#SET(AIX_CXX_OPT_FLAGS "-O3 -Q -qlargepage -qprefetch")
#SET(AIX_CXX_FLAGS "-O3 -Q -qlargepage -qprefetch")

SET(CMAKE_CXX_FLAGS "${AIX_ARCH_FLAGS} ${AIX_CXX_COMMON_FLAGS} ${AIX_OPT_FLAGS}")
SET(CMAKE_C_FLAGS "${AIX_ARCH_FLAGS} ${AIX_OPT_FLAGS}")

#link lapack, essl, mass
link_libraries(/soft/apps/LAPACK/lapack_3.3_BGP.a 
    /soft/apps/ESSL/lib/libesslbg.a
    -L/soft/apps/ibmcmp/xlf/bg/11.1/lib -lxlf90_r -lxlsmp  
    -L/soft/apps/ibmcmp/xlmass/bg/4.4/lib -lmass  
    )

