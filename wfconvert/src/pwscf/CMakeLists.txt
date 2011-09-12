set (QE_SRCS
  esh5_interfaces.c
  openfil_pp.f90
  stop_pp.f90
  )

set(QE_HOME /uf/ncsa/jnkim/svnwork/espresso-4.2)
INCLUDE_DIRECTORIES(${QE_HOME}/include)

link_libraries(
  ${QE_HOME}/PW/libpw.a
  ${QE_HOME}/EE/libee.a
  ${QE_HOME}/Modules/libqemod.a
  ${QE_HOME}/flib/ptools.a 
  ${QE_HOME}/flib/flib.a ${QE_HOME}/clib/clib.a
  ${QE_HOME}/iotk/src/libiotk.a
  )

ADD_EXECUTABLE(pw2qmcpack pw2qmcpack.f90 ${QE_SRCS})        
#ADD_EXECUTABLE(fit_sk  FitSk.cc)
#ADD_EXECUTABLE(test_tile  TestTile.cc)
#ADD_EXECUTABLE(test_spline  TestSpline.cc)
#
