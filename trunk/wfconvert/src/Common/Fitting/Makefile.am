AM_LDFLAGS = $(ALL_STATIC)

check_PROGRAMS = TestFitting TestBirch
bin_PROGRAMS   = LiFit VinetFit VinetFit2 BirchFit cBNModel MGFit DatchiRamanModel GoncharovRamanModel VinetFitP MyRamanModel

LDADD = libfitting.la ../MatrixOps/libmatops.la ../IO/libIO.la $(LAPACK_LIBS) $(BLITZ_LIBS) $(HDF5_LIBS) $(GSL_LIBS) $(SPRNG2_LIBS) $(FLIBS)

TestFitting_SOURCES = TestFitting.cc
LiFit_SOURCES       = TestEOS.cc LiEOS.h
VinetFit_SOURCES    = VinetFit.cc VinetEOS.h
VinetFitP_SOURCES   = VinetFitP.cc VinetEOS2.h
VinetFit2_SOURCES   = VinetFit2.cc PhononFreeEnergy.h PhononFreeEnergy.cc \
                      DebyeModel.h DebyeModel.cc VinetEOS2.h              \
                      RamanFreq.h RamanFreq.cc
BirchFit_SOURCES    = BirchFit.cc PhononFreeEnergy.h PhononFreeEnergy.cc  \
                      DebyeModel.h DebyeModel.cc VinetEOS2.h              \
                      RamanFreq.h RamanFreq.cc
TestBirch_SOURCES   = TestBirch.cc
cBNModel_SOURCES    = cBNModel.cc cBNModel.h Raman.h Raman.cc             \
	              PhononFreeEnergy.h PhononFreeEnergy.cc              \
                      DebyeModel.h DebyeModel.cc VinetEOS2.h
MGFit_SOURCES       = MGFit.cc MieGruineisen.h DebyeModel.h DebyeModel.cc
DatchiRamanModel_SOURCES = DatchiRamanModel.h DatchiRamanModel.cc
GoncharovRamanModel_SOURCES = GoncharovRamanModel.h GoncharovRamanModel.cc
MyRamanModel_SOURCES = MyRamanModel.h MyRamanModel.cc

cBNModel_LDADD = libfitting.la ../MatrixOps/libmatops.la ../Splines/libspline.la ../IO/libIO.la $(LAPACK_LIBS) $(BLITZ_LIBS) $(HDF5_LIBS) $(GSL_LIBS) $(SPRNG2_LIBS) $(FLIBS)

AM_CXXFLAGS = $(BLITZ_CFLAGS) $(HDF5_CFLAGS) $(XML2_CFLAGS) $(GSL_CFLAGS)

noinst_LTLIBRARIES = libfitting.la

libfitting_la_SOURCES = Fitting.cc        \
                        Fitting.h         \
                        NonlinearFit.h    \
                        NonlinearFitTemp.h

#libfitting_la_CXXFLAGS = $(BLITZ_CFLAGS) $(AM_CXXFLAGS)
include_HEADERS = Fitting.h NonlinearFit.h
includedir=$(prefix)/include/Common/Fitting

