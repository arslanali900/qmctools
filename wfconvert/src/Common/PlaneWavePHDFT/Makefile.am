AM_LDFLAGS = $(ALL_STATIC)

check_PROGRAMS = TestPW2 BandStructure MPIBandStructure TestPWMPI TestLDA

TestPW2_SOURCES = TestPW2.cc
TestPWMPI_SOURCES = TestPWMPI.cc
BandStructure_SOURCES = BandStructure.cc BandStructure.h
MPIBandStructure_SOURCES = MPIBandStructure.cc MPIBandStructure.h
TestLDA_SOURCES = TestLDA.cc

LDADD = libpw.la                  ../IO/libIO.la                          \
        ../FFT/libfft.la          ../PH/libPH.la                          \
        ../MatrixOps/libmatops.la ../Splines/libspline.la                 \
        ../Fitting/libfitting.la  ../Integration/libquad.la               \
        ../MPI/libcomm.la         ../Ewald/libbreakup.la                  \
	../Atom/libatom.la        ../DFT/libdft.la                        \
        $(HDF5_LIBS) $(BLITZ_LIBS) $(G2C_LIBS) $(FFTW3_LIBS)              \
        $(LAPACK_LIBS) $(BLAS_LIBS) $(GSL_LIBS) $(FLIBS) $(MPI_LIBS)      \
        $(FLIBS)

TestLDA_LDADD = $(LDADD) ../Atom/libatom.la ../DFT/libdft.la $(FLIBS)

AM_CXXFLAGS = $(FFTW3_CFLAGS) $(BLITZ_CFLAGS)                             \
              $(HDF5_CFLAGS)  $(BLAS_CFLAGS) $(MPI_DEF) $(MPI_CFLAGS)     \
              $(FFT_DEF) $(CBLAS_DEF)

noinst_LTLIBRARIES = libpw.la
libpw_la_SOURCES = ChargeMixer.h      \
                   ChargeMixer.cc     \
	           ConjGrad2.h        \
                   ConjGrad2.cc       \
                   ConjGradMPI.h      \
                   ConjGradMPI.cc     \
                   Coulomb.h          \
                   Coulomb.cc         \
                   CoulombFFT.h       \
                   CoulombFFT.cc      \
                   FermiSmear.h       \
                   FermiSmear.cc      \
                   FFTBox.h           \
                   FFTBox.cc          \
                   GVecs.h            \
                   GVecs.cc           \
                   HamiltonianBase.h  \
                   HamiltonianBase.cc \
                   Hamiltonians.h     \
                   Hamiltonians.cc    \
                   LDA.cc             \
                   LocalPotFFT.h      \
                   LocalPotFFT.cc     \
                   NLPP_FFT.h         \
                   NLPP_FFT.cc        \
                   PHPot.h            \
                   PHPot.cc           \
                   PHPotFFT.h         \
                   PHPotFFT.cc        \
                   PlaneWaves.h       \
                   PlaneWaves.cc      \
                   PlaneWavesMPI.h    \
                   PlaneWavesMPI.cc   \
                   VectorOps.h

includedir=$(prefix)/include/Common/PlaneWavePHDFT


#if WITH_F77
#     libsprng_la_SOURCES += common/fwrap_mpi.c
#endif
libpw_la_CXXFLAGS = $(FFTW3_CFLAGS) $(BLITZ_CFLAGS) $(HDF5_CFLAGS)  \
                    $(CBLAS_CFLAGS) $(AM_CXXFLAGS) $(FFT_DEF)
libpw_la_LIBADD   =

#if HAVE_GMP
#     libsprng_la_SOURCES += pmlcg/pmlcg.h pmlcg/pmlcg.c pmlcg/info.h pmlcg/basic.h 
#endif
 
include_HEADERS = ConjGrad2.h FFTBox.h GVecs.h Hamiltonian2.h PlaneWaves.h \
                  VectorOps.h PlaneWavesMPI.h ConjGradMPI.h Hamiltonians.h \
                  Coulomb.h CoulombFFT.h HamiltonianBase.h LocalPotFFT.h   \
                  PHPot.h PHPotFFT.h ChargeMixer.h FermiSmear.h NLPP_FFT.h

