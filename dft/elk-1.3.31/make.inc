MAKE = make
F90 = ifort
F90_OPTS = -O3 -ip -unroll -no-prec-div -parallel -openmp -I/home/jnkim/share/intel11/hdf5-1.8.7/include
F77 = ifort
F77_OPTS = -O3 -ip -unroll -no-prec-div -parallel -openmp
AR = ar
LIB_SYS = 
LIB_LPK = -L/usr/local/intel/11.1.069/mkl/lib/em64t  -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lmkl_solver_lp64_sequential
LIB_FFT = fftlib.a -L/home/jnkim/share/intel11/hdf5-1.8.7/lib -lhdf5_fortran -lhdf5hl_fortran -lhdf5_hl -lhdf5