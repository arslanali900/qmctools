A collection of utilities and tools to run QMC codes released under GPL v2.

## 2012-07-16 ##

### QE 4.2. with pw2qmcpack.x PP tool ###
pw2qmcpack.x : a utility to generate ESHDF file that can be used by QMCPACK
```
svn co http://qmctools.googlecode.com/svn/dft/espresso-4.2/
```

It is a fork of QE 4.2. See how to compile and build at [wiki](http://qmcpack.cmscc.org/how-to-guides/pwscf-converter).

### Tools to handle PPs and wavefunctions ###
wfconv/ppconvert are now available as a package and can be download by subversion as:
```
svn co http://qmctools.googlecode.com/svn/trunk/wfconvert
```

The distribution includes blitz header files and Common library. It uses the same build system based on cmake and needs these libraries to compile the tools
  * libxml2
  * hdf5
  * fftw3

Change: einspline is included with the distribution and there is no need to build the library separately.

With XYZ\_HOME, simply run
```
cd build
cmake ..
make
```

Consult QMCPACK wiki http://qmcpack.cmscc.org/  for the details.