from distutils.core import setup, Extension
import numpy
import numpy.numarray as nn
numpyincludedirs = numpy.get_include()
numarrayincludedirs = nn.get_numarray_include_dirs()
module1 = Extension('IOSection',
                    sources = ['IONumPyWrapper.cc','IO.cc',\
                               'IOHDF5.cc', 'IOASCII.cc',\
                               'IOVarHDF5.cc' ],\
                    library_dirs=['/home/esler/lib',
                                  '/usr/include',\
                                  numpyincludedirs,\
                                  numarrayincludedirs[0]],
                    libraries =  ['blitz', 'xml2', 'hdf5'])

setup (name = 'IO',
       version = '1.0',
       description = 'This is a demo package',
       ext_modules = [module1],
       py_modules=['IO'])
