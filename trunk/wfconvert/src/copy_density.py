#!/bin/env python
from IO import *
import sys


infile = IOSectionClass()
infile.OpenFile (sys.argv[1])

outfile = IOSectionClass()
outfile.OpenFile(sys.argv[2])

infile.OpenSection('electrons')
outfile.OpenSection('electrons')
nspins   = infile.ReadVar('number_of_spins')
infile.OpenSection('density')
if (outfile.OpenSection('density') != 0):
    print 'There is already a density section in ' + sys.argv[1]
    print 'Exitting.'
    sys.exit(1)

outfile.NewSection('density')
gvectors = infile.ReadVar('gvectors')
mesh     = infile.ReadVar('mesh')
outfile.WriteVar('gvectors', gvectors)
outfile.WriteVar('mesh', mesh)
outfile.SetUnderscores(True)

for spin in range(0,nspins):
    infile.OpenSection2('spin', spin)
    success = outfile.OpenSection2('spin', spin)
    if (success == 0):
        outfile.NewSection('spin')
    density_r = infile.ReadVar('density_r')
    density_g = infile.ReadVar('density_g')
    outfile.WriteVar('density_r', density_r)
    outfile.WriteVar('density_g', density_g)
    
    infile.CloseSection()  # 'spin'
    outfile.CloseSection() # 'spin'

outfile.CloseFile()

    

