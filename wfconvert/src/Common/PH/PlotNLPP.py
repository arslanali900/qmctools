#!/bin/env python
from pylab import *

zeta_r = load ('zeta_r.dat');
zeta_q = load ('zeta_q.dat');
chi_r  = load ('chi_r.dat');
chi_q  = load ('chi_q.dat');

plot (zeta_r[:,0], zeta_r[:,1], 'b',\
      chi_r[:,0], chi_r[:,1], 'b--',\
      zeta_r[:,0], zeta_r[:,2], 'r',\
      chi_r[:,0], chi_r[:,2], 'r--')
axis ((1.5, 1.58, -1.0e-3, 1.0e-3))

figure(2);
plot (zeta_q[:,0], zeta_q[:,1], 'b',\
      chi_q[:,0], chi_q[:,1], 'b--',\
      zeta_q[:,0], zeta_q[:,2], 'r',\
      chi_q[:,0], chi_q[:,2], 'r--')

show()
