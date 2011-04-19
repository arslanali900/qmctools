#!/bin/env python
from pylab import *

c = load ('Chebyshev.dat')
N = c.shape[1]
for i in range(1,N):
    plot (c[:,0], c[:,i]);

figure(2)
c = load ('dChebyshev.dat')
N = c.shape[1]
for i in range(1,N):
    plot (c[:,0], c[:,i]);

figure(3)
c = load ('d2Chebyshev.dat')
N = c.shape[1]
for i in range(1,N):
    plot (c[:,0], c[:,i]);



show()
