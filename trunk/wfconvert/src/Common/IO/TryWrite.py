#!/bin/env python

from IO import *
from numpy import *

io = IOSectionClass()

io.NewFile ("test.h5")
io.WriteVar ("testfloat", 1.5)
io.WriteVar ("testbool", False)
io.WriteVar ("testint", 3)
io.WriteVar ("testcomplex", 1.0+2.0j);
io.WriteVar ("teststring", "Thomas");

I = array ([ 1, 2, 3, 4] )
M = array([1.1, 2.2, 3.3, 4.4])
B = array([True, True, False])
D = array([1.0 + 1.0j])
mat = array ([ [ 1, 2, 3, 4], [5, 6, 7, 8] ])
S = array(["Hello1", "world!!!!"])
S2 = array([["Hello1", "world!!!!"], ["abc", "123"]] )
S3 = array([[["Hello1", "world!!!!"], ["abc", "123"]], [["Hello1", "world!!!!"], ["abc", "123"]]] )
S4 = array([[[["Hello1", "world!!!!"], ["abc", "123"]], [["Hello1", "world!!!!"], ["abc", "123"]]]] )
io.WriteVar("ArrayDouble", M)
io.WriteVar("ArrayInt", I)
io.WriteVar("ArrayDouble2", mat)
io.WriteVar("ArrayString2", S2)
io.WriteVar("ArrayString3", S3)
io.WriteVar("ArrayString4", S4)

io.CloseFile()
