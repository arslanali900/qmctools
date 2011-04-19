from IOSection import *

handle = New()
OpenFile (handle, "junk.dat")

int_data = ReadVar (handle, "int_data")
print int_data

double_data = ReadVar (handle, "double_data")
print double_data

bool_data = ReadVar (handle, "bool_data")
print bool_data

string_data=ReadVar(handle,"string_data")
print string_data
