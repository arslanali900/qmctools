/////////////////////////////////////////////////////////////
// Copyright (C) 2003-2006 Bryan Clark and Kenneth Esler   //
//                                                         //
// This program is free software; you can redistribute it  //
// and/or modify it under the terms of the GNU General     //
// Public License as published by the Free Software        //
// Foundation; either version 2 of the License, or         //
// (at your option) any later version.  This program is    //
// distributed in the hope that it will be useful, but     //
// WITHOUT ANY WARRANTY; without even the implied warranty //
// of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. //  
// See the GNU General Public License for more details.    //
// For more information, please see the PIMC++ Home Page:  //
//           http://pathintegrals.info                     //
/////////////////////////////////////////////////////////////

#include <Python.h>
#include <numpy/arrayobject.h>
#include "IO.h"

using namespace IO;

extern "C" PyObject* IOSection_SetUnderscores(PyObject *self, PyObject *args)
{
  void *IOSectionPtr;
  int use;
  char *str = (char*)( (sizeof(int) == sizeof(void*)) ? "ii" : "li");


  if (!PyArg_ParseTuple (args, str ,&IOSectionPtr,&use)) {
    cerr << "Error in IOSection_OpenFile.\n";
    return NULL;
  }
  else
    ((IOSectionClass*)IOSectionPtr)->SetUnderscores(use);
  return Py_None;
}


extern "C" PyObject* IOSection_New(PyObject *self, PyObject *args)
{
  IOSectionClass *newSect=new IOSectionClass();
  if (sizeof(void*) == sizeof(int))
    return Py_BuildValue("i", newSect);
  else
    return Py_BuildValue("l", newSect);
}


extern "C" PyObject*
IOSection_CountSectionsName (PyObject *self, PyObject *args)
{
  char *name;
  void *IOSectionPtr;
  char *str = (char*)( (sizeof(int) == sizeof(void*)) ? "is" : "ls");

  if (!PyArg_ParseTuple (args, str ,&IOSectionPtr,&name))
    return NULL;
  else {
    int num = ((IOSectionClass*)IOSectionPtr)->CountSections(name);
    return Py_BuildValue("i",(int)num);
  }
}




extern "C" PyObject*
IOSection_CountSections (PyObject *self, PyObject *args)
{
  void *IOSectionPtr;
  char *str = (char*)((sizeof(int) == sizeof(void*)) ? "i" : "l");

  if (!PyArg_ParseTuple (args, str ,&IOSectionPtr))
    return NULL;
  else {
    int num = ((IOSectionClass*)IOSectionPtr)->CountSections();
    return Py_BuildValue("i",num);
  }
}



extern "C" PyObject*
IOSection_CountVars (PyObject *self, PyObject *args)
{
  void *IOSectionPtr;
  char *str = (char*)((sizeof(int) == sizeof(void*)) ? "i" : "l");

  if (!PyArg_ParseTuple (args, str ,&IOSectionPtr))
    return NULL;
  else {
    int num = ((IOSectionClass*)IOSectionPtr)->CountVars();
    return Py_BuildValue("i",num);
  }
}



extern "C" PyObject*
IOSection_OpenFile (PyObject *self, PyObject *args)
{
  char *fileName;
  void *IOSectionPtr;
  char *str = (char*)((sizeof(int) == sizeof(void*)) ? "is" : "ls");

  if (!PyArg_ParseTuple (args, str ,&IOSectionPtr,&fileName)) {
    cerr << "Error in IOSection_OpenFile.\n";
    return NULL;
  }
  else {
    bool success = ((IOSectionClass*)IOSectionPtr)->OpenFile(fileName);
    return Py_BuildValue("i",(int)success);
  }
}


extern "C" PyObject*
IOSection_GetName(PyObject *self, PyObject *args)
{
  void *IOSectionPtr;
  char *str = (char*)((sizeof(int) == sizeof(void*)) ? "i" : "l");

  if (!PyArg_ParseTuple (args, str ,&IOSectionPtr))
    return NULL;
  else {
    string myName;
    myName = ((IOSectionClass*)IOSectionPtr)->GetName();
    return Py_BuildValue("s",myName.c_str());
  }
}



extern "C" PyObject*
IOSection_NewFile(PyObject *self, PyObject *args)
{
  char *fileName;
  void *IOSectionPtr;
  char *str = (char*)((sizeof(int) == sizeof(void*)) ? "is" : "ls");

  if (!PyArg_ParseTuple (args, str ,&IOSectionPtr,&fileName))
    return NULL;
  else {
    bool success;
    success = ((IOSectionClass*)IOSectionPtr)->NewFile(fileName);
    return Py_BuildValue("i",(int)success);
  }
}


extern "C" PyObject*
IOSection_CloseFile(PyObject *self, PyObject *args)
{
  void *IOSectionPtr;
  char *str = (char*)((sizeof(int) == sizeof(void*)) ? "i" : "l");

  if (!PyArg_ParseTuple (args, str ,&IOSectionPtr))
    return NULL;
  else {
    ((IOSectionClass*)IOSectionPtr)->CloseFile();
    return Py_True;
  }
}


extern "C" PyObject*
IOSection_FlushFile(PyObject *self, PyObject *args)
{
  void *IOSectionPtr;
  char *str = (char*)((sizeof(int) == sizeof(void*)) ? "i" : "l");

  if (!PyArg_ParseTuple (args, str ,&IOSectionPtr))
    return NULL;
  else {
    ((IOSectionClass*)IOSectionPtr)->FlushFile();
    return Py_True;
    //    return NULL;
  }
}

extern "C" PyObject*
IOSection_OpenSectionName (PyObject *self, PyObject *args)
{
  char *sectionName;
  void *IOSectionPtr;
  char *str = (char*)((sizeof(int) == sizeof(void*)) ? "is" : "ls");

  if (!PyArg_ParseTuple (args, str ,&IOSectionPtr,&sectionName))
    return NULL;
  else {
    bool success = ((IOSectionClass*)IOSectionPtr)->OpenSection(sectionName);
    return Py_BuildValue("i",(int)success);
  }
}


extern "C" PyObject*
IOSection_OpenSectionNameNum (PyObject *self, PyObject *args)
{
  char *sectionName;
  void *IOSectionPtr;
  int num;
  char *str = (char*)((sizeof(int) == sizeof(void*)) ? "isi" : "lsi");

  if (!PyArg_ParseTuple (args, str, &IOSectionPtr,&sectionName,&num))
    return NULL;
  else {
    bool success = ((IOSectionClass*)IOSectionPtr)->OpenSection(sectionName,num);
    return Py_BuildValue("i",(int)success);
  }
}


extern "C" PyObject*
IOSection_OpenSectionNum (PyObject *self, PyObject *args)
{
  void *IOSectionPtr;
  int num;
  char *str = (char*)((sizeof(int) == sizeof(void*)) ? "ii" : "li");

  if (!PyArg_ParseTuple (args, str, &IOSectionPtr,&num))
    return NULL;
  else {
    bool success = ((IOSectionClass*)IOSectionPtr)->OpenSection(num);
    return Py_BuildValue("i",(int)success);
  }
}

extern "C" PyObject*
IOSection_IncludeSection (PyObject *self, PyObject *args)
{
  void *IOSectionPtr;
  char *name;
  char *fileName;
  char *str = (char*)((sizeof(int) == sizeof(void*)) ? "iss" : "lss");

  if (!PyArg_ParseTuple (args, str, &IOSectionPtr,&name,&fileName))
    return NULL;
  else {
   bool success = ((IOSectionClass*)IOSectionPtr)->IncludeSection(name,fileName);
   return Py_BuildValue("i",(int)success);
  }
}


extern "C" PyObject*
IOSection_NewSectionName(PyObject *self, PyObject *args)
{
  char *sectionName;
  void *IOSectionPtr;
  char *str = (char*)((sizeof(int) == sizeof(void*)) ? "is" : "ls");

  if (PyArg_ParseTuple (args, str, &IOSectionPtr,&sectionName))
    ((IOSectionClass*)IOSectionPtr)->NewSection(sectionName);
  bool success = true;
   return Py_BuildValue("i",(int)success);
   //  return NULL;
}

extern "C" PyObject*
IOSection_NewSectionFile (PyObject *self, PyObject *args)
{
  void *IOSectionPtr;
  char *name;
  char *fileName;
  char *str = (char*)((sizeof(int) == sizeof(void*)) ? "iss" : "lss");

  if (!PyArg_ParseTuple (args, str, &IOSectionPtr,&name,&fileName))
    return NULL;
  else {
   bool success = ((IOSectionClass*)IOSectionPtr)->NewSection(name,fileName);
   return Py_BuildValue("i",(int)success);
  }
}

extern "C" PyObject*
IOSection_CloseSection(PyObject *self, PyObject *args)
{
  void *IOSectionPtr;
  char *str = (char*)((sizeof(int) == sizeof(void*)) ? "i" : "l");

  if (!PyArg_ParseTuple (args, str, &IOSectionPtr))
    return NULL;
  else {
    ((IOSectionClass*)IOSectionPtr)->CloseSection();
    return Py_None;
  }
}

extern "C" PyObject*
IOSection_GetVarName(PyObject *self, PyObject *args)
{
  int num;
  void *IOSectionPtr;
  char *str = (char*)((sizeof(int) == sizeof(void*)) ? "ii" : "li");

  if (!PyArg_ParseTuple (args, str, &IOSectionPtr,&num))
    return NULL;
  else {
    string name=((IOSectionClass*)IOSectionPtr)->GetVarName(num);
    return Py_BuildValue("s",(name.c_str()));
  }
}


extern "C" PyObject*
IOSection_WriteVar (PyObject *self, PyObject *args)
{
  // Get the name from the python arguments
  char *str = (char*)((sizeof(int) == sizeof(void*)) ? "isO" : "lsO");
  const char *name;
  void *IOSectionPtr;
  PyObject *dataObject;
  if (!PyArg_ParseTuple (args, str, &IOSectionPtr, &name, &dataObject))
    return Py_None;

  IOSectionClass &io = *((IOSectionClass*)IOSectionPtr);

  /////////////////////////////////////////////////
  //              Atomic Writes                  //
  /////////////////////////////////////////////////
  if (PyBool_Check (dataObject)) {
    bool val = (dataObject == Py_True);
    io.WriteVar (name, val);
  }
  else if (PyInt_Check(dataObject)) {
    int val = PyInt_AsLong (dataObject);
    io.WriteVar (name, val);
  }
  else if (PyFloat_Check(dataObject)) {
    double val = PyFloat_AsDouble (dataObject);
    io.WriteVar (name, val);
  }
  else if (PyComplex_Check (dataObject)) {
    complex<double> val (PyComplex_RealAsDouble(dataObject), 
			 PyComplex_ImagAsDouble(dataObject));
    io.WriteVar (name, val);
  }
  else if (PyString_Check (dataObject)) {
    string str = PyString_AS_STRING (dataObject);
    io.WriteVar (name, str);
  }
  else if (PyList_Check (dataObject)) {
  }
  ////////////////////////////////////////////////////
  //                  Array Writes                  //
  ////////////////////////////////////////////////////
  else if (PyArray_Check (dataObject)) {
    PyArrayObject *array = (PyArrayObject*) dataObject;
    if (array->descr->type == 'l') {
      int safe = PyArray_CanCastSafely(PyArray_LONG, PyArray_INT);
      array = (PyArrayObject*) PyArray_Cast(array, PyArray_INT);
    }
    PyArray_Descr *descr = array->descr;
    char type = descr->type;
    int ndim = array->nd;
    int typesize;
    switch (type) {
      case NPY_INTLTR:
	typesize = sizeof(int);              break;
      case NPY_DOUBLELTR:
	typesize = sizeof(double);           break;
      case NPY_BOOLLTR:
	typesize = sizeof(bool);             break;
      case NPY_CDOUBLE:
	typesize = sizeof(complex<double>);  break;
      case NPY_STRING:
        typesize = sizeof(char);             break;
      default:
        typesize = 1;                        break;
    }
    //    cerr << "typecode = " << type << endl;
    ////////
    // 1D //
    ////////
    if (ndim == 1) {
      TinyVector<int,1> tvdims, tvstrides; 
      tvdims[0]     = array->dimensions[0];
      tvstrides[0]  = array->strides[0] / typesize;

      if (type == NPY_INTLTR) {
	int *data = (int*) array->data;
	Array<int,1> blitzArray(data, tvdims, tvstrides, blitz::neverDeleteData);
	io.WriteVar (name, blitzArray);
      }
      else if (type == NPY_DOUBLELTR) {
	double *data = (double*) array->data;
	Array<double,1> blitzArray(data, tvdims, tvstrides, blitz::neverDeleteData);
	io.WriteVar (name, blitzArray);
      }
      else if (type == NPY_BOOLLTR){
	double *data = (double*) array->data;
	Array<double,1> blitzArray(data, tvdims, tvstrides, blitz::neverDeleteData);
	io.WriteVar (name, blitzArray);
      }
      else if (type == NPY_CDOUBLELTR) {
	complex<double> *data = (complex<double>*) array->data;
	Array<complex<double>,1> blitzArray(data, tvdims, tvstrides, blitz::neverDeleteData);
	io.WriteVar (name, blitzArray);
      }
      else if (type == NPY_STRINGLTR) {
	char tmp[tvstrides[0]+1];
	tmp[tvstrides[0]] = '\0';
	Array<string,1> strArray(tvdims[0]);
	for (int i=0; i<tvdims[0]; i++) {
	  for (int j=0; j<tvstrides[0]; j++)
	    tmp[j] = ((char*)array->data)[i*tvstrides[0]+j];
	  strArray(i) = tmp;
	}
	io.WriteVar(name, strArray);
      }
    }
    ////////
    // 2D //
    ////////
    else if (ndim == 2) {
      TinyVector<int,2> tvdims, tvstrides; 
      for (int i=0; i<2; i++) {
	tvdims[i]     = array->dimensions[i];
	tvstrides[i]  = array->strides[i] / typesize;
      }
      
      if (type == NPY_INTLTR) {
	int *data = (int*) array->data;
	Array<int,2> blitzArray(data, tvdims, tvstrides, blitz::neverDeleteData);
	io.WriteVar (name, blitzArray);
      }
      else if (type == NPY_DOUBLELTR) {
	double *data = (double*) array->data;
	Array<double,2> blitzArray(data, tvdims, tvstrides, blitz::neverDeleteData);
	io.WriteVar (name, blitzArray);
      }
      else if (type == NPY_BOOLLTR){
	double *data = (double*) array->data;
	Array<double,2> blitzArray(data, tvdims, tvstrides, blitz::neverDeleteData);
	io.WriteVar (name, blitzArray);
      }
      else if (type == NPY_CDOUBLELTR) {
	complex<double> *data = (complex<double>*) array->data;
	Array<complex<double>,2> blitzArray(data, tvdims, tvstrides, blitz::neverDeleteData);
	io.WriteVar (name, blitzArray);
      }
      else if (type == NPY_STRINGLTR) {
	char tmp[tvstrides[1]+1];
	tmp[tvstrides[1]] = '\0';
	Array<string,2> strArray(tvdims[0], tvdims[1]);
	for (int i=0; i<tvdims[0]; i++) 
	  for (int j=0; j<tvdims[1]; j++) {
	    for (int k=0; k<tvstrides[1]; k++)
	      tmp[k] = ((char*)array->data)[i*tvstrides[0]+j*tvstrides[1]+k];
	    strArray(i,j) = tmp;
	  }
	io.WriteVar(name, strArray);
      }

    }
    ////////
    // 3D //
    ////////
    else if (ndim == 3) {
      TinyVector<int,3> tvdims, tvstrides; 
      for (int i=0; i<3; i++) {
	tvdims[i]     = array->dimensions[i];
	tvstrides[i]  = array->strides[i] / typesize;
      }
      
      if (type == NPY_INTLTR) {
	int *data = (int*) array->data;
	Array<int,3> blitzArray(data, tvdims, tvstrides, blitz::neverDeleteData);
	io.WriteVar (name, blitzArray);
      }
      else if (type == NPY_DOUBLELTR) {
	double *data = (double*) array->data;
	Array<double,3> blitzArray(data, tvdims, tvstrides, blitz::neverDeleteData);
	io.WriteVar (name, blitzArray);
      }
      else if (type == NPY_BOOLLTR){
	double *data = (double*) array->data;
	Array<double,3> blitzArray(data, tvdims, tvstrides, blitz::neverDeleteData);
	io.WriteVar (name, blitzArray);
      }
      else if (type == NPY_CDOUBLELTR) {
	complex<double> *data = (complex<double>*) array->data;
	Array<complex<double>,3> blitzArray(data, tvdims, tvstrides, blitz::neverDeleteData);
	io.WriteVar (name, blitzArray);
      }
      else if (type == NPY_STRINGLTR) {
	char tmp[tvstrides[2]+1];
	tmp[tvstrides[2]] = '\0';
	Array<string,3> strArray(tvdims[0], tvdims[1], tvdims[2]);
	for (int i=0; i<tvdims[0]; i++) 
	  for (int j=0; j<tvdims[1]; j++) 
	    for (int k=0; k<tvdims[2]; k++) {
	      for (int l=0; l<tvstrides[2]; l++) 
		tmp[l] = 
		  ((char*)array->data) [i*tvstrides[0]+j*tvstrides[1]+k*tvstrides[2]+l];
	      strArray(i,j,k) = tmp;
	    }
	io.WriteVar(name, strArray);
      }

    }
    ////////
    // 4D //
    ////////
    else if (ndim == 4) {
      TinyVector<int,4> tvdims, tvstrides; 
      for (int i=0; i<4; i++) {
	tvdims[i]     = array->dimensions[i];
	tvstrides[i]  = array->strides[i] / typesize;
      }
      
      if (type == NPY_INTLTR) {
	int *data = (int*) array->data;
	Array<int,4> blitzArray(data, tvdims, tvstrides, blitz::neverDeleteData);
	io.WriteVar (name, blitzArray);
      }
      else if (type == NPY_DOUBLELTR) {
	double *data = (double*) array->data;
	Array<double,4> blitzArray(data, tvdims, tvstrides, blitz::neverDeleteData);
	io.WriteVar (name, blitzArray);
      }
      else if (type == NPY_BOOLLTR){
	double *data = (double*) array->data;
	Array<double,4> blitzArray(data, tvdims, tvstrides, blitz::neverDeleteData);
	io.WriteVar (name, blitzArray);
      }
      else if (type == NPY_CDOUBLELTR) {
	complex<double> *data = (complex<double>*) array->data;
	Array<complex<double>,4> blitzArray(data, tvdims, tvstrides, blitz::neverDeleteData);
	io.WriteVar (name, blitzArray);
      }
      else if (type == NPY_STRINGLTR) {
	char tmp[tvstrides[3]+1];
	tmp[tvstrides[3]] = '\0';
	Array<string,4> strArray(tvdims[0], tvdims[1], tvdims[2], tvdims[3]);
	for (int i=0; i<tvdims[0]; i++) 
	  for (int j=0; j<tvdims[1]; j++) 
	    for (int k=0; k<tvdims[2]; k++) 
	      for (int l=0; l<tvdims[3]; l++) {
		for (int m=0; m<tvstrides[3]; m++) 
		  tmp[m] = ((char*)array->data) 
		    [i*tvstrides[0]+j*tvstrides[1]+k*tvstrides[2]+l*tvstrides[3]+m];
		strArray(i,j,k,l) = tmp;
	      }
	io.WriteVar(name, strArray);
      }
    }
  }
  else {
    cerr << "Error:  unrecognized object type passed to WriteVar.\n";
    return Py_False;
  }
  return Py_True;
}


extern "C" PyObject*
IOSection_ReadVar(PyObject *self, PyObject *args)
{
  const char *name;    
  void *IOSectionPtr;
  IOVarBase *varPtr;
  bool success;
  //Checking to see if we are passed a string or an int
  PyObject* handle;
  PyObject* varToRead;
  handle=PyTuple_GetItem(args,0);
  varToRead=PyTuple_GetItem(args,1);
  if (PyInt_Check(varToRead)){
    char *str = (char*)((sizeof(int) == sizeof(void*)) ? "ii" : "li");
    int num;
    if (!PyArg_ParseTuple (args, str, &IOSectionPtr,&num)){
      return Py_None;
    }
    varPtr=((IOSectionClass*)IOSectionPtr)->GetVarPtr(num);
    name=varPtr->GetName().c_str();
  }
  else {
    char *str = (char*)((sizeof(int) == sizeof(void*)) ? "is" : "ls");
    if (!PyArg_ParseTuple (args, str, &IOSectionPtr,&name)){
      return Py_None;
    }
    varPtr=((IOSectionClass*)IOSectionPtr)->GetVarPtr(name);
  }
  varPtr=((IOSectionClass*)IOSectionPtr)->GetVarPtr(name);
  if (varPtr==NULL)
    return Py_None;
  IODataType type=varPtr->GetType();
  int dim=varPtr->GetRank();
  if (dim == 0) {
    if (type==INT_TYPE){
      int val;    
      success = ((IOSectionClass*)IOSectionPtr)->ReadVar(name,val);
      if (success)
	return Py_BuildValue("i",val);
      else
	return Py_None;
    }
    else if (type==DOUBLE_TYPE){
      double val;    
      success = ((IOSectionClass*)IOSectionPtr)->ReadVar(name,val);
      if (success)
	return Py_BuildValue("d",val);
      else
	return Py_None;
    }
    else if (type==COMPLEX_TYPE){
      complex<double> val;    
      success = ((IOSectionClass*)IOSectionPtr)->ReadVar(name,val);
      if (success) {
	Py_complex z;
	z.real = val.real();
	z.imag = val.imag();
	return Py_BuildValue("D",&z);
      }
      else
	return Py_None;
    }
    else if (type==STRING_TYPE){
      string val;    
      success = ((IOSectionClass*)IOSectionPtr)->ReadVar(name,val);
      if (success)
	return Py_BuildValue("s",val.c_str());
      else
	return Py_None;
    }
    else if (type==BOOL_TYPE){
      bool val;    
      success = ((IOSectionClass*)IOSectionPtr)->ReadVar(name,val);
      if (success)
	return Py_BuildValue("i",(int)val);
      else
	return Py_None;
    }
  }
  // 1D arrays
  else if (dim == 1) {
    if (type==INT_TYPE){
      blitz::Array<int,1> val;
      success = ((IOSectionClass*)IOSectionPtr)->ReadVar(name,val);
      npy_intp dims = val.size();
      if (success) {
	PyArrayObject *array;
	//array = NA_NewArray (NULL, NPY_INT, 1, len);
	array = (PyArrayObject*)PyArray_SimpleNew (1, &dims, NPY_INT);
	// Now copy data into new array
	for (int i=0; i<dims; i++)
	  *(((int*)array->data)+i) = val(i);
	return (PyObject *) array;
      }
      else
	return Py_None;
    }
    if (type==DOUBLE_TYPE){
      blitz::Array<double,1> val;
      success = ((IOSectionClass*)IOSectionPtr)->ReadVar(name,val);
      npy_intp dims = val.size();
      if (success) {
	PyArrayObject *array;
	//	array = NA_NewArray (NULL, NPY_DOUBLE, 1, len);
	array = (PyArrayObject*)PyArray_SimpleNew (1, &dims, NPY_DOUBLE);
	// Now copy data into new array
	for (int i=0; i<dims; i++)
	  *(((double*)array->data)+i) = val(i);
	return (PyObject *) array;
      }
      else
	return Py_None;
    }
    if (type==COMPLEX_TYPE){
      blitz::Array<complex<double>,1> val;
      success = ((IOSectionClass*)IOSectionPtr)->ReadVar(name,val);
      int len = val.size();
      npy_intp dims = val.size();
      if (success) {
	PyArrayObject *array;
	//	array = NA_NewArray (NULL, NPY_CDOUBLE, 1, len);
	array = (PyArrayObject*)PyArray_SimpleNew (1, &dims, NPY_CDOUBLE);
	// Now copy data into new array
	for (int i=0; i<dims; i++)
	  *(((complex<double>*)array->data)+i) = val(i);
	return (PyObject *) array;
      }
      else
	return Py_None;
    }
    if (type==BOOL_TYPE){
      blitz::Array<bool,1> val;
      success = ((IOSectionClass*)IOSectionPtr)->ReadVar(name,val);
      int len = val.size();
      npy_intp dims = val.size();
      if (success) {
	PyArrayObject *array;
	// array = NA_NewArray (NULL, NPY_BOOL, 1, len);
	array = (PyArrayObject*)PyArray_SimpleNew (1, &dims, NPY_BOOL);
	// Now copy data into new array
	for (int i=0; i<dims; i++)
	  *(((bool*)array->data)+i) = val(i);
	return (PyObject *) array;
      }
      else
	return Py_None;
    }
    else if (type==STRING_TYPE){
      blitz::Array<string,1> val;    
      success = ((IOSectionClass*)IOSectionPtr)->ReadVar(name,val);
      if (success){
	int len=val.size();
	PyObject* array=PyList_New(0);
	for (int i=0;i<len;i++){
	  PyList_Append(array,Py_BuildValue("s",val(i).c_str()));
	}
	return (PyObject*)array;
      }
      else
	return Py_None;
    }
  }
  // 2D arrays
  else if (dim == 2) {
    if (type==INT_TYPE){
      blitz::Array<int,2> val;
      success = ((IOSectionClass*)IOSectionPtr)->ReadVar(name,val);
      npy_intp dims[2] = { val.extent(0), val.extent(1) };
      if (success) {
	PyArrayObject *array;
	// array = NA_NewArray (NULL, NPY_INT, 2, d0, d1);
	array = (PyArrayObject*)PyArray_SimpleNew (2, dims, NPY_INT);
	// Now copy data into new array
	for (int i=0; i<dims[0]; i++)
	  for (int j=0; j<dims[1]; j++)
	    *(((int*)array->data)+(i*dims[1])+j) = val(i,j);
	return (PyObject *) array;
      }
      else
	return Py_None;
    }
    if (type==DOUBLE_TYPE) {
      blitz::Array<double,2> val;
      success = ((IOSectionClass*)IOSectionPtr)->ReadVar(name,val);
      npy_intp dims[2] = { val.extent(0), val.extent(1) };
      if (success) {
	PyArrayObject *array;
	// array = NA_NewArray (NULL, NPY_DOUBLE, 2, d0, d1);
	array = (PyArrayObject*)PyArray_SimpleNew (2, dims, NPY_DOUBLE);
	// Now copy data into new array
	for (int i=0; i<dims[0]; i++)
	  for (int j=0; j<dims[1]; j++)
	    *(((double*)array->data)+(i*dims[1])+j) = val(i,j);
	return (PyObject *) array;
      }
      else
	return Py_None;
    }
    if (type==COMPLEX_TYPE) {
      blitz::Array<complex<double>,2> val;
      npy_intp dims[2] = { val.extent(0), val.extent(1) };
      success = ((IOSectionClass*)IOSectionPtr)->ReadVar(name,val);
      if (success) {
	PyArrayObject *array;
	// array = NA_NewArray (NULL, NPY_CDOUBLE, 2, d0, d1);
	array = (PyArrayObject*)PyArray_SimpleNew (2, dims, NPY_CDOUBLE);
	// Now copy data into new array
	for (int i=0; i<dims[0]; i++)
	  for (int j=0; j<dims[1]; j++)
	    *(((complex<double>*)array->data)+(i*dims[1])+j) = val(i,j);
	return (PyObject *) array;
      }
      else
	return Py_None;
    }
    if (type==BOOL_TYPE){
      blitz::Array<bool,2> val;
      success = ((IOSectionClass*)IOSectionPtr)->ReadVar(name,val);
      npy_intp dims[2] = { val.extent(0), val.extent(1) };
      if (success) {
	PyArrayObject *array;
	//array = NA_NewArray (NULL, NPY_BOOL, 2, d0, d1);
	array = (PyArrayObject*)PyArray_SimpleNew (2, dims, NPY_BOOL);
	// Now copy data into new array
	for (int i=0; i<dims[0]; i++)
	  for (int j=0; j<dims[1]; j++)
	    *(((bool*)array->data)+(i*dims[1])+j) = val(i,j);
	return (PyObject *) array;
      }
      else
	return Py_None;
    }
    else if (type==STRING_TYPE){
      blitz::Array<string,2> val;    
      success = ((IOSectionClass*)IOSectionPtr)->ReadVar(name,val);
      if (success){
	int d0=val.extent(0);
	int d1=val.extent(1);
	PyObject* totalArray=PyList_New(0);
	for (int i=0;i<d0;i++){
	  PyObject* array=PyList_New(0);
	  for (int j=0;j<d1;j++){
	    PyList_Append(array,Py_BuildValue("s",val(i,j).c_str()));
	  }
	  PyList_Append(totalArray,array);
	}
	return (PyObject*)totalArray;
      }
      else
	return Py_None;
    }
  }
  // 3D arrays
  else if (dim == 3) {
    if (type==INT_TYPE){
      blitz::Array<int,3> val;
      success = ((IOSectionClass*)IOSectionPtr)->ReadVar(name,val);
      npy_intp dims[3] = { val.extent(0), val.extent(1), val.extent(2) };
      if (success) {
	PyArrayObject *array;
	//array = NA_NewArray (NULL, NPY_INT, 3, d0, d1, d2);
	array = (PyArrayObject*)PyArray_SimpleNew (3, dims, NPY_INT);
	// Now copy data into new array
	for (int i=0; i<dims[0]; i++)
	  for (int j=0; j<dims[1]; j++)
	    for (int k=0; k<dims[2]; k++)
	      *(((int*)array->data)+(i*dims[1]*dims[2])+(j*dims[2])+k) = val(i,j,k);
	return (PyObject *) array;
      }
      else
	return Py_None;
    }
    if (type==DOUBLE_TYPE) {
      blitz::Array<double,3> val;
      success = ((IOSectionClass*)IOSectionPtr)->ReadVar(name,val);
      npy_intp dims[3] = { val.extent(0), val.extent(1), val.extent(2) };
      if (success) {
	PyArrayObject *array;
	//array = NA_NewArray (NULL, NPY_DOUBLE, 3, d0, d1, d2);
	array = (PyArrayObject*)PyArray_SimpleNew (3, dims, NPY_DOUBLE);
	// Now copy data into new array
	for (int i=0; i<dims[0]; i++)
	  for (int j=0; j<dims[1]; j++)
	    for (int k=0; k<dims[2]; k++)
	      *(((double*)array->data)+(i*dims[1]*dims[2])+j*dims[2]+k) = val(i,j,k);
	return (PyObject *) array;
      }
      else
	return Py_None;
    }
    if (type==COMPLEX_TYPE) {
      blitz::Array<complex<double>,3> val;
      success = ((IOSectionClass*)IOSectionPtr)->ReadVar(name,val);
      npy_intp dims[3] = { val.extent(0), val.extent(1), val.extent(2) };
      if (success) {
	PyArrayObject *array;
	//array = NA_NewArray (NULL, NPY_CDOUBLE, 3, d0, d1, d2);
	array = (PyArrayObject*)PyArray_SimpleNew (3, dims, NPY_CDOUBLE);
	// Now copy data into new array
	for (int i=0; i<dims[0]; i++)
	  for (int j=0; j<dims[1]; j++)
	    for (int k=0; k<dims[2]; k++)
	      *(((complex<double>*)array->data)+(i*dims[1]*dims[2])+j*dims[2]+k) = val(i,j,k);
	return (PyObject *) array;
      }
      else
	return Py_None;
    }
    if (type==BOOL_TYPE){
      blitz::Array<bool,3> val;
      success = ((IOSectionClass*)IOSectionPtr)->ReadVar(name,val);
      npy_intp dims[3] = { val.extent(0), val.extent(1), val.extent(2) };
      if (success) {
	PyArrayObject *array;
	//array = NA_NewArray (NULL, NPY_BOOL, 3, d0, d1, d2);
	array = (PyArrayObject*)PyArray_SimpleNew (3, dims, NPY_BOOL);
	// Now copy data into new array
	for (int i=0; i<dims[0]; i++)
	  for (int j=0; j<dims[1]; j++)
	    for (int k=0; k<dims[2]; k++)
	      *(((bool*)array->data)+(i*dims[1]*dims[2])+j*dims[2]+k) = val(i,j,k);
	return (PyObject *) array;
      }
      else
	return Py_None;
    }
    else if (type==STRING_TYPE){
      blitz::Array<string,3> val;    
      success = ((IOSectionClass*)IOSectionPtr)->ReadVar(name,val);
      int d0 = val.extent(0); int d1 = val.extent(1); int d2 = val.extent(2);
      if (success){
	PyObject* l0=PyList_New(0);
	for (int i=0; i<d0; i++) {
	  PyObject* l1=PyList_New(0);
	  for (int j=0;j<d1;j++) {
	    PyObject* l2 = PyList_New(0);
	    for (int k=0; k<d2; k++)
	      PyList_Append(l2,Py_BuildValue("s",val(i,j,k).c_str()));
	    PyList_Append(l1,l2);
	  }
	  PyList_Append (l0, l1);
	}
	return (PyObject*)l0;
      }
      else
	return Py_None;
    }
  }
  // 4D arrays
  else if (dim == 4) {
    if (type==INT_TYPE){
      blitz::Array<int,4> val;
      success = ((IOSectionClass*)IOSectionPtr)->ReadVar(name,val);
      npy_intp dims[4] = { val.extent(0), val.extent(1), 
			   val.extent(2), val.extent(3) };
      if (success) {
	PyArrayObject *array;
	// array = NA_NewArray (NULL, NPY_INT, 4, d0, d1, d2, d3);
	array = (PyArrayObject*)PyArray_SimpleNew (4, dims, NPY_INT);
	// Now copy data into new array
	for (int i=0; i<dims[0]; i++)
	  for (int j=0; j<dims[1]; j++)
	    for (int k=0; k<dims[2]; k++)
	      for (int m=0; m<dims[3]; m++)
		*(((int*)array->data)+(i*dims[1]*dims[2]*dims[3])+(j*dims[2]*dims[3])+k*dims[3]+m) = 
		  val(i,j,k,m);
	return (PyObject *) array;
      }
      else
	return Py_None;
    }
    if (type==DOUBLE_TYPE) {
      blitz::Array<double,4> val;
      success = ((IOSectionClass*)IOSectionPtr)->ReadVar(name,val);
      npy_intp dims[4] = { val.extent(0), val.extent(1), 
			   val.extent(2), val.extent(3) };
      if (success) {
	PyArrayObject *array;
	// array = NA_NewArray (NULL, NPY_DOUBLE, 4, d0, d1, d2, d3);
	array = (PyArrayObject*)PyArray_SimpleNew (4, dims, NPY_DOUBLE);
	// Now copy data into new array
	for (int i=0; i<dims[0]; i++)
	  for (int j=0; j<dims[1]; j++)
	    for (int k=0; k<dims[2]; k++)
	      for (int m=0; m<dims[3]; m++)
		*(((double*)array->data)+(i*dims[1]*dims[2]*dims[3])+j*dims[2]*dims[3]+k*dims[3]+m) = 
		  val(i,j,k,m);
	return (PyObject *) array;
      }
      else
	return Py_None;
    }
    if (type==COMPLEX_TYPE) {
      blitz::Array<complex<double>,4> val;
      success = ((IOSectionClass*)IOSectionPtr)->ReadVar(name,val);
      npy_intp dims[4] = { val.extent(0), val.extent(1), 
			   val.extent(2), val.extent(3) };
      if (success) {
	PyArrayObject *array;
	// array = NA_NewArray (NULL, NPY_CDOUBLE, 4, d0, d1, d2, d3);
	array = (PyArrayObject*)PyArray_SimpleNew (4, dims, NPY_CDOUBLE);
	// Now copy data into new array
	for (int i=0; i<dims[0]; i++)
	  for (int j=0; j<dims[1]; j++)
	    for (int k=0; k<dims[2]; k++)
	      for (int m=0; m<dims[3]; m++)
		*(((complex<double>*)array->data)+(i*dims[1]*dims[2]*dims[3])+j*dims[2]*dims[3]+k*dims[3]+m) 
		  = val(i,j,k,m);
	return (PyObject *) array;
      }
      else
	return Py_None;
    }
    if (type==BOOL_TYPE){
      blitz::Array<bool,4> val;
      success = ((IOSectionClass*)IOSectionPtr)->ReadVar(name,val);
      npy_intp dims[4] = { val.extent(0), val.extent(1), 
			   val.extent(2), val.extent(3) };
      if (success) {
	PyArrayObject *array;
	// array = NA_NewArray (NULL, NPY_BOOL, 4, d0, d1, d2, d3);
	array = (PyArrayObject*)PyArray_SimpleNew (4, dims, NPY_BOOL);
	// Now copy data into new array
	for (int i=0; i<dims[0]; i++)
	  for (int j=0; j<dims[1]; j++)
	    for (int k=0; k<dims[2]; k++)
	      for (int m=0; m<dims[3]; m++)
		*(((bool*)array->data)+(i*dims[1]*dims[2]*dims[3])+j*dims[2]*dims[3]+k*dims[3]+m) = 
		  val(i,j,k,m);
	return (PyObject *) array;
      }
      else
	return Py_None;
    }
    else if (type==STRING_TYPE){
      blitz::Array<string,4> val;    
      success = ((IOSectionClass*)IOSectionPtr)->ReadVar(name,val);
      int d0 = val.extent(0); int d1 = val.extent(1); 
      int d2 = val.extent(2); int d3 = val.extent(3);
      if (success){
	PyObject* l0=PyList_New(0);
	for (int i=0; i<d0; i++) {
	  PyObject* l1=PyList_New(0);
	  for (int j=0;j<d1;j++) {
	    PyObject* l2 = PyList_New(0);
	    for (int k=0; k<d2; k++) {
	      PyObject* l3 = PyList_New(0);
	      for (int m=0; m<d3; m++)
		PyList_Append(l3,Py_BuildValue("s",val(i,j,k,m).c_str()));
	      PyList_Append(l2,l3);
	    }
	    PyList_Append (l1, l2);
	  }
	  PyList_Append(l0, l1);
	}
	return (PyObject*)l0;
      }
      else
	return Py_None;
    }
  }
  return (PyObject*)NULL;
}





static PyMethodDef IOSectionMethods[] = {
    {"OpenFile",  IOSection_OpenFile, METH_VARARGS,
     "Open a file for reading/writing."},
    {"New", IOSection_New, METH_VARARGS,
     "Create a new IOSection object."},
    {"ReadVar", IOSection_ReadVar, METH_VARARGS,
     "Reads the named variable, returning appropriate object with data."},
    {"WriteVar", IOSection_WriteVar, METH_VARARGS,
     "Writes the named variable with the data pased."},
    {"GetName", IOSection_GetName, METH_VARARGS,
     "Gets the name of the current section"},
    {"NewFile", IOSection_NewFile, METH_VARARGS,
     "Creates a new file"},
    {"CloseFile", IOSection_CloseFile, METH_VARARGS,
     "Closes an open file"},
    {"FlushFile", IOSection_FlushFile, METH_VARARGS,
     "Flushes the file buffer"},
    {"OpenSectionName", IOSection_OpenSectionName, METH_VARARGS,
     "Opens the current section given the name"},
    {"OpenSectionNum", IOSection_OpenSectionNum, METH_VARARGS,
     "Opens the nth section"},
    {"OpenSectionNameNum", IOSection_OpenSectionNameNum, METH_VARARGS,
     "Opens the nth section"},
    {"IncludeSection", IOSection_IncludeSection, METH_VARARGS,
     "Includes a section"},
    {"NewSectionName", IOSection_NewSectionName, METH_VARARGS,
     "Creates a new section given its name"},
    {"NewSectionFile", IOSection_NewSectionFile, METH_VARARGS,
     "Creates a new section in a new file given the section and file name"},
    {"CloseSection", IOSection_CloseSection, METH_VARARGS,
     "Closes a section"},
    {"CountSections", IOSection_CountSections, METH_VARARGS,
     "Counts the total number of sections in the current section"},
    {"CountSectionsName", IOSection_CountSectionsName, METH_VARARGS,
     "Counts the total number of sections in the current section"},
    {"CountVars", IOSection_CountVars, METH_VARARGS,
     "Counts the total number of variables in the current section"},
    {"GetVarName", IOSection_GetVarName, METH_VARARGS,
     "Returns the name of the num'th variable"},
    {"SetUnderscores", IOSection_SetUnderscores, METH_VARARGS,
     "Sets whether to always include underscores or not."},
    {NULL, NULL, 0, NULL}        /* Sentinel */
};

PyMODINIT_FUNC
initIOSection(void)
{
  (void) Py_InitModule("IOSection", IOSectionMethods);
  import_array();
}
