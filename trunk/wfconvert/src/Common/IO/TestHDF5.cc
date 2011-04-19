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

#include "InputOutput.h"

using namespace blitz;

void TestHDF5Bool()
{
  IOSectionClass out;
  out.NewFile("BoolTest.h5");
  out.WriteVar ("MyBool", true);
  Array<bool,1> boolArray(5), boolArray2;
  boolArray(0) = true;
  boolArray(1) = false;
  boolArray(2) = true;
  boolArray(3) = false;
  boolArray(4) = true;
  out.WriteVar ("BoolArray", boolArray);
  VarClass *arrayPtr = out.GetVarPtr("BoolArray");
  arrayPtr->Append (true);
  out.CloseFile();

  bool isTrue;
  IOSectionClass in;
  in.OpenFile ("BoolTest.h5");
  assert (in.ReadVar ("MyBool", isTrue));
  cerr << "MyBool = " << (isTrue ? "true" : "false") << endl;
  assert (in.ReadVar ("BoolArray", boolArray2));
  cerr << "BoolArray2 = " << boolArray2 << endl;
}

main()
{
  TestHDF5Bool();


}
