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

#ifndef INPUT_OUTPUT_XML_H
#define INPUT_OUTPUT_XML_H
#include "InputOutputBase.h"
#include <libxml/parser.h>
#include <iostream>
#include <stack>
#include <string>
#include <list>



class VarXMLClass : public VarClass
{  
protected:
  inline void ComplainReadInto()
  { 
    cerr << "Trying to read into wrong type in ReadInto.\n";
    exit(1);
  }
  inline void ComplainAppend()
  {
    cerr << "Trying to Append to a variable of wrong type.\n";
  }

public:
  /// Hack to hold different types.  We cast to to whatever type is
  /// indicated by the Type and Dim variables in the parent class.

  /// The ReadIntos copy their contents into val if the type
  /// that's passed to them matches the type that's actually stored in
  /// it's void *Value, as indicated by Type and Dim.
  virtual bool ReadInto (double &val);
  virtual bool ReadInto (int &val);
  virtual bool ReadInto (string &val);
  virtual bool ReadInto (bool &val);
  virtual bool ReadInto (blitz::Array<double,1> &v);
  virtual bool ReadInto (blitz::Array<double,2> &v);
  virtual bool ReadInto (blitz::Array<double,3> &v);
  virtual bool ReadInto (blitz::Array<double,4> &v);
  virtual bool ReadInto (blitz::Array<int,1> &v);
  virtual bool ReadInto (blitz::Array<int,2> &v);
  virtual bool ReadInto (blitz::Array<int,3> &v);
  virtual bool ReadInto (blitz::Array<int,4> &v);
  virtual bool ReadInto (blitz::Array<string,1> &v);
  virtual bool ReadInto (blitz::Array<string,2> &v);
  virtual bool ReadInto (blitz::Array<string,3> &v);
  virtual bool ReadInto (blitz::Array<string,4> &v);
  virtual bool ReadInto (blitz::Array<bool,1> &v);
  virtual bool ReadInto (blitz::Array<bool,2> &v);
  virtual bool ReadInto (blitz::Array<bool,3> &v);
  virtual bool ReadInto (blitz::Array<bool,4> &v);

  virtual bool Append (double val);
  virtual bool Append (blitz::Array<double,1> &val);
  virtual bool Append (blitz::Array<double,2> &val);
  virtual bool Append (blitz::Array<double,3> &val);
  virtual bool Append (int val);
  virtual bool Append (blitz::Array<int,1> &val);
  virtual bool Append (blitz::Array<int,2> &val);
  virtual bool Append (blitz::Array<int,3> &val);
  virtual bool Append (string val);
  virtual bool Append (blitz::Array<string,1> &val);
  virtual bool Append (blitz::Array<string,2> &val);
  virtual bool Append (blitz::Array<string,3> &val);
  virtual bool Append (bool val);
  virtual bool Append (blitz::Array<bool,1> &val);
  virtual bool Append (blitz::Array<bool,2> &val);
  virtual bool Append (blitz::Array<bool,3> &val);
  
  virtual void Print(ofstream &outFile) = 0;

  virtual void ReadVals(list<string> &vals) = 0;
};

class VarXMLdouble0Class : public VarXMLClass
{
public:
  double Value;
  bool ReadInto (double &val);
  void Print(ofstream &outFile);
  void ReadVals(list<string> &vals);
};

class VarXMLdouble1Class : public VarXMLClass
{
public:
  blitz::Array<double,1> Value;
  bool ReadInto (blitz::Array<double,1> &val);
  bool Append (double val);
  void Print(ofstream &outFile);
  void ReadVals(list<string> &vals);
};

class VarXMLdouble2Class : public VarXMLClass
{
public:
  blitz::Array<double,2> Value;
  bool ReadInto (blitz::Array<double,2> &val);
  bool Append (blitz::Array<double,1> &val);
  void Print(ofstream &outFile);
  void ReadVals(list<string> &vals);
};

class VarXMLdouble3Class : public VarXMLClass
{
public:
  blitz::Array<double,3> Value;
  bool ReadInto (blitz::Array<double,3> &val);
  bool Append (blitz::Array<double,2> &val);
  void Print(ofstream &outFile);
  void ReadVals(list<string> &vals);
};

class VarXMLdouble4Class : public VarXMLClass
{
public:
  blitz::Array<double,4> Value;
  bool ReadInto (blitz::Array<double,4> &val);
  bool Append (blitz::Array<double,3> &val);
  void Print(ofstream &outFile);
  void ReadVals(list<string> &vals);
};


class VarXMLint0Class : public VarXMLClass
{
public:
  int Value;
  bool ReadInto (int &val);
  void Print(ofstream &outFile);
  void ReadVals(list<string> &vals);
};

class VarXMLint1Class : public VarXMLClass
{
public:
  blitz::Array<int,1> Value;
  bool ReadInto (blitz::Array<int,1> &val);
  bool Append (int val);
  void Print(ofstream &outFile);
  void ReadVals(list<string> &vals);
};

class VarXMLint2Class : public VarXMLClass
{
public:
  blitz::Array<int,2> Value;
  bool ReadInto (blitz::Array<int,2> &val);
  bool Append (blitz::Array<int,1> &val);
  void Print(ofstream &outFile);
  void ReadVals(list<string> &vals);
};

class VarXMLint3Class : public VarXMLClass
{
public:
  blitz::Array<int,3> Value;
  bool ReadInto (blitz::Array<int,3> &val);
  bool Append (blitz::Array<int,2> &val);
  void Print(ofstream &outFile);
  void ReadVals(list<string> &vals);
};


class VarXMLint4Class : public VarXMLClass
{
public:
  blitz::Array<int,4> Value;
  bool ReadInto (blitz::Array<int,4> &val);
  bool Append (blitz::Array<int,3> &val);
  void Print(ofstream &outFile);
  void ReadVals(list<string> &vals);
};

 
class VarXMLstring0Class : public VarXMLClass
{
public:
  string Value;
  bool ReadInto (string &val);
  void Print(ofstream &outFile);
  void ReadVals(list<string> &vals);
};

class VarXMLstring1Class : public VarXMLClass
{
public:
  blitz::Array<string,1> Value;
  bool ReadInto (blitz::Array<string,1> &val);
  bool Append (string val);
  void Print(ofstream &outFile);
  void ReadVals(list<string> &vals);
};

class VarXMLstring2Class : public VarXMLClass
{
public:
  blitz::Array<string,2> Value;
  bool ReadInto (blitz::Array<string,2> &val);
  bool Append (blitz::Array<string,1> &val);
  void Print(ofstream &outFile);
  void ReadVals(list<string> &vals);
};

class VarXMLstring3Class : public VarXMLClass
{
public:
  blitz::Array<string,3> Value;
  bool ReadInto (blitz::Array<string,3> &val);
  bool Append (blitz::Array<string,2> &val);
  void Print(ofstream &outFile);
  void ReadVals(list<string> &vals);
};

class VarXMLstring4Class : public VarXMLClass
{
public:
  blitz::Array<string,4> Value;
  bool ReadInto (blitz::Array<string,4> &val);
  bool Append (blitz::Array<string,3> &val);
  void Print(ofstream &outFile);
  void ReadVals(list<string> &vals);
};

 

class VarXMLbool0Class : public VarXMLClass
{
public:
  bool Value;
  bool ReadInto (bool &val);
  void Print(ofstream &outFile);
  void ReadVals(list<string> &vals);
};

class VarXMLbool1Class : public VarXMLClass
{
public:
  blitz::Array<bool,1> Value;
  bool ReadInto (blitz::Array<bool,1> &val);
  bool Append (bool val);
  void Print(ofstream &outFile);
  void ReadVals(list<string> &vals);
};

class VarXMLbool2Class : public VarXMLClass
{
public:
  blitz::Array<bool,2> Value;
  bool ReadInto (blitz::Array<bool,2> &val);
  bool Append (blitz::Array<bool,1> &val);
  void Print(ofstream &outFile);
  void ReadVals(list<string> &vals);
};

class VarXMLbool3Class : public VarXMLClass
{
public:
  blitz::Array<bool,3> Value;
  bool ReadInto (blitz::Array<bool,3> &val);
  bool Append (blitz::Array<bool,2> &val);
  void Print(ofstream &outFile);
  void ReadVals(list<string> &vals);
};

class VarXMLbool4Class : public VarXMLClass
{
public:
  blitz::Array<bool,4> Value;
  bool ReadInto (blitz::Array<bool,4> &val);
  bool Append (blitz::Array<bool,3> &val);
  void Print(ofstream &outFile);
  void ReadVals(list<string> &vals);
};


class XMLattribute
{
public:
  string Name;
  string Value;
};


/// This is the XML specialization of IOTreeClass for XML text
/// files.  It's syntax is as follows:
/// <SectionName>
///   <x type="double"> 3 </x>
///   <y type="int" dim="3"> 1, 2, 3 </y>
///   <z type="int" dim="2,2,1"> 1, 2, 3, 4 </z>
///   <Species file="species1.h5/> <!-- Include syntax -->
/// </SectionName>                             
///
class IOTreeXMLClass : public IOTreeClass
{
  /// Reads a text file into a buffer eliminating c++ and c-style
  /// comments.  
  void ParseFile(string fileName);
  bool ReadSection (IOTreeClass *parent, string name,
		    list<TokenClass>::iterator &iter,
		    list<TokenClass> &tokenList,
		    bool wantEndBrace);
 public:
  void AddSection (IOTreeClass *newSection);
  void AddVar (VarClass *newVar);
  void WriteSection(ofstream &outFile,int indent);
  /// Print an indented tree of section variable names.
  void PrintTree(int level);
  /// Same thing, just calls above with level 0;
  void PrintTree();

  IOTreeClass* NewSection(string name);
  void IncludeSection (IOTreeClass *);
  /// Takes the name of a file to read, the name of my section and a
  /// pointer to my parent.  Reads the file into a tree of
  /// IOTreeClass's.
  bool OpenFile (string filename, string myName, 
		 IOTreeClass *parent);
  bool NewFile (string fileName, string mySectionName,
		IOTreeClass *parent);
  /// Do any file handling necessary and delete the whole tree of data.
  void CloseFile();
  void FlushFile();

  void WriteVar(string name, double val);
  void WriteVar(string name, blitz::Array<double,1> &val);
  void WriteVar(string name, blitz::Array<double,2> &val);
  void WriteVar(string name, blitz::Array<double,3> &val);
  void WriteVar(string name, blitz::Array<double,4> &val);

  void WriteVar(string name, int val);
  void WriteVar(string name, blitz::Array<int,1> &val);
  void WriteVar(string name, blitz::Array<int,2> &val);
  void WriteVar(string name, blitz::Array<int,3> &val);
  void WriteVar(string name, blitz::Array<int,4> &val);

  void WriteVar(string name, bool val);
  void WriteVar(string name, blitz::Array<bool,1> &val);
  void WriteVar(string name, blitz::Array<bool,2> &val);
  void WriteVar(string name, blitz::Array<bool,3> &val);
  void WriteVar(string name, blitz::Array<bool,4> &val);

  void WriteVar(string name, string val);
  void WriteVar(string name, blitz::Array<string,1> &val);
  void WriteVar(string name, blitz::Array<string,2> &val);
  void WriteVar(string name, blitz::Array<string,3> &val);
  void WriteVar(string name, blitz::Array<string,4> &val);
  IOTreeXMLClass()
  { IsModified = false; }
};

class XMLparserClass
{
private:
  IOTreeClass *CurrTree;
  VarXMLClass *CurrVar;
  xmlSAXHandler handler;
  stack<string> charBuffers;
  void SetHandler();
  bool ElementIsSection;
  string FileName;
public:
  bool ParseFile (string fileName, IOTreeXMLClass *rootNode);
  void startElement(string &name, list<XMLattribute> &attributes);
  void endElement(string &name);
  void characters(string &newChars);
};



#endif
