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

#include <libxml/parser.h>




void startElementWrapper(void *parserPtr, const xmlChar *name,
			 const xmlChar **attr)
{
  list<XMLattribute> attributes;
  string nameStr = (const char *)name;
  int i=0;
  XMLattribute myAttr;
  if (attr != NULL)
    while (attr[i] != NULL) {
      myAttr.Name = (const char *)(attr[i++]);
      myAttr.Value = (const char *)(attr[i++]);
      attributes.push_back(myAttr);
    }
  XMLparserClass &parser = *((XMLparserClass *) parserPtr);
  parser.startElement(nameStr, attributes);
}


void endElementWrapper(void *parserPtr, const xmlChar *name)
{
  string nameStr = (const char *)name;
  XMLparserClass &parser = *((XMLparserClass *) parserPtr);
  parser.endElement(nameStr);
}


void charactersWrapper(void *parserPtr, const xmlChar *ch, int len)
{
  string str((const char *)ch, len);
  XMLparserClass &parser = *((XMLparserClass *) parserPtr);
  parser.characters(str);
}



void XMLparserClass::SetHandler()
{
  handler.internalSubset        = NULL;
  handler.isStandalone          = NULL;
  handler.hasInternalSubset     = NULL;
  handler.hasExternalSubset     = NULL;
  handler.resolveEntity         = NULL;
  handler.getEntity             = NULL;
  handler.entityDecl            = NULL;
  handler.notationDecl          = NULL;
  handler.attributeDecl         = NULL;
  handler.elementDecl           = NULL;
  handler.unparsedEntityDecl    = NULL;
  handler.setDocumentLocator    = NULL;
  handler.startDocument         = NULL;
  handler.endDocument           = NULL;
  handler.startElement          = startElementWrapper;
  handler.endElement            = endElementWrapper;
  handler.reference             = NULL;
  handler.characters            = charactersWrapper;
  handler.ignorableWhitespace   = NULL;
  handler.processingInstruction = NULL;
  handler.comment               = NULL;
  handler.warning               = NULL;
  handler.error                 = NULL;
  handler.fatalError            = NULL;
}


bool XMLparserClass::ParseFile(string fileName, IOTreeXMLClass *rootNode)
{
  FileName = fileName;
  CurrTree = rootNode;
  // set up the handler functions
  SetHandler();
  // store our filename
  CurrTree->FileName = fileName;
  
  int retval = xmlSAXUserParseFile(&handler, this, fileName.c_str());
  // don't know if this is right.
  return (retval == 0);
}

inline bool IsWhiteSpace(char ch)
{
  return ((ch==' ') || (ch=='\t') || (ch=='\n'));
}

void StripWhiteSpace(string &str)
{
  int i=0;
  while (i<str.length()) {
    if (IsWhiteSpace(str[i]))
      str.erase(i,1);
    else
      i++;
  }
}

void CommaSplit(string &str, list<string> &strList)
{
  int len=0;
  int start=0;
  string tmpString;
  while ((len+start)<str.length()) {
    if (str[len+start] == ',') {
      tmpString = str.substr(start,len);
      strList.push_back(tmpString);
      start+=len+1;
      len=0;
    }
    else
      len++;
  }
  tmpString = str.substr(start,len);
  strList.push_back(tmpString);
}

void TestCommaSplit()
{
  string s = "abc,def, 123";
  list<string> strList;
  CommaSplit(s, strList);
  list<string>::iterator iter = strList.begin();
  while (iter != strList.end()) {
    cerr << "\"" << *iter << "\"\n";
    iter++;
  }
}

void TestStripWhiteSpace()
{
  string s = " abcdef 123456 ";
  StripWhiteSpace(s);
  cerr << s <<endl;
}


inline int StrToInt (string &str)
{
  char *endPtr;
  int d = strtol(str.c_str(), &endPtr,10);
  if (*endPtr!='\0') {
    cerr << "Error reading integer in XML parser.\n";
    abort();
  }
  return (d);
}

inline double StrToDouble(string str)
{
  double val;
  char *endPtr;
  val = strtod(str.c_str(), &endPtr);
  if (*endPtr!='\0') {
    cerr << "Badly formed double in XML parser.\n";
    abort();
  }
  return (val);
}

inline string StrToString(string str)
{
  string val;
  int len = str.length();
  if ((str[0] == '\"') && (str[len-1] == '\"')) {
    val = str;
    val.erase(len-1,1);
    val.erase(0,1);
    return (val);
  }
  else {
    cerr << "Badly formed string in XML parser.\n";
    abort();
  } 
}

void Lower (string &str)
{
  for (int i=0; i<str.length(); i++)
    str[i] = tolower(str[i]);
}

inline bool StrToBool (string str)
{
  Lower (str);
  if (str == "true")
    return (true);
  else if (str == "false")
    return (false);
  else {
    cerr << "Badly formed bool in XML parser.\n";
    abort();
  }
}



void ParseDim(string dimStr, blitz::Array<int,1> &dimensions)
{
  StripWhiteSpace(dimStr);
  list<string> strList;
  CommaSplit(dimStr, strList);
  list<string>::iterator iter = strList.begin();
  int numDim = 0;
  while (iter != strList.end()) {
    iter++;
    numDim++;
  }
  dimensions.resize(numDim);
  iter = strList.begin();
  for (int i=0; i<numDim; i++)
    dimensions(i) = StrToInt(*iter);
}



VarXMLClass *NewXMLVar (AtomicType newType, int ndim,
			    blitz::Array<int,1> dims)
{
  if (ndim == 0) {
    if (newType == DOUBLE_TYPE)
      return new VarXMLdouble0Class;
    else if (newType == INT_TYPE)
      return new VarXMLint0Class;
    else if (newType == STRING_TYPE)
      return new VarXMLstring0Class;
    else if (newType == BOOL_TYPE)
      return new VarXMLbool0Class;
  }
  else if (ndim == 1) {
    if (newType == DOUBLE_TYPE)
      {
	VarXMLdouble1Class *newVar = new VarXMLdouble1Class;
	newVar->Value.resize(dims(0));
	return newVar;
      }
    else if (newType == INT_TYPE)
      {
	VarXMLint1Class *newVar = new VarXMLint1Class;
	newVar->Value.resize(dims(0));
	return newVar;
      }
    else if (newType == STRING_TYPE)
      {
	VarXMLstring1Class *newVar = new VarXMLstring1Class;
	newVar->Value.resize(dims(0));
	return newVar;
      }
    else if (newType == BOOL_TYPE)
      {
	VarXMLbool1Class *newVar = new VarXMLbool1Class;
	newVar->Value.resize(dims(0));
	return newVar;
      }
  }
  else if (ndim == 2) {
    if (newType == DOUBLE_TYPE)
      {
	VarXMLdouble2Class *newVar = new VarXMLdouble2Class;
	newVar->Value.resize(dims(0), dims(1));
	return newVar;
      }
    else if (newType == INT_TYPE)
      {
	VarXMLint2Class *newVar = new VarXMLint2Class;
	newVar->Value.resize(dims(0), dims(1));
	return newVar;
      }
    else if (newType == STRING_TYPE)
      {
	VarXMLstring2Class *newVar = new VarXMLstring2Class;
	newVar->Value.resize(dims(0), dims(1));
	return newVar;
      }
    else if (newType == BOOL_TYPE)
      {
	VarXMLbool2Class *newVar = new VarXMLbool2Class;
	newVar->Value.resize(dims(0), dims(1));
	return newVar;
      }
  }  
  else if (ndim == 3) {
    if (newType == DOUBLE_TYPE)
      {
	VarXMLdouble3Class *newVar = new VarXMLdouble3Class;
	newVar->Value.resize(dims(0), dims(1), dims(2));
	return newVar;
      }
    else if (newType == INT_TYPE)
      {
	VarXMLint3Class *newVar = new VarXMLint3Class;
	newVar->Value.resize(dims(0), dims(1), dims(2));
	return newVar;
      }
    else if (newType == STRING_TYPE)
      {
	VarXMLstring3Class *newVar = new VarXMLstring3Class;
	newVar->Value.resize(dims(0), dims(1), dims(2));
	return newVar;
      }
    else if (newType == BOOL_TYPE)
      {
	VarXMLbool3Class *newVar = new VarXMLbool3Class;
	newVar->Value.resize(dims(0), dims(1), dims(2));
	return newVar;
      }
  }  
  else if (ndim == 4) {
    if (newType == DOUBLE_TYPE)
      {
	VarXMLdouble4Class *newVar = new VarXMLdouble4Class;
	newVar->Value.resize(dims(0), dims(1), dims(2), dims(3));
	return newVar;
      }
    else if (newType == INT_TYPE)
      {
	VarXMLint4Class *newVar = new VarXMLint4Class;
	newVar->Value.resize(dims(0), dims(1), dims(2), dims(3));
	return newVar;
      }
    else if (newType == STRING_TYPE)
      {
	VarXMLstring4Class *newVar = new VarXMLstring4Class;
	newVar->Value.resize(dims(0), dims(1), dims(2), dims(3));
	return newVar;
      }
    else if (newType == BOOL_TYPE)
      {
	VarXMLbool4Class *newVar = new VarXMLbool4Class;
	newVar->Value.resize(dims(0), dims(1), dims(2), dims(3));
	return newVar;
      }
  }  
}


AtomicType StrToType (string &str)
{
  if (str == "double")
    return DOUBLE_TYPE;
  else if (str == "int")
    return INT_TYPE;
  else if (str == "string")
    return STRING_TYPE;
  else if (str == "bool")
    return BOOL_TYPE;
  else { 
    cerr << "Unknown Atomic type " << str << endl;
    abort();
  }
}

void XMLparserClass::startElement(string &name, list<XMLattribute> &attributes)
{
  // First, check to see if we have a section or variable
  ElementIsSection = true;
  // Push a new character buffer onto the stack
  string s;
  charBuffers.push(s);

  string typeStr, dimStr, file;
  list<XMLattribute>::iterator iter = attributes.begin();
  while (iter != attributes.end()) {
    cerr << "Atrribute = " << iter->Name << endl;
    if (iter->Name == "type") {
      typeStr = iter->Value;
      ElementIsSection = false;
    }
    else if (iter->Name == "dim")
      dimStr = iter->Value;
    else if (iter->Name == "file")
      file = iter->Value;
    else {
      cerr << "Unrecognized attribute " << iter->Name 
	   << " in file " << FileName << ".  Exitting\n";
      abort();
    }
    iter++;
  }

  if (ElementIsSection) {
    if (file != "") { // We're including another file
      // DO SOMETHING!!!!
      cerr << "Section (" << name << ") is including "
	   << file << ".\n";
    }
    else {  // This is a new section
      IOTreeXMLClass *newTree = new IOTreeXMLClass;
      newTree->Name = name;
      newTree->Parent = CurrTree;
      CurrTree->SectionList.push_back(newTree);
      CurrTree = newTree;
      cerr << "New Section (" << name << ")\n";
    }
  }
  else { // This is a variable
    int ndim;
    AtomicType type;
    blitz::Array<int,1> dim;
    type = StrToType(typeStr);
    if (dimStr != "") { // This is an array
      ParseDim(dimStr, dim);
      ndim = dim.size();
      cerr << "blitz::Array variable " << name << " of type "
	   << typeStr << " and dimensions (" << dimStr
	   << ").\n";
    }
    else { // This is a scalar variable
      ndim = 0;
      cerr << "Scalar variable " << name << " of type " 
	   << typeStr << ".\n";
    }
    CurrVar = NewXMLVar(type, ndim, dim);
    CurrVar->Name = name;
    CurrVar->Type = type;
    CurrVar->Dim = ndim;
  }


//   cerr << "Element name:  " << name << endl;
//   cerr << "Attributes:  " << endl;
//   list<XMLattribute>::iterator iter = attributes.begin();
//   while (iter != attributes.end()) {
//     cerr << iter->Name << "=\"" << iter->Value << "\"\n";
//     iter++;
//   }
}

template <class T>
void ReadArrayData(list<string> numbers, blitz::Array<T,1> valArray)
{
  


}


void VarXMLdouble0Class::ReadVals(list<string> &vals)
{
  if (vals.size()!=1) {
    cerr << "Wrong number of doubles in XML parser for variable "
	 << Name << endl;
    abort();
  }
  Value = StrToDouble (vals.front());
}

void VarXMLdouble1Class::ReadVals(list<string> &vals)
{
  if (vals.size()!=Value.size()) {
    cerr << "Wrong number of doubles in XML parser for variable "
	 << Name << endl;
    cerr << "Expected " << Value.size() << ", but got " 
	 << vals.size() << ".\n";
    abort();
  }
  list<string>::iterator iter = vals.begin();
  for (int i=0; i<Value.extent(0); i++) {
    Value(i) = StrToDouble(*iter);
    iter++;
  }
}

void VarXMLdouble2Class::ReadVals(list<string> &vals)
{
  if (vals.size()!=Value.size()) {
    cerr << "Wrong number of doubles in XML parser for variable "
	 << Name << endl;
    cerr << "Expected " << Value.size() << ", but got " 
	 << vals.size() << ".\n";
    abort();
  }
  list<string>::iterator iter = vals.begin();
  for (int i=0; i<Value.extent(0); i++) 
    for (int j=0; j<Value.extent(1); j++) {
      Value(i,j) = StrToDouble(*iter);
      iter++;
    }
}

void VarXMLdouble3Class::ReadVals(list<string> &vals)
{
  if (vals.size()!=Value.size()) {
    cerr << "Wrong number of ints in XML parser for variable "
	 << Name << endl;
    cerr << "Expected " << Value.size() << ", but got " 
	 << vals.size() << ".\n";
    abort();
  }
  list<string>::iterator iter = vals.begin();
  for (int i=0; i<Value.extent(0); i++) 
    for (int j=0; j<Value.extent(1); j++) 
      for (int k=0; k<Value.extent(2); k++){
	Value(i,j,k) = StrToDouble(*iter);
	iter++;
      }
}

void VarXMLdouble4Class::ReadVals(list<string> &vals)
{
  if (vals.size()!=Value.size()) {
    cerr << "Wrong number of ints in XML parser for variable "
	 << Name << endl;
    cerr << "Expected " << Value.size() << ", but got " 
	 << vals.size() << ".\n";
    abort();
  }
  list<string>::iterator iter = vals.begin();
  for (int i=0; i<Value.extent(0); i++) 
    for (int j=0; j<Value.extent(1); j++) 
      for (int k=0; k<Value.extent(2); k++)
	for (int l=0; l<Value.extent(3); l++) {
	  Value(i,j,k,l) = StrToDouble(*iter);
	  iter++;
	}
}


void VarXMLint0Class::ReadVals(list<string> &vals)
{
  if (vals.size()!=1) {
    cerr << "Wrong number of ints in XML parser for variable "
	 << Name << endl;
    abort();
  }
  Value = StrToInt (vals.front());
}

void VarXMLint1Class::ReadVals(list<string> &vals)
{
  if (vals.size()!=Value.size()) {
    cerr << "Wrong number of ints in XML parser for variable "
	 << Name << endl;
    cerr << "Expected " << Value.size() << ", but got " 
	 << vals.size() << ".\n";
    abort();
  }
  list<string>::iterator iter = vals.begin();
  for (int i=0; i<Value.extent(0); i++) {
    Value(i) = StrToInt(*iter);
    iter++;
  }
}

void VarXMLint2Class::ReadVals(list<string> &vals)
{
  if (vals.size()!=Value.size()) {
    cerr << "Wrong number of ints in XML parser for variable "
	 << Name << endl;
    cerr << "Expected " << Value.size() << ", but got " 
	 << vals.size() << ".\n";
    abort();
  }
  list<string>::iterator iter = vals.begin();
  for (int i=0; i<Value.extent(0); i++) 
    for (int j=0; j<Value.extent(1); j++) {
      Value(i,j) = StrToInt(*iter);
      iter++;
    }
}

void VarXMLint3Class::ReadVals(list<string> &vals)
{
  if (vals.size()!=Value.size()) {
    cerr << "Wrong number of ints in XML parser for variable "
	 << Name << endl;
    cerr << "Expected " << Value.size() << ", but got " 
	 << vals.size() << ".\n";
    abort();
  }
  list<string>::iterator iter = vals.begin();
  for (int i=0; i<Value.extent(0); i++) 
    for (int j=0; j<Value.extent(1); j++) 
      for (int k=0; k<Value.extent(2); k++){
	Value(i,j,k) = StrToInt(*iter);
	iter++;
      }
}

void VarXMLint4Class::ReadVals(list<string> &vals)
{
  if (vals.size()!=Value.size()) {
    cerr << "Wrong number of ints in XML parser for variable "
	 << Name << endl;
    cerr << "Expected " << Value.size() << ", but got " 
	 << vals.size() << ".\n";
    abort();
  }
  list<string>::iterator iter = vals.begin();
  for (int i=0; i<Value.extent(0); i++) 
    for (int j=0; j<Value.extent(1); j++) 
      for (int k=0; k<Value.extent(2); k++)
	for (int l=0; l<Value.extent(3); l++){
	  Value(i,j,k,l) = StrToInt(*iter);
	  iter++;
	}
}



void VarXMLstring0Class::ReadVals(list<string> &vals)
{
  if (vals.size()!=1) {
    cerr << "Wrong number of ints in XML parser for variable "
	 << Name << endl;
    abort();
  }
  Value = StrToString (vals.front());
}

void VarXMLstring1Class::ReadVals(list<string> &vals)
{
  if (vals.size()!=Value.size()) {
    cerr << "Wrong number of ints in XML parser for variable "
	 << Name << endl;
    cerr << "Expected " << Value.size() << ", but got " 
	 << vals.size() << ".\n";
    abort();
  }
  list<string>::iterator iter = vals.begin();
  for (int i=0; i<Value.extent(0); i++) {
    Value(i) = StrToString(*iter);
    iter++;
  }
}

void VarXMLstring2Class::ReadVals(list<string> &vals)
{
  if (vals.size()!=Value.size()) {
    cerr << "Wrong number of ints in XML parser for variable "
	 << Name << endl;
    cerr << "Expected " << Value.size() << ", but got " 
	 << vals.size() << ".\n";
    abort();
  }
  list<string>::iterator iter = vals.begin();
  for (int i=0; i<Value.extent(0); i++) 
    for (int j=0; j<Value.extent(1); j++) {
      Value(i,j) = StrToString(*iter);
      iter++;
    }
}

void VarXMLstring3Class::ReadVals(list<string> &vals)
{
  if (vals.size()!=Value.size()) {
    cerr << "Wrong number of ints in XML parser for variable "
	 << Name << endl;
    cerr << "Expected " << Value.size() << ", but got " 
	 << vals.size() << ".\n";
    abort();
  }
  list<string>::iterator iter = vals.begin();
  for (int i=0; i<Value.extent(0); i++) 
    for (int j=0; j<Value.extent(1); j++) 
      for (int k=0; k<Value.extent(2); k++){
	Value(i,j,k) = StrToString(*iter);
	iter++;
      }
}


void VarXMLstring4Class::ReadVals(list<string> &vals)
{
  if (vals.size()!=Value.size()) {
    cerr << "Wrong number of ints in XML parser for variable "
	 << Name << endl;
    cerr << "Expected " << Value.size() << ", but got " 
	 << vals.size() << ".\n";
    abort();
  }
  list<string>::iterator iter = vals.begin();
  for (int i=0; i<Value.extent(0); i++) 
    for (int j=0; j<Value.extent(1); j++) 
      for (int k=0; k<Value.extent(2); k++)
	for (int l=0; l<Value.extent(3); l++){
	  Value(i,j,k,l) = StrToString(*iter);
	  iter++;
	}
}



void VarXMLbool0Class::ReadVals(list<string> &vals)
{
  if (vals.size()!=1) {
    cerr << "Wrong number of bools in XML parser for variable "
	 << Name << endl;
    abort();
  }
  Value = StrToBool (vals.front());
}

void VarXMLbool1Class::ReadVals(list<string> &vals)
{
  if (vals.size()!=Value.size()) {
    cerr << "Wrong number of bools in XML parser for variable "
	 << Name << endl;
    cerr << "Expected " << Value.size() << ", but got " 
	 << vals.size() << ".\n";
    abort();
  }
  list<string>::iterator iter = vals.begin();
  for (int i=0; i<Value.extent(0); i++) {
    Value(i) = StrToBool(*iter);
    iter++;
  }
}

void VarXMLbool2Class::ReadVals(list<string> &vals)
{
  if (vals.size()!=Value.size()) {
    cerr << "Wrong number of bools in XML parser for variable "
	 << Name << endl;
    cerr << "Expected " << Value.size() << ", but got " 
	 << vals.size() << ".\n";
    abort();
  }
  list<string>::iterator iter = vals.begin();
  for (int i=0; i<Value.extent(0); i++) 
    for (int j=0; j<Value.extent(1); j++) {
      Value(i,j) = StrToBool(*iter);
      iter++;
    }
}

void VarXMLbool3Class::ReadVals(list<string> &vals)
{
  if (vals.size()!=Value.size()) {
    cerr << "Wrong number of bools in XML parser for variable "
	 << Name << endl;
    cerr << "Expected " << Value.size() << ", but got " 
	 << vals.size() << ".\n";
    abort();
  }
  list<string>::iterator iter = vals.begin();
  for (int i=0; i<Value.extent(0); i++) 
    for (int j=0; j<Value.extent(1); j++) 
      for (int k=0; k<Value.extent(2); k++) {
	Value(i,j,k) = StrToBool(*iter);
	iter++;
      }
}

void VarXMLbool4Class::ReadVals(list<string> &vals)
{
  if (vals.size()!=Value.size()) {
    cerr << "Wrong number of ints in XML parser for variable "
	 << Name << endl;
    cerr << "Expected " << Value.size() << ", but got " 
	 << vals.size() << ".\n";
    abort();
  }
  list<string>::iterator iter = vals.begin();
  for (int i=0; i<Value.extent(0); i++) 
    for (int j=0; j<Value.extent(1); j++) 
      for (int k=0; k<Value.extent(2); k++)
	for (int l=0; l<Value.extent(3); l++){
	  Value(i,j,k,l) = StrToBool(*iter);
	  iter++;
	}
}






void XMLparserClass::endElement(string &name)
{
  string buffer=charBuffers.top();
  charBuffers.pop();
  if (ElementIsSection)
    CurrTree = CurrTree->Parent;
  else {
    // Parse variable character data here.
    list<string> numbers;
    StripWhiteSpace(buffer);
    CommaSplit(buffer,numbers);
    CurrVar->ReadVals(numbers);
    CurrTree->VarList.push_back(CurrVar);
    
    ElementIsSection = true;
  }
  cerr << "End element " << name << endl;
}

void XMLparserClass::characters(string &newChars)
{
  charBuffers.top() += newChars;
}



bool IOTreeXMLClass::OpenFile(string fileName, string myName,
			      IOTreeClass *parent)
{
  Parent = parent;
  Name = myName;
  XMLparserClass parser;
  return parser.ParseFile (fileName, this);
}




/// Simply prints 3*num spaces
inline void XMLPrintIndent(int num)
{
  for (int counter=0;counter<num*2;counter++)
    cout<<' ';
}

/// Simply prints 3*num spaces
inline void XMLPrintIndent(int num,ofstream &outFile)
{
  for (int counter=0;counter<num*2;counter++)
    outFile<<' ';
}


/// Prints an indented hierarchy of sections and variable names to
/// cout. 
void IOTreeXMLClass::PrintTree(int indentNum)
{
  XMLPrintIndent(indentNum);
  cout<<"Section: "<<Name<<endl;
  list<VarClass*>::iterator varIter=VarList.begin();
  while (varIter!=VarList.end()){
    XMLPrintIndent(indentNum+1);
    cout<<"Variable: "<<(*varIter)->Name<<" "<<endl;
    varIter++;
  }
  list<IOTreeClass*>::iterator secIter=SectionList.begin();
  while (secIter!=SectionList.end()){
    //    cout<<"Section: "<<(*secIter)->Name<<endl;
    (*secIter)->PrintTree(indentNum+1);
    secIter++;
  }
}

/// Calls PrintTree(0)
void IOTreeXMLClass::PrintTree()
{
  PrintTree(0);
}






IOTreeClass* IOTreeXMLClass::NewSection(string name)
{
  IOTreeClass* tempSection=new IOTreeXMLClass();
  tempSection->Name=name;
  tempSection->Parent=this;
  tempSection->MyNumber=CurrSecNum;
  CurrSecNum++;
  SectionList.push_back(tempSection);
  MarkModified();
  return tempSection;
}

void IOTreeXMLClass::IncludeSection(IOTreeClass *newSection)
{
  newSection->MyNumber=CurrSecNum++;
  SectionList.push_back(newSection);
  MarkModified();
}



bool IOTreeXMLClass::NewFile (string fileName,
			      string mySectionName,
			      IOTreeClass *parent)
{
  FileName=fileName;
  Parent=parent;
  Name=mySectionName;
}


void IOTreeXMLClass::WriteSection(ofstream &outFile,int indentNum)
{
  list<VarClass*>::iterator varIter=VarList.begin();
  while (varIter!=VarList.end()){
    XMLPrintIndent(indentNum,outFile);
    ((VarXMLClass*)(*varIter))->Print(outFile);    
    varIter++;
  }
  list<IOTreeClass*>::iterator secIter=SectionList.begin();
  while (secIter!=SectionList.end()){
    if ((*secIter)->FileName==""){
      XMLPrintIndent(indentNum,outFile);
      outFile<<"Section ("<<(*secIter)->Name<<")\n";
      XMLPrintIndent(indentNum,outFile);
      outFile<<"{\n";
      ((IOTreeXMLClass*)(*secIter))->WriteSection(outFile,indentNum+1);
      XMLPrintIndent(indentNum,outFile);
      outFile<<"}\n\n";
    }
    else {
      XMLPrintIndent(indentNum,outFile);
      outFile<<"Section ("<<(*secIter)->Name<<", \"";
      outFile<<(*secIter)->FileName<<"\");"<<endl;
      (*secIter)->FlushFile();
    }
    secIter++;
  }
}



void IOTreeXMLClass::FlushFile()
{
  ofstream outfile;
  if ((FileName!="") && IsModified){
    outfile.open(FileName.c_str());
    WriteSection(outfile,0);
  }

  list<IOTreeClass*>::iterator iter = SectionList.begin();
  while (iter != SectionList.end()) {
    (*iter)->FlushFile();
    iter++;
  }
}

/// CloseFile recursively destroys the tree of data we constructed.
void IOTreeXMLClass::CloseFile()
{  
  // First, free all the variables in the list
  while (!VarList.empty()) {
    delete(VarList.front());
    VarList.pop_front();
  }
  
  // Now, call all closes recursively and delete all sections
  while (!SectionList.empty())
    {
      SectionList.front()->CloseFile();
      delete SectionList.front();
      SectionList.pop_front();
    }
}



//////////////////////////////////////////////////////////////////////
//                             ReadInto's                           //
//////////////////////////////////////////////////////////////////////

bool VarXMLClass::ReadInto (double &val)
{ ComplainReadInto(); return false; }
bool VarXMLClass::ReadInto (int &val)
{ ComplainReadInto(); return false; }
bool VarXMLClass::ReadInto (string &val)
{ ComplainReadInto(); return false; }
bool VarXMLClass::ReadInto (bool &val)
{ ComplainReadInto(); return false; }
bool VarXMLClass::ReadInto (blitz::Array<double,1> &val)
{ ComplainReadInto(); return false; }
bool VarXMLClass::ReadInto (blitz::Array<double,2> &val)
{ ComplainReadInto(); return false; }
bool VarXMLClass::ReadInto (blitz::Array<double,3> &val)
{ ComplainReadInto(); return false; }
bool VarXMLClass::ReadInto (blitz::Array<double,4> &val)
{ ComplainReadInto(); return false; }
bool VarXMLClass::ReadInto (blitz::Array<int,1> &val)
{ ComplainReadInto(); return false; }
bool VarXMLClass::ReadInto (blitz::Array<int,2> &val)
{ ComplainReadInto(); return false; }
bool VarXMLClass::ReadInto (blitz::Array<int,3> &val)
{ ComplainReadInto(); return false; }
bool VarXMLClass::ReadInto (blitz::Array<int,4> &val)
{ ComplainReadInto(); return false; }
bool VarXMLClass::ReadInto (blitz::Array<string,1> &val)
{ ComplainReadInto(); return false; }
bool VarXMLClass::ReadInto (blitz::Array<string,2> &val)
{ ComplainReadInto(); return false; }
bool VarXMLClass::ReadInto (blitz::Array<string,3> &val)
{ ComplainReadInto(); return false; }
bool VarXMLClass::ReadInto (blitz::Array<string,4> &val)
{ ComplainReadInto(); return false; }
bool VarXMLClass::ReadInto (blitz::Array<bool,1> &val)
{ ComplainReadInto(); return false; }
bool VarXMLClass::ReadInto (blitz::Array<bool,2> &val)
{ ComplainReadInto(); return false; }
bool VarXMLClass::ReadInto (blitz::Array<bool,3> &val)
{ ComplainReadInto(); return false; }
bool VarXMLClass::ReadInto (blitz::Array<bool,4> &val)
{ ComplainReadInto(); return false; }

bool VarXMLClass::Append (double val)
{ ComplainAppend(); return false; }
bool VarXMLClass::Append (int val)
{ ComplainAppend(); return false; }
bool VarXMLClass::Append (string val)
{ ComplainAppend(); return false; }
bool VarXMLClass::Append (bool val)
{ ComplainAppend(); return false; }
bool VarXMLClass::Append (blitz::Array<double,1> &val)
{ ComplainAppend(); return false; }
bool VarXMLClass::Append (blitz::Array<double,2> &val)
{ ComplainAppend(); return false; }
bool VarXMLClass::Append (blitz::Array<double,3> &val)
{ ComplainAppend(); return false; }
bool VarXMLClass::Append (blitz::Array<int,1> &val)
{ ComplainAppend(); return false; }
bool VarXMLClass::Append (blitz::Array<int,2> &val)
{ ComplainAppend(); return false; }
bool VarXMLClass::Append (blitz::Array<int,3> &val)
{ ComplainAppend(); return false; }
bool VarXMLClass::Append (blitz::Array<string,1> &val)
{ ComplainAppend(); return false; }
bool VarXMLClass::Append (blitz::Array<string,2> &val)
{ ComplainAppend(); return false; }
bool VarXMLClass::Append (blitz::Array<string,3> &val)
{ ComplainAppend(); return false; }
bool VarXMLClass::Append (blitz::Array<bool,1> &val)
{ ComplainAppend(); return false; }
bool VarXMLClass::Append (blitz::Array<bool,2> &val)
{ ComplainAppend(); return false; }
bool VarXMLClass::Append (blitz::Array<bool,3> &val)
{ ComplainAppend(); return false; }


bool VarXMLdouble0Class::ReadInto (double &val)
{ val = Value; return true;}
bool VarXMLint0Class::ReadInto (int &val)
{  val = Value; return true; }
bool VarXMLstring0Class::ReadInto (string &val)
{  val = Value; return true; }
bool VarXMLbool0Class::ReadInto (bool &val)
{ val = Value; return true; }
bool VarXMLdouble1Class::ReadInto (blitz::Array<double,1> &val)
{ val.resize(Value.extent(0)); val = Value; return true; }
bool VarXMLdouble2Class::ReadInto (blitz::Array<double,2> &val)
{ val.resize(Value.extent(0),Value.extent(1)); val = Value; return true; }
bool VarXMLdouble3Class::ReadInto (blitz::Array<double,3> &val)
{ val.resize(Value.extent(0),Value.extent(1),Value.extent(2)); 
  val = Value; return true; }
bool VarXMLdouble4Class::ReadInto (blitz::Array<double,4> &val)
{ val.resize(Value.extent(0),Value.extent(1),Value.extent(2),Value.extent(3)); 
  val = Value; return true; }
bool VarXMLint1Class::ReadInto (blitz::Array<int,1> &val)
{ val.resize(Value.extent(0)); val = Value; return true; }
bool VarXMLint2Class::ReadInto (blitz::Array<int,2> &val)
{ val.resize(Value.extent(0),Value.extent(1)); val = Value; return true; }
bool VarXMLint3Class::ReadInto (blitz::Array<int,3> &val)
{ val.resize(Value.extent(0),Value.extent(1),Value.extent(2)); 
  val = Value; return true; }
bool VarXMLint4Class::ReadInto (blitz::Array<int,4> &val)
{ val.resize(Value.extent(0),Value.extent(1),Value.extent(2),Value.extent(3)); 
  val = Value; return true; }
bool VarXMLstring1Class::ReadInto (blitz::Array<string,1> &val)
{ val.resize(Value.extent(0)); val = Value; return true; }
bool VarXMLstring2Class::ReadInto (blitz::Array<string,2> &val)
{ val.resize(Value.extent(0),Value.extent(1)); val = Value; return true; }
bool VarXMLstring3Class::ReadInto (blitz::Array<string,3> &val)
{ val.resize(Value.extent(0),Value.extent(1),Value.extent(2)); 
 val = Value; return true; }
bool VarXMLstring4Class::ReadInto (blitz::Array<string,4> &val)
{ val.resize(Value.extent(0),Value.extent(1),Value.extent(2),Value.extent(3)); 
 val = Value; return true; }
bool VarXMLbool1Class::ReadInto (blitz::Array<bool,1> &val)
{ val.resize(Value.extent(0)); val = Value; return true; }
bool VarXMLbool2Class::ReadInto (blitz::Array<bool,2> &val)
{ val.resize(Value.extent(0),Value.extent(1)); val = Value; return true; }
bool VarXMLbool3Class::ReadInto (blitz::Array<bool,3> &val)
{ val.resize(Value.extent(0),Value.extent(1),Value.extent(2)); 
 val = Value; return true; }
bool VarXMLbool4Class::ReadInto (blitz::Array<bool,4> &val)
{ val.resize(Value.extent(0),Value.extent(1),Value.extent(2),Value.extent(3)); 
 val = Value; return true; }


/*****************************************************************
 *                         Variable  Appends                     *
 *****************************************************************/
bool VarXMLdouble1Class::Append (double val)
{
  int n = Value.extent(0);
  Value.resizeAndPreserve(n+1);
  Value(n) = val;
  return(true);
}
bool VarXMLint1Class::Append (int val)
{
  int n = Value.extent(0);
  Value.resizeAndPreserve(n+1);
  Value(n) = val;
  return(true);
}
bool VarXMLstring1Class::Append (string val)
{
  int n = Value.extent(0);
  Value.resizeAndPreserve(n+1);
  Value(n) = val;
  return(true);
}
bool VarXMLbool1Class::Append (bool val)
{
  int n = Value.extent(0);
  Value.resizeAndPreserve(n+1);
  Value(n) = val;
  return(true);
}
bool VarXMLdouble2Class::Append (blitz::Array<double,1> &val)
{
  int n = Value.extent(0);  int m = Value.extent(1);
  assert(val.extent(0) == m);
  Value.resizeAndPreserve(n+1,m);
  Value(n,blitz::Range::all()) = val;
  return(true);
}
bool VarXMLdouble3Class::Append (blitz::Array<double,2> &val)
{
  int n=Value.extent(0); int m=Value.extent(1); int o=Value.extent(2);
  assert(val.extent(0) == m);
  assert(val.extent(1) == o);
  Value.resizeAndPreserve(n+1,m,o);
  Value(n,blitz::Range::all(),blitz::Range::all()) = val;
  return(true);
}
bool VarXMLdouble4Class::Append (blitz::Array<double,3> &val)
{
  int n=Value.extent(0); int m=Value.extent(1); 
  int o=Value.extent(2); int p=Value.extent(3);
  assert(val.extent(0) == m);
  assert(val.extent(1) == o);
  assert(val.extent(2) == p);
  Value.resizeAndPreserve(n+1,m,o,p);
  Value(n,blitz::Range::all(),blitz::Range::all(),blitz::Range::all()) = val;
  return(true);
}
bool VarXMLint2Class::Append (blitz::Array<int,1> &val)
{
  int n = Value.extent(0);  int m = Value.extent(1);
  assert(val.extent(0) == m);
  Value.resizeAndPreserve(n+1,m);
  Value(n,blitz::Range::all()) = val;
  return(true);
}
bool VarXMLint3Class::Append (blitz::Array<int,2> &val)
{
  int n=Value.extent(0); int m=Value.extent(1); int o=Value.extent(2);
  assert(val.extent(0) == m);
  assert(val.extent(1) == o);
  Value.resizeAndPreserve(n+1,m,o);
  Value(n,blitz::Range::all(),blitz::Range::all()) = val;
  return(true);
}
bool VarXMLint4Class::Append (blitz::Array<int,3> &val)
{
  int n=Value.extent(0); int m=Value.extent(1); 
  int o=Value.extent(2); int p=Value.extent(3);
  assert(val.extent(0) == m);
  assert(val.extent(1) == o);
  assert(val.extent(2) == p);
  Value.resizeAndPreserve(n+1,m,o,p);
  Value(n,blitz::Range::all(),blitz::Range::all(),blitz::Range::all()) = val;
  return(true);
}
bool VarXMLstring2Class::Append (blitz::Array<string,1> &val)
{
  int n = Value.extent(0);  int m = Value.extent(1);
  assert(val.extent(0) == m);
  Value.resizeAndPreserve(n+1,m);
  Value(n,blitz::Range::all()) = val;
  return(true);
}
bool VarXMLstring3Class::Append (blitz::Array<string,2> &val)
{
  int n=Value.extent(0); int m=Value.extent(1); int o=Value.extent(2);
  assert(val.extent(0) == m);
  assert(val.extent(1) == o);
  Value.resizeAndPreserve(n+1,m,o);
  Value(n,blitz::Range::all(),blitz::Range::all()) = val;
  return(true);
}
bool VarXMLstring4Class::Append (blitz::Array<string,3> &val)
{
  int n=Value.extent(0); int m=Value.extent(1); 
  int o=Value.extent(2); int p=Value.extent(3);
  assert(val.extent(0) == m);
  assert(val.extent(1) == o);
  assert(val.extent(1) == p);
  Value.resizeAndPreserve(n+1,m,o,p);
  Value(n,blitz::Range::all(),blitz::Range::all(),blitz::Range::all()) = val;
  return(true);
}
bool VarXMLbool2Class::Append (blitz::Array<bool,1> &val)
{
  int n = Value.extent(0);  int m = Value.extent(1);
  assert(val.extent(0) == m);
  Value.resizeAndPreserve(n+1,m);
  Value(n,blitz::Range::all()) = val;
  return(true);
}
bool VarXMLbool3Class::Append (blitz::Array<bool,2> &val)
{
  int n=Value.extent(0); int m=Value.extent(1); int o=Value.extent(2);
  assert(val.extent(0) == m);
  assert(val.extent(1) == o);
  Value.resizeAndPreserve(n+1,m,o);
  Value(n,blitz::Range::all(),blitz::Range::all()) = val;
  return(true);
}
bool VarXMLbool4Class::Append (blitz::Array<bool,3> &val)
{
  int n=Value.extent(0); int m=Value.extent(1); 
  int o=Value.extent(2); int p=Value.extent(3);
  assert(val.extent(0) == m);
  assert(val.extent(1) == o);
  assert(val.extent(2) == p);
  Value.resizeAndPreserve(n+1,m,o,p);
  Value(n,blitz::Range::all(),blitz::Range::all(),blitz::Range::all()) = val;
  return(true);
}

void VarXMLdouble0Class::Print (ofstream &outFile)
{ outFile << "double " << Name << " = " << Value << ";\n";}
void VarXMLint0Class::Print (ofstream &outFile)
{ outFile << "int " << Name << " = " << Value << ";\n";} 
void VarXMLstring0Class::Print (ofstream &outFile)
{ outFile << "string " << Name << " = " << "\"" << Value << "\";\n";} 
void VarXMLbool0Class::Print (ofstream &outFile)
{ outFile << "bool " << Name <<" = "<< (Value ? "true;\n" : "false;\n"); }
void VarXMLdouble1Class::Print (ofstream &outFile)
{ 
  outFile << "blitz::Array<double,1> " << Name 
	  << "(" << Value.extent(0) << ") = [";
  for (int i=0; i<(Value.extent(0)-1); i++)
    outFile << Value(i) << ", ";
  outFile << Value(Value.extent(0)-1) << "];\n";
}
void VarXMLdouble2Class::Print (ofstream &outFile)
{
  outFile << "blitz::Array<double,2> " << Name
	  << "(" << Value.extent(0) << "," 
	  << Value.extent(1) << ") = [";
  for (int i=0; i<(Value.extent(0)); i++)
    for (int j=0; j<(Value.extent(1)); j++)
      if ((i < (Value.extent(0)-1)) || (j<(Value.extent(1)-1)))
	outFile << Value(i,j) << ", ";
  outFile << Value(Value.extent(0)-1,Value.extent(1)-1) << "];\n";
}
void VarXMLdouble3Class::Print (ofstream &outFile)
{
  outFile << "blitz::Array<double,3> " << Name
	  << "(" << Value.extent(0) << "," 
	  << Value.extent(1) << "," << Value.extent(2) << ") = [";
  for (int i=0; i<(Value.extent(0)); i++)
    for (int j=0; j<(Value.extent(1)); j++)
      for (int k=0; k<(Value.extent(2)); k++)
      if ((i < (Value.extent(0)-1)) || (j<(Value.extent(1)-1))
	  || (k<(Value.extent(2)-1)))
	outFile << Value(i,j,k) << ", ";
  outFile << Value(Value.extent(0)-1,Value.extent(1)-1,Value.extent(2)-1) 
	  << "];\n";
}
void VarXMLdouble4Class::Print (ofstream &outFile)
{
  outFile << "blitz::Array<double,4> " << Name
	  << "(" << Value.extent(0) << "," 
	  << Value.extent(1) << "," 
	  << Value.extent(2) << ","
	  << Value.extent(3) <<") = [";
  for (int i=0; i<(Value.extent(0)); i++)
    for (int j=0; j<(Value.extent(1)); j++)
      for (int k=0; k<(Value.extent(2)); k++)
	for (int l=0; l<(Value.extent(3)); l++)
      if ((i < (Value.extent(0)-1)) || (j<(Value.extent(1)-1))
	  || (k<(Value.extent(2)-1)) || (l<(Value.extent(3)-1)))
	outFile << Value(i,j,k,l) << ", ";
  outFile << Value(Value.extent(0)-1,Value.extent(1)-1,
		   Value.extent(2)-1,Value.extent(3)-1) 
	  << "];\n";
}
void VarXMLint1Class::Print (ofstream &outFile)
{
  outFile << "blitz::Array<int,1> " << Name 
	  << "(" << Value.extent(0) << ") = [";
  for (int i=0; i<(Value.extent(0)-1); i++)
    outFile << Value(i) << ", ";
  outFile << Value(Value.extent(0)-1) << "];\n";
}
void VarXMLint2Class::Print (ofstream &outFile)
{
  outFile << "blitz::Array<int,2> " << Name 
	  << "(" << Value.extent(0) << "," 
	  << Value.extent(1) << ") = [";
  for (int i=0; i<(Value.extent(0)); i++)
    for (int j=0; j<(Value.extent(1)); j++)
      if ((i < (Value.extent(0)-1)) || (j<(Value.extent(1)-1)))
	outFile << Value(i,j) << ", ";
  outFile << Value(Value.extent(0)-1,Value.extent(1)-1) << "];\n";
}
void VarXMLint3Class::Print (ofstream &outFile)
{
  outFile << "blitz::Array<int,3> " << Name
	  << "(" << Value.extent(0) << "," 
	  << Value.extent(1) << "," << Value.extent(2) << ") = [";
  for (int i=0; i<(Value.extent(0)); i++)
    for (int j=0; j<(Value.extent(1)); j++)
      for (int k=0; k<(Value.extent(2)); k++)
      if ((i < (Value.extent(0)-1)) || (j<(Value.extent(1)-1))
	  || (k<(Value.extent(2)-1)))
	outFile << Value(i,j,k) << ", ";
  outFile << Value(Value.extent(0)-1,Value.extent(1)-1,Value.extent(2)-1) 
	  << "];\n";
}
void VarXMLint4Class::Print (ofstream &outFile)
{
  outFile << "blitz::Array<int,4> " << Name
	  << "(" 
	  << Value.extent(0) << "," 
	  << Value.extent(1) << "," 
	  << Value.extent(2) << ","
	  << Value.extent(3) <<") = [";
  for (int i=0; i<(Value.extent(0)); i++)
    for (int j=0; j<(Value.extent(1)); j++)
      for (int k=0; k<(Value.extent(2)); k++)
	for (int l=0; l<(Value.extent(3)); l++)
	  if ((i < (Value.extent(0)-1)) || (j<(Value.extent(1)-1))
	      || (k<(Value.extent(2)-1)) || (l<(Value.extent(3)-1)))
	    outFile << Value(i,j,k,l) << ", ";
  outFile << Value(Value.extent(0)-1,Value.extent(1)-1,
		   Value.extent(2)-1,Value.extent(3)-1) 
	  << "];\n";
}
void VarXMLstring1Class::Print (ofstream &outFile)
{
  outFile << "blitz::Array<string,1> " << Name 
	  << "(" << Value.extent(0) << ") = [";
  for (int i=0; i<(Value.extent(0)-1); i++)
    outFile << "\"" << Value(i) << "\"" << ", ";
  outFile << "\"" << Value(Value.extent(0)-1) << "\"];\n";
}
void VarXMLstring2Class::Print (ofstream &outFile)
{
  outFile << "blitz::Array<string,2> " << Name 
	  << "(" << Value.extent(0) << "," 
	  << Value.extent(1) << ") = [";
  for (int i=0; i<(Value.extent(0)); i++)
    for (int j=0; j<(Value.extent(1)); j++)
      if ((i < (Value.extent(0)-1)) || (j<(Value.extent(1)-1)))
	outFile << "\"" << Value(i,j) << "\", ";
  outFile << Value(Value.extent(0)-1,Value.extent(1)-1) << "];\n";
}
void VarXMLstring3Class::Print (ofstream &outFile)
{
  outFile << "blitz::Array<string,3> " << Name
	  << "(" << Value.extent(0) << "," 
	  << Value.extent(1) << "," << Value.extent(2) << ") = [";
  for (int i=0; i<(Value.extent(0)); i++)
    for (int j=0; j<(Value.extent(1)); j++)
      for (int k=0; k<(Value.extent(2)); k++)
      if ((i < (Value.extent(0)-1)) || (j<(Value.extent(1)-1))
	  || (k<(Value.extent(2)-1)))
	outFile << "\"" << Value(i,j,k) << "\", ";
  outFile << "\"" 
	  << Value(Value.extent(0)-1,Value.extent(1)-1,Value.extent(2)-1) 
	  << "\"];\n";
}
void VarXMLstring4Class::Print (ofstream &outFile)
{
  outFile << "blitz::Array<string,4> " << Name
	  << "(" 
	  << Value.extent(0) << "," 
	  << Value.extent(1) << "," 
	  << Value.extent(2) << ","
	  << Value.extent(3) <<") = [";
  for (int i=0; i<(Value.extent(0)); i++)
    for (int j=0; j<(Value.extent(1)); j++)
      for (int k=0; k<(Value.extent(2)); k++)
	for (int l=0; l<(Value.extent(3)); l++)
	  if ((i < (Value.extent(0)-1)) || (j<(Value.extent(1)-1))
	      || (k<(Value.extent(2)-1)) || (l<(Value.extent(3)-1)))
	    outFile << "\"" << Value(i,j,k,l) << "\", ";
  outFile << "\"" << Value(Value.extent(0)-1,Value.extent(1)-1,
		   Value.extent(2)-1,Value.extent(3)-1) 
	  << "\"];\n";
}

void VarXMLbool1Class::Print (ofstream &outFile)
{
  outFile << "blitz::Array<bool,1> " << Name 
	  << "(" << Value.extent(0) << ") = [";
  for (int i=0; i<(Value.extent(0)-1); i++)
    outFile << (Value(i) ? "true" : "false") << ", ";
  outFile << (Value(Value.extent(0)-1) ? "true" : "false") << "];\n";
}
void VarXMLbool2Class::Print (ofstream &outFile)
{
  outFile << "blitz::Array<bool,2> " << Name 
	  << "(" << Value.extent(0) << "," 
	  << Value.extent(1) << ") = [";
  for (int i=0; i<(Value.extent(0)); i++)
    for (int j=0; j<(Value.extent(1)); j++)
      if ((i < (Value.extent(0)-1)) || (j<(Value.extent(1)-1)))
	outFile << (Value(i,j) ? "true" : "false") << ", ";
  outFile << 
    (Value(Value.extent(0)-1,Value.extent(1)-1) ? "true" : "false") << "];\n";
}
void VarXMLbool3Class::Print (ofstream &outFile)
{
  outFile << "blitz::Array<bool,3> " << Name
	  << "(" << Value.extent(0) << "," 
	  << Value.extent(1) << "," << Value.extent(2) << ") = [";
  for (int i=0; i<(Value.extent(0)); i++)
    for (int j=0; j<(Value.extent(1)); j++)
      for (int k=0; k<(Value.extent(2)); k++)
      if ((i < (Value.extent(0)-1)) || (j<(Value.extent(1)-1))
	  || (k<(Value.extent(2)-1)))
	outFile << (Value(i,j,k) ? "true" : "false") << ", ";
  outFile << (Value(Value.extent(0)-1,Value.extent(1)-1,Value.extent(2)-1) ?
	      "true" : "false") << "];\n";
}
void VarXMLbool4Class::Print (ofstream &outFile)
{
  outFile << "blitz::Array<bool,4> " << Name
	  << "(" 
	  << Value.extent(0) << "," 
	  << Value.extent(1) << "," 
	  << Value.extent(2) << ","
	  << Value.extent(3) << ") = [";
  for (int i=0; i<(Value.extent(0)); i++)
    for (int j=0; j<(Value.extent(1)); j++)
      for (int k=0; k<(Value.extent(2)); k++)
	for (int l=0; l<(Value.extent(3)); l++)
	  if ((i < (Value.extent(0)-1)) || (j<(Value.extent(1)-1))
	      || (k<(Value.extent(2)-1)) || (l<(Value.extent(3)-1)))
	    outFile << (Value(i,j,k,l) ? "true" : "false") << ", ";
  outFile << (Value(Value.extent(0)-1,Value.extent(1)-1,Value.extent(2)-1) ?
	      "true" : "false") << "];\n";
}




//////////////////////////////////////////////////////////////////////
//                             Append's                             //
//////////////////////////////////////////////////////////////////////


void IOTreeXMLClass::WriteVar(string name, double val)
{
  VarXMLdouble0Class *newVar = new VarXMLdouble0Class;
  newVar->Name=name;
  newVar->Value=val;
  VarList.push_back(newVar);
  MarkModified();
}
void IOTreeXMLClass::WriteVar(string name, blitz::Array<double,1> &val)
{
  VarXMLdouble1Class *newVar = new VarXMLdouble1Class;
  newVar->Name=name;
  newVar->Value.resize(val.extent(0));
  newVar->Value=val;
  VarList.push_back(newVar);
  MarkModified();
}
void IOTreeXMLClass::WriteVar(string name, blitz::Array<double,2> &val)
{
  VarXMLdouble2Class *newVar = new VarXMLdouble2Class;
  newVar->Name=name;
  newVar->Value.resize(val.extent(0), val.extent(1));
  newVar->Value=val;
  VarList.push_back(newVar);
  MarkModified();
}
void IOTreeXMLClass::WriteVar(string name, blitz::Array<double,3> &val)
{
  VarXMLdouble3Class *newVar = new VarXMLdouble3Class;
  newVar->Name=name;
  newVar->Value.resize(val.extent(0), val.extent(1), val.extent(2));
  newVar->Value=val;
  VarList.push_back(newVar);
  MarkModified();
}
void IOTreeXMLClass::WriteVar(string name, blitz::Array<double,4> &val)
{
  VarXMLdouble4Class *newVar = new VarXMLdouble4Class;
  newVar->Name=name;
  newVar->Value.resize(val.extent(0), val.extent(1), 
		       val.extent(2), val.extent(3));
  newVar->Value=val;
  VarList.push_back(newVar);
  MarkModified();
}

void IOTreeXMLClass::WriteVar(string name, int val)
{
  VarXMLint0Class *newVar = new VarXMLint0Class;
  newVar->Name=name;
  newVar->Value=val;
  VarList.push_back(newVar);
  MarkModified();
}
void IOTreeXMLClass::WriteVar(string name, blitz::Array<int,1> &val)
{
  VarXMLint1Class *newVar = new VarXMLint1Class;
  newVar->Name=name;
  newVar->Value.resize(val.extent(0));
  newVar->Value=val;
  VarList.push_back(newVar);
  MarkModified();
}
void IOTreeXMLClass::WriteVar(string name, blitz::Array<int,2> &val)
{
  VarXMLint2Class *newVar = new VarXMLint2Class;
  newVar->Name=name;
  newVar->Value.resize(val.extent(0), val.extent(1));
  newVar->Value=val;
  VarList.push_back(newVar);
  MarkModified();
}
void IOTreeXMLClass::WriteVar(string name, blitz::Array<int,3> &val)
{
  VarXMLint3Class *newVar = new VarXMLint3Class;
  newVar->Name=name;
  newVar->Value.resize(val.extent(0), val.extent(1), val.extent(2));
  newVar->Value=val;
  VarList.push_back(newVar);
  MarkModified();
}
void IOTreeXMLClass::WriteVar(string name, blitz::Array<int,4> &val)
{
  VarXMLint4Class *newVar = new VarXMLint4Class;
  newVar->Name=name;
  newVar->Value.resize(val.extent(0), val.extent(1), 
		       val.extent(2), val.extent(3));
  newVar->Value=val;
  VarList.push_back(newVar);
  MarkModified();
}

void IOTreeXMLClass::WriteVar(string name, string val)
{
  VarXMLstring0Class *newVar = new VarXMLstring0Class;
  newVar->Name=name;
  newVar->Value=val;
  VarList.push_back(newVar);
  MarkModified();
}
void IOTreeXMLClass::WriteVar(string name, blitz::Array<string,1> &val)
{
  VarXMLstring1Class *newVar = new VarXMLstring1Class;
  newVar->Name=name;
  newVar->Value.resize(val.extent(0));
  newVar->Value=val;
  VarList.push_back(newVar);
  MarkModified();
}
void IOTreeXMLClass::WriteVar(string name, blitz::Array<string,2> &val)
{
  VarXMLstring2Class *newVar = new VarXMLstring2Class;
  newVar->Name=name;
  newVar->Value.resize(val.extent(0), val.extent(1));
  newVar->Value=val;
  VarList.push_back(newVar);
  MarkModified();
}
void IOTreeXMLClass::WriteVar(string name, blitz::Array<string,3> &val)
{
  VarXMLstring3Class *newVar = new VarXMLstring3Class;
  newVar->Name=name;
  newVar->Value.resize(val.extent(0), val.extent(1), val.extent(2));
  newVar->Value=val;
  VarList.push_back(newVar);
  MarkModified();
}
void IOTreeXMLClass::WriteVar(string name, blitz::Array<string,4> &val)
{
  VarXMLstring4Class *newVar = new VarXMLstring4Class;
  newVar->Name=name;
  newVar->Value.resize(val.extent(0), val.extent(1), 
		       val.extent(2), val.extent(3));
  newVar->Value=val;
  VarList.push_back(newVar);
  MarkModified();
}
 
void IOTreeXMLClass::WriteVar(string name, bool val)
{
  VarXMLbool0Class *newVar = new VarXMLbool0Class;
  newVar->Name=name;
  newVar->Value=val;
  VarList.push_back(newVar);
  MarkModified();
}
void IOTreeXMLClass::WriteVar(string name, blitz::Array<bool,1> &val)
{
  VarXMLbool1Class *newVar = new VarXMLbool1Class;
  newVar->Name=name;
  newVar->Value.resize(val.extent(0));
  newVar->Value=val;
  VarList.push_back(newVar);
  MarkModified();
}
void IOTreeXMLClass::WriteVar(string name, blitz::Array<bool,2> &val)
{
  VarXMLbool2Class *newVar = new VarXMLbool2Class;
  newVar->Name=name;
  newVar->Value.resize(val.extent(0), val.extent(1));
  newVar->Value=val;
  VarList.push_back(newVar);
  MarkModified();
}
void IOTreeXMLClass::WriteVar(string name, blitz::Array<bool,3> &val)
{
  VarXMLbool3Class *newVar = new VarXMLbool3Class;
  newVar->Name=name;
  newVar->Value.resize(val.extent(0), val.extent(1), val.extent(2));
  newVar->Value=val;
  VarList.push_back(newVar);
  MarkModified();
}
void IOTreeXMLClass::WriteVar(string name, blitz::Array<bool,4> &val)
{
  VarXMLbool4Class *newVar = new VarXMLbool4Class;
  newVar->Name=name;
  newVar->Value.resize(val.extent(0), val.extent(1), 
		       val.extent(2), val.extent(3));
  newVar->Value=val;
  VarList.push_back(newVar);
  MarkModified();
}
