#include "ParserClass.h"

#include <fstream>
#include <cstring>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <cassert>
#include <cstdlib>

streamsize ParserClass::FileSize(string fname)
{
  ifstream infile;
  infile.open(fname.c_str(), ifstream::in);
  if (!infile.is_open())
    return false;
  streampos fileSize = 0;
  infile.seekg(0, ios_base::end);
  fileSize = infile.tellg();
  infile.close();
  return fileSize;
}

bool 
MemParserClass::OpenFile(string fname)
{
  ifstream infile;
  infile.open(fname.c_str(), ifstream::in);
  if (!infile.is_open())
    return false;
  streampos fileSize = 0;
  infile.seekg(fileSize, ios_base::end);
  fileSize = infile.tellg();
  infile.seekg(0, ios_base::beg);
  if (fileSize != (streampos)-1) {
    Buffer.resize(fileSize);
    infile.read(&(Buffer[0]), fileSize);
    infile.close();
  }
  else {
    int fd = open (fname.c_str(), O_RDONLY);
    char ch[100000];
    int len = read (fd, ch, 99999);
    ch[len] = '\0';
    close (fd);
    Buffer = ch;
  }
  Pos = 0;
  return true;
}

void
MemParserClass::CloseFile()
{
  Buffer.clear();
}

bool
MemParserClass::FindToken(string token)
{
  bool found = false;
  int toklen = token.size();
  int tokenPos = Buffer.find(token, Pos);
  if (tokenPos == -1)
    return false;
  Pos = tokenPos+token.size();
  return true;
}


bool 
MemParserClass::ReadInt (int &val)
{
//   int numChars;
//   int success = sscanf (&(Buffer[Pos]), " %d %n", &val, &numChars);
//   if (success) 
//     Pos += numChars;
//   return (success == 1);
  char * endptr;
  val = strtol (&(Buffer[Pos]), &endptr, 10);
  if (endptr == &(Buffer[Pos]))
    return false;
  Pos += endptr - &(Buffer[Pos]);
  return (true);  
}

bool 
MemParserClass::ReadLong (long &val)
{
  char * endptr;
  val = strtol (&(Buffer[Pos]), &endptr, 10);
  if (endptr == &(Buffer[Pos]))
    return false;
  Pos += endptr - &(Buffer[Pos]);
  return (true);  
}


bool 
MemParserClass::ReadDouble (double &val)
{
  char *endptr;
  val =strtod (&(Buffer[Pos]), &endptr);
  if (endptr == &(Buffer[Pos]))
    return false;
  Pos += endptr - &(Buffer[Pos]);
  return (true);
}

void
MemParserClass::SavePos()
{  saved = Pos;  }


void
MemParserClass::RestorePos()
{  Pos = saved;  }

void
FileParserClass::SavePos()
{  saved = Pos;  }

void
FileParserClass::RestorePos()
{  
  Pos = saved;
  Infile.seekg(Pos, ios_base::end);
}


bool
ParserClass::ReadComplex (complex<double> &val)
{
  double re, im;
  if (FindToken ("("))
    if (ReadDouble(re))
      if (FindToken(","))
	if (ReadDouble(im))
	  if (FindToken(")")) {
	    val = complex<double>(re,im);
	    return true;
	  }
  return false;
}

bool isWhiteSpace (char c)
{
  return ((c==' ') || (c=='\t') || (c=='\n') || (c=='\r'));
}

bool
MemParserClass::ReadWord (string &word)
{
  word = "";
  char str[2];
  str[1] = '\0';
  while (isWhiteSpace (Buffer[Pos]) && (Pos<(Buffer.size()-1)))
    Pos++;
  while (!isWhiteSpace(Buffer[Pos]) && (Pos<Buffer.size()-1)) {
    str[0] = Buffer[Pos];
    word.append(str);
    Pos++;
  }
  return true;
}

bool
MemParserClass::ReadLine (string &word)
{
  word = "";
  char str[2];
  str[1] = '\0';
  while ((Buffer[Pos]!='\n') && (Pos<Buffer.size()-1)) {
    str[0] = Buffer[Pos];
    word.append(str);
    Pos++;
  }
  return true;
}

bool
MemParserClass::NextLine ()
{
  while ((Buffer[Pos]!='\n') && (Pos<Buffer.size()-1)) 
    Pos++;
  if (Pos < Buffer.size()) {
    Pos++;
    return true;
  }
  else
    return false;
}


void
MemParserClass::Reset()
{
  Pos = 0;
}



bool
FileParserClass::OpenFile (string fname)
{
  Infile.open(fname.c_str(), ifstream::in);
  if (!Infile.is_open())
    return false;
  FileSize = 0;
  Infile.seekg(FileSize, ios_base::end);
  FileSize = Infile.tellg();
  Infile.seekg((streampos) 0, ios_base::beg);
  Pos = 0;
  return true;
}

void
FileParserClass::CloseFile()
{
  if (!Infile.is_open()) 
    cerr << "Tried to close a FileParserClass that's not open.\n";
  else
    Infile.close();
}

bool
FileParserClass::FindToken(string token)
{
  assert (Infile.is_open());
  char compare[token.size()+1];
  Pos = Infile.tellg();
  bool found = false;
  while (!found) {
    Infile.seekg(Pos, ios_base::beg);
    Infile.get(compare, (streamsize)token.size()+1, '\0');
    if (token == compare) {
      Pos += token.size();
      Infile.seekg(Pos, ios_base::beg);
      return true;
    }
    if (Pos >= FileSize)
      return false;
    Pos+=1;
  }
  return false;
}

bool
FileParserClass::ReadInt (int &val)
{
  Infile >> val;
  return !Infile.fail();
}

bool
FileParserClass::ReadLong(long &val)
{
  Infile >> val;
  return !Infile.fail();
}

bool
FileParserClass::ReadDouble(double &val)
{
  Infile >> val;
  return !Infile.fail();
}

bool
FileParserClass::ReadComplex (complex<double> &val)
{
  double re, im;
  if (FindToken ("("))
    if (ReadDouble(re))
      if (FindToken(","))
	if (ReadDouble(im))
	  if (FindToken(")")) {
	    val = complex<double>(re,im);
	    return true;
	  }
  return false;
}

bool
FileParserClass::ReadWord (string &word)
{
  word = "";
  char ch;
  char str[2];
  str[1] = '\0';
  while (isWhiteSpace(ch = Infile.get()) && !Infile.eof());
  if (Infile.eof())
    return false;
  while (!isWhiteSpace(ch = Infile.get()) && !Infile.eof()) {
    str[0] = ch;
    word.append (str);
  }
  if (isWhiteSpace(ch))
    Infile.unget();
  return true;
}

bool
FileParserClass::ReadLine (string &line)
{
  line = "";
  char ch;
  char str[2];
  str[1] = '\0';
  if (Infile.eof())
    return false;
  while (((ch = Infile.get())!='\n') && !Infile.eof()) {
    str[0] = ch;
    line.append (str);
  }
  return true;
}

bool
FileParserClass::NextLine ()
{
  if (Infile.eof())
    return false;
  while (Infile.get()!='\n' && !Infile.eof());

  return true;
}



void
FileParserClass::Reset()
{
  assert (Infile.is_open());
  Pos = 0;
  Infile.seekg (Pos, ios_base::beg);
}




////////////////////////////////////////////////////////////
//                   FileParserClass2                     //
////////////////////////////////////////////////////////////
bool
FileParserClass2::OpenFile(string fname)
{
  Infile.open(fname.c_str(), ifstream::in);
  if (!Infile.is_open())
    return false;
  FileSize = 0;
  Infile.seekg(FileSize, ios_base::end);
  FileSize = Infile.tellg();
  Infile.seekg((streampos) 0, ios_base::beg);
  Pos = 0;
  ReadChunk (0);
  return true;
}

void
FileParserClass2::Reset()
{
  assert (Infile.is_open());
  Pos = 0;
  Infile.seekg (Pos, ios_base::beg);
  ReadChunk(0);
}


void 
FileParserClass2::CloseFile()
{
  if (!Infile.is_open()) 
    cerr << "Tried to close a FileParserClass that's not open.\n";
  else
    Infile.close();
  Buffer.resize(0);
}



void
FileParserClass2::ReadChunk (long start)
{
  long n = min(MaxBufferSize, FileSize-start);
  if (Buffer.size() != n)
    Buffer.resize(n);
  Infile.seekg(start, ios_base::beg);
  Infile.read(&Buffer[0], n);
  BufferStart = start;
  Pos = start;
}

bool
FileParserClass2::FindToken (string token)
{
  bool fileEnd = false;
  if (Pos < BufferStart)
    ReadChunk (Pos);
  do {
    long tokenPos = Buffer.find (token, Pos-BufferStart);
    if (tokenPos != -1) {
      Pos = tokenPos + BufferStart + token.size();
      return true;
    }
    else if (BufferStart + Buffer.size() >= FileSize) {
      return false;
    }
    else {
      ReadChunk (BufferEnd()- (token.size()+1));
      Pos = BufferStart;
    }
  } while (true);
}

void
FileParserClass2::SavePos()
{  saved = Pos;  }

void
FileParserClass2::RestorePos()
{  
  Pos = saved;
  if (Pos < BufferStart)
    ReadChunk (Pos);
  else if (BufferStart + Buffer.size() < Pos)
    ReadChunk (Pos);
}



bool 
FileParserClass2::ReadInt(int &val)
{
  if ((BufferEnd() - Pos) < 100)
    ReadChunk (Pos);
  char *endptr;
  long offset = Pos - BufferStart;
  val = strtol (&(Buffer[offset]), &endptr, 10);
  if (endptr == &(Buffer[offset]))
    return false;
  Pos += endptr - &(Buffer[offset]);
  return (true);  
}

bool 
FileParserClass2::ReadLong(long &val)
{
  if (BufferEnd() - Pos < 100)
    ReadChunk (Pos);
  char *endptr;
  long offset = Pos - BufferStart;
  val = strtol (&(Buffer[offset]), &endptr, 10);
  if (endptr == &(Buffer[offset]))
    return false;
  Pos += endptr - &(Buffer[offset]);
  return (true);  
}


bool 
FileParserClass2::ReadDouble(double &val)
{
  if (BufferEnd() - Pos < 100)
    ReadChunk (Pos);
  char *endptr;
  long offset = Pos - BufferStart;
  val = strtod (&(Buffer[offset]), &endptr);
  if (endptr == &(Buffer[offset])) {
    cerr << "Couldn't file double.\n";
    return false;
  }
  Pos += endptr - &(Buffer[offset]);
  return (true);  
}

bool
FileParserClass2::ReadWord (string &word)
{
  if (BufferEnd() - Pos < 100)
    ReadChunk (Pos);
   char str[2];
  str[1] = '\0';
  long offset = Pos - BufferStart;
  while (isWhiteSpace (Buffer[offset]) && (offset<(Buffer.size()-1))) {
    offset++; 
    Pos++;
  }

  while (!isWhiteSpace(Buffer[offset]) && (offset<Buffer.size()-1)) {
    str[0] = Buffer[offset];
    word.append(str);
    offset++;
    Pos++;
  }
  return true;
}

bool
FileParserClass2::ReadLine (string &line)
{
  if (BufferEnd() - Pos < 100)
    ReadChunk (Pos);
  char str[2];
  str[1] = '\0';
  long offset = Pos - BufferStart;

  while ((Buffer[offset] != '\n') && (offset<Buffer.size()-1)) {
    str[0] = Buffer[offset];
    line.append(str);
    offset++;
    Pos++;
  }
  return true;
}


bool
FileParserClass2::NextLine ()
{
  if (BufferEnd() - Pos < 100)
    ReadChunk (Pos);
  long offset = Pos - BufferStart;

  while ((Buffer[offset] != '\n') && (offset<Buffer.size()-1)) {
    offset++;
    Pos++;
  }
  return true;
}

