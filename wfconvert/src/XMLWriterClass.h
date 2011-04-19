#ifndef XML_WRITER_CLASS
#define XML_WRITER_CLASS

#include <vector>
#include <string>
#include <libxml/xmlwriter.h>

using namespace std;

class XMLWriterClass
{
private:
  xmlTextWriterPtr Writer;
public:
  bool StartDocument(string fname, string version="",
		     string encoding="", string standalone="");
  bool EndDocument();
  bool StartElement (string name);
  bool EndElement();
  bool FullEndElement();
  bool WriteAttribute (string name, string content);
  bool WriteAttribute (string name, double val, bool scientific=false);
  bool WriteAttribute (string name, int val);
  bool WriteData(vector<double> data);
  bool WriteElement(string name, vector<double> data);
};

#endif
