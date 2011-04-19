#include <XMLWriterClass.h>
#include <sstream>
#include <iomanip>

bool
XMLWriterClass::StartDocument(string fname, string version,
			      string encoding, string standalone)
{
  const char *xversion, *xencoding, *xstandalone;
  xversion    = (version=="")    ? NULL : version.c_str();
  xencoding   = (encoding=="")   ? NULL : encoding.c_str();
  xstandalone = (standalone=="") ? NULL : standalone.c_str();
  Writer = xmlNewTextWriterFilename(fname.c_str(), 0);
  xmlTextWriterStartDocument (Writer, xversion, xencoding, xstandalone);
  xmlTextWriterSetIndent (Writer, 1);
  return (Writer != NULL);
}

bool
XMLWriterClass::EndDocument()
{
  bool success = (Writer != NULL);
  success = (xmlTextWriterEndDocument (Writer) != 0);
  xmlFreeTextWriter(Writer);
  return success;
}

bool
XMLWriterClass::StartElement(string name)
{
  return (xmlTextWriterStartElement(Writer, (xmlChar*)name.c_str()) != 0);
}

bool
XMLWriterClass::EndElement()
{
  //  xmlTextWriterWriteRaw (Writer, (const xmlChar*)"\n");
  bool success =  (xmlTextWriterEndElement(Writer) != 0);
  //  xmlTextWriterWriteRaw (Writer, (const xmlChar*)"\n");
  return success;
}

bool
XMLWriterClass::FullEndElement()
{
  xmlTextWriterWriteRaw (Writer, (const xmlChar*)"\n");
  bool success =  (xmlTextWriterFullEndElement(Writer) != 0);
  xmlTextWriterWriteRaw (Writer, (const xmlChar*)"\n");
  return success;
}

bool
XMLWriterClass::WriteAttribute(string name, string content)
{
  int rc = xmlTextWriterWriteAttribute
    (Writer, (xmlChar*)name.c_str(), (xmlChar*)content.c_str());
  return rc != 0;
}

bool
XMLWriterClass::WriteAttribute(string name, double val, bool scientific)
{
  stringstream content;
  if (scientific)
    content.setf(ios_base::scientific, ios_base::floatfield);
  else
    content.setf(ios_base::fixed, ios_base::floatfield);
  content << setprecision(14) << val;
  return WriteAttribute (name, content.str());
}

bool
XMLWriterClass::WriteAttribute(string name, int val)
{
  stringstream content;
  content << val;
  return WriteAttribute (name, content.str());
}

bool 
XMLWriterClass::WriteData(vector<double> data)
{
  stringstream content;
  content.setf(ios_base::scientific, ios_base::floatfield);
  content << setprecision(14);
  for (int i=0; i<data.size(); i++) {
    if ((i % 3) == 0)
      content << "\n";
    content.width(24);
    content << data[i];
  }
  content << "\n";
  xmlTextWriterWriteRaw(Writer, (xmlChar*)content.str().c_str());
}

bool 
XMLWriterClass::WriteElement(string name, vector<double> data)
{
  stringstream content;
  content.setf(ios_base::scientific, ios_base::floatfield);
  content << setprecision(14);
  for (int i=0; i<data.size(); i++) {
    if ((i % 3) == 0)
      content << "\n";
    content.width(24);
    content << data[i];
  }
  content << "\n";
  xmlTextWriterWriteElement(Writer, (xmlChar*)name.c_str(), 
			    (xmlChar*)content.str().c_str());
}

