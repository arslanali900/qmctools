#include "NLPP.h"
#include "../IO/FileExpand.h"

void 
TestRead()
{
  IOSectionClass in;
  NLPPClass nlpp;

  cerr << "Na pseudopotential:\n";
  assert (in.OpenFile(ExpandFileName("~/NLPP/Na/CASINO/Na_CASINO.h5")));
  nlpp.Read (in);
  in.CloseFile();
  nlpp.SetupProjectors (4.0, 16.0);

  cerr << "N pseudopotential:\n";
  assert (in.OpenFile(ExpandFileName("~/NLPP/N/CASINO/N_CASINO.h5")));
  nlpp.Read (in);
  in.CloseFile();
  nlpp.SetupProjectors(10.0, 40.0);

  cerr << "B pseudopotential:\n";
  assert (in.OpenFile(ExpandFileName("~/NLPP/B/CASINO/B_CASINO.h5")));
  nlpp.Read (in);
  in.CloseFile();
  nlpp.SetupProjectors(10.0, 40.0);

}

main()
{
  TestRead();
}
