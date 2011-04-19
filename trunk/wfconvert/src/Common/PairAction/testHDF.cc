#include "DavidPAClass.h"

int main(int argc, char** argv)
{
  DavidPAClass pa;
  pa.ReadDavidSquarerFileHDF5("LJ.h5");
  return 0;
}
