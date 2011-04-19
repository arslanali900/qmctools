#ifndef ONE_PATH_H
#define ONE_PATH_H

#include <Common/Blitz.h>
#include <vector>

class OnePath
{
public:
  vector<Vec3> Path;
  vector<Vec3> Color;
  bool Closed, Winding;
};


#endif
