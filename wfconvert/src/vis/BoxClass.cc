#include "BoxClass.h"


// Returns true if a segment walks across a box boundary.  In this
// case, it "breaks" the segment into two disjoint pieces, one from
// r1 to wall1 and one from wall2 to r2.  Note that to account for
// segments that cross more than one wall (i.e. a corner-crosser), we
// have to call this 3 times to ensure we get it right.
bool BoxClass::BreakSegment (Vec3 &r1, Vec3 &r2,
			     Vec3 &wall1, Vec3 &wall2)
{
  PutInBox (r1);
  PutInBox (r2);
  Vec3 disp = r2-r1;
  double eps = 1.0e-12;
  for (int dim=0; dim<3; dim++) 
    if (disp[dim]<-0.5*Box[dim]) { // path wraps in + direction
      double d1 = 0.5*Box[dim]-r1[dim];
      double d2 = r2[dim] + 0.5*Box[dim];

      double s1 = d2/(d1+d2);
      double s2 = 1.0-s1;
      Vec3 rshift = r2;
      rshift[dim] += Box[dim];
      wall1 = s1*r1 + s2*rshift;
      wall2 = wall1;
      wall2[dim] -= Box[dim];
      wall1[dim] -= eps;
      wall2[dim] += eps;
      return true;
    }
    else if (disp[dim]>0.5*Box[dim]) {
      double d2 = 0.5*Box[dim] - r2[dim];
      double d1 = r1[dim] + 0.5*Box[dim];
      double s1 = d2/(d1+d2);
      double s2 = 1.0-s1;
      Vec3 rshift = r2;
      rshift[dim] -= Box[dim];
      wall1 = s1*r1 + s2*rshift;
      wall2 = wall1;
      wall2[dim] += Box[dim];
      wall1[dim] += eps;
      wall2[dim] -= eps;
      return true;
    }
  return false;
}


bool BoxClass::BreakSegment (Vec3 &r1, Vec3 &r2,
			     Vec3 &wall1, Vec3 &wall2,
			     int dim)
{
  PutInBox (r1, dim);
  PutInBox (r2, dim);
  Vec3 disp = r2-r1;
  if (disp[dim]<-0.5*Box[dim]) { // path wraps in + direction
    double d1 = 0.5*Box[dim]-r1[dim];
    double d2 = r2[dim] + 0.5*Box[dim];
    
    double s1 = d2/(d1+d2);
    double s2 = 1.0-s1;
    Vec3 rshift = r2;
    rshift[dim] += Box[dim];
    wall1 = s1*r1 + s2*rshift;
    wall2 = wall1;
    wall2[dim] -= Box[dim];
    return true;
  }
  else if (disp[dim]>0.5*Box[dim]) {
    double d2 = 0.5*Box[dim] - r2[dim];
    double d1 = r1[dim] + 0.5*Box[dim];
    double s1 = d2/(d1+d2);
    double s2 = 1.0-s1;
    Vec3 rshift = r2;
    rshift[dim] -= Box[dim];
    wall1 = s1*r1 + s2*rshift;
    wall2 = wall1;
    wall2[dim] += Box[dim];
    return true;
  }
  return false;
}



void BoxClass::PutPathsInBox (vector<OnePath*>& inList)
{
  // To be perfectly correct, we need to repeat this process three
  // times to ensure that segments that cross the box corners get
  // properly broken up into two or three segments.
  for (int dim=0; dim<3; dim++) {
    vector<OnePath*> outList;
    outList.resize(0);
    for (int in=0; in<inList.size(); in++) {
      OnePath *newPath = new OnePath;
      OnePath *oldPath = inList[in];
      bool closed = oldPath->Closed;
      for (int slice=0; slice<(oldPath->Path.size()-1); slice++) {
	Vec3 r1, r2, wall1, wall2;
	r1 = oldPath->Path[slice];
	r2 = oldPath->Path[slice+1];
	if (BreakSegment (r1, r2, wall1, wall2, dim)) {
	  newPath->Path.push_back(r1);
	  newPath->Path.push_back(wall1);
	  newPath->Winding = oldPath->Winding;
	  newPath->Closed = false;
	  outList.push_back(newPath);
	  newPath = new OnePath;
	  newPath->Path.push_back(wall2);
	} 
	else 
	  newPath->Path.push_back(r1);
      }
      // Push last coordinate onto new path
      Vec3 rlast = oldPath->Path[oldPath->Path.size()-1];
      PutInBox(rlast, dim);
      newPath->Path.push_back(rlast);
      // Close the path
//       rlast = oldPath->Path[0];
//       PutInBox(rlast, dim);
//       newPath->Path.push_back(rlast);
      newPath->Winding = oldPath->Winding;
      newPath->Closed = oldPath->Closed;
      outList.push_back (newPath);
      delete oldPath;
    }
    inList.resize(outList.size());
    for (int i=0; i<outList.size(); i++)
      inList[i] = outList[i];
  }
}
