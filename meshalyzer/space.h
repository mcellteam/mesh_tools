// Author: Justin Kinney
// Date: Feb 2009

#ifndef SPACE_H
#define SPACE_H 1

#include "meshalyzer.h"

class Space
{
private:
  int num_space[3];
  int num_boxes;	// WHY DO I NEED THIS?
  double space_length;
  vec_b b;		// vector of boxes
  double world[6];      // minX, maxX, minY, maxY, minZ, maxZ of world
public:
  ~Space(void);
  Space (void)
    :num_boxes(0),space_length(0),b()
  {
  }
  void initBoxes(double);
  void clearBoxes(void);
  void deleteBoxes(void);
  void recordFace(vec_b&,Face*);
  void computeBoxesToCheck(Face*,vec_b&);
  void index2Range(int,double[2],char*);
  //int location2Index(double,char*);
  int location2Index(double,const char *);
  int indices2Index(int,int,int);
  void getBoxesFor3DLocations(double[6],vec_b&);
  void getBoxesFor3DIndices(int[6],vec_b&);
  int screenIndex(int,const char*);
  c_b_iterator first_box (void)
  {
    return b.begin();
  }
  c_b_iterator one_past_last_box (void)
  {
    return b.end();
  }
  void setWorld (int i,double a)
  {
    world[i] = a;
  }
  double getWorld (int i)
  {
    return world[i];
  }
  double getSpaceLength (void)
  {
    return space_length;
  }
};

#endif
