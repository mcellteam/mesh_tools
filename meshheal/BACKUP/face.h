// Author: Justin Kinney
// Date: Sep 2008

#ifndef FACE_H
#define FACE_H 1

class Face
{
public:
  int index;	// Face index
  int v1,v2,v3;	// vertex indices
  Face(char *triplet);
};

#endif    

