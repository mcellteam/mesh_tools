#ifndef GROUP_H
#define GROUP_H 1

#include "contour.h"

class Group
{
private:
  Group (Group const &);
  Group & operator = (Group const &);
public:
  int count;	// number of contours in group
  Contour *c;
  Group(void);
  ~Group(void);
  void printOtherContourData(char *str);
  void printContourFacesVertices(char *str);
  void getContours(std::list<Face> & evfsp,Vertex **vert_array,int max_vert,int z_value,std::list<Vertex> & S);
  void addFacesVertices(std::list<Face> & F,std::list<Vertex> & V,int max_vertex,int max_face);
  void gatherThirdDeleteFace(std::list<Face> & F,int val,Vertex **vert_array);
  void convertLateral(std::list<Vertex> & V,double epsilon);
  int faceContainsLaterals(Vertex *v1,Vertex *v2,Face *q,int val,Vertex **vert_array);
  void addFaces(std::list<Face> & F,int v1,int v2, int th,int orient,int &MF); 
  void addVertex(std::list<Vertex> & V,Vertex *EV,int MV); 
  void lateralVertices(std::list<Face> & S,int z_value,Vertex **vert_array); 
  void groupExtraVertices(void);
  void initializeSites(void);
  void orderSites(void);
  void extraVertices (std::list<Vertex> & S);
};

#endif
