// Author: Justin Kinney
// Date: Feb 2009

#ifndef CONTAINER_H
#define CONTAINER_H 1

#include <vector>

#include "point.h"
#include "stats.h"
#include "vertex.h"

class Container
{
private:
  Container &  operator =            (Container const &);
  Container                          (Container const &);
public:
  Vertex *eg;		            // example vertex
  vec_o o;			    //vector of object pointers
  std::vector<std::string> files;   // array of input file names
  int num_files;		    // number of input files
  int pairs[3][2];
  Stats area;
  Stats aspect_ratio;
  Stats edge_length;
  Stats edge_angle;
  Stats adjacent_face;
  ///////////
  ~Container(void);
  Container(void);
  ///////////
  void update(Object*);
  void clear(void);
  int getScore(void);
  void printBatch(Controls&);
  void scanDir(const char*);
  //	Object* processFile(char[1024]);
  Object* processFile(std::string);
  //	void scanFile(Object*,char*);
  void scanFile(Object*,std::string);
  void createEdges(void);
  void addVertexPointers(void);
  void findVertexAdjacencies(void);
  void writeDistances(void);
  ///////////
  void printCumulative(void);
  void printAttr(void);
  void printChars(void);
  void analyzeBatch(Space &s);
  void printIntegrity(void);
  void vertexAdjacentFaces(void);
  void areaAspectRatio(void);
  void processEdgeLengths(void);
  void computeEdgeAngles(void);
  //	void analyzeCumulative(Container*);
  void analyzeCumulative(void);
  ///////////
  void boundWorld(Space&,Controls&);
  void getNiceSet(Space&,Monitor&);
  void findNice(Space&);
  bool checkNiceness(Space&,Vertex*);
  void collectCrossed(Space&,Vertex*,vec_o&);
  bool updateNiceness(Vertex*,vec_o&);	
  void findClosestAxis(Space&,Vertex*,double[2][3]);
  int findExtraPoint(Space&,Vertex*,double[3],int);
  void findCrossed1(Space&,Vertex*,double[2][3],vec_o&);
  void findCrossed2(Space&,double[2][3],vec_o&);
  void getExtraRay(Vertex*,double[2][3],int);
  void collectNiceFaces(Space&,double[2][3],std::vector<Face*>&);
  void getBoxIndexRange(Space&,double[2][3],int[6]);
  void findIntersectedFaces(double[2][3],std::vector<Face*>&,std::vector<Face*>&,std::vector<int>&);
  void findOddMeshes(std::vector<Face*>&,std::vector<int>&,int&,vec_o&);
  ///////////
  void getSeparationDistances(Space&);
  bool findClosest(Space&,Vertex*);
  bool computeClosest(Face*,Vertex*,double&,double[3]);
  bool getPlaneIntersection(Face*,Vertex*,double*,double,double,Point&);
  void getEdgeIntersection(Vertex*,double*[3],Point&);
  void getBoxes(std::vector<Box*>&,Vertex*,int,Space&);
  void getCandidateFaces(std::vector<Box*>&,Vertex*,hset_f&);
  bool faceInNeighborhood(Face*,Vertex*);
  ///////////
  bool checkFaceFaceIntersections(Face*,Face*);
  int checkEdgeEdgeIntersection(Face*,Face*,bool);
  int checkFaceEdgeIntersection(Face*,Face*);
  bool facesParallel(Face*,Face*);
  bool facesColinear(Face*,Face*);
  int numUniqueVertices(Face*,Face*,int[2]);
  void assignFacesToBoxes(Space&);
  ///////////
  int num_orph;
  int countOrphan(void);
  int num_mis;
  int countMissing(void);
  int num_deg;
  int countDegen(void);
  int num_bor;
  int countBorder(void);
  int num_flip;
  int countFlipped(void);
  int num_nonman_e;
  int countNonmanE(void);
  int num_nonman_v;
  int countNonmanV(void);
  int num_obj;
  int countObject(void);
  int num_vert;
  int countVertex(void);
  int num_face;
  int countFace(void);
  int num_edge;
  int countEdge(void);
  double num_vol;
  double countVol(void);
  double countArea(void);
  int num_sep;	// yep, sep for separate components, not a mistake
  int countComponents(void);
  int num_bou;
  int countBoundaries(void);
  int num_indistin;
  int countIndistin(void);
  int countIntFace(void);
  int num_dupl_v;
  int countDuplV(void);
  int num_dupl_f;
  int countDuplF(void);
  ///////////
  int num_clo_cc,num_clo_nn;
  std::pair<int,int> countClosed(void);
  //	int num_man_cc,num_man_nn;
  //	std::pair<int,int> countManifold(void);
  int num_man[3];
  void countManifold(int[3]);
  int num_cons[3];
  void countConsistent(int[3]);
  int num_out[3];
  void countOutward(int[3]);
};

#endif
