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
  bool good_integrity;	// integrity of mesh
  double setbb[6];	// bounding box [xmin,ymin,zmin,xmax,ymax,zmax]
  int num_orph;
  int num_mis;
  int num_deg;
  int num_bor;
  int num_flip;
  int num_nonman_e;
  int num_nonman_v;
  int num_obj;
  int num_vert;
  int num_face;
  int num_edge;
  double num_vol;
  int num_sep;	// yep, sep for separate components, not a mistake
  int num_bou;
  int num_indistin;
  int num_dupl_v;
  int num_dupl_f;
  int num_clo_cc,num_clo_nn;
  int num_man[3];
  int num_cons[3];
  int num_out[3];
public:
  ~Container(void);
  Container(void);
  ///////////
  void setIntegrity (bool b)
  {
    good_integrity = b;
  }
  void addAreaRange (c_d_iterator i, c_d_iterator j)
  {
    area.insertElementRange(i,j);
  }
  void addAspectRatioRange (c_d_iterator i, c_d_iterator j)
  {
    aspect_ratio.insertElementRange(i,j);
  }
  void addEdgeLengthRange (c_d_iterator i, c_d_iterator j)
  {
    edge_length.insertElementRange(i,j);
  }
  void addEdgeAngleRange (c_d_iterator i, c_d_iterator j)
  {
    edge_angle.insertElementRange(i,j);
  }
  void addAdjacentFaceRange (c_d_iterator i, c_d_iterator j)
  {
    adjacent_face.insertElementRange(i,j);
  }
  void clearObjects (void)
  {
    o.clear();
  }
  void addFile (std::string s)
  {
    files.push_back(s);
  }
  std::string getFile (int i)
  {
    return files[i];
  }
  void incrementNumFiles (void)
  {
    num_files++;
  }
  int getNumFiles (void)
  {
    return num_files;
  }
  ///////////
  int getScore(void);
  int findExtraPoint(Space&,Vertex*,double[3],int);
  int checkEdgeEdgeIntersection(Face*,Face*,bool);
  int checkFaceEdgeIntersection(Face*,Face*);
  int numUniqueVertices(Face*,Face*,int[2]);
  int countOrphan(void);
  int countMissing(void);
  int countDegen(void);
  int countBorder(void);
  int countFlipped(void);
  int countNonmanE(void);
  int countNonmanV(void);
  int countObject(void);
  int countVertex(void);
  int countFace(void);
  int countEdge(void);
  int countComponents(void);
  int countBoundaries(void);
  int countIndistin(void);
  int countIntFace(void);
  int countDuplV(void);
  int countDuplF(void);
  bool checkNiceness(Space&,Vertex*);
  bool updateNiceness(Vertex*,vec_o&);	
  bool findClosest(Space&,Vertex*);
  bool computeClosest(Face*,Vertex*,double&,double[3]);
  bool getPlaneIntersection(Face*,Vertex*,double*,double,double,Point&);
  bool faceInNeighborhood(Face*,Vertex*);
  bool checkFaceFaceIntersections(Face*,Face*);
  bool facesParallel(Face*,Face*);
  bool facesColinear(Face*,Face*);
  void update(Object*);
  void clear(void);
  void printBatch(Controls&);
  void scanDir(const char*);
  void scanFile(Object*,std::string);
  void createEdges(void);
  void addVertexPointers(void);
  void findVertexAdjacencies(void);
  void writeDistances(void);
  void printCumulative(void);
  void printAttr(void);
  void printChars(void);
  void analyzeBatch(Space &s);
  void printIntegrity(void);
  void vertexAdjacentFaces(void);
  void areaAspectRatio(void);
  void processEdgeLengths(void);
  void computeEdgeAngles(void);
  void analyzeCumulative(void);
  void boundWorld(Space&);
  void getNiceSet(Space&,Monitor&);
  void findNice(Space&);
  void collectCrossed(Space&,Vertex*,vec_o&);
  void findClosestAxis(Space&,Vertex*,double[2][3]);
  void findCrossed1(Space&,Vertex*,double[2][3],vec_o&);
  void findCrossed2(Space&,double[2][3],vec_o&);
  void getExtraRay(Vertex*,double[2][3],int);
  void collectNiceFaces(Space&,double[2][3],std::vector<Face*>&);
  void getBoxIndexRange(Space&,double[2][3],int[6]);
  void findIntersectedFaces(double[2][3],std::vector<Face*>&,std::vector<Face*>&,std::vector<int>&);
  void findOddMeshes(std::vector<Face*>&,std::vector<int>&,int&,vec_o&);
  void getSeparationDistances(Space&);
  void getEdgeIntersection(Vertex*,double*[3],Point&);
  void getBoxes(std::vector<Box*>&,Vertex*,int,Space&);
  void getCandidateFaces(std::vector<Box*>&,Vertex*,hset_f&);
  void assignFacesToBoxes(Space&);
  void countManifold(int[3]);
  void countConsistent(int[3]);
  void countOutward(int[3]);
  double countVol(void);
  double countArea(void);
  Object* processFile(std::string);
  std::pair<int,int> countClosed(void);
};

#endif
