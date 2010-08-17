// Author: Justin Kinney
// Date: Feb 2009

#ifndef OBJECT_H
#define OBJECT_H 1

#include "face.h"
#include "stats.h"

class Object
{
private:
  std::string name;		// object name
  std::vector<Vertex*> v;	// container of pointers to all vertices in object
  std::vector<Face*> f;		// container of pointers to all faces in object
  std::vector<Edge*> e;		// container of pointers to all edges in object
  bool outward;			// true=outwardly-oriented face normals
  hmap_fdp iv;			// store intersection force
  hmap_ff intf;			// store intersecting faces
  hmap_vi nice;			// vertex* is key to int nice value
  hmap_e flip;			// collection of flip values for edges
  int num_sep;			// number of components (separate meshes in object)
  int num_bou;			// number of separate boundaries in object, same as # holes
  double vol;			// object volume
  int genus;			// object genus
  mmap_iv vp;
  mmap_ib found;
  bool contig_f;		// true if vface indexing is contiguous
  int vec_cf[5];		// bad face index +/- 2
  bool contig_v;		// true if vertex indexing is contiguous
  int vec_cv[5];		// bad vertex index +/- 2
  double bb[6];		// bounding box [xmin,ymin,zmin,xmax,ymax,zmax]
  Stats area;
  Stats aspect_ratio;
  Stats edge_length;
  Stats edge_angle;
  Stats adjacent_face;
  vec_e border;		        // border edges (edge with a single adjacent face)
  vec_e nonman_e;	        // nonmanifold edges (edge with three or more adjacent faces)
  vec_e flipped;		// flipped edges (edges traversed more than once)
  vec_v indistin_v;             // indistinguishable vertices 
  vec_e indistin;	        // indistinguishable edges (edges with vertices 
                                // that are indistinguishbale using double precision)
  vec_v dupl_v_index;           // vertices with duplicate vertex indices
  vec_f dupl_f_index;           // faces with duplicate face indices
  vec_v nonman_v;	        // nonmanifold vertices 
  vec_v orphan;	                // vertices not referenced by any face
  vec_f missing_f;	        // faces with one or more nonexistant referenced vertices
  vec_i missing_v;	        // nonexistant referenced vertices
  vec_f degen;	                // faces with two or more identical vertex references
  hmap_fd bad_aspect;	        // faces with aspect ratio larger than threshold
  hmap_ed bad_angle;	        // edges with angle > max threshold or < min threshold
  hmap_ed bad_length;	        // edges with length > max threshold or < min threshold
  vec_f zero_area;
public:
  Object(std::string);
  ~Object(void);
  int getNumDegenerateFaces (void)
  {
    return degen.size();
  }
  void addMissingVertex (int i)
  {
    missing_v.push_back(i);
  }
  void addMissingFace (Face * ff)
  {
    missing_f.push_back(ff);
  }
  int getNumMissingVertices (void)
  {
    return missing_v.size();
  }
  int getNumOrphanVertices (void)
  {
    return orphan.size();
  }
  void addNonmanVertices (Vertex * vv)
  {
    nonman_v.push_back(vv);
  }
  int getNumNonmanVertices (void)
  {
    return nonman_v.size();
  }
  int getNumDuplFaceIndices (void)
  {
    return dupl_f_index.size();
  }
  int getNumDuplVertexIndices (void)
  {
    return dupl_v_index.size();
  }
  int getNumIndistinVertices (void)
  {
    return indistin_v.size();
  }
  int getNumFlippedEdges (void)
  {
    return flipped.size();
  }
  int getNumNonmanEdges (void)
  {
    return nonman_e.size();
  }
  int getNumBorderEdges (void)
  {
    return border.size();
  }
  double getAreaSum (void)
  {
    return area.getSum();
  }
  std::pair<ib_iterator,ib_iterator> findIndexInMap2 (int i)
  {
    return found.equal_range(i);
  }
  std::pair<iv_iterator,iv_iterator>  findIndexInMap (int i)
  {
    return vp.equal_range(i);
  }
  void addIndexBoolPair (int i,bool b)
  {
    found.insert(std::make_pair(i,b));
  }
  void addVertexIndexPair (int i,Vertex * vv)
  {
    vp.insert(std::make_pair(i,vv));
  }
  double getVolume (void)
  {
    return vol;
  }
  int getNumBoundaries (void)
  {
    return num_bou;
  }
  int getNumSep (void)
  {
    return num_sep;
  }
  void removeIntFaceLHS (Face * lhs)
  {
    // delete vector<face*>*
    delete intf[lhs];
    // remove element from table
    intf.erase(lhs);
  }
  void addIntFaceRHS (Face * lhs,Face * rhs)
  {
     (*intf[lhs]).push_back(rhs);
  }
  void addIntFaceLHS (Face * ff,vec_f * i)
  {
     intf[ff]=i;
  }
  int getNumIntFace (void)
  {
    return intf.size();
  }
  ff_iterator getFirstIntFace (void)
  {
    return intf.begin();
  }
  ff_iterator getOnePastLastIntFace (void)
  {
    return intf.end();
  }
  bool getOutward (void)
  {
    return outward;
  }
  bool noFaces (void)
  {
    return f.empty();
  }
  bool noVertices (void)
  {
    return v.empty();
  }
  void addVertex (Vertex * vv)
  {
    v.push_back(vv);
  }
  void addFace (Face * ff)
  {
    f.push_back(ff);
  }
  int getNumFaces (void)
  {
    return f.size();
  }
  int getNumVertices (void)
  {
    return v.size();
  }
  int getNumEdges (void)
  {
    return e.size();
  }
  Vertex * getFrontVertex (void)
  {
    return v.front();
  }
  v_iterator getFirstVertex (void)
  {
    return v.begin();
  }
  v_iterator getOnePastLastVertex (void)
  {
    return v.end();
  }
  f_iterator getFirstFace (void)
  {
    return f.begin();
  }
  f_iterator getOnePastLastFace (void)
  {
    return f.end();
  }
  e_iterator getFirstEdge (void)
  {
    return e.begin();
  }
  e_iterator getOnePastLastEdge (void)
  {
    return e.end();
  }
  // *********************
  void orphanMissingContig(void);
  void degenContig(void);
  void vertexAdjacentFaces(void);
  void areaAspectRatio(Controls&);
  void findIntersectingFaces(Container*,Space&);
  void processEdgeLengths(Controls&);
  void computeEdgeAngles(void);
  void computeVolume(void);
  void printIntegrity(Controls&);
  // *********************
  void analyze(Container*,Controls&,Space&);
  void print(Controls&);
  void printChars(Controls&);
  void checkIntegrity(void);
  void setAll(Vertex*,hmap_fi&,int&);
  void setZero(Edge*,hmap_fi&,int);
  void getGroups(Vertex*,hmap_fi&,set_i&);
  int getLowest(set_i&);
  void replaceGroups(Vertex*,hmap_fi&,int);
  void printAttr(Controls&);
  void vertexDistin(void);
  // *********************
  void addOriginal(int,double*);
  void createEdges(void);
  int setNumDigits(void);
  void buildEdge(Face*,Vertex*,Vertex*,map_se&,int);
  Edge* findEdge(Vertex*,Vertex*,map_se&,int);
  void createEdge(Face*,Vertex*,Vertex*,map_se&,int);
  void addVertexPointers(void);
  int getMaxVertex(void);
  void fixFaces(Vertex**);
  void fixEdges(Vertex**);
  void findVertexAdjacencies(void);
  void findNeighborhoods(void);
  void boundObject(double*);
  double getMeanEdgeLength(void);
  //////////////
  bool intersectingFacesExist(Face*);
  bool faceInTable_intf(Face*f);
  ff_iterator findFaceInTable_intf(Face*f);
  //////////////
  bool vertexIsNice(Vertex*);
  int getVertexNiceness(Vertex*);
  void setVertexNiceness(Vertex*,int);
  //////////////
  void addEdgeFlip(Edge*,int);
  void clearFlipTable(void);
  int getEdgeFlip(Edge*);
  //////////////
  bool ifFrozen(hmap_vd&,Vertex*);
  void newFindNeighborhoods(void);
  bool thawedAndAble(hmap_vd&,set_v&);
  void collectFaces(hmap_vd&,set_v&,vec_f&);
  bool processEdge(Edge*,hmap_vd&,vec_e&,Vertex*);	
  //////////////
  void evalCharacteristics(Container*,Controls&,Space&);
  void evalAttributes(Space&);
  bool isClosed(void);
  bool isManifold(void);
  bool isConsistent(void);
  bool isOutward(Space&);
  bool rayIntersectsBB(double[2][3],Face*,double[3]);
  bool rayIntersectsSide(const char*,double[2][3],double[6],double[3]);
  void countBoundaries(void);
  int countComponents(void);
  void computeGenus(void);
  bool edgesManifold(bool);
  bool verticesManifold(bool);
  bool removeSelectedFaces(vec_f&,std::vector<Face*>&);
  void getSelectedFaces(vec_f&,std::vector<Face*>);
  void verifyEdges(void);
  bool goodIntegrity(void);
  //
  void store(Container &);
  //
  std::string getName (void)
  {
    return name;
  }
};

#endif
