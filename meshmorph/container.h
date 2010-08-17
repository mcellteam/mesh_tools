// Author: Justin Kinney
// Date: Sep 2008

#ifndef CONTAINER_H
#define CONTAINER_H 1

#include "meshmorph.h"

#include "octree_visitor_face.h"

typedef std::multimap<Object*,int,lto>            mmap_oi;
typedef std::multimap<Object*,int,lto>::iterator  oi_it;

struct Complex
{
  Object * o;
  int index;
  Complex (Object *oo,int i)
        :o(oo),index(i)
  {
  }
};

class Container
{
private:
  static Container * only_one;
  Object     *fo; // pointer to first object to detect vector copies
  vec_s    files; // array of input file names
  vec_vp  frozen; // sorted vector of Vertex* to frozen vertices
  vec_d    world; // minX, maxX, minY, maxY, minZ, maxZ of world
  vec_vp   cleft; // sorted vector of Vertex* to cleft vertices
  vec_vp    peri; // sorted vector of Vertex* to perisynaptic vertices
  Container                        (void);
  Container                        (Container const &);
  Container & operator =           (Container const &);
  ~Container                       (void);
public:
  hxa7241_graphics::Octree<Face> * octree;
  vec_o        o; //vector of object pointers
  static Container & instance                  (void);
  int      getVertexCount                      (void) const;
  int      getFaceCount                        (void) const;
  int      getEdgeCount                        (void) const;
  void     boundWorld                          (void);
  void     writeSummary                        (std::ostream &);
  void     scanDir                             (void);
  void     scanFiles                           (void);
  void     assessFile                          (Object * const,char const *);
  void     scanFile                            (Object * const,char const *);
  void     computeVertexNormals                (void);
  void     createEdges                         (void);
  void     findVertAdj                         (void);
  void     readFrozen                          (char const *);
  void     readCleft                           (char const *);
  void     findPerisynaptic                    (void);
  void     sortAdjacentFaces                   (void);
  void     deleteme_checkClosestFace           (int const &,int const &,std::string);
  void     checkClosestFace                    (int const &,std::string);
  void     checkClosestFace2                   (int const &,std::string);
  void     checkClosestFace3                   (int const &,std::string);
  void     checkFaces                          (std::string);
  void     checkFacesInOctree                  (void);
  void     printRegionInOctree                 (void);
  void     writeMeshData                       (int const &)  const;
  int      countObjects                        (fp_cit begin,fp_cit end) const;
  //int      calculateVertexNeighborCount        (Vertex const * const v) const;
  int      calculateVertexNeighborCount (Vertex const * const v,vector3 const * pos,vector3 const & normal) const;
  void     findClosestFaceToEachVertex         (void);
  void     findClosestPtInFaceToLocation       (vector3 const & pt,
                                                Face const * const face,
                                                vector3 & p,
                                                double & squareD) const;
  void     findPtInFaceWhereNormalInt          (vector3 const & pt,
                                                vector3 n,
                                                double const &,
                                                Face const * const face,
                                                vector3 & p,
                                                double & squareD) const;
  bool     faceLiesOppositeToNormal            (vector3 const &,
                                                Face const * const face,
                                                vector3 const & n) const;
  bool     boundingBoxFullyInSearchRegion      (vec_d const & sr,
                                                vec_d const & bb) const;
  bool     vertexOutsideOctreeBounds           (vector3 const * const new_pos);
  bool     closestPtIsInSearchCone             (vector3 const & pt,
                                                vector3 const & p,
                                                vector3 const & n,
                                                double const &) const;
  bool     findClosestPtToBarycenterAmongFaces (vector3 const &,Face * const f,
                                                fp_cit,fp_cit,
                                                double const & cone_radius,
                                                vector3 &,
                                                double &,
                                                Face * &,int &) const;
  bool     findClosestPtToVertexAmongFaces     (Vertex const * const,
                                                fp_cit,
                                                fp_cit,
                                                double const &,
                                                vector3 &,
                                                double &,
                                                Face * &,int &) const;
  bool     findClosestPtToVertex               (Vertex const * const,
                                                vector3 &,double &,
                                                Face * &) const;
  bool     findClosestPtToBarycenter           (vector3 const & pt,Face * const f,
                                                vector3 &,double &,
                                                Face * &) const;
  double   getMinEdgeAngle                     (void) const;
  double   getWorld                            (int const &) const;
  double   calculateCleftArea                  (void) const;
  vec_d    getSphericalConeBoundingBox         (vector3 const & pt,
                                                vector3 n,
                                                double const & cone_radius) const;
  void     getSphericalConeRays                (vector3 const & pt,
                                                vector3 n,
                                                double const & cone_radius,
                                                std::vector<vector3> & rays) const;
  mmap_oi  loadMap                       (char const *,s_set &);
  Object * getObjectPointer                    (std::string) const;
  void     findVertexNeighbors (const int & group);
//  vec_d    getTruncatedConeBoundingBox         (vector3 const & pt,
//                                                vector3 n,
//                                                double const &) const;
  void findOrthogonalRegions                   (void);
  void writeOrthVerts2Recon3D                  (std::vector<const Vertex*> & orthogonal_vertices);
  void writeOrthVerts2Dreamm                   (std::vector<const Vertex*> & orthogonal_vertices);
  void openContourFile (std::ofstream & contour_handle,double z);
  void addVertexToFile (std::ofstream & contour_handle,const Vertex * vp);
  //void initializeContour (std::ofstream & contour_handle,const std::string & name);
  void initializeContour (std::ofstream & contour_handle,const std::string & name,const double &,const double &,const double &);
  void finalizeContour (std::ofstream & contour_handle);
  void writeContoursSingleSection (std::vector<const Vertex*> & single_section, double z);
  void binVertexNormalAngle (double * bins,int * count,const double & cos_angle);

  std::vector<Complex> loadVector              (const char *filename,s_set & not_found);

  /** Check if vertex is located in a synaptic cleft
    * as specified by user input data file.
   * \param[in] v Vertex of interest.
   * \return True if vertex is synaptic cleft; otherwise false;
   */

  bool vertexIsCleft (Vertex * const v) const
  {
    return binary_search(cleft.begin(),cleft.end(),v);
  }

  /** Check if vertex is located in perisynaptic regions
    * (surrounding cleft regions as specified by user input data file).
   * \param[in] v Vertex of interest.
   * \return True if vertex is synaptic cleft; otherwise false;
   */

  bool vertexIsPeri (Vertex * const v) const
  {
    return binary_search(peri.begin(),peri.end(),v);
  }

  /** Check if vertex is frozen.
   * \param[in] v Vertex of interest.
   * \return True if vertex is frozen; otherwise false;
   */

  bool vertexIsFrozen (Vertex * const v) const
  {
    return binary_search(frozen.begin(),frozen.end(),v);
  }

  /** Get number of input mesh files.
   * \return Number of input files.
   */

  int getFileCount (void) const
  {
    return files.size();
  }
};

#endif
