// Author: Justin Kinney
// Date: Sep 2008

#ifndef LOG_H
#define LOG_H 1

#include <fstream>

#include <cstdlib>
#include <cstdarg>
#include <cstdio>
#include "meshmorph.h"

#include "state.h"

typedef std::map<std::string,double,lts>                             map_s_d;
typedef std::map<std::string,double,lts>::iterator                   sd_it;

class Log
{
private:
  Log                         (void);
  Log                         (Log const &);
  Log& operator =             (Log const &);
  ~Log                        (void);
  static Log * only_one;
  std::ofstream  Cfile;
  std::ofstream  Mfile;
  // closest point search
  int        num_verts; // # vertices for which closest point has been searched
  double   leaves_mean; // mean number of octree leaves visited in search for closest point
  int       leaves_min; // min  number of octree leaves visited in search for closest point
  int       leaves_max; // max  number of octree leaves visited in search for closest point
  double     face_mean; // mean number of faces returned from octree in search for closest point
  int         face_min; // min  number of faces returned from octree in search for closest point
  int         face_max; // max  number of faces returned from octree in search for closest point
  double  f_check_mean; // mean number of faces checked in search for closest point
  int      f_check_min; // min  number of faces checked in search for closest point
  int      f_check_max; // max  number of faces checked in search for closest point
  // vertex scheduling
  int         num_sets; // number of vertex sets used in statistics
  double      vps_mean; // mean number of vertices per set
  int          vps_min; // min  number of vertices per set
  int          vps_max; // max  number of vertices per set
  // (note not all verts necessarily moved)
  double    moved_mean; // mean number of moved vertices per set
  int        moved_min; // min  number of moved vertices per set
  int        moved_max; // max  number of moved vertices per set
  // moved vertices
  int                N; // count of moved vertices per group
  double         d_min; // min actual displacement
  double         d_max; // max actual displacement
  vec_d             md; // mean vertex displacement for {N-1,N} vertices of group
  double       fs_mean; // mean number of vertices undergoing full closest point search
  int           fs_min; // min number of vertices undergoing full closest point search
  int           fs_max; // max number of vertices undergoing full closest point search
  double       ps_mean; // mean number of vertices undergoing partial closest point search
  double        ps_min; // min number of vertices undergoing partial closest point search
  double        ps_max; // max number of vertices undergoing partial closest point search

  double        fs_face_change_mean; // number of vertices whose closest faces changed
  double        fs_face_change_min;  // number of vertices whose closest faces changed
  double        fs_face_change_max;  // number of vertices whose closest faces changed
  double        ps_face_change_mean; // number of vertices whose closest faces changed
  double        ps_face_change_min;  // number of vertices whose closest faces changed
  double        ps_face_change_max;  // number of vertices whose closest faces changed
//  double        fs_mean_face_change_is_adj; // number of new closest faces that are adjacent to moved vertex
//  double        ps_mean_face_change_is_adj; // number of new closest faces that are adjacent to moved vertex
  // topN_val[0]: stored virtual displacement as estimated
  // before the move for the last moved vertex.
  // topN_val[1]: stored virtual displacement as estimated
  // after the move for the last moved vertex.
  vec_d       topN_val;
  // vd_val[0]: stored virtual displacement as estimated
  // before the move for the last moved vertex.
  // vd_val[1]: stored virtual displacement as estimated
  // after the move for the last moved vertex.
  vec_d         vd_val;
  // sepdis[0]: stored extracellular width as measured
  // before the move for the last moved vertex.
  // sepdis[1]: stored extracellular width as measured
  // after the move for the last moved vertex.
  vec_d         sepdis;
  //vec_d         p_orig; 
  vector3         p_orig; 
  vector3 const *  p_new; 
  vector3        cp_orig; 
  vector3         cp_new; 
  vec_i        trouble;
  // space partitioning
  int        num_faces; // number of faces used in statistics
  double      bpf_mean; // mean number of boxes per face
  int          bpf_min; // min  number of boxes per face
  int          bpf_max; // max  number of boxes per face
  int        num_boxes; // number of boxes used in statistics
  int  num_empty_boxes; // number of boxes with no faces
  double      fpb_mean; // mean number of faces per box
  int          fpb_min; // min  number of faces per box
  int          fpb_max; // max  number of faces per box
  // time
  double            t1; // track elapsed time
  timeval          tim; // time the main loop
  time_t      currtime; // store time
public:
  static Log & instance               (void);

  int getN (void)
  {
    return N;
  }

  void    groupInit                   (void);
  void    writeVertMoveDistribution   (int const &);
  void    writeFiles                  (int const & group);
  void    writeSepDistances           (int const &);
  void    updateHauss_1               (std::string const &,double const &,map_s_d &);
  void    updateHauss_2               (std::string const &,o_it const &,double const &,map_s_d &);
  void    updateHauss_1_noself        (std::string const &,double const &,map_s_d &);
  void    updateHauss_2_noself        (std::string const &,o_it const &,double const &,map_s_d &);
  void    writeCommandSettings        (void);
  void    closeFiles                  (void);
  void    updateFile                  (int const &,bool,double const &);
  void    statusFileInit              (void);
  void    writeObjectList             (void);
  void    writeObjectData             (void);
  //void    updateVertDisplStats        (double const &);
  void    updateVertDisplStats        (Vertex * const,vector3 const &);
  void    clearVals                   (void);
  void    writeMoveSummary            (int const &,int const &);
  void    writeRefracted              (int const &);
  void    writeRefractedNow           (int const & group,int code);
  void    writeIntersected            (int const &);
  void    writeNonnice                (int const &);
  void    writeVD                     (int const &);
  void    writeAD                     (int const & group);
  void    printNumNonnice             (std::ostream &) const;
  void    printNumInt                 (std::ostream &) const;
  void    setTime                     (time_t);
  void    openMainFile                (void);
  void    recordTime                  (std::string const &);
  void    updateBoxesPerFaceStats     (int const &);
  void    calculateFacesPerBoxStats   (void);
  void    printPartitioningStats      (std::ostream &) const;
  void    updateClosestPtStats        (int const &,int const &,int const &);
  void    printClosestPtStats         (std::ostream &) const;
  void    updateVertexSchedulingStats (int const &);
  void    updateMovedVertsFromSet     (void);
  void    printVertexSchedulingStats  (std::ostream &) const;
  void    updateClosestPtSearchStats  (int const &,int const &,Search_Stats const &);
  void    printClosestPtSearchStats   (std::ostream &) const;
  void    openOrDie                   (std::ofstream * const handle,std::string str,const int & group);
  // originally vtrack
  bool    isGood                      (Vertex const * const) const;
  void    writeDetailedMoveInfo       (Vertex *,int,int);
  void    printBad                    (Vertex const * const) const;
  void    setDetailedInfoPreMove      (Vertex * const); 
  void    setDetailedInfoPostMove     (Vertex *); 
  //
//  std::string formatv(char const *fmt, va_list args)
//  {
//    char buffer[512];
//    uint needLen;
//    //va_list argsCopy;
//    //va_copy(argsCopy, args);
//    va_list argsCopy = args;
//  
//    needLen = vsnprintf(buffer, sizeof(buffer), fmt, args);
//    if (needLen <= sizeof(buffer))
//    {
//      va_end(argsCopy);
//      return std::string(buffer);
//    }
//    else
//    {
//      char *big_buffer = new char[needLen + 1];
//      vsnprintf(big_buffer, needLen + 1, fmt, argsCopy);
//      std::string ret(big_buffer);
//      delete[] big_buffer;
//      va_end(argsCopy);
//      return ret;
//    }
//  };
//  
//  std::string format(char const *fmt, ...)
//  {
//    va_list args;
//    va_start(args, fmt);
//    std::string ret(formatv(fmt, args));
//    va_end(args);
//    return ret;
//  };

  std::string format(char const *fmt, ...)
  {
    va_list args;
    va_start(args, fmt);
    std::string ret;
    //std::string formatv(char const *fmt, va_list args)
    //{
    //uint needLen;
    //va_list argsCopy;
    //va_copy(argsCopy, args);
    //va_list argsCopy = args;
  
    char buffer[512];
    uint needLen = vsnprintf(buffer, sizeof(buffer), fmt, args);
    if (needLen <= sizeof(buffer))
    {
      //va_end(argsCopy);
      //return std::string(buffer);
      ret = buffer;
    }
    else
    {
      char *big_buffer = new char[needLen + 1];
      vsnprintf(big_buffer, needLen + 1, fmt, args);
      //std::string ret(big_buffer);
      ret = big_buffer;
      delete[] big_buffer;
      //va_end(argsCopy);
      //return ret;
    }
    //};
  
    va_end(args);
    return ret;
  };

};


//struct Grid
//{
//  vector3 unit_u;
//  vector3 unit_v;
//  double uv_vert1_u;
//  double uv_vert2[2];
//  Grid(Face*);
//  void computeBarycenter(vector3 &,int,int,Face*);
//};

#endif
