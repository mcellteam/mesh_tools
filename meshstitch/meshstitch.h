#ifndef MESHSTITCH_H
#define MESHSTITCH_H 1

#include <list>

class Group;

int distinguishable(double a,double b,double eps);
int maxVert (std::list<Vertex> & v);
void matchZ (std::list<Vertex> & v,std::list<Vertex> & th,int z_value);
void candidateBadVertices (std::list<Vertex> & v, std::list<Vertex> & th,int val);
void printVertices (std::list<Vertex> & v,char *str);
void removeBadVertices (std::list<Vertex> & V,std::list<Vertex> & B);
void findSharedVertices (std::list<Vertex> & v1_match_h,
                         std::list<Vertex> & v2_match_h,
                         std::list<Vertex> & v1_shared_h,
                         std::list<Vertex> & v2_shared_h,
                         double epsilon);
void identifyBadVerticesAndFaces (std::list<Face> & cbf,
                                  std::list<Vertex> & cbv,
                                  std::list<Vertex> & v_shared,
                                  Group &g,
                                  std::list<Face> & bfh,
                                  std::list<Vertex> &bvh);
void candidateBadFaces (std::list<Face> & F,std::list<Face> & th,std::list<Vertex> & V);
void findBadFaces (std::list<Face> & CBF,std::list<Vertex> & V,std::list<Face> & BF);
int maxFace (std::list<Face> & L);
void getFaceSearchPool (std::list<Face> & F,std::list<Face> & qh,int z_value,Vertex** vert_array);
void printFaces (std::list<Face> & L,char *str);
void removeBadFaces (std::list<Face> & F,std::list<Face> & B);
int compare (const void* a, const void* b);

typedef std::vector<ExtraVertex> vec_ev;
typedef std::vector<ExtraVertex>::iterator iter_vec_ev;
typedef std::list<Vertex> list_v;
typedef std::list<Vertex>::iterator iter_list_v;
typedef std::list<Face> list_f;
typedef std::list<Face>::iterator iter_list_f;

#endif
