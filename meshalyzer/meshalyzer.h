// Author: Justin Kinney
// Date: Feb 2009

#ifndef MESHALYZER_H
#define MESHALYZER_H 1

#include <string>
#include <vector>
#include <map>
#include <unordered_map>
#include <set>
#include <unordered_set>

class Boundary;
class Box;
class Container;
class Controls;
class Edge;
class Face;
class Face_Pair;
class Manip;
class Monitor;
class Object;
class Neighbor;
class Space;
class Vertex;

typedef  unsigned long int  u4;   /* unsigned 4-byte type */
typedef  unsigned     char  u1;   /* unsigned 1-byte type */

struct eqf
{
  bool operator()(const Face* s1, const Face* s2) const
  {
    return s1==s2;
  }
};

struct eqv
{
  bool operator()(const Vertex* s1, const Vertex* s2) const
  {
    return s1==s2;
  }
};

struct eqi
{
  bool operator()(const int s1, const int s2) const
  {
    return s1==s2;
  }
};

struct eqe
{
  bool operator()(const Edge* s1, const Edge* s2) const
  {
    return s1==s2;
  }
};

struct f_hash
{
  u4 operator()(Face* i) const { return (u4) i; }
};

struct v_hash
{
  u4 operator()(Vertex* i) const { return (u4) i; }
};

struct e_hash
{
  u4 operator()(Edge* i) const { return (u4) i; }
};

struct i_hash
{
  u4 operator()(int i) const { return (u4) i; }
};

struct lti
{
  bool operator()(const int s1, const int s2) const
  {
    return s1 < s2;
  }
};

struct gti
{
  bool operator()(const int s1, const int s2) const
  {
    return s1 > s2;
  }
};

struct gtd
{
  bool operator()(const double s1, const double s2) const
  {
    return s1 > s2;
  }
};

struct ltv
{
  bool operator()(const Vertex* s1, const Vertex* s2) const
  {
    return s1 < s2;
  }
};

struct lts
{
  bool operator()(const std::string s1, const std::string s2) const
  {
    return s1 < s2;
  }
};

struct vSortComp
{
  bool operator()(const std::pair<double,Vertex*> lhs,const std::pair<double,Vertex*> rhs)
  {
    //return (lhs->index < rhs->index);
    return (lhs.second > rhs.second);
  }
};

typedef std::vector<Vertex*>			vec_v;
typedef std::vector<Vertex*>::iterator	v_iterator;
typedef std::vector<Face*>				vec_f;
typedef std::vector<Face*>::iterator	f_iterator;
typedef std::vector<Face*>::const_iterator	c_f_iterator;
typedef std::vector<Edge*>				vec_e;
typedef std::vector<Edge*>::iterator	e_iterator;
typedef std::vector<Object*>			vec_o;
typedef std::vector<Object*>::iterator	o_iterator;
typedef std::vector<Box*>				vec_b;
typedef std::vector<Box*>::iterator		b_iterator;
typedef std::vector<Box*>::const_iterator		c_b_iterator;
typedef std::vector<int>				vec_i;
typedef std::vector<int>::iterator		i_iterator;
typedef std::vector<double>				vec_d;
typedef std::vector<double>::iterator	d_iterator;
typedef std::vector<double>::const_iterator	c_d_iterator;

typedef std::map<std::string,Edge*,lts,std::allocator<Edge*> >					map_se;

typedef std::multimap<int,bool,lti,std::allocator<bool> >						mmap_ib;
typedef std::multimap<int,bool,lti,std::allocator<bool> >::iterator				ib_iterator;
typedef std::multimap<double,Vertex*,gtd,std::allocator<Vertex*> >				mmap_dv;
typedef std::multimap<double,Vertex*,gtd,std::allocator<Vertex*> >::iterator	dv_iterator;
typedef std::multimap<int,Face*,lti,std::allocator<Face*> >						mmap_if;
typedef std::multimap<int,Face*,lti,std::allocator<Face*> >::iterator			if_iterator;
typedef std::multimap<int,Vertex*,gti,std::allocator<Vertex*> >					mmap_iv;
typedef std::multimap<int,Vertex*,gti,std::allocator<Vertex*> >::iterator		iv_iterator;

typedef std::unordered_map<Vertex*,int,v_hash,eqv,std::allocator<int> >								hmap_vi;
typedef std::unordered_map<Face*,double*,f_hash,eqf,std::allocator<double*> >							hmap_fdp;
typedef std::unordered_map<Face*,double,f_hash,eqf,std::allocator<double*> >							hmap_fd;
typedef std::unordered_map<Face*,double,f_hash,eqf,std::allocator<double*> >::iterator					fd_iterator;
typedef std::unordered_map<Face*,int,f_hash,eqf,std::allocator<int> >									hmap_fi;
typedef std::unordered_map<Face*,int,f_hash,eqf,std::allocator<int> >::iterator						fi_iterator;
typedef std::unordered_map<Vertex*,double,v_hash,eqv,std::allocator<double> >							hmap_vd;
typedef std::unordered_map<Vertex*,double,v_hash,eqv,std::allocator<double> >::iterator				vd_iterator;
typedef std::unordered_map<Edge*,double,e_hash,eqe,std::allocator<double> >							hmap_ed;
typedef std::unordered_map<Edge*,double,e_hash,eqe,std::allocator<double> >::iterator					ed_iterator;
typedef std::unordered_map<Face*,vec_f*,f_hash,eqf,std::allocator<std::vector<Face*>* > >				hmap_ff;
typedef std::unordered_map<Face*,vec_f*,f_hash,eqf,std::allocator<std::vector<Face*>* > >::iterator	ff_iterator;
typedef std::unordered_map<Edge*,int,e_hash,eqe,std::allocator<int> >									hmap_e;
typedef std::unordered_map<int,Vertex*,i_hash,eqi,std::allocator<Vertex*> >							hmap_iv;

typedef std::set<int,lti>									set_i;
typedef std::set<Vertex*,ltv>								set_v;

typedef std::unordered_set<Face*,f_hash,eqf>				hset_f;
typedef std::unordered_set<Face*,f_hash,eqf>::iterator		hf_iterator;

int distinguishable(double a,double b,double epsilon);
int distinguishable(double a,double b);
double dot(double a[3],double b[3]);
bool getPointEdgeDistance(double p[3],Face *f);
void biggest(double *x, int& big); 
void threeValueSort(double p1,double p2,double p3, double &biggest, double &smallest);
void checkLineFaceIntersection(Face *f,double lp[2][3],bool &line_flag,
                               bool &poly_flag, bool &poly_edge_flag);
std::string keyPair(int a,int b,int num_digits);
bool edgeMatch(Edge *e,int va,int vb);
bool overlap(double a1,double a2,double b1,double b2);
bool faceBBsOverlap(double bb[6],Face *f);

#endif
