// Author: Justin Kinney
// Date: Sep 2008

#ifndef MESHMORPH_H
#define MESHMORPH_H 1

#include <cmath>
#include <iostream>

//#include <ext/hash_map>
#include <unordered_map>
//#include <ext/hash_set>
#include <unordered_set>
#include <map>
#include <set>
#include <string>
#include <vector>

class Log;
class Face;
class Edge;
class Nice;
class State;
class Object;
class Vertex;
class Container;
class Controls;
class Refractory;
class Virtual_Disp;
class Gain_Schedule;
class Vertex_Schedule;
class Intersecting_Faces;

typedef  unsigned long int  u4;   /* unsigned 4-byte type */
typedef  unsigned     char  u1;   /* unsigned 1-byte type */

struct f_hash
{
  u4 operator () (Face* i) const { return (u4) i; }
};

struct v_hash
{
  u4 operator () (Vertex* i) const { return (u4) i; }
};

struct e_hash
{
  u4 operator () (Edge* i) const { return (u4) i; }
};

typedef std::equal_to<const Face *>     eqf;
typedef std::equal_to<const Edge *>     eqe;
typedef std::equal_to<const Vertex *>   eqv;

typedef std::less<const Edge *>         lte;
typedef std::less<const Object *>       lto;
typedef std::less<const Vertex *>       ltv;
typedef std::less<const Face *>         ltf;
typedef std::less<const std::string>    lts;
typedef std::less<const double>         ltd;

typedef std::vector<Vertex>		     vec_v;
typedef std::vector<Vertex>::iterator	     v_it;
typedef std::vector<Vertex>::const_iterator  v_cit;
typedef std::vector<Face>		     vec_f;
typedef std::vector<Face>::iterator	     f_it;
typedef std::vector<Face>::const_iterator    f_cit;
typedef std::vector<Edge>		     vec_e;
typedef std::vector<Edge>::iterator	     e_it;
typedef std::vector<Edge>::const_iterator    e_cit;
typedef std::vector<Object>		     vec_o;
typedef std::vector<Object>::iterator	     o_it;
typedef std::vector<Object>::const_iterator  o_cit;

typedef std::vector<Vertex*>		     vec_vp;
typedef std::vector<const Vertex*>  	     vec_cvp;
typedef std::vector<Vertex*>::iterator	     vp_it;
typedef std::vector<Vertex*>::const_iterator vp_cit;
typedef std::vector<const Vertex*>::const_iterator cvp_cit;
typedef std::vector<Face*>		     vec_fp;
typedef std::vector<Face*>::iterator	     fp_it;
typedef std::vector<Face*>::const_iterator   fp_cit;
typedef std::vector<Edge*>		     vec_ep;
typedef std::vector<Edge*>::iterator	     ep_it;
typedef std::vector<Edge*>::const_iterator   ep_cit;
typedef std::vector<Object*>		     vec_op;
typedef std::vector<Object*>::iterator	     op_it;

typedef std::vector<int>		     vec_i;
typedef std::vector<int>::iterator	     i_it;
typedef std::vector<double>		     vec_d;
typedef std::vector<double>::iterator	     d_it;
typedef std::vector<double>::const_iterator  d_cit;
typedef std::vector<std::string>             vec_s;
typedef std::vector<std::string>::const_iterator s_cit;

typedef std::set<std::string,lts>            s_set;
typedef std::set<std::string,lts>::iterator  ss_it;
typedef std::set<Edge*,lte>                  e_set;
typedef std::set<Edge*,lte>::iterator        es_it;
typedef std::set<Vertex*,ltv>                v_set;
typedef std::set<Vertex*,ltv>::iterator      vs_it;
typedef std::set<Face*,ltf>                  f_set;
typedef std::set<Face*,ltf>::iterator        fs_it;
typedef std::set<Object const*,lto>                o_set;
typedef std::set<Object const*,lto>::iterator      os_it;

//typedef __gnu_cxx::hash_set<Vertex*,v_hash,eqv>                      hashset_v;
typedef std::unordered_set<Vertex*,v_hash,eqv>                      hashset_v;
typedef std::map<std::string,Edge*,lts>                              map_s_ep;
//typedef __gnu_cxx::hash_map<Vertex*,int,v_hash,eqv>                  hmap_v;
//typedef __gnu_cxx::hash_map<Vertex*,int,v_hash,eqv>::const_iterator  vhm_cit;
typedef std::unordered_map<Vertex*,int,v_hash,eqv>                  hmap_v;
typedef std::unordered_map<Vertex*,int,v_hash,eqv>::const_iterator  vhm_cit;

struct vector3
{
  double p[3];
  vector3 (void)
  {
    p[0] = 0.0;
    p[1] = 0.0;
    p[2] = 0.0;
  }
  vector3 ( double x, double y, double z)
  {
    p[0] = x;
    p[1] = y;
    p[2] = z;
  }
  vector3& operator= ( const vector3 & v)
  {
    if (&v != this)
    {
      p[0] = v.p[0];
      p[1] = v.p[1];
      p[2] = v.p[2];
    }
    return *this;
  }
  vector3& operator+= ( double a )
  {
    p[0] += a;
    p[1] += a;
    p[2] += a;
    return *this;
  }
  vector3& operator-= ( double a )
  {
    p[0] -= a;
    p[1] -= a;
    p[2] -= a;
    return *this;
  }
  vector3& operator+= ( const vector3 & v )
  {
    if (&v != this)
    {
      p[0] += v.p[0];
      p[1] += v.p[1];
      p[2] += v.p[2];
    }
    return *this;
  }
  vector3& operator*= ( double a )
  {
    p[0] *= a;
    p[1] *= a;
    p[2] *= a;
    return *this;
  }
  vector3 operator- ( const vector3 & v) const
  {
    return vector3( p[0] - v.p[0],
                    p[1] - v.p[1],
                    p[2] - v.p[2] );
  }
  vector3 operator* ( double a) const
  {
    return vector3( p[0] * a,
                    p[1] * a,
                    p[2] * a );
  }
  vector3 operator+ ( double a) const
  {
    return vector3( p[0] + a,
                    p[1] + a,
                    p[2] + a );
  }
  vector3 operator+ ( const vector3 & v) const
  {
    return vector3( p[0] + v.p[0],
                    p[1] + v.p[1],
                    p[2] + v.p[2] );
  }
  vector3 operator/ ( double a) const
  {
    vector3 kQuot;

    if (a != 0.0)
    {
      double fInvScalar = 1.0/a;
      kQuot.p[0] = fInvScalar*p[0];
      kQuot.p[1] = fInvScalar*p[1];
      kQuot.p[2] = fInvScalar*p[2];
    }
    else
    {
      kQuot.p[0] = 1E300;
      kQuot.p[1] = 1E300;
      kQuot.p[2] = 1E300;
    }

    return kQuot;
  }
  double dot ( const vector3 & v) const
  {
    return p[0]*v.p[0]+p[1]*v.p[1]+p[2]*v.p[2];
  }
  double length (void) const
  {
    return sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]);
  }
  vector3 cross ( const vector3 & v) const
  {
    return vector3(p[1]*v.p[2]-p[2]*v.p[1],
                   p[2]*v.p[0]-p[0]*v.p[2],
                   p[0]*v.p[1]-p[1]*v.p[0]);
  }
  void print ( std::ostream & target) const
  {
    target << "["
          << p[0] << " "
          << p[1] << " "
          << p[2] << "]";
  }
};

struct Minmax
{
  double const * min;
  double const * max;
  Minmax (double const * a,double const * b,double const * c)
        :min(NULL),max(NULL)
  {
    if ( *a < *b )
    {
      // a < b
      if ( *a < *c )
      {
        // a < b
        // a < c
        if (*b < *c)
        {
          // a < b
          // b < c
          // a is smallest
          // c is largest
          min=a;
          max=c;
          return;
        }
        else
        {
          // a < b
          // c <= b
          // a is smallest
          // b is largest or tied with c
          min=a;
          max=b;
          return;
        }
      }
      else
      {
        // a < b
        // c <= a
        // c is smallest or tied with a
        // b is largest
        min=c;
        max=b;
        return;
      }
    }
    else
    {
      // b <= a
      if ( *a < *c )
      {
        // b <= a
        // a < c
        // b is smallest or tied with a
        // c is largest
        min=b;
        max=c;
        return;
      }
      else
      {
        // b <= a
        // c <= a
        if (*b < *c)
        {
          // b < c
          // c <= a
          // b is smallest
          // a is largest
          min=b;
          max=a;
          return;
        }
        else
        {
          // b <= a
          // c <= b
          // c is smallest
          // a is largest
          min=c;
          max=a;
          return;
        }
      }
    }
  }
};

struct result
{
  bool line_flag;
  bool poly_flag;
  bool poly_edge_flag;
  result (void)
        :line_flag(false),poly_flag(false),poly_edge_flag(false)
  {
  }
};

bool distinguishable (double a,double b,double epsilon);
bool distinguishable (double a,double b);
bool checkIntSize   (void);

#endif
