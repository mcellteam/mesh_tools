// Author: Justin Kinney
// Date: June 2010

#ifndef MESHHEAL_H
#define MESHHEAL_H 1

#include <iostream>
#include <cmath>
#include <string>
#include <vector>
#include <list>
#include <map>

class Controls;
class Edge;
class Face;
class Object;
class Vertex;

struct Closest
{
  double distance;
  Vertex * good_vertex;
  Vertex * bad_vertex;
  void clear (void)
  {
    distance = 1e300;
    good_vertex=NULL;
    bad_vertex=NULL;
  }
  void add (const double & d,Vertex * vA,Vertex * vB)
  {
    if (d<distance)
    {
      distance = d;
      good_vertex = vA;
      bad_vertex = vB;
    }
  }
};

struct ltd
{
  bool operator()(const double s1, const double s2) const
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

struct lts
{
  bool operator()(const std::string s1, const std::string s2) const
  {
    return s1 < s2;
  }
};

typedef std::list<Vertex*>                                                 list_v;
typedef std::list<Vertex*>::iterator                                  vl_iterator;
typedef std::vector<Vertex*>                                                vec_v;
typedef std::vector<Vertex*>::iterator                                 v_iterator;
typedef std::vector<Vertex*>::const_iterator                          cv_iterator;
typedef std::vector<Face*>                                                  vec_f;
typedef std::vector<Face*>::iterator                                   f_iterator;
typedef std::vector<Edge*>                                                  vec_e;
typedef std::vector<Edge*>::iterator                                   e_iterator;
typedef std::vector<Edge*>::const_iterator                            ce_iterator;
typedef std::map<std::string,Edge*,lts,std::allocator<Edge*> >             map_se;
typedef std::map<int,Vertex*,gti,std::allocator<Vertex*> >                mmap_iv;
typedef std::map<int,Vertex*,gti,std::allocator<Vertex*> >::iterator  iv_iterator;
typedef std::multimap<double,Vertex*,ltd,std::allocator<Vertex*> >                         map_dv;
typedef std::multimap<double,Vertex*,ltd,std::allocator<Vertex*> >::iterator          dv_iterator;
typedef std::multimap<double,Vertex*,ltd,std::allocator<Vertex*> >::const_iterator  c_dv_iterator;

std::string format(char const *fmt, ...);

struct Box
{
  double p[6];
  std::string name;
  Box (void)
        :name()
  {
    p[0] = 0.0;
    p[1] = 0.0;
    p[2] = 0.0;
    p[3] = 0.0;
    p[4] = 0.0;
    p[5] = 0.0;
  }
  void init (std::string str, 
            double xMin,
            double xMax,
            double yMin,
            double yMax,
            double zMin,
            double zMax
             )
  {
    name = str;
    p[0] = xMin;
    p[1] = xMax;
    p[2] = yMin;
    p[3] = yMax;
    p[4] = zMin;
    p[5] = zMax;
  }
  double xmin (void) const
  {
    return p[0];
  }
  double xmax (void) const
  {
    return p[1];
  }
  double ymin (void) const
  {
    return p[2];
  }
  double ymax (void) const
  {
    return p[3];
  }
  double zmin (void) const
  {
    return p[4];
  }
  double zmax (void) const
  {
    return p[5];
  }
};

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
  void init ( double x, double y, double z)
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

#endif
