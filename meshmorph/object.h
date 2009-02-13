// Author: Justin Kinney
// Date: Sep 2008

#ifndef OBJECT_H
#define OBJECT_H 1

#include "meshmorph.h"

typedef std::map<std::string,Edge*,lts>::const_iterator              msep_cit;

class Object
{
private:
  std::string  name; // object name, derived from file name
  Vertex       * fv; // pointer to first Vertex to detect vector copies
  Face         * ff; // pointer to first Face to detect vector copies
  Edge         * fe; // pointer to first Edge to detect vector copies
public:
  vec_v     v; // container of all vertices in object
  vec_f     f; // container of all faces in object
  vec_e     e; // container of all edges in object
  Object                       (std::string);
  Object                       (Object const &);
  Object &     operator =      (Object const &);
  void         findVertAdj     (void);
  void         boundObject     (double * const &) const;
  void         setFF           (Face *);
  void         setFV           (Vertex *);
  Edge   *     getFE           (void) const;
  Face   *     getFF           (void) const;
  Vertex *     getFV           (void) const;
  std::string const & getName  (void) const;
  // build edges
  int          setNumDigits    (void) const;
  void         createEdges     (void);
  void         processEdge     (Face * const,
                                Vertex * const,
                                Vertex * const,
                                Vertex * const,
                                map_s_ep &,
                                int const &);
  void         createEdge      (Face * const,
                                Vertex * const,
                                Vertex * const,
                                Vertex * const,
                                map_s_ep &,
                                int const &);
  Edge *       findEdge        (Vertex const * const,
                                Vertex const * const,
                                map_s_ep const &,
                                int const &) const;
  std::string  keyPair         (int const & a,
                                int const & b,
                                int const & num_digits) const;
};

#endif
