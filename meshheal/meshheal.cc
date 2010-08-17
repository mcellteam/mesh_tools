#include "meshheal.h"

#include <cstdlib>
#include <cstdarg>
#include <iostream>
#include <fstream>

#include "controls.h"
#include "object.h"

using std::cout;
using std::cerr;
using std::endl;

int main(int argc,char **argv)
{
  // instantiate controls class
  Controls & cs(Controls::instance());

  // parse command line 
  cs.parse(argc,argv,cs.getUsageMessage());

  // create object with vertices and faces
  Object o(cs.get_inpath());

  // create edges
  o.createEdges();

  // gather free vertices
  vec_v free_vertices;
  o.gatherFreeVertices(free_vertices);

  // project free_vertices onto random vector
  o.projectVerts(free_vertices);

  while(free_vertices.size()>0)
  {
    // compute distances between free vertices
    Closest distances;
    o.computeDistances(distances,free_vertices);

    // report progress
    cerr << "\nnumber of free vertices = " << free_vertices.size() << endl
          << "  squared separation distance"
          << "  between closest free vertices = " << distances.distance << endl
          << "  good vertex = ";
    distances.good_vertex->print(cerr);
    cerr << endl;
    cerr << "   bad vertex = ";
    distances.bad_vertex->print(cerr);
    cerr << endl;

    // enforce threshold
    if (distances.distance > cs.get_distance_threshold())
    {
      cerr << "\n\nError: minimum vertex separation distance ("
            << distances.distance << ") is greater than threshold ("
            << cs.get_distance_threshold() << ")\n";
      cerr << "\nnumber of free vertices = " << free_vertices.size() << endl;
      o.printMesh(cout);
      exit(1);
    }

    // remove second vertex in pair from vertex list
    o.removeVertex(distances.bad_vertex);

    // replace all instances of bad vertex in face list with good vertex
    o.fixFaces(distances.bad_vertex,distances.good_vertex);

    // rebuild edges
    o.rebuildEdges();

    // gather free vertices
    o.gatherFreeVertices(free_vertices);
   
    // detect end game
    if (free_vertices.size()==3)
    {
      o.addFinalFace(free_vertices);
      break;
    }
  }
  cerr << "\nnumber of free vertices = " << free_vertices.size() << endl;
  o.printMesh(cout);
  cerr << "Finished.\n";
}

std::string format (char const *fmt, ...)
{
  va_list args;
  va_start(args, fmt);
  std::string ret;
  char buffer[512];
  uint needLen = vsnprintf(buffer, sizeof(buffer), fmt, args);
  if (needLen <= sizeof(buffer))
  {
    ret = buffer;
  }
  else
  {
    char *big_buffer = new char[needLen + 1];
    vsnprintf(big_buffer, needLen + 1, fmt, args);
    ret = big_buffer;
    delete[] big_buffer;
  }
  va_end(args);
  return ret;
};

