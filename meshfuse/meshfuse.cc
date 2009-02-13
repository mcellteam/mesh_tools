#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <string>
#include <iostream>
#include <vector>
#include <map>
#include <list>

using std::cout;
using std::endl;

#define EPSILON 1E-10
#define PRINT false

#include "classes.cc"
#include "functions.cc"

const double random_vector[3] = {0.48565478568609,0.01968678962794,0.87392897849559};


int main(int argc,char *argv[])
{
  std::string message = "\n";
  message=message+
        "NAME\n"+
        "       meshfuse - fill holes in mesh by merging vertices\n"+
        "\nSYNOPSIS\n"+
        "       meshfuse FILE [OPTIONS]\n"+
        "\nDESCRIPTION\n"+
        "       Meshfuse merges pairs of border vertices to fill in.\n"+
        "       holes in the mesh. Closest vertices are merged first.\n"+
        "       Vertices must be on an open edge, i.e. border, to be merged.\n"+
        "       If an optional maximum vertex separation distance threshold\n"+
        "       is specfied, then vertices separated by more than the threshold\n"+
        "       will not be merged.\n\n"+
        "       Input mesh must pass meshalyzer integrity check,\n"+
        "       but mesh is allowed to be nonmanifold.\n"+
        "\nEXAMPLES\n"+
        "       meshfuse in_file.mesh > out_file.mesh\n"+
        "              Fuse all candidate border vertices in in_file.mesh.\n"+
        "       meshfuse in_file.mesh -t 1E-8 > out_file.mesh\n"+
        "              Fuse candidate border vertices in in_file.mesh\n"+
        "              separated by less than 1E-8 units.\n"+
        "\nOPTIONS\n"+
        "       -t NUM\n"+
        "              Do not fuse vertices separated by a distance larger than NUM.\n"+
        "              Units are same as FILE.\n\n"+
        "       -h\n"+
        "              Print meshfuse man page.\n\n"+
        "\nJustin Kinney				2008/05/20\n";

  std::string filename;
  const double user_threshold = parse(argc,argv,message,filename);

  // build object
  Object o(filename);

  // gather free vertices
  o.gatherFreeVertices();

  // get max edge length squared
  double max_edge_length_sq = o.findLongestEdge();

  // sort vertices by dot product with random vector

  // optimization
  // purge map of all edges containing no free vertices
  o.purgeMap();

  // optimization
  // purge vector of all edges containing no free vertices
  o.purgeEdges();

  // compute distances between free vertices
  o.computeDistances(user_threshold,max_edge_length_sq);

  ///// merge vertices /////
  while(o.distances.empty()==false)
  {
    // identify free vertex pair separated by the smallest distance
    Distance d = o.distances.back();
    o.distances.pop_back();
    fprintf(stderr,"Keeping vertex %d and removing vertex %d, separation distance %g\n",d.vA->index,d.vB->index,sqrt(d.d));
    fflush(stderr);

    // replace all instances of second vertex in face list with first vertex
    o.fixFaces(d.vA,d.vB);

    // update edges
    o.updateEdges(d.vA,d.vB);

    // remove second vertex in pair from vertex list
    o.fixVertices(d.vB);

    // gather free vertices
    o.gatherFreeVertices();

    // purge nonfree vertices from distances vector
    o.updateDistances();
  }
  // print to stdout
  o.printVerticesFaces();

  return 0;
}
