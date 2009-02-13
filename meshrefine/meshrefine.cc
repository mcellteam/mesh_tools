#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <iostream>
#include <fstream>

#include <vector>
#include <algorithm>
#include <ext/hash_map>
#include <ext/hash_set>
#include <map>
#include <float.h>

using std::cout;
using std::cerr;
using std::endl;

double THRESHOLD = 20;
std::string INPUT_FILE = "";
std::string FROZEN_FILE = "";
std::string MESH_SUFFIX = ".mesh";
#define LINE_SIZE       2048
#define FROZEN_LINE_SIZE  128
#define SUFFIX ".UPDATED"

typedef std::vector<std::string>::iterator vs;

#include "classes.cc"
#include "functions.cc"

int main(int argc,char *argv[])
{
  std::string message = "\n";
  message=message+
        "NAME\n"+
        "       meshrefine - shorten edges by adding vertices and faces\n"+
        "\nSYNOPSIS\n"+
        "       meshrefine input_file [options]\n"+
        "\nDESCRIPTION\n"+
        "       Polygons are iteratively bisected so that no edge is longer than threshold.\n"+
        "       If threshold is negative each edge is bisected once if edge length\n"+
        "       is longer than absolute value of threshold.\n"+
        "       New mesh is written to stdout. Input mesh must be fully closed and\n"+
        "       consistently oriented with vertex and face indices sequentially numbered.\n"+
        "\nEXAMPLES\n"+
        "       meshrefine input_file -t threshold -f frozen_vertices.dat\n"+
        "              Read mesh from file 'input_file' and write new refined\n"+
        "              mesh to stdout. Compare edge lengths to threshold.\n"+
        "              Use frozen_vertices.dat to freeze new vertices when appropriate.\n"+
        "\nOPTIONS\n"+
        "       -t NUM\n"+
        "              Edge length threshold:\n"+
        "              If threshold is positive then edges are iteratively bisected until\n"+
        "              no edge is longer than threshold. If threshold is negative then\n"+
        "              each edge is bisected once if edge length is longer than\n"+
        "              absolute value of threshold. Units are same as meshes in input mesh.\n"+
        "              Default is 20.\n\n"+
        "       -f FILE\n"+
        "              Frozen vertices file:\n"+
        "              Any vertices from input mesh found in FILE are considered frozen.\n"+
        "              An edge defined by two frozen vertices generates new frozen\n"+
        "              vertices when bisected.\n"+
        "              The format of FILE must be 'object_name vertex_index'.\n"+
        "              For example, 'd000 134' indicates that vertex number 134\n"+
        "              in object d000 should be frozen.\n"+
        "              Default behavior in the case that '-f' option is not used\n"+
        "              is to not freeze any vertices.\n\n"+
        "       -h\n"+
        "              Print meshrefine man page.\n"+
        "\nJustin Kinney				2008/05/01\n";

  bool freeze = parse(argc,argv,message);

  // build object
  Object o(INPUT_FILE);
  //cerr << "object name = " << o.name << endl;

  // identify frozen vertices
  //uint orig_frozen_size=0;
  if(freeze==true)
  {
    o.processFrozenFile(FROZEN_FILE);
  //  orig_frozen_size=o.frozen.size();
  //  cerr << "Original # frozen vertices = " << orig_frozen_size << endl;
  }

  int j=0;
  do {
    int a = static_cast<int>(o.v.size());
    int b = static_cast<int>(o.f.size());
    fprintf(stderr,"iteration %i: ",j++);
    fprintf(stderr,"max vertex %i max face %i\n",a,b);
    fprintf(stderr,"Building edges..................");
    fflush(stderr);
    o.createEdges();
    fprintf(stderr,"complete.\n");
    fflush(stderr);
  } while(refineMesh(o,THRESHOLD)==true);
  // print updated frozen vertices to new file
  if(freeze==true)
  {
  //  if(o.frozen.size()==orig_frozen_size)
  //  {
  //    cerr << "Error. No new frozen vertices identified!\n";
  //    exit(0);
  //  }
    o.printFrozen(FROZEN_FILE);
  }
  // print new mesh to stdout
  printMesh(o);
  return 0;
}
