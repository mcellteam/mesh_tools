#include <math.h>
#include <iostream>
#include <vector>
#include <set>
#include <list>
#include <algorithm>
#include <fstream>
#include <string>
#include <dirent.h>
#include <errno.h>
#include <time.h>
#include <sys/timeb.h>
#include <sys/time.h>
#include <ext/hash_map>
#include <ext/hash_set>
#include <numeric>
#include <map>
#include <cassert>

using std::map;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using std::left;
using std::right;
using std::set;

std::string INPUT_DATA_DIR  = "./";
std::string OUTPUT_DATA_DIR = "./";
std::string FROZEN_VERTICES_FILE = "";

const std::string VERTS_SUFFIX = "_verts.dat";
const std::string CONNECTIONS_SUFFIX = "_connections.dat";
const std::string COLORS_SUFFIX = "_colors.dat";
const std::string CP_FILE = "custom_points.dat";

//const double MAX_VD = 42;
//const double MID_VD = 21;
const double MAX_VD = 22;
const double MID_VD = 11;
const double MIN_VD = 0;
const double MAX_COLOR[3] = {1.0,1.0,0.0};
const double MID_COLOR[3] = {1.0,0.0,0.0};
const double MIN_COLOR[3] = {0.0,0.0,1.0};
const std::string DX_FILE = "vd.dx";

const double MAX_ERROR = 0;
const double MID_ERROR = -10;
const double MIN_ERROR = -20;
//const double MAX_COLOR[3] = {0.0,0.0,1.0};
//const double MID_COLOR[3] = {1.0,0.0,0.0};
//const double MIN_COLOR[3] = {1.0,1.0,0.0};
//const std::string DX_FILE = "ecw_error.dx";

const bool MEASURE_VD = true;

const double NULL_COLOR[3] = {0.1,0.1,0.1};

//####################################################
//#################### includes ######################
//####################################################

#include "controls.cc"
#include "classes.cc"
#include "subroutines.cc"

//####################################################
//###################### main  #######################
//####################################################

int main(int argc,char **argv)
{

  std::string message = "\n";
  message=message+
    "\nSyntax: measure_ecw_error -i input_dir -o output_dir -t target_ecw -f frozen_vertices_file\n\n"+
    "The extracellular width of the model in input_dir is\n"+
    "computed at each vertex and compared to target_ecw to\n"+
    "generate an ecw error which is mapped to vertex color.\n"+
    "A dx file (ecw_error.dx) containing the model with colored vertices is\n"+
    "written to the output directory.\n\n";

  bool freeze = parse(argc,argv,message);

  // check that assumption of 32 bit int is correct
  assert(checkIntSize());

  // save input directory
  //char filename[FILENAME_SIZE];
  //strcpy(filename,argv[1]);
  //char *temp=strrchr(filename,'/');
  //if(!temp) {strcat(filename,"/");}
  //else if(*++temp) {strcat(filename,"/");}
  //INPUT_DATA_DIR = filename;

  // read target ecw
  //char *eptr;
  //double t = strtod(argv[2],&eptr);

  // create container, objects, vertices, faces, edges, and find adjacencies
  Container c;
  
  // read frozen vertices
  if(freeze==true){ c.readFrozen(FROZEN_VERTICES_FILE.c_str()); }

  // print frozen vertices
  std::ofstream CP;
  std::string strcp = OUTPUT_DATA_DIR + CP_FILE;
  CP.open(strcp.c_str());
  if(CP.is_open()==false){
    fprintf(stderr,"\nCouldn't open output file %s\n",strcp.c_str());
    exit(0);
  }
  CP << "x_coordinate y_coordinate z_coordinate state_value x_normal y_normal z_normal\n";
  // for each object in container
  for(std::vector<Object*>::iterator i=c.o.begin();i!=c.o.end();i++)
  {
    // for each vertex in object
    for(std::vector<Vertex*>::iterator j=(*i)->v.begin();j!=(*i)->v.end();j++)
    {
      if(binary_search(c.frozen.begin(),c.frozen.end(),*j)==true)
      {
        CP << (*j)->pN[0] << " " << (*j)->pN[1] << " " << (*j)->pN[2] << " 1 0 0 1\n";
      }
    }
  }
  CP.close();

  // initialize space data structure
  Space s(c);

  // assign faces to boxes
  c.assignFacesToBoxes(s);

  // instantiate statistics class
  Monitor stats;

  // identify vertices that lie inside of another object
  c.findNice(s);

  if (MEASURE_VD==false)
  {
    // identify the closest point on a mesh to each vertex
    c.getSeparationDistancesAndECW(s,stats);
  }
  else
  {
    // identify the closest point on a mesh to each vertex
    c.getSeparationDistances(s,stats);
    
    // compute the vertex forces due to face intersections
    c.computeFaceIntersectionForce(s);

    // measure vd
    stats.groupInit(c);

    // print vd
    c.printVD(stats);
  }
}
