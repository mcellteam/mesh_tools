#include <iostream>
#include <fstream>
#include <string>
#include <dirent.h>
#include <map>

using std::cout;
using std::endl;

#define FILENAME_SIZE 1024

////////////////////////////////////////////
////////////////////////////////////////////

struct lti
{
  bool operator()(const int s1, const int s2) const
  {
    return s1 < s2;
  }
};

typedef std::map<int,int,lti,std::allocator<int> > map_ii;
typedef std::map<int,int,lti,std::allocator<int> >::iterator ii_iterator;

void parseFace(char* triplet,int index[4])
{
  char val[80];
  char *eptr;
  char *cp=triplet;
  // get past 'Face'
  while (strchr("Face",*triplet)!=NULL) {triplet++;}
  // grab Face index
  while (strchr(" \t,",*triplet)!=NULL) { triplet++; }
  int i=0;
  while (strchr("0123456789+-eE.",*triplet)!=NULL)
  {
    val[i++] = *triplet++;
  }
  val[i]=0;
  index[0] = (int) strtod(val,&eptr);
  if (val==eptr)
  {
    printf("Error in reading face index\n");
    cout << "from line: " << cp << endl;
    exit(0);
  }
  // grab first vertex index
  while (strchr(" \t,",*triplet)!=NULL) { triplet++; }
  i=0;
  while (strchr("0123456789+-eE.",*triplet)!=NULL)
  {
    val[i++] = *triplet++;
  }
  val[i]=0;
  index[1] = (int)strtod(val,&eptr);
  if (val==eptr)
  {
    printf("Error in reading vertex index\n");
    cout << "from line: " << cp << endl;
    exit(0);
  }
  // grab second vertex index
  while (strchr(" \t,",*triplet)!=NULL) triplet++;
  i=0;
  while (strchr("0123456789+-eE.",*triplet))
  {
    val[i++] = *triplet++;
  }
  val[i]=0;
  index[2] = (int)strtod(val,&eptr);
  if (val==eptr)
  {
    printf("Error in reading vertex index\n");
    cout << "from line: " << cp << endl;
    exit(0);
  }
  // grab third vertex index
  while (strchr(" \t,",*triplet)!=NULL) triplet++;
  i=0;
  while (strchr("0123456789+-eE.",*triplet))
  {
    val[i++] = *triplet++;
  }
  val[i]=0;
  index[3] = (int)strtod(val,&eptr);
  if (val==eptr)
  {
    printf("Error in reading vertex index\n");
    cout << "from line: " << cp << endl;
    exit(0);
  }
}

int parseVertex(char* triplet)
{
  char val[80];
  char *eptr;
  char *cp=triplet;
  // skip leading whitespace
  while (strchr(" \t,",*triplet)!=NULL) { triplet++;}
  // get past 'Vertex'
  while (strchr("Vertx",*triplet)!=NULL) {triplet++;}
  // skip leading whitespace
  while (strchr(" \t,",*triplet)!=NULL) { triplet++;}
  // grab vertex index
  int i=0;
  while (strchr("0123456789+-eE.",*triplet)!=NULL)
  {
    val[i++] = *triplet++;
  }
  val[i]=0;
  int index = (int) strtod(val,&eptr);
  if (val==eptr)
  {
    index=0;
    printf("Error in reading vertex index\n");
    printf("Error in reading vertex: string %s\n",cp);
    exit(0);
  }
  return index;
}

int parseMeshalyzer(char* triplet)
{
  // example
  //Vertex <obj>d000.mesh<ind>25 [2.8 4.7 6.6] 
  char val[80];
  char *eptr;
  char *cp=triplet;
  // get to first '>'
  while (strchr(">",*triplet)==NULL) {triplet++;}
  triplet++;
  // get to second '>'
  while (strchr(">",*triplet)==NULL) {triplet++;}
  triplet++;
  // grab vertex index
  int i=0;
  while (strchr("0123456789+-eE.",*triplet)!=NULL)
  {
    val[i++] = *triplet++;
  }
  val[i]=0;
  int index = (int) strtod(val,&eptr);
  if (val==eptr)
  {
    index=0;
    printf("Error in reading vertex index\n");
    printf("line: %s\n",cp);
    exit(0);
  }
  return index;
}

////////////////////////////////////////////
////////////////////////////////////////////

int main(int argc, const char* argv[]){

  if(argc!=3)
  {
    cout << "\nError. Usage: remove_duplicate_vertices meshalyzer_output_file input_mesh_file\n"
          << "Before using remove_duplicate_vertices, first run meshalyzer with -p option on\n"
          << "input_mesh_file to generate meshalyzer_output_file. meshalyzer_output_file is\n"
          << "then scanned for keyphrase \"#  indistinguishable vertices: vertex \" after\n"
          << "which follows a list of indistinguishable vertices in input_mesh_file.\n"
          << "For each indistinguishable pair of vertices, all instances of the first vertex\n"
          << "are replaced with the second vertex. The new mesh file is written to STDOUT.\n"
          << "Note that indistinguishable vertices are assumed to come in pairs, and three or\n"
          << "more indistinguishable vertices will not be handled appropriately.\n";
    exit(0);
  }

  // open file
  std::ifstream inFile(argv[1]);
  if (inFile.fail()) // if stream cannot be opened
  {
    cout << "\n\nCan not open input meshalyzer file ["
          << argv[1] << "]\n\n";
    exit(1); 
  }

  //	map: old vertex index -> new vertex index
  map_ii id;
  bool found=false;
  bool o=true;
  int old=0;
  std::string start_pattern ="#  indistinguishable vertices: vertex";
  std::string line;
  while(getline(inFile,line))
  {
    if(found==false)
    {
      std::string::size_type loc = line.find( start_pattern, 0 );
      if( loc != std::string::npos )
      {
        found=true;
      }
    }
    if(found==true)
    {
      char *str=const_cast<char*>(line.c_str());
      // if first character is V for Vertex
      if ( line.empty()==false && strchr("V",*str)!=NULL)
      {
        int index = parseMeshalyzer(str);
        if(o==true){
          old=index;
          o=false;
        } else {
          id[old]=index;
          o=true;
        }
      }
    }
  }
  inFile.close();


  //	cout << "# vertex pairs = " << id.size() << endl;

  // open file
  inFile.open(argv[2]);
  if (inFile.fail()) // if stream cannot be opened
  {
    cout << "\n\nCan not open input mesh file ["
          << argv[2] << "]\n\n";
    exit(1); 
  }

  bool vert=true;

  while(getline(inFile,line))
  {
    if(vert==true)
    {
      // skip leading whitespace
      char *str=const_cast<char*>(line.c_str());
      while (strchr(" \t,",*str)!=NULL) { str++;}
      // if first character is V for Vertex
      if (strchr("V",*str)!=NULL)
      {
        int index = parseVertex(str);
        if(id.find(index)==id.end()){cout << line << endl;}
      } else{vert=false;}
    }
    if(vert==false)
    {
      // skip leading whitespace
      //			char *str=&line[0];
      char *str=const_cast<char*>(line.c_str());
      while (strchr(" \t,",*str)!=NULL) { str++;}
      // if first character is F for Face
      if (strchr("F",*str)!=NULL)
      {
        int index[4] = {0,0,0,0};
        parseFace(str,index);
        //				$faces_pattern="/Face\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+/";
        //				preg_match($faces_pattern, $i, $j);
        //				if ($j) {
        //					$a=$j[2];
        //					$b=$j[3];
        //					$c=$j[4];
        if(id.find(index[1])!=id.end()){index[1]=id[index[1]];}
        if(id.find(index[2])!=id.end()){index[2]=id[index[2]];}
        if(id.find(index[3])!=id.end()){index[3]=id[index[3]];}
        cout << "Face "
              << index[0] << " "
              << index[1] << " "
              << index[2] << " "
              << index[3] << endl;
        //				}
      }
    }
  }
}
