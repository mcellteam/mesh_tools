#include<stxxl.h>
#include<iostream>
#include <float.h>
#include <dirent.h>
#include <set>
#include <time.h>
#include <sys/timeb.h>
#include <sys/time.h>

using std::cout;
using std::endl;

//note: 2^64=18446744073709551616
typedef long long int         key;

// pi
#define PI 3.14159265358979

#define DOUBLE_EPSILON 1e-12

// number of bins in histograms
#define NUM_BINS 16

// default filename size
#define FILENAME_SIZE 1024

// max number of vertex adjacent faces
#define ADJACENT_FACES 40

// max number of faces per box
#define NUMBER_OF_FACES 24

// bytes
#define NODE_CACHE_SIZE (128*1024*1024) // 128MB
#define LEAF_CACHE_SIZE (128*1024*1024)

struct ltk
{
  bool operator()(const key s1, const key s2) const
  {
    return s1 < s2;
  }
};

struct pSort
{
  bool operator () (const std::pair<key,key> & lhs,
                    const std::pair<key,key> & rhs) const
  {
    if(lhs.first != rhs.first)
    {
      return lhs.first < rhs.first;
    }
    else
    {
      return lhs.second < rhs.second;
    }
  }
  static std::pair<key,key> max_value()
  {
    return std::make_pair(std::numeric_limits<key>::max(),
                          std::numeric_limits<key>::max());
  }
};

struct kSortComp
{
  bool operator () (const key & lhs,const key & rhs) const
  {
    return lhs < rhs;
  }
  static key max_value()
  {
    return std::numeric_limits<key>::max();
  }
};

typedef stxxl::vector<key>                  vec_key;
typedef stxxl::vector<key>::iterator        vk;
typedef stxxl::vector<double>               vec_d;
typedef stxxl::vector<double>::iterator     vd;
typedef std::set<key,ltk>                   set_key;
typedef std::set<key,ltk>::iterator         sk;
typedef std::pair<key,key>                  pkk;
typedef std::set<pkk,pSort>                 ne_set;
typedef std::set<pkk,pSort>::iterator       nes;
typedef	stxxl::map<pkk,key,pSort>           map_edge;
typedef	stxxl::map<pkk,key,pSort>::iterator m_e;
typedef stxxl::map<key,key,kSortComp>::const_iterator m_kk;

class Container;
class E;
class F;
class FE;
class NE;
class Object;
class Stats;
class V;
class VE;

int distinguishable(double a,double b)
{
  double c;
  c=a-b;
  if (c<0) c=-c;
  if (a<0) a=-a;
  if (a<1) a=1;
  if (b<0) b=-b;
  if (b<a) return (c>a*DOUBLE_EPSILON);
  else return (c>b*DOUBLE_EPSILON);
}

std::string message1 = "\n";
std::string message=message1+
"NAME\n"+

"       meshalyzer - mesh quality analyzer\n"+
"\nSYNOPSIS\n"+
"       meshalyzer [options] FILE|DIR\n"+
"\nDESCRIPTION\n"+
"       Meshalyzer is a general purpose mesh analyzer useful for\n"+
"       generating a complete summary of the current state of a mesh.\n"+
"       Meshalyzer assesses mesh integrity (e.g. missing data),\n"+
"       mesh attributes (e.g. closed, manifold, oriented), and\n"+
"       mesh characteristics (e.g. number of vertices, faces, edges).\n"+
"       Batch processing is easy by passing a directory name\n"+
"       as input on the command line.\n"+
"\nEXAMPLES\n"+
"       meshalyzer filename\n"+
"              Evaluate mesh integrity, attributes,\n"+
"              and characteristics for the single mesh file.\n\n"+
"       meshalyzer directoryname\n"+
"              Evaluate mesh integrity, attributes,\n"+
"              and characteristics for each single mesh file in directory.\n\n"+
"       meshalyzer -p filename\n"+
"              Evaluate mesh integrity, attributes, and characteristics \n"+
"              for the single mesh file, print the results, and print \n"+
"              the mesh elements preventing the mesh from being good with\n"+
"              regards to the mesh characteristics and the attributes, if any.\n\n"+
"       meshalyzer -a -p filename\n"+
"              Evaluate the five mesh attributes for the single mesh file,\n"+
"              print the state of each attribute, and print the mesh elements\n"+
"              preventing the mesh from being good with regards to the attributes, if any.\n\n"+
"       meshalyzer -b 10.0 -c 1.0 -p filename\n"+
"              Evaluate mesh integrity, attributes, and characteristics \n"+
"              for the single mesh file, print the results, and print \n"+
"              the mesh elements preventing the mesh from being good with\n"+
"              regards to the mesh characteristics and the attributes, if any.\n"+
"              Additionally, screen faces with aspect ratios larger than 10.0 and\n"+
"              screen edges with lengths larger than 1.0, and print detailed\n"+
"              information about the offending mesh elements.\n"+
"\nOPTIONS\n"+
"       -a\n"+
"              Evaluate the attributes of the mesh and report the results.\n"+
"              Skip the evaluation of mesh characteristics.\n\n"+
"       -b NUM\n"+
"              Detect edges with length smaller than NUM.\n"+
"              Units are same as FILE.\n\n"+
"       -c NUM\n"+
"              Detect edges with length greater than NUM.\n"+
"              Units are same as FILE.\n\n"+
"       -d NUM\n"+
"              Detect edges with angle between two adjacent faces\n"+
"              smaller than NUM degrees.\n\n"+
"       -e NUM\n"+
"              Detect edges with angle between two adjacent faces\n"+
"              greater than NUM degrees.\n\n"+
"       -f NUM\n"+
"              Detect faces with aspect ratio greater than NUM.\n\n"+
"       -h\n"+
"              Print meshalyzer man page.\n\n"+
"       -i\n"+
"              Detect intersections between faces from different objects.\n"+
"              Faceintersection detection is performed once all objects\n"+
"              are loaded into memory. Single object intersection detection\n"+
"              is omitted.\n\n"+
"       -p\n"+
"              Print detailed information about offending mesh elements \n"+
"              (i.e. flipped faces, borders, nonmanifold edges,\n"+
"              nonmanifold vertices, intersecting faces).\n\n"+
"       -q\n"+
"              Same as '-p' option, but prints vertex information\n"+
"              in dreamm custom points format.\n"+
"\nJustin Kinney				2007/10/01\n";

class Stats{
public:
  double total;
  double min;
  double max;
  double sum;
  double sum2;
  int count[NUM_BINS];
  double bins[NUM_BINS+1];
  vec_d x;
  //vd x_begin;
  key x_begin;
  key x_size;
  Stats(void);
  double mean(void);
  double median(void);
  double variance(void);
  void createHistogram(void);
  void createAdjacentFaceHistogram(void);
  void createAspectRatioHistogram(void);
  void printHistogram(void);
  void printAspectRatioHistogram(void);
  void printAdjacentFaceHistogram(void);
  void printStats(void);
  void clear();
  void process();
};

Stats::Stats(void)
{
  total=sum=sum2=0;
  min=DBL_MAX;
  max=-DBL_MAX;
  for(int i=0;i<NUM_BINS;i++)
  {
    count[i]=0;
    bins[i]=0;
  }
  bins[NUM_BINS]=0;
  x_begin=0;
}

void Stats::clear(void)
{
  total=sum=sum2=0;
  min=DBL_MAX;
  max=-DBL_MAX;
  for(int i=0;i<NUM_BINS;i++)
  {
    count[i]=0;
    bins[i]=0;
  }
  bins[NUM_BINS]=0;
}

void Stats::process()
{
  clear();
  vd j = x.begin();
  for(key i=0;i<x_begin;i++){j++;}
  for(vd i=j;i!=x.end();i++)
  {
    sum+=*i;
    sum2+=(*i)*(*i);
    total+=*i;
    if(*i<min) {min=*i;}
    if(*i>max) {max=*i;}
  }
}

void Stats::printHistogram(void)
{
  printf("  %9.4g - %-9.4g    : %9d  | %9.4g - %-9.4g   : %9d\n",
         bins[0], bins[1], count[0], bins[8], bins[9], count[8]);
  printf("  %9.4g - %-9.4g    : %9d  | %9.4g - %-9.4g   : %9d\n",
         bins[1], bins[2], count[1], bins[9], bins[10], count[9]);
  printf("  %9.4g - %-9.4g    : %9d  | %9.4g - %-9.4g   : %9d\n",
         bins[2], bins[3], count[2], bins[10], bins[11], count[10]);
  printf("  %9.4g - %-9.4g    : %9d  | %9.4g - %-9.4g   : %9d\n",
         bins[3], bins[4], count[3], bins[11], bins[12], count[11]);
  printf("  %9.4g - %-9.4g    : %9d  | %9.4g - %-9.4g   : %9d\n",
         bins[4], bins[5], count[4], bins[12], bins[13], count[12]);
  printf("  %9.4g - %-9.4g    : %9d  | %9.4g - %-9.4g   : %9d\n",
         bins[5], bins[6], count[5], bins[13], bins[14], count[13]);
  printf("  %9.4g - %-9.4g    : %9d  | %9.4g - %-9.4g   : %9d\n",
         bins[6], bins[7], count[6], bins[14], bins[15], count[14]);
  printf("  %9.4g - %-9.4g    : %9d  | %9.4g - %-9.4g   : %9d\n",
         bins[7], bins[8], count[7], bins[15], bins[16], count[15]);
}

void Stats::printAdjacentFaceHistogram(void)
{
  int L[16],U[16];
  for(int i=0;i<16;i++)
  {
    L[i]=1+static_cast<int>(bins[i]);
    U[i]=static_cast<int>(bins[i+1]);
  }
  printf("  %9d - %-9d    : %9d  | %9d - %-9d   : %9d\n",
         0, U[0], count[0], L[8], U[8], count[8]);
  printf("  %9d - %-9d    : %9d  | %9d - %-9d   : %9d\n",
         L[1], U[1], count[1], L[9], U[9], count[9]);
  printf("  %9d - %-9d    : %9d  | %9d - %-9d   : %9d\n",
         L[2], U[2], count[2], L[10], U[10], count[10]);
  printf("  %9d - %-9d    : %9d  | %9d - %-9d   : %9d\n",
         L[3], U[3], count[3], L[11], U[11], count[11]);
  printf("  %9d - %-9d    : %9d  | %9d - %-9d   : %9d\n",
         L[4], U[4], count[4], L[12], U[12], count[12]);
  printf("  %9d - %-9d    : %9d  | %9d - %-9d   : %9d\n",
         L[5], U[5], count[5], L[13], U[13], count[13]);
  printf("  %9d - %-9d    : %9d  | %9d - %-9d   : %9d\n",
         L[6], U[6], count[6], L[14], U[14], count[14]);
  printf("  %9d - %-9d    : %9d  | %9d - %-9d   : %9d\n",
         L[7], U[7], count[7], L[15], U[15], count[15]);
}

void Stats::printAspectRatioHistogram(void)
{
  printf("  %9.4g - %-9.4g    : %9d  | %9.4g - %-9.4g   : %9d\n",
         bins[0], bins[1], count[0], bins[8], bins[9], count[8]);
  printf("  %9.4g - %-9.4g    : %9d  | %9.4g - %-9.4g   : %9d\n",
         bins[1], bins[2], count[1], bins[9], bins[10], count[9]);
  printf("  %9.4g - %-9.4g    : %9d  | %9.4g - %-9.4g   : %9d\n",
         bins[2], bins[3], count[2], bins[10], bins[11], count[10]);
  printf("  %9.4g - %-9.4g    : %9d  | %9.4g - %-9.4g   : %9d\n",
         bins[3], bins[4], count[3], bins[11], bins[12], count[11]);
  printf("  %9.4g - %-9.4g    : %9d  | %9.4g - %-9.4g   : %9d\n",
         bins[4], bins[5], count[4], bins[12], bins[13], count[12]);
  printf("  %9.4g - %-9.4g    : %9d  | %9.4g - 10000       : %9d\n",
         bins[5], bins[6], count[5], bins[13], count[13]);
  printf("  %9.4g - %-9.4g    : %9d  |     10000 - 100000      : %9d\n",
         bins[6], bins[7], count[6], count[14]);
  printf("  %9.4g - %-9.4g    : %9d  |    100000 -             : %9d\n",
         bins[7], bins[8], count[7], count[15]);
}

void Stats::printStats(void)
{
  cout << "       min      " << min << endl;
  cout << "       max      " << max << endl;
  cout << "       median   " << median() << endl;
  cout << "       mean     " << mean() << endl;
  cout << "       variance " << variance() << endl << endl;
}

void Stats::createAdjacentFaceHistogram(void)
{
  // create bins
  bins[0]=0;
  bins[1]=0.99;
  bins[2]=1.99;
  bins[3]=2.99;
  bins[4]=3.99;
  bins[5]=4.99;
  bins[6]=5.99;
  bins[7]=6.99;
  bins[8]=7.99;
  bins[9]=8.99;
  bins[10]=9.99;
  bins[11]=10.99;
  bins[12]=11.99;
  bins[13]=12.99;
  bins[14]=13.99;
  bins[15]=14.99;
  if(max>15.99){bins[16]=max;}
  else{bins[16]=15.99;}
  // bin data
  unsigned int start,mid,end;
  x_size=0;
  // for each element in x from x_begin to x.end
  vd i = x.begin();
  for(key j=0;j<x_begin;j++){i++;}
  for(vd j=i;j!=x.end();j++)
  {
    start=0;
    end=NUM_BINS-1;
    while (end-start>1){
      mid=(start+end)/2;
      if (*j >= bins[start] && *j <= bins[mid]) { end=mid;}
      else { start=mid; }
    }
    if (*j >= bins[start] && *j <= bins[end]) { count[start]++;}
    else { count[end]++; }
    x_size++;
  }
}

void Stats::createAspectRatioHistogram(void)
{
  // create bins
  bins[0]=1.1547;
  bins[1]=1.5;
  bins[2]=2;
  bins[3]=2.5;
  bins[4]=3;
  bins[5]=4;
  bins[6]=6;
  bins[7]=10;
  bins[8]=15;
  bins[9]=25;
  bins[10]=50;
  bins[11]=100;
  bins[12]=300;
  bins[13]=1000;
  bins[14]=10000;
  bins[15]=100000;
  // bin data
  unsigned int start,mid,end;
  // for each element in x from x_begin to x.end
  x_size=0;
  vd ii = x.begin();
  for(key j=0;j<x_begin;j++){ii++;}
  for(vd i=ii;i!=x.end();i++)
  {
    start=0;
    end=NUM_BINS-1;
    while (end-start>1)
    {
      mid=(start+end)/2;
      if (*i >= bins[start] && *i <= bins[mid]) { end=mid;}
      else { start=mid; }
    }
    if (*i >= bins[start] && *i <= bins[end]) { count[start]++;}
    else { count[end]++; }
    x_size++;
  }
}

struct SortDec // comparison function
{
  bool operator() ( const double & a,
                    const double & b) const
  {
    return a<b;
  }
  static double min_value()
  {
    return std::numeric_limits<double>::min();
  }
  static double max_value()
  {
    return std::numeric_limits<double>::max();
  }
};

double Stats::mean(void)
{
  return sum/x_size;
}

double Stats::median(void)
{
  // bound the main memory consumption by M
  // during sorting
  const unsigned M = 1024*1024*1024; // bytes
  // sort records by callers number
  vd ii = x.begin();
  for(key j=0;j<x_begin;j++){ii++;}
  stxxl::sort(ii,x.end(),SortDec(),M);
  int i= x_size/2;
  return x[i];
}

double Stats::variance(void)
{
  if(x_size>1)
  {
    return (sum2/x_size-(mean()*mean()));
  }
  else
  {
    return 0;
  }
}

void Stats::createHistogram()
{
  // compute x size
  x_size=0;
  vd ii = x.begin();
  for(key j=0;j<x_begin;j++){ii++;}
  for(vd i=ii;i!=x.end();i++){x_size++;}
  double bin_width=(4.0*sqrt(variance()))/(NUM_BINS-1.0);
  bool nonneg=false;
  double mm = mean();
  while(nonneg==false)
  {
    // create bins
    bins[0]=0;
    bins[1]=0;
    for(int i=2;i<NUM_BINS;i++)
    {
      bins[i]=mm+(((i-1)-(NUM_BINS/2.0))*bin_width);	
    }
    bins[NUM_BINS]=max;
    nonneg=true;
    //check
    for(int i=2;i<NUM_BINS+1;i++)
    {
      if(bins[i]<0)
      {
        nonneg=false;
      }
    }
    if(nonneg==false)
    {
      // increase mm
      mm = mm*1.1;
    }
  }
  // bin data
  unsigned int start,mid,end;
  // for each element in x from first to last
  ii = x.begin();
  for(key j=0;j<x_begin;j++){ii++;}
  for(vd i=ii;i!=x.end();i++)
  {
    start=0;
    end=NUM_BINS-1;
    while (end-start>1)
    {
      mid=(start+end)/2;
      if (*i >= bins[start] && *i <= bins[mid]) { end=mid;}
      else { start=mid; }
    }
    if (*i >= bins[start] && *i <= bins[end]) { count[start]++;}
    else { count[end]++; }
  }
}

/* *********************************************** */
/* *********************************************** */
/* *********************************************** */

struct Edata // edges
{
  key v1,v2; // relative
  key f1,f2; // relative
};

struct Fdata // faces
{
  key index;    // native, i.e. from file
  key v[3];	// relative
};

struct FEdata
{
  key e[3];	// always three edges per face, relative
};

struct ObjectData
{
  // absolute indices
  key v_start, v_end;
  key f_start, f_end;
  key e_start, e_end;
  key ve_start, ve_end;
  key fe_start, fe_end;
};

struct Vdata // vertices
{
  key index;    // native, i.e. from file
  double x, y, z;
};

struct VEdata
{
  key adj_f[ADJACENT_FACES]; // indices of adjacent faces to vertex, relative
                             // define first element of array
                             // to store the array index of
                             // the next available (i.e. empty)
                             // element. so adj_f[0] must be
                             // larger than 0.
};

vec_key border;     // keys to edges with only one face, relative
key     degen_f;    // # degenerate faces with duplicate vertex indices
key     disc_v;     // first discontiguous Vertex index 
key     dupl_f;     // number of duplicate Face indices
key     dupl_v;	    // number of duplicate Vertex indices
vec_key flipped;    // keys to flipped edges, relative
set_key found;      // detect duplicate face indices, relative
vec_key indistin_e; // keys to edges with indistinguishable vertices, relative
vec_key indistin_v; // keys to indistinguishable vertices, relative
key     miss_v;	    // number of missing vertices discovered
ne_set  nonm_e;     // nonmanifold edges: pair<Edge key,face key>, relative
set_key orphan;     // vertices referenced by no face, relative
bool    seq_v;      // true=sequential Vertex indexing,false=otherwise
std::vector<Object>       objects;

//WORKStypedef stxxl::VECTOR_GENERATOR<key,16,32,8388608>::result new_vec_key;
//WORKStypedef stxxl::VECTOR_GENERATOR<key,16,32,8388608>::result::iterator new_vk;
typedef stxxl::VECTOR_GENERATOR<key,4,8,8388608>::result new_vec_key;
typedef stxxl::VECTOR_GENERATOR<key,4,8,8388608>::result::iterator new_vk;
new_vec_key nonman_e;   // keys to edges with more than two faces, relative

typedef stxxl::VECTOR_GENERATOR<key,4,8,2097152>::result new_vec_key2;
typedef stxxl::VECTOR_GENERATOR<key,4,8,2097152>::result::iterator new_vk2;
new_vec_key2 nonman_v;   // keys to nonmanifold vertices, relative

//typedef stxxl::VECTOR_GENERATOR<Edata,4,8,16777216>::result new_edata;
typedef stxxl::VECTOR_GENERATOR<Edata,4,8,33554432>::result new_edata;
new_edata e_vector;
const new_edata & cev = e_vector;

//typedef stxxl::VECTOR_GENERATOR<Fdata,4,8,16777216>::result new_fdata;
typedef stxxl::VECTOR_GENERATOR<Fdata,4,8,33554432>::result new_fdata;
new_fdata f_vector;
const new_fdata & cfv = f_vector;

typedef stxxl::VECTOR_GENERATOR<FEdata,4,8,33554432>::result new_fedata;
new_fedata fe_vector;
const new_fedata & cfev = fe_vector;

stxxl::vector<ObjectData> object_vector;
const stxxl::vector<ObjectData> & cov = object_vector;

//stxxl::vector<Vdata>      v_vector;
typedef stxxl::VECTOR_GENERATOR<Vdata,4,8,33554432>::result new_vdata;
new_vdata v_vector;

//typedef stxxl::VECTOR_GENERATOR<VEdata,4,16,16777216>::result new_vedata;
typedef stxxl::VECTOR_GENERATOR<VEdata,4,8,33554432>::result new_vedata;
new_vedata ve_vector;

key vertex_index; // relative
key face_index;   // relative
key edge_index;   // relative

// command line globals
std::string inpath,outpath;

bool distvert;
// true if '-w', false otherwise
// true=end execution after printing indistinguishable vertices
// false=do nothing

bool folder;
// true if folder, false otherwise

bool attr;
// true if '-a', false otherwise
// true=eval attributes only
// false=allow a full analysis, i.e. evaluate characteristics

bool my_print;
// true if '-p', false otherwise
// true=print detailed information about offending mesh elements
// false=print summary information only

char style[32];
// output style for offending mesh elements (-p option)
// cp=dreamm custom points format, i.e. x y z 1 0 0 1
// everything else=detailed face,vertex,edge information

bool interf;
// true if '-i', false otherwise
// true=collect all objects,detect intersections between objects,
//	in contrast to the intra-object search that is the default
// false=do not perform inter-object face intersection check

bool my_signal[5];
// true if threshold value found on command line

double thresholds[5];
// user-defined thresholds found on command line
// [0]=aspect ratio threshold set
// [1]=min edge angle threshold set
// [2]=max edge angle threshold set
// [3]=min edge length threshold set
// [4]=max edge length threshold set
// TODO add more parameters to threshold, e.g. face area

/* *********************************************** */
/* *********************************************** */
/* *********************************************** */

class Container
{
public:
  void scanDir(const char *filename);
  void update();
  void printCumulative();
  void printIntegrity();
  void printAttr();
  void printChars();
  void vertexAdjacentFaces();
  void areaAspectRatio();
  void processEdgeLengths();
  void computeEdgeAngles();
  Object* processFile(std::string filename);

  void addFile(std::string s)
  {
    files.push_back(s);
    num_files++;
  }
  void a_create ()
  {
    area.createHistogram();
  }

  void af_create ()
  {
    adjacent_face.createAdjacentFaceHistogram();
  }

  void ar_create ()
  {
    aspect_ratio.createAspectRatioHistogram();
  }

  void ea_create ()
  {
    edge_angle.createHistogram();
  }

  void el_create ()
  {
    edge_length.createHistogram();
  }

  double area_sum()
  {
    return area.sum;
  }
  
  void clearData()
  {
    next_v=next_f=next_e=0;
    next_ve=next_fe=0;
  }

  Container()
  {
    good_integrity=true;
    attr=distvert=folder=my_print=interf=false;
    bb[0]=bb[1]=bb[2]=1E30;
    bb[3]=bb[4]=bb[5]=-1E30;
    for(int i=0;i<5;i++){my_signal[i]=false;thresholds[i]=0.0;}
    next_e = 0; // absolute
    next_f = 0; // absolute
    next_fe = 0; // absolute
    next_v = 0; // absolute
    next_ve = 0; // absolute
    num_bor=0;
    num_bou=0;
    num_clo=0;
    num_op=0;
    num_cons[0]=num_cons[1]=num_cons[2]=0;
    num_deg=0;
    num_dupl_f=0;
    num_dupl_v=0;
    num_e=0;
    num_files=0;
    num_flip=0;
    num_f=0;
    num_indistin=0;
    num_man[0]=num_man[1]=num_man[2]=0;
    num_mis=0;
    num_nonman_e=0;
    num_nonman_v=0;
    num_objects=0;
    num_orph=0;
    num_out[0]=num_out[1]=num_out[2]=0;
    num_v=0;
    object_id = 0;
    vol=0.0;
  }

  void countObj()
  {
    num_objects++;
  }
  
  std::string getFile(key i)
  {
    return files[i];
  }

  key getNextE()
  {
    return next_e; // absolute
  }

  key getNextF()
  {
    return next_f; // absolute
  }

  key getNextFE()
  {
    return next_fe; // absolute
  }

  key getNextV()
  {
    return next_v; // absolute
  }

  key getNextVE()
  {
    return next_ve; // absolute
  }

  key getNumFiles()
  {
    return num_files;
  }

  key getObjectId()
  {
    return object_id;
  }

  void incObjectId()
  {
    object_id++;
   }

  void prepFEdata(key oid)
  {
    object_vector[oid].fe_start = getNextFE(); // absolute
    // for each face in object
    key num=object_vector[oid].f_end-object_vector[oid].f_start+1;
    for(key i=0;i<num;i++) // relative
    {
      FEdata fed;
      fed.e[0]=fed.e[1]=fed.e[2]=-2;
      fe_vector.push_back(fed);
    }
    setNextFE(object_vector[oid].fe_start + num); // absolute+relative
    object_vector[oid].fe_end = getNextFE()-1; // absolute
  }

  void print_adjacent_face_histo()
  {
    adjacent_face.printAdjacentFaceHistogram();
  }

  void print_adjacent_face_stats()
  {
    adjacent_face.printStats();
  }

  void print_area_histo()
  {
    area.printHistogram();
  }

  void print_area_stats()
  {
    area.printStats();
  }

  void print_aspect_ratio_histo()
  {
    aspect_ratio.printAspectRatioHistogram();
  }

  void print_aspect_ratio_stats()
  {
    aspect_ratio.printStats();
  }

  void print_edge_angle_histo()
  {
    edge_angle.printHistogram();
  }

  void print_edge_angle_stats()
  {
    edge_angle.printStats();
  }

  void print_edge_length_histo()
  {
    edge_length.printHistogram();
  }

  void print_edge_length_stats()
  {
    edge_length.printStats();
  }
  
  void setBB(double d[6])
  {
    for(int i=0;i<6;i++){ bb[i]=d[i]; }
  }

  void setIntegrity(bool b)
  {
    good_integrity = b;
  }

  void setNextE(key i)
  {
    next_e=i; // absolute
  }

  void setNextF(key i)
  {
    next_f=i; // absolute
  }

  void setNextFE(key i)
  {
    next_fe=i; // absolute
  }

  void setNextV(key i)
  {
    next_v=i; // absolute
  }

  void setNextVE(key i)
  {
    next_ve=i; // absolute
  }

  void setNumE(key i)
  {
    num_e+=i;
  }

  void setNumF(key i)
  {
    num_f+=i;
  }

  void setNumV(key i)
  {
    num_v+=i;
  }

  void update_x_begin()
  {
    adjacent_face.x_begin=adjacent_face.x.size();
    area.x_begin=area.x.size();
    aspect_ratio.x_begin=aspect_ratio.x.size();
    edge_angle.x_begin=edge_angle.x.size();
    edge_length.x_begin=edge_length.x.size();
  }

private:
  double bb[6];     // bounding box of single mesh object
                    // [xmin,ymin,zmin,xmax,ymax,zmax]
  std::vector<std::string> files;	// array of input file names
  key next_e;       // absolute
  key next_f;       // absolute
  key next_fe;      // absolute
  key next_v;       // absolute
  key next_ve;      // absolute
  key num_bor;      // # edges with only one adjacent face, i.e. boundary edges
  key num_bou;      // # separate boundaries in object, i.e. # holes
  key num_clo;      // # closed objects
  key num_op;       // # open objects
  key num_cons[3];  // 0: # manifold objects with consistently oriented faces
                    // 1: # manifold objects with INconsistently oriented faces
                    // 2: # NONmanifold objects, irrespective of face orientation
                    //    i.e. # objects with undefined orientation
  key num_deg;      // # degenerate faces, i.e. faces without 3 unique vertex indices
  key num_dupl_v;   // # vertices with duplicate indices
  key num_dupl_f;   // # faces with duplicate indices
  key num_e;        // # edges in data set
  key num_files;    // # input files
  key num_flip;     // # edges with at least one adjacent face oriented incorrectly
  key num_f;        // # faces in data set
  key num_indistin; // # vertices indistinguishable from another vertex in x,y,z
  key num_man[3];   // 0: # manifold objects
                    // 1: # NONmanifold, closed objects
                    // 2: # NONmanifold, open objects
  key num_mis;      // # missing vertices, i.e. face references missing vertex index
  key num_nonman_e; // # edges with three or more adjacent faces
  key num_nonman_v; // # vertices whose adjacent faces are nontraversable
  key num_objects;  // # objects in data set
  key num_orph;     // # vertices not referenced by any face
  key num_out[3];   // 0: # objects with consistently-, outward-oriented faces
                    // 1: # obejcts with consistently-, inward-oriented faces
                    // 2: # objects with INconsistently oriented faces
  key num_v;        // # vertices in data set
  key object_id;    // current object id
  double vol;       // volume inside data set objects
  bool good_integrity;
  // cumulative statistics
  Stats adjacent_face;
  Stats area;
  Stats aspect_ratio;
  Stats edge_angle;
  Stats edge_length;
};

class E
{
public:
  void updateEdge(key f,key object_id);
  void print(key object_id);
  void printCP(key object_id);
  key getNewFace(key,key,key);
  bool isConsistent(key);
  bool getOrientation(key,key);
  double getSqLength(key object_id);
  double getAngle(key);
  void getVertices(key&,key&,key&,key&);

  E(key i)
  {
    e_id=i; //  absolute
  }

  key get_v1() // relative
  {
    //return e_vector[e_id].v1; 
    return cev[e_id].v1; 
  }

  key get_v2() // relative
  {
    //return e_vector[e_id].v2;
    return cev[e_id].v2;
  }

  key get_f1(key vector_key) // relative
  {
    //return e_vector[e_id].f1;
    if(get_v1()==vector_key)
    {
      return cev[e_id].f1;
    }
    else
    {
      return cev[e_id].f2;
    }
  }

  key get_f2(key vector_key) // relative
  {
    //return e_vector[e_id].f2;
    if(get_v1()==vector_key)
    {
      return cev[e_id].f2;
    }
    else
    {
      return cev[e_id].f1;
    }
  }

  key get_rel_key(key object_id) // relative
  { // absolute-absolute=relative
    //return e_id-object_vector[object_id].e_start;
    return e_id-cov[object_id].e_start;
  }

  bool isManifold(key object_id)
  {
    // look for 3rd and higher edges stored in array
    std::pair<key,key> p;
    p = std::make_pair(get_rel_key(object_id),get_rel_key(object_id));
    return nonm_e.find(p)==nonm_e.end();
  }

private:
  key e_id; // absolute
};

class F
{
public:
  key getNewEdge(key,key,key);
  void print(key object_id);
  void getNormal(key,double[3]);
  double getAspectRatio(key);
  key get_v1_file(key object_id); // from file
  key get_v2_file(key object_id); // from file
  key get_v3_file(key object_id); // from file

  F(key i)
  {
    f_id=i; // absolute
  }

  F()
  {
    f_id=0; // absolute
  }

  void assign(key i) // absolute
  {
    f_id=i;
  }

  key get_o(key v1,key v2)
  {
    if(get_v1_rel()!=v1 && get_v1_rel()!=v2) return get_v1_rel();
    if(get_v2_rel()!=v1 && get_v2_rel()!=v2) return get_v2_rel();
    if(get_v3_rel()!=v1 && get_v3_rel()!=v2) return get_v3_rel();
    cout << "\nError: other vertex not found!\n"; exit(0);
  }

  double get_index() // from file
  {
    return f_vector[f_id].index;
  }

  key get_v1_rel() // relative in v_vector
  {
    //return f_vector[f_id].v[0];
    return cfv[f_id].v[0];
  }

  key get_v2_rel() // relative in v_vector
  {
    //return f_vector[f_id].v[1];
    return cfv[f_id].v[1];
  }

  key get_v3_rel() // relative in v_vector
  {
    //return f_vector[f_id].v[2];
    return cfv[f_id].v[2];
  }

  key get_abs_key() // absolute
  {
    return f_id;
  }

  key get_rel_key(key object_id) // relative
  {
    return f_id-object_vector[object_id].f_start;
  }

private:
  key f_id; // absolute
};

class FE
{
public:
  FE(key i) // absolute
  {
    fe_id=i;
  }

  void addEdge(key edge) // relative
  {
    for(int i=0;i<3;i++)
    {
      if(fe_vector[fe_id].e[i]<0)
      {
        fe_vector[fe_id].e[i]=edge;
        break;
      }
    }
  }

  key* getEdges()
  {
    return &fe_vector[fe_id].e[0];
  }

  key get_edge(int i)
  {
    //return fe_vector[fe_id].e[i];
    return cfev[fe_id].e[i];
  }

private:
  key fe_id; // absolute
};

class NE
{
public:

  void addPair(key edge,key face) // relative,relative
  {
    nonm_e.insert(std::make_pair(edge,face));
  }

  nes findMatch(key edge) // relative
  { 
    return nonm_e.find(std::make_pair(edge,edge));
  }

};

class V
{
public:
  bool isManifold(key,bool);
  void getAdjacentEdges(key,std::vector<key>&);
  bool scanAdjFaces(key,key,key,bool&);
  void printCP();
  void print(key object_id);
  void printAdjacent(key object_id);

  V(key i)
  {
    v_id=i; // absolute
  }

  key get_index()
  {
    return v_vector[v_id].index;
  }

  double get_x()
  {
    return v_vector[v_id].x;
  }

  double get_y()
  {
    return v_vector[v_id].y;
  }

  double get_z()
  {
    return v_vector[v_id].z;
  }

  key get_rel_key(key object_id) // absolute
  { // absolute-absolute=relative
    //return v_id-object_vector[object_id].v_start;
    return v_id-cov[object_id].v_start;
  }

  key get_abs_key()
  {
    return v_id;
  }

private:
  key v_id; // absolute
};

class VE
{
// access arrays of adjacent faces for vertices
public:

  VE(key i) // absolute
  {
    ve_id=i;
  }

  VE()
  {
    ve_id=0;
  }

  void assign(key i) // absolute 
  {
    ve_id=i;
  }

  void print(key object_id);

  key* get_array()
  {
    return &ve_vector[ve_id].adj_f[0];
  }

  void getAf(std::vector<key> &vv)
  {
    //for(int i=0;i<ADJACENT_FACES;i++)
    for(int i=1;i<ADJACENT_FACES;i++)
    {
      if(ve_vector[ve_id].adj_f[i]<0) break;
      vv.push_back(ve_vector[ve_id].adj_f[i]);
    }
  }

/*  unsigned int size()
  {
    for(int i=0;i<ADJACENT_FACES;i++)
    {
      if(ve_vector[ve_id].adj_f[i]<0) return i;
    }
    return ADJACENT_FACES;
  }*/

  unsigned int size()
  {
      return ve_vector[ve_id].adj_f[0]-1;
  }

/*  int get_next()
  {
    int i;
    for(i=0;i<ADJACENT_FACES;i++)
    {
      if(ve_vector[ve_id].adj_f[i]<0){break;}
    }
    return i;
  }
  */

  int get_next()
  {
    return ve_vector[ve_id].adj_f[0];
  }

  void set_next(key f,int next)
  {
    ve_vector[ve_id].adj_f[next] = f;
    ve_vector[ve_id].adj_f[0]++;
  }

  void addFace(key f)
  {
    VEdata* ved = &ve_vector[ve_id];
    int next = ved->adj_f[0]++;
    // if array is full
    if(next==ADJACENT_FACES)
    {
      cout << "\nError: increase ADJACENT_FACES (current value = "
            << ADJACENT_FACES << ").\n";
      cout << "face relative key = " << f << endl;
      exit(0);
    }
    else
    {
      // the next spot is empty
      ved->adj_f[next] = f;
      //ved->adj_f[0]++;
    }
  }

private:
  key ve_id; // absolute
};

class Object
{
public:
  bool isClosed();
  void evalAttributes(Container&);
  void evalCharacteristics(Container&);
  bool verticesManifold(bool);
  bool edgesManifold(bool);
  bool isManifold();
  bool isConsistent();
  void vertexDistin();
  void boundObject();
  void computeVolume();
  void computeGenus();
  void print(Container&);
  void printAttr(Container&);
  void printIntegrity(Container&);
  void printChars(Container&);
  void analyze(Container&);
  VE *GetVE(key id); // relative
  V *GetV(key id); // relative
  void findVertexAdjacencies(Container &c);
  VE ve;
  F f;
  
  Object (key i,std::string filename)
  {
    object_id = i;
    name = filename;
    bb[0]=bb[1]=bb[2]=1E30;
    bb[3]=bb[4]=bb[5]=-1E30;
    closed=consistent=outward=manifold=false;
  }
  
  E *GetE(key id); // relative

  NE *GetNE() // relative
  {
    return new NE();
  }

  FE *GetFE(key id) // relative
  {
    key fe_id = object_vector[object_id].fe_start + id; // absolute
    return new FE(fe_id);
  }

  F *GetF(key id) // relative
  {
    key f_id = object_vector[object_id].f_start + id; // absolute
    return new F(f_id);
  }

  std::string GetName()
  {
    return name;
  }

  bool getManifold()
  {
    return manifold;
  }

  bool getClosed()
  {
    return closed;
  }

  key getBoundary()
  {
    return num_bou;
  }

  bool getConsistent()
  {
    return consistent;
  }

  bool getOutward()
  {
    return outward;
  }

  double getVol()
  {
    return vol;
  }

private:
  double bb[6];         // bounding box [xmin,ymin,zmin,xmax,ymax,zmax]
  bool closed;		// true=closed mesh
  bool consistent;	// true=consistently-oriented face normals
  int genus;		// object genus
  bool manifold;	// true=2d manifold in R3
  std::string name;
  int num_sep;		// number of components (separate meshes in object)
  int num_bou;		// number of separate boundaries in object, same as # holes
  key object_id;
  bool outward;		// true=outwardly-oriented face normals
  double vol;		// object volume
};

/* *********************************************** */
/* *********************************************** */
/* *********************************************** */

void V::printAdjacent(key object_id)
{
  VE* ve_ptr = objects[object_id].GetVE(get_rel_key(object_id));
  ve_ptr->print(object_id);
  delete ve_ptr;
}

void E::print(key object_id)
{
  cout.precision(12);
  cout << "Object " << objects[object_id].GetName() << endl
        << "Edge :\n";
  V* v_ptr = objects[object_id].GetV(get_v1());
  cout  << "Vertex " << v_ptr->get_index()
        << " " << v_ptr->get_x()
        << " " << v_ptr->get_y()
        << " " << v_ptr->get_z()
        << endl;
  delete v_ptr;
  v_ptr = objects[object_id].GetV(get_v2());
  cout  << "Vertex " << v_ptr->get_index()
        << " " << v_ptr->get_x()
        << " " << v_ptr->get_y()
        << " " << v_ptr->get_z()
        << endl;
  delete v_ptr;
  //F* f_ptr = objects[object_id].GetF(get_f1());
  objects[object_id].f.assign(get_f1(get_v1()));
  F* f_ptr = &objects[object_id].f;
  cout  << "Face " << f_ptr->get_index()
        << " " << f_ptr->get_v1_file(object_id)
        << " " << f_ptr->get_v2_file(object_id)
        << " " << f_ptr->get_v3_file(object_id)
        << endl;
  //delete f_ptr;
  //f_ptr = objects[object_id].GetF(get_f2());
  objects[object_id].f.assign(get_f2(get_v1()));
  f_ptr = &objects[object_id].f;
  cout  << "Face " << f_ptr->get_index()
        << " " << f_ptr->get_v1_file(object_id)
        << " " << f_ptr->get_v2_file(object_id)
        << " " << f_ptr->get_v3_file(object_id)
        << endl;
  //delete f_ptr;
}

void E::printCP(key object_id)
{
  cout.precision(12);
  V* v_ptr = objects[object_id].GetV(get_v1());
  v_ptr->printCP();
  delete v_ptr;
  v_ptr = objects[object_id].GetV(get_v2());
  v_ptr->printCP();
  delete v_ptr;
}

void E::updateEdge(key f,key object_id)
{
  if (e_vector[e_id].f2<0)
  {
    // add 2nd face to edge
    e_vector[e_id].f2=f;
  }
  else
  {
    // save 3rd or greater edge as edge,face pair
    NE* ne_ptr = objects[object_id].GetNE();
    ne_ptr->addPair(get_rel_key(object_id),f);
    delete ne_ptr;
  }
  // add edge pointer to face
  FE* fe_ptr = objects[object_id].GetFE(f);
  fe_ptr->addEdge(get_rel_key(object_id));
  delete fe_ptr;
}

key F::get_v1_file(key object_id) // from file
{
  V *v_ptr = objects[object_id].GetV(f_vector[f_id].v[0]);
  key i = v_ptr->get_index();
  delete v_ptr;
  return i;
}

key F::get_v2_file(key object_id) // from file
{
  V *v_ptr = objects[object_id].GetV(f_vector[f_id].v[1]);
  key i = v_ptr->get_index();
  delete v_ptr;
  return i;
}

key F::get_v3_file(key object_id) // from file
{
  V *v_ptr = objects[object_id].GetV(f_vector[f_id].v[2]);
  key i = v_ptr->get_index();
  delete v_ptr;
  return i;
}

void F::print(key object_id)
{
  cout.precision(12);
  cout << "Object " << objects[object_id].GetName() << endl
        << "Face " << get_index()
        << " " << get_v1_file(object_id)
        << " " << get_v2_file(object_id)
        << " " << get_v3_file(object_id)
        << endl;
}

E* Object::GetE(key id) // relative
{
  key e_id = object_vector[object_id].e_start + id; // absolute
  return new E(e_id);
}

V* Object::GetV(key id) // relative
{
  // DEBUG
  //cout << "111"; cout.flush();
  //cout << ", object_id="<<object_id<<", "; cout.flush();
  //ObjectData* pp = &object_vector[object_id];
  //cout << "222"; cout.flush();
  //cout << "f_start=" << object_vector[object_id].f_start<< ", ";
  //cout << "v_start=" << object_vector[object_id].v_start<< ", ";
  //cout.flush();
  key v_id = object_vector[object_id].v_start + id; // absolute
  //cout << "333"; cout.flush();
  // DEBUG
  return new V(v_id);
}

VE *Object::GetVE(key id) // relative
{
  key ve_id = object_vector[object_id].v_start + id; // absolute
  return new VE(ve_id);
}

void V::print(key object_id)
{
  cout.precision(12);
  cout << "Object " << objects[object_id].GetName() << endl
        << "Vertex " << get_index()
        << " " << get_x()
        << " " << get_y()
        << " " << get_z()
        << endl;
}

void V::printCP()
{
  cout.precision(12);
  cout <<  get_x() << " "
        << get_y() << " " 
        << get_z() << " 1 0 0 1"
        << endl;
}

void VE::print(key object_id)
{
  cout.precision(12);
  // for each adjacent face of vertex
  //for(int i=0;i<ADJACENT_FACES;i++)
  for(int i=1;i<ADJACENT_FACES;i++)
  {
    if(ve_vector[ve_id].adj_f[i]<0){break;}
    //F* f_ptr = objects[object_id].GetF(ve_vector[ve_id].adj_f[i]);
    objects[object_id].f.assign(ve_vector[ve_id].adj_f[i]);
    F* f_ptr = &objects[object_id].f;
    f_ptr->print(object_id);
    //delete f_ptr;
  }
}

/* *********************************************** */
/* *********************************************** */
/* *********************************************** */

void clearData(Container &c)
{
  c.clearData();
  e_vector.clear();
  f_vector.clear();
  fe_vector.clear();
  v_vector.clear();
  ve_vector.clear();
  objects.clear();
  object_vector.clear();

  border.clear();
  degen_f=disc_v=dupl_f=dupl_v=0;
  flipped.clear();
  found.clear();
  indistin_e.clear();
  indistin_v.clear();
  miss_v=0;
  nonm_e.clear();
  nonman_e.clear();
  nonman_v.clear();
  orphan.clear();
  seq_v=true;

  vertex_index=face_index=edge_index=0;
}

double getCoord(char* &triplet,char *cp)
{
  char val[80];
  char *eptr;
  int i;
  // get past white space
  while (strchr(" \t,",*triplet)!=NULL) { triplet++; }
  i=0;
  while (strchr("0123456789+-eE.",*triplet)!=NULL)
  {
    val[i++] = *triplet++;
  }
  val[i]=0;
  double a=strtod(val,&eptr);
  if (val==eptr)
  {
    printf("\nError in reading vertex: string %s\n",cp);
    exit(0);
  }
  return a;
}

void ReadVertex(char* triplet,Vdata &vd)
{
  char *cp=triplet;
  // get past 'Vertex'
  while (strchr("Vertx",*triplet)!=NULL) {triplet++;}
  // grab vertex index
  // TODO find conversion which supports long long int
  key v = static_cast<key>(getCoord(triplet,cp));
  if(v > std::numeric_limits<key>::max())
  {
    cout << "\nError: vertex index = "
          << v << " is larger than max allowable integer value = "
          << std::numeric_limits<key>::max() << endl;
    exit(0);
  }
  vd.index = v;
  // grab coords
  vd.x=getCoord(triplet,cp);
  vd.y=getCoord(triplet,cp);
  vd.z=getCoord(triplet,cp);
}


key getIndex(char* &triplet,char *line)
{
  char val[80];
  char *eptr;
  int i;
  // get past white space
  while (strchr(" \t,",*triplet)!=NULL) { triplet++; }
  i=0;
  // scan number
  while (strchr("0123456789+-eE.",*triplet)!=NULL)
  {
    val[i++] = *triplet++;
  }
  val[i]=0;
  // TODO find conversion which supports long long int
  key v = static_cast<key>(strtod(val,&eptr));
  // check for string length == 0
  if (val==eptr)
  {
    printf("Error in reading index from line: %s\n",line);
    exit(0);
  }
  // check for index > max int
  if(v > std::numeric_limits<key>::max())
  {
    cout << "\nError: face index = "
          << v << " is larger than max allowable integer value = "
          << std::numeric_limits<key>::max() << endl;
    exit(0);
  }
  return v;
}

key processIndex(char* &triplet,char* line,m_kk p,bool flag)
{
  // map<vertex index from file, relative vertex index>
  //stxxl::map<key,key,kSortComp>::const_iterator p;
  // grab first vertex index
  // search for matching vertex index in map
  // p=vp.find(getIndex(triplet,line));
  // if match found
  if(flag)
  {
    // look for vertex index in orphan set 
    // orphan set contains vertex index from file
    const std::set<key>::iterator q = orphan.find((*p).first);
    // if found
    if(q!=orphan.end())
    {
      // remove matching element from orphan set
      orphan.erase(q);
    }
    // set Face vertex key, first
    //fd.v[i] = (*p).second;
    return (*p).second;
  }
  else
  {
    // missing vertex!
    miss_v++;
    // TODO add -p option dependency
    // print details
    cout << "\nError: Face references missing vertex " 
          << (*p).first << endl
          << "line: " << line << endl;
    // null variable
    //fd.v[0] = 0;
    return 0;
  }
}


void ReadFace(char *triplet,Fdata &fd,stxxl::map<key,key,kSortComp> &vp)
{
  // TODO use strtod return values
  // 1) If no valid conversion could be performed, 
  // a zero value (0.0) is returned.
  // 2) If the correct value is out of 
  // the range of representable values, 
  // a positive or negative HUGE_VAL is returned, 
  // and the global variable errno is set to ERANGE.
  // 3) If the correct value would cause underflow, 
  // zero is returned and errno is set to ERANGE.

  char *line = triplet;
  // map<vertex index from file, relative vertex index>
  stxxl::map<key,key,kSortComp>::const_iterator p;

  // get past 'Face'
  while (strchr("Face",*triplet)!=NULL) {triplet++;}
  // grab Face index
  fd.index = getIndex(triplet,line);

  // grab first vertex index
  // search for matching vertex index in map
  p=vp.find(getIndex(triplet,line));
  fd.v[0] = processIndex(triplet,line,p,p!=vp.end());
  // grab second vertex index
  // search for matching vertex index in map
  p=vp.find(getIndex(triplet,line));
  fd.v[1] = processIndex(triplet,line,p,p!=vp.end());
  // grab third vertex index
  // search for matching vertex index in map
  p=vp.find(getIndex(triplet,line));
  fd.v[2] = processIndex(triplet,line,p,p!=vp.end());
  // check for vertex indices used in Face
  // require indices by nonzero
  if(fd.v[0]>0 && fd.v[1]>0 && fd.v[2]>0)
  {
    if(	fd.v[0]==fd.v[1] ||
        fd.v[1]==fd.v[2] ||
        fd.v[2]==fd.v[0] )
    {
      // duplicate!
      degen_f++;
      // TODO add -p option dependency
      // print details
      cout << "\nError: Face " << fd.index
            << " does NOT reference three unique vertex indices.\n" 
            << line;
    }
  }
}

void scanFile (key object_id,std::string filename)
{
  char line[2048];
  FILE *FF = fopen(filename.c_str(),"r");
  if(!FF)
  {
    cout <<"Couldn't open input file " << filename << endl;
    exit(0);
  }
  else
  {
    cout << "\n\n" << "/* ********************** "
          << "OBJECT ********************** */\n";
    cout << "name: " << filename << endl;
    cout.flush();
  }
  // for tracking
  cout << "Counting lines in input file .........................";
  cout.flush();
  int count=0;
  // for every line in file
  for(char *str=fgets(line,2048,FF);str!=NULL;str=fgets(line,2048,FF))
  {
    count++;
  }
  cout << count;
  cout << "...complete.\n";cout.flush();
  rewind(FF);
  cout << "Reading object into memory .........................";
  cout.flush();
  int iii=1;
  double goal = 0.2;
  printf("0%%..");
  fflush(stdout);

  // map<vertex index from file, relative vertex index>
  stxxl::map<key,key,kSortComp> vp(NODE_CACHE_SIZE,LEAF_CACHE_SIZE);

  // for every line in file
  for(char *str=fgets(line,2048,FF);str!=NULL;str=fgets(line,2048,FF))
  {
    //cout << line;cout.flush();
    // skip leading whitespace
    while (strchr(" \t,",*str)!=NULL) { str++;}
    // if first character is V for Vertex
    if (strchr("V",*str)!=NULL)
    {
      Vdata vd;
      ReadVertex(str,vd);
      v_vector.push_back(vd);
      // load orphan vector for later pruning
      orphan.insert(vd.index);
      // if index found in map
      // map<vertex index from file, relative vertex index>
      stxxl::map<key,key,kSortComp>::const_iterator p=vp.find(vd.index);
      if(p!=vp.end())
      {
        // duplicate!
        dupl_v++;
        // TODO add -p option dependency
        // print details
        Vdata vv = v_vector[(*p).second];
        cout << "\nError: "
              << "Two Vertices found with the same index ("
              << vd.index << ")\n";
        cout << "Current line: " << line;
        cout << "Matching line: Vertex "
              << vv.index << " "
              << vv.x << " "
              << vv.y << " "
              << vv.z << endl << endl;
      }
      else
      {
        // add vertex to map
        // map<vertex index from file, relative vertex index>
        vp.insert(std::make_pair(vd.index,vertex_index));
      }
      // check vertex index sequentiality
      vertex_index++; // relative 
      if(vertex_index!=vd.index) 
      {
        seq_v=false;
        disc_v = vd.index;
        // TODO add -p option dependency
        // print details
        cout << "\nWarning: Nonsequential vertex index:\n"
              << line;
      }
    }
    // if first character is F for Face
    else if (strchr("F",*str)!=NULL)
    {
      // require no duplicate vertex indices
      if(dupl_v>0){
        cout << "\nDuplicate vertex indices found. Skipping Face scan.\n";
        break;
      }
      Fdata fd;
      ReadFace(str,fd,vp);
      f_vector.push_back(fd);
      face_index++; // relative
      // look for face index in found set
      std::set<key>::iterator q=found.find(fd.index);
      // if match discovered
      if(q!=found.end())
      {
        // duplicate!
        dupl_f++;
        // TODO add -p option dependency
        // print details
        cout << "\nError: "
              << "Two Faces found with the same index ("
              << fd.index << ")\n"; 
        cout << "Current line: " << line << endl;
      }
      else
      {
        // add index to set
        found.insert(fd.index);
      }
    }
    // track progress
    double progress = static_cast<double>(iii++)/count;
    if(progress>goal){
      printf("%d%%..",static_cast<int>(goal*100));
      fflush(stdout);
      goal+=0.2;
    }
  }
  vp.clear();
  fclose(FF);
  printf("100%%..");
  fflush(stdout);
  cout << "complete.\n";cout.flush();
}

Object* Container::processFile(std::string filename)
{
  // create new object
  Object o(object_id,filename);
  //o.clearStats();
  objects.push_back(o);
  // initialize vector bounds
  ObjectData od = {0};
  od.v_start = getNextV(); // absolute
  od.f_start = getNextF(); // absolute
  object_vector.push_back(od);
  // vertex and face count in object
  vertex_index=face_index=0; // relative
  seq_v = true;
  // track degeneracies
  miss_v=dupl_v=dupl_f=degen_f=0;
  orphan.clear();
  found.clear();
  // scan file
  scanFile(object_id,filename);
  // check object contents
  if(vertex_index==0 || face_index==0) // relative
  { 
    // delete object
    cout << "\n Object::processFile: "
          << "no valid mesh object found in "
          << filename << ". Skipping file.\n";
    return NULL;
  }
  else
  {
    setNextV(od.v_start + vertex_index); // absolute+relative
    setNextF(od.f_start + face_index); // absolute+relative
    object_vector[object_id].v_end = getNextV()-1; // absolute
    object_vector[object_id].f_end = getNextF()-1; // absolute
    countObj();
    setNumV(vertex_index);
    setNumF(face_index);
    // return
    return &objects.back();
  }
}

void Container::scanDir(const char *filename)
{
  struct dirent *pent;			// pointer to dirent structure
  DIR *pdir = opendir(filename);	// pointer to a directory data structure
  if (!pdir) {printf("Error. Could not open %s.\n",filename);exit(1);}
  else { cout << "\nFolder found " << filename << endl << endl;}
  while((pent=readdir(pdir)))
  {
    // copy char array to string
    std::string str = pent->d_name;
    // if file of typ *.mesh
    std::string::size_type found = str.find(".mesh",0);
    // if found
    if(found != std::string::npos)
    {
      // save filename
      files.push_back(str);
      // update index
      num_files++;
    }
  }
  closedir(pdir);
}

void parse(int argc,char **argv)
{
  // if no arguments passed
  if(argc==1){
    cout << message << endl;
    exit(0);
  }

  int c;
  opterr=0;
  char *eptr=NULL;
  while((c=getopt(argc,argv,"ab:c:d:e:f:hipq")) != -1)
    switch(c)
    {
      case 'a':
        // determine attributes only
        attr=true;
        break;
      case 'b':
        // specify max edge length threshold (units of input file)
        my_signal[4]=true;
        thresholds[4]= strtod(optarg,&eptr);
        break;
      case 'c':
        // specify min edge length threshold (units of input file)
        my_signal[3]=true;
        thresholds[3]= strtod(optarg,&eptr);
        break;
      case 'd':
        // specify max edge angle threshold (degrees)
        my_signal[2]=true;
        thresholds[2]= strtod(optarg,&eptr);
        break;
      case 'e':
        // specify min edge angle threshold (degrees)
        my_signal[1]=true;
        thresholds[1]= strtod(optarg,&eptr);
        break;
      case 'f':
        // specify max aspect ratio threshold
        my_signal[0]=true;
        thresholds[0]= strtod(optarg,&eptr);
        break;
      case 'h':
        cout << message << endl;
        abort ();
      case 'i':
        interf=true;
        break;
      case 'p':
        my_print=true;
        sprintf(style,"detail");
        break;
      case 'q':
        my_print=true;
        sprintf(style,"cp");
        break;
      case '?':
        if (optopt == 'c')
          fprintf (stderr,"Option -%c requires an argument.\n", optopt);
        else if (isprint (optopt))
          fprintf (stderr, "Unknown option '-%c.\n", optopt);
        else
          fprintf (stderr,
                   "Unknown option character '\\x%x'.\n",
                   optopt);
        exit(0);
      default:
        cout << message << endl;
        abort ();
    }

  if(interf==true && folder==false){
    cout << "\n\nError. The '-i' option applies only to folders.\n\n";
    exit(0);
  }

  if(optind<argc)
  {
    for (int index=optind;index<argc;index++)
    {
      // determine if file or folder was passed
      struct stat stat_p;
      stat(argv[index],&stat_p);
      if(stat_p.st_mode & S_IFDIR){ folder=true;}
      else{ folder=false;}
      // adjust input directory
      char filename[FILENAME_SIZE];
      strcpy(filename,argv[index]);
      if(folder)
      {
        char *temp=strrchr(filename,'/');
        if(!temp) {strcat(filename,"/");}
        else if(*++temp) {strcat(filename,"/");}
      }
      if(inpath.empty()==true){ inpath=filename; }
      else if(outpath.empty()==true){ outpath=filename; }
    }
  }
  else
  {
    fprintf (stderr,"No input file argument found on command line.\n");
    exit(0);
  }
}

void PrintObject(key object_id)
{
  cout << endl << endl;
  int num=object_vector[object_id].v_end
        -object_vector[object_id].v_start+1;
  for(int i=0;i<num;i++)
  {
    V* v_ptr = objects[object_id].GetV(i);
    v_ptr->print(object_id);
    delete v_ptr;
  }
  num=object_vector[object_id].f_end
        -object_vector[object_id].f_start+1;
  for(int i=0;i<num;i++)
  {
    //F* f_ptr = objects[object_id].GetF(i);
    objects[object_id].f.assign(i);
    F* f_ptr = &objects[object_id].f;
    f_ptr->print(object_id);
    //delete f_ptr;
  }
  cout << endl << endl;
}

void PrintObjectSum(key object_id)
{
  int nv=object_vector[object_id].v_end
        -object_vector[object_id].v_start+1;
  int nf=object_vector[object_id].f_end
        -object_vector[object_id].f_start+1;
  int ne=object_vector[object_id].e_end
        -object_vector[object_id].e_start+1;
  cout << "# vertices = " << nv
        << ", # faces = " << nf
        << ", # edges = " << ne << endl;
}
/*
void addFace(VE* ve_ptr,key f)
{
  int next = ve_ptr->get_next();
  // if array is full
  if(next==ADJACENT_FACES)
  {
    cout << "\nError: increase ADJACENT_FACES (current value = "
          << ADJACENT_FACES << ").\n";
    cout << "face relative key = " << f << endl;
    exit(0);
  }
  else
  {
    // the next spot is empty
    ve_ptr->set_next(f,next);
  }
}*/

void addFace(VE* ve_ptr,key f)
{
  ve_ptr->addFace(f);
/*  int next = ve_ptr->get_next();
  // if array is full
  if(next==ADJACENT_FACES)
  {
    cout << "\nError: increase ADJACENT_FACES (current value = "
          << ADJACENT_FACES << ").\n";
    cout << "face relative key = " << f << endl;
    exit(0);
  }
  else
  {
    // the next spot is empty
    ve_ptr->set_next(f,next);
  }
  */
}

void Object::boundObject()
{
  //initialize bounds
  V* v_ptr = GetV(0);
  bb[0] = v_ptr->get_x();
  bb[1] = v_ptr->get_y();
  bb[2] = v_ptr->get_z();
  bb[3] = v_ptr->get_x();
  bb[4] = v_ptr->get_y();
  bb[5] = v_ptr->get_z();
  delete v_ptr;
  // for each vertex in object
  int num=object_vector[object_id].v_end
        -object_vector[object_id].v_start+1;
  for(int i=0;i<num;i++) // relative
  {
    V* v_ptr = objects[object_id].GetV(i);
    if(v_ptr->get_x()<bb[0]){bb[0]=v_ptr->get_x();}
    if(v_ptr->get_x()>bb[3]){bb[3]=v_ptr->get_x();}
    if(v_ptr->get_y()<bb[1]){bb[1]=v_ptr->get_y();}
    if(v_ptr->get_y()>bb[4]){bb[4]=v_ptr->get_y();}
    if(v_ptr->get_z()<bb[2]){bb[2]=v_ptr->get_z();}
    if(v_ptr->get_z()>bb[5]){bb[5]=v_ptr->get_z();}
    delete v_ptr;
  }
}
/*
void PrintEdges(Container &c)
{
  key oid = c.getObjectId();
  // for each edge in object
  int num=object_vector[oid].e_end
        -object_vector[oid].e_start+1;
  for(int i=0;i<num;i++) // relative
  {
    E* e_ptr = objects[oid].GetE(i);
    e_ptr->print(oid);
    delete e_ptr;
  }
}

void PrintAdjacencies(Container &c)
{
  key oid = c.getObjectId();
  // for each extra vertex data in object
  int num=object_vector[oid].ve_end
        -object_vector[oid].ve_start+1;
  for(int i=0;i<num;i++) // relative
  {
    V* v_ptr = objects[oid].GetV(i);
    v_ptr->print(oid);
    VE* ve_ptr = objects[oid].GetVE(i);
    ve_ptr->print(oid);
    delete ve_ptr;
  }
}*/

double getTime()
{
  //cout << "getTime(): begin\n"; 
  cout.flush();
  timeval tv;
  gettimeofday(&tv,NULL);
  //cout << "time = " << tv.tv_sec*1E6+tv.tv_usec << endl;
  return tv.tv_sec*1E6+tv.tv_usec;
  //return static_cast<int>(tv.tv_sec*1E6+tv.tv_usec);
}

void Object::findVertexAdjacencies(Container &c)
{
  //cout << "Finding vertex adjacencies .........................";
  //cout.flush();
  cout << "Finding vertex adjacencies .........................\n";
  cout.flush();
  ObjectData od = object_vector[object_id];
  od.ve_start = c.getNextVE(); // absolute
  cout << "Create vertex adjacency containers .........................\n";
  cout.flush();
  // vertices in object
  int num=object_vector[object_id].v_end
        -object_vector[object_id].v_start+1;
  // tracking
  double trackingInterval = 0.01;
  int goal = static_cast<int>(trackingInterval*num);
  printf("0%%..");
  fflush(stdout);
  // declare time variable
  double begintime = getTime();
  // for each vertex in object
  for(int i=0;i<num;i++) // relative
  {
    // create matching ve_vector entry
    VEdata ved;
    ved.adj_f[0]=1;
    for(int j=1;j<ADJACENT_FACES;j++)
    {
      ved.adj_f[j]=-2;
    }
    ve_vector.push_back(ved);
/**/    // DEBUG
    // track progress
    if(i>goal)
    {
      //double progress = static_cast<double>(i)/num;
      printf("%g%%(%g verts/sec)..\n",static_cast<double>(i)/static_cast<double>(num)*100.0,num*trackingInterval/(getTime()-begintime));
      fflush(stdout);
      //goal+=trackingInterval;
      //begintime = time(NULL);
      begintime = getTime();
      goal+=static_cast<int>(trackingInterval*num);
    }
    // DEBUG /**/
  }
  printf("100%%..complete.\n");
 // printf("100%%(%g verts/sec)..\n",num*trackingInterval/(getTime()-begintime));
  fflush(stdout);
  cout << "Fill vertex adjacency containers .........................\n";
  cout.flush();
  printf("0%%..");
  fflush(stdout);
  // faces in object
  num=object_vector[object_id].f_end
        -object_vector[object_id].f_start+1;
  //tracking
  //goal = trackingInterval;
  goal = static_cast<int>(trackingInterval*num);
  // declare time variable
  begintime = getTime();
  double starttime = begintime;
  // for each face in object
  for(int i=0;i<num;i++) // relative
  {
    //cout << i << endl;
    //F* f_ptr = GetF(i);
    f.assign(i);
    F* f_ptr = &f;
    key th = f_ptr->get_rel_key(object_id);
    //VE* ve_ptr = GetVE(f_ptr->get_v1_rel());
    ve.assign(f_ptr->get_v1_rel());
    //addFace(ve_ptr,th);
    addFace(&ve,th);
    //delete ve_ptr;
    //ve_ptr = GetVE(f_ptr->get_v2_rel());
    ve.assign(f_ptr->get_v2_rel());
    //addFace(ve_ptr,th);
    addFace(&ve,th);
    //delete ve_ptr;
    //ve_ptr = GetVE(f_ptr->get_v3_rel());
    ve.assign(f_ptr->get_v3_rel());
    //addFace(ve_ptr,th);
    addFace(&ve,th);
    //delete f_ptr;
    //delete ve_ptr;
/**/    // DEBUG
    // track progress
    //double progress = static_cast<double>(i)/static_cast<double>(num);
    if(i>goal)
    {
      printf("%g%%(%g faces/sec)..\n",static_cast<double>(i)/static_cast<double>(num)*100.0,num*trackingInterval/(getTime()-begintime));
      //printf("%d%%(%g faces/sec)..\n",static_cast<int>(goal*100),num*0.05/(time(NULL)-begintime));
      fflush(stdout);
      //goal+=trackingInterval;
      begintime = getTime();
      //begintime = time(NULL);
      goal+=static_cast<int>(trackingInterval*num);
    }
    // DEBUG /**/
  }
  //printf("100%%(%g faces/sec)..\n",num*trackingInterval/(getTime()-begintime));
  printf("100%%..complete.\n");
  fflush(stdout);
  c.setNextVE(od.ve_start + face_index); // absolute+relative=absolute
  object_vector[object_id].ve_end = c.getNextV()-1; // relative
  //cout << "complete.\n";cout.flush();
  printf("...complete.\n");
  fflush(stdout);
  cout << "elapsed face processing time (usec): " << getTime()-starttime << endl;
}

bool checkIntegrity(void)
{
  cout << "Checking object integrity .........................";
  cout.flush();
  bool good=true;
  // if duplicate vertex indices found
  if(dupl_v>0)
  {
    // print number of duplicate vertices
    cout << "\n    Error: " << dupl_v
          << " pair(s) of duplicate vertex indices were found.\n";
    good=false;
  }
  // if duplicate face indices found
  if(dupl_f>0)
  {
    // print number of duplicate face indices
    cout << "\n    Error: " << dupl_f
          << " pair(s) of duplicate face indices were found.\n";
    good=false;
  }
  // if missing vertices encountered
  if(miss_v>0)
  {
    // print number of missing vertices
    cout << "\n    Error: " << miss_v
          << " missing vertices were found.\n";
    good=false;
  }
  // if vertex indexing is not contiguous
  if(seq_v==false)
  {
    // not an error, simply a warning
    cout << "Warning: First nonsequential vertex index found = "
          << disc_v << endl;
  }
  // if any face did not contain three unique vertex indices
  if(degen_f>0)
  {
    // print number of such faces
    cout << "\n    Error: " << degen_f
          << " faces did not contain three unique vertex indices.\n";
    good=false;
  }
  // check for unused vertices, orphan data from file
  if(orphan.empty()==false)
  {
    cout << "\n    Error: The following " << orphan.size()
          << " vertex indices were not referenced by any face.\n";
    for(std::set<key>::iterator i=orphan.begin();i!=orphan.end();i++)
    {
      cout << (*i) << ", ";
    }
    cout << endl;
    good=false;
  }
  // TODO check for contiguous face indexing, warning only
  cout << "complete.\n";cout.flush();
  return good;
}

std::pair<key,key> getPair(key va,key vb)
{
  // create pair of vertices
  // introduce bias [smallest first]
  if(va<vb)
  {
    return std::make_pair(va,vb);
  }
  else
  {
    return std::make_pair(vb,va);
  }
}

key look4Edge(Container &c,key va,key vb,map_edge &hm)
{
  // map<pair<vertex key,vertex key>,edge key>
  //           from file, from file ,relative
  //
  // if element exists given key, then get Edge pointer
  m_e i = hm.find(getPair(va,vb));
  if (i!=hm.end())
  {
    return (*i).second; // relative
  }
  return -2;
}

void printFaceEdges(Container &c)
{
  key oid = c.getObjectId();
  // for each extra face data in object
  int num=object_vector[oid].fe_end
        -object_vector[oid].fe_start+1;
  for(int i=0;i<num;i++) // relative
  {
    cout << endl << endl;
    //F* f_ptr = objects[oid].GetF(i);
    objects[oid].f.assign(i);
    F* f_ptr = &objects[oid].f;
    f_ptr->print(oid);
    //delete f_ptr;
    FE* fe_ptr = objects[oid].GetFE(i);
    key triplet[3] = {fe_ptr->get_edge(0),
                      fe_ptr->get_edge(1),
                      fe_ptr->get_edge(2)};
    delete fe_ptr;
    E *e_ptr = objects[oid].GetE(triplet[0]);
    e_ptr->print(oid);
    delete e_ptr;
    e_ptr = objects[oid].GetE(triplet[1]);
    e_ptr->print(oid);
    delete e_ptr;
    e_ptr = objects[oid].GetE(triplet[2]);
    e_ptr->print(oid);
    delete e_ptr;
  }
}

void halfEdge(Container &c,key f,key va,key vb,map_edge &hm)
{
  // map<pair<vertex key,vertex key>,edge key>
  //           from file, from file ,relative
  //
  // create new Edata
  Edata ed;
  ed.v1=va; // relative
  ed.v2=vb; // relative
  ed.f1=f;
  ed.f2=-2;
  // store edge pointer in hash table
  hm[getPair(ed.v1,ed.v2)]=edge_index; // relative
  // add edge pointer to face
  FE* fe_ptr = objects[c.getObjectId()].GetFE(f);
  fe_ptr->addEdge(edge_index); // relative
  delete fe_ptr;
  // save new Edata
  e_vector.push_back(ed);
  edge_index++;
}

void buildEdge(Container &c,key f,key va,key vb,map_edge &hm)
{
  // va and vb passed as v_vector relative indices
  key e = look4Edge(c,va,vb,hm);
  if(e<0)
  { 
    halfEdge(c,f,va,vb,hm); 
  }
  else
  {
    E *e_ptr = objects[c.getObjectId()].GetE(e);
    e_ptr->updateEdge(f,c.getObjectId());
    delete e_ptr;
  }
}

void createEdges(Container &c)
{
  key oid = c.getObjectId();
  cout << "Create edges for object ["
        << objects[oid].GetName()
        << "]...............................";
  cout.flush();
  // prep FEdata
  c.prepFEdata(oid);
  // prepare edge data
  object_vector[oid].e_start = c.getNextE();
  edge_index=0;
  // create map for finding edges
  // map<pair<vertex key,vertex key>,edge key>
  //           for v_vector, for v_vector ,for e_vector; ALL relative 
  stxxl::map<std::pair<key,key>,key,pSort> hm(NODE_CACHE_SIZE,LEAF_CACHE_SIZE);
  // faces in object
  key num=object_vector[oid].f_end-object_vector[oid].f_start+1;
  // tracking
  double trackingInterval = 0.2;
  int goal = static_cast<int>(trackingInterval*num);
  printf("0%%..");
  fflush(stdout);
  // for each face in object
  for(key i=0;i<num;i++)
  {
    //F* f_ptr = objects[oid].GetF(i);
    objects[oid].f.assign(i);
    F* f_ptr = &objects[oid].f;
    // TODO pass relative index as key
    buildEdge(c,f_ptr->get_rel_key(oid),f_ptr->get_v1_rel(),f_ptr->get_v2_rel(),hm);
    buildEdge(c,f_ptr->get_rel_key(oid),f_ptr->get_v2_rel(),f_ptr->get_v3_rel(),hm);
    buildEdge(c,f_ptr->get_rel_key(oid),f_ptr->get_v3_rel(),f_ptr->get_v1_rel(),hm);
    //delete f_ptr;
    // track progress
    if(i>goal)
    {
      printf("%g%%..",static_cast<double>(i)/static_cast<double>(num)*100.0);
      fflush(stdout);
      goal+=static_cast<int>(trackingInterval*num);
    }
  }
  printf("100%%..complete.\n");
  fflush(stdout);
  c.setNumE(edge_index);
  c.setNextE(object_vector[oid].e_start + edge_index);
  object_vector[oid].e_end = c.getNextE()-1;
  //cout << "complete.\n";cout.flush();
}

bool Object::isClosed(void)
{
  bool flag=true;
  // for each edge in object
  int num=object_vector[object_id].e_end
        -object_vector[object_id].e_start+1;
  for(int i=0;i<num;i++) // relative
  {
    E* e_ptr = objects[object_id].GetE(i);
    // if edge has no adjacent face f2
    if(e_ptr->get_f2(e_ptr->get_v1())<0)
    {
      flag=false;
      // record offending edge
      border.push_back(e_ptr->get_rel_key(object_id));
    }
  }
  return flag;
}

void V::getAdjacentEdges(key object_id,std::vector<key> &e)
{
  //double firsttime = getTime();
  key vid = get_rel_key(object_id);
  // gather adjacent faces of vertex
  Object* o_ptr = &objects[object_id];
  //objects[object_id].ve.assign(vid);
  o_ptr->ve.assign(vid);
  std::vector<key> af;
  //objects[object_id].ve.getAf(af);
  o_ptr->ve.getAf(af);
  //double secondtime = getTime();
  //double preamble=0,secondloop=0;
  //double ifclause=0,follow=0;
  //double inner1=0,inner2=0,inner3=0,inner4=0;
  // for each adjacent face of vertex
  for (std::vector<key>::iterator i=af.begin();i!=af.end();i++)
  {
    //double p1 = getTime();
    // gather edges of face
    //FE *fe_ptr = objects[object_id].GetFE(*i);
    FE *fe_ptr = o_ptr->GetFE(*i);
    // ep[i] == edge index, relative
    key triplet[3] = {fe_ptr->get_edge(0),
                      fe_ptr->get_edge(1),
                      fe_ptr->get_edge(2)};
    delete fe_ptr;
    //double p2 = getTime();
    // for each edge of face
    for(int j=0;j<3;j++)
    {
      //double pa = getTime();
      // if edge is undefined
      if(triplet[j]<0)
      {
        cout << "\nV::getAdjacentEdges: Error.\n"
              << "Face edge is missing.\n";
        //objects[object_id].f.assign(*i);
        o_ptr->f.assign(*i);
        //F *f_ptr = &objects[object_id].f;
        F *f_ptr = &o_ptr->f;
        f_ptr->print(object_id);
        for(int k=0;k<3;k++)
        {
          if((triplet[k]<0)==false)
          {
            //E *e_ptr = objects[object_id].GetE(triplet[k]);
            E *e_ptr = o_ptr->GetE(triplet[k]);
            cout << "\nep[" << k << "]:" << endl; 
            e_ptr->print(object_id);
            delete e_ptr;
          }
        }
        exit(0);
      }
      //double pa = getTime();
      // edge is ok
      //E *e_ptr = objects[object_id].GetE(triplet[j]);
      E *e_ptr = o_ptr->GetE(triplet[j]);
      // if edge contains vertex
      //double pb = getTime();
      //double pc = pb,pd = pb;
      if(e_ptr->get_v1()==vid || e_ptr->get_v2()==vid)
      {
        //pc = getTime();
        // add edge to vector
        e.push_back(triplet[j]);
        //pd = getTime();
      }
      //double pe = getTime();
      delete e_ptr;
      //double pf = getTime();
      //ifclause+=pb-pa;
      //follow+=pc-pb;
      //inner1+=pb-pa;
      //inner2+=pe-pc;
      //inner3+=pd-pc;
      //inner4+=pf-pe;
    }
    //double p3 = getTime();
    //preamble+=p2-p1;
    //secondloop+=p3-p2;
  }
  //double thirdtime = getTime();
  // keep unique edges
  sort(e.begin(),e.end());
  std::vector<key>::iterator new_end = unique(e.begin(),e.end());
  e.assign(e.begin(),new_end);
  //double fourthtime = getTime();
/*  cout << "V::getAdjacentEdges: "
       <<" delt1 = " << secondtime-firsttime
       <<", firstloop = " << thirdtime-secondtime
       <<" (= " << preamble << " + " << secondloop << ")"
       <<", secondloop = " << secondloop
       <<" (= " << ifclause << " + " << follow << ")"
       <<", delt3 = " << fourthtime-thirdtime << endl;
       */
  //cout << "V::getAdjacentEdges: "
  //     <<" pre_if = " << inner1
  //     <<", if = " << inner2
  //     <<", inside_if = " << inner3
  //     <<", delete = " << inner4 << endl;
  //cout.flush();
}

key F::getNewEdge(key object_id,key old,key vv)
{
  //double firsttime = getTime();
  Object* o_ptr = &objects[object_id];
  // return face edge that contains vertex vv
  // and is different from old edge
  // 
  // find old edge in Face
  bool efound=false;
  FE *fe_ptr = o_ptr->GetFE(get_rel_key(object_id));
  //FE *fe_ptr = objects[object_id].GetFE(get_rel_key(object_id));
  key triplet[3] = {fe_ptr->get_edge(0),
                    fe_ptr->get_edge(1),
                    fe_ptr->get_edge(2)};
  delete fe_ptr;
  //double secondtime = getTime();
  //cout << "\nF::getNewEdge: " << "delt1 = " << secondtime-firsttime << endl; cout.flush();
  // for each edge of face
  for(int i=0;i<3;i++)
  {
    // if edge is undefined
    if(triplet[i]<0)
    {
      cout << "\nF::getNewEdge: Error.\n"
            << "Face edge is missing.\n";
      exit(0);
    }
    if(triplet[i]==old) {efound=true;}
  }
  //double thirdtime = getTime();
  //cout << "F::getNewEdge: " << "delt2 = " << thirdtime-secondtime << endl; cout.flush();
  // if old edge is not found in face
  if(efound==false)
  {
    cout << "\n\nOld edge does not match any on face.\n";
    cout << "Current face:\n";
    print(object_id);
    cout << endl;
    cout << "Old edge:\n";
    E *e_ptr = o_ptr->GetE(old);
    //E *e_ptr = objects[object_id].GetE(old);
    e_ptr->print(object_id);
    delete e_ptr;
    cout << endl;
    exit(0);
  }
  //double fourthtime = getTime();
  //cout << "F::getNewEdge: " << "delt3 = " << fourthtime-thirdtime << endl; cout.flush();
  // if current vertex not found in face
  if(get_v1_rel()!=vv && get_v2_rel()!=vv && get_v3_rel()!=vv)
  {
    cout << "\n\nCurrent vertex does not match any on face.\n";
    cout << "Current face:\n";
    print(object_id);
    cout << endl;
    cout << "Current vertex:\n";
    V *v_ptr = o_ptr->GetV(vv);
    //V *v_ptr = objects[object_id].GetV(vv);
    v_ptr->print(object_id);
    delete v_ptr;
    cout << endl << endl;
    exit(0);
  }
  //double fifthtime = getTime();
  //cout << "F::getNewEdge: " << "delt4 = " << fifthtime-fourthtime << endl; cout.flush();
  for(int i=0;i<3;i++)
  {
    if(triplet[i]!=old)
    {
      E *e_ptr = o_ptr->GetE(triplet[i]);
      //E *e_ptr = objects[object_id].GetE(triplet[i]);
      key a = e_ptr->get_v1();
      key b = e_ptr->get_v2();
      delete e_ptr;
      if(a==vv || b==vv)
      {
        //cout << "F::getNewEdge: " << "delt5 = " 
        //      << getTime()-fifthtime << endl; cout.flush();
        return triplet[i];
      }
    }
  }
  cout << "\n\nFace::getNewEdge: Error. No edge on face contains current vertex \n"
        << "and is different from old vertex.\n\n";
  cout << "current vertex:\n";
  V *v_ptr = o_ptr->GetV(vv);
  //V *v_ptr = objects[object_id].GetV(vv);
  v_ptr->print(object_id);
  delete v_ptr;
  cout << endl << endl;
  cout << "old edge:\n";
  E *e_ptr = o_ptr->GetE(old);
  //E *e_ptr = objects[object_id].GetE(old);
  e_ptr->print(object_id);
  delete e_ptr;
  cout << endl << endl;
  cout << "current face:\n";
  print(object_id);
  cout << endl << endl;
  exit(0);
}

key E::getNewFace(key object_id,key old_face,key v)
{
  // return the adjacent face that is
  // different from old adjacent face
  if(old_face!=get_f1(v) && old_face!=get_f2(v)){
    cout << "Error. Neither edge adjacent face matches old face.\n";
    exit(0);
  }
  if(old_face==get_f1(v)){return get_f2(v);}
  else {return get_f1(v);}
}

bool V::scanAdjFaces(key object_id,key se,key sf,bool &nonman)
{
  //double starttime = getTime();
  Object* o_ptr = &objects[object_id];
  // collect touched edges
  std::vector<key> te;
  te.push_back(se);
  // collect touched faces
  std::vector<key> tf;
  tf.push_back(sf);
  // get new edge
  o_ptr->f.assign(sf);
  F* sf_ptr = &o_ptr->f;
  //double zerothtime = getTime();
  se=sf_ptr->getNewEdge(object_id,se,get_rel_key(object_id));
  //double firsttime = getTime();
  E* se_ptr = o_ptr->GetE(se);
  bool man = se_ptr->isManifold(object_id);
  delete se_ptr;
  // if new edge is not manifold
  if(man==false)
  {
    // vertex manifoldness will not be determined
    // because it's complicated, so just return flag
    nonman=true;
    return false;
  }
  //double secondtime = getTime();
  //double thirdtime = secondtime, fourthtime=secondtime, fifthtime=secondtime;
  //double sixthtime=secondtime, seventhtime=secondtime, eighthtime=secondtime;
  //double ninethtime=secondtime, tenthtime=secondtime;
  //double inner1=0, inner2=0, inner3=0, inner4=0, inner5=0, inner6=0, inner7=0;
  // while new edge has not already been touched
  while(find(te.begin(),te.end(),se)==te.end())
  {
    //thirdtime = getTime();
    // keep new edge
    te.push_back(se);
    // get new face
    se_ptr = o_ptr->GetE(se);
    //fourthtime = getTime();
    sf=se_ptr->getNewFace(object_id,sf,get_rel_key(object_id));
    //fifthtime = getTime();
    delete se_ptr;
    if (sf<0) {break;}
    // keep new face
    tf.push_back(sf);
    // get new edge
    o_ptr->f.assign(sf);
    sf_ptr = &o_ptr->f;
    //sixthtime = getTime();
    se=sf_ptr->getNewEdge(object_id,se,get_rel_key(object_id));
    //seventhtime = getTime();
    // if new edge is not manifold
    se_ptr = o_ptr->GetE(se);
    //eighthtime = getTime();
    bool man = se_ptr->isManifold(object_id);
    //ninethtime = getTime();
    delete se_ptr;
    if(man==false)
    {
      // vertex manifoldness cannot be determined
      // because it's complicated, so just return flag
      nonman=true;
      return false;
    }
    //tenthtime = getTime();
    //inner1+=fourthtime-thirdtime;
    //inner2+=fifthtime-fourthtime;
    //inner3+=sixthtime-fifthtime;
    //inner4+=seventhtime-sixthtime;
    //inner5+=eighthtime-seventhtime;
    //inner6+=ninethtime-eighthtime;
    //inner7+=tenthtime-ninethtime;
  }
  //double eleventhtime = getTime();
  // if number of touched faces != number of vertex adjacent faces
  o_ptr->ve.assign(get_rel_key(object_id));
  bool a = (tf.size()==o_ptr->ve.size());
  /*cout << "\nV::scanAdjFaces: "
        << "delt0 = " << zerothtime-starttime
        << "delt1 = " << firsttime-zerothtime
        << "delt2 = " << secondtime-firsttime
        << ", whileouter = " << eleventhtime-secondtime
        << ", inner1 = " << inner1
        << ", inner2 = " << inner2
        << ", inner3 = " << inner3
        << ", inner4 = " << inner4
        << ", inner5 = " << inner5
        << ", inner6 = " << inner6
        << ", inner7 = " << inner7
        << endl;
        */
  return a;
}

bool V::isManifold(key object_id,bool flag)
{
  //cout << "\nV:isManifold: starting\n"; cout.flush();
  //print(object_id);
  //printAdjacent(object_id);
  //cout << "\n\n"; cout.flush();
  //ObjectData* aa = &object_vector[object_id];
  Object* o_ptr = &objects[object_id];
  // try to "walk" around vertex adjacent faces
  // using adjacent face edges
  // if an edge is nonmanifold then abort mission
  // stop when the walk returns to starting adjacent face

  std::vector<key> e; // relative edge indices
  //double firsttime = getTime();
  getAdjacentEdges(object_id,e);
  //double secondtime = getTime();
  //cout << "V::isManifold: delt = " << secondtime-firsttime << endl;
  //cout.flush();
  // DEBUG
  //for(std::vector<key>::iterator i=e.begin();i!=e.end();i++)
  //{
  //  E* e_ptr = objects[object_id].GetE(*i);
  //  e_ptr->print(object_id);
  //  delete e_ptr;
  //}
  // DEBUG
  //cout << "bbb"; cout.flush();
  //ObjectData* bb = &object_vector[object_id];
  //cout << "ccc"; cout.flush();
  if(e.empty()==true){
    cout << "\n\nVertex::isManifold: Error."
          << " Vertex was not an 'orphan', but has no edges.\n";
    cout << "\n\nVertex::isManifold: Confused vertex:\n";
    print(object_id);
    cout << endl;
    exit(0);
  }
  // grab starting edge
  //E* se = objects[object_id].GetE(e.front());
  E* se = o_ptr->GetE(e.front());
  key cw=-2,ccw=-2;
  // if edge is manifold
  //double thirdtime = getTime();
  if(se->isManifold(object_id)==true)
  {
    // DEBUG
    //cout << "\n\nstarting edge is manifold:\n";cout.flush();
    //se->print(object_id);
    //cout << "\n";cout.flush();
    // DEBUG
    // grab starting face from edge
    // cw face should be clockwise from starting edge
    // relative to current vertex
    //
    // analyze f1
    //F *f_ptr = objects[object_id].GetF(se->get_f1());
    //objects[object_id].f.assign(se->get_f1());
    cw = se->get_f1(get_rel_key(object_id));
    ccw = se->get_f2(get_rel_key(object_id));
    // DEBUG
    //cout << "\ncw face key = " << cw
    //      << ", ccw face key = " << ccw << endl;
    //cout.flush();
    //o_ptr->f.assign(se->get_f1(get_rel_key(object_id)));
    //F *f_ptr = &o_ptr->f;
    //f_ptr->print(object_id);
    //cout << "\n";cout.flush();
    // DEBUG
    /*
    //F *f_ptr = &objects[object_id].f;
    key triplet[3] = {f_ptr->get_v1_rel(),
                      f_ptr->get_v2_rel(),
                      f_ptr->get_v3_rel()};
    //delete f_ptr;
    key e1 = se->get_v1();
    key e2 = se->get_v2();
    key a=-2,b=-2;
    if(triplet[0]==e1) a=0;
    else if(triplet[1]==e1) a=1;
    else if(triplet[2]==e1) a=2;
    else {cout << "Error. What the 1?\n";exit(0);}
    if(triplet[0]==e2) b=0;
    else if(triplet[1]==e2) b=1;
    else if(triplet[2]==e2) b=2;
    else {cout << "Error. What the 2?\n";exit(0);}
    cout << "a = " << a << ", b = " << b << endl << endl;
    if(
       (a==0 && b==1) ||
       (a==1 && b==2) ||
       (a==2 && b==1))
    {
      if(triplet[a]==get_rel_key(object_id))
      {
        ccw = se->get_f1(get_rel_key(object_id));
        cw = se->get_f2(get_rel_key(object_id));
      }
      else if(triplet[b]==get_rel_key(object_id))
      {
        cw = se->get_f1(get_rel_key(object_id));
        ccw = se->get_f2(get_rel_key(object_id));
      }
      else
      {
        cout << "\nV::isManifold: neither triplet[a]"
              << " nor triplet[b] match current vertex.\n";
        exit(0);
      }
    }
    else if(
       (b==0 && a==1) ||
       (b==1 && a==2) ||
       (b==2 && a==1))
    {
      if(triplet[b]==get_rel_key(object_id))
      {
        ccw = se->get_f1(get_rel_key(object_id));
        cw = se->get_f2(get_rel_key(object_id));
      }
      else if(triplet[a]==get_rel_key(object_id))
      {
        cw = se->get_f1(get_rel_key(object_id));
        ccw = se->get_f2(get_rel_key(object_id));
      }
      else
      {
        cout << "\nV::isManifold: neither triplet[a]"
              << " nor triplet[b] match current vertex.\n";
        exit(0);
      }
    }
    else
    {
      cout << "\nSomething weird happened.\n";exit(0);
    }
    */
  }
  else
  {
    // vertex manifoldness cannot be determined
    // because the starting edge is not manifold
    // thus determining manifoldness of vertex
    // is complicated, so just return flag
    return flag;
    // TODO: report number of vertices for which manifoldness
    // was not determined and alert user
  }
  //ObjectData* cc = &object_vector[object_id];
  //cout << "ddd"; cout.flush();

  // try clockwise
  //objects[object_id].f.assign(cw); // cw=relative
  o_ptr->f.assign(cw); // cw=relative
  //F *sf=&objects[object_id].f; // cw=relative
  F *sf=&o_ptr->f; // cw=relative
  bool nonman=false;
  //cout << "eee"; cout.flush();
  //ObjectData* dd = &object_vector[object_id];
  //cout << "fff"; cout.flush();
  // if fail i.e. not all adjacent faces touched
  //double fourthtime = getTime();
  //double fifthtime = fourthtime, sixthtime = fourthtime, seventhtime = fourthtime;
  //double eighthtime = fourthtime, ninethtime = fourthtime, tenthtime = fourthtime;
  if(scanAdjFaces(object_id,se->get_rel_key(object_id),sf->get_rel_key(object_id),nonman)==false){
    //cout << "ggg"; cout.flush();
    //ObjectData* ee = &object_vector[object_id];
    //cout << "hhh"; cout.flush();
    // if nonmanifold edge found, then bail
    //fifthtime = getTime();
    if(nonman==true){return flag;}
    // if ccw face was set
    if((ccw<0)==false)
    {
      //cout << "iii"; cout.flush();
      //ObjectData* ff = &object_vector[object_id];
      //cout << "jjj"; cout.flush();
      // try counter-clockwise
      //objects[object_id].f.assign(ccw); // ccw=relative
      o_ptr->f.assign(ccw); // ccw=relative
      sf=&o_ptr->f; // ccw=relative
      //sf=&objects[object_id].f; // ccw=relative
      // if fail i.e. all adjacent faces still not touched
      //cout << "kkk"; cout.flush();
      //ObjectData* gg = &object_vector[object_id];
      //cout << "lll"; cout.flush();
      //sixthtime = getTime();
      if(scanAdjFaces(object_id,se->get_rel_key(object_id),sf->get_rel_key(object_id),nonman)==false)
      {
        //cout << "\n\nvertex is NOT manifold. Are you sure?\n\n";cout.flush();exit(0);
        //cout << "mmm"; cout.flush();
        //ObjectData* hh = &object_vector[object_id];
        //cout << "nnn"; cout.flush();
        // record offending vertex
        //seventhtime = getTime();
        nonman_v.push_back(get_rel_key(object_id));
        flag=false;
        //eighthtime = getTime();
        //cout << "ooo"; cout.flush();
        //ObjectData* ii = &object_vector[object_id];
        //cout << "ppp"; cout.flush();
      }
      //ninethtime = getTime();
    }
    else
    {
      // record offending vertex
      nonman_v.push_back((get_rel_key(object_id)));
      flag=false;
    }
    //tenthtime = getTime();
  } // else all adjacent faces touched
  //double eleventhtime = getTime();
  //cout << "qqq"; cout.flush();
  //ObjectData* jj = &object_vector[object_id];
  //cout << "rrr\n"; cout.flush();
  //cout << "V:isManifold: finished\n"; cout.flush();
  //delete sf;
  delete se;
  /*cout << "V::isManifold: "
        << "delt1 = " << secondtime-firsttime
        << ", delt2 = " << thirdtime-secondtime
        << ", delt3 = " << fourthtime-thirdtime
        << ", outerif = " << eleventhtime-fourthtime
        << ", outerifbody = " << tenthtime-fifthtime
        << ", innerif = " << ninethtime-sixthtime
        << ", innerifbody = " << eighthtime-seventhtime
        << endl;
  cout.flush();
  */
  return flag;
}

bool Object::verticesManifold(bool flag)
{
  cout << "\nObject::verticesManifold starting\n";
  cout.flush();
  ///// if vertices are manifold /////
  ////// confirm that all adjacent faces /////
  ///// are consecutively reachable by edge hops /////

  // for each vertex in object
  //int num=object_vector[object_id].v_end
  //      -object_vector[object_id].v_start+1;
  int num=cov[object_id].v_end
        -cov[object_id].v_start+1;
  for(int i=0;i<num;i++) // relative
  {
    //double firsttime = getTime();
    //V* v_ptr = objects[object_id].GetV(i);
    //cout << "aaa"; cout.flush();
    Object* o_ptr = &objects[object_id];
    //cout << "zzz"; cout.flush();
    //cout << "object name="<<o_ptr->GetName()<<endl;cout.flush();
    V* v_ptr = o_ptr->GetV(i);
    //cout << "bbb"; cout.flush();
    // if vertex is not an orphan
    sk j = orphan.find(v_ptr->get_rel_key(object_id));
    //ObjectData* vv = &object_vector[object_id];
    //cout << "ccc"; cout.flush();
    //double secondtime = getTime();
    //double thirdtime=secondtime,fourthtime=secondtime;
    if(j==orphan.end())
    {
      //double thirdtime = getTime();
      //ObjectData* zz = &object_vector[object_id];
      //cout << "ddd"; cout.flush();
      flag = v_ptr->isManifold(object_id,flag);
      //cout << "DDD"; cout.flush();
      //ObjectData* tt = &object_vector[object_id];
      //double fourthtime = getTime();
      //cout << "Object::verticesManifold: vertex " << i << " of " << num
      //      <<", delt = " << fourthtime-thirdtime << endl << endl << endl;
      //cout.flush();
    }
    //cout << "eee"; cout.flush();
    //ObjectData* qq = &object_vector[object_id];
    delete v_ptr;
    //cout << "fff"; cout.flush();
    //ObjectData* pp = &object_vector[object_id];
    //cout << "ggg\n"; cout.flush();
  }
  cout << "\nObject::verticesManifold finished\n";
  cout.flush();
  return flag;
}

bool Object::edgesManifold(bool flag)
{
  cout << "\nObject::edgesManifold starting\n";
  cout.flush();
  // for each edge in object
  //int num=object_vector[object_id].e_end
  //      -object_vector[object_id].e_start+1;
  int num=cov[object_id].e_end
        -cov[object_id].e_start+1;
  for(int i=0;i<num;i++) // relative
  {
    E* e_ptr = objects[object_id].GetE(i);
    if(e_ptr->isManifold(object_id)==false)
    {
      flag=false;
      // record offending edge
      nonman_e.push_back(e_ptr->get_rel_key(object_id));
    }
    delete e_ptr;
  }
  cout << "\nObject::edgesManifold finished\n";
  cout.flush();
  return flag;
}

bool Object::isManifold(void)
{
  bool flag=true;
  return verticesManifold(edgesManifold(flag));
}

bool E::getOrientation(key fid,key object_id)
{
  // gather f1 vertices
  //F *f_ptr = objects[object_id].GetF(fid);
  objects[object_id].f.assign(fid);
  F *f_ptr = &objects[object_id].f;
  key f[3] = {f_ptr->get_v1_rel(),
    f_ptr->get_v2_rel(),
    f_ptr->get_v3_rel()};
  //delete f_ptr;
  key v1 = get_v1();
  key v2 = get_v2();
  if((f[0]==v1 && f[1]==v2) ||
     (f[1]==v1 && f[2]==v2) ||
     (f[2]==v1 && f[0]==v2))
  {
    return true;
  }
  return false;
}

bool E::isConsistent(key object_id)
{
  // border edge
  if(get_f2(get_v1())<0) { return true; }
  // not border
  bool a = getOrientation(get_f1(get_v1()),object_id);
  bool b = getOrientation(get_f2(get_v1()),object_id);
  if(a==b){ return false; }
  else {return true;}
}

bool Object::isConsistent(void)
{
  // assuming object mesh is manifold...
  bool flag=true;
  // for each edge in object
  int num=object_vector[object_id].e_end
        -object_vector[object_id].e_start+1;
  for(int i=0;i<num;i++) // relative
  {
    // if the one or two adjacent edge faces traverse the edge
    // in the same direction, either both v1->v2 or both v2->v1,
    // then the entire set of mesh faces is declared inconsistent
    // i.e. the first edge with a flipped face will trigger signal return
    // and no more edge checking is performed, but -p option
    // needs ALL flipped edges.
    E* e_ptr = objects[object_id].GetE(i);
    if(e_ptr->isConsistent(object_id)==false)
    {
      flag=false;
      // record offending edge
      flipped.push_back(e_ptr->get_rel_key(object_id));
    }
    delete e_ptr;
  }
  return flag;
}
/* TODO port this old code
bool Object::isOutward(Space &s){
  // assuming object mesh is closed, manifold, and consistent...
  //
  bool line_flag=false, poly_flag=false, poly_edge_flag=false;
  f_iterator fff=f.begin();
  double count;
  do{
    // set origin face
    Face *ff = *fff;
    count=0;
    // compute normal of origin face
    double n[3];
    ff->getNormal(n);
    double L=sqrt( dot(n,n) );
    n[0]=n[0]/L;
    n[1]=n[1]/L;
    n[2]=n[2]/L;
    // ray origin = centroid of origin face
    double lp[2][3];
    lp[0][0] = (ff->v[0]->pN[0]+
                ff->v[1]->pN[0]+
                ff->v[2]->pN[0])/3.0;
    lp[0][1] = (ff->v[0]->pN[1]+
                ff->v[1]->pN[1]+
                ff->v[2]->pN[1])/3.0;
    lp[0][2] = (ff->v[0]->pN[2]+
                ff->v[1]->pN[2]+
                ff->v[2]->pN[2])/3.0;
    // ray end = point on normal advanced from origin
    //			 a distance equal to 2*world width along each axis
    double del[3]={fabs(s.world[1]-s.world[0]),fabs(s.world[3]-s.world[2]),fabs(s.world[5]-s.world[4])};
    int big;
    biggest(del,big);
    lp[1][0] = lp[0][0]+n[0]*2*del[big];
    lp[1][1] = lp[0][1]+n[1]*2*del[big];
    lp[1][2] = lp[0][2]+n[2]*2*del[big];
    // for each face in object
    for(f_iterator i=f.begin();i!=f.end();i++){
      // if face in Object not same as origin face 
      if(*i!=ff){
        // if ray intersects face bounding box
        if(rayIntersectsBB(lp,*i,n)){
          // if hit, determine face-line intersection
          checkLineFaceIntersection(*i,lp,line_flag,poly_flag,poly_edge_flag);
          // check
          if(poly_edge_flag==true) {
            fff++;
            if(fff==f.end()){
              cout << "Object::isOutward: "
                    << "Error. Ray-tracing failed from all object <"
                    << name << "> faces.\n";
              exit(0);
            }
            break;
          }
          // does point intersect polygon
          if (poly_flag && line_flag) {count++;}
        }
      }
    }
  }while(poly_edge_flag==true);

  count = ceil(count);
  // if odd hit count then inward normal
  if(static_cast<int>(count)%2){return false;}
  // if even hit count then outward normal
  else {return true;}
}*/

void Object::evalAttributes(Container &c)
{
  cout << "Checking if object [" << name << "] is closed..........................";
  cout.flush();
  closed=isClosed();
  cout << "complete.\n";
  cout.flush();
  cout << "Checking if object [" << name << "] is manifold........................";
  cout.flush();
  manifold=isManifold();
  cout << "complete.\n";cout.flush();
  if(manifold==true)
  {
    cout << "Checking if object [" << name << "] faces are consistently oriented....";
    cout.flush();
    consistent=isConsistent();
    cout << "complete.\n";cout.flush();
    /* TODO port this old code
    if(consistent && closed)
    {
      cout << "Checking if object [" << name << "] faces are oriented outward.........";
      cout.flush();
      outward=isOutward(s);
      cout << "complete.\n";cout.flush();
    }
    */
  }
}

double dot(const double a[3],const double b[3]){
  return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
}

const double rand_vec[3]={
  0.236416584579274058342,
  0.927225593011826276779,
  0.389099507126957178116};

struct SortByRndVec // comparison function
{
  bool operator () (const key & a,
                    const key & b) const
  {
    // sort by random vector
    return  (v_vector[a].x*rand_vec[0]
            +v_vector[a].y*rand_vec[1]
            +v_vector[a].z*rand_vec[2]) < 
            (v_vector[b].x*rand_vec[0]
            +v_vector[b].y*rand_vec[1]
            +v_vector[b].z*rand_vec[2]);
  }
  static key min_value()
  {
    return std::numeric_limits<long long int>::min();
  }
  static key max_value()
  {
    return std::numeric_limits<long long int>::max();
  }
};

void Object::vertexDistin(void)
{
  cout << "\nObject::vertexDistin: aaa"; cout.flush();
  ///// check vertex distinguishability /////
  vec_key verts;
  key start = object_vector[object_id].v_start;
  key finish = object_vector[object_id].v_end;
  for(key i=start;i<finish+1;i++) { verts.push_back(i); }
  cout << "bbb"; cout.flush();
  // bound the main memory consumption by M
  // during sorting
  const unsigned M = 1024*1024*1024; // bytes
  // sort records by callers number
  stxxl::sort(verts.begin(),verts.end(),SortByRndVec(),M);
  cout << "ccc\n"; cout.flush();
  // for each vector element
  vk j;
  // DEBUG
  key num = finish-start+1;
  int k=1;
  // DEBUG
  for(vk i=verts.begin();i!=verts.end();i++)
  {
    //cout << "Object::vertexDistin: "<<k++<<" of "<<num<<endl;cout.flush();
    j=i;j++;
    if(j!=verts.end())
    {
      // if sequential pair of vector elements is not distinguishable
      if( !distinguishable(v_vector[*i].x,v_vector[*j].x) &&
          !distinguishable(v_vector[*i].y,v_vector[*j].y) &&
          !distinguishable(v_vector[*i].z,v_vector[*j].z) )
      {
        indistin_v.push_back(*i-object_vector[object_id].v_start);
        indistin_v.push_back(*j-object_vector[object_id].v_start);
      }
    }
  }
  cout << "Object::vertexDistin: finish\n"; cout.flush();
}

void Container::vertexAdjacentFaces()
{
  // DEBUG
  // sanityCheck();
  // DEBUG
  // for each vertex in object
  int num=object_vector[object_id].v_end
        -object_vector[object_id].v_start+1;
  for(int i=0;i<num;i++) // relative
  {
    //VE* ve_ptr = objects[object_id].GetVE(i);
    objects[object_id].ve.assign(i);
    //int n=static_cast<int>(ve_ptr->size());
    int n=static_cast<int>(objects[object_id].ve.size());
    //delete ve_ptr;
    // analyze vertex adjacent faces
    adjacent_face.x.push_back(n);
  }
  // process adjacent faces
  adjacent_face.process();
  // build adjacent face histogram
  af_create();
}

/* TODO port this old code
void Object::setAll(Vertex *vv,hmap_fi &group,int &g){
  // for each vertex adjacent face
  for(f_iterator i=vv->f.begin();i!=vv->f.end();i++){
    group[*i]=g;
  }
  g++;
}

void Object::getGroups(Vertex *vv,hmap_fi &group,set_i &s){
  s.clear();
  // for each vertex adjacent face
  for(f_iterator i=vv->f.begin();i!=vv->f.end();i++){
    s.insert(group[*i]);
  }
}

int Object::getLowest(set_i &s){
  int i = 100000;
  // for each element in set
  for(std::set<int,lti>::iterator j=s.begin();j!=s.end();j++){
    // if element is lower than i and not 0
    if(*j<i && *j){i=*j;}
  }
  return i;
}

void Object::replaceGroups(Vertex *vv,hmap_fi &group,int z){
  // for each vertex adjacent face
  for(f_iterator i=vv->f.begin();i!=vv->f.end();i++){
    group[*i]=z;
  }
}

void Object::setZero(Edge *ee,hmap_fi &group,int z){
  fi_iterator i=group.find(ee->f1);
  if((*i).second==0){group[ee->f1]=z;}
  // if second adjacent face
  if(ee->f2!=NULL){
    i=group.find(ee->f2);
    if((*i).second==0){group[ee->f2]=z;}
  }
  // look for edge in table
  ef_iterator k = findEdgeInTable_nonm_e(ee);
  // if found
  if(k!=nonm_e.end()){
    // grab face vector pointer
    vec_f *nfv = (*k).second;
    // for each adjacent face
    for(f_iterator j=nfv->begin();j!=nfv->end();j++){
      // if adjacent face has no group
      i=group.find(*j);
      if((*i).second==0){group[*j]=z;}
    }
  }
}

bool Object::removeSelectedFaces(vec_f &sf,vec_f &fv){
  bool flag=false;
  // for each selected face
  f_iterator i=sf.begin();
  while(i!=sf.end()){
    f_iterator j= find(fv.begin(),fv.end(),*i);
    // if selected face found in face vector
    if(j!=fv.end()){
      // remove face from face vector
      fv.erase(j);
      flag=true;
      i++;
    } else {
      // remove face from selected face vector
      //  so that the face is not used next round
      sf.erase(i);
    }
  }
  return flag;
}

void Object::getSelectedFaces(vec_f &sf,vec_f fv){
  // if selected faces vector is empty
  if(sf.empty()==true){
    // then grab first face from vector
    sf.push_back(fv.front());
  } else {
    // adjacent faces
    vec_f af;
    // for each selected face
    for(f_iterator i=sf.begin();i!=sf.end();i++){
      // grab adjacent faces of selected face
      if((*i)->e[0]->f1!=*i){af.push_back((*i)->e[0]->f1);}
      else 				  {af.push_back((*i)->e[0]->f2);}
      if((*i)->e[1]->f1!=*i){af.push_back((*i)->e[1]->f1);}
      else 				  {af.push_back((*i)->e[1]->f2);}
      if((*i)->e[2]->f1!=*i){af.push_back((*i)->e[2]->f1);}
      else 				  {af.push_back((*i)->e[2]->f2);}
    }
    // keep unique faces
    sort(af.begin(),af.end());
    f_iterator i=unique(af.begin(),af.end());
    sf.assign(af.begin(),i);
  }
}
int Object::countComponents(void){
  // initialize group number
  int g=0;
  ///// create hashed map face*->integer (group #) /////
  // hash_map<Face*,int>
  hmap_fi group;
  // for each face
  for(f_iterator i=f.begin();i!=f.end();i++){
    group[*i]=g;
  }
  // increment group number
  g++;
  // set of integers
  set_i s;
  bool changes = true;
  while(changes==true){
    changes=false;
    // for each vertex in object
    for(v_iterator i=v.begin();i!=v.end();i++){
      // get vertex adjacent face groups
      getGroups(*i,group,s);
      std::set<int,lti>::iterator j=s.begin();
      // if no adjacent face group has been set
      if(s.size()==1 && *j==0){
        // set all adjacent faces to next available group #
        setAll(*i,group,g);
        changes=true;
      } else if(s.size()>1){ // more than one group
        // identify lowest group # larger than 0 in set
        int z=getLowest(s);
        // replace the larger group # with lowest group # in all faces
        replaceGroups(*i,group,z);
        changes=true;
      } 
    }
  }
  // analyze map
  s.clear();
  // for each element in group
  for(fi_iterator i=group.begin();i!=group.end();i++){
    s.insert((*i).second);
  }
  return s.size();
}*/

void F::getNormal(key object_id,double n[3])
{
  V *v1_ptr = objects[object_id].GetV(get_v1_rel());
  V *v2_ptr = objects[object_id].GetV(get_v2_rel());
  V *v3_ptr = objects[object_id].GetV(get_v3_rel());
  double uX, uY, uZ, vX, vY, vZ;
  uX = v2_ptr->get_x()-v1_ptr->get_x();
  uY = v2_ptr->get_y()-v1_ptr->get_y();
  uZ = v2_ptr->get_z()-v1_ptr->get_z();
  vX = v3_ptr->get_x()-v1_ptr->get_x();
  vY = v3_ptr->get_y()-v1_ptr->get_y();
  vZ = v3_ptr->get_z()-v1_ptr->get_z();
  delete v1_ptr;
  delete v2_ptr;
  delete v3_ptr;
  // compute cross product (u x v)
  n[0] = uY*vZ-uZ*vY;
  n[1] = uZ*vX-uX*vZ;
  n[2] = uX*vY-uY*vX;
}

double F::getAspectRatio(key object_id)
{
  V *v1_ptr = objects[object_id].GetV(get_v1_rel());
  V *v2_ptr = objects[object_id].GetV(get_v2_rel());
  V *v3_ptr = objects[object_id].GetV(get_v3_rel());
  /* Make triangle edge vectors */
  double va[3]={v2_ptr->get_x()-v1_ptr->get_x(),v2_ptr->get_y()-v1_ptr->get_y(),v2_ptr->get_z()-v1_ptr->get_z()};
  double vb[3]={v3_ptr->get_x()-v2_ptr->get_x(),v3_ptr->get_y()-v2_ptr->get_y(),v3_ptr->get_z()-v2_ptr->get_z()};
  double vc[3]={v1_ptr->get_x()-v3_ptr->get_x(),v1_ptr->get_y()-v3_ptr->get_y(),v1_ptr->get_z()-v3_ptr->get_z()};
  double vbase[3]={0,0,0};
  double vopp[3]={0,0,0};
  /* Find length of longest edge */
  double lmax=-DBL_MAX;
  double la=sqrt(dot(va,va));
  double lb=sqrt(dot(vb,vb));
  double lc=sqrt(dot(vc,vc));
  if (la>lmax)
  {
    lmax=la;
    vbase[0]=va[0];
    vbase[1]=va[1];
    vbase[2]=va[2];
    vc[0]=v3_ptr->get_x()-v1_ptr->get_x();
    vc[1]=v3_ptr->get_y()-v1_ptr->get_y();
    vc[2]=v3_ptr->get_z()-v1_ptr->get_z();
    vopp[0]=vc[0];
    vopp[1]=vc[1];
    vopp[2]=vc[2];
  }
  if (lb>lmax)
  {
    lmax=lb;
    vbase[0]=vb[0];
    vbase[1]=vb[1];
    vbase[2]=vb[2];
    va[0]=v1_ptr->get_x()-v2_ptr->get_x();
    va[1]=v1_ptr->get_y()-v2_ptr->get_y();
    va[2]=v1_ptr->get_z()-v2_ptr->get_z();
    vopp[0]=va[0];
    vopp[1]=va[1];
    vopp[2]=va[2];
  }
  if (lc>lmax)
  {
    lmax=lc;
    vbase[0]=vc[0];
    vbase[1]=vc[1];
    vbase[2]=vc[2];
    vb[0]=v2_ptr->get_x()-v3_ptr->get_x();
    vb[1]=v2_ptr->get_y()-v3_ptr->get_y();
    vb[2]=v2_ptr->get_z()-v3_ptr->get_z();
    vopp[0]=vb[0];
    vopp[1]=vb[1];
    vopp[2]=vb[2];
  }
  delete v1_ptr;
  delete v2_ptr;
  delete v3_ptr;

  /* Find shortest altitude */
  double ll = sqrt(dot(vbase,vbase));
  vbase[0]=vbase[0]/ll;
  vbase[1]=vbase[1]/ll;
  vbase[2]=vbase[2]/ll;
  double dot_prod = dot(vbase,vopp);
  double alt[3]={vopp[0]-(dot_prod*vbase[0]),
    vopp[1]-(dot_prod*vbase[1]),
    vopp[2]-(dot_prod*vbase[2])};
  double amin=sqrt(dot(alt,alt));

  return lmax/amin;
}

void Container::areaAspectRatio()
{
  // for each face in object
  int num=object_vector[object_id].f_end
        -object_vector[object_id].f_start+1;
  for(int i=0;i<num;i++) // relative
  {
    // compute face normal vector
    //F* f_ptr = objects[object_id].GetF(i);
    objects[object_id].f.assign(i);
    F* f_ptr = &objects[object_id].f;
    double n[3];
    f_ptr->getNormal(object_id,n);
    double ar = f_ptr->getAspectRatio(object_id);
    //delete f_ptr;
    // compute face area = half normal vector length
    double aa=sqrt(dot(n,n))/2.0;
    ///// face area /////
    area.x.push_back(aa);
    ///// aspect ratio /////
    aspect_ratio.x.push_back(ar);
  }
  // process_area
  area.process();
  // build face area histogram
  a_create();
  aspect_ratio.process();
  // build aspect ratio histogram
  ar_create();
}

double E::getSqLength(key object_id)
{
  V *v1_ptr = objects[object_id].GetV(get_v1());
  V *v2_ptr = objects[object_id].GetV(get_v2());
  double a = (v1_ptr->get_x()-v2_ptr->get_x())*
              (v1_ptr->get_x()-v2_ptr->get_x())+
              (v1_ptr->get_y()-v2_ptr->get_y())*
              (v1_ptr->get_y()-v2_ptr->get_y())+
              (v1_ptr->get_z()-v2_ptr->get_z())*
              (v1_ptr->get_z()-v2_ptr->get_z());
  delete v1_ptr;
  delete v2_ptr;
  return a;
}

void Container::processEdgeLengths()
{
  // for each edge in object
  int num=object_vector[object_id].e_end
        -object_vector[object_id].e_start+1;
  for(int i=0;i<num;i++) // relative
  {
    E* e_ptr = objects[object_id].GetE(i);
    //e_ptr->print(object_id);
    double l = sqrt(e_ptr->getSqLength(object_id));
    delete e_ptr;
    // process_edge_lengths(i,l);
    // check distinguishability
    e_ptr = objects[object_id].GetE(i);
    V *v1_ptr = objects[object_id].GetV(e_ptr->get_v1());
    V *v2_ptr = objects[object_id].GetV(e_ptr->get_v2());
    if( !distinguishable(v1_ptr->get_x(),v2_ptr->get_x()) &&
        !distinguishable(v1_ptr->get_y(),v2_ptr->get_y()) &&
        !distinguishable(v1_ptr->get_z(),v2_ptr->get_z()) 
      ){ indistin_e.push_back(i); }
    delete e_ptr;
    delete v1_ptr;
    delete v2_ptr;
    // add to vector
    edge_length.x.push_back(l);
  }
  edge_length.process();
  el_create();
}

double E::getAngle(key object_id)
{
  // get outward normals of edge faces
  //F* f1_ptr = objects[object_id].GetF(get_f1());
  F* f1_ptr = objects[object_id].GetF(get_f1(get_v1()));
  //objects[object_id].f.assign(get_f1(get_v1()));
  //F* f1_ptr = &objects[object_id].f;
  //F* f2_ptr = objects[object_id].GetF(get_f2());
  F* f2_ptr = objects[object_id].GetF(get_f2(get_v1()));
  //objects[object_id].f.assign(get_f2(get_v1()));
  //F* f2_ptr = &objects[object_id].f;
  key o2 = f2_ptr->get_o(get_v1(),get_v2());
  double n1[3],n2[3];
  f1_ptr->getNormal(object_id,n1);
  f2_ptr->getNormal(object_id,n2);
  //delete f1_ptr;
  //delete f2_ptr;
  // compute the cosine of angle between normals
  double normal_angle_cosine = dot(n1,n2)/sqrt(dot(n1,n1))/sqrt(dot(n2,n2));
  if(normal_angle_cosine >= 1) { return PI; }
  if (normal_angle_cosine <= -1) { return 0; } // normal_angle = PI, gamma == 0 or 2PI
  // normal_angle = acos(normal_angle_cosine);
  // use the edge itself as a reference vector
  V *o2_ptr = objects[object_id].GetV(o2);
  V *v2_ptr = objects[object_id].GetV(get_v2()); // note: get_v1() could be v2 instead
  double refvec[3] = {v2_ptr->get_x()-o2_ptr->get_x(),
                      v2_ptr->get_y()-o2_ptr->get_y(),
                      v2_ptr->get_z()-o2_ptr->get_z()};
  delete o2_ptr;
  delete v2_ptr;
  // dot product of refvec and n1
  double d = dot(refvec,n1);
  if(d==0.0) { return PI; }
  //cout << "this one is legit:";
  return PI+d/fabs(d)*acos(normal_angle_cosine);
  // I EXPECT 0 <= angle < 2*PI
}

void Container::computeEdgeAngles()
{
  // for each edge in object
  key num=object_vector[object_id].e_end
        -object_vector[object_id].e_start+1;
  for(key i=0;i<num;i++) // relative
  {
    E* e_ptr = objects[object_id].GetE(i);
    // if edge has exactly two adjacent faces
    // look for edge in nonmanifold vector
    if(std::find(nonman_e.begin(),nonman_e.end(),i)==nonman_e.end())
    {
      if((e_ptr->get_f2(e_ptr->get_v1())<0)==false)
      {
        double angle = e_ptr->getAngle(object_id)*180/PI; // degrees
        // DEBUG
        //cout.precision(12);
        //cout << angle << endl;
        // DEBUG
        edge_angle.x.push_back(angle);
      }
    }
    delete e_ptr;
  }
  edge_angle.process();
  ea_create();
}

void Object::computeVolume()
{
  vol = 0.0;
  // for each face in object
  int num=object_vector[object_id].f_end
        -object_vector[object_id].f_start+1;
  for(int i=0;i<num;i++) // relative
  {
    //F* f_ptr = GetF(i);
    f.assign(i);
    F* f_ptr = &f;
    V *v1_ptr = GetV(f_ptr->get_v1_rel());
    V *v2_ptr = GetV(f_ptr->get_v2_rel());
    V *v3_ptr = GetV(f_ptr->get_v3_rel());
    double x1=v1_ptr->get_x();
    double y1=v1_ptr->get_y();
    double z1=v1_ptr->get_z();
    double x2=v2_ptr->get_x();
    double y2=v2_ptr->get_y();
    double z2=v2_ptr->get_z();
    double x3=v3_ptr->get_x();
    double y3=v3_ptr->get_y();
    double z3=v3_ptr->get_z();
    //delete f_ptr;
    delete v1_ptr;
    delete v2_ptr;
    delete v3_ptr;
    /* compute determinant of oriented triangle */
    double det=x1*(y2*z3-y3*z2)+x2*(y3*z1-y1*z3)+x3*(y1*z2-y2*z1);
    vol+=det;
  }
  vol=vol/6.0;
}
/* TODO port this old code
void Object::computeGenus(void)
{
  int v_size = object_vector[object_id].v_end
              -object_vector[object_id].v_start+1;
  int e_size = object_vector[object_id].e_end
              -object_vector[object_id].e_start+1;
  int f_size = object_vector[object_id].f_end
              -object_vector[object_id].f_start+1;
  int num = v_size-e_size+f_size;
  if (num%2) {
    cout << "Round off error in genus computation. "
          << "(#vertices + #faces - #edges) is not evenly divisible by 2.\n";
    cout << "v-e+f = " << num << endl;
    cout << v_size << "-" << e_size << "+" << f_size << endl;
    exit(0);
  } else {
    genus = num_sep-num/2;
  }
}*/

void Object::evalCharacteristics(Container& c)
{
  // DEBUG
  //sanityCheck();
  // DEBUG
  if(distvert==true){
    cout << "Check if vertices are distinguishable for object ["
          << name << "]......";
    cout.flush();
    vertexDistin();
    cout << "complete.\n";cout.flush();
    return;
  }
  // vertices
  cout << "Bound object [" << name << "]..........................................";
  cout.flush();
  boundObject();
  cout << "complete.\n";cout.flush();
  cout << "Check if vertices are distinguishable for object [" << name << "]......";
  cout.flush();
  // DEBUG
  // sanityCheck();
  // DEBUG
  vertexDistin();
  cout << "complete.\n";cout.flush();
  cout << "Analyze vertex adjacent faces for object [" << name << "]..............";
  cout.flush();
  c.vertexAdjacentFaces();
  cout << "complete.\n";cout.flush();
  // faces
  /* TODO port this old code
  cout << "Identify separate components of object [" << name << "]................";
  cout.flush();
  num_sep=countComponents();
  cout << "complete.\n";cout.flush();
  */
  cout << "Compute face area and aspect ratio for object [" << name << "].........";
  cout.flush();
  c.areaAspectRatio();
  cout << "complete.\n";cout.flush();
  /* TODO port this old code
  // if not batch mode
  if(cs.interf==false){
  cout << "Find intersecting faces for object [" << name << "]....................";
  cout.flush();
  findIntersectingFaces(c,s);
  cout << "complete.\n";cout.flush();
  }
  // edges
  cout << "Identify boundaries for object [" << name << "]........................";
  cout.flush();
  countBoundaries();
  cout << "complete.\n";cout.flush();
  */
  cout << "Analyze edge lengths for object [" << name << "].......................";
  cout.flush();
  c.processEdgeLengths();
  cout << "complete.\n";cout.flush();
  if(manifold && consistent)
  {
    // manifold and consistently oriented face normals
    cout << "Analyze edge angles for object [" << name << "]........................";
    cout.flush();
    c.computeEdgeAngles();
    cout << "complete.\n";cout.flush();
    if(closed){
      // closed, manifold, and consistently oriented face normals
      cout << "Compute volume of object [" << name << "]..............................";
      cout.flush();
      computeVolume();
      cout << "complete.\n";cout.flush();
      /* TODO add this when countComponents is added
      // if number of components == 1
      // FUTURE IMPROVEMENT: only require orientable, not oriented
      // FUTURE IMPROVEMENT: separate components and compute genus of each component
      if(num_sep==1 && orphan.empty()){
        cout << "Compute genus of object [" << name << "]...............................";
        cout.flush();
        computeGenus();
        cout << "complete.\n";cout.flush();
      }
      */
    }
  }
}

void Object::analyze(Container &c)
{
  // DEBUG
  //sanityCheck();
  // DEBUG
  evalAttributes(c);
  if(attr==false) evalCharacteristics(c);
}

void Container::update()
{
  key oid = getObjectId();
  // TODO num_sep+=oo->num_sep;
  if(objects[oid].getManifold()==true){num_man[0]++;}
  else if(objects[oid].getManifold()==false && 
          objects[oid].getClosed()==true) {num_man[1]++;}
  else {num_man[2]++;}
  num_bou+=objects[object_id].getBoundary();
  num_indistin+=indistin_v.size();
  if     (objects[oid].getManifold()==true &&
  objects[oid].getConsistent()==true) {num_cons[0]++;}
  else if(objects[oid].getManifold()==true &&
          objects[oid].getConsistent()==false){num_cons[1]++;}
  else {num_cons[2]++;}
  vol+=objects[oid].getVol();
  if(objects[oid].getClosed()==true){num_clo++;}
  else {num_op++;}
  num_bor+=border.size();
  num_nonman_v+=nonman_v.size();
  num_nonman_e+=nonman_e.size();
  num_flip+=flipped.size();
  if     (objects[oid].getConsistent()==true && 
          objects[oid].getOutward()==true) {num_out[0]++;}
  else if(objects[oid].getConsistent()==true && 
          objects[oid].getOutward()==false){num_out[1]++;}
  else {num_out[2]++;}
  num_orph+=orphan.size();
  num_mis+=miss_v;
  num_deg+=degen_f;
  num_dupl_v+=dupl_v;
  num_dupl_f+=dupl_f;
}

void Object::printChars(Container &c)
{
  if(distvert==true)
  {
    ///////////////////// indistinguishable vertices
    if(indistin_v.empty())
    {
      cout << "    # indistinguishable vertices: none\n";
    }
    else
    {
      cout << "    # indistinguishable vertices: "
            << indistin_v.size() << endl;
      // if -p option, print offending
      if(my_print==true)
      {
        int j=1;
        // for each indistinguishable vertex
        for(vk i=indistin_v.begin();i!=indistin_v.end();i++)
        {
          cout << "    #  indistinguishable vertices: vertex "
                << j++ << endl;
          V *v_ptr = GetV(*i);
          if(strcmp(style,"cp")==false)
          {
            v_ptr->printCP();
          }
          else
          {
            v_ptr->print(object_id);
          }
          delete v_ptr;
        }
      }
    }
    return;
  }
  int v_size = object_vector[object_id].v_end
        -object_vector[object_id].v_start+1;
  int e_size = object_vector[object_id].e_end
        -object_vector[object_id].e_start+1;
  int f_size = object_vector[object_id].f_end
        -object_vector[object_id].f_start+1;
  cout << "\nMESH CHARACTERISTICS\n\n";
  cout << "    # vertices: " << v_size << endl
        << "    # faces: " << f_size << endl
        << "    # edges: " << e_size << endl;
  // TODO add component counting
  //      << "    # components: " << num_sep << endl;
  if (manifold==false)
  {
    cout << "    # boundaries: Since object is nonmanifold,\n"
          << "    # boundaries: the number of boundaries may be underestimated.\n";
  }
  //////////////////// borders
  if (border.empty())
  {
    cout << "    # boundaries: none\n";
  }
  else
  {
    cout << "    # boundaries: " << border.size() << endl;
    // if -p option, print offending
    if(my_print==true)
    {
      int j=1;
      // for each border edge
      for(vk i=border.begin();i!=border.end();i++)
      {
        cout << "    # boundaries: boundary edge " << j++ << endl;
        E *e_ptr = GetE(*i);
        if(strcmp(style,"cp")==false)
        {
          e_ptr->printCP(object_id);
        }
        else
        {
          e_ptr->print(object_id);
        }
        delete e_ptr;
      }
    }
  }
  ///////////////////// indistinguishable vertices
  if(indistin_v.empty()==true)
  {
    cout << "    # indistinguishable vertices: none\n";
  }
  else
  {
    cout << "    # indistinguishable vertices: " << indistin_v.size() << endl;
    //	if -p option, print offending
    if(my_print==true)
    {
      int j=1;
      // for each indistinguishable vertex
      for(vk i=indistin_v.begin();i!=indistin_v.end();i++)
      {
        cout << "    #  indistinguishable vertices: vertex " << j++ << endl;
        V *v_ptr = GetV(*i);
        if(strcmp(style,"cp")==false)
        {
          v_ptr->printCP();
        }
        else
        {
          v_ptr->print(object_id);
        }
        delete v_ptr;
      }
    }
  }
  //////////////// if volume computed
  if(closed==true && manifold==true && consistent==true)
  {
    cout << "    object volume: [(data units)^3]" << endl;
    cout << "    object volume: " << vol << endl;
  }
  else
  {
    cout << "    object volume: not computed, since ";
    if(closed==false){cout << "not closed,";}
    if(consistent==false){cout << "not consistent,";}
    if(manifold==false){cout << "not manifold";}
    cout << endl;
  }
  /* TODO add this after adding countComponents
  /////////////// if genus computed
  if(closed==true && manifold==true && consistent==true && num_sep==1 && orphan.empty())
  {
    cout << "    object genus: " << genus << endl;
  }
  else
  {
    cout << "    object genus: not computed, since ";
    if(closed==false){cout << "not closed,";}
    if(consistent==false){cout << "not consistent,";}
    if(manifold==false){cout << "not manifold,";}
    if(num_sep>1){cout << "#components=" << num_sep << ",";}
    if(!orphan.empty()){cout << "orphan vertices were found";}
    cout << endl;
  }*/
  //////////////// bounding box
  cout << "    bounding box: [data units]\n";
  cout << "    bounding box: [xmin,ymin,zmin][xmax,ymax,zmax]\n";
  cout << "    bounding box: ["
        << bb[0] << ","
        << bb[1] << ","
        << bb[2] << "]["
        << bb[3] << ","
        << bb[4] << ","
        << bb[5] << "]" << endl;
  ////////////////// edges with indistinguishable vertices
  if(indistin_e.empty()==true)
  {
    cout << "    # edges with indistinguishable vertices: none\n";
  }
  else
  {
    cout << "    # edges with indistinguishable vertices: "
          << indistin_e.size() << endl;
  }
  //	if -p option, print offending
  if(my_print==true)
  {
    if(indistin_e.empty()==false)
    {
      // for each afflicted edge
      for(vk i=indistin_e.begin();i!=indistin_e.end();i++)
      {
        E *e_ptr = GetE(*i);
        if(strcmp(style,"cp")==false)
        {
          e_ptr->printCP(object_id);
        }
        else
        {
          e_ptr->print(object_id);
        }
        delete e_ptr;
      }
    }
  }
  /* TODO add this after adding findIntersectingFaces
  /////////////// intersecting faces
  // if not batch mode
  if(cs.interf==false){
    // intersecting faces
    if (intf.empty()){
      cout << "    # intersecting faces: none\n\n";
    } else {
      cout << "    # intersecting faces: " << intf.size() << endl;
      //	if -p option, print offending
      if(my_print==true){
        int j=1;
        // for each intersected face
        for(ff_iterator i=intf.begin();i!=intf.end();i++){
          cout << "    # intersecting faces: intersected face " << j++ << endl;
          // print intersected face
          if(strcmp(style,"cp")==false){
            (*i).first->printFaceCP();
          } else {
            (*i).first->printFace((*i).first->v[0]->o);
          }
          // keep unique list of intersecting faces
          sort((*(*i).second).begin(),(*(*i).second).end());
          f_iterator new_end = unique((*(*i).second).begin(),(*(*i).second).end());
          (*(*i).second).assign((*(*i).second).begin(),new_end);
          // print intersecting faces
          for(f_iterator k=(*(*i).second).begin();k!=(*(*i).second).end();k++){
            if(strcmp(style,"cp")==false){
              (*k)->printFaceCP();
            } else {
              (*k)->printFace((*k)->v[0]->o);
              cout << endl;
            }
          }
        }
        cout << endl << endl;
      }
    }
  }*/

  //////////////// vertex adjacent faces
  cout << "    Vertex adjacent face statistics [faces]:" << endl;
  c.print_adjacent_face_stats();
  cout << "    Vertex adjacent face histogram [faces]:" << endl;
  c.print_adjacent_face_histo();
  cout << endl;
  ///////////////// face area
  cout << "    Face area statistics [(data units)^2]:" << endl;
  cout << "       total    " << c.area_sum() << endl;
  c.print_area_stats();
  cout << "    Face area histogram [(data units)^2]:" << endl;
  c.print_area_histo();
  cout << endl;
  /////////////////// face aspect ratio
  cout << "    Face aspect ratio statistics [unitless]:" << endl;
  c.print_aspect_ratio_stats();
  cout << "    Face aspect ratio histogram [unitless]:" << endl;
  c.print_aspect_ratio_histo();
  cout << "      (Aspect ratio is longest edge "
        << "divided by shortest altitude)\n";
  /* TODO port this old code
  // if face aspect ratio threshold specified
  if(my_signal[0]==true)
  {
    if(bad_aspect.empty()==true){
      cout << "    # faces with bad aspect ratios: none" << endl;
    } else {
      cout << "    # faces with bad aspect ratios: "
            << bad_aspect.size() << endl;
    }
  }
  //	if -p option, print offending
  if(my_print==true){
    if(bad_aspect.empty()==false){
      // print faces with aspect ratios that violate threshold
      Face_Pair fp(this);
      for(fd_iterator i=bad_aspect.begin();i!=bad_aspect.end();i++){
        cout << "    face aspect ratio: " << (*i).second << endl;
        if(strcmp(style,"cp")==false){
          (*i).first->printFaceCP();
        } else {
          (*i).first->printFace((*i).first->v[0]->o);
          cout << endl;
        }
        fp.processBadFace((*i).first);
      }
    }
  }
  cout << "\n";
  */
  ////////////// edge length
  cout << "    Edge length statistics [data units]:" << endl;
  c.print_edge_length_stats();
  cout << "    Edge length histogram [data units]:" << endl;
  c.print_edge_length_histo();
  cout << endl;
  /* TODO port this old code
  // if edge length threshold specified
  if(cs.signal[3]==true || cs.signal[4]==true){
    if(bad_length.empty()==true){
      cout << "    # edges with bad lengths: none" << endl;
    } else {
      cout << "    # edges with bad lengths: "
            << bad_length.size() << endl;
    }
  }
  //	if -p option, print offending
  if(my_print==true){
    if(bad_length.empty()==false){
      // print edges with lengthss that violate threshold
      for(ed_iterator i=bad_length.begin();i!=bad_length.end();i++){
        cout << "    edge length: " << (*i).second << endl;
        if(strcmp(style,"cp")==false){
          (*i).first->printEdgeCP();
        } else {
          (*i).first->printEdge((*i).first->vv1->o->name);
          cout << endl;
        }
      }
    }
  }*/
  //////////////////// if edge angles computed
  if(manifold==true && consistent==true){
    cout << "    Edge angle statistics [degrees]:" << endl;
    c.print_edge_angle_stats();
    cout << "    Edge angle histogram [degrees]:" << endl;
    c.print_edge_angle_histo();
    cout << endl;
    /* TODO port this old code
    // if edge angle threshold specified
    if(cs.signal[1]==true || cs.signal[2]==true){
      if(bad_angle.empty()==true){
        cout << "    # edges with bad angles: none" << endl;
      } else {
        cout << "    # edges with bad angles: "
              << bad_angle.size() << endl;
      }
    }
    //	if -p option, print offending
    if(my_print==true){
      if(bad_angle.empty()==false){
        // print edges with angles that violate threshold
        for(ed_iterator i=bad_angle.begin();i!=bad_angle.end();i++){
          cout << "    edge angle: " << (*i).second << endl;
          if(strcmp(style,"cp")==false){
            (*i).first->printEdgeCP();
          } else {
            (*i).first->printEdge((*i).first->vv1->o->name);
            cout << endl;
          }
        }
      }
    }
    */
  }
  else
  {
    cout << "    edge angles: not computed, since ";
    if(consistent==false){cout << "not consistent,";}
    if(manifold==false){cout << "not manifold";}
    cout << endl;
  }
}

void Object::printIntegrity(Container &c)
{
  cout << "\nMESH FILE INTEGRITY\n\n";
  // orphan vertices
  if(orphan.empty()==true)
  {
    cout << "    # orphan vertices: none\n";
  }
  else
  {
    cout << "    # orphan vertices: " << orphan.size() << endl;
    //	if -p option, print offending
    if(my_print==true)
    {
      int j=1;
      // for each orphan vertex
      for(sk i=orphan.begin();i!=orphan.end();i++)
      {
        cout << "    # orphan vertices: orphan vertex " << j++ << endl;
        V* v_ptr = objects[object_id].GetV(*i);
        if(strcmp(style,"cp")==false)
        {
          v_ptr->printCP();
        }
        else
        {
          v_ptr->print(object_id);
        }
        delete v_ptr;
      }
    }
  }
  // missing vertices
  if(miss_v==0)
  {
    cout << "    # missing vertices: none\n";
  }
  else
  {
    cout << "    # missing vertices: " << miss_v << endl;
  }
  // degenerate faces
  if(degen_f==0)
  {
    cout << "    # degenerate faces: none\n";
  }
  else
  {
    cout << "    # degenerate faces: " << degen_f << endl;
  }
  // duplicate indices
  if(dupl_v==0)
  {
    cout << "    # duplicate vertex indices: none\n";
  }
  else
  {
    cout << "    # duplicat vertex indices: " << dupl_v << endl;
  }
  if(dupl_f==0)
  {
    cout << "    # duplicate face indices: none\n";
  }
  else
  {
    cout << "    # duplicat face indices: " << dupl_f << endl;
  }
  // contiguous numbering
  if(seq_v==true)
  {
    cout << "    contiguous vertex indexing from 1: yes\n";
  }
  else
  {
    cout << "    contiguous vertex indexing from 1: no\n";
    cout << "    contiguous vertex indexing from 1: first nonsequential index\n";
    cout << "    contiguous vertex indexing from 1: " << disc_v << endl;
  }
}

void Object::printAttr(Container &c)
{
  cout << "\nMESH ATTRIBUTES\n\n";
  // closed
  if(closed==true)
  {
    cout << "    mesh is closed: yes\n";
  }
  else
  {
    cout << "    mesh is closed: no\n";
    //	if -p option, print offending
    if(my_print==true)
    {
      int j=1;
      cout << "    mesh is closed: # border edges - " << border.size() << endl;
      // for each border edge
      for(vk i=border.begin();i!=border.end();i++)
      {
        cout << "    mesh is closed: border edge # " << j++ << endl;
        E *e_ptr = objects[object_id].GetE(*i);
        if(strcmp(style,"cp")==false)
        {
          e_ptr->printCP(object_id);
        }
        else
        {
          e_ptr->print(object_id);
        }
        delete e_ptr;
      }
    }
  }
  // manifold
  if(manifold==true)
  {
    cout << "    mesh is manifold: yes\n";
  }
  else if(manifold==false && closed==true)
  {
    cout << "    mesh is manifold: no\n";
    //	if -p option, print offending
    if(my_print==true)
    {
      if(nonman_v.empty()==true)
      {
        cout << "    mesh is manifold: # nonmanifold vertices - none\n";
      } else {
        int j=1;
        cout << "    mesh is manifold: # nonmanifold vertices - " << nonman_v.size() << endl;
        // for each nonmanifold vertex
        //for(vk i=nonman_v.begin();i!=nonman_v.end();i++){
        for(new_vk2 i=nonman_v.begin();i!=nonman_v.end();i++){
          cout << "    mesh is manifold: nonmanifold vertex # " << j++ << endl;
          V *v_ptr = objects[object_id].GetV(*i);
          if(strcmp(style,"cp")==false)
          {
            v_ptr->printCP();
          }
          else
          {
            v_ptr->print(object_id);
          }
          delete v_ptr;
        }
      }
      if(nonman_e.empty()==true)
      {
        cout << "    mesh is manifold: # nonmanifold edges - none\n";
      }
      else
      {
        int j=1;
        cout << "    mesh is manifold: # nonmanifold edges - " << nonman_e.size() << endl;
        // for each nonmanifold edge
        //for(vk i=nonman_e.begin();i!=nonman_e.end();i++){
        for(new_vk i=nonman_e.begin();i!=nonman_e.end();i++){
          cout << "    mesh is manifold: nonmanifold edge # " << j++ << endl;
          E *e_ptr = objects[object_id].GetE(*i);
          if(strcmp(style,"cp")==false)
          {
            e_ptr->printCP(object_id);
          }
          else
          {
            e_ptr->print(object_id);
          }
          delete e_ptr;
        }
      }
    }
  }
  else
  {
    cout << "    mesh is manifold: undefined since mesh is open\n";
  }
  // consistent
  if(manifold==false)
  {
    cout << "    mesh has consistently oriented face normals: undefined since not manifold\n";
  }
  else
  {
    if(consistent==true)
    {
      cout << "    mesh has consistently oriented face normals: yes\n";
    }
    else
    {
      cout << "    mesh has consistently oriented face normals: no\n";
      //	if -p option, print offending
      if(my_print==true)
      {
        int j=1;
        cout << "    mesh has consistently oriented face normals: # flipped edges - " << flipped.size() << endl;
        // for each flipped edge
        // TODO
        // FUTURE IMPROVEMENT: if edge is nonmanifold then exclude from flipped list
        for(vk i=flipped.begin();i!=flipped.end();i++){
          cout << "    mesh has consistently oriented face normals: flipped edge # " << j++ << endl;
          E *e_ptr = objects[object_id].GetE(*i);
          if(strcmp(style,"cp")==false)
          {
            e_ptr->printCP(object_id);
          }
          else
          {
            e_ptr->print(object_id);
          }
          delete e_ptr;
        }
      }
    }
  }
  // outward
  /* TODO port this old code
     if(manifold==false || consistent==false || closed==false)
     {
     cout << "    mesh has outward oriented face normals: uncomputable since ";
     if(closed==false){cout << "not closed,";}
     if(consistent==false){cout << "not consistent,";}
     if(manifold==false){cout << "not manifold";}
     cout << endl;
     }
     else
     {
     if(outward==true)
     {
     cout << "    mesh has outward oriented face normals: yes\n";
     }
     else
     {
     cout << "    mesh has outward oriented face normals: no\n";
     }
     }*/
}

void Object::print(Container &c)
{
  if(distvert==true){
    printChars(c);
    return;
  }
  printIntegrity(c);
  if(checkIntegrity()==false)
  {
    cout << "\n\nWarning: Attributes and "
          << "characteristics were not evaluated,\n"
          << " since mesh file failed the integrity check.\n\n";
    return;
  }
  printAttr(c);
  if(attr==false)
  { 
    printChars(c);
  }
}

////////////////////////////////////////////////////
////////////////////////////////////////////////////
////////////////////////////////////////////////////
/*
void Container::processEdgeLengths()
{
  // for each element in area vector
  for(d_iterator i=edge_length.x.begin();i!=edge_length.x.end();i++){
    //		double l = (*i)->l;
    //		edge_length.n++;
    edge_length.sum+=*i;
    edge_length.sum2+=(*i)*(*i);
    edge_length.total+=*i;
    if(*i<edge_length.min) edge_length.min=*i;
    if(*i>edge_length.max) edge_length.max=*i;
  }
  edge_length.createHistogram();
}

//void Container::areaAspectRatio()
void Container::areaAspectRatio(key first,key last)
{
  // for each element in area vector from first to last
  key gate=0;
  for(d_iterator i=area.x.begin();i!=area.x.end();i++){
    if(gate<x_first==false && gate>x_last==false)
    {
      area.sum+=*i;
      area.sum2+=(*i)*(*i);
      area.total+=*i;
      if(*i<area.min) {area.min=*i;}
      if(*i>area.max) {area.max=*i;}
    }
    gate++;
  }

  // for each element in aspect ratio vector
  gate=0;
  for(d_iterator i=aspect_ratio.x.begin();i!=aspect_ratio.x.end();i++){
    if(gate<x_first==false && gate>x_last==false)
    {
      aspect_ratio.sum+=*i;
      aspect_ratio.sum2+=(*i)*(*i);
      aspect_ratio.total+=*i;
      if(*i<aspect_ratio.min) aspect_ratio.min=*i;
      if(*i>aspect_ratio.max) aspect_ratio.max=*i;
    }
    gate++;
  }
  // build face area histogram
  area.createHistogram();
  // build aspect ratio histogram
  aspect_ratio.createAspectRatioHistogram();
}

void Container::vertexAdjacentFaces()
{
  //	mmap_iv af;
  // for each element in vector
  for(d_iterator i=adjacent_face.x.begin();i!=adjacent_face.x.end();i++){
    //		int c=(*i)->f.size();
    //		af.insert(std::make_pair(c,*i));
    //		adjacent_face.n++;
    adjacent_face.sum+=*i;
    adjacent_face.sum2+=(*i)*(*i);
    adjacent_face.total+=*i;
    if(*i<adjacent_face.min) {adjacent_face.min=*i;}
    if(*i>adjacent_face.max) {adjacent_face.max=*i;}
    // add to vector
    //		adjacent_face.x.push_back(c);
  }
  // build adjacent face histogram
  adjacent_face.createAdjacentFaceHistogram();
}

void Container::computeEdgeAngles()
{
  if(edge_angle.x.empty()==false)
  {
    // for each element in area vector
    for(d_iterator i=edge_angle.x.begin();i!=edge_angle.x.end();i++){
      edge_angle.sum+=*i;
      edge_angle.sum2+=(*i)*(*i);
      edge_angle.total+=*i;
      if(*i<edge_angle.min) edge_angle.min=*i;
      if(*i>edge_angle.max) edge_angle.max=*i;
    }
    edge_angle.createHistogram();
  }
}*/

void Container::printCumulative()
{
  cout << "\n" << "/* ********************** "
        << "SET OF OBJECTS ********************** */\n\n";\
  // NOTE CUMULATIVE VOLUME ASSUMES ALL 
  // MESHES HAVE SAME ORIENTATION
  //	print Integrity
  printIntegrity();
  if(good_integrity==false)
  {
    cout << "\n\nWarning: Some attributes and "
          << "characteristics were not evaluated,\n"
          << " since some mesh files failed the integrity check.\n\n";
  }
  else
  {
    //	print attributes
    printAttr();
    //	print characteristics
    adjacent_face.x_begin=0;
    area.x_begin=0;
    aspect_ratio.x_begin=0;
    edge_angle.x_begin=0;
    edge_length.x_begin=0;
    if(attr==false) { printChars(); }
  }
  cout << "/* ********************** "
        << "END ********************** */\n\n";
}

void Container::printIntegrity()
{
  cout << "\nMESH SET INTEGRITY:\n\n";
  //int a=c.countOrphan();
  if(num_orph==0){
    cout << "    # orphan vertices: none\n";
  } else {
    cout << "    # orphan vertices: " << num_orph << endl;
  }
  //a=c.countMissing();
  if(num_mis==0){
    cout << "    # missing vertices: none\n";
  } else {
    cout << "    # missin vertices: " << num_mis << endl;
  }
  //a=c.countDegen();
  if(num_deg==0){
    cout << "    # degenerate faces: none\n";
  } else {
    cout << "    # degenerate faces: " << num_deg << endl;
  }
  //a=c.countDuplV();
  if(num_dupl_v==0){
    cout << "    # duplicate vertex indices: none\n";
  } else {
    cout << "    # duplicate vertex indices: " << num_dupl_v << endl;
  }
  //a=c.countDuplF();
  if(num_dupl_f==0){
    cout << "    # duplicate face indices: none\n";
  } else {
    cout << "    # duplicate face indices: " << num_dupl_f << endl;
  }
  cout << endl;
}

void Container::printAttr()
{
  cout << "MESH SET ATTRIBUTES:\n\n";
  if(good_integrity==false){
    cout << "    Warning: These attribute summaries may be inaccurate,\n"
          << "    since some mesh files failed the integrity check.\n";
  }
  // closed
  //std::pair<int,int> pp=c.countClosed();
  if(num_clo==0){
    cout << "    # closed mesh objects: none\n";
  } else {
    cout << "    # closed mesh objects: " << num_clo << endl;
  }
  if(num_op==0){
    cout << "    # open mesh objects: none\n";
    cout << "    # border edges: none\n";
  } else {
    cout << "    # open mesh objects: " << num_op << endl;
    cout << "    # border edges: " << num_bor << endl;
  }
  // manifold
  //int val[3]={0,0,0};
  //c.countManifold(val);
  if(num_man[0]==0){
    cout << "    # manifold mesh objects: none\n";
  } else {
    cout << "    # manifold mesh objects: " << num_man[0] << endl;
  }
  /*
  if(val[2]==0){
    cout << "    # mesh objects with undefined manifoldness: none\n";
  } else {
    cout << "    # mesh objects with undefined manifoldness: " << val[2] << endl;
  }
  */
  if(num_man[1]==0 && num_man[2]==0){
    cout << "    # nonmanifold mesh objects: none\n";
    cout << "    # nonmanifold vertices: none\n";
    cout << "    # nonmanifold edges: none\n";
  } else {
    cout << "    # nonmanifold mesh objects: " << num_man[1]+num_man[2] << endl;
    cout << "    # nonmanifold vertices: " << num_nonman_v << endl;
    cout << "    # nonmanifold edges: " << num_nonman_e << endl;
  }
  // consistent
  //c.countConsistent(val);
  if(num_cons[0]==0){
    cout << "    # mesh objects with consistently "
          << "oriented face normals: none\n";
  } else {
    cout << "    # mesh objects with consistently "
          << "oriented face normals: " << num_cons[0] << endl;
  }
  if(num_cons[1]==0){
    cout << "    # mesh objects with inconsistently "
          << "oriented face normals: none\n";
  } else {
    cout << "    # mesh objects with inconsistently "
          << "oriented face normals: " << num_cons[1] << endl;
    cout << "       # flipped edges: " << num_flip << endl;
  }
  if(num_cons[2]==0){
    cout << "    # mesh objects whose face normal "
          << "orientation is undefined: none\n";
  } else {
    cout << "    # mesh objects whose face normal "
          << "orientation is undefined: " << num_cons[2] << endl;
  }
  // outward
  //c.countOutward(val);
  if(num_out[0]==0){
    cout << "    # mesh objects with outward "
          << "oriented face normals: none\n";
  } else {
    cout << "    # mesh objects with outward "
          << "oriented face normals: " << num_out[0] << endl;
  }
  if(num_out[1]==0){
    cout << "    # mesh objects with inward "
          << "oriented face normals: none\n";
  } else {
    cout << "    # mesh objects with inward "
          << "oriented face normals: " << num_out[1] << endl;
  }
  if(num_out[2]==0){
    cout << "    # mesh objects whose face normal "
          << "orientation is undefined: none\n";
  } else {
    cout << "    # mesh objects whose face normal "
          << "orientation is undefined: " << num_out[2] << endl;
  }
  cout << endl;
}

void Container::printChars()
{
  cout << "MESH SET CHARACTERISTICS:\n\n";
  //	print characteristics
  if(good_integrity==false){
    cout << "    Warning: These characteristics "
          << "summaries may be inaccurate,\n"
          << "    since some mesh files failed the integrity check.\n";
  }
  cout << "    # objects: " << num_objects << endl;
  cout << "    # vertices: " << num_v << endl;
  cout << "    # faces: " << num_f << endl;
  cout << "    # edges: " << num_e << endl;
  // TODO add this functionality
  //cout << "    # components: " << c.countComponents() << endl;
  //int val[3]={0,0,0};
  //c.countManifold(val);
  if (num_man[1]>0 || num_man[2]>0){
    cout << "    # boundaries: Since at aleast one object is nonmanifold,\n"
          << "    # boundaries: the number of boundaries may be inaccurate.\n";
  }
  //int a=c.countBoundaries();
  if(num_bou==0){
    cout << "    # boundaries: none\n";
  } else {
    cout << "    # boundaries: " << num_bou << endl;
  }
  //a=c.countIndistin();
  if(num_indistin==0){
    cout << "    # indistinguishable vertices: none\n";
  } else {
    cout << "    # indistinguishable vertices: " << num_indistin << endl;
  }
  // volume
  //c.countConsistent(val);
  if(vol>0){
    cout << "    object volume: [(data units)^3]" << endl;
    cout << "    object volume: " << vol << endl;
  } else {
    cout << "    object volume: not computed, since no mesh"
          << "    object with consistently oriented face normals was found.\n";
  }
  // TODO
  // ADD GENUS STATS (HISTOGRAM TOO)
  // TODO
  // bounding box
  cout << "    bounding box: [data units]\n";
  cout << "    bounding box: [xmin,ymin,zmin][xmax,ymax,zmax]\n";
  cout << "    bounding box: ["
        << bb[0] << ","
        << bb[1] << ","
        << bb[2] << "]["
        << bb[3] << ","
        << bb[4] << ","
        << bb[5] << "]" << endl << endl;
  // vertex adjacent faces
  adjacent_face.process();
  adjacent_face.createAdjacentFaceHistogram();
  cout << "    Vertex adjacent face statistics [faces]:" << endl;
  adjacent_face.printStats();
  cout << "    Vertex adjacent face histogram [faces]:" << endl;
  adjacent_face.printAdjacentFaceHistogram();
  cout << endl;
  // face area
  area.process();
  area.createHistogram();
  cout << "    Face area statistics [(data units)^2]:" << endl;
  cout << "       total    " << area.sum << endl;
  area.printStats();
  cout << "    Face area histogram [(data units)^2]:" << endl;
  area.printHistogram();
  cout << endl;
  // face aspect ratio
  aspect_ratio.process();
  aspect_ratio.createAspectRatioHistogram();
  cout << "    Face aspect ratio statistics [unitless]:" << endl;
  aspect_ratio.printStats();
  cout << "    Face aspect ratio histogram [unitless]:" << endl;
  aspect_ratio.printAspectRatioHistogram();
  printf("      (Aspect ratio is longest edge divided by shortest altitude)\n");
  cout << endl;
  // edge length
  edge_length.process();
  edge_length.createHistogram();
  cout << "    Edge length statistics [data units]:" << endl;
  edge_length.printStats();
  cout << "    Edge length histogram [data units]:" << endl;
  edge_length.printHistogram();
  cout << endl;
  // if edge angles computed, i.e. if at least one mesh was conistent
  if(num_cons[0]>0){
    edge_angle.process();
    edge_angle.createHistogram();
    cout << "    Edge angle statistics [degrees]:" << endl;
    edge_angle.printStats();
    cout << "    Edge angle histogram [degrees]:" << endl;
    edge_angle.printHistogram();
    cout << endl;
  }
}

int main(int argc,char **argv)
{
  // TODO store less by printing immediately if -p option
  // TODO use 'const' where sensible
  // TODO use 'std::for_each' from STL where sensible
  Container c;

  // parse command line 
  parse(argc,argv);

  // if single input file was found
  if(folder==false)
  {
    // save filename
    c.addFile(inpath);
    // build data structure of mesh
    cout << "input filename = " << inpath << endl;
    // TODO make following functions members of Object class
    // instead of relying on c.getObjectId()
    Object *o = c.processFile(inpath);
    if(o!=NULL)
    {
      if(checkIntegrity())
      {
        o->findVertexAdjacencies(c);
        //exit(0);
        createEdges(c);
        o->analyze(c);
      }
      o->print(c);
    }
  } 
  else 
  { 
    // else scan folder
    c.scanDir(inpath.c_str());
    // for each file in folder
    for (key count=0;count<c.getNumFiles();count++)
    {
      // build data structure of mesh
      Object* o=c.processFile(inpath+c.getFile(count));
      // DEBUG
      //cout << "\npostprocess:\n";cout.flush();
      //o->sanityCheck();
      //cout << "\npostprocess:\n";cout.flush();
      // DEBUG
      if(o!=NULL)
      {
        if(checkIntegrity())
        {
          //cout << "\nPROCESSING: " << c.getFile(c.getObjectId()) << endl;
          o->findVertexAdjacencies(c);
          createEdges(c);
          o->analyze(c);
          c.update();
        }
        else
        {
          c.setIntegrity(false);
        }
        o->print(c);
        //PrintObjectSum(c.getObjectId());
        c.update_x_begin();
      }
      if(true)
      {
        // clear Object data
        clearData(c);
      }
      else
      {
        c.incObjectId();
      }
    }
    // remember to set first and last to x limits
    // analyze data as set
    //cs.analyzeCumulative();
    // print cumulative surface area, volume, 
    c.printCumulative();
  }
}

