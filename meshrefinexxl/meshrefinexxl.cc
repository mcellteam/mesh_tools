#include <iostream>
#include <float.h>
#include <stxxl.h>

using std::cout;
using std::cerr;
using std::endl;


//note: 2^64=18446744073709551616
typedef long long int         key;

// bytes
#define LINE_SIZE       2048
#define FROZEN_LINE_SIZE  128
#define NODE_CACHE_SIZE (128*1024*1024) // 128MB
#define LEAF_CACHE_SIZE (128*1024*1024)
#define NUMBER_OF_BLOCKS 4
#define SUFFIX ".UPDATED"

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

typedef std::pair<key,key>                  pkk;
typedef	stxxl::map<pkk,key,pSort>           map_edge;
typedef	stxxl::map<pkk,key,pSort>::iterator m_e;

class E;
class F;
class FE;
class Object;
class V;
class VE;

struct Edata // edges
{
  key v1,v2; // relative
  key f1,f2; // relative
};

struct Fdata // faces
{
  key index;    // from file
  key v[3];	// index of vertex in context 
                // of other vertices in object
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
  key index;    // from file
  double x, y, z;
};

typedef std::vector<std::string> vec_s;
typedef std::vector<std::string>::iterator vs;
vec_s storage; // lines from frozen file
//stxxl::vector<key> frozen; // frozen vertices, file
stxxl::map<key,key,kSortComp> frozen(NODE_CACHE_SIZE,LEAF_CACHE_SIZE);
stxxl::vector<key> bisected; // bisected edges, relative
stxxl::vector<key> subdivided; // subdivided faces, relative
// edge key -> midpoint vertex key
stxxl::map<key,key,kSortComp> midpoint(NODE_CACHE_SIZE,LEAF_CACHE_SIZE);

stxxl::vector<Edata> e_vector;
stxxl::vector<Fdata> f_vector;
stxxl::vector<FEdata> fe_vector;
std::vector<Object> objects;
stxxl::vector<ObjectData> object_vector;
stxxl::vector<Vdata> v_vector;

key max_vertex = -1;   // absolute
key max_face = -1;   // absolute
key vertex_index; // relative
key face_index;   // relative
key edge_index;   // relative

class E
{
public:
  void updateEdge(key f,key object_id);
  double getSqLength(key object_id);
  bool threshold(double t);
  void print(key object_id);
  /*
  void printCP(key object_id);
  key getNewFace(key,key);
  bool isConsistent(key);
  bool getOrientation(key,key);
  double getAngle(key);
  void getVertices(key&,key&,key&,key&);
  */
  E(key i)
  {
    e_id=i; //  absolute
  }

  key get_v1() // relative
  {
    return e_vector[e_id].v1; 
  }

  key get_v2() // relative
  {
    return e_vector[e_id].v2;
  }

  key get_f1() // relative
  {
    return e_vector[e_id].f1;
  }

  key get_f2() // relative
  {
    return e_vector[e_id].f2;
  }

  key get_rel_key(key object_id) // relative
  { // absolute-absolute=relative
    return e_id-object_vector[object_id].e_start;
  }
  /*
  bool isManifold(key object_id)
  {
    // look for 3rd and higher edges stored in array
    std::pair<key,key> p;
    p = std::make_pair(get_rel_key(object_id),get_rel_key(object_id));
    return nonm_e.find(p)==nonm_e.end();
  }
  */
private:
  key e_id; // absolute
};

class F
{
public:
  key get_v1_file(key object_id); // from file
  key get_v2_file(key object_id); // from file
  key get_v3_file(key object_id); // from file
  void print(key object_id);
  /*
  key getNewEdge(key,key,key);
  void getNormal(key,double[3]);
  double getAspectRatio(key);
  */
  F(key i)
  {
    f_id=i; // absolute
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
    return f_vector[f_id].v[0];
  }

  key get_v2_rel() // relative in v_vector
  {
    return f_vector[f_id].v[1];
  }

  key get_v3_rel() // relative in v_vector
  {
    return f_vector[f_id].v[2];
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
    return fe_vector[fe_id].e[i]; // edge key, relative
  }

private:
  key fe_id; // absolute
};

class Object
{
public:
  /*
  bool isClosed();
  void evalAttributes(Container&);
  void evalCharacteristics(Container&);
  bool verticesManifold(bool);
  bool edgesManifold(bool);
  bool isManifold();
  bool isConsistent();
  void vertexDistin();
  void boundObject();
  void vertexAdjacentFaces(Container&);
  void areaAspectRatio(Container&);
  void processEdgeLengths(Container&);
  void computeEdgeAngles(Container&);
  void computeVolume();
  void computeGenus();
  void print(Container&);
  void printAttr(Container&);
  void printIntegrity(Container&);
  void printChars(Container&);
  void analyze(Container&);
  void process_edge_lengths(key e,double l);

  void clearStats()
  {
    adjacent_face.clear();
    area.clear();
    aspect_ratio.clear();
    edge_angle.clear();
    edge_length.clear();
  }
  */
  Object (key i,std::string filename)
  {
    object_id = i;
    name = filename;
    bb[0]=bb[1]=bb[2]=1E30;
    bb[3]=bb[4]=bb[5]=-1E30;
  }
  /*
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
  */

  E *GetE(key id); // relative
  /*
  NE *GetNE() // relative
  {
    return new NE();
  }
  */
  FE *GetFE(key id) // relative
  {
    key fe_id = object_vector[object_id].fe_start + id; // absolute

    return new FE(fe_id);
  }
  /*
  VE *GetVE(key id); // relative
  */
  V *GetV(key id); // relative

  F *GetF(key id) // relative
  {
    key f_id = object_vector[object_id].f_start + id; // absolute

    return new F(f_id);
  }

  std::string GetName()
  {
    return name;
  }

  /*
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

  vd area_begin()
  {
    return area.x.begin();
  }
  vd area_end()
  {
    return area.x.end();
  }
  vd aspect_ratio_begin()
  {
    return aspect_ratio.x.begin();
  }
  vd aspect_ratio_end()
  {
    return aspect_ratio.x.end();
  }
  vd edge_length_begin()
  {
    return edge_length.x.begin();
  }
  vd edge_length_end()
  {
    return edge_length.x.end();
  }
  vd edge_angle_begin()
  {
    return edge_angle.x.begin();
  }
  vd edge_angle_end()
  {
    return edge_angle.x.end();
  }
  vd adjacent_face_begin()
  {
    return adjacent_face.x.begin();
  }
  vd adjacent_face_end()
  {
    return adjacent_face.x.end();
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

  void process_adjacent_face (int n)
  {
    //cout << "adjacent_face.sum = " << adjacent_face.sum << endl;
    adjacent_face.sum+=n;
    adjacent_face.sum2+=n*n;
    adjacent_face.total+=n;
    if(n<adjacent_face.min) {adjacent_face.min=n;}
    if(n>adjacent_face.max) {adjacent_face.max=n;}
    adjacent_face.x.push_back(n);
  }

  void process_area (double n)
  {
    area.sum+=n;
    area.sum2+=n*n;
    area.total+=n;
    if(n<area.min) {area.min=n;}
    if(n>area.max) {area.max=n;}
    area.x.push_back(n);
  }

  void process_aspect_ratio (double n)
  {
    aspect_ratio.sum+=n;
    aspect_ratio.sum2+=n*n;
    aspect_ratio.total+=n;
    if(n<aspect_ratio.min) {aspect_ratio.min=n;}
    if(n>aspect_ratio.max) {aspect_ratio.max=n;}
    aspect_ratio.x.push_back(n);
    // TODO
    // compare face aspect ratio to user-defined threshold
    // if(ar>cs.thresholds[0]) bad_aspect[*i]=ar;
  }

  void process_edge_angle (double angle)
  {
    edge_angle.sum+=angle;
    edge_angle.sum2+=angle*angle;
    edge_angle.total+=angle;
    if(angle<edge_angle.min) edge_angle.min=angle;
    if(angle>edge_angle.max) edge_angle.max=angle;
    // add to vector
    edge_angle.x.push_back(angle);
  }

  void sanityCheck()
  {
     cout << "sanityCheck:adjacent_face.x.size() = "
           << adjacent_face.x.size() << endl;cout.flush();
     if(adjacent_face.x.begin()!=adjacent_face.x.end())
     {
       cout << "sanityCheck:not equal\n";
     }
     else if (adjacent_face.x.begin()==adjacent_face.x.end())
     {
       cout << "sanityCheck:equal\n";
     }
  }
  */

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

class V
{
public:
  void print(key object_id);
  /*
  bool isManifold(key,bool);
  void getAdjacentEdges(key,std::vector<key>&);
  bool scanAdjFaces(key,key,key,bool&);
  void printCP();
  */
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
    return v_id-object_vector[object_id].v_start;
  }

  key get_abs_key()
  {
    return v_id;
  }

private:
  key v_id; // absolute
};

E* Object::GetE(key id) // relative
{
  key e_id = object_vector[object_id].e_start + id; // absolute
  return new E(e_id);
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
  F* f_ptr = objects[object_id].GetF(get_f1());
  cout  << "Face " << f_ptr->get_index()
        << " " << f_ptr->get_v1_file(object_id)
        << " " << f_ptr->get_v2_file(object_id)
        << " " << f_ptr->get_v3_file(object_id)
        << endl;
  delete f_ptr;
  f_ptr = objects[object_id].GetF(get_f2());
  cout  << "Face " << f_ptr->get_index()
        << " " << f_ptr->get_v1_file(object_id)
        << " " << f_ptr->get_v2_file(object_id)
        << " " << f_ptr->get_v3_file(object_id)
        << endl;
  delete f_ptr;
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
    cout << "\nE::updateEdge: Error: edge is nonmanifold.\n"
          << "f key = " << f << endl;
          print(object_id);
    exit(0);
    // save 3rd or greater edge as edge,face pair
    //NE* ne_ptr = objects[object_id].GetNE();
    //ne_ptr->addPair(get_rel_key(object_id),f);
    //delete ne_ptr;
  }
  // add edge pointer to face
  FE* fe_ptr = objects[object_id].GetFE(f);
  fe_ptr->addEdge(get_rel_key(object_id));
  delete fe_ptr;
}

V *Object::GetV(key id) // relative
{
  key v_id = object_vector[object_id].v_start + id; // absolute
  return new V(v_id);
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

bool E::threshold(double t)
{
  double l = getSqLength(0);
  if (l>t)
  {
    return true;
  }
  else
  {
    return false;
  }
}

void F::print(key object_id)
{
  cout.precision(12);
  cout << "Face " << get_index()
        << " " << get_v1_file(object_id)
        << " " << get_v2_file(object_id)
        << " " << get_v3_file(object_id)
        << endl;
}

void V::print(key object_id)
{
  cout.precision(12);
  cout << "Vertex " << get_index()
        << " " << get_x()
        << " " << get_y()
        << " " << get_z()
        << endl;
}

void printMesh()
{
  key oid = 0;
  // for each vertex in object
  int num=object_vector[oid].v_end
        -object_vector[oid].v_start+1;
  for(int i=0;i<num;i++) // relative
  {
    // access vertex file index
    V* v_ptr = objects[oid].GetV(i);
    v_ptr->print(oid);
    delete v_ptr;
  }
  // for each face in object
  num=object_vector[oid].f_end
        -object_vector[oid].f_start+1;
  for(int i=0;i<num;i++) // relative
  {
    // screen psuedo-erased faces, i.e. inactive
    F* f_ptr = objects[oid].GetF(i);
    bool good = (f_ptr->get_index()>0);
    delete f_ptr;
    if(good==true)
    {
      // face is active
      // access vertex file index
      F* f_ptr = objects[oid].GetF(i);
      f_ptr->print(oid);
      delete f_ptr;
    }
  }
}

void printFrozen(std::string filename)
{
  key oid = 0;
  char fileout[LINE_SIZE];
  // open output data file
  sprintf(fileout,"%s%s",filename.c_str(),SUFFIX);
  std::ofstream outFile(fileout);
  if (outFile.fail()) // if stream cannot be opened
  {
    cout << "\n\nprintFrozen: Can not open output file ["
          << fileout << "]\n\n";
    exit(1); 
  }
  // print storage
  for(vs i=storage.begin();i!=storage.end();i++)
  {
    outFile << *i;
  }
  // for each vertex in object
  stxxl::map<key,key,kSortComp>::iterator m_f;
  for(m_f=frozen.begin();m_f!=frozen.end();m_f++)
  {
    // if vertex if frozen
    if((*m_f).second==true)
    {
      // print
      outFile << objects[oid].GetName()
             << " " << (*m_f).first << endl;
    }
  }
  outFile.close();
}

struct SortByKey // comparison function
{
  bool operator () (const key & a,
                    const key & b) const
  {
    return  a < b;
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

void validateThresholding()
{
  key oid = 0;
  // for each face in object
  int num=object_vector[oid].f_end
        -object_vector[oid].f_start+1;
  for(int i=0;i<num;i++) // relative
  {
    // screen psuedo-erased faces, i.e. inactive
    F* f_ptr = objects[oid].GetF(i);
    bool good = (f_ptr->get_index()>0);
    delete f_ptr;
    if(good==true)
    {
      // face is active
      F* f_ptr = objects[oid].GetF(i);
      FE* fe_ptr = objects[oid].GetFE(f_ptr->get_rel_key(oid));
      delete f_ptr;
      key edges[3] = {fe_ptr->get_edge(0),
        fe_ptr->get_edge(1),
        fe_ptr->get_edge(2)};
      delete fe_ptr;
      // if any edge is bisected, then face must be subdivided
      if((bisected[edges[0]]>0) ||
         (bisected[edges[1]]>0) ||
         (bisected[edges[2]]>0))
      {
        // then face must be subdivided
        // so if face is not subdivided
        if(subdivided[i]<0)
        {
          cout << "\nvalidateThresholding: Error."
                << "Inconsistent thresholding.\n"
                << "bisected=["
                << bisected[edges[0]] << " "
                << bisected[edges[1]] << " "
                << bisected[edges[2]] << "]"
                << ", subdivided = " << subdivided[i] << endl;
          exit(0);
        }
      }
      else
      {
        // no edge is bisected
        // then face must not be subdivided
        // so if face is subdivided
        if(subdivided[i]>0)
        {
          cout << "\nvalidateThresholding: Error."
                << "Inconsistent thresholding.\n"
                << "bisected=["
                << bisected[edges[0]] << " "
                << bisected[edges[1]] << " "
                << bisected[edges[2]] << "]"
                << ", subdivided = " << subdivided[i] << endl;
          exit(0);
        }
      }
    }
  }
}

bool thresholdEdges(double t,int &n)
{
  bisected.clear();
  subdivided.clear();
  key oid = 0;
  n=0;
  bool flag = false;
  // prep subdivided
  // for each face in object
  int num=object_vector[oid].f_end
        -object_vector[oid].f_start+1;
  for(int i=0;i<num;i++) // relative
  {
    subdivided.push_back(-2);
  }
  // prep bisected
  // for each edge in object
  num=object_vector[oid].e_end
        -object_vector[oid].e_start+1;
  for(int i=0;i<num;i++) // relative
  {
    bisected.push_back(-2);
  }
  // for each face in object
  num=object_vector[oid].f_end
        -object_vector[oid].f_start+1;
  for(int i=0;i<num;i++) // relative
  {
    F* f_ptr = objects[oid].GetF(i);
    bool good = (f_ptr->get_index()>0);
    delete f_ptr;
    if(good==true)
    {
      F* f_ptr = objects[oid].GetF(i);
      FE* fe_ptr = objects[oid].GetFE(f_ptr->get_rel_key(oid));
      delete f_ptr;
      key edges[3] = {fe_ptr->get_edge(0),
                      fe_ptr->get_edge(1),
                      fe_ptr->get_edge(2)};
      delete fe_ptr;
      // for each edge of face
      for(int j=0;j<3;j++)
      {
        E* e_ptr = objects[oid].GetE(edges[j]);
        bool too_long = e_ptr->threshold(t);
        delete e_ptr;
        // threshold criteria
        if(too_long==true)
        {
          // mark face as subdivided
          //subdivided.push_back(i);
          subdivided[i] = 1;
          // mark edge as bisected
          //bisected.push_back(triplet[j]);
          bisected[edges[j]] = 1;
          // set flag
          flag = true;
          n++;
        }
      }
    }
  }
  // bound the main memory consumption by M during sorting
  //const unsigned M = 1024*1024*1024; // bytes
  //stxxl::sort(bisected.begin(),bisected.end(),SortByKey(),M);
  //stxxl::sort(subdivided.begin(),subdivided.end(),SortByKey(),M);
  return flag;
}

int createNewVertices(void)
{
  int n=0;
  key oid = 0;
  // for each edge in object
  int num=object_vector[oid].e_end
        -object_vector[oid].e_start+1;
  //cout << "\ncreateNewVertices: #edges = " << num << endl;
  for(int i=0;i<num;i++) // relative
  {
    //cout << "createNewVertices: edge " << i;
    // if edge is bisected, i.e. edge key found in bisected vector
    //stxxl::vector<key>::iterator j;
    //j = stxxl::find(bisected.begin(),bisected.end(),i,NUMBER_OF_BLOCKS);
    // if match found
    //if(j!=bisected.end())
    if(bisected[i]>0)
    {
      // new vertex index
      key newi = ++max_vertex;
      //cout << " is bisected\n";
      E* e_ptr = objects[oid].GetE(i);
      V *v1_ptr = objects[oid].GetV(e_ptr->get_v1());
      V *v2_ptr = objects[oid].GetV(e_ptr->get_v2());
      delete e_ptr;
      // compute new vertex
      double x=(v1_ptr->get_x()+v2_ptr->get_x())/2.0;
      double y=(v1_ptr->get_y()+v2_ptr->get_y())/2.0;
      double z=(v1_ptr->get_z()+v2_ptr->get_z())/2.0;
      // determine if frozen
      if(frozen[v1_ptr->get_index()]==true &&
         frozen[v2_ptr->get_index()]==true)
      {
        frozen[newi]=true;
      }
      else
      {
        frozen[newi]=false;
      }
      delete v1_ptr;
      delete v2_ptr;
      // create new vertex
      Vdata vd;
      vd.index = newi;
      vd.x=x;
      vd.y=y;
      vd.z=z;
      // save vertex in edge
      midpoint[i]=vertex_index++;
      // save vertex in object
      v_vector.push_back(vd);
      object_vector[oid].v_end++; // absolute
      // count new vertices
      n++;
    }
    //else
    //{
    //  cout << " is not bisected\n";
    //}
  }
  return n;
}

int countBisected(key f,key mid_verts[3])
{
  // collect face edge keys
  key oid = 0;
  FE* fe_ptr = objects[oid].GetFE(f);
  key edges[3] = {fe_ptr->get_edge(0),
                    fe_ptr->get_edge(1),
                    fe_ptr->get_edge(2)};
  delete fe_ptr;
  /*
  cout << "\ncountBisected: "
        << "edges ["
        << edges[0] << " "
        << edges[1] << " "
        << edges[2] << "]"
        << ", bisected ["
        << bisected[0] << " "
        << bisected[1] << " "
        << bisected[2] << "]\n";
  */
  int count=0;
  // for each face edge
  for(int i=0;i<3;i++)
  {
    // if edge is bisected, i.e. edge key found in bisected vector
    if(bisected[edges[i]]>0)
    {
      // add new vertex key, relative
      // NOT edge key !!!!
      mid_verts[i]=midpoint[edges[i]];
      count++;
    }
    else
    {
      mid_verts[i]=-2;
    }
  }
  return count;
}

bool addFace(stxxl::vector<Fdata> &nf,int new_index,key va,key vb,key vc)
{
  Fdata fd;
  fd.index = new_index;
  fd.v[0] = va;
  fd.v[1] = vb;
  fd.v[2] = vc;
  nf.push_back(fd);
  if(va==vc)
  {
    cout << "\naddFace: Error. duplicate indices.\n"
          << "Face " << new_index
          << " " << va
          << " " << vb
          << " " << vc << endl;
    return true;
    //exit(0);
  }
  return false;
}

double dot(double a[3],double b[3])
{
  return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
}

double getAspectRatio(key v_a,key v_b,key v_c)
{
  //cout << "\ngetAspectRatio: v_a="<<v_a
  //  << ", v_b="<<v_b
  //  << ", v_c="<<v_c
  //  << ", v_vector.size="<<v_vector.size()
  //  << ", a";cout.flush();
  key oid = 0;
  V *v0_ptr = objects[oid].GetV(v_a);
  V *v1_ptr = objects[oid].GetV(v_b);
  V *v2_ptr = objects[oid].GetV(v_c);
  //cout << "b1";cout.flush();
  double v0[3]={v0_ptr->get_x(),v0_ptr->get_y(),v0_ptr->get_z()};
  //cout << "b2";cout.flush();
  double v1[3]={v1_ptr->get_x(),v1_ptr->get_y(),v1_ptr->get_z()};
  //cout << "b3";cout.flush();
  double v2[3]={v2_ptr->get_x(),v2_ptr->get_y(),v2_ptr->get_z()};
  //cout << "c";cout.flush();
  delete v0_ptr;
  delete v1_ptr;
  delete v2_ptr;
  //cout << "d";cout.flush();
  // Make triangle edge vectors
  double va[3]={v1[0]-v0[0],v1[1]-v0[1],v1[2]-v0[2]};
  double vb[3]={v2[0]-v1[0],v2[1]-v1[1],v2[2]-v1[2]};
  double vc[3]={v0[0]-v2[0],v0[1]-v2[1],v0[2]-v2[2]};
  //cout << "e";cout.flush();
  double vbase[3]={0,0,0};
  double vopp[3]={0,0,0};
  // Find length of longest edge
  double lmax=-DBL_MAX;
  double la=sqrt(dot(va,va));
  double lb=sqrt(dot(vb,vb));
  double lc=sqrt(dot(vc,vc));
  //cout << "f";cout.flush();
  if(la>lmax)
  {
    lmax=la;
    vbase[0]=va[0];
    vbase[1]=va[1];
    vbase[2]=va[2];
    vc[0]=v2[0]-v0[0];
    vc[1]=v2[1]-v0[1];
    vc[2]=v2[2]-v0[2];
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
    va[0]=v0[0]-v1[0];
    va[1]=v0[1]-v1[1];
    va[2]=v0[2]-v1[2];
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
    vb[0]=v1[0]-v2[0];
    vb[1]=v1[1]-v2[1];
    vb[2]=v1[2]-v2[2];
    vopp[0]=vb[0];
    vopp[1]=vb[1];
    vopp[2]=vb[2];
  }
  // Find shortest altitude
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

double getMaxAR(key v0,key v1,key v2,key v3)
{
  // get max aspect ratio connecting v0 and v2
  // get aspect ratio of face v0,v1,v2
  double ar012 = getAspectRatio(v0,v1,v2);
  // get aspect ratio of face v0,v2,v3
  double ar023 = getAspectRatio(v0,v2,v3);
  // identify max aspect ratios
  if( ar012 > ar023 ){ return ar012;}
  else { return ar023;}
}

void validateVertices()
{
  key oid = 0;
  key mid_verts[3];
  // for each face in object
  int num=object_vector[oid].f_end
        -object_vector[oid].f_start+1;
  for(int i=0;i<num;i++) // relative
  {
    // screen psuedo-erased faces, i.e. inactive
    F* f_ptr = objects[oid].GetF(i);
    bool good = (f_ptr->get_index()>0);
    delete f_ptr;
    if(good==true)
    {
      // face is active
      // if face is subdivided, i.e. face key found in subdivided vector
      if(subdivided[i]>0)
      {
        // identify new vertices
        countBisected(i,mid_verts);
        key f_verts[3];
        F* f_ptr = objects[oid].GetF(i);
        f_verts[0] = f_ptr->get_v1_rel();
        f_verts[1] = f_ptr->get_v2_rel();
        f_verts[2] = f_ptr->get_v3_rel();
        delete f_ptr;
        key maxval = v_vector.size()-1;
        for(int j=0;j<3;j++)
        {
          if(mid_verts[j]>maxval || f_verts[j]>maxval)
          {
            cout << "\nvalidateVertices: Error."
                  << "vertex key too large.\n"
                  << "v_vector.size=" << v_vector.size()
                  << ", mid_verts=["
                  << mid_verts[0] << " "
                  << mid_verts[1] << " "
                  << mid_verts[2] << "]"
                  << ", f_verts=["
                  << f_verts[0] << " "
                  << f_verts[1] << " "
                  << f_verts[2] << "]"
                  << endl;
            exit(0);
          }
        }
      }
    }
  }
  cout << "\nVertices are valid.\n";cout.flush();
}

int createNewSubdividedFaces()
{
  // DEBUG
  bool xxx = false;
  // DEBUG
  //Vertex* mid_verts[3];
  key mid_verts[3];
  //std::vector<Face*> nf; // new faces
  stxxl::vector<Fdata> nf; // new faces
  key oid=0;
  // for each face in object
  key num=object_vector[oid].f_end
        -object_vector[oid].f_start+1;
  for(int i=0;i<num;i++) // relative
  {
    F* f_ptr = objects[oid].GetF(i);
    bool good = (f_ptr->get_index()>0);
    delete f_ptr;
    if(good==true)
    {
      // if face is subdivided, i.e. face key found in subdivided vector
      //stxxl::vector<key>::iterator j;
      //j = stxxl::find(subdivided.begin(),subdivided.end(),i,NUMBER_OF_BLOCKS);
      // if match found
      //if(j!=subdivided.end())
      if(subdivided[i]>0)
      {
        // identify new vertices
        int num_subdivided_edges = countBisected(i,mid_verts);
        ///// if one edge subdivided /////
        if(num_subdivided_edges==1)
        {
          F* f_ptr = objects[oid].GetF(i);
          // add new faces
          if(mid_verts[0]>0)
          {
            xxx=addFace(nf,++max_face,f_ptr->get_v1_rel(),mid_verts[0],f_ptr->get_v3_rel());
            if(xxx==true){cout << "\ncreateNewSubdividedFaces: a.\n";exit(0);}
            xxx=addFace(nf,++max_face,f_ptr->get_v3_rel(),mid_verts[0],f_ptr->get_v2_rel());
            if(xxx==true){cout << "\ncreateNewSubdividedFaces: b.\n";exit(0);}
          }
          else if(mid_verts[1]>0)
          {
            xxx=addFace(nf,++max_face,f_ptr->get_v2_rel(),mid_verts[1],f_ptr->get_v1_rel());
            if(xxx==true){cout << "\ncreateNewSubdividedFaces: c.\n";exit(0);}
            xxx=addFace(nf,++max_face,f_ptr->get_v1_rel(),mid_verts[1],f_ptr->get_v3_rel());
            if(xxx==true){cout << "\ncreateNewSubdividedFaces: d.\n";exit(0);}
          }
          else if(mid_verts[2]>0)
          {
            xxx=addFace(nf,++max_face,f_ptr->get_v3_rel(),mid_verts[2],f_ptr->get_v2_rel());
            if(xxx==true){cout << "\ncreateNewSubdividedFaces: e.\n";exit(0);}
            xxx=addFace(nf,++max_face,f_ptr->get_v2_rel(),mid_verts[2],f_ptr->get_v1_rel());
            if(xxx==true){cout << "\ncreateNewSubdividedFaces: f.\n";exit(0);}
          }
          delete f_ptr;
        }
        else if(num_subdivided_edges==2)
        {
          key f_verts[3];
          F* f_ptr = objects[oid].GetF(i);
          f_verts[0] = f_ptr->get_v1_rel();
          f_verts[1] = f_ptr->get_v2_rel();
          f_verts[2] = f_ptr->get_v3_rel();
          delete f_ptr;
          if(mid_verts[0]<0)
          {
            // get max aspect ratio connecting mid_verts[2] and f->v2
            double max1 = getMaxAR(mid_verts[2],f_verts[0],f_verts[1],mid_verts[1]);
            // get max aspect ratio connecting mid_verts[1] and f->v1
            double max2 = getMaxAR(mid_verts[1],mid_verts[2],f_verts[0],f_verts[1]);
            // if connecting mid_verts[2],f->v2 generates
            // faces with smallest aspect ratio
            if(max1<max2)
            {
              // add new faces
              xxx=addFace(nf,++max_face,f_verts[2],mid_verts[2],mid_verts[1]);
              if(xxx==true){cout << "\ncreateNewSubdividedFaces: g.\n";exit(0);}
              xxx=addFace(nf,++max_face,mid_verts[2],f_verts[0],f_verts[1]);
              if(xxx==true){cout << "\ncreateNewSubdividedFaces: h.\n";exit(0);}
              xxx=addFace(nf,++max_face,mid_verts[2],f_verts[1],mid_verts[1]);
              if(xxx==true){cout << "\ncreateNewSubdividedFaces: i.\n";exit(0);}
            }
            // else connecting mid_verts[1],f->v1 generates
            // faces with smallest aspect ratio
            else
            {
              // add new faces
              xxx=addFace(nf,++max_face,f_verts[2],mid_verts[2],mid_verts[1]);
              if(xxx==true){cout << "\ncreateNewSubdividedFaces: j.\n";exit(0);}
              xxx=addFace(nf,++max_face,mid_verts[1],f_verts[0],f_verts[1]);
              if(xxx==true){cout << "\ncreateNewSubdividedFaces: k.\n";exit(0);}
              xxx=addFace(nf,++max_face,mid_verts[2],f_verts[0],mid_verts[1]);
              if(xxx==true){cout << "\ncreateNewSubdividedFaces: l.\n";exit(0);}
            }
          }
          else if(mid_verts[1]<0)
          {
            double max1 = getMaxAR(mid_verts[0],f_verts[1],f_verts[2],mid_verts[2]);
            double max2 = getMaxAR(mid_verts[2],mid_verts[0],f_verts[1],f_verts[2]);
            // if connecting mid_verts[0],v3 generates
            // faces with smallest aspect ratio
            if(max1<max2)
            {
              // add new faces
              xxx=addFace(nf,++max_face,f_verts[0],mid_verts[0],mid_verts[2]);
              if(xxx==true){cout << "\ncreateNewSubdividedFaces: m.\n";exit(0);}
              xxx=addFace(nf,++max_face,mid_verts[0],f_verts[1],f_verts[2]);
              if(xxx==true){cout << "\ncreateNewSubdividedFaces: n.\n";exit(0);}
              xxx=addFace(nf,++max_face,mid_verts[0],f_verts[2],mid_verts[2]);
              if(xxx==true){cout << "\ncreateNewSubdividedFaces: o.\n";exit(0);}
            }
            // else connecting mid_verts[2],v2 generates
            // faces with smallest aspect ratio
            else
            {
              // add new faces
              xxx=addFace(nf,++max_face,f_verts[0],mid_verts[0],mid_verts[2]);
              if(xxx==true){cout << "\ncreateNewSubdividedFaces: p.\n";exit(0);}
              xxx=addFace(nf,++max_face,mid_verts[0],f_verts[1],mid_verts[2]);
              if(xxx==true){cout << "\ncreateNewSubdividedFaces: q.\n";exit(0);}
              xxx=addFace(nf,++max_face,f_verts[1],f_verts[2],mid_verts[2]);
              if(xxx==true){cout << "\ncreateNewSubdividedFaces: s.\n";exit(0);}
            }
          }
          else if(mid_verts[2]<0)
          {
            double max1 = getMaxAR(mid_verts[1],f_verts[2],f_verts[0],mid_verts[0]);
            double max2 = getMaxAR(mid_verts[0],mid_verts[1],f_verts[2],f_verts[0]);
            // if connecting mid_verts[1],v1 generates
            // faces with smallest aspect ratio
            if(max1<max2)
            {
              // add new faces
              xxx=addFace(nf,++max_face,f_verts[1],mid_verts[1],mid_verts[0]);
              if(xxx==true){cout << "\ncreateNewSubdividedFaces: t.\n";exit(0);}
              xxx=addFace(nf,++max_face,mid_verts[1],f_verts[2],f_verts[0]);
              if(xxx==true){cout << "\ncreateNewSubdividedFaces: u.\n";exit(0);}
              xxx=addFace(nf,++max_face,mid_verts[1],f_verts[0],mid_verts[0]);
              if(xxx==true){cout << "\ncreateNewSubdividedFaces: v.\n";exit(0);}
            }
            // else connecting mid_verts[0],v3 generates
            // faces with smallest aspect ratio
            else
            {
              // add new faces
              xxx=addFace(nf,++max_face,f_verts[1],mid_verts[1],mid_verts[0]);
              if(xxx==true){cout << "\ncreateNewSubdividedFaces: w.\n";exit(0);}
              xxx=addFace(nf,++max_face,mid_verts[1],f_verts[2],mid_verts[0]);
              if(xxx==true)
              {
                cout << "\ncreateNewSubdividedFaces: x.\n"
                      << "nf.size=" << nf.size()
                      << ", max_face=" << max_face
                      << ", mid_verts=["
                      << mid_verts[0] << " "
                      << mid_verts[1] << " "
                      << mid_verts[2] << "]"
                      <<"\nf_verts=[" 
                      << f_verts[0] << " "
                      << f_verts[1] << " "
                      << f_verts[2] << "]\n\n";
                exit(0);
              }
              xxx=addFace(nf,++max_face,f_verts[2],f_verts[0],mid_verts[0]);
              if(xxx==true){cout << "\ncreateNewSubdividedFaces: y.\n";exit(0);}
            }
          }
        }
        else if(num_subdivided_edges==3)
        {
          key f_verts[3];
          F* f_ptr = objects[oid].GetF(i);
          f_verts[0] = f_ptr->get_v1_rel();
          f_verts[1] = f_ptr->get_v2_rel();
          f_verts[2] = f_ptr->get_v3_rel();
          delete f_ptr;
          // add new faces
          xxx=addFace(nf,++max_face,f_verts[0],mid_verts[0],mid_verts[2]);
          if(xxx==true){cout << "\ncreateNewSubdividedFaces: z.\n";exit(0);}
          xxx=addFace(nf,++max_face,mid_verts[0],f_verts[1],mid_verts[1]);
          if(xxx==true){cout << "\ncreateNewSubdividedFaces: aa.\n";exit(0);}
          xxx=addFace(nf,++max_face,mid_verts[0],mid_verts[1],mid_verts[2]);
          if(xxx==true){cout << "\ncreateNewSubdividedFaces: bb.\n";exit(0);}
          xxx=addFace(nf,++max_face,mid_verts[2],mid_verts[1],f_verts[2]);
          if(xxx==true){cout << "\ncreateNewSubdividedFaces: cc.\n";exit(0);}
        }
        else
        {
          printf("\n\ncreateNewSubdividedFaces: Error: Face is subdivided, ");
          printf("but num_subdivided_edges = %d\n\n",num_subdivided_edges);
          exit(0);
        }
        // mark face as erased
        f_vector[i].index=-1;
      }
    }
  }
  // for each new face
  stxxl::vector<Fdata>::iterator i;
  for(i=nf.begin();i!=nf.end();i++)
  {
    // add new faces to old faces
    f_vector.push_back(*i);
    face_index++; // relative
    object_vector[oid].f_end++; // absolute
  }
  return nf.size();
}

bool refineMesh(double t)
{
  bool flag = true;
  fprintf(stderr,"Thresholding edges..............");
  fflush(stderr);
  int n;
  flag = thresholdEdges(t,n);
  cerr << "complete. " << n
        << " divided edges of " << edge_index << endl;
  fflush(stderr);
  //validateThresholding();
  if(flag==true)
  {
    fprintf(stderr,"Creating new vertices...........");
    fflush(stderr);
    n = createNewVertices();
    fprintf(stderr,"complete. %d new verts.\n",n);
    fflush(stderr);
    //validateVertices();
    //validateThresholding();
    fprintf(stderr,"Creating new faces..............");
    fflush(stderr);
    n = createNewSubdividedFaces();
    //validateThresholding();
    //validateVertices();
    fprintf(stderr,"complete. %d new faces.\n",n);
    fflush(stderr);
  }
  return flag;
}

void ReadVertex(char* triplet,Vdata &vd)
{
  char val[80];
  char *eptr;
  int i;

  char *cp=triplet;

  // get past 'Vertex'
  while (strchr("Vertx",*triplet)!=NULL) {triplet++;}

  // grab vertex index
  while (strchr(" \t,",*triplet)!=NULL) { triplet++; }
  i=0;
  while (strchr("0123456789+-eE.",*triplet)!=NULL)
  {
    val[i++] = *triplet++;
  }
  val[i]=0;
  // TODO find conversion which supports long long int
  key v = static_cast<key>(strtod(val,&eptr));
  if (val==eptr)
  {
    vd.index=0;
    vd.x=vd.y=vd.z=0;
    printf("Error in reading vertex index\n");
    return;
  }
  if(v > std::numeric_limits<key>::max())
  {
    cout << "\nError: vertex index = "
          << v << " is larger than max allowable integer value = "
          << std::numeric_limits<key>::max() << endl;
    exit(0);
  }
  vd.index = v;
  if(v>max_vertex){max_vertex=v;}

  // grab x coord
  while (strchr(" \t,",*triplet)!=NULL) { triplet++; }
  i=0;
  while (strchr("0123456789+-eE.",*triplet)!=NULL)
  {
    val[i++] = *triplet++;
  }
  val[i]=0;
  vd.x=strtod(val,&eptr);
  if (val==eptr)
  {
    vd.x=vd.y=vd.z=0;
    printf("Error in reading vertex\n");
    printf("Error in reading vertex: string %s\n",cp);
    return;
  }

  // grab y coord
  while (strchr(" \t,",*triplet)!=NULL) triplet++;
  i=0;
  while (strchr("0123456789+-eE.",*triplet))
  {
    val[i++] = *triplet++;
  }
  val[i]=0;
  vd.y=strtod(val,&eptr);
  if (val==eptr)
  {
    vd.x=vd.y=vd.z=0;
    printf("Error in reading vertex\n");
    printf("Error in reading vertex: string %s\n",cp);
    return;
  }

  // grab z coord
  while (strchr(" \t,",*triplet)!=NULL) triplet++;
  i=0;
  while (strchr("0123456789+-eE.",*triplet))
  {
    val[i++] = *triplet++;
  }
  val[i]=0;
  vd.z=strtod(val,&eptr);
  if (val==eptr)
  {
    vd.x=vd.y=vd.z=0;
    printf("Error in reading vertex\n");
    printf("Error in reading vertex: string %s\n",cp);
    return;
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
  char val[80];
  char *eptr;
  int i;

  // get past 'Face'
  while (strchr("Face",*triplet)!=NULL) {triplet++;}

  /////////// grab Face index ///////////////////
  while (strchr(" \t,",*triplet)!=NULL) { triplet++; }
  i=0;
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
    fd.index=0;
    fd.v[0]=fd.v[1]=fd.v[2]=0;
    printf("Error in reading face index\n");
    return;
  }
  // check for index > max int
  if(v > std::numeric_limits<key>::max())
  {
    cout << "\nError: face index = "
          << v << " is larger than max allowable integer value = "
          << std::numeric_limits<key>::max() << endl;
    exit(0);
  }
  fd.index = v;
  if(v>max_face){max_face=v;}

  //////////// grab first vertex index ////////////////
  while (strchr(" \t,",*triplet)!=NULL) { triplet++; }
  i=0;
  while (strchr("0123456789+-eE.",*triplet)!=NULL)
  {
    val[i++] = *triplet++;
  }
  val[i]=0;
  // TODO find conversion which supports long long int
  v = static_cast<key>(strtod(val,&eptr));
  // check for string length == 0
  if (val==eptr)
  {
    fd.v[0]=fd.v[1]=fd.v[2]=0;
    printf("Error in reading vertex index\n");
    return;
  }
  // check for index > max int
  if(v > std::numeric_limits<key>::max())
  {
    cout << "\nError: first vertex index = "
          << v << " is larger than max allowable integer value = "
          << std::numeric_limits<key>::max() << endl;
    exit(0);
  }
  // search for matching vertex index in map
  // map<vertex index from file, relative vertex index>
  stxxl::map<key,key,kSortComp>::const_iterator p=vp.find(v);
  // if match found
  if(p!=vp.end())
  {
    // set Face vertex key, first
    fd.v[0] = (*p).second;
  }
  else
  {
    // print details
    cout << "\nError: Face " << fd.index
          << " references missing vertex " 
          << v << endl << line;
          exit(0);
  }

  ////////////// grab second vertex index ///////////////
  while (strchr(" \t,",*triplet)!=NULL) triplet++;
  i=0;
  while (strchr("0123456789+-eE.",*triplet))
  {
    val[i++] = *triplet++;
  }
  val[i]=0;
  // TODO find conversion which supports long long int
  v = static_cast<key>(strtod(val,&eptr));
  // check for string length == 0
  if (val==eptr)
  {
    fd.v[0]=fd.v[1]=fd.v[2]=0;
    printf("Error in reading vertex index\n");
    return;
  }
  // check for index > max int
  if(v > std::numeric_limits<key>::max())
  {
    cout << "\nError: second vertex index = "
          << v << " is larger than max allowable integer value = "
          << std::numeric_limits<key>::max() << endl;
    exit(0);
  }
  // search for matching vertex index in map
  // map<vertex index from file, relative vertex index>
 p=vp.find(v);
  // if match found
  if(p!=vp.end())
  {
    // set Face vertex key, second
    fd.v[1] = (*p).second;
  }
  else
  {
    // print details
    cout << "\nError: Face " << fd.index
          << " references missing vertex " 
          << v << endl << line;
          exit(0);
  }

  //////////// grab third vertex index //////////////
  while (strchr(" \t,",*triplet)!=NULL) triplet++;
  i=0;
  while (strchr("0123456789+-eE.",*triplet))
  {
    val[i++] = *triplet++;
  }
  val[i]=0;
  // TODO find conversion which supports long long int
  v = static_cast<key>(strtod(val,&eptr));
  // check for string length == 0
  if (val==eptr)
  {
    fd.v[0]=fd.v[1]=fd.v[2]=0;
    printf("Error in reading vertex index\n");
    return;
  }
  // check for index > max int
  if(v > std::numeric_limits<key>::max())
  {
    cout << "\nError: third vertex index = "
          << v << " is larger than max allowable integer value = "
          << std::numeric_limits<key>::max() << endl;
    exit(0);
  }
  // search for matching vertex index in map
  // map<vertex index from file, relative vertex index>
  p=vp.find(v);
  // if match found
  if(p!=vp.end())
  {
    // set Face vertex key, second
    fd.v[2] = (*p).second;
  }
  else
  {
    // TODO add -p option dependency
    // print details
    cout << "\nError: Face " << fd.index
          << " references missing vertex " 
          << v << endl << line;
    exit(0);
  }
  // check for vertex indices used in Face
  // require indices by nonzero
  if(fd.v[0]>0 && fd.v[1]>0 && fd.v[2]>0)
  {
    if(	fd.v[0]==fd.v[1] ||
        fd.v[1]==fd.v[2] ||
        fd.v[2]==fd.v[0] )
    {
      // print details
      cout << "\nError: Face " << fd.index
            << " does NOT reference three unique vertex indices.\n" 
            << line;
      exit(0);
    }
  }
}

void scanFile (key object_id,std::string filename)
{
  char line[LINE_SIZE];
  FILE *FF = fopen(filename.c_str(),"r");
  if(!FF)
  {
    cout <<"Couldn't open input file " << filename << endl;
    exit(0);
  }
  else
  {
    cerr << "\n\n" << "/* ********************** "
          << "OBJECT ********************** */\n";
    cerr << "name: " << filename << endl;
    cerr.flush();
  }

  // map<vertex index from file, relative vertex index>
  stxxl::map<key,key,kSortComp> vp(NODE_CACHE_SIZE,LEAF_CACHE_SIZE);

  // for every line in file
  for(char *str=fgets(line,LINE_SIZE,FF);str!=NULL;str=fgets(line,LINE_SIZE,FF))
  {
    // skip leading whitespace
    while (strchr(" \t,",*str)!=NULL) { str++;}
    // if first character is V for Vertex
    if (strchr("V",*str)!=NULL)
    {
      Vdata vd;
      ReadVertex(str,vd);
      v_vector.push_back(vd);
      // if index found in map
      // map<vertex index from file, relative vertex index>
      stxxl::map<key,key,kSortComp>::const_iterator p=vp.find(vd.index);
      if(p!=vp.end())
      {
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
        exit(0);
      }
      else
      {
        // ASSUME std::make_pair is compatible with stxxl
        // add vertex to map
        // map<vertex index from file, relative vertex index>
        vp.insert(std::make_pair(vd.index,vertex_index));
      }
      vertex_index++; // relative 
    }
    // if first character is F for Face
    else if (strchr("F",*str)!=NULL)
    {
      Fdata fd;
      ReadFace(str,fd,vp);
      f_vector.push_back(fd);
      face_index++; // relative
    }
  }
  vp.clear();
  fclose(FF);
}

// frozen packet
struct FP
{
  std::string name; // frozen object name
  key index;        // frozen vertex index
};

void parseFrozen(char* str,FP &myfp)
{
  char val[80];
  char *eptr;
  int i;

  char *cp=str;

  // skip leading whitespace
  while (strchr(" \t,",*str)!=NULL) { str++;}
  // grab object name
  i=0;
  while (strchr(" \t,",*str)==NULL)
  {
    val[i++] = *str++;
  }
  val[i]=0;
  if(i==0)
  {
    myfp.name="";
    myfp.index=-1;
    printf("\nError in reading frozen object name: string %s\n",cp);
    exit(0);
  }
  myfp.name=val;

  // skip leading whitespace
  while (strchr(" \t,",*str)!=NULL) { str++;}
  // grab frozen vertex index
  i=0;
  while (strchr("0123456789",*str)!=NULL)
  {
    val[i++] = *str++;
  }
  val[i]=0;
  myfp.index=static_cast<key>(strtod(val,&eptr));
  if (val==eptr)
  {
    myfp.name="";
    myfp.index=-1;
    printf("\nError in reading frozen vertex index: string %s\n",cp);
    exit(0);
  }
}

void processFrozenFile(std::string filename)
{
  key oid = 0;
  bool match=false;
  char line[LINE_SIZE];
  // prep fmap
  // for each vertex in object
  int num=object_vector[oid].v_end
        -object_vector[oid].v_start+1;
  for(int i=0;i<num;i++) // relative
  {
    // access vertex file index
    V* v_ptr = objects[oid].GetV(i);
    frozen[v_ptr->get_index()]=false;
    delete v_ptr;
  }
  // open file
  FILE *FF = fopen(filename.c_str(),"r");
  if(!FF)
  {
    cout <<"Couldn't open input frozen file " << filename << endl;
    exit(0);
  }
  // for every line in file
  for(char *str=fgets(line,LINE_SIZE,FF);str!=NULL;str=fgets(line,LINE_SIZE,FF))
  {
    FP myfp;
    parseFrozen(str,myfp);
    // if object name matches this object's name
    if(myfp.name==objects[oid].GetName())
    {
      // add vertex index to frozen vector
      //frozen.push_back(myfp.index);
      frozen[myfp.index]=true;
      match=true;
    }
    else
    {
      /*
      // if line length is larger than FROZEN_LINE_SIZE
      if(strlen(str)>FROZEN_LINE_SIZE)
      {
        cout << "\nprocessFrozenFile: Error. "
              << "Increase FROZEN_LINE_SIZE=="
              << FROZEN_LINE_SIZE << endl;
        exit(0);
      }
      else
      {
        // add line to storage vector
        storage.push_back(str);
      }
      */
      // add line to storage vector
      storage.push_back(str);
    }
  }
  // if no lines successfully parsed and matches this object
  if(match==false)
  {
    cout << "\n\nprocessFrozenFile: Warning no lines in the frozen"
          << " vertex file\n    " << filename
          << "\nreferenced this object.\n\n";
  }
}

Object* processFile(std::string filename)
{
  key oid = 0;
  // create new object
  Object o(oid,filename);
  objects.push_back(o);
  // initialize vector bounds
  ObjectData od = {0};
  od.v_start = 0; // absolute
  od.f_start = 0; // absolute
  object_vector.push_back(od);
  // vertex and face count in object
  vertex_index=face_index=0; // relative
  // scan file
  scanFile(oid,filename);
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
    object_vector[oid].v_end = od.v_start + vertex_index-1; // absolute
    object_vector[oid].f_end = od.f_start + face_index-1; // absolute
    return &objects.back();
  }
}

void prepFEdata(key oid)
{
  object_vector[oid].fe_start = 0; // absolute
  // for each face in object
  key num=object_vector[oid].f_end-object_vector[oid].f_start+1;
  for(key i=0;i<num;i++) // relative
  {
    FEdata fed;
    fed.e[0]=fed.e[1]=fed.e[2]=-2;
    fe_vector.push_back(fed);
  }
  //setNextFE(object_vector[oid].fe_start + num); // absolute+relative
  object_vector[oid].fe_end = object_vector[oid].fe_start + num-1; // absolute
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

key look4Edge(key va,key vb,map_edge &hm)
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

void halfEdge(key f,key va,key vb,map_edge &hm)
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
  FE* fe_ptr = objects[0].GetFE(f);
  fe_ptr->addEdge(edge_index); // relative
  delete fe_ptr;
  // save new Edata
  e_vector.push_back(ed);
  edge_index++;
}

void buildEdge(key f,key va,key vb,map_edge &hm)
{
  // va and vb passed as v_vector relative indices
  key e = look4Edge(va,vb,hm);
  if(e<0)
  { 
    halfEdge(f,va,vb,hm); 
  }
  else
  {
    E *e_ptr = objects[0].GetE(e);
    e_ptr->updateEdge(f,0);
    delete e_ptr;
  }
}

/*
void cleanEdges()
{
  fe_vector.clear();
  e_vector.clear();
}
*/

void createEdges()
{
  key oid = 0;
  // prep FEdata
  prepFEdata(oid);
  // prepare edge data
  object_vector[oid].e_start = 0;
  edge_index=0;
  // create map for finding edges
  stxxl::map<std::pair<key,key>,key,pSort> hm(NODE_CACHE_SIZE,LEAF_CACHE_SIZE);
  // for each face in object
  key num=object_vector[oid].f_end-object_vector[oid].f_start+1;
  for(key i=0;i<num;i++)
  {
    F* f_ptr = objects[oid].GetF(i);
    if(f_ptr->get_index()>0)
    {
      // TODO pass relative index as key
      buildEdge(f_ptr->get_rel_key(oid),f_ptr->get_v1_rel(),f_ptr->get_v2_rel(),hm);
      buildEdge(f_ptr->get_rel_key(oid),f_ptr->get_v2_rel(),f_ptr->get_v3_rel(),hm);
      buildEdge(f_ptr->get_rel_key(oid),f_ptr->get_v3_rel(),f_ptr->get_v1_rel(),hm);
    }
    delete f_ptr;
  }
  object_vector[oid].e_end = object_vector[oid].e_start + edge_index-1;
}

int main(int argc,char *argv[])
{
  if (argc != 3 && argc != 4)
  {
    printf("\nSyntax: meshrefine input_file threshold [frozen_file]\n\n");
    printf("Detail: Polygons are subdivided so that");
    printf(" no edge is longer than threshold.\n");
    printf("New mesh is written to stdout.\n\n");
    printf("Input mesh must be fully closed, consistently ");
    printf("oriented with vertex and face indices sequentially numbered.\n\n");
    return 1;
  }
  char *eptr;
  double threshold = strtod(argv[2],&eptr)*strtod(argv[2],&eptr);
  // build object
  Object *o = processFile(argv[1]);
  if(o!=NULL)
  {
    if(argc==4)
    {
      processFrozenFile(argv[3]);
    }
    int j=0;
    do {
      fprintf(stderr,"iteration %i: ",j++);
      cerr << "max vertex index " << vertex_index
            << ", max face index " << face_index << endl;
      fprintf(stderr,"Building edges..................");
      fflush(stderr);
      fe_vector.clear();
      e_vector.clear();
      createEdges();
      fprintf(stderr,"complete.\n");
      fflush(stderr);
    } while(refineMesh(threshold)==true);
    // print new mesh to stdout
    printMesh();
    // print updated frozen vertices to new file
    printFrozen(argv[3]);
  }
}
