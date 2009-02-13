
typedef  unsigned long int  u4;   /* unsigned 4-byte type */
typedef  unsigned     char  u1;   /* unsigned 1-byte type */

// ######################################
// ######################################
// declarations - so compiler doesn't complain
class Vertex;
class Face;
class Edge;
class Object;
class Container;
class Box;
class Space;
class Manip;
class void_list;
class Bit;
class Monitor;
class Neighbor;
class Triplet;
class VTrack;
u4 computeHashValue(int);

// ######################################
// ######################################

struct equ4
{
  bool operator()(const u4 s1, const u4 s2) const
  {
    return s1==s2;
  }
};

struct eqf
{
  bool operator()(const Face* s1, const Face* s2) const
  {
    return s1==s2;
  }
};

struct eqv
{
  bool operator()(const Vertex* s1, const Vertex* s2) const
  {
    return s1==s2;
  }
};

struct eqe
{
  bool operator()(const Edge* s1, const Edge* s2) const
  {
    return s1==s2;
  }
};

struct u4_hash
{
  u4 operator()(u4 i) const { return i; }
};

struct f_hash
{
  u4 operator()(Face* i) const { return (u4) i; }
};

struct v_hash
{
  u4 operator()(Vertex* i) const { return (u4) i; }
};

struct e_hash
{
  u4 operator()(Edge* i) const { return (u4) i; }
};

struct ltd
{
  bool operator()(const double s1, const double s2) const
  {
    return s1 < s2;
  }
};

struct gtd
{
  bool operator()(const double s1, const double s2) const
  {
    return s1 > s2;
  }
};

struct lte
{
  bool operator()(const Edge* s1, const Edge* s2) const
  {
    return s1 < s2;
  }
};

struct lto
{
  bool operator()(const Object* s1, const Object* s2) const
  {
    return s1 < s2;
  }
};

struct ltv
{
  bool operator()(const Vertex* s1, const Vertex* s2) const
  {
    return s1 < s2;
  }
};

struct ltf
{
  bool operator()(const Face* s1, const Face* s2) const
  {
    return s1 < s2;
  }
};

struct lts
{
  bool operator()(const std::string s1, const std::string s2) const
  {
    return s1 < s2;
  }
};

typedef std::map<double,int,ltd,std::allocator<int> > map_di;
typedef std::map<double,int,ltd,std::allocator<int> >::iterator di_iterator;
typedef std::map<std::string,Edge*,lts,std::allocator<Edge*> > hashtable_t;
typedef std::map<std::string,Edge*,lts,std::allocator<Edge*> >::iterator ht_iterator;
typedef std::map<std::string,double,lts,std::allocator<double> > map_sd;
typedef std::map<std::string,double,lts,std::allocator<double> >::iterator sd_iterator;
typedef std::multimap<Face*,Box*,ltf,std::allocator<Box*> > hashtable_f;
typedef std::multimap<Face*,Box*,ltf,std::allocator<Box*> >::iterator tf_iterator;
typedef __gnu_cxx::hash_map<Vertex*,int,v_hash,eqv,std::allocator<int> > hashtable_v;
typedef __gnu_cxx::hash_map<Vertex*,int,v_hash,eqv,std::allocator<int> >::iterator vhm_iterator;
typedef __gnu_cxx::hash_map<Face*,double*,f_hash,eqf,std::allocator<double*> > hashtable_f_double;
typedef __gnu_cxx::hash_map<Vertex*,double,v_hash,eqv,std::allocator<double> > hashtable_v_double;
typedef __gnu_cxx::hash_map<Vertex*,double,v_hash,eqv,std::allocator<double> >::iterator vdhm_iterator;
typedef __gnu_cxx::hash_map<Edge*,double,e_hash,eqe,std::allocator<double> > hashtable_e_double;
typedef __gnu_cxx::hash_map<Edge*,double,e_hash,eqe,std::allocator<double> >::iterator edhm_iterator;
typedef __gnu_cxx::hash_map<Face*,std::vector<Face*>*,f_hash,eqf,std::allocator<std::vector<Face*>* > > hashtable_f_face;
typedef __gnu_cxx::hash_map<Face*,std::vector<Face*>*,f_hash,eqf,std::allocator<std::vector<Face*>* > >::iterator htff_iterator;
typedef __gnu_cxx::hash_map<Edge*,int,e_hash,eqe,std::allocator<int> > hashtable_e;
typedef std::map<Vertex*,double,ltv,std::allocator<double> > table_d;
typedef std::map<Vertex*,double,ltv,std::allocator<double> >::iterator td_iterator;
typedef std::map<Face*,double*,ltf,std::allocator<double*> > table_fd;
typedef std::map<Face*,double*,ltf,std::allocator<double*> >::iterator fd_iterator;
typedef std::multimap<double,Vertex*,gtd,std::allocator<Vertex*> > table_v;
typedef std::multimap<double,Vertex*,gtd,std::allocator<Vertex*> >::iterator tv_iterator;
typedef std::multimap<double,Face*,ltd,std::allocator<Face*> > mmap_d_f;
typedef std::multimap<double,Face*,ltd,std::allocator<Face*> >::iterator df_iterator;
typedef std::pair<double,Vertex*> vd_pair;
typedef std::pair<double,double> dd_pair;
typedef std::set<Edge*,lte> e_set;
typedef std::set<Edge*,lte>::iterator es_iterator;
typedef std::set<Vertex*,ltv> v_set;
typedef std::set<Vertex*,ltv>::iterator vs_iterator;
typedef std::set<Face*,ltf> f_set;
typedef std::set<Face*,ltf>::iterator fs_iterator;
typedef std::list<Vertex*> v_list;
typedef std::list<Vertex*>::iterator vl_iterator;
typedef __gnu_cxx::hash_set<Face*,f_hash,eqf> hashset_f;
typedef __gnu_cxx::hash_set<Face*,f_hash,eqf>::iterator hf_iterator;
typedef __gnu_cxx::hash_set<Vertex*,v_hash,eqv> hashset_v;
typedef __gnu_cxx::hash_set<Vertex*,v_hash,eqv>::iterator hv_iterator;
typedef __gnu_cxx::hash_set<Edge*,e_hash,eqe> hashset_e;
typedef __gnu_cxx::hash_set<Edge*,e_hash,eqe>::iterator he_iterator;
typedef std::map<double,Face*,ltd,std::allocator<Face*> > map_df;
typedef std::map<double,Face*,ltd,std::allocator<Face*> >::iterator mdf_iterator;
typedef std::map<double,Triplet*,ltd,std::allocator<Triplet*> > map_dt;
typedef std::map<double,Triplet*,ltd,std::allocator<Triplet*> >::iterator mdt_iterator;
typedef std::multimap<Object*,int,lto,std::allocator<int> > mmap_oi;
typedef std::multimap<Object*,int,lto,std::allocator<int> >::iterator oi_iterator;

// ######################################
// ######################################

int distinguishable(double a,double b,double epsilon) {
  double c;
  c=a-b;
  if (c<0) c=-c;
  if (a<0) a=-a;
  if (a<1) a=1;
  if (b<0) b=-b;
  if (b<a) return (c>a*epsilon);
  else return (c>b*epsilon);
}

int distinguishable(double a,double b) {
  double c;
  c=a-b;
  if (c<0) c=-c;
  if (a<0) a=-a;
  if (a<1) a=1;
  if (b<0) b=-b;
  if (b<a) return (c>a*DOUBLE_EPSILON);
  else return (c>b*DOUBLE_EPSILON);
}

// #####################################################
// #####################################################

class void_list
{
public:
  void_list *previous;
  void_list *next;
  void *data;
};

// ######################################
// ######################################

class Triplet{
public:
  double x,y,z;
  Triplet(double,double,double);
};

Triplet::Triplet(double j,double k,double l){
  x=j;
  y=k;
  z=l;
}

class Point{
public:
  double x,y,z; // current vertex
  double a,b,c,L; // closest point and squared distance
  void add(double,double,double);
  Point(double,double,double);
  void clear(void);
};

void Point::add(double j,double k,double l){
  double t=(x-j)*(x-j)+(y-k)*(y-k)+(z-l)*(z-l);
  if (t<L){a=j;b=k;c=l;L=t;}
}

Point::Point(double j,double k,double l){
  x=j;
  y=k;
  z=l;
  a=b=c=0.0;
  L=1e300;
}

void Point::clear(void){
  a=b=c=0.0;
  L=1e300;
}

// ######################################
// ######################################

class Vertex {
public:
  int index;
  double pN[3];		// current position coordinates (x,y,z)
  Face *cl;			// pointer to face on which closest mesh position lies
  Object *o;			// pointer to parent object
  vector<Face*> f;	// pointers to adjacent faces
  vector<Face*> nf;	// pointers to neighborhood faces
  // NOTE FULL HOOD = neighborhood+adjacent
  Vertex(char* triplet,Object*);
  double energy;
  void getAdjacentVertices(vector<Vertex*>&);
  void getAdjacentFaces(hashset_f&);
  void getAdjacentEdges(std::vector<Edge*>&);
  double getSqSepDist(Container*); 
  void getNormal(double*);
  void getForceEnergy(double[3],Container*);
  double getSeparationForceEnergy(double[3],bool,Container*);
  double getEdgeStretchForceEnergy(double[3],bool,Container*);
  double getEdgeAngleFaceIntersectionForceEnergy(double[3],bool);
  void computeNewCoords(Container*,double[3],double);
  void assignHolding(double[3]);
  void printNeighborhood(void);
  void printVertex(std::string);
  void printVertexCP(void);
  void getBoundingBox(double[6]);
  //////////////
  bool vertexIsNice(void);
  int getVertexNiceness(void);
  void setVertexNiceness(int,Container*);
  void computeEdgeFlip(void);
  bool match(int,std::string);
};

void Vertex::printVertex(std::string s){
  cout.precision(12);
  cout << "Vertex <obj>" << s << "<ind>" << index << " "
        << "["
        << pN[0] << " "
        << pN[1] << " "
        << pN[2] << "] ";
}

void Vertex::printVertexCP(void){
  cout.precision(12);
  cout
        << pN[0] << " "
        << pN[1] << " "
        << pN[2] << " 1 0 0 1";
}

Vertex::Vertex(char* triplet,Object *q) {

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
  index = (int) strtod(val,&eptr);
  if (val==eptr)
  {
    index=0;
    printf("Error in reading vertex index\n");
    return;
  }

  // grab x coord
  while (strchr(" \t,",*triplet)!=NULL) { triplet++; }
  i=0;
  while (strchr("0123456789+-eE.",*triplet)!=NULL)
  {
    val[i++] = *triplet++;
  }
  val[i]=0;
  pN[0]=strtod(val,&eptr);
  if (val==eptr)
  {
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
  pN[1]=strtod(val,&eptr);
  if (val==eptr)
  {
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
  pN[2]=strtod(val,&eptr);
  if (val==eptr)
  {
    printf("Error in reading vertex\n");
    printf("Error in reading vertex: string %s\n",cp);
    return;
  }
  o=q;
  cl=NULL;
  energy=0;
}

// ######################################
// ######################################

class Face {
public:
  int index;		// face index
  Vertex *v[3];	// pointers to vertices
  Edge   *e[3];	// pointers to edges
  double bb[6];	// bounding box [xmin xmax ymin ymax zmin zmax]
  Face(char*,std::vector<Vertex*>&); 
  void addEdge(Edge*);
  void getNormal(double[3]);	// compute outward normal vector
  void getVertexCoordinates(double *[3]);
  bool computeIntersectionForce(Container*,Space&);
  void calculateIntersectionForce(Container*);
  //	bool getFaceIntersection(Container*,std::vector<Face*>&,Space&);
  bool getFaceIntersection(Container*,Space&);
  bool getFaceIntersectionCheck(Container*,Space&,hashtable_f&);
  double getAngle(Vertex *v);
  void printFace(std::string);
  void printFaceCP(void);
  void updateBoxes(hashtable_f&,hashtable_f&);
  void getBoundingBox(void);
  /////////
  std::vector<Face*> * getIntersectingFaces(void);
  void addFaceToTable_iv(void);
  void addFaceToTable_intf(void);
  bool faceInTable_iv(void);
  bool faceInTable_intf(void);
  void addForceToFace(double[3]);
  void addFaceToVector(Face*,Container*);
  void clearFaceFromTable_iv(void);
  //	void clearFaceFromTable_intf(void);
  void clearFaceFromTable_intf(Container*);
  void getForceFromTable(double[3]);
  bool noMoreIntersectingFaces(void);
  //	void removeFaceFromVector(Face*);
  void removeFaceFromVector(Face*,Container*);
  bool match(int,std::string);
  /////////
  bool faceInVector(Face*);
  double getAspectRatio(void);
};

// ######################################
// ######################################

class Edge {
public:
  Face *f1,*f2;	// pointers to adjacent faces (i.e. faces that contain edge)
  double l;		// original edge length
  void update(Face*);
  Edge(Face*,Vertex*,Vertex*);
  double getSqLength(void);
  int computeFlip(void);
  double getAngle(void);
  //double getForceEnergy(int,double[3],bool);  
  double getForceEnergy(int,double[3],bool,bool);  
  double getReactionForceEnergy(double[3],bool,bool);  
  double getStretchForceEnergy(Vertex*,double[3],bool,Container*);
  void printEdge(std::string);
  void getVertices(Vertex*&,Vertex*&,Vertex*&,Vertex*&);
  double energy;
};

void Edge::printEdge(std::string s){
  Vertex *v1=NULL,*v2=NULL,*o1=NULL,*o2=NULL;
  getVertices(v1,v2,o1,o2);
  cout.precision(12);
  if (v1==NULL || v2==NULL ||o1==NULL ||o2==NULL ||f1==NULL ||f2==NULL){
    cout << "\nNULL POINTERS IN EDGE!\n";
    if(v1==NULL){cout << "v1 is NULL\n";}
    if(v2==NULL){cout << "v2 is NULL\n";}
    if(o1==NULL){cout << "o1 is NULL\n";}
    if(o2==NULL){cout << "o2 is NULL\n";}
    if(f1==NULL){cout << "f1 is NULL\n";}
    if(f2==NULL){cout << "f2 is NULL\n";}
    exit(0);
  }
  cout << "\nprintEdge: <obj>" << s
        << " v1 "<< v1->index << " ["
        << v1->pN[0] << " "
        << v1->pN[1] << " "
        << v1->pN[2] << "],"
        << " v2 "<< v2->index << " ["
        << v2->pN[0] << " "
        << v2->pN[1] << " "
        << v2->pN[2] << "]\n"
        << "printEdge: <obj>" << s
        << " o1 "<< o1->index << " ["
        << o1->pN[0] << " "
        << o1->pN[1] << " "
        << o1->pN[2] << "],"
        << " o2 "<< o2->index << " ["
        << o2->pN[0] << " "
        << o2->pN[1] << " "
        << o2->pN[2] << "]\n"
        << "printEdge: <obj>" << s
        << " f1 "<< f1->index << " ["
        << f1->v[0]->index << " "
        << f1->v[1]->index << " "
        << f1->v[2]->index << "],"
        << " f2 "<< f2->index << " ["
        << f2->v[0]->index << " "
        << f2->v[1]->index << " "
        << f2->v[2]->index << "]\n"
        << "printEdge: <obj>" << s
        << " orignal length "<< l
        << endl;
}

Edge::Edge(Face *f,Vertex *va,Vertex *vb){
  f1=f;
  f2=NULL;
  // compute original edge length
  l=sqrt((va->pN[0]-vb->pN[0])*(va->pN[0]-vb->pN[0])+
         (va->pN[1]-vb->pN[1])*(va->pN[1]-vb->pN[1])+
         (va->pN[2]-vb->pN[2])*(va->pN[2]-vb->pN[2]));
  energy=0;
}

double Edge::getSqLength(void){
  Vertex *v1=NULL,*v2=NULL,*o1=NULL,*o2=NULL;
  getVertices(v1,v2,o1,o2);
  return (v1->pN[0]-v2->pN[0])*(v1->pN[0]-v2->pN[0])
        +(v1->pN[1]-v2->pN[1])*(v1->pN[1]-v2->pN[1])
        +(v1->pN[2]-v2->pN[2])*(v1->pN[2]-v2->pN[2]);
}

// ######################################
// ######################################

class Object
{
public:
  std::string name;	// object name
  vector<Vertex*> v;		// container of pointers to all vertices in object
  vector<Face*> f;		// container of pointers to all faces in object
  vector<Edge*> e;		// container of pointers to all edges in object
  Object(std::string);
  double mean_edge_length;
  ~Object(void);
  void addOriginal(int,double*);
  void createEdges(void);
  int setNumDigits(void);
  void checkEdge(Face*,Vertex*,Vertex*,hashtable_t&,int);
  Edge* findEdge(Vertex*,Vertex*,hashtable_t&,int);
  void createEdge(Face*,Vertex*,Vertex*,hashtable_t&,int);
  void addVertexPointers(void);
  int getMaxVertex(void);
  void fixFaces(Vertex**);
  void fixEdges(Vertex**);
  void findVertexAdjacencies(void);
  void findNeighborhoods(void);
  void boundObject(double*);
  double getMeanEdgeLength(void);
  //////////////
  hashtable_f_double iv; // store intersection force
  hashtable_f_face intf; // store intersecting faces
  bool intersectingFacesExist(Face*);
  bool faceInTable_intf(Face*f);
  bool faceInTable_iv(Face*f);	
  //////////////
  //	assume all vertices are initially nice, so nice is empty
  //	nice=0;// since initial # nonnice==0 (e.g. s_nonnice=nonnice=0;)
  hashtable_v nice; // vertex* is key to int nice value
  bool vertexIsNice(Vertex*);
  int getVertexNiceness(Vertex*);
  void setVertexNiceness(Vertex*,int);
  //////////////
  hashtable_e flip;	// collection of flip values for edges
  void addEdgeFlip(Edge*,int);
  void clearFlipTable(void);
  int getEdgeFlip(Edge*);
  //////////////
  bool ifFrozen(hashtable_v_double&,Vertex*);
  void newFindNeighborhoods(void);
  bool thawedAndAble(hashtable_v_double&,v_set&);
  void collectFaces(hashtable_v_double&,v_set&,std::vector<Face*>&);
  bool processEdge(Edge*,hashtable_v_double&,std::vector<Edge*>&,Vertex*);
  void buildNeighborhood(Vertex*);
  //////////////
};

bool Vertex::match(int i, std::string str){
  return i==index && str==o->name;
}

bool Face::match(int i, std::string str){
  return i==index && str==v[0]->o->name;
}

void Object::clearFlipTable(void){
  flip.clear();
}

void Object::addEdgeFlip(Edge* ee,int i){
  // if no element in hashtable has this key
  if(flip.find(ee)==flip.end()){flip[ee]=i;}
  else {cout << "Edge element already in hashtable.\n";exit(0);}
}

int Object::getEdgeFlip(Edge* ee){
  // if element with this edge key exists in hashtable
  if(flip.find(ee)!=flip.end()){return flip[ee];}
  else {cout << "Edge element does not exist in hashtable.\n";exit(0);}
}

bool Vertex::vertexIsNice(void){
  return o->vertexIsNice(this);	
}

bool Object::vertexIsNice(Vertex *vv){
  return nice.find(vv)==nice.end();	
}

int Vertex::getVertexNiceness(void){
  return o->getVertexNiceness(this);
}

int Object::getVertexNiceness(Vertex *vv){
  if (vertexIsNice(vv)){
    return 0;
  } else {
    return nice[vv];
  }
}

void Object::setVertexNiceness(Vertex *vv,int val){
  if (val==0){
    nice.erase(vv);
  } else {
    nice[vv]=val;
  }
}

Object::~Object(void){
  std::vector<Vertex*>::iterator i;
  std::vector<Face*>::iterator j;
  std::vector<Edge*>::iterator k;
  for(i=v.begin();i!=v.end();i++){ delete *i; }
  for(j=f.begin();j!=f.end();j++){ delete *j; }
  for(k=e.begin();k!=e.end();k++){ delete *k; }
}


Object::Object(std::string s) {
  name=s;
}

bool Object::faceInTable_intf(Face *ff){
  return intf.find(ff)!=intf.end();
}

bool Object::faceInTable_iv(Face *ff){
  return iv.find(ff)!=iv.end();
}

void Face::addFaceToTable_iv(void){
  // create new force*
  double* fp = new double[3];
  fp[0]=0;
  fp[1]=0;
  fp[2]=0;
  // create element in table
  v[0]->o->iv[this]=fp;
}

void Face::addFaceToTable_intf(void){
  // create new face vector in table
  std::vector<Face*> *nfv = new std::vector<Face*>();
  v[0]->o->intf[this]=nfv;
}

bool Face::faceInTable_iv(void){
  return v[0]->o->faceInTable_iv(this);
}

bool Face::faceInTable_intf(void){
  return v[0]->o->faceInTable_intf(this);
}

void Face::addForceToFace(double fvec[3]){
  v[0]->o->iv[this][0]+=fvec[0];
  v[0]->o->iv[this][1]+=fvec[1];
  v[0]->o->iv[this][2]+=fvec[2];	
}
/*
    void Face::addFaceToVector(Face* f){
    (*v[0]->o->intf[this]).push_back(f);
    }*/

void Face::getForceFromTable(double fvec[3]){
  fvec[0]+=v[0]->o->iv[this][0];
  fvec[1]+=v[0]->o->iv[this][1];
  fvec[2]+=v[0]->o->iv[this][2];	
}

void Face::clearFaceFromTable_iv(void){
  // if this face is in iv table 
  if(faceInTable_iv()){
    // delete force*
    delete[] v[0]->o->iv[this];
    // remove element from table
    v[0]->o->iv.erase(this);
  }
}

/*
    void Face::removeFaceFromVector(Face *f){
// if intersecting face has intersecting faces
if(faceInTable_intf()){
// remove face from intersecting faces vector
std::vector<Face*> *ifv=this->getIntersectingFaces();
(*ifv).erase(remove((*ifv).begin(),(*ifv).end(),f),(*ifv).end());
// if intersecting face vector is now empty
if (noMoreIntersectingFaces()){
// then remove intersecting face vector from hashtable
clearFaceFromTable_intf();
// then update face intersect force, i.e. set force to zero
clearFaceFromTable_iv();
}
}
}*/

bool Face::noMoreIntersectingFaces(void){
  return (*v[0]->o->intf[this]).empty();
}

std::vector< Face* > * Face::getIntersectingFaces(void){
  return v[0]->o->intf[this];
}

// THE FOLLOWING CODE IS USEFUL. IT PRINTS IN COMPREHENSIVE FORMAT.
void Face::printFace(std::string s){
  cout.precision(12);
  cout << "Face <obj>" << s << "<ind>" << index << endl
        << "[v0 "
        << v[0]->index << " "
        << v[0]->pN[0] << " "
        << v[0]->pN[1] << " "
        << v[0]->pN[2] << "]\n"
        << "[v1 "
        << v[1]->index << " "
        << v[1]->pN[0] << " "
        << v[1]->pN[1] << " "
        << v[1]->pN[2] << "]\n"
        << "[v2 "
        << v[2]->index << " "
        << v[2]->pN[0] << " "
        << v[2]->pN[1] << " "
        << v[2]->pN[2] << "]";
}

// THE FOLLOWING CODE IS USEFUL. IT PRINTS IN DReAMM CUSTOM POINTS FORMAT.
void Face::printFaceCP(void){
  cout << v[0]->pN[0] << " "
        << v[0]->pN[1] << " "
        << v[0]->pN[2] << " 1 0 0 1\n"
        << v[1]->pN[0] << " "
        << v[1]->pN[1] << " "
        << v[1]->pN[2] << " 1 0 0 1\n"
        << v[2]->pN[0] << " "
        << v[2]->pN[1] << " "
        << v[2]->pN[2] << " 1 0 0 1\n";
}

void Face::getVertexCoordinates(double *cpvc[3]){
  cpvc[0]=v[0]->pN;
  cpvc[1]=v[1]->pN;
  cpvc[2]=v[2]->pN;
}

Face::Face(char *triplet,std::vector<Vertex*> &vp){

  e[0]=e[1]=e[2]=NULL;

  char val[80];
  char *eptr;
  int i;

  // get past 'Face'
  while (strchr("Face",*triplet)!=NULL) {triplet++;}

  // grab Face index
  while (strchr(" \t,",*triplet)!=NULL) { triplet++; }
  i=0;
  while (strchr("0123456789+-eE.",*triplet)!=NULL)
  {
    val[i++] = *triplet++;
  }
  val[i]=0;
  index = (int) strtod(val,&eptr);
  if (val==eptr)
  {
    index=0;
    v[0]=v[1]=v[2]=0;
    printf("Error in reading face index\n");
    return;
  }

  // grab first vertex index
  while (strchr(" \t,",*triplet)!=NULL) { triplet++; }
  i=0;
  while (strchr("0123456789+-eE.",*triplet)!=NULL)
  {
    val[i++] = *triplet++;
  }
  val[i]=0;
  v[0] = vp[(int)strtod(val,&eptr)-1];
  if (val==eptr)
  {
    v[0]=v[1]=v[2]=0;
    printf("Error in reading vertex index\n");
    return;
  }

  // grab second vertex index
  while (strchr(" \t,",*triplet)!=NULL) triplet++;
  i=0;
  while (strchr("0123456789+-eE.",*triplet))
  {
    val[i++] = *triplet++;
  }
  val[i]=0;
  v[1] = vp[(int)strtod(val,&eptr)-1];
  if (val==eptr)
  {
    v[0]=v[1]=v[2]=0;
    printf("Error in reading vertex index\n");
    return;
  }

  // grab third vertex index
  while (strchr(" \t,",*triplet)!=NULL) triplet++;
  i=0;
  while (strchr("0123456789+-eE.",*triplet))
  {
    val[i++] = *triplet++;
  }
  val[i]=0;
  v[2] = vp[(int)strtod(val,&eptr)-1];
  if (val==eptr)
  {
    v[0]=v[1]=v[2]=0;
    printf("Error in reading vertex index\n");
    return;
  }

  getBoundingBox();

}


// ######################################
// ######################################

class Container
{
public:
  std::ofstream Sfile;
  std::ofstream Cfile;
  std::ofstream Olist;
  vector<Object*> o;	//vector of object pointers
  int num_files;		// number of input files
  std::vector<std::string> files;	// array of input file names
  double energy;		// total energy in all objects (vertices, faces, edges)
  double force;		// total force in all objects
  int nonnice;		// total # nonnice vertices in all objects
  int s_nonnice;		// total # nonnice vertices that are inside self in all objects
  /////
  int si;				// twice total # self-intersections in all objects
  int ti;				// twice total # intersections in all objects
  // because each intersecting pair
  // of faces will be discovered twice
  double md[2];		///// mean displacement of vertices in object during each group
  // md[0] = mean displacement of N-1 vertices
  // md[1] = mean displacement of N vertices
  double d_min,d_max;	// min and max actual displacement
  std::vector<Vertex*> frozen;	// sorted vector of Vertex* to frozen vertices
  double min_edge_angle;
  int N;
  void clear(void);
  void scanDir(void);
  void scanFiles(void);
  void scanFile(Object*,char*);
  void createEdges(void);
  void testData(void);
  void addVertexPointers(void);
  void findVertexAdjacencies(void);
  void writeObjectList(void);
  void statusFileInit(void);
  void sepFileInit(int);
  void sepOrigLog(Monitor&);
  void sepSnapshotLog(Monitor&,int);
  void buildMeshAfter(int const);
  void writeDistances(int);
  void writeDistances50(int);
  void writeDistancesNOCP(int);
  void writeDistancesOneSidedHausdorff_noself(void);
  void writeDistancesTwoSidedHausdorff_noself(void);
  void writeDistancesOneSidedHausdorff(void);
  void writeDistancesTwoSidedHausdorff(void);
  void computePC(Face*,Vertex*,double[3]);
  void fileOutit(void);
  void sepFileOutit(void);
  //	void updateFile(int,bool,double,double);
  void updateFile(int,bool,double);
  void updateSepFile(double,std::string,int,double);
  void checkAdjacentFaces(void);
  void checkAngle(double);
  void cleanup(void);
  ~Container(void);
  Container(void);
  void updateStats(double);
  void removeOldIntersections(Vertex*,hashset_v&);
  void updateNewIntersections(Vertex*,hashset_v&,Space&);
  //// from Manip
  std::ofstream Mfile;
  int object_count;				// total number of objects in model
  int vertex_count;				// total number of vertices in model
  int face_count;					// total number of faces in model
  int edge_count;					// total number of edges in model
  bool self_intersection;			// flag signifies if self-intersecting object was found
  int pairs[3][2];
  ///////////
  ///////////
  ///////////
  void getNonnice(hashset_v&);
  void getIntersectedVertices(hashset_v&);
  void checkEdgeAngles(void);
  double computeSeparationError(Vertex*,double);
  ///////////
  void getNiceSet(Space&,Monitor&);
  void findNice(Space&);
  bool checkNiceness(Space&,Vertex*);
  void collectCrossed(Space&,Vertex*,std::vector<Object*>&);
  bool updateNiceness(Vertex*,std::vector<Object*>&);	
  void findClosestAxis(Space&,Vertex*,double[2][3],double [6][3]);
  int findExtraPoint(Space&,Vertex*,double[3],int);
  //	void findCrossed1(Space&,Vertex*,double[2][3],std::vector<Object*>&);
  bool findCrossed1(Space&,Vertex*,double[2][3],std::vector<Object*>&);
  //	void findCrossed2(Space&,double[2][3],std::vector<Object*>&);
  bool findCrossed2(Space&,double[2][3],std::vector<Object*>&);
  void getExtraRay(Vertex*,double[2][3],int);
  void collectNiceFaces(Space&,double[2][3],vector<Face*>&);
  void getBoxIndexRange(Space&,double[2][3],int[6]);
  void findIntersectedFaces(double[2][3],vector<Face*>&,
                            std::vector<Face*>&,std::vector<Face*>&);
  void findOddMeshes(std::vector<Face*>&,std::vector<Face*>&,std::vector<Object*>&);
  void findOddObjects(std::vector<Object*>&,std::vector<Object*>&);
  void findEvenObjects(std::vector<Object*>&,std::vector<Object*>&);
  ///////////
  void getSeparationDistancesAndECW(Space&,Monitor&);
  void getSeparationDistances(Space&,Monitor&);
  bool findClosest(Space&,Vertex*,Monitor&,bool);
  bool computeClosest(Face*,Vertex*,double&,double[3]);
  bool computeClosestColl(std::vector<Vertex*>&,std::vector<Edge*>&,std::vector<Face*>&,
                          Vertex*,double&,double[3]);
  bool getPlaneIntersection(Face*,Vertex*,double*,double,double,Point&);
  void getEdgeIntersection(Vertex*,double*[3],Point&);
  void getBoxes(std::vector<Box*>&,Vertex*,Space&);
  void getCandidateFaces(std::vector<Box*>&,Vertex*,std::vector<Face*>&);
  bool faceInNeighborhood(Face*,Vertex*);
  bool faceIsAdjacent(Face*,Vertex*);
  ///////////
  //	bool assignNewVertexCoords(Space&,Vertex*,double[3],Monitor&); 
  bool assignNewVertexCoords(Space&,Vertex*,double[3],Monitor&,bool&,bool&); 
  void getAffectedVerticesAndEdgesBefore(Space&,Vertex*,double[3],Monitor&);
  void getAffectedVerticesAndEdgesAfter(Space&,Vertex*,double[3],Monitor&);
  bool boxesOverlap(double[6],double[6]);
  void computeNewVertexCoords(double); 
  void computeFaceIntersectionForce(Space&); 
  void setFaceIntersection(Face*,Space&); 
  //	bool checkForIntersections(Vertex*,Space&,bool,hashtable_f&);
  bool checkForIntersections(Vertex*,Space&,hashtable_f&);
  bool checkForSmallAngles(Monitor&);
  bool angleChangeIsWrong(double,double);
  void collectEdgeAngles(Vertex*,Monitor&);
  void updateClosest(Space&,Vertex*,Monitor&);
  ///////////
  bool checkFaceFaceIntersections(Face*,Face*);
  int checkEdgeEdgeIntersection(Face*,Face*,bool);
  int checkFaceEdgeIntersection(Face*,Face*);
  bool facesParallel(Face*,Face*);
  bool facesColinear(Face*,Face*);
  int numUniqueVertices(Face*,Face*,int[2]);
  ///////////
  void assignFacesToBoxes(Space&);
  void identifyBoxes(Space&,Face*,vector<Box*>&);
  /////////
  void updateAdjacentFaceBoxes(Vertex*,Monitor&);
  void updateAdjacentFaceBoundingBoxes(Vertex*);
  void getNiceCheckSet(Vertex*,Monitor&,Space&);
  double getVertexSqD(Vertex*,double);
  void collectAdjacentFaceNormals(table_fd&,Vertex*);
  void freeAdjacentFaceNormals(table_fd&);
  void getVertexAndEdgeEnergy(Monitor&);
  void computeGlobalEnergy(void);
  void updateMovedVertexEnergy(Vertex*,Monitor&);
  void updateEdgeAngles(Vertex*);
  void updateVertexVD(Vertex*,Monitor&);
  bool searchCP(Vertex*,std::vector<Face*>&,double&);
  void getPD2(double[6],double[6],Vertex*);
  //
  Object* getObjectPointer(char[FILENAME_SIZE]);
  void loadFrozenMap(mmap_oi&,const char*);
  void readFrozen (const char*);
  void printFrozenCP(void);
  void checkFrozenAndNoCP(void);
  void writeObjectData(void);
  void writeSeparationDistances(void);
  void computeMeanEdgeLengths(void);
  void sortAdjacentFaces(void);
  Face* getNonAdjacentFaceOfVertex(Vertex*,Vertex*);
  Face* getNonAdjacentFaceOfEdge(Vertex*,Edge*);
  void summarize1(void);
  void summarize2(void);
  void countNonnice(int&,int&);
  void countIntersections(int&,int&);
  void reportVertexNiceness(int,std::string,Space&);
  void gatherDiagnostics(void);
  void printVD(Monitor&);

};

void Container::summarize1(void){
  // for each object, accumulate number of vertices and faces
  for (std::vector<Object*>::iterator i=o.begin();i!=o.end();i++) {
    object_count++;
    vertex_count += (*i)->v.size();
    face_count += (*i)->f.size();
  }
}

void Container::summarize2(void){
  // for each object, accumulate number of vertices and faces
  for (std::vector<Object*>::iterator i=o.begin();i!=o.end();i++) {
    edge_count += (*i)->e.size();
  }
}

Container::Container(void){
  si=ti=0;				// twice number of intersecting faces
  // in model, self and total, respectively
  min_edge_angle=500;		// minimum edge angle in model
  s_nonnice=nonnice=0;	// total number of nonnice vertices
  // in model, self and total, respectively
  energy=0; // NOT USED
  force=0; // NOT USED
  pairs[0][0] = 0;
  pairs[0][1] = 1;
  pairs[1][0] = 1;
  pairs[1][1] = 2;
  pairs[2][0] = 2;
  pairs[2][1] = 0;
  object_count=vertex_count=face_count=edge_count=0;
  scanDir();
  scanFiles();
  summarize1();
  createEdges();
  summarize2();
  findVertexAdjacencies();
  sortAdjacentFaces();
  checkEdgeAngles();
  computeMeanEdgeLengths();
  N=0;			// count of moved vertices
  md[0]=md[1]=0;	// mean displacement of N moved vertices
  d_min=1E30;		// min displacement of N moved vertices
  d_max=-1E30;	// max displacement of N moved vertices
}

Container::~Container(void){
  std::vector<Object*>::iterator i;
  for(i=o.begin();i!=o.end();i++){ delete *i; }
}

void Container::checkAngle(double a){
  // a is edge angle
  if (a<min_edge_angle){min_edge_angle=a;}
}

void Container::clear(void){
  N=0;			// count of moved vertices
  md[0]=md[1]=0;	// mean displacement of N moved vertices
  d_min=1E30;		// min displacement of N moved vertices
  d_max=-1E30;	// max displacement of N moved vertices
}

// ######################################
// ######################################

class Box {
public:
  vector<Face*> f;			// vector of face pointers
  int x,y,z;					// indices of box in space
  Box(int,int,int);
  double xmin(double);
  double xmax(double);
  double ymin(double);
  double ymax(double);
  double zmin(double);
  double zmax(double);
  void printBox(double,double,double);
};

void Box::printBox(double a,double b,double c){
  cout << "Box::printBox: world 0,2,4 ["
        << a << " "
        << b << " "
        << c << "]";
  cout.flush();
  cout << ", Box [";
  cout.flush();
  cout << x << " ";
  cout.flush();
  cout << y << " ";
  cout.flush();
  cout << z << "]\n";
  cout.flush();
  cout << "Box::printBox: "
        << "limits ["
        << xmin(a) << " "
        << xmax(a) << " "
        << ymin(b) << " "
        << ymax(b) << " "
        << zmin(c) << " "
        << zmax(c) << "]\n";
  cout.flush();
}

double Box::xmin(double i){return x*SPACE_LENGTH*SCALE+i;}
double Box::xmax(double i){return (x+1)*SPACE_LENGTH*SCALE+i;}
double Box::ymin(double i){return y*SPACE_LENGTH*SCALE+i;}
double Box::ymax(double i){return (y+1)*SPACE_LENGTH*SCALE+i;}
double Box::zmin(double i){return z*SPACE_LENGTH*SCALE+i;}
double Box::zmax(double i){return (z+1)*SPACE_LENGTH*SCALE+i;}

Box::Box(int a,int b,int c){
  x=a;
  y=b;
  z=c;
}

// ######################################
// ######################################

class Space {
public:
  std::vector<double> xv,yv,zv;
  std::ofstream Sfile;
  int num_space[3];
  int num_boxes;
  vector<Box*> b;	// vector of boxes
  double world[6];	// minX, maxX, minY, maxY, minZ, maxZ of world
  void boundWorld(Container&);
  void initBoxes(void);
  void clearBoxes(void);
  void recordFace(std::vector<Box*>&,Face*);
  void computeBoxesToCheck(Face*,std::vector<Box*>&);
  void fileInit(void);
  void index2Range(int,double[2],char*);
  int location2Index(double,char*);
  int indices2Index(int,int,int);
  void getBoxesFor3DLocations(double[6],std::vector<Box*>&);
  void getBoxesFor3DIndices(int[6],std::vector<Box*>&);
  int screenIndex(int,char*);
  Space(Container&);
  ~Space(void);
  void getSR(double[6],int[6]);
  bool expandSR(bool,double,double[6],int[6],int[6]);
};

Space::Space(Container& c){
  cout << "Initializing space data structure..............";
  cout.flush();	
  xv.reserve(3);
  yv.reserve(3);
  zv.reserve(3);
  boundWorld(c);
  initBoxes();
  fileInit();
  cout << "complete.\n";
  cout.flush();	
}

Space::~Space(void){
  std::vector<Box*>::iterator i;
  for(i=b.begin();i!=b.end();i++){
    delete *i;
  }
}

int Space::screenIndex(int i,char *c){
  int a=0;
  if		(!strcmp(c,"x")){a=num_space[0]-1;}
  else if (!strcmp(c,"y")){a=num_space[1]-1;}
  else if (!strcmp(c,"z")){a=num_space[2]-1;}
  else {cout << "Received unexpected string.\n";exit(0);}
  //
  if (i<0){i=0;}
  if (i>a){i=a;}
  return i;
}

void Space::getBoxesFor3DLocations(double p[6],std::vector<Box*>& bp){
  // assume p = [xmin,xmax,ymin,ymax,zmin,zmax]
  int br[6];
  // compute 3D index range that includes 3D location range
  br[0] = location2Index(p[0],"x");   // -x
  br[1] = location2Index(p[1],"x");	//  x
  br[2] = location2Index(p[2],"y");	// -y
  br[3] = location2Index(p[3],"y");	//  y
  br[4] = location2Index(p[4],"z");   // -z
  br[5] = location2Index(p[5],"z");	//  z
  /*	if(!distinguishable(p[0],3687.94) &&
        !distinguishable(p[1],3723.04) &&
        !distinguishable(p[2],7044.62) &&
        !distinguishable(p[3],7063.26) &&
        !distinguishable(p[4],5685.73) &&
        !distinguishable(p[5],5721.23)){
        cout << "\n\n index range ["
        << br[0] << " "
        << br[1] << " "
        << br[2] << " "
        << br[3] << " "
        << br[4] << " "
        << br[5] << "]\n";
        }
        */
  // collect boxes to check
  getBoxesFor3DIndices(br,bp);
}

void Space::getBoxesFor3DIndices(int br[6],std::vector<Box*>& bp)
{
  bp.reserve(VECTOR_RESERVE);
  // assume br = [xmin,xmax,ymin,ymax,zmin,zmax]
  int x0=br[0],x1=br[1],y0=br[2],y1=br[3],z0=br[4];
  std::vector<Box*>::iterator i;
  i = b.begin()+num_space[0]*((z0-1)*num_space[1]+(y0-1));
  for (int z = br[4];z<br[5]+1;z++){
    i += num_space[0]*num_space[1];
    for (int y = br[2];y<br[3]+1;y++){
      i += num_space[0];
      if(x1==x0){
        bp.push_back(*(i+x0));
      } else {
        bp.insert(bp.end(),i+x0,i+x1+1);
      }
    }
    i-=(y1-y0+1)*num_space[0];
  }
  if(bp.empty()==true){
    cout << "br ["
          << br[0] << ","
          << br[1] << ","
          << br[2] << ","
          << br[3] << ","
          << br[4] << ","
          << br[5] << "]\n";
    exit(0);
  }
}

int Space::indices2Index(int x,int y,int z){
  return num_space[0]*(z*num_space[1]+y)+x;
}

void Space::index2Range(int i,double r[2],char *c){
  r[0] = i*SPACE_LENGTH*SCALE;
  r[1] = (i+1)*SPACE_LENGTH*SCALE;
  if (!strcmp(c,"x")){
    r[0] += world[0];
    r[1] += world[0];
  } else if (!strcmp(c,"y")){
    r[0] += world[2];
    r[1] += world[2];
  } else if (!strcmp(c,"z")){
    r[0] += world[4];
    r[1] += world[4];
  }
}

int Space::location2Index(double ss,char *c){
  // since world[] is calculated only once, face limits may be less than world[0,2,4]
  // and greater than world[1,3,5]. Consequently, box_range as computed above may
  // yield negative indices or indices beyond declared range, respectively.
  // Therefore, constrain box_range to known index range.
  // Effectively, border boxes are understood to extend to infinity

  // DEBUG
  /*	bool flag=false;
        if(!distinguishable(ss,7044.62)){flag=true;}
        if(flag){
        cout << "Space::location2Index: "
        << " ss=" << ss
        << ", c=" << c << endl;
        cout << "Space::location2Index: "
        << "world ["
        << world[0] << " "
        << world[1] << " "
        << world[2] << " "
        << world[3] << " "
        << world[4] << " "
        << world[5] << "]\n";
        }*/
  // DEBUG
  int a=0;
  if (!strcmp(c,"x")){
    a = (int) floor((ss-world[0])/(SPACE_LENGTH*SCALE));
    a = screenIndex(a,"x");
  } else if (!strcmp(c,"y")){
    a = (int) floor( (ss-world[2])/(SPACE_LENGTH*SCALE) );
    //		if(flag){
    //			cout << "Space::location2Index: "
    //			<< " a=" << a << endl;
    //		}
    a = screenIndex(a,"y");
    //		if(flag){
    //			cout << "Space::location2Index: "
    //			<< " a=" << a << endl;
    //		}
  } else if (!strcmp(c,"z")){
    a = (int) floor( (ss-world[4])/(SPACE_LENGTH*SCALE) );
    a = screenIndex(a,"z");
  } else {cout << "Received unexpected string.\n"; exit(0);}
  return a;
}

//############################################################################
//############################################################################

class Bit
{
public:
  int *series;
  unsigned int reg;
  vector<unsigned int> adjacent,mask,newbits;
  int num_ints,nv,count,max; //nv = number of vertices
  void init(int);
  void addToAdjacent(std::vector<Vertex*>&);
  void adjacentAndNotMask(void);
  void getNewVertices(std::vector<Vertex*>&,std::vector<Vertex*>&);
  bool get_bit(unsigned int,int);
  void unset_bit(unsigned int &,int);
  void set_bit(unsigned int &,int);
  void clear(void);
};

void Bit::clear(void){
  delete[] series;
}

void Bit::set_bit(unsigned int &i, int bitno)
{
  i |= 1u << bitno;
}

bool Bit::get_bit(unsigned int i, int bitno)
{
  return (i & (1u << bitno)) != 0;
}

void Bit::getNewVertices(std::vector<Vertex*> &t,std::vector<Vertex*> &v){
  std::vector<unsigned int>::iterator j;
  int k=0;
  count=0;
  adjacentAndNotMask();
  mask.assign(adjacent.begin(),adjacent.end());
  adjacent.assign(num_ints,0);
  t.clear();
  ///// read out new bits with 1 /////
  // for each int in newbits
  for (j=newbits.begin();j!=newbits.end();j++){
    reg=*j;
    count = series[k];
    // if int has any bits with 1s
    while (reg!=0) {
      // if bit is active
      if (count<max){
        // if bit is 1, then save new vertex index
        if(reg&1){
          t.push_back(v[count]);
        }
      }
      // bit shift right 1
      reg = reg >> 1;
      count++;
    }
    k++;
  }
}

struct mask_oper
{
  inline int operator()(int i1, int i2) const { return i1 & ~i2; }
};

void Bit::adjacentAndNotMask(void){
  std::transform(adjacent.begin(), adjacent.end(), mask.begin(), newbits.begin(), mask_oper());
}

void Bit::addToAdjacent(std::vector<Vertex*> &b){
  std::vector<Vertex*>::iterator k;
  int i,j;
  // for each adjacent vertex of new vertex
  for(k=b.begin();k!=b.end();k++){
    // compute which integer (byte) the vertex index is located
    i = ((*k)->index-1)/32;		// 0-based
    // compute which bit in byte the vertex index is located
    // -1 accounts for 0 base of vector compared to 1 base of vertex index
    j = (*k)->index-1-32*i;
    // add 1 in corresponding bit location
    set_bit(adjacent[i],j);
  }
}

void Bit::init(int num_vertices){
  nv=num_vertices;
  max=nv+1;
  num_ints = (int) ceil((double)num_vertices/32.0);
  adjacent.assign(num_ints,0);
  mask.assign(num_ints,0);
  newbits.assign(num_ints,0);
  series = new int[num_ints+1];
  for(int i=0;i<(num_ints+1);i++){
    series[i]=32*i;
  }
}

class Monitor {
public:
  int old_rank;
  int used_rank;
  double used_vd;
  double old_vd;
  double old_position[3];
  double new_position[3];
  int touches;
  table_v topN;					// double f -> vertex *v
  table_d old;					// vertex *v -> double se
  hashtable_v_double v_energy;	// vertex* -> double potential energy
  hashtable_e_double e_energy;	// edge* -> double potential energy
  hashtable_v touch_map;			// vertex selection hash map
  hashtable_f ob;					// build old face* hashtable
  hashtable_f nb;					// build new face* hashtable
  hashset_v set_nice_check;		// store vertices for which niceness may have changed
  hashset_v set_nice_changed;		// store vertices for which niceness did change
  hashset_v fs;					// store all vertices requiring a full search for closest point
  //hashset_v ps;					// store all vertices requiring a partial search for closest point
  std::vector<Vertex*> ps;					// store all vertices requiring a partial search for closest point
  //hashset_e ae;					// store adjacent face edges
  std::vector<Edge*> ae;					// store adjacent face edges
  hashtable_e_double e_angle;		// edge* -> double edge angle in radians
  double sum,avg_new,avg_old;		// sliding sum and average of global energy
  int num;						// size of sliding window
  double *window,*begin,*end;		// double pointers to elements of sliding
  // collection of global energies
  v_set refrac_s;					// set of last N vertices moved
  v_list refrac_l;				// set of last N vertices moved

  Vertex* set_seed;

  timeval tim;					// time the main loop
  int interval;					// used to track number of times each vertex is moved
  int gain_period;				// used to schedule gain changes
  double gain_step;				// gain step size
  double max_gain;				// maximum allowed gain
  double ref_gain;				// reference gain
  double gain;					// active gain, all the time
  std::vector<Vertex*> vset;		// unique set of vertices of all faces in box
  Vertex *cv;						// Vertex* of current vertex
  //	hashtable_v gate; 				// list of vertices moved small distance and punished
  // list of vertices moved small distance and punished
  hashtable_v pun_n;				// moved max # times
  hashtable_v pun_int;			// moves result in face intersections
  hashtable_v pun_ang;			// moves result in really small or large angle
  hashtable_v pun_com;			// neither of above cases hold
  tv_iterator tvi; 				// multimap iterator to pair with largest separation error
  double t1; 						// track elapsed time
  int count;						// number of moved vertices in group
  double store;					// square of actual virtual displacement of current vertex
  double pH[3];					// destination of vertex move
  int count_ref;					// to detect if any vertices from set were moved

  Box *bp;						// Box* of worst vertex
  double sep_dis;					// square of separation distance of current vertex before move

  void initAvg(void);
  void freeAvg(void);
  void clearAvg(void);
  void updateAvg(double);
  void initTable(Vertex*,Space&);
  void saveOld(void);
  void updateTopN(Vertex*,double,double,bool);
  bool entryInTopN(Vertex*,double,tv_iterator&);
  void updateOld(Vertex*,double);
  void loadTopN(Container*);
  void validateOld(void);
  void validateTopN(Container*);
  void validateVertex(int,Container&);
  void validateMultimap(void);
  void printVertexSelect(Container&,int);
  void printVertexSelect(Container&,char*);
  void updateTouchMap(Vertex*);
  void updateEnergyMap(Vertex*,double);
  void updateSets(Vertex*,double,bool);
  void cleanTopN(Vertex*,double);
  void prep(Container*);
  double bpf_mean, bpf_min, bpf_max;	// #boxes per face
  double fpb_mean, fpb_min, fpb_max;	// #faces per box
  void getBoxesPerFace(Container&);
  void getFacesPerBox(Space&);
  void printVerticesWithCP(void);
  void initRefrac(void);
  bool Refracted(Vertex*);
  void updateRefractoryWindow(Vertex*);
  int getMaxVertexMoves(void);
  Monitor(void);
  void groupInit(Container&);
  void identifyMeshRegionToUpdate(Space&,Container&);
  bool vertexIsMoveCandidate(Vertex *v,Container&);
  void computeVertex(Vertex*,Container*);
  void updateStatsAndPrint(Container&,int,VTrack&);
  //std::vector<Vertex*>::iterator detectPunishableVertex(std::vector<Vertex*>::iterator,bool,bool,VTrack&,Container&);
  std::vector<Vertex*>::iterator detectPunishableVertex(std::vector<Vertex*>::iterator,bool,bool);
  bool noSetVerticesMoved(void);
  void updateGain(void);
  void writeFilesAndUpdateMaxGain(Container&,int);
  void writeFile(Container&,int);
  //	void punishGate(Vertex*,int);
  void punishN(Vertex*,int);
  void vertexSearchInTopN(Container*,int,std::string);
  bool gainReachedSteadyState(void);
  bool findTopN(Vertex*,tv_iterator&,int&);
  void punishIfOverTouched(Vertex*);
  void writePunished(int);
  void writeIntersected(int,Container*);
  void writeNonnice(int,Container*);

  void screenDisp(Vertex*,double[3]);
  bool isPunished(Vertex*);
  void punishInt(Vertex*,int);
  void punishAng(Vertex*,int);
  void punishCom(Vertex*,int);
};

bool Monitor::isPunished(Vertex *v){
  return (pun_n.find(v)!=pun_n.end() ||
          pun_int.find(v)!=pun_int.end() ||
          pun_com.find(v)!=pun_com.end() || 
          pun_ang.find(v)!=pun_ang.end());
}

void Monitor::writeFile(Container &c,const int group){
  if (WRITE_MESH_EVERY_ITERATION) {
    if(group==7 && count>2749 && count<2801){
      //		printVertexSelect(c,group);
      cout << "Iteration " << group << ": ";
      cout << "Build Mesh after..................";
      cout.flush();
      //c.buildMeshAfter(group);
      c.buildMeshAfter(count);
      cout << "complete.\n";
      cout.flush();
    }
  }
}
void Monitor::writeFilesAndUpdateMaxGain(Container &c,const int group){
  gettimeofday(&tim,NULL);
  double t2=tim.tv_sec+(tim.tv_usec/1000000.0);
  //	c.updateFile(group,count>1,t2-t1,max_gain);
  c.updateFile(group,count>1,t2-t1);

  if(DISABLE_GAIN_SCHEDULING==false)
  {
    //	if (!distinguishable(avg_new,avg_old,ENERGY_EPSILON)){
    //   if(group>NUM_GROUPS/2)
    //   {
    if(max_gain>gain_step){max_gain-=gain_step;}
    if(gain>max_gain){gain=max_gain;}
    /*     }
           else
           {
           max_gain+=gain_step;
           if(gain<max_gain){gain=max_gain;}
           }
           */
    //	}
  }

  /**/	if (WRITE_MESH_EVERY_ITERATION) {
    //		printVertexSelect(c,group);
    cout << "Iteration " << group << ": ";
    cout << "Build Mesh after..................";
    cout.flush();
    c.buildMeshAfter(group);
    cout << "complete.\n";
    cout.flush();
  }
  /*     */
  if (WRITE_DISTANCES_EVERY_ITERATION) {
    cout << "Iteration " << group << ": ";
    cout << "Write closest point distances.....";
    cout.flush();
    c.writeDistances(group);
    cout << "complete.\n";
    cout.flush();
  }
  if(WRITE_PUNISHED_VERTICES_TO_FILE){
    cout << "Iteration " << group << ": ";
    cout << "Write punished vertices...........";
    cout.flush();
    writePunished(group);
    cout << "complete.\n";
    cout.flush();
  }
  if(WRITE_INTERSECTED_FACES_TO_FILE){
    cout << "Iteration " << group << ": ";
    cout << "Write intersected faces...........";
    cout.flush();
    writeIntersected(group,&c);
    cout << "complete.\n";
    cout.flush();
  }
  if(WRITE_NONNICE_VERTICES_TO_FILE){
    cout << "Iteration " << group << ": ";
    cout << "Write nonnice vertices...........";
    cout.flush();
    writeNonnice(group,&c);
    cout << "complete.\n";
    cout.flush();
  }


  /*	//DEBUG
  // close separation distance log file
  c.sepFileOutit();
  // separation distance snapshot log file
  c.sepSnapshotLog((*this),group);
  */	// DEBUG
}

void Monitor::writeNonnice(const int group,Container *c){
  char file[FILENAME_SIZE];
  // create output filename
  if (APPEND_ITERATION_NUMBER_TO_MESH) {
    sprintf(file,"%s%s.%i",OUTPUT_DATA_DIR.c_str(),NONNICE_FILE,group);
  } else {
    sprintf(file,"%s%s",OUTPUT_DATA_DIR.c_str(),NONNICE_FILE);
  }
  // open output file
  std::ofstream newfile (file,std::ios::out);
  if(newfile.is_open()){
    // create vector set
    v_set vs;
    // for each object in container
    for(std::vector<Object*>::iterator i=c->o.begin();i!=c->o.end();i++) {
      // for each nonnice hashtable element in object
      for(vhm_iterator j=(*i)->nice.begin();j!=(*i)->nice.end();j++){
        // add vertex* to set
        vs.insert((*j).first);
      }
    }
    // vector set was collected
    // if cp
    if("FORMAT_NONNICE_VERTICES"=="cp"){
      // print header
      newfile << "x_coordinate y_coordinate z_coordinate "
            << "state_value x_normal y_normal z_normal\n";
      // for each vertex in set
      newfile.precision(12);
      for(vs_iterator i=vs.begin();i!=vs.end();i++){
        // print vertex
        newfile << (*i)->pN[0] << " "
              << (*i)->pN[1] << " "
              << (*i)->pN[2] << " 1 0 0 1\n";
      }
    } else { // else detail
      newfile.precision(12);
      // for each vertex in set
      for(vs_iterator i=vs.begin();i!=vs.end();i++){
        // print vertex
        newfile << "Vertex <obj>" << (*i)->o->name
              << "<ind>" << (*i)->index << " "
              << "["
              << (*i)->pN[0] << " "
              << (*i)->pN[1] << " "
              << (*i)->pN[2] << "]\n";
      }
    }
    newfile.close();
  }
}

void Monitor::writeIntersected(const int group,Container *c){
  char file[FILENAME_SIZE];
  // create output filename
  if (APPEND_ITERATION_NUMBER_TO_MESH) {
    sprintf(file,"%s%s.%i",OUTPUT_DATA_DIR.c_str(),INTERSECTED_FILE,group);
  } else {
    sprintf(file,"%s%s",OUTPUT_DATA_DIR.c_str(),INTERSECTED_FILE);
  }
  // open output file
  std::ofstream newfile (file,std::ios::out);
  if(newfile.is_open()){
    // if cp
    if("FORMAT_INTERSECTED_FACES"=="cp"){
      // create vector set
      v_set vs;
      // for each object in container
      for(std::vector<Object*>::iterator i=c->o.begin();i!=c->o.end();i++) {
        // for each intersected hashtable element in object
        for(htff_iterator j=(*i)->intf.begin();j!=(*i)->intf.end();j++){
          // grab Face *key
          Face *ff=(*j).first;
          // for each vertex of Face* key
          for(int k=0;k<3;k++){
            // add vertex* to set
            vs.insert(ff->v[k]);
          }
          // for each Face* in vector
          std::vector<Face*> *fv=(*j).second;
          for(std::vector<Face*>::iterator k=fv->begin();k!=fv->end();k++){
            // for each vertex of Face*
            for(int nn=0;nn<3;nn++){
              // add vertex* to set
              vs.insert((*k)->v[nn]);
            }
          }
        }
      }
      // vector set was collected
      // print header
      newfile << "x_coordinate y_coordinate z_coordinate "
            << "state_value x_normal y_normal z_normal\n";
      // for each vertex in set
      newfile.precision(12);
      for(vs_iterator i=vs.begin();i!=vs.end();i++){
        // print vertex
        newfile << (*i)->pN[0] << " "
              << (*i)->pN[1] << " "
              << (*i)->pN[2] << " 1 0 0 1\n";
      }

    } else { // else detail
      // create face set
      f_set ffs;
      // for each object in container
      for(std::vector<Object*>::iterator i=c->o.begin();i!=c->o.end();i++) {
        // for each intersected hashtable element in object
        for(htff_iterator j=(*i)->intf.begin();j!=(*i)->intf.end();j++){
          // add Face* to set
          ffs.insert((*j).first);
          // for each Face* in vector
          std::vector<Face*> *fv=(*j).second;
          for(std::vector<Face*>::iterator k=fv->begin();k!=fv->end();k++){
            // add Face* to set
            ffs.insert(*k);
          }
        }
      }
      // face set was collected
      newfile.precision(12);
      // for each face in set
      for(fs_iterator i=ffs.begin();i!=ffs.end();i++){
        // print face
        newfile << "Face <obj>" << (*i)->v[0]->o->name
              << "<ind>" << (*i)->index << endl
              << "[v0 "
              << (*i)->v[0]->index << " "
              << (*i)->v[0]->pN[0] << " "
              << (*i)->v[0]->pN[1] << " "
              << (*i)->v[0]->pN[2] << "]\n"
              << "[v1 "
              << (*i)->v[1]->index << " "
              << (*i)->v[1]->pN[0] << " "
              << (*i)->v[1]->pN[1] << " "
              << (*i)->v[1]->pN[2] << "]\n"
              << "[v2 "
              << (*i)->v[2]->index << " "
              << (*i)->v[2]->pN[0] << " "
              << (*i)->v[2]->pN[1] << " "
              << (*i)->v[2]->pN[2] << "]\n";
      }
    }
    newfile.close();
  }
}

void Monitor::writePunished(const int group){
  char file[FILENAME_SIZE];
  // create output filename
  if (APPEND_ITERATION_NUMBER_TO_MESH) {
    sprintf(file,"%s%s.%i",OUTPUT_DATA_DIR.c_str(),PUNISHED_FILE,group);
  } else {
    sprintf(file,"%s%s",OUTPUT_DATA_DIR.c_str(),PUNISHED_FILE);
  }
  // open output file
  std::ofstream newfile (file,std::ios::out);
  if(newfile.is_open()){
    newfile.precision(12);
    if(pun_n.empty()==false){
      // for each vertex moved too many time
      newfile << "# vertices punished for moving MAX_TOUCHES="
            << MAX_TOUCHES << "times in group=" << group << endl;
      for(vhm_iterator i=pun_n.begin();i!=pun_n.end();i++) {
        Vertex *vv=(*i).first;
        newfile << "Vertex <obj>" << vv->o->name
              << "<ind>" << vv->index << " "
              << "["
              << vv->pN[0] << " "
              << vv->pN[1] << " "
              << vv->pN[2] << "]\n";
      }
    }
    if(pun_int.empty()==false){
      // for each punished vertex
      newfile << "# vertices punished for causing intersecting faces"
            << " in group=" << group << endl;
      for(vhm_iterator i=pun_int.begin();i!=pun_int.end();i++) {
        Vertex *vv=(*i).first;
        newfile << "Vertex <obj>" << vv->o->name
              << "<ind>" << vv->index << " "
              << "["
              << vv->pN[0] << " "
              << vv->pN[1] << " "
              << vv->pN[2] << "]\n";
      }
    }
    if(pun_ang.empty()==false){
      // for each punished vertex
      newfile << "# vertices punished for causing very small or very large edge angles"
            << " in group=" << group << endl;
      for(vhm_iterator i=pun_ang.begin();i!=pun_ang.end();i++) {
        Vertex *vv=(*i).first;
        newfile << "Vertex <obj>" << vv->o->name
              << "<ind>" << vv->index << " "
              << "["
              << vv->pN[0] << " "
              << vv->pN[1] << " "
              << vv->pN[2] << "]\n";
      }
    }
    if(pun_com.empty()==false){
      // for each punished vertex
      newfile << "# vertices whose punishment is not "
            << "attributable to either of the above categories" << endl;
      for(vhm_iterator i=pun_com.begin();i!=pun_com.end();i++) {
        Vertex *vv=(*i).first;
        newfile << "Vertex <obj>" << vv->o->name
              << "<ind>" << vv->index << " "
              << "["
              << vv->pN[0] << " "
              << vv->pN[1] << " "
              << vv->pN[2] << "]\n";
      }
    }
    newfile.close();
  } else {
    cout << "Monitor::writePunished: Error.\n"
          << "Could not open " << file << " for writing punished vertices.\n";
  }
}
void Monitor::updateGain(void){
  // if gain changeable?
  if(gain_period==0){
    // let cc equal the max number of moves
    // made by a single vertex in the last gain period
    //		int cc = getMaxVertexMoves();
    //		// if oscillating, i.e. max number of vertex moves > upper threshold
    //		// and gain is larger than gain step size
    //		if(cc>UPPER_MOVES_THRESHOLD && gain>gain_step){
    //
    // if change in mean energy < epsilon
    // and gain is larger than gain step size
    if( gainReachedSteadyState()==true && gain>gain_step){
      //			cout << "cc " << cc
      //			<< ", gain " << gain
      //			<< ", gain_step " << gain_step
      //			<< ", max_gain-gain_step " << max_gain-gain_step
      //			<< endl;
      // decrease gain
      gain -= gain_step;
      // reset gain period
      gain_period = REFRACTORY_PERIOD;
      // update reference gain
      ref_gain = gain;
    } 
    // else if system not active enough
    // i.e. max number of vertex moves < lower threshold
    // and gain is small enough to be incremented without exceeding max_gain
    /*		else if (cc<LOWER_MOVES_THRESHOLD && gain<=(max_gain-gain_step)){
    // increase gain
    gain += gain_step;
    // reset gain period
    gain_period = REFRACTORY_PERIOD;
    // update reference gain
    ref_gain=gain;
    }
    */
  }
  }

  bool Monitor::noSetVerticesMoved(void){
    return count_ref==count;
  }

  void Monitor::screenDisp(Vertex *ccvv,double PT[3]){
    // square of displacement
    double d[3]={PT[0]-ccvv->pN[0],PT[1]-ccvv->pN[1],PT[2]-ccvv->pN[2]};
    //double sqdisp = d[0]*d[0]+d[1]*d[1]+d[2]*d[2];
    double disp = sqrt(d[0]*d[0]+d[1]*d[1]+d[2]*d[2]);
    // compute max allowed displacement
    double min_edge_length = 1e30;
    // for each adjacent face
    for(std::vector<Face*>::iterator i=ccvv->f.begin();i!=ccvv->f.end();i++)
    {
      // find shortest edge
      for(int j=0;j<3;j++)
      {
        if((*i)->e[j]->l < min_edge_length)
        {
          min_edge_length = (*i)->e[j]->l;
        }
      }
    }
    double MAX_ACTUAL_DISPL = min_edge_length*MAX_ACTUAL_DISPL_FRACTION;
    // if too big
    if(disp>MAX_ACTUAL_DISPL){
      for(int k=0;k<3;k++){
        PT[k] = ccvv->pN[k]+MAX_ACTUAL_DISPL*d[k]/disp;
      }
    }
    /*	// if too big
        if(sqdisp>MAX_ACTUAL_DISPL_SQ){
        double disp=sqrt(sqdisp);
        double newgain=sqrt(MAX_ACTUAL_DISPL_SQ);
        for(int k=0;k<3;k++){
        PT[k] = ccvv->pN[k]+newgain*d[k]/disp;
        }
        }
        */
  }

  void Monitor::computeVertex(Vertex *v,Container *c){
    cv=v;
    // compute new holding position coordinates (x,y,z)
    cv->computeNewCoords(c,pH,gain);
    // impose max displacement policy
    screenDisp(cv,pH);
    // displacement along x,y,z
    double a=pH[0]-cv->pN[0];
    double b=pH[1]-cv->pN[1];
    double d=pH[2]-cv->pN[2];
    // compute and store squared displacement 
    store=a*a+b*b+d*d;
    // compute and store squared separation distance 
    sep_dis = cv->getSqSepDist(c);

    /*	// DEBUG
        if(v==set_seed){
        old_rank=0;
        old_vd=0;
        tv_iterator	target;
        if(findTopN(v,target,old_rank)==true){
        old_vd = sqrt((*target).first);
        } else {
        cout << "vertex not found in topN!" << endl;
        v->printVertex(v->o->name);cout << endl;
        exit(0);
        }

        old_position[0] = cv->pN[0];
        old_position[1] = cv->pN[1];
        old_position[2] = cv->pN[2];
        new_position[0] = pH[0];
        new_position[1] = pH[1];
        new_position[2] = pH[2];
        }
        */	// DEBUG

  }


  bool Monitor::vertexIsMoveCandidate(Vertex *v,Container &c){
    // if vertex is candidate, i.e. closest point was found for vertex
    // and vertex not found in refractory list of vertices
    // i.e. list of last REFRACTORY_PERIOD (e.g. 1000) vertices moved
    // and vertex is not in gate, i.e. failed to move and punished for extent of w group
    return (v->cl!=NULL) && (isPunished(v)==false) && !Refracted(v)
          && (binary_search(c.frozen.begin(),c.frozen.end(),v)==false);
    //	return (v->cl!=NULL) && (gate.find(v)==gate.end()) && !Refracted(v);
  }



  void Monitor::identifyMeshRegionToUpdate(Space &s,Container &c){
    set_seed=NULL;
    if(MOVE_NONNICE_AND_INTERSECTED_FIRST){
      // try to pick nonnice vertex as set seed
      // for each object
      for(std::vector<Object*>::iterator i=c.o.begin();i!=c.o.end();i++){
        // for each element of nice hashmap
        for(vhm_iterator j=(*i)->nice.begin();j!=(*i)->nice.end();j++){
          Vertex *vv=(*j).first;
          // if vertex is not punished and not refracted
          // and not frozen and has a closest face
          if ( isPunished(vv)==false && Refracted(vv)==false && vv->cl!=NULL
               && (binary_search(c.frozen.begin(),c.frozen.end(),vv)==false))
          {
            set_seed = vv;break;
          }
        }
        if(set_seed!=NULL){break;}
      }

      // else try to pick intersected face vertex as set seed
      if(set_seed==NULL){
        // for each object
        for(std::vector<Object*>::iterator i=c.o.begin();i!=c.o.end();i++){
          // for each element of intf hashmap
          for(htff_iterator j=(*i)->intf.begin();j!=(*i)->intf.end();j++){
            Face *ff=(*j).first;
            // for each vertex of face
            for(int k=0;k<3;k++){
              Vertex *vv=ff->v[k];
              // if vertex is not punished and not refracted
              // and not frozen and has a closest face
              if (isPunished(vv)==false && Refracted(vv)==false && vv->cl!=NULL
                  && (binary_search(c.frozen.begin(),c.frozen.end(),vv)==false))
              {
                set_seed = vv;break;
              }
            }
            if(set_seed!=NULL){break;}
          }
          if(set_seed!=NULL){break;}
        }
      }
    }

    // else pick vertex from sorted virtual displacement list
    if(set_seed==NULL){
      // starting at current location of tvi in topN
      // find vertex with largest virtual displacement that is not punished
      // and not refracted and not frozen and has a closest face
      while(isPunished((*tvi).second)==true || Refracted((*tvi).second)==true
            || (binary_search(c.frozen.begin(),c.frozen.end(),(*tvi).second)==true)
            || (*tvi).second->cl==NULL){
        tvi++;
      }
      set_seed = (*tvi).second;
    }

    // DEBUG
    /*
        if(set_seed!=NULL){
        cout << "\n>>>>> SET_SEED VERTEX: " << set_seed->o->name << "->" << set_seed->index << endl;
        } else {
        cout << "\n>>>>> SET_SEED VERTEX: NULL" << endl;
        }*/
    // DEBUG

    // grab Box* of worst vertex
    bp=s.b[s.indices2Index(
                           //		s.location2Index((*tvi).second->pN[0],"x"),
                           //		s.location2Index((*tvi).second->pN[1],"y"),
                           //		s.location2Index((*tvi).second->pN[2],"z"))];
          s.location2Index(set_seed->pN[0],"x"),
          s.location2Index(set_seed->pN[1],"y"),
          s.location2Index(set_seed->pN[2],"z"))];
    // grab unique set of vertices of all faces in box
    std::vector<Vertex*> bin;
    for(std::vector<Face*>::iterator z=bp->f.begin();z!=bp->f.end();z++){
      for(int k=0;k<3;k++){
        if((*z)->v[k]!=set_seed){
          bin.push_back((*z)->v[k]);
        }
      }
    }
    // sort and keep unique Vertex*
    sort(bin.begin(),bin.end());
    std::vector<Vertex*>::iterator new_end = unique(bin.begin(),bin.end());
    bin.assign(bin.begin(),new_end);
    // build vset
    vset.clear();	
    vset.push_back(set_seed);
    for(std::vector<Vertex*>::iterator i=bin.begin();i!=bin.end();i++){
      vset.push_back(*i);
    }
    // to detect if any vertices from set were moved
    count_ref = count;
  }

  Monitor::Monitor(void){
    // used to track number of times each vertex is moved
    interval=0;
    // used to schedule gain changes
    gain_period=0;
    // inital max gain
    max_gain = TIME_STEP/DAMPING;
    // initialize gain
    gain = max_gain;
    // gain step size
    gain_step = max_gain/NUM_GROUPS;
    // Box* of worst vertex
    bp=NULL;
    // Vertex* of current vertex
    cv=NULL;
    // reserve capacity for vectors
    ae.reserve(VECTOR_RESERVE);
    ps.reserve(VECTOR_RESERVE);
  }
  /*
      Monitor::Monitor(void){
// used to track number of times each vertex is moved
interval=0;
// used to schedule gain changes
gain_period=0;
// inital max gain
max_gain = TIME_STEP/DAMPING;
// gain step size
gain_step = max_gain/NUM_GROUPS*2;
// initialize gain
gain = gain_step;
max_gain=gain;
// Box* of worst vertex
bp=NULL;
// Vertex* of current vertex
cv=NULL;
}
*/
void Monitor::groupInit(Container &c){
  //	if(group-1){
  // load topN, save old, and initialize avg
  prep(&c);
  //	}
  // instantiate list of vertices moved small distance (gate)
  // and instantiate refractory period list
  initRefrac();
  // clear punished vertex lists
  pun_n.clear();
  pun_int.clear();
  pun_com.clear();
  pun_ang.clear();
  // initialize reference gain
  ref_gain = gain;
  // intialize sliding window average of global energy
  clearAvg();
  updateAvg(c.energy);
  c.clear();
  // point multimap iterator to pair with largest separation error
  tvi = topN.begin();

  // initialize variable used to track how many times each vertex was moved
  touches = 0;
  touch_map.clear();
  // track elapsed time
  gettimeofday(&tim,NULL);
  t1=tim.tv_sec+(tim.tv_usec/1000000.0);
  // initialize count of vertices moved in group
  count = 1;
}

bool Face::faceInVector(Face *ff){
  // get this face's intersecting faces vector
  std::vector<Face*> *ifv=this->getIntersectingFaces();
  // look for face f in this face's intersecting faces vector
  std::vector<Face*>::iterator i=find((*ifv).begin(),(*ifv).end(),ff);
  // return true if found
  // return false if no found
  return (i!=(*ifv).end());
}

void Face::addFaceToVector(Face* f,Container *c){
  // if this face not in intf, then add new element to intf
  if(!faceInTable_intf()){addFaceToTable_intf();}
  // face f not already in vector
  if(faceInVector(f)==false){
    // add face f to this face's intersecting faces vector
    (*v[0]->o->intf[this]).push_back(f);
  }
  // update global intersections count
  if (v[0]->o==f->v[0]->o){c->si++;}
  c->ti++;
}

void Face::removeFaceFromVector(Face *f,Container *c){
  // update si for removal of both adjacent and intersecting faces
  if(v[0]->o==f->v[0]->o){c->si-=1;}
  // update ti for removal of both adjacent and intersecting faces
  c->ti-=1;
  // if this face has intersecting faces in hashmap intf
  if(faceInTable_intf()==true){
    // get this face's intersecting faces vector
    std::vector<Face*> *ifv=this->getIntersectingFaces();
    // look for face f in this face's intersecting faces vector
    std::vector<Face*>::iterator i=find((*ifv).begin(),(*ifv).end(),f);
    // if found
    if (i!=(*ifv).end()){
      // remove face f from this face's intersecting faces vector
      (*ifv).erase(remove((*ifv).begin(),(*ifv).end(),f),(*ifv).end());
      // if intersecting face vector is now empty
      if (noMoreIntersectingFaces()){
        // then remove intersecting face vector from hashtable
        clearFaceFromTable_intf(c);
        // then update face intersect force, i.e. set force to zero
        clearFaceFromTable_iv();
      }
    } else {
      /*			cout << "\nFace::removeFaceFromVector: "
                                << "Face " << f->v[0]->o->name
                                << "->" << f->index
                                << " NOT found in intersecting face vector of Face "
                                << v[0]->o->name << "->" << index
                                << " intersecting face vector.\n";
                                */
    }
  } else {
    /*		cout << "\nFace::removeFaceFromVector: "
                << "Face " << v[0]->o->name
                << "->" << index
                << " NOT found in intf table.\n";
                */
  }
}

void Face::clearFaceFromTable_intf(Container *c){
  // if this face is in intf table 
  if(faceInTable_intf()==true){
    // get this face's intersecting faces vector
    std::vector<Face*> *ifv=this->getIntersectingFaces();
    // if intersecting faces vector is not empty
    if((*ifv).empty()==false){
      // for each face in intersection vector
      for(std::vector<Face*>::iterator i=(*ifv).begin();i!=(*ifv).end();i++){
        /*				// DEBUG
                                        cout << "Container::clearFaceFromTable_intf: "
                                        << "adjface=" << index
                                        << ", intface=" << (*i)->index
                                        << endl;
                                        */				// DEBUG
        // update si for removal of both adjacent and intersecting faces
        if(v[0]->o==(*i)->v[0]->o){c->si--;}
        // update ti for removal of both adjacent and intersecting faces
        c->ti--;
      }
      (*ifv).clear();
    }
    // delete vector<face*>*
    delete v[0]->o->intf[this];
    // remove element from table
    v[0]->o->intf.erase(this);
  }
}

void Vertex::setVertexNiceness(int newval,Container *c){
  int oldval=0;
  vhm_iterator i = o->nice.find(this);
  // if vertex in hashtable
  if(i!=o->nice.end()){oldval=o->nice[this];}
  // process
  // not old==new needs no action
  switch (oldval+newval*3) {
    case 1 :
      // old==1 + new==0*3
      // remove this vertex from hashtable
      o->nice.erase(i);
      c->nonnice--;
      break;

    case 2 :
      // old==2 + new==0*3
      o->nice.erase(i);
      c->nonnice--;
      c->s_nonnice--;
      break;

    case 3 :
      // old==0 + new==1*3
      o->nice[this]=newval;
      c->nonnice++;
      break;

    case 5 :
      // old==2 + new==1*3
      o->nice[this]=newval;
      c->s_nonnice--;
      break;

    case 6 :
      // old==0 + new==2*3
      o->nice[this]=newval;
      c->nonnice++;
      c->s_nonnice++;
      break;

    case 7 :
      // old==1 + new==2*3
      o->nice[this]=newval;
      c->s_nonnice++;
      break;

  }	
}


//############################################################################
//############################################################################

class VTrack{
public:
  Vertex *vv;
  // ..[0]==original, ..[1]==new
  double p[2][3];
  double topN[2];
  double vd[2];
  double sepdis[2];
  double cp[2][3];
  void clear(void);
  bool isGood(void);
  void print(void);
  void printBad(void);
  void addNewSepDis(double);
  double getNewSepDis(void);
  void addNewVD(double);
  double getNewVD(void);
  void addNewTopN(double);
  double getNewTopN(void);
  void addOrigSepDis(double);
  double getOrigSepDis(void);
  void addOrigVD(double);
  double getOrigVD(void);
  void addOrigTopN(double);
  double getOrigTopN(void);
  void addOrigP(double[3]);
  void getOrigP(double[3]);
  void getNewP(double[3]);
  void addNewP(double[3]);
  void addOrigCP(double[3]);
  void getOrigCP(double[3]);
  void getNewCP(double[3]);
  void addNewCP(double[3]);
  void addVertex(Vertex*);
  Vertex* getVertex(void);
  std::vector<int> trouble;
};

void VTrack::printBad(void){
  if (vv==NULL){cout << "vertex* is NULL.\n";}
  if (p[0][0]<0){cout << "orig pos not set.\n";}
  if (p[1][0]<0){cout << "new pos not set.\n";}
  if (cp[0][0]<0){cout << "orig clos pos not set.\n";}
  if (cp[1][0]<0){cout << "new clos pos not set.\n";}
  if (topN[0]<0){cout << "orig topN not set.\n";}
  if (topN[1]<0){cout << "new topN not set.\n";}
  if (vd[0]<0){cout << "orig vd not set.\n";}
  if (vd[1]<0){cout << "new vd not set.\n";}
  if (sepdis[0]<0){cout << "orig sepdis not set.\n";}
  if (sepdis[1]<0){cout << "new sepdis not set.\n";}
  cout << "p[0][0] = " << p[0][0] << endl;
  cout << "p[1][0] = " << p[1][0] << endl;
  cout << "cp[0][0] = " << cp[0][0] << endl;
  cout << "cp[1][0] = " << cp[1][0] << endl;
  cout << "topN[0] = " << topN[0] << endl;
  cout << "topN[1] = " << topN[1] << endl;
  cout << "vd[0] = " << vd[0] << endl;
  cout << "vd[1] = " << vd[1] << endl;
  cout << "sepdis[0] = " << sepdis[0] << endl;
  cout << "sepdis[1] = " << sepdis[1] << endl;
}

void VTrack::print(void){
  cout << "pos ["
        << p[0][0] << " "
        << p[0][1] << " "
        << p[0][2] << "]->["
        << p[1][0] << " "
        << p[1][1] << " "
        << p[1][2] << "]\n"
        << "closest pos ["
        << cp[0][0] << " "
        << cp[0][1] << " "
        << cp[0][2] << "]->["
        << cp[1][0] << " "
        << cp[1][1] << " "
        << cp[1][2] << "]\n"
        << "topN rank ["
        << topN[0] << "]->["
        << topN[1] << "], "
        << "vd ["
        << vd[0] << "]->["
        << vd[1] << "], "
        << "sepdis ["
        << sepdis[0] << "]->["
        << sepdis[1] << "]\n";
}

Vertex* VTrack::getVertex(void){
  return vv;
}

void VTrack::addVertex(Vertex *uu){
  vv=uu;
}

bool VTrack::isGood(void){
  return	vv!=NULL &&
        (p[0][0]>0) &&
        (p[1][0]>0) &&
        (cp[0][0]>0) &&
        (cp[1][0]>0) &&
        (topN[0]>0) &&
        (topN[1]>0) &&
        (vd[0]>=0) &&
        (vd[1]>=0) &&
        (sepdis[0]>=0) &&
        (sepdis[1]>0);
}

void VTrack::addNewSepDis(double i){
  sepdis[1]=i;
}

double VTrack::getNewSepDis(void){
  return sepdis[1];
}

void VTrack::addNewVD(double i){
  vd[1]=i;
}

double VTrack::getNewVD(void){
  return vd[1];
}

void VTrack::addNewTopN(double i){
  topN[1]=i;
}

double VTrack::getNewTopN(void){
  return topN[1];
}

void VTrack::addOrigSepDis(double i){
  sepdis[0]=i;
}

double VTrack::getOrigSepDis(void){
  return sepdis[0];
}

void VTrack::addOrigVD(double i){
  vd[0]=i;
}

double VTrack::getOrigVD(void){
  return vd[0];
}

void VTrack::addOrigTopN(double i){
  topN[0]=i;
}

double VTrack::getOrigTopN(void){
  return topN[0];
}

void VTrack::clear(void){
  for(int i=0;i<2;i++){
    p[i][0]=p[i][1]=p[i][2]=-1.0;
    cp[i][0]=cp[i][1]=cp[i][2]=-1.0;
    topN[i]=vd[i]=sepdis[i]=-1.0;
  }
}

void VTrack::addOrigP(double q[3]){
  p[0][0]=q[0];
  p[0][1]=q[1];
  p[0][2]=q[2];
}

void VTrack::getOrigP(double q[3]){
  q[0]=p[0][0];
  q[1]=p[0][1];
  q[2]=p[0][2];
}

void VTrack::getNewP(double q[3]){
  q[0]=p[1][0];
  q[1]=p[1][1];
  q[2]=p[1][2];
}

void VTrack::addNewP(double q[3]){
  p[1][0]=q[0];
  p[1][1]=q[1];
  p[1][2]=q[2];
}

void VTrack::addOrigCP(double q[3]){
  cp[0][0]=q[0];
  cp[0][1]=q[1];
  cp[0][2]=q[2];
}

void VTrack::getOrigCP(double q[3]){
  q[0]=cp[0][0];
  q[1]=cp[0][1];
  q[2]=cp[0][2];
}

void VTrack::getNewCP(double q[3]){
  q[0]=cp[1][0];
  q[1]=cp[1][1];
  q[2]=cp[1][2];
}

void VTrack::addNewCP(double q[3]){
  cp[1][0]=q[0];
  cp[1][1]=q[1];
  cp[1][2]=q[2];
}

std::vector<Vertex*>::iterator Monitor::detectPunishableVertex(std::vector<Vertex*>::iterator i,bool int_flag,bool ang_flag){
  if(store<MIN_DISPLACEMENT_SQ*SCALE*SCALE){
    /**/	// DEBUG
    //		pod.addNewP(stats.cv->pN);
    //		pod.addNewSepDis(sqrt(stats.cv->getSqSepDist(&c)));
    /*
        pod.addNewP(pod.p[0]);
        pod.addNewSepDis(pod.sepdis[0]);
        if(true){
        int rrank;
        tv_iterator ttt;
        findTopN(cv,ttt,rrank);
    // get closest point
    double pC[3];
    c.computePC(cv->cl,cv,pC);
    // set
    pod.addNewCP(pC);
    pod.addNewVD(sqrt((*ttt).first));
    pod.addNewTopN(rrank);
    }
    if(pod.isGood()==true){pod.print();}
    else { pod.printBad();exit(0);}
    */
    /**/	// DEBUG
    cout << "\nvertex: " << cv->o->name << "->" << cv->index << " was punished.\n"
          << "virtual displacement sqd(" << store << ")"
          << " < threshold(" << MIN_DISPLACEMENT_SQ*SCALE*SCALE << ")\n";
    cv->printVertex(cv->o->name);
    cout << endl << endl;
    // vertex is effectively immovable, so punish vertex
    if(int_flag==true){	punishInt(cv,static_cast<int>(GROUP_SIZE));}
    else if(ang_flag==true){	punishAng(cv,static_cast<int>(GROUP_SIZE));}
    else {	punishCom(cv,static_cast<int>(GROUP_SIZE));}
    // reset gain
    gain=ref_gain;
    // move on to next vertex in set
    i++;
  } else {
    // try halving gain to halve attempted virtual displacement and try again
    gain=gain/2.0;
  }
  return i;
}

//void Monitor::updateStatsAndPrint(Container &c,const int group){
void Monitor::updateStatsAndPrint(Container &c,const int group,VTrack &pod){
  // automatically strict face intersection prevention
  if(STRICT_FACE_INTERSECTION_PREVENTION==false){
    if(DETECT_COMPLETE_SEPARATION==true){
      if(c.ti==0){STRICT_FACE_INTERSECTION_PREVENTION=true;}
    }
  }
  // update gain_period
  if(gain_period>0){gain_period--;}
  // increment number of times current vertex has been moved
  updateTouchMap(cv);
  // punish vertex if moved too many time
  punishIfOverTouched(cv);
  // update sliding energy average
  updateAvg(c.energy);
  // print number of times each vertex has been moved
  if(!WRITE_MESH_EVERY_ITERATION && touches==1000){
    touches = 0;
    printVertexSelect(c,++interval);
    touch_map.clear();
  }
  // update container stats
  c.updateStats(sqrt(store));
  // update list of vertices in refractory period
  updateRefractoryWindow(cv);

  if(false)
  {cout << pod.trouble.size() << endl;}

  if(ENABLE_VTRACK==true)
  {
    //if(binary_search(pod.trouble.begin(),pod.trouble.end(),cv->index)==true && cv->o->name=="g_clipped")
    if(
       ((cv->index)==330 && cv->o->name=="g_clipped") ||
       ((cv->index)==158 && cv->o->name=="d000_clipped")
      )
    {
      cout << "TROUBLE VERTEX MOVED!\n";
    }
  }
  // print to stdout
  if(!(count++%10000)){
    //	if(count++){
    char name[1024];
    sprintf(name,"count %d, ",count-1);
    cout << name;
    sprintf(name,"group %d, ",group);
    cout << name;
    cout.width(45);
    sprintf(name,"%s->%d, ",cv->o->name.c_str(),cv->index);
    cout << left << name;
    cout.width(24);
    sprintf(name,"vd %.12g, ",sqrt(store));
    cout << left << name;
    cout.width(24);
    sprintf(name,"e %.12g, ",c.energy);
    cout << left << name;
    cout.width(28);
    sprintf(name,"e_avg %.12g, ",avg_new);
    cout << left << name;
    cout.width(12);
    sprintf(name,"gain %.6g, ",gain);
    cout << left << name;
    //        cout.width(14);
    //        sprintf(name,"max_gain %.6g, ",max_gain);
    //        cout << left << name;
    //        sprintf(name,"min edge angle %.6g, ",c.min_edge_angle);
    //        cout << left << name;
    if(PRINT_FLAG){
      sprintf(name,"bpf (%.6g,%.6g,%.6g), ",bpf_min,bpf_mean,bpf_max);
      cout << left << name;
      sprintf(name,"fpb (%.6g,%.6g,%.6g), ",fpb_min,fpb_mean,fpb_max);
      cout << left << name;
    }
    cout << endl << endl;
  }

  /*	// DEBUG
        bool dec=false;
  //	if(count && group==3){dec=true;}
  //	if(count++ && cv==set_seed){dec=true;}
  if(dec){
  char name[1024];
  sprintf(name,"count %d, ",count-1);
  cout << name;
  sprintf(name,"iter %d, ",group);
  cout << name;
  if(cv==set_seed){cout << " seed, ";}
  else {cout << "member, ";}
  cout.width(45);
  sprintf(name,"%s->%d, ",cv->o->name.c_str(),cv->index);
  cout << left << name;

  tv_iterator	target;
  int rank=0;
  if(findTopN(cv,target,rank)==true){
  cout.width(28);
  sprintf(name,"rank [%d->%d], ",old_rank,rank);
  cout << left << name;
  cout.width(38);
  sprintf(name,"vd [%.12g->%.12g], ",old_vd,sqrt((*target).first));
  cout << left << name;
  } else {
  cout << "vertex not found in topN!" << endl;
  cv->printVertex(cv->o->name);cout << endl;
  exit(0);
  }
  cout.width(24);
  sprintf(name,"actual disp %.12g, ",sqrt(store));
  cout << left << name << endl;

  cout.width(88);
  sprintf(name,"pos [%.12g,%.12g,%.12g]->[%.12g,%.12g,%.12g], "
  ,old_position[0],old_position[1],old_position[2]
  ,cv->pN[0],cv->pN[1],cv->pN[2]);
  cout << left << name;

  // get closest point
  double pC[3];
  c.computePC(cv->cl,cv,pC);
  cout.width(48);
  sprintf(name,"new clos pt [%.12g,%.12g,%.12g]",pC[0],pC[1],pC[2]);
  cout << left << name;

  cout << endl;
  }*/
  // update iterator
  gain=ref_gain;
  /*	if(group==3){
  // print debug data to file
  char file[256];	
  sprintf(file,"%s%s_%d.dat",OUTPUT_DATA_DIR,cv->o->name.c_str(),cv->index);
  std::ofstream File(file,std::ofstream::app);
  if (File.fail()) // if stream cannot be opened
  { cout << "Can't open " << file ; exit(1); }
  File.precision(12);
  File << endl;
  File << "group " << group << endl;
  File << "count " << count-1 << endl;
  if(cv==set_seed){File << "seed\n";}
  else {File << "member\n";}
  File << "old topN rank " << old_rank << endl;
  File << "old virtual displacement " << old_vd << endl;
  tv_iterator	target;
  int rank=0;
  if(findTopN(cv,target,rank)==true){
  File << "new topN rank " << rank << endl;
  File << "new virtual displacement "
  << sqrt((*target).first) << endl;
  } else {
  cout << "vertex not found in topN!" << endl;
  cv->printVertex(cv->o->name);cout << endl;
  exit(0);
  }
  File << "actual displacement " << sqrt(store) << endl;
  File << "old position ["
  << old_position[0] << ", "
  << old_position[1] << ", "
  << old_position[2] << "]\n";
  File << "new position [" << cv->pN[0] << ", "
  << cv->pN[1] << ", " << cv->pN[2] << "]\n";
  // get closest point
  double pC[3];
  c.computePC(cv->cl,cv,pC);
  File << "new closest point [" << pC[0] << ", "
  << pC[1] << ", " << pC[2] << "]\n";
  // close files
  File.close();
  }
  */	// DEBUG

  }

  //############################################################################
  //############################################################################
