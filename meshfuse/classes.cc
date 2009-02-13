class Object;
class Face;
class Vertex;
class Edge;
class Distance;

struct lts
{
  bool operator()(const std::string s1, const std::string s2) const
  {
    return s1 < s2;
  }
};

typedef std::map<std::string,Edge*,lts,std::allocator<Edge*> > hashtable_t;
typedef std::map<std::string,Edge*,lts,std::allocator<Edge*> >::iterator ht_iterator;
typedef std::vector<Vertex*>		vec_v;
typedef std::vector<Vertex*>::iterator	v_iterator;
typedef std::vector<Face*>		vec_f;
typedef std::vector<Face*>::iterator	f_iterator;
typedef std::list<Edge*>		list_e;
typedef std::list<Edge*>::iterator	e_iterator;
typedef std::vector<Distance>		vec_d;
typedef std::vector<Distance>::iterator	d_iterator;

class Object
{
public:
  std::vector<Vertex*> v;		// container of pointers to all vertices in object
  std::vector<Face*> f;		// container of pointers to all faces in object
  list_e e;		// container of pointers to all edges in object
  hashtable_t hm;               // create map for finding edges
  int num_digits;               // determine number of digits in largest vertex index
  vec_v free_vertices;          // vertices that are on border
  vec_d distances;
  Object(std::string const &);
  void scanFile(const char*);
  void createEdges(void);
  void createEdge(Face*,Vertex*,Vertex*);
  void buildEdge(Face*,Vertex*,Vertex*);
  int setNumDigits(void);
  void gatherFreeVertices(void);
  double findLongestEdge(void);
  void computeDistances(const double,double);
  void fixVertices(Vertex*);
  void fixFaces(Vertex*,Vertex*);
  void updateEdges(Vertex*,Vertex*);
  void purgeEdgeFromMap(Edge*);
  void purgeEdgeFromFaces(Edge*);
  void updateDistances(void);
  void printVerticesFaces(void);
  void purgeMap(void);
  void purgeEdges(void);
};

Object::Object(std::string const &file)
{
  fprintf(stderr,"\n\nBuilding vertices and faces.....");
  fflush(stderr);
  scanFile(file.c_str());
  sort(v.begin(),v.end());
  sort(f.begin(),f.end());
  fprintf(stderr,"complete.\n");
  fflush(stderr);
  fprintf(stderr,"Building edges..................");
  fflush(stderr);
  createEdges();
  fprintf(stderr,"complete.\n");
  fflush(stderr);
}

class Vertex
{
public:
  double p[3];
  int index;
  Vertex(char*);
  void printVertex(void);
};

void Vertex::printVertex(void)
{
  cout.precision(12);
  cout << "Vertex " << index << " "
        << p[0] << " "
        << p[1] << " "
        << p[2] << "\n";
  cout.flush();
}

Vertex::Vertex(char *triplet)
{
  char val[80];
  char *eptr;
  int i;

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
    p[0]=p[1]=p[2]=0;
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
  p[0] = strtod(val,&eptr);
  if (val==eptr)
  {
    p[0]=p[1]=p[2]=0;
    printf("Error in reading vertex\n");
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
  p[1] = strtod(val,&eptr);
  if (val==eptr)
  {
    p[0]=p[1]=p[2]=0;
    printf("Error in reading vertex\n");
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
  p[2] = strtod(val,&eptr);
  if (val==eptr)
  {
    p[0]=p[1]=p[2]=0;
    printf("Error in reading vertex\n");
    return;
  }
}

class Face
{
public:
  int index;	// Face index
  Vertex *v[3];	// vertex indices
  Edge *e[3];	// edge pointers
  Face(char*,std::vector<Vertex*>&);
  void addEdge(Edge*);
  void printFace(void);
};

Face::Face(char *triplet,std::vector<Vertex*> &vp)
{
  v[0]=v[1]=v[2]=NULL;
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
    v[0]=v[1]=v[2]=NULL;
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
    v[0]=v[1]=v[2]=NULL;
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
    v[0]=v[1]=v[2]=NULL;
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
    v[0]=v[1]=v[2]=NULL;
    printf("Error in reading vertex index\n");
    return;
  }
}

void Face::printFace(void){
  cout.precision(12);
  cout << "Face <ind>" << index << endl
        << "[v0 "
        << v[0]->index << " "
        << v[0]->p[0] << " "
        << v[0]->p[1] << " "
        << v[0]->p[2] << "]\n"
        << "[v1 "
        << v[1]->index << " "
        << v[1]->p[0] << " "
        << v[1]->p[1] << " "
        << v[1]->p[2] << "]\n"
        << "[v2 "
        << v[2]->index << " "
        << v[2]->p[0] << " "
        << v[2]->p[1] << " "
        << v[2]->p[2] << "]";
}

class Edge
{
public:
  Vertex *v1,*v2; 	// pointers to vertices on edge
  Face *f1,*f2;	        // pointers to adjacent faces (i.e. faces that contain edge)
  vec_f fvec;	        // pointers to additional adjacent faces
  bool bisected;
  Edge(Face*,Vertex*,Vertex*);
  void update(Face*);
  double getLengthSq(void);
};

Edge::Edge(Face *f,Vertex *va,Vertex *vb)
{
  v1=va;
  v2=vb;
  f1=f;
  f2=NULL;
}

double Edge::getLengthSq(void)
{
  double x = v1->p[0]-v2->p[0];
  double y = v1->p[1]-v2->p[1];
  double z = v1->p[2]-v2->p[2];
  return x*x+y*y+z*z;
}

class Distance
{
public:
  double d; // squared distance between vertices with indices vA and vB
  Vertex *vA,*vB;
  Distance(double,Vertex*,Vertex*);
};

Distance::Distance(double dist,Vertex *p,Vertex *q)
{
  d = dist;
  vA = p;
  vB = q;
}
