
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
/*
struct eqs
{
  bool operator()(const std::string s1, const std::string s2) const
  {
    return s1==s2;
  }
};*/

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

/*
struct s_hash
{
    u4 operator()(std::string i) const { return (u4) i; }
};*/

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

typedef std::map<std::string,Edge*,lts,std::allocator<Edge*> > hashtable_t;
typedef std::map<std::string,Edge*,lts,std::allocator<Edge*> >::iterator ht_iterator;
//typedef __gnu_cxx::hash_map<Face*,std::vector<Box*>*,f_hash,eqf,std::allocator<std::vector<Box*>* > > hashtable_f;
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
typedef __gnu_cxx::hash_map<Edge*,int,e_hash,eqe,std::allocator<int> > hashtable_e;
typedef std::map<Vertex*,double,ltv,std::allocator<double> > table_d;
typedef std::map<Vertex*,double,ltv,std::allocator<double> >::iterator td_iterator;
typedef std::map<Face*,double*,ltf,std::allocator<double*> > table_fd;
typedef std::map<Face*,double*,ltf,std::allocator<double*> >::iterator fd_iterator;
typedef std::multimap<double,Vertex*,gtd,std::allocator<Vertex*> > table_v;
typedef std::multimap<double,Vertex*,gtd,std::allocator<Vertex*> >::iterator tv_iterator;
typedef std::pair<double,Vertex*> vd_pair;
typedef std::pair<double,double> dd_pair;
typedef std::set<Vertex*,ltv> v_set;
typedef std::set<Vertex*,ltv>::iterator vs_iterator;
typedef std::set<Face*,ltf> f_set;
typedef std::set<Face*,ltf>::iterator fs_iterator;
typedef std::list<Vertex*> v_list;
typedef __gnu_cxx::hash_set<Face*,f_hash,eqf> hashset_f;
typedef __gnu_cxx::hash_set<Face*,f_hash,eqf>::iterator hf_iterator;
typedef __gnu_cxx::hash_set<Vertex*,v_hash,eqv> hashset_v;
typedef __gnu_cxx::hash_set<Vertex*,v_hash,eqv>::iterator hv_iterator;
typedef __gnu_cxx::hash_set<Edge*,e_hash,eqe> hashset_e;
typedef __gnu_cxx::hash_set<Edge*,e_hash,eqe>::iterator he_iterator;

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

/*
// structure for neighborhood element
class neighbor
{
public:
	Vertex *v;	// vertex* to vertex in neighborhood
	double l;	// minimum length between v and origin vertex
	neighbor(Vertex*,double);
	bool ifFrozen(void);
};

bool neighbor::ifFrozen(void){
	if (l<NEIGHBORHOOD_RADIUS){return true;}
	else {return false;}
}

neighbor::neighbor(Vertex *a, double b){
	v=a;
	l=b;
}
*/

class Point{
public:
	double x,y,z; // current vertex
	double a,b,c,L; // closest point and squared distance
	void add(double,double,double);
	Point(double,double,double);
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

// ######################################
// ######################################

class Vertex {
public:
	int index;
	double pN[3];		// current position coordinates (x,y,z)
	double pC[3];		// closest mesh position coordinates (x,y,z)
	Face *cl;			// pointer to face on which closest mesh position lies
	Object *o;			// pointer to parent object
	vector<Edge*> e;	// pointers to adjacent edges
	vector<Face*> f;	// pointers to adjacent faces
	vector<Face*> nf;	// pointers to neighborhood faces
	Vertex(char* triplet,Object*);
	void getAdjacentVertices(vector<Vertex*>&);
	void getAdjacentFaces(hashset_f&);
	double getSqSepDist(void); 
	void getNormal(double*);
	void getForceEnergy(double[3]);
    double getSeparationForceEnergy(double[3],bool);
    double getEdgeStretchForceEnergy(double[3],bool);
    double getEdgeAngleFaceIntersectionForceEnergy(double[3],bool);
//	void computeNewCoords(Container*,double[3]);
	void computeNewCoords(Container*,double[3],double);
	void assignHolding(double[3]);
	void printNeighborhood(void);
	void printVertex(std::string);
	void printVertexCP(void);
	void getBoundingBox(double[6]);
	//////////////
	bool vertexIsNice(void);
	int getVertexNiceness(void);
	void setVertexNiceness(int);
	void computeEdgeFlip(void);
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
	pC[0]=pN[0];
	pC[1]=pN[1];
	pC[2]=pN[2];
	o=q;
	cl=NULL;
}

// ######################################
// ######################################

class Face {
public:
	int index;		// face index
	Vertex *v[3];	// pointers to vertices
	Edge   *e[3];	// pointers to edges
	std::vector<Box*> b;	// pointers to boxes
	Face(char*,std::vector<Vertex*>&); 
	void addEdge(Edge*);
	void recordBoxes(vector<Box*>&); 
	void getNormal(double[3]);	// compute outward normal vector
	void getVertexCoordinates(double *[3]);
	bool computeIntersectionForce(Container*);
	void calculateIntersectionForce(Container*);
	bool getFaceIntersection(Container*,bool,std::vector<Face*>&);
	bool getFaceIntersectionCheck(Container*,Space&,hashtable_f&);
	double getAngle(Vertex *v);
	void printFace(std::string);
	void updateBoxes(hashtable_f&);
	void getBoundingBox(double br[6]);
	/////////
	std::vector<Face*> * getIntersectingFaces(void);
	void addFaceToTable_iv(void);
	void addFaceToTable_intf(void);
	bool faceInTable_iv(void);
	bool faceInTable_intf(void);
	void addForceToFace(double[3]);
	void addFaceToVector(Face*);
	void clearFaceFromTable_iv(void);
	void clearFaceFromTable_intf(void);
	void getForceFromTable(double[3]);
	bool noMoreIntersectingFaces(void);
	void removeFaceFromVector(Face*);
};

// ######################################
// ######################################

class Edge {
public:
	Object *o;		// pointer to parent object
	Vertex *v1,*v2; // pointers to vertices on edge
	Vertex *o1,*o2; // pointers to vertices on adjacent faces not on edge
	Face *f1,*f2;	// pointers to adjacent faces (i.e. faces that contain edge)
	double l;		// original edge length
	void update(Face*,Vertex*);
	Edge(Face*,Vertex*,Vertex*,Vertex*,Object*);
	double getSqLength(void);
	int computeFlip(void);
	double getAngle(void);
	double getForceEnergy(int,double[3],bool); // force[0][] = fx,fy,fz for o1 
															// force[1][] = fx,fy,fz for o2
	double getStretchForceEnergy(Vertex*,double[3],bool);
	void printEdge(std::string);
};

void Edge::printEdge(std::string s){
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

Edge::Edge(Face *f,Vertex *va,Vertex *vb,Vertex *vc,Object *op){
	v1=va;
	v2=vb;
	f1=f;
	f2=NULL;
	o=op;
	o1=vc;
	o2=NULL;
	// compute original edge length
	l=sqrt((va->pN[0]-vb->pN[0])*(va->pN[0]-vb->pN[0])+
			(va->pN[1]-vb->pN[1])*(va->pN[1]-vb->pN[1])+
			(va->pN[2]-vb->pN[2])*(va->pN[2]-vb->pN[2]));
}

double Edge::getSqLength(void){
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
	~Object(void);
	void addOriginal(int,double*);
	void createEdges(void);
	int setNumDigits(void);
	void checkEdge(Face*,Vertex*,Vertex*,Vertex*,hashtable_t&,int);
	Edge* findEdge(Vertex*,Vertex*,hashtable_t&,int);
	void createEdge(Face*,Vertex*,Vertex*,Vertex*,hashtable_t&,int);
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
//	bool processEdge(Edge*,hashtable_v_double&,std::vector<Edge*>&);	
	bool processEdge(Edge*,hashtable_v_double&,std::vector<Edge*>&,Vertex*);	
};

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

void Vertex::setVertexNiceness(int val){
	o->setVertexNiceness(this,val);
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
//	frozen=false;
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

void Face::addFaceToVector(Face* f){
	(*v[0]->o->intf[this]).push_back(f);
}

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

void Face::clearFaceFromTable_intf(void){
	// if this face is in intf table 
	if(faceInTable_intf()){
		// delete vector<face*>*
		delete v[0]->o->intf[this];
		// remove element from table
		v[0]->o->intf.erase(this);
	}
}

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
}

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
/*void Face::printFace(std::string s){
	cout << v[0]->pN[0] << " "
	<< v[0]->pN[1] << " "
	<< v[0]->pN[2] << " 1 0 0 1\n"
	<< v[1]->pN[0] << " "
	<< v[1]->pN[1] << " "
	<< v[1]->pN[2] << " 1 0 0 1\n"
	<< v[2]->pN[0] << " "
	<< v[2]->pN[1] << " "
	<< v[2]->pN[2] << " 1 0 0 1\n";
}*/

void Face::getVertexCoordinates(double *cpvc[3]){
	cpvc[0]=v[0]->pN;
	cpvc[1]=v[1]->pN;
	cpvc[2]=v[2]->pN;
}

Face::Face(char *triplet,std::vector<Vertex*> &vp){

	e[0]=e[1]=e[2]=NULL;
//	iv[0]=iv[1]=iv[2]=0;

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
}


// ######################################
// ######################################

class Container
{
public:
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
						///// because each intersecting pair of faces will be discovered twice
	double md[2];		///// mean displacement of vertices in object during each iteration
						// md[0] = mean displacement of N-1 vertices
						// md[1] = mean displacement of N vertices
	double d_min,d_max;	// min and max actual displacement
//	double se_min,se_max,se_mean;	// min and max separation error
//	double stor_min,stor_max;	// storage for printing
	double min_edge_angle;
	double gain;
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
	void buildMeshAfter(int);
	void writeDistances(int);
	void writeDistances50(int);
	void writeDistancesNOCP(int);
	void fileOutit(void);
	void updateFile(int,bool,double);
	void checkAdjacentFaces(void);
	void checkAngle(double);
	void cleanup(void);
	~Container(void);
	Container(void);
	void updateStats(double);
	void removeOldIntersections(Vertex*,hashset_v&);
	void updateNewIntersections(Vertex*,hashset_v&);
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
	void findClosestAxis(Space&,Vertex*,double[2][3]);
	int findExtraPoint(Space&,Vertex*,double[3],int);
	void findCrossed1(Space&,Vertex*,double[2][3],std::vector<Object*>&);
	void findCrossed2(Space&,double[2][3],std::vector<Object*>&);
	void getExtraRay(Vertex*,double[2][3],int);
	void collectNiceFaces(Space&,double[2][3],vector<Face*>&);
	void getBoxIndexRange(Space&,double[2][3],int[6]);
	void findIntersectedFaces(double[2][3],vector<Face*>&,
								vector<Face*>&,vector<int>&);
	void findOddMeshes(vector<Face*>&,vector<int>&,int&,std::vector<Object*>&);
	///////////
	void getSeparationDistances(Space&,Monitor&);
	bool findClosest(Space&,Vertex*,Monitor&,bool);
//	bool computeClosest(Face*,Vertex*,double&,double[3],bool);
	bool computeClosest(Face*,Vertex*,double&,double[3]);
	bool getPlaneIntersection(Face*,Vertex*,double*,double,double,Point&);
	void getEdgeIntersection(Vertex*,double*[3],Point&);
	void getBoxes(vector<Box*>&,Vertex*,int,Space&);
//	void getBoxes(vector<Box*>&,Vertex*,Space&);
	void getCandidateFaces(vector<Box*>&,Vertex*,hashset_f&);
	bool faceInNeighborhood(Face*,Vertex*);
	///////////
	bool assignNewVertexCoords(Space&,Vertex*,double[3],Monitor&); 
	void getAffectedVerticesAndEdgesBefore(Space&,Vertex*,double[3],Monitor&);
	void getAffectedVerticesAndEdgesAfter(Space&,Vertex*,double[3],Monitor&);
	bool boxesOverlap(double[6],double[6]);
	void computeNewVertexCoords(double); 
	void computeFaceIntersectionForce(void); 
	void setFaceIntersection(Face*,Space&); 
	bool checkForIntersections(Vertex*,Space&,bool,hashtable_f&);
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
	void getNiceCheckSet(Vertex*,Monitor&);
	double getVertexSqD(Vertex*);
	void collectAdjacentFaceNormals(table_fd&,Vertex*);
	void freeAdjacentFaceNormals(table_fd&);
	void getVertexAndEdgeEnergy(Monitor&);
	void computeGlobalEnergy(void);
	void updateMovedVertexEnergy(Vertex*,Monitor&);
	void updateEdgeAngles(Vertex*);
	void updateVertexVD(Vertex*,Monitor&);
};

Container::Container(void){
	si=ti=0; // twice number of intersecting faces in model, self and total, respectively
	min_edge_angle=500; // minimum edge angle in model
	s_nonnice=nonnice=0; // total number of nonnice vertices in model, self and total, respectively
	energy=0; // NOT USED
	force=0; // NOT USED
	gain = TIME_STEP/DAMPING;
	pairs[0][0] = 0;
	pairs[0][1] = 1;
	pairs[1][0] = 1;
	pairs[1][1] = 2;
	pairs[2][0] = 2;
	pairs[2][1] = 0;
	object_count=vertex_count=face_count=edge_count=0;
	scanDir();
	scanFiles();
	createEdges();
	findVertexAdjacencies();
	checkEdgeAngles();
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
	cout << "Box::printBox world 0,2,4 ["
	<< a << " "
	<< b << " "
	<< c << "]" << endl;
	cout.flush();
	cout << "Box ";
	cout.flush();
	cout << x << " ";
	cout.flush();
	cout << y << " ";
	cout.flush();
	cout << z << "]" << endl;
	cout.flush();
	cout << " limits ["
	<< xmin(a) << " "
	<< xmax(a) << " "
	<< ymin(b) << " "
	<< ymax(b) << " "
	<< zmin(c) << " "
	<< zmax(c) << "] ";
	cout.flush();
}

double Box::xmin(double i){return x*SPACE_LENGTH+i;}
double Box::xmax(double i){return (x+1)*SPACE_LENGTH+i;}
double Box::ymin(double i){return y*SPACE_LENGTH+i;}
double Box::ymax(double i){return (y+1)*SPACE_LENGTH+i;}
double Box::zmin(double i){return z*SPACE_LENGTH+i;}
double Box::zmax(double i){return (z+1)*SPACE_LENGTH+i;}

Box::Box(int a,int b,int c){
	x=a;
	y=b;
	z=c;
}

// ######################################
// ######################################

class Space {
public:
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
	void getBoxesFor3DIndices(int[6],std::vector<Box*>&,bool);
	int screenIndex(int,char*);
	Space(Container&);
	~Space(void);
};

Space::Space(Container& c){
	cout << "initializing space data structure..............";
	cout.flush();	
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
	// collect boxes to check
	getBoxesFor3DIndices(br,bp,false);
}

void Space::getBoxesFor3DIndices(int br[6],std::vector<Box*>& bp,bool shell){
	// assume br = [xmin,xmax,ymin,ymax,zmin,zmax]
	std::vector<Box*>::iterator i;
	for (int z = br[4];z<br[5]+1;z++){
		for (int y = br[2];y<br[3]+1;y++){
			for (int x = br[0];x<br[1]+1;x++){
				if (shell) {
					if (z==br[4]||z==br[5]||x==br[0]||x==br[1]||y==br[2]||y==br[3]){
						bp.push_back(b[indices2Index(x,y,z)]);
					}
				} else {
					bp.push_back(b[indices2Index(x,y,z)]);
				}
			}
		}
	}
}

int Space::indices2Index(int x,int y,int z){
	return z*num_space[0]*num_space[1]+y*num_space[0]+x;
}

void Space::index2Range(int i,double r[2],char *c){
		r[0] = i*SPACE_LENGTH;
		r[1] = (i+1)*SPACE_LENGTH;
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
	int a=0;
	if (!strcmp(c,"x")){
		a = (int) floor((ss-world[0])/SPACE_LENGTH);
		a = screenIndex(a,"x");
	} else if (!strcmp(c,"y")){
		a = (int) floor( (ss-world[2])/SPACE_LENGTH );
		a = screenIndex(a,"y");
	} else if (!strcmp(c,"z")){
		a = (int) floor( (ss-world[4])/SPACE_LENGTH );
		a = screenIndex(a,"z");
	} else {cout << "Received unexpected string.\n"; exit(0);}
	return a;
}

// ######################################
// ######################################

/*class Manip
{
public:
	///////////
	Manip(void);
};

Manip::Manip(void){
}*/

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
	int touchs;
	table_v topN;					// double f -> vertex *v
	table_d old;					// vertex *v -> double se
	hashtable_v_double v_energy;	// vertex* -> double potential energy
	hashtable_e_double e_energy;	// edge* -> double potential energy
	hashtable_v touch_map;			// vertex selection hash map
	hashtable_f nb;					// build face* hashtable
	hashset_v set_nice_check;		// store vertices for which niceness may have changed
	hashset_v set_nice_changed;		// store vertices for which niceness did change
	hashset_v av;					// store all vertices within NUM_ADJACENT_BOXES of adjacent faces
	hashset_e ae;					// store adjacent face edges
	hashtable_e_double e_angle;		// edge* -> double edge angle in radians
	double sum,avg_new,avg_old;		// sliding sum and average of global energy
	int num;						// size of sliding window
	double *window,*begin,*end;		// double pointers to elements of sliding
									// collection of global energies
//	std::vector<double> window;		// sliding collection of global energies
//	std::vector<double>::iterator next; // iterator pointing to next element in vector to replace
	v_set refrac_s;					// set of last N vertices moved
	v_list refrac_l;					// set of last N vertices moved
//	vs_iterator next;				// iterator pointing to oldest element in refractory window
	void initAvg(void);
	void freeAvg(void);
	void clearAvg(void);
	void updateAvg(double);
	void initTable(void);
	void clearTable(void);
	void saveOld(void);
	void updateTopN(Vertex*,double,double,bool);
	bool entryInTopN(Vertex*,double,tv_iterator&);
	void updateOld(Vertex*,double);
	void loadTopN(Container*);
	void validateOld(void);
	void validateTopN(char*);
	void validateVertex(int,Container&);
	void validateMultimap(void);
	void printVertexSelect(Container&,int);
	void printVertexSelect(Container&,char*);
	void updateTouchMap(Vertex*);
	void updateEnergyMap(Vertex*,double);
	void updateSets(Vertex*,double,bool);
	void prep(Container*);
	double bpf_mean, bpf_min, bpf_max;	// #boxes per face
	double fpb_mean, fpb_min, fpb_max;	// #faces per box
	void getBoxesPerFace(Container&);
	void getFacesPerBox(Space&);
	void printVerticesWithCP(void);
//	std::vector<Vertex*> added;
	void initRefrac(void);
	bool Refracted(Vertex*);
	void updateRefractoryWindow(Vertex*);
};

//############################################################################
//############################################################################








