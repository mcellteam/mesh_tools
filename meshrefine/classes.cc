// ######################################
// ######################################
class Vertex;
class Face;
class Edge;
class Object;

typedef  unsigned long int  u4;   /* unsigned 4-byte type */
typedef  unsigned     char  u1;   /* unsigned 1-byte type */

struct eqf
{
  bool operator()(const Face* s1, const Face* s2) const
  {
    return s1==s2;
  }
};

struct equ4
{
  bool operator()(const u4 s1, const u4 s2) const
  {
    return s1==s2;
  }
};

struct pair_hash
{
    u4 operator()(u4 i) const { return i; }
};

struct f_hash
{
    u4 operator()(Face* i) const { return (u4) i; }
};

struct lts
{
  bool operator()(const std::string s1, const std::string s2) const
  {
    return s1 < s2;
  }
};

typedef __gnu_cxx::hash_map<int,Vertex*> hashtable_v;
typedef __gnu_cxx::hash_set<Face*,f_hash,eqf> hashset_f;
typedef __gnu_cxx::hash_set<Face*,f_hash,eqf>::iterator hf_iterator;
typedef std::map<std::string,Edge*,lts,std::allocator<Edge*> > hashtable_t;
typedef std::map<std::string,Edge*,lts,std::allocator<Edge*> >::iterator ht_iterator;

// ######################################
// ######################################

class Vertex
{
public:
	double p[3];	// vertex position
	int index;
	Vertex(char*);
	Vertex(int,double,double,double);
	void getAdjacentVertices(std::vector<Vertex*>&);
	void getAdjacentFaces(hashset_f&);
	void getAdjacentEdges(std::vector<Edge*>&);
	void printVertex(std::string);
	void printVertexCP(void);
	bool match(int,std::string);
};

void Vertex::printVertex(std::string s){
	cout.precision(12);
	cout << "Vertex <obj>" << s << "<ind>" << index << " "
	<< "["
	<< p[0] << " "
	<< p[1] << " "
	<< p[2] << "] ";
}

void Vertex::printVertexCP(void){
	cout.precision(12);
	cout
	<< p[0] << " "
	<< p[1] << " "
	<< p[2] << " 1 0 0 1";
}

Vertex::Vertex(int i,double xval,double yval,double zval)
{
    index=i;
    p[0]=xval;
    p[1]=yval;
    p[2]=zval;
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
    Face(int,Vertex*,Vertex*,Vertex*);
	void clearEdges(void);
	bool subdivided;
	bool bisubdivided;
	void printFace(void);
	void printFaceCP(void);
	void addEdge(Edge*);
	int matchEdges(Vertex*[3]);
};

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

// THE FOLLOWING CODE IS USEFUL. IT PRINTS IN DReAMM CUSTOM POINTS FORMAT.
void Face::printFaceCP(void){
	cout << v[0]->p[0] << " "
	<< v[0]->p[1] << " "
	<< v[0]->p[2] << " 1 0 0 1\n"
	<< v[1]->p[0] << " "
	<< v[1]->p[1] << " "
	<< v[1]->p[2] << " 1 0 0 1\n"
	<< v[2]->p[0] << " "
	<< v[2]->p[1] << " "
	<< v[2]->p[2] << " 1 0 0 1\n";
}

void Face::clearEdges(void){
	e[0]=e[1]=e[2]=NULL;
	subdivided=bisubdivided=false;
}

Face::Face(int i,Vertex *va,Vertex *vb,Vertex *vc)
{
    subdivided = false;
    bisubdivided = false;
    e[0]=e[1]=e[2]=NULL;
    index=i;
    v[0]=va;
    v[1]=vb;
    v[2]=vc;
}

Face::Face(char *triplet,std::vector<Vertex*> &vp)
{
	v[0]=v[1]=v[2]=NULL;
	e[0]=e[1]=e[2]=NULL;
    subdivided = false;
    bisubdivided = false;

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

class Edge
{
public:
	Face *f1,*f2;	// pointers to adjacent faces (i.e. faces that contain edge)
	Vertex* v; // pointer to new vertex if edge is bisected
	bool bisected;
	Edge(Face*);
	void update(Face*);
	double getSqLength(void);
	void getVertices(Vertex*&,Vertex*&,Vertex*&,Vertex*&);
	bool threshold(double);
};

double Edge::getSqLength(void){
	Vertex *v1=NULL,*v2=NULL,*o1=NULL,*o2=NULL;
	getVertices(v1,v2,o1,o2);
	return (v1->p[0]-v2->p[0])*(v1->p[0]-v2->p[0])
			+(v1->p[1]-v2->p[1])*(v1->p[1]-v2->p[1])
			+(v1->p[2]-v2->p[2])*(v1->p[2]-v2->p[2]);
}

Edge::Edge(Face *f){
	f1=f;
	f2=NULL;
	v=NULL;
	bisected=false;
}

void Edge::getVertices(Vertex *&v1,Vertex *&v2,Vertex *&o1,Vertex *&o2){
	// find pair of vertices va and vb in common between f1 and f2
	Vertex *va=NULL,*vb=NULL;
	if(f1->v[0]==f2->v[0] || f1->v[0]==f2->v[1] || f1->v[0]==f2->v[2]){va=f1->v[0];}
	if(f1->v[1]==f2->v[0] || f1->v[1]==f2->v[1] || f1->v[1]==f2->v[2]){
		if(va==NULL){va=f1->v[1];}
		else {vb=f1->v[1];}
	} 
	if(f1->v[2]==f2->v[0] || f1->v[2]==f2->v[1] || f1->v[2]==f2->v[2]){
		vb=f1->v[2];
	}
	if(va==NULL || vb==NULL){
		cout << "\n\nEdge::getVertices: "
		<< "common face vertices were not identified.\n\n";
		if(f1!=NULL){f1->printFace();cout << endl;}
		else{cout << "f1 is NULL.\n";}
		if(f2!=NULL){f2->printFace();cout << endl;}
		else{cout << "f2 is NULL.\n";}
		if(va!=NULL){cout << "va index = " << va->index << endl;}
		else{cout << "va is NULL.\n";}
		if(vb!=NULL){cout << "vb index = " << vb->index << endl;}
		else{cout << "vb is NULL.\n";}
		exit(0);
	}
	// identify v1 and v2 using f1
	if( (f1->v[0]==va && f1->v[1]==vb) || (f1->v[0]==vb && f1->v[1]==va)){
		v1=f1->v[0];
		v2=f1->v[1];
		o1=f1->v[2];
	} else if( (f1->v[1]==va && f1->v[2]==vb) || (f1->v[1]==vb && f1->v[2]==va)){
		v1=f1->v[1];
		v2=f1->v[2];
		o1=f1->v[0];
	} else if( (f1->v[2]==va && f1->v[0]==vb) || (f1->v[2]==vb && f1->v[0]==va)){
		v1=f1->v[2];
		v2=f1->v[0];
		o1=f1->v[1];
	} else {
		cout << "\n\nEdge::getVertices: "
		<< "v1 and v2 were not successfully found.\n\n";
		exit(0);
	}
	// identify o2
	if( (f2->v[0]==va && f2->v[1]==vb) || (f2->v[0]==vb && f2->v[1]==va)){
		o2=f2->v[2];
	} else if( (f2->v[1]==va && f2->v[2]==vb) || (f2->v[1]==vb && f2->v[2]==va)){
		o2=f2->v[0];
	} else if( (f2->v[2]==va && f2->v[0]==vb) || (f2->v[2]==vb && f2->v[0]==va)){
		o2=f2->v[1];
	} else {
		cout << "\n\nEdge::getVertices: "
		<< "o2 was not successfully found.\n\n";
		exit(0);
	}
}

class Object
{
public:
	std::vector<Vertex*> v;		// container of pointers to all vertices in object
	std::vector<Face*> f;		// container of pointers to all faces in object
	std::vector<Edge*> e;		// container of pointers to all edges in object
        std::vector<std::string> storage; // store unamtched lines from frozen vert file
        std::vector<int> frozen;        // indices of frozen vertices
        std::string name;
	int max_faces;
	Object(std::string const &);
	void scanFile(const char*);
	void clearEdges(void);
	void createEdges(void);
	void createEdge(Face*,Vertex*,Vertex*,hashtable_t&,int);
	void checkEdge(Face*,Vertex*,Vertex*,hashtable_t&,int);
	int setNumDigits(void);
	bool thresholdEdges(double);
	int createNewVertices(void);
	int createNewSubdividedFaces(void);
        void processFrozenFile(std::string filename);
        void printFrozen(std::string filename);
        std::string strip_tail(std::string const &,uint);
        std::string getName(std::string const &);
};

std::string Object::strip_tail(std::string const &str, uint len)
{
  if (str.size() < len)
    return str;

  return str.substr(0, str.size() - len);
}

std::string Object::getName(std::string const &str)
{
  // extract name (== everything before last '_')
  uint start = 0;
  if (str.find_last_of('/') != std::string::npos)
  {
    start = str.find_last_of('/')+1;
  }
  uint end = str.size()-MESH_SUFFIX.size();
  uint len = end-start;
  return str.substr(start,len);
}


Object::Object(std::string const &file)
{
  name = getName(file);
  fprintf(stderr,"\n\nBuilding vertices and faces.....");
  fflush(stderr);
  scanFile(file.c_str());
  fprintf(stderr,"complete.\n");
  fflush(stderr);
  max_faces = f.size();
}

