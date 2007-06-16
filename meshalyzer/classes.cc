
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
class Controls;
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

struct eqi
{
  bool operator()(const int s1, const int s2) const
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

struct i_hash
{
    u4 operator()(int i) const { return (u4) i; }
};

/*
struct s_hash
{
    u4 operator()(std::string i) const { return (u4) i; }
};*/

struct lti
{
  bool operator()(const int s1, const int s2) const
  {
    return s1 < s2;
  }
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

struct gtf
{
  bool operator()(const Face* s1, const Face* s2) const
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
typedef std::multimap<Face*,Face*,ltf,std::allocator<Face*> > multimap_ff;
typedef std::multimap<Face*,Face*,ltf,std::allocator<Face*> >::iterator mmff_iterator;
typedef std::multimap<Face*,Box*,ltf,std::allocator<Box*> > hashtable_f;
typedef std::multimap<Face*,Box*,ltf,std::allocator<Box*> >::iterator tf_iterator;
typedef __gnu_cxx::hash_map<int,int,i_hash,eqi,std::allocator<int> > hashtable_ii;
typedef __gnu_cxx::hash_map<Vertex*,int,v_hash,eqv,std::allocator<int> > hashtable_v;
typedef __gnu_cxx::hash_map<Vertex*,int,v_hash,eqv,std::allocator<int> >::iterator vhm_iterator;
typedef __gnu_cxx::hash_map<Face*,double*,f_hash,eqf,std::allocator<double*> > hashtable_f_double;
typedef __gnu_cxx::hash_map<Face*,int,f_hash,eqf,std::allocator<int> > hashtable_fi;
typedef __gnu_cxx::hash_map<Face*,int,f_hash,eqf,std::allocator<int> >::iterator fi_iterator;
typedef __gnu_cxx::hash_map<Vertex*,double,v_hash,eqv,std::allocator<double> > hashtable_v_double;
typedef __gnu_cxx::hash_map<Vertex*,double,v_hash,eqv,std::allocator<double> >::iterator vdhm_iterator;
typedef __gnu_cxx::hash_map<Edge*,double,e_hash,eqe,std::allocator<double> > hashtable_e_double;
typedef __gnu_cxx::hash_map<Edge*,double,e_hash,eqe,std::allocator<double> >::iterator edhm_iterator;
typedef __gnu_cxx::hash_map<Face*,std::vector<Face*>*,f_hash,eqf,std::allocator<std::vector<Face*>* > > hashtable_f_face;
typedef __gnu_cxx::hash_map<Face*,std::vector<Face*>*,f_hash,eqf,std::allocator<std::vector<Face*>* > >::iterator ff_iterator;
typedef __gnu_cxx::hash_map<Edge*,int,e_hash,eqe,std::allocator<int> > hashtable_e;
typedef __gnu_cxx::hash_map<int,Vertex*,i_hash,eqi,std::allocator<Vertex*> > hashtable_iv;
typedef std::map<int,bool,lti,std::allocator<bool> > table_ib;
typedef std::map<int,bool,lti,std::allocator<bool> >::iterator ib_iterator;
typedef std::map<int,Vertex*,lti,std::allocator<Vertex*> > table_iv;
typedef std::map<Vertex*,double,ltv,std::allocator<double> > table_d;
typedef std::map<Vertex*,double,ltv,std::allocator<double> >::iterator td_iterator;
typedef std::map<Face*,double*,ltf,std::allocator<double*> > table_fd;
typedef std::map<Face*,double*,ltf,std::allocator<double*> >::iterator fd_iterator;
typedef std::multimap<double,Vertex*,gtd,std::allocator<Vertex*> > table_v;
typedef std::multimap<double,Vertex*,gtd,std::allocator<Vertex*> >::iterator tv_iterator;
typedef std::multimap<int,Face*,lti,std::allocator<Face*> > table_if;
typedef std::multimap<int,Face*,lti,std::allocator<Face*> >::iterator if_iterator;
typedef std::pair<double,Vertex*> vd_pair;
typedef std::pair<double,double> dd_pair;
typedef std::set<int,lti> i_set;
typedef std::set<int,lti>::iterator is_iterator;
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
class Stats{
public:
	int n;
	double total;
	double min;
	double max;
	double sum;
	double sum2;
	int count[NUM_BINS];
	double bins[NUM_BINS+1];
	std::vector<double> x;
	void clear(void);
	Stats(void);
	double mean(void);
	double median(void);
	double variance(void);
	void createHistogram(void);
	void createAspectRatioHistogram(void);
	void printHistogram(void);
	void printAspectRatioHistogram(void);
	void printStats(void);
};

void Stats::printHistogram(void){
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

void Stats::printAspectRatioHistogram(void){
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

void Stats::printStats(void){
	cout << "   min      " << min << endl;
	cout << "   max      " << max << endl;
	cout << "   median   " << median() << endl;
	cout << "   mean     " << mean() << endl;
	cout << "   variance " << variance() << endl << endl;
}

double Stats::median(void){
	sort(x.begin(),x.end());
	int i= x.size()/2;
	return x[i];
}

void Stats::createHistogram(void){
	// create bins
	double range=max-min;
	double step = range/NUM_BINS;
	bins[0]=min;
	for(int i=1;i<NUM_BINS;i++){
		bins[i]=min+step*i;
	}
	bins[NUM_BINS]=max;

	// bin data
	unsigned int start,mid,end;

	// for each element in x
	for(std::vector<double>::iterator i=x.begin();i!=x.end();i++){
		start=0;
		end=NUM_BINS-1;
		while (end-start>1){
			mid=(start+end)/2;
			if (*i >= bins[start] && *i <= bins[mid]) { end=mid;}
			else { start=mid; }
		}
		if (*i >= bins[start] && *i <= bins[end]) { count[start]++;}
		else { count[end]++; }
	}
}
void Stats::createAspectRatioHistogram(void){
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

	// for each element in x
	for(std::vector<double>::iterator i=x.begin();i!=x.end();i++){
		start=0;
		end=NUM_BINS-1;
		while (end-start>1){
			mid=(start+end)/2;
			if (*i >= bins[start] && *i <= bins[mid]) { end=mid;}
			else { start=mid; }
		}
		if (*i >= bins[start] && *i <= bins[end]) { count[start]++;}
		else { count[end]++; }
	}
}

double Stats::mean(void){
	return sum/n;
}

double Stats::variance(void){
 	return sum2/n;
}

Stats::Stats(void){
	n=0;
	total=sum=sum2=0;
	min=DBL_MAX;
	max=-DBL_MAX;
	for(int i=0;i<NUM_BINS;i++){
		count[i]=0;
		bins[i]=0;
	}
	bins[NUM_BINS]=0;
}

void Stats::clear(void){
	n=0;
	sum=sum2=0;
	min=DBL_MAX;
	max=-DBL_MAX;
	for(int i=0;i<NUM_BINS;i++){
		count[i]=0;
		bins[i]=0;
	}
	bins[NUM_BINS]=0;
}

class Controls{
public:
    bool folder;	// true if folder, false if file
    bool attr;		// true=eval attributes only
    bool print;     // true=print offenders
    bool report;    // true=report goodness
	bool interf;	// true=detect intersections between objects (aka batch-mode, false=serial mode)
	bool sepdist;	// true=compute separation distances
    double t[3];	// thresholds aspect ratio, edge angle, edge length
	// ***************** cumulative stats ******************
	double bb[6];	// data set bounding box [xmin,xmaxm,ymin,ymax,zmin,zmax]
	double a;		// cumulative area
	int num_obj;	// number of objects
	int num_sep;	// number of components (separate meshes in object)
	int num_bou;	// number of separate boundaries in object, same as # holes
	int num_bor;	// number of border edges
	int num_nonman_e;	// number of nonmanifold edges
	int num_nonman_v;	// number of nonmanifold vertices
	int num_flip;	// number of flipped edges
	int num_intf;	// number of intersecting faces
	int num_vert;	// number of vertices
	int num_faces;	// number of faces
	int num_edges;	// number of edges
	double vol;		// cumulative mesh volume
	int num_orphan;	// number of orphan vertices
	int num_missing;// number of missing vertices
	int num_degen;	// number of degenerate faces
	int num_score;	// cumulative goodness score
	int num_indistin;// cumulative number of indistinguishable vertices
	// ***************** individual stats ******************
	bool contig_v;	// true if vertex indexing is contiguous
	int vec_cv[5];	// bad vertex index +/- 2
	bool contig_f;	// true if vface indexing is contiguous
	int vec_cf[5];	// bad face index +/- 2
//	multimap_ff ff;	// multimap of intersecting faces in entire container
	// multimap Face*->Face*
	Stats area;
	Stats aspect_ratio;
	Stats edge_length;
	Stats edge_angle;
	std::vector<Edge*> border;		// border edges (edge with a single adjacent face)
	std::vector<Edge*> nonman_e;	// nonmanifold edges (edge with three or more adjacent faces)
	std::vector<Edge*> flipped;		// flipped edges (edges traversed more than once)
	std::vector<Vertex*> indistin_v;// indistinguishable vertices 
	std::vector<Edge*> indistin;	// indistinguishable edges (edges with vertices 
									// that are indistinguishbale using double precision)
	std::vector<Vertex*> nonman_v;	// nonmanifold vertices 
	std::vector<Vertex*> orphan;	// vertices not referenced by any face
	std::vector<Face*> missing_f;	// faces with one or more nonexistant referenced vertices
	std::vector<int> missing_v;	// nonexistant referenced vertices
	std::vector<Face*> degen;	// faces with two or more identical vertex references
	// *****************
    Controls(void);
	void clear(void);
    void printClass(void);
	void parse(int,char**,char[2048]);
	bool gParse(char*);
	void computeVolume(Object*);
	void findIntersectingFaces(Container*,Object*);
	void processEdgeLengths(Object*);
	void orphanMissingContig(Object*);
	void areaAspectRatio(Object*);
	void DegenContig(Object*);
	int isGood(Object*);
	void computeEdgeAngles(Object*);
	void printIntegrity(void);
	void printCumulative(bool,Object*);
	bool binaryOutput(void);
};

Controls::Controls(void){
    folder=attr=print=report=interf=sepdist=contig_f=contig_v=false;
    t[0]=t[1]=t[2]=0.0;
	bb[0]=bb[1]=bb[2]=bb[3]=bb[4]=bb[5]=0;
	a=vol=0.0;
	num_obj=num_sep=num_bou=num_bor=num_nonman_e=num_nonman_v=num_flip=num_vert=0;
	num_faces=num_edges=num_orphan=num_missing=num_degen=num_intf=num_score=num_indistin=0;
	vec_cv[0]=vec_cv[1]=vec_cv[2]=vec_cv[3]=vec_cv[4]=0;
	vec_cf[0]=vec_cf[1]=vec_cf[2]=vec_cf[3]=vec_cf[4]=0;
	area.clear();
	aspect_ratio.clear();
	edge_length.clear();
	edge_angle.clear();
	border.clear();
	nonman_e.clear();
	flipped.clear();
	indistin.clear();
	nonman_v.clear();
	orphan.clear();
	missing_f.clear();
	missing_v.clear();
	degen.clear();
	indistin_v.clear();
}

void Controls::clear(void){
	contig_v=contig_f=false;
	vec_cv[0]=vec_cv[1]=vec_cv[2]=vec_cv[3]=vec_cv[4]=0;
	vec_cf[0]=vec_cf[1]=vec_cf[2]=vec_cf[3]=vec_cf[4]=0;
	area.clear();
	aspect_ratio.clear();
	edge_length.clear();
	edge_angle.clear();
	border.clear();
	nonman_e.clear();
	flipped.clear();
	indistin.clear();
	nonman_v.clear();
	orphan.clear();
	missing_f.clear();
	missing_v.clear();
	degen.clear();
}

void Controls::printClass(void){
    cout << "\n\nControls class\n"
    << "folder " << folder << endl
    << "attr " << attr << endl
    << "print " << print << endl
    << "report " << report << endl
    << "art " << t[0] << endl
    << "eat " << t[1] << endl
    << "elt " << t[2] << endl << endl;
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
	void getNormal(double*);
	void printVertex(std::string);
	bool vertexIsNice(void);
	int getVertexNiceness(void);
	void setVertexNiceness(int);
};

void Vertex::printVertex(std::string s){
	cout.precision(12);
	cout << "Vertex <obj>" << s << "<ind>" << index << " "
	<< "["
	<< pN[0] << " "
	<< pN[1] << " "
	<< pN[2] << "] ";
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
	double getAspectRatio(void);
	int index;		// face index
	Vertex *v[3];	// pointers to vertices
	Edge   *e[3];	// pointers to edges
	std::vector<Box*> b;	// pointers to boxes
	Face(char*,hashtable_iv&,Object*); 
	void addEdge(Edge*);
	void recordBoxes(vector<Box*>&); 
	void getNormal(double[3]);	// compute outward normal vector
	void getVertexCoordinates(double *[3]);
	bool getFaceIntersection(Container*,bool,std::vector<Face*>&);
	bool getFaceIntersectionCheck(Container*,Space&,hashtable_f&);
	double getAngle(Vertex *v);
	void printFace(Object*);
	void addFaceToTable_intf(void);
	bool faceInTable_intf(void);
	void addFaceToVector(Face*);
	void clearFaceFromTable_intf(void);
};

// ######################################
// ######################################

class Edge {
public:
	Object *o;					// pointer to parent object
	Vertex *v1,*v2; 			// pointers to vertices on edge
	Vertex *o1,*o2; 			// pointers to vertices on adjacent faces not on edge
	Face *f1,*f2;				// pointers to adjacent faces (i.e. faces that contain edge)
	std::vector<Face*> fvec;	// pointers to additional adjacent faces
	double l;					// original edge length
	Edge(Face*,Vertex*,Vertex*,Vertex*,Object*);
	void update(Face*,Vertex*);
	double getSqLength(void);
	double getAngle(void);
	void printEdge(std::string);
	bool isConsistent(void);
};

void Edge::printEdge(std::string s){
	cout.precision(12);
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
	<< " f1 "<< f1->index << " ["
	<< f1->v[0]->index << " "
	<< f1->v[1]->index << " "
	<< f1->v[2]->index << "],";
	if (f2!=NULL){
		cout << " f2 "<< f2->index << " ["
		<< f2->v[0]->index << " "
		<< f2->v[1]->index << " "
		<< f2->v[2]->index << "]\n";
	} else {
		cout << " f2 is NULL\n";
	}
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
	hashtable_ii ii;		// face index -> missing vertex index
	std::string name;	// object name
	vector<Vertex*> v;		// container of pointers to all vertices in object
	vector<Face*> f;		// container of pointers to all faces in object
	vector<Edge*> e;		// container of pointers to all edges in object
	bool closed;			// true=closed mesh
	bool consistent;		// true=consistently-oriented face normals
	bool outward;			// true=outwardly-oriented face normals
	bool manifold;			// true=2d manifold in R3
	// *********************
	int num_sep;		// number of components (separate meshes in object)
	int num_bou;		// number of separate boundaries in object, same as # holes
	int score;
	double vol;				// object volume
	int genus;				// object genus
	void analyze(Container*,Controls&);
	void analyzeBatch(Container*,Controls&,Space&);
	void print(Controls&,bool);
	void printChars(Controls&);
	int getScore(Controls&);
	void checkIntegrity(Controls&);
	void setAll(Edge*,hashtable_fi&,int&);
	void getGroups(Edge*,hashtable_fi&,i_set&);
	int getLowest(i_set&);
	void replaceGroups(hashtable_fi&,i_set&,int);
    void printAttr(Controls&);
	void vertexDistin(Controls&);
	// *********************
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
	void findVertexAdjacencies(Controls&);
	void findNeighborhoods(void);
	void boundObject(double*);
	double getMeanEdgeLength(void);
	//////////////
	hashtable_f_double iv; // store intersection force
	hashtable_f_face intf; // store intersecting faces
	bool intersectingFacesExist(Face*);
	bool faceInTable_intf(Face*f);
//	bool faceInTable_iv(Face*f);	
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
	//////////////
	void evalCharacteristics(Container*,Controls&);
	void evalAttributes(Container*,Controls&);
	bool isClosed(Controls&);
	bool isManifold(Controls&);
	bool isConsistent(Controls&);
	bool isOutward(Container*);
	bool rayIntersectsBB(double[3][3],Face*,double[3]);
	bool rayIntersectsSide(char*,double[3][3],double[6],double[3]);
	void countBoundaries(Controls&);
	int countComponents(void);
	void computeGenus(void);
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
	v.clear();
	f.clear();
	e.clear();
}


Object::Object(std::string s) {
	name=s;
	closed=consistent=outward=manifold=false;
	score = 0;
}

bool Object::faceInTable_intf(Face *ff){
	return intf.find(ff)!=intf.end();
}

//bool Object::faceInTable_iv(Face *ff){
//	return iv.find(ff)!=iv.end();
//}

void Face::addFaceToTable_intf(void){
	// create new face vector in table
	std::vector<Face*> *nfv = new std::vector<Face*>();
	v[0]->o->intf[this]=nfv;
}

bool Face::faceInTable_intf(void){
	return v[0]->o->faceInTable_intf(this);
}

void Face::addFaceToVector(Face* f){
	(*v[0]->o->intf[this]).push_back(f);
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

// THE FOLLOWING CODE IS USEFUL. IT PRINTS IN COMPREHENSIVE FORMAT.
void Face::printFace(Object *o){
	cout.precision(12);
	cout << "Face <obj>" << o->name << "<ind>" << index << endl;
	if(v[0]!=NULL){
		cout << "[v0 "
		<< v[0]->index << " "
		<< v[0]->pN[0] << " "
		<< v[0]->pN[1] << " "
		<< v[0]->pN[2] << "]\n";
	} else {
		cout << "[v0 "
		<< o->ii[index] << "]\n";
	}

	if(v[1]!=NULL){
		cout << "[v1 "
		<< v[1]->index << " "
		<< v[1]->pN[0] << " "
		<< v[1]->pN[1] << " "
		<< v[1]->pN[2] << "]\n";
	} else {
		cout << "[v1 "
		<< o->ii[index] << "]\n";
	}

	if(v[2]!=NULL){
		cout << "[v2 "
		<< v[2]->index << " "
		<< v[2]->pN[0] << " "
		<< v[2]->pN[1] << " "
		<< v[2]->pN[2] << "]";
	} else {
		cout << "[v2 "
		<< o->ii[index] << "]";
	}
}

void Face::getVertexCoordinates(double *cpvc[3]){
	cpvc[0]=v[0]->pN;
	cpvc[1]=v[1]->pN;
	cpvc[2]=v[2]->pN;
}

Face::Face(char *triplet,hashtable_iv &vp,Object *o){

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
	int zz = (int)strtod(val,&eptr);
	if(vp.find(zz)!=vp.end()){
	    v[0] = vp[zz];
	} else {
	    v[0] = NULL;
		o->ii[index]=zz;
	}
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
	zz = (int)strtod(val,&eptr);
	if(vp.find(zz)!=vp.end()){
	    v[1] = vp[zz];
	} else {
	    v[1] = NULL;
		o->ii[index]=zz;
	}
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
	zz = (int)strtod(val,&eptr);
	if(vp.find(zz)!=vp.end()){
	    v[2] = vp[zz];
	} else {
	    v[2] = NULL;
		o->ii[index]=zz;
	}
    if (val==eptr)
    {
		v[0]=v[1]=v[2]=NULL;
        printf("Error in reading vertex index\n");
        return;
    }
}


// ######################################
// ######################################

class Container
{
public:
	std::vector<Object*> o;				//vector of object pointers
	std::vector<std::string> files;	// array of input file names
	int num_files;					// number of input files
	int object_count;				// total number of objects in model
	int vertex_count;				// total number of vertices in model
	int face_count;					// total number of faces in model
	int edge_count;					// total number of edges in model
	double world[6];				// minX, maxX, minY, maxY, minZ, maxZ of world
	int pairs[3][2];
	///////////
	~Container(void);
	Container(void);
	///////////
	void clear(void);
	int getScore(void);
	void printBatch(Controls&);
	void scanDir(Controls&,char*);
	void scanFiles();
	void scanFile(Object*,char*);
	void createEdges(void);
	void addVertexPointers(void);
	void findVertexAdjacencies(Controls&);
	void writeDistances(void);
	///////////
	void boundWorld(void);
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
	void findIntersectedFaces(double[2][3],vector<Face*>&,vector<Face*>&,vector<int>&);
	void findOddMeshes(vector<Face*>&,vector<int>&,int&,std::vector<Object*>&);
	///////////
	void getSeparationDistances(Space&);
	bool findClosest(Space&,Vertex*);
	bool computeClosest(Face*,Vertex*,double&,double[3]);
	bool getPlaneIntersection(Face*,Vertex*,double*,double,double,Point&);
	void getEdgeIntersection(Vertex*,double*[3],Point&);
	void getBoxes(vector<Box*>&,Vertex*,int,Space&);
	void getCandidateFaces(vector<Box*>&,Vertex*,hashset_f&);
	bool faceInNeighborhood(Face*,Vertex*);
	///////////
	bool checkFaceFaceIntersections(Face*,Face*);
	int checkEdgeEdgeIntersection(Face*,Face*,bool);
	int checkFaceEdgeIntersection(Face*,Face*);
	bool facesParallel(Face*,Face*);
	bool facesColinear(Face*,Face*);
	int numUniqueVertices(Face*,Face*,int[2]);
	void assignFacesToBoxes(Space&);
};

Container::Container(void){
	pairs[0][0] = 0;
	pairs[0][1] = 1;
	pairs[1][0] = 1;
	pairs[1][1] = 2;
	pairs[2][0] = 2;
	pairs[2][1] = 0;
	object_count=vertex_count=face_count=edge_count=0;
	world[0]=world[1]=world[2]=world[3]=world[4]=world[5]=0.0;
}

Container::~Container(void){
	std::vector<Object*>::iterator i;
	for(i=o.begin();i!=o.end();i++){ delete *i; }
}

void Container::clear(void){
	o.clear();
//	files.clear();
	object_count=vertex_count=face_count=edge_count=0;
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
	~Box(void);
};

Box::~Box(void){
	f.clear();
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
	int num_space[3];
	int num_boxes;
	vector<Box*> b;	// vector of boxes
	double world[6];	// minX, maxX, minY, maxY, minZ, maxZ of world
	void boundWorld(Container&);
	void initBoxes(Container&);
	void clearBoxes(void);
	void deleteBoxes(void);
	void recordFace(std::vector<Box*>&,Face*);
	void computeBoxesToCheck(Face*,std::vector<Box*>&);
	void index2Range(int,double[2],char*);
	int location2Index(double,char*);
	int indices2Index(int,int,int);
	void getBoxesFor3DLocations(double[6],std::vector<Box*>&);
	void getBoxesFor3DIndices(int[6],std::vector<Box*>&,bool);
	int screenIndex(int,char*);
	~Space(void);
};

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

//############################################################################
//############################################################################
