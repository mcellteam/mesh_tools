// ######################################
// ######################################

class Boundary;
class Box;
class Container;
class Controls;
class Edge;
class Face;
class Face_Pair;
class Manip;
class Monitor;
class Object;
class Neighbor;
class Space;
class Vertex;
class Contour;

// ######################################
// ######################################

typedef  unsigned long int  u4;   /* unsigned 4-byte type */
typedef  unsigned     char  u1;   /* unsigned 1-byte type */

// ######################################
// ######################################

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

struct lti
{
  bool operator()(const int s1, const int s2) const
  {
    return s1 < s2;
  }
};

struct gti
{
  bool operator()(const int s1, const int s2) const
  {
    return s1 > s2;
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


struct lts
{
  bool operator()(const std::string s1, const std::string s2) const
  {
    return s1 < s2;
  }
};

struct vSortComp
{
    bool operator()(const std::pair<double,Vertex*> lhs,const std::pair<double,Vertex*> rhs)
    {
        //return (lhs->index < rhs->index);
        return (lhs.second > rhs.second);
    }
};

// ######################################
// ######################################

typedef std::vector<Vertex*>			vec_v;
typedef std::vector<Vertex*>::iterator	v_iterator;
typedef std::vector<Face*>				vec_f;
typedef std::vector<Face*>::iterator	f_iterator;
typedef std::vector<Edge*>				vec_e;
typedef std::vector<Edge*>::iterator	e_iterator;
typedef std::vector<Object*>			vec_o;
typedef std::vector<Object*>::iterator	o_iterator;
typedef std::vector<Box*>				vec_b;
typedef std::vector<Box*>::iterator		b_iterator;
typedef std::vector<int>				vec_i;
typedef std::vector<int>::iterator		i_iterator;
typedef std::vector<double>				vec_d;
typedef std::vector<double>::iterator	d_iterator;

typedef std::map<std::string,Edge*,lts,std::allocator<Edge*> >					map_se;

typedef std::multimap<int,bool,lti,std::allocator<bool> >						mmap_ib;
typedef std::multimap<int,bool,lti,std::allocator<bool> >::iterator				ib_iterator;
typedef std::multimap<double,Vertex*,gtd,std::allocator<Vertex*> >				mmap_dv;
typedef std::multimap<double,Vertex*,gtd,std::allocator<Vertex*> >::iterator	dv_iterator;
typedef std::multimap<int,Face*,lti,std::allocator<Face*> >						mmap_if;
typedef std::multimap<int,Face*,lti,std::allocator<Face*> >::iterator			if_iterator;
typedef std::multimap<int,Vertex*,gti,std::allocator<Vertex*> >					mmap_iv;
typedef std::multimap<int,Vertex*,gti,std::allocator<Vertex*> >::iterator		iv_iterator;

typedef __gnu_cxx::hash_map<Vertex*,int,v_hash,eqv,std::allocator<int> >								hmap_vi;
typedef __gnu_cxx::hash_map<Face*,double*,f_hash,eqf,std::allocator<double*> >							hmap_fdp;
typedef __gnu_cxx::hash_map<Face*,double,f_hash,eqf,std::allocator<double*> >							hmap_fd;
typedef __gnu_cxx::hash_map<Face*,double,f_hash,eqf,std::allocator<double*> >::iterator					fd_iterator;
typedef __gnu_cxx::hash_map<Face*,int,f_hash,eqf,std::allocator<int> >									hmap_fi;
typedef __gnu_cxx::hash_map<Face*,int,f_hash,eqf,std::allocator<int> >::iterator						fi_iterator;
typedef __gnu_cxx::hash_map<Vertex*,double,v_hash,eqv,std::allocator<double> >							hmap_vd;
typedef __gnu_cxx::hash_map<Vertex*,double,v_hash,eqv,std::allocator<double> >::iterator				vd_iterator;
typedef __gnu_cxx::hash_map<Edge*,double,e_hash,eqe,std::allocator<double> >							hmap_ed;
typedef __gnu_cxx::hash_map<Edge*,double,e_hash,eqe,std::allocator<double> >::iterator					ed_iterator;
typedef __gnu_cxx::hash_map<Face*,vec_f*,f_hash,eqf,std::allocator<std::vector<Face*>* > >				hmap_ff;
typedef __gnu_cxx::hash_map<Face*,vec_f*,f_hash,eqf,std::allocator<std::vector<Face*>* > >::iterator	ff_iterator;
typedef __gnu_cxx::hash_map<Edge*,int,e_hash,eqe,std::allocator<int> >									hmap_e;
typedef __gnu_cxx::hash_map<int,Vertex*,i_hash,eqi,std::allocator<Vertex*> >							hmap_iv;

typedef std::set<int,lti>									set_i;
typedef std::set<Vertex*,ltv>								set_v;

typedef __gnu_cxx::hash_set<Face*,f_hash,eqf>				hset_f;
typedef __gnu_cxx::hash_set<Face*,f_hash,eqf>::iterator		hf_iterator;

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

class Contour
{
public:
	vec_v p;			// Vertex* points along contour
	void print(Controls &cs,std::string);
};

// ######################################
// ######################################

class Stats{
public:
//	int n;
	double total;
	double min;
	double max;
	double sum;
	double sum2;
	int count[NUM_BINS];
	double bins[NUM_BINS+1];
	vec_d x;
	void clear(void);
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
};

Stats::Stats(void){
//	n=0;
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
//	n=0;
	sum=sum2=0;
	min=DBL_MAX;
	max=-DBL_MAX;
	for(int i=0;i<NUM_BINS;i++){
		count[i]=0;
		bins[i]=0;
	}
	bins[NUM_BINS]=0;
	x.clear();
}

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

void Stats::printAdjacentFaceHistogram(void){
	int L[16],U[16];
	for(int i=0;i<16;i++){
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
	cout << "       min      " << min << endl;
	cout << "       max      " << max << endl;
	cout << "       median   " << median() << endl;
	cout << "       mean     " << mean() << endl;
	cout << "       variance " << variance() << endl << endl;
}

double Stats::median(void){
	sort(x.begin(),x.end());
	int i= x.size()/2;
//	cout << "\nStats::median: "
//	<< "x.size() = " << x.size()
//	<< ", i = " << i << endl;
	return x[i];
}

void Stats::createHistogram(void){

	double bin_width=(4.0*sqrt(variance()))/(NUM_BINS-1.0);

	bool nonneg=false;

	double mm = mean();

	while(nonneg==false){
		// create bins
		bins[0]=0;
		bins[1]=0;
	
		for(int i=2;i<NUM_BINS;i++){
			bins[i]=mm+(((i-1)-(NUM_BINS/2.0))*bin_width);	
		}
		bins[NUM_BINS]=max;

		nonneg=true;
		//check
		for(int i=2;i<NUM_BINS+1;i++){
			if(bins[i]<0){
				nonneg=false;
//				cout << "Stats::createHistogram: " << bins[i] << " < 0\n";
			}
		}
		if(nonneg==false){
			// increase mm
			mm = mm*1.1;
		}
	}

	// bin data
	unsigned int start,mid,end;

	// for each element in x
	for(d_iterator i=x.begin();i!=x.end();i++){
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

void Stats::createAdjacentFaceHistogram(void){
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

	// for each element in x
	for(d_iterator i=x.begin();i!=x.end();i++){
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
	for(d_iterator i=x.begin();i!=x.end();i++){
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
	return sum/x.size();
}

double Stats::variance(void){
	if(x.size()>1){
		return (sum2-(sum*sum/x.size()))/(x.size()-1);
	} else {
//		cout << "variance not computed since sample size==1.\n";
//		return -666;
		return 0;
	}
}

// ######################################
// ######################################

class Controls{
public:
	double zEM;
	std::string inpath,outpath;
	// command line info
    bool folder;			// true if folder, false otherwise
    bool attr;				// true if '-a', false otherwise
							// true=eval attributes only
							// false=allow a full analysis, i.e. evaluate characteristics
    bool print;				// true if '-p', false otherwise
							// true=print detailed information about offending mesh elements
							// false=print summary information only
	char style[32];			// output style for offending mesh elements (-p option)
							// cp=dreamm custom points format, i.e. x y z 1 0 0 1
							// everything else=detailed face,vertex,edge information
	bool interf;			// true if '-i', false otherwise
							// true=collect all objects,
							//	detect intersections between objects,
							//	in contrast to the intra-object search that is the default
							// false=do not perform inter-object face intersection check
	bool signal[5];			// true if threshold value found on command line
    double thresholds[5];	// user-defined thresholds found on command line
							// [0]=aspect ratio threshold set
							// [1]=min edge angle threshold set
							// [2]=max edge angle threshold set
							// [3]=min edge length threshold set
							// [4]=max edge length threshold set
	// cumulative statistics
	double bb[6];			// bounding box [xmin,xmaxm,ymin,ymax,zmin,zmax]
//	double a;				// area
	bool good_integrity;	// integrity of mesh
	Stats area;
	Stats aspect_ratio;
	Stats edge_length;
	Stats edge_angle;
	Stats adjacent_face;
	// *****************
    Controls(void);
//	void parse(int,char**,std::string,char[FILENAME_SIZE]);
	void parse(int,char**,std::string);
	void printCumulative(Container&);
	void analyzeBatch(Container*,Space&);
	void store(Object*);
	void printIntegrity(Container&);
    void printAttr(Container&);
    void printChars(Container&);
	void vertexAdjacentFaces(void);
	void areaAspectRatio(void);
	void processEdgeLengths(void);
	void computeEdgeAngles(void);
//	void analyzeCumulative(Container*);
	void analyzeCumulative(void);
	
};

Controls::Controls(void){
	// command line arguments
	inpath.clear();
	outpath.clear();
    folder=attr=print=interf=false;
	signal[0]=signal[1]=signal[2]=signal[3]=signal[4]=false;
    thresholds[0]=thresholds[1]=thresholds[2]=thresholds[3]=thresholds[4]=0.0;
	// cumulative statistics
	good_integrity=true;
	bb[0]=bb[1]=bb[2]=bb[3]=bb[4]=bb[5]=0;
	// cumulative Stats
	area.clear();
	aspect_ratio.clear();
	edge_length.clear();
	edge_angle.clear();
	adjacent_face.clear();
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
	void getAdjacentFaces(hset_f&);
	void getNormal(double*);
	void printVertex(std::string);
	void printVertexCP(void);
	bool vertexIsNice(void);
	int getVertexNiceness(void);
	void setVertexNiceness(int);
	bool isManifold(bool);
	bool scanAdjFaces(Edge*,Face*,bool&);
};

void Vertex::printVertexCP(void){
	cout.precision(12);
	cout << pN[0] << " "
	<< pN[1] << " "
	<< pN[2] << " 1 0 0 1\n";
}

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
	int index;		// face index
	int vi[3];		// vertex indices
	Vertex *v[3];	// pointers to vertices
	Edge   *e[3];	// pointers to edges
	vec_b b;		// pointers to boxes
	Face(char*,Object*); 
	Face(Vertex*,Vertex*,Vertex*); 
	void addEdge(Edge*);
	void recordBoxes(vector<Box*>&); 
	void getNormal(double[3]);
	void getVertexCoordinates(double *[3]);
	double getAngle(Vertex *v);
	void printFace(Object*);
	void printFaceCP(void);
	void addFaceToTable_intf(void);
	bool faceInTable_intf(void);
	ff_iterator findFaceInTable_intf(void);
	void addFaceToVector(Face*);
	void clearFaceFromTable_intf(void);
	Edge* getNewEdge(Edge*,Vertex*);
	bool match(int);
	double getAspectRatio(void);
	void getzEM(Vertex*&,Vertex*&,double);
};

bool Face::match(int i){
	return i==index;
}

Face::Face(Vertex *v1,Vertex *v2,Vertex *v3){
	v[0]=v1;
	v[1]=v2;
	v[2]=v3;
}

// ######################################
// ######################################

class Edge {
public:
	Vertex *vv1,*vv2; 			// pointers to vertices on edge
	Face *f1,*f2;				// pointers to adjacent faces (i.e. faces that contain edge)
	vec_f fvec;					// pointers to additional adjacent faces
	double l;					// original edge length
	Edge(Face*,Vertex*,Vertex*);
	void update(Face*);
	double getSqLength(void);
	double getAngle(void);
	void printEdge(std::string);
	void printEdgeCP();
	bool isConsistent(void);
	bool isManifold(void);
	bool getStartingFace(Face*&);
	Face* getNewFace(Face*);
	void getVertices(Vertex*&,Vertex*&,Vertex*&,Vertex*&);
	bool valid(void);
};

Edge::Edge(Face *f,Vertex *va,Vertex *vb){
	vv1=va;
	vv2=vb;
	f1=f;
	f2=NULL;
	// compute original edge length
	l=sqrt((va->pN[0]-vb->pN[0])*(va->pN[0]-vb->pN[0])+
			(va->pN[1]-vb->pN[1])*(va->pN[1]-vb->pN[1])+
			(va->pN[2]-vb->pN[2])*(va->pN[2]-vb->pN[2]));
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
	std::string name;		// object name
	vector<Vertex*> v;		// container of pointers to all vertices in object
	vector<Face*> f;		// container of pointers to all faces in object
	vector<Edge*> e;		// container of pointers to all edges in object
	bool closed;			// true=closed mesh
	bool consistent;		// true=consistently-oriented face normals
	bool outward;			// true=outwardly-oriented face normals
	bool manifold;			// true=2d manifold in R3
	hmap_fdp iv;			// store intersection force
	hmap_ff intf;			// store intersecting faces
	hmap_vi nice;			// vertex* is key to int nice value
	hmap_e flip;			// collection of flip values for edges
	int num_sep;			// number of components (separate meshes in object)
	int num_bou;			// number of separate boundaries in object, same as # holes
	double vol;				// object volume
	int genus;				// object genus
	mmap_iv vp;
	mmap_ib found;
	bool contig_v;			// true if vertex indexing is contiguous
	int vec_cv[5];			// bad vertex index +/- 2
	bool contig_f;			// true if vface indexing is contiguous
	int vec_cf[5];			// bad face index +/- 2
	Stats area;
	Stats aspect_ratio;
	Stats edge_length;
	Stats edge_angle;
	Stats adjacent_face;
	vec_e border;		// border edges (edge with a single adjacent face)
	vec_e nonman_e;	// nonmanifold edges (edge with three or more adjacent faces)
	vec_e flipped;		// flipped edges (edges traversed more than once)
	vec_v indistin_v;// indistinguishable vertices 
	vec_e indistin;	// indistinguishable edges (edges with vertices 
									// that are indistinguishbale using double precision)
	vec_v dupl_v_index;// vertices with duplicate vertex indices
	vec_f dupl_f_index;// faces with duplicate face indices
	vec_v nonman_v;	// nonmanifold vertices 
	vec_v orphan;	// vertices not referenced by any face
	vec_f missing_f;	// faces with one or more nonexistant referenced vertices
	vec_i missing_v;	// nonexistant referenced vertices
	vec_f degen;	// faces with two or more identical vertex references
	hmap_fd bad_aspect;	// faces with aspect ratio larger than threshold
	hmap_ed bad_angle;	// edges with angle > max threshold or < min threshold
	hmap_ed bad_length;	// edges with length > max threshold or < min threshold
	// *********************
	void orphanMissingContig(void);
	void degenContig(void);
	void vertexAdjacentFaces(void);
	void areaAspectRatio(Controls&);
	void findIntersectingFaces(Container*,Space&);
	void processEdgeLengths(Controls&);
	void computeEdgeAngles(Controls&);
	void computeVolume(void);
	void printIntegrity(Controls&);
	// *********************
	void analyze(Container*,Controls&,Space&);
	void print(Controls&);
	void printChars(Controls&);
	void checkIntegrity(void);
	void setAll(Vertex*,hmap_fi&,int&);
	void setZero(Edge*,hmap_fi&,int);
	void getGroups(Vertex*,hmap_fi&,set_i&);
	int getLowest(set_i&);
	void replaceGroups(Vertex*,hmap_fi&,int);
    void printAttr(Controls&);
	void vertexDistin(void);
	// *********************
	Object(std::string);
	~Object(void);
	void addOriginal(int,double*);
	void createEdges(void);
	int setNumDigits(void);
	void buildEdge(Face*,Vertex*,Vertex*,map_se&,int);
	Edge* findEdge(Vertex*,Vertex*,map_se&,int);
	void createEdge(Face*,Vertex*,Vertex*,map_se&,int);
	void addVertexPointers(void);
	int getMaxVertex(void);
	void fixFaces(Vertex**);
	void fixEdges(Vertex**);
	void findVertexAdjacencies(void);
	void findNeighborhoods(void);
	void boundObject(double*);
	double getMeanEdgeLength(void);
	//////////////
	bool intersectingFacesExist(Face*);
	bool faceInTable_intf(Face*f);
	ff_iterator findFaceInTable_intf(Face*f);
	//////////////
	bool vertexIsNice(Vertex*);
	int getVertexNiceness(Vertex*);
	void setVertexNiceness(Vertex*,int);
	//////////////
	void addEdgeFlip(Edge*,int);
	void clearFlipTable(void);
	int getEdgeFlip(Edge*);
	//////////////
	bool ifFrozen(hmap_vd&,Vertex*);
	void newFindNeighborhoods(void);
	bool thawedAndAble(hmap_vd&,set_v&);
	void collectFaces(hmap_vd&,set_v&,vec_f&);
	bool processEdge(Edge*,hmap_vd&,vec_e&,Vertex*);	
	//////////////
	void evalCharacteristics(Container*,Controls&,Space&);
	void evalAttributes(Space&);
	bool isClosed(void);
	bool isManifold(void);
	bool isConsistent(void);
	bool isOutward(Space&);
	bool rayIntersectsBB(double[2][3],Face*,double[3]);
	bool rayIntersectsSide(char*,double[2][3],double[6],double[3]);
	void countBoundaries(void);
	int countComponents(void);
	void computeGenus(void);
	bool edgesManifold(bool);
	bool verticesManifold(bool);
	bool removeSelectedFaces(vec_f&,std::vector<Face*>&);
	void getSelectedFaces(vec_f&,std::vector<Face*>);
	void verifyEdges(void);
	bool goodIntegrity(void);
	void findContours(double);
	void addContour(vec_f&,double);
	std::vector<Contour*> contours;
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
	v_iterator i;
	f_iterator j;
	e_iterator k;
	for(i=v.begin();i!=v.end();i++){ delete *i; }
	for(j=f.begin();j!=f.end();j++){ delete *j; }
	for(k=e.begin();k!=e.end();k++){ delete *k; }
	v.clear();
	f.clear();
	e.clear();
}

Object::Object(std::string s) {
	name=s;
    contig_f=contig_v=false;
	closed=consistent=outward=manifold=false;
	vec_cv[0]=vec_cv[1]=vec_cv[2]=vec_cv[3]=vec_cv[4]=0;
	vec_cf[0]=vec_cf[1]=vec_cf[2]=vec_cf[3]=vec_cf[4]=0;
	area.clear();
	aspect_ratio.clear();
	edge_length.clear();
	edge_angle.clear();
	adjacent_face.clear();
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
	bad_aspect.clear();
	bad_angle.clear();
	bad_length.clear();
	dupl_v_index.clear();
	dupl_f_index.clear();
}

bool Object::faceInTable_intf(Face *ff){
	return intf.find(ff)!=intf.end();
}

ff_iterator Object::findFaceInTable_intf(Face *ff){
	return intf.find(ff);
}

void Face::addFaceToTable_intf(void){
	// create new face vector in table
	vec_f *nfv = new std::vector<Face*>();
	v[0]->o->intf[this]=nfv;
}

bool Face::faceInTable_intf(void){
	return v[0]->o->faceInTable_intf(this);
}

ff_iterator Face::findFaceInTable_intf(void){
	return v[0]->o->findFaceInTable_intf(this);
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
		<< vi[0] << "]\n";
	}

	if(v[1]!=NULL){
		cout << "[v1 "
		<< v[1]->index << " "
		<< v[1]->pN[0] << " "
		<< v[1]->pN[1] << " "
		<< v[1]->pN[2] << "]\n";
	} else {
		cout << "[v1 "
		<< vi[1] << "]\n";
	}

	if(v[2]!=NULL){
		cout << "[v2 "
		<< v[2]->index << " "
		<< v[2]->pN[0] << " "
		<< v[2]->pN[1] << " "
		<< v[2]->pN[2] << "]\n";
	} else {
		cout << "[v2 "
		<< vi[2] << "]\n";
	}

	if(e[0]!=NULL){ cout << "[e0 " << e[0] << "]\n";}
	else { cout << "[e0 NULL]\n"; }
	if(e[1]!=NULL){ cout << "[e1 " << e[1] << "]\n";}
	else { cout << "[e1 NULL]\n"; }
	if(e[2]!=NULL){ cout << "[e2 " << e[2] << "]\n";}
	else { cout << "[e2 NULL]\n"; }

}

void Face::printFaceCP(void){
	cout.precision(12);
	if(v[0]!=NULL){
		cout << v[0]->pN[0] << " "
		<< v[0]->pN[1] << " "
		<< v[0]->pN[2] << " 1 0 0 1\n";
	}

	if(v[1]!=NULL){
		cout << v[1]->pN[0] << " "
		<< v[1]->pN[1] << " "
		<< v[1]->pN[2] << " 1 0 0 1\n";
	}

	if(v[2]!=NULL){
		cout << v[2]->pN[0] << " "
		<< v[2]->pN[1] << " "
		<< v[2]->pN[2] << " 1 0 0 1\n";
	}
}

void Face::getVertexCoordinates(double *cpvc[3]){
	cpvc[0]=v[0]->pN;
	cpvc[1]=v[1]->pN;
	cpvc[2]=v[2]->pN;
}

Face::Face(char *triplet,Object *obj){
	
	std::pair<iv_iterator,iv_iterator> pp;
	std::pair<ib_iterator,ib_iterator> qq;

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
	// search for index in multimap
	qq=obj->found.equal_range(zz);
	if(qq.first!=qq.second){
		// set all matching elements to true
		for(ib_iterator k=qq.first;k!=qq.second;k++){
			// set flag to true;
			(*k).second=true;
		}
	}
	pp=obj->vp.equal_range(zz);
	if(pp.first!=pp.second){
	    v[0] = (*(pp.first)).second;
	} else {
	    v[0] = NULL;
		vi[0]=zz;
		// add index to cs.missing_v
		obj->missing_v.push_back(zz);
		// add face to cs.missing_f
		obj->missing_f.push_back(this);
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
	// search for index in multimap
	qq=obj->found.equal_range(zz);
	if(qq.first!=qq.second){
		// set all matching elements to true
		for(ib_iterator k=qq.first;k!=qq.second;k++){
			// set flag to true;
			(*k).second=true;
		}
	}
	pp=obj->vp.equal_range(zz);
	if(pp.first!=pp.second){
	    v[1] = (*(pp.first)).second;
	} else {
	    v[1] = NULL;
		vi[1]=zz;
		// add index to cs.missing_v
		obj->missing_v.push_back(zz);
		// add face to cs.missing_f
		obj->missing_f.push_back(this);
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
	// search for index in multimap
	qq=obj->found.equal_range(zz);
	if(qq.first!=qq.second){
		// set all matching elements to true
		for(ib_iterator k=qq.first;k!=qq.second;k++){
			// set flag to true;
			(*k).second=true;
		}
	}
	pp=obj->vp.equal_range(zz);
	if(pp.first!=pp.second){
	    v[2] = (*(pp.first)).second;
	} else {
	    v[2] = NULL;
		vi[2]=zz;
		// add index to cs.missing_v
		obj->missing_v.push_back(zz);
		// add face to cs.missing_f
		obj->missing_f.push_back(this);
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
	Vertex *eg;		// example vertex
	vec_o o;						//vector of object pointers
	std::vector<std::string> files;	// array of input file names
	int num_files;					// number of input files
	int pairs[3][2];
	///////////
	~Container(void);
	Container(void);
	///////////
	void update(Object*);
	void clear(void);
	int getScore(void);
	void printBatch(Controls&);
	void scanDir(const char*);
//	Object* processFile(char[1024]);
	Object* processFile(std::string);
//	void scanFile(Object*,char*);
	void scanFile(Object*,std::string);
	void createEdges(void);
	void addVertexPointers(void);
	void findVertexAdjacencies(void);
	void writeDistances(void);
	///////////
	void boundWorld(Space&);
	void getNiceSet(Space&,Monitor&);
	void findNice(Space&);
	bool checkNiceness(Space&,Vertex*);
	void collectCrossed(Space&,Vertex*,vec_o&);
	bool updateNiceness(Vertex*,vec_o&);	
	void findClosestAxis(Space&,Vertex*,double[2][3]);
	int findExtraPoint(Space&,Vertex*,double[3],int);
	void findCrossed1(Space&,Vertex*,double[2][3],vec_o&);
	void findCrossed2(Space&,double[2][3],vec_o&);
	void getExtraRay(Vertex*,double[2][3],int);
	void collectNiceFaces(Space&,double[2][3],vector<Face*>&);
	void getBoxIndexRange(Space&,double[2][3],int[6]);
	void findIntersectedFaces(double[2][3],vector<Face*>&,vector<Face*>&,vector<int>&);
	void findOddMeshes(vector<Face*>&,vector<int>&,int&,vec_o&);
	///////////
	void getSeparationDistances(Space&);
	bool findClosest(Space&,Vertex*);
	bool computeClosest(Face*,Vertex*,double&,double[3]);
	bool getPlaneIntersection(Face*,Vertex*,double*,double,double,Point&);
	void getEdgeIntersection(Vertex*,double*[3],Point&);
	void getBoxes(vector<Box*>&,Vertex*,int,Space&);
	void getCandidateFaces(vector<Box*>&,Vertex*,hset_f&);
	bool faceInNeighborhood(Face*,Vertex*);
	///////////
	bool checkFaceFaceIntersections(Face*,Face*);
	int checkEdgeEdgeIntersection(Face*,Face*,bool);
	int checkFaceEdgeIntersection(Face*,Face*);
	bool facesParallel(Face*,Face*);
	bool facesColinear(Face*,Face*);
	int numUniqueVertices(Face*,Face*,int[2]);
	void assignFacesToBoxes(Space&);
	///////////
	int num_orph;
	int countOrphan(void);
	int num_mis;
	int countMissing(void);
	int num_deg;
	int countDegen(void);
	int num_bor;
	int countBorder(void);
	int num_flip;
	int countFlipped(void);
	int num_nonman_e;
	int countNonmanE(void);
	int num_nonman_v;
	int countNonmanV(void);
	int num_obj;
	int countObject(void);
	int num_vert;
	int countVertex(void);
	int num_face;
	int countFace(void);
	int num_edge;
	int countEdge(void);
	double num_vol;
	double countVol(void);
	double countArea(void);
	int num_sep;	// yep, sep for separate components, not a mistake
	int countComponents(void);
	int num_bou;
	int countBoundaries(void);
	int num_indistin;
	int countIndistin(void);
	int countIntFace(void);
	int num_dupl_v;
	int countDuplV(void);
	int num_dupl_f;
	int countDuplF(void);
	///////////
	int num_clo_cc,num_clo_nn;
	std::pair<int,int> countClosed(void);
//	int num_man_cc,num_man_nn;
//	std::pair<int,int> countManifold(void);
	int num_man[3];
	void countManifold(int[3]);
	int num_cons[3];
	void countConsistent(int[3]);
	int num_out[3];
	void countOutward(int[3]);
};

Container::Container(void){
	num_obj=num_vert=num_face=num_edge=num_sep=num_flip=0;
//	num_man_nn=num_man_cc=num_clo_nn=num_clo_cc=0;
	num_clo_nn=num_clo_cc=0;
	num_orph=num_mis=num_deg=num_dupl_v=num_dupl_f=0;
	num_bou=num_indistin=num_bor=num_nonman_e=num_nonman_v=0;
	num_vol=0.0;
	num_man[0]=num_man[1]=num_man[2]=0;
	num_cons[0]=num_cons[1]=num_cons[2]=0;
	num_out[0]=num_out[1]=num_out[2]=0;
	pairs[0][0] = 0;
	pairs[0][1] = 1;
	pairs[1][0] = 1;
	pairs[1][1] = 2;
	pairs[2][0] = 2;
	pairs[2][1] = 0;
	num_files=0;
}

int Container::countIntFace(void){
	int a=0;
	for(o_iterator i=o.begin();i!=o.end();i++){
		a+=(*i)->intf.size();
	}
	return a;
}

void Container::countOutward(int val[3]){
	//			consistent	outward
	// val[0] 	true		true
	// val[1] 	true		false
	// val[2] 	false		true or false
//	val[0]=0; val[1]=0; val[2]=0;
//	for(o_iterator i=o.begin();i!=o.end();i++){
//		if     ((*i)->consistent==true && (*i)->outward==true) {val[0]++;}
//		else if((*i)->consistent==true && (*i)->outward==false){val[1]++;}
//		else {val[2]++;}
//	}
	val[0]=num_out[0];
	val[1]=num_out[1];
	val[2]=num_out[2];
}

void Container::countConsistent(int val[3]){
	//			manifold	consistent
	// val[0] 	true		true
	// val[1] 	true		false
	// val[2] 	false		true or false
//	val[0]=0; val[1]=0; val[2]=0;
//	for(o_iterator i=o.begin();i!=o.end();i++){
//		if     ((*i)->manifold==true && (*i)->consistent==true) {val[0]++;}
//		else if((*i)->manifold==true && (*i)->consistent==false){val[1]++;}
//		else {val[2]++;}
//	}
	val[0]=num_cons[0];
	val[1]=num_cons[1];
	val[2]=num_cons[2];
}

void Container::countManifold(int val[3]){
//	int cc=0,nn=0; // cc=manifold, nn=not manifold
//	for(o_iterator i=o.begin();i!=o.end();i++){
//		if((*i)->manifold==true){cc++;}
//		else {nn++;}
//	}
//	return std::make_pair(cc,nn);
	//			close		manifold
	// val[0] 	true		true			'manifold'
	// val[1] 	true		false			'nonmanifold'
	// val[2] 	false		true or false	'undefined'
	val[0]=num_man[0];
	val[1]=num_man[1];
	val[2]=num_man[2];
}
/*
std::pair<int,int> Container::countManifold(void){
//	int cc=0,nn=0; // cc=manifold, nn=not manifold
//	for(o_iterator i=o.begin();i!=o.end();i++){
//		if((*i)->manifold==true){cc++;}
//		else {nn++;}
//	}
//	return std::make_pair(cc,nn);
	return std::make_pair(num_man_cc,num_man_nn);
}*/

std::pair<int,int> Container::countClosed(void){
//	int cc=0,nn=0; // cc=closed, nn=not closed
//	for(o_iterator i=o.begin();i!=o.end();i++){
//		if((*i)->closed==true){cc++;}
//		else {nn++;}
//	}
//	return std::make_pair(cc,nn);
	return std::make_pair(num_clo_cc,num_clo_nn);
}

int Container::countDuplV(void){
//	int a=0;
//	for(o_iterator i=o.begin();i!=o.end();i++){
//		a+=(*i)->dupl_v_index.size();
//	}
//	return a;
	return num_dupl_v;
}

int Container::countDuplF(void){
//	int a=0;
//	for(o_iterator i=o.begin();i!=o.end();i++){
//		a+=(*i)->dupl_f_index.size();
//	}
//	return a;
	return num_dupl_f;
}

int Container::countIndistin(void){
//	int a=0;
//	for(o_iterator i=o.begin();i!=o.end();i++){
//		a+=(*i)->indistin_v.size();
//	}
//	return a;
	return num_indistin;
}

int Container::countBoundaries(void){
//	int a=0;
//	for(o_iterator i=o.begin();i!=o.end();i++){
//		a+=(*i)->num_bou;
//	}
//	return a;
	return num_bou;
}

int Container::countComponents(void){
//	int a=0;
//	for(o_iterator i=o.begin();i!=o.end();i++){
//		a+=(*i)->num_sep;
//	}
//	return a;
	return num_sep;
}

double Container::countArea(void){
	double a=0.0;
	for(o_iterator i=o.begin();i!=o.end();i++){
		a+=(*i)->area.sum;
	}
	return a;
}

double Container::countVol(void){
//	double a=0.0;
//	for(o_iterator i=o.begin();i!=o.end();i++){
//		a+=(*i)->vol;
//	}
//	return a;
	return num_vol;
}

int Container::countEdge(void){
//	int a=0;
//	for(o_iterator i=o.begin();i!=o.end();i++){
//		a+=(*i)->e.size();
//	}
//	return a;
	return num_edge;
}

int Container::countFace(void){
//	int a=0;
//	for(o_iterator i=o.begin();i!=o.end();i++){
//		a+=(*i)->f.size();
//	}
//	return a;
	return num_face;
}

int Container::countVertex(void){
//	int a=0;
//	for(o_iterator i=o.begin();i!=o.end();i++){
//		a+=(*i)->v.size();
//	}
//	return a;
	return num_vert;
}

int Container::countObject(void){
//	return o.size();
	return num_obj;
}

int Container::countNonmanV(void){
//	int a=0;
//	for(o_iterator i=o.begin();i!=o.end();i++){
//		a+=(*i)->nonman_v.size();
//	}
//	return a;
	return num_nonman_v;
}

int Container::countNonmanE(void){
//	int a=0;
//	for(o_iterator i=o.begin();i!=o.end();i++){
//		a+=(*i)->nonman_e.size();
//	}
//	return a;
	return num_nonman_e;
}

int Container::countFlipped(void){
//	int a=0;
//	for(o_iterator i=o.begin();i!=o.end();i++){
//		a+=(*i)->flipped.size();
//	}
//	return a;
	return num_flip;
}

int Container::countBorder(void){
//	int a=0;
//	for(o_iterator i=o.begin();i!=o.end();i++){
//		a+=(*i)->border.size();
//	}
//	return a;
	return num_bor;
}

int Container::countDegen(void){
//	int a=0;
//	for(o_iterator i=o.begin();i!=o.end();i++){
//		a+=(*i)->degen.size();
//	}
//	return a;
	return num_deg;
}

int Container::countMissing(void){
//	int a=0;
//	for(o_iterator i=o.begin();i!=o.end();i++){
//		a+=(*i)->missing_v.size();
//	}
//	return a;
	return num_mis;
}

int Container::countOrphan(void){
//	int a=0;
//	for(o_iterator i=o.begin();i!=o.end();i++){
//		a+=(*i)->orphan.size();
//	}
//	return a;
	return num_orph;
}

Container::~Container(void){
	o_iterator i;
	for(i=o.begin();i!=o.end();i++){ delete *i; }
}

void Container::clear(void){
	o.clear();
}

// ######################################
// ######################################

class Box {
public:
	vector<Face*> f;			// vector of face pointers
	int x,y,z;					// indices of box in space
	Box(int,int,int);
	~Box(void);
	double xmin(double,double);
	double xmax(double,double);
	double ymin(double,double);
	double ymax(double,double);
	double zmin(double,double);
	double zmax(double,double);
	void getFaceIntersection(Container*);
	bool faceIntersectionAlreadyKnown(Face*,Face*);
	void printBox(Space*);
};

Box::~Box(void){
	f.clear();
}

double Box::xmin(double i,double sl){return x*sl+i;}
double Box::xmax(double i,double sl){return (x+1)*sl+i;}
double Box::ymin(double i,double sl){return y*sl+i;}
double Box::ymax(double i,double sl){return (y+1)*sl+i;}
double Box::zmin(double i,double sl){return z*sl+i;}
double Box::zmax(double i,double sl){return (z+1)*sl+i;}

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
	int num_boxes;			// WHY DO I NEED THIS?
	double space_length;
	vec_b b;				// vector of boxes
	double world[6];		// minX, maxX, minY, maxY, minZ, maxZ of world
	// *****************
	void initBoxes(double);
	void clearBoxes(void);
	void deleteBoxes(void);
	void recordFace(vec_b&,Face*);
	void computeBoxesToCheck(Face*,vec_b&);
	void index2Range(int,double[2],char*);
	int location2Index(double,char*);
	int indices2Index(int,int,int);
	void getBoxesFor3DLocations(double[6],vec_b&);
	void getBoxesFor3DIndices(int[6],vec_b&);
	int screenIndex(int,char*);
	~Space(void);
};

Space::~Space(void){
	b_iterator i;
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

void Space::getBoxesFor3DLocations(double dr[6],vec_b& bp){
	// assume p = [xmin,xmax,ymin,ymax,zmin,zmax]
	int ir[6]; // index range
	// compute 3D index range that includes 3D location range
	ir[0] = location2Index(dr[0],"x");   // -x
	ir[1] = location2Index(dr[1],"x");	//  x
	ir[2] = location2Index(dr[2],"y");	// -y
	ir[3] = location2Index(dr[3],"y");	//  y
	ir[4] = location2Index(dr[4],"z");   // -z
	ir[5] = location2Index(dr[5],"z");	//  z
/*	// DEBUG
	cout << "\nSpace::getBoxesFor3DLocations: "
	<< "data range ["
	<< dr[0] << ","
	<< dr[1] << ","
	<< dr[2] << ","
	<< dr[3] << ","
	<< dr[4] << ","
	<< dr[5] << "]\n";
	cout << "\nSpace::getBoxesFor3DLocations: "
	<< "num_space ["
	<< num_space[0] << ","
	<< num_space[1] << ","
	<< num_space[2] << "]\n";
	cout << "\nSpace::getBoxesFor3DLocations: "
	<< "world ["
	<< world[0] << ","
	<< world[1] << ","
	<< world[2] << ","
	<< world[3] << ","
	<< world[4] << ","
	<< world[5] << "]\n";
	cout << "\nSpace::getBoxesFor3DLocations: "
	<< "index range ["
	<< ir[0] << ","
	<< ir[1] << ","
	<< ir[2] << ","
	<< ir[3] << ","
	<< ir[4] << ","
	<< ir[5] << "]\n";
*/	// DEBUG
	getBoxesFor3DIndices(ir,bp);
}

void Space::getBoxesFor3DIndices(int br[6],vec_b& bp){
	// assume br = [xmin,xmax,ymin,ymax,zmin,zmax]
	int x0=br[0],x1=br[1],y0=br[2],y1=br[3],z0=br[4];
	b_iterator i;
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
//	return z*num_space[0]*num_space[1]+y*num_space[0]+x;
	return num_space[0]*(z*num_space[1]+y)+x;
}

void Space::index2Range(int i,double r[2],char *c){
		r[0] = i*space_length;
		r[1] = (i+1)*space_length;
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
		a = (int) floor((ss-world[0])/space_length);
//		cout << "\nSpace::location2Index: "
//		<< "x index (proposed,accepted) = " << a;
		a = screenIndex(a,"x");
//		cout << "," << a << endl;
	} else if (!strcmp(c,"y")){
		a = (int) floor( (ss-world[2])/space_length );
//		cout << "\nSpace::location2Index: "
//		<< "y index (proposed,accepted) = " << a;
		a = screenIndex(a,"y");
//		cout << "," << a << endl;
	} else if (!strcmp(c,"z")){
		a = (int) floor( (ss-world[4])/space_length );
//		cout << "\nSpace::location2Index: "
//		<< "z index (proposed,accepted) = " << a;
		a = screenIndex(a,"z");
//		cout << "," << a << endl;
	} else {cout << "Received unexpected string.\n"; exit(0);}
	return a;
}

void Box::printBox(Space *s){
	cout << "Box indices ["
	<< x << " "
	<< y << " "
	<< z << "]\n"
	<< "Box range ["
	<< xmin(s->world[0],s->space_length) << " "
	<< xmax(s->world[0],s->space_length) << " "
	<< ymin(s->world[2],s->space_length) << " "
	<< ymax(s->world[2],s->space_length) << " "
	<< zmin(s->world[4],s->space_length) << " "
	<< zmax(s->world[4],s->space_length) << "]";
}

//############################################################################
//############################################################################
class Face_Pair
{
public:
	Vertex *a,*b,*c,*d;	// vertex*s
	Edge *e;			// shared edge
	Face *f1,*f2;		// bad faces
	Object *o;			// parent object of faces
	int next_i;			// next available face index
	Face_Pair(Object*);
	void clear(void);
	void processBadFace(Face*);
	void print(void);
	void analyzeF1(void);
	void findF2(void);
	bool aspectRatiosImprove(void);
	bool existingEdge(void);
};

bool Face_Pair::existingEdge(void){
	// for each edge in object
	for(e_iterator i=o->e.begin();i!=o->e.end();i++){
		// if edge vertices are b and d
		if( ((*i)->vv1==b && (*i)->vv2==d) || 
			((*i)->vv1==d && (*i)->vv2==b)){return true;}
	}
	return false;
}

bool Face_Pair::aspectRatiosImprove(void){
	double arf1 = f1->getAspectRatio();
	double arf2 = f2->getAspectRatio();
	Face nf1(a,b,d);
	Face nf2(b,c,d);
	double arnf1 = nf1.getAspectRatio();
	double arnf2 = nf2.getAspectRatio();

	double old_max,new_max;
	if(arf1>arf2){old_max=arf1;}
	else {old_max=arf2;}
	if(arnf1>arnf2){new_max=arnf1;}
	else {new_max=arnf2;}

	if(new_max<old_max) {return true;}
	else {return false;}
}

void Face_Pair::print(void){
	cout << "Follow these instructions to improve face aspect ratio.\n"
	<< "open " << o->name << endl
	<< "remove Face "
	<< f1->index << " "
	<< f1->v[0]->index << " "
	<< f1->v[1]->index << " "
	<< f1->v[2]->index << endl
	<< "remove Face "
	<< f2->index << " "
	<< f2->v[0]->index << " "
	<< f2->v[1]->index << " "
	<< f2->v[2]->index << endl
	<< "add Face "
	<< next_i++ << " "
	<< a->index << " "
	<< b->index << " "
	<< d->index << endl
	<< "add Face "
	<< next_i++ << " "
	<< b->index << " "
	<< c->index << " "
	<< d->index << endl;
}

void Face_Pair::findF2(void){
	// for each face in object
	for(f_iterator i=o->f.begin();i!=o->f.end();i++){
		// if face not f1
		if(*i!=f1){
			// if face contains edge e
			if((*i)->e[0]==e||(*i)->e[1]==e||(*i)->e[2]==e){
				f2=*i;
				break;
			}
		}
	}
	//check
	if(f2==NULL){
		cout << "Face_Pair::findF2: Unable to identify f2.\n";
		exit(0);
	}
	// identify d
	if(f2->v[0]!=a && f2->v[0]!=c){d=f2->v[0];}
	else if(f2->v[1]!=a && f2->v[1]!=c){d=f2->v[1];}
	else if(f2->v[2]!=a && f2->v[2]!=c){d=f2->v[2];}
	else {cout << "Face_Pair::findF2: Error: Unable to identify d.\n";
		exit(0);
	}
}

void Face_Pair::analyzeF1(void){
	double ll = -1.0;
	// find longest edge
	for(int i=0;i<3;i++){
		if(f1->e[i]->l>ll){
			e=f1->e[i];
			ll=f1->e[i]->l;
		}
	}
	// identify a,b,c
	if		 (f1->v[0]==e->vv1 && f1->v[1]==e->vv2){
		a=e->vv2;
		c=e->vv1;
		b=f1->v[2];
	} else if(f1->v[1]==e->vv1 && f1->v[2]==e->vv2){
		a=e->vv2;
		c=e->vv1;
		b=f1->v[0];
	} else if(f1->v[2]==e->vv1 && f1->v[0]==e->vv2){
		a=e->vv2;
		c=e->vv1;
		b=f1->v[1];
	} else if  (f1->v[0]==e->vv2 && f1->v[1]==e->vv1){
		a=e->vv1;
		c=e->vv2;
		b=f1->v[2];
	} else if(f1->v[1]==e->vv2 && f1->v[2]==e->vv1){
		a=e->vv1;
		c=e->vv2;
		b=f1->v[0];
	} else if(f1->v[2]==e->vv2 && f1->v[0]==e->vv1){
		a=e->vv1;
		c=e->vv2;
		b=f1->v[1];
	} else {
		cout << "Face_Pair::analyzeF1: Error: Unable to identify a,b,c.\n";
		exit(0);
	}
}

Face_Pair::Face_Pair(Object *oo){
	int max=-1;
	o=oo;
	// for each face in object
	for(f_iterator i=o->f.begin();i!=o->f.end();i++){
		if((*i)->index > max){max=(*i)->index;}
	}
	next_i=max+1;
	//
	a=b=c=d=NULL;
	f1=f2=NULL;
	e=NULL;
}

void Face_Pair::clear(void){
	a=b=c=d=NULL;
	f1=f2=NULL;
	e=NULL;
}

void Face_Pair::processBadFace(Face *f){
	clear();
	f1=f;
	analyzeF1();
	findF2();
	// if reconfiguring the faces would improve the face aspect ratios
	// i.e. the largest aspect ratio would get smaller
	if(aspectRatiosImprove()==true && existingEdge()==false){
		print();
	}
}

bool Edge::valid(void){
	// if edge has two adjacent faces
	if(f1!=NULL && f2!=NULL){
		// and the faces are of interest
		if(	(f1->index==91824 && f2->index==91825 ) || 
			(f1->index==91825 && f2->index==91824 ) ){
			cout << "Edge::valid: <" << this << "> vv1=" << vv1->index 
			<< ", vv2=" << vv2->index << endl;
			printEdge(vv1->o->name);
			// if the edge vertex indices match the mesh
			if ((vv1->index==9430 && vv2->index==9632)==false &&
				(vv1->index==9632 && vv2->index==9430)==false 
				){return false;}
		}

	}
	return true;
}

void Edge::printEdgeCP(void){
	Vertex *v1=NULL,*v2=NULL,*o1=NULL,*o2=NULL;
	getVertices(v1,v2,o1,o2);
	if((v1!=vv1) || (v2!=vv2)){
		cout << "Edge::printEdge: "
		<< "vertices don't match:\n"
		<< "	v1:\n";
		v1->printVertex(v1->o->name);
		cout << endl << "	v2:\n";
		v2->printVertex(v2->o->name);
		cout << endl << "	o1:\n";
		o1->printVertex(o1->o->name);
		cout << endl << "	o2:\n";
		o2->printVertex(o2->o->name);
		cout << endl << "	vv1:\n";
		vv1->printVertex(vv1->o->name);
		cout << endl << "	vv2:\n";
		vv2->printVertex(vv2->o->name);
		cout << endl;
		exit(0);
	}
	cout.precision(12);
	cout << "printEdge: " << this << endl;
	cout << "printEdge: <obj>" << v1->o->name << endl;
	cout << v1->pN[0] << " "
	<< v1->pN[1] << " "
	<< v1->pN[2] << " 1 0 0 1\n";
	cout << v2->pN[0] << " "
	<< v2->pN[1] << " "
	<< v2->pN[2] << " 1 0 0 1\n";
}

void Edge::printEdge(std::string s){
	Vertex *v1=NULL,*v2=NULL,*o1=NULL,*o2=NULL;
	getVertices(v1,v2,o1,o2);
	if((v1!=vv1) || (v2!=vv2)){
		cout << "Edge::printEdge: "
		<< "vertices don't match:\n"
		<< "	v1:\n";
		v1->printVertex(v1->o->name);
		cout << endl << "	v2:\n";
		v2->printVertex(v2->o->name);
		cout << endl << "	o1:\n";
		o1->printVertex(o1->o->name);
		cout << endl << "	o2:\n";
		o2->printVertex(o2->o->name);
		cout << endl << "	vv1:\n";
		vv1->printVertex(vv1->o->name);
		cout << endl << "	vv2:\n";
		vv2->printVertex(vv2->o->name);
		cout << endl;
		exit(0);
	}
	cout.precision(12);
	cout << "printEdge: " << this << endl;
	cout << "printEdge: <obj>" << s << endl;
	cout << "printEdge:"
	<< " v1 "<< v1->index << " ["
	<< v1->pN[0] << " "
	<< v1->pN[1] << " "
	<< v1->pN[2] << "]\n";
	cout << "printEdge:"
	<< " v2 "<< v2->index << " ["
	<< v2->pN[0] << " "
	<< v2->pN[1] << " "
	<< v2->pN[2] << "]\n"
	<< "printEdge:"
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

class Boundary
{
public:
	vec_e e;
	Vertex *end,*begin;
	bool open;
	Boundary(void);
	void init(Edge*);
	bool edgeExtendsBoundary(Edge*);
	bool add(Edge*);
	void print(void);
	bool closed(void);
};

bool Boundary::closed(void){
	if(begin==end){
		open=false;
		return true;
	} else {
		return false;
	}
}

void Boundary::print(void){
	cout << "\nBEGIN VERTEX:\n";
	begin->printVertex(begin->o->name);
	cout << endl;
	cout << "END VERTEX:\n";
	end->printVertex(end->o->name);
	cout << endl;
	for(e_iterator i=e.begin();i!=e.end();i++){
		(*i)->printEdge((*i)->f1->v[0]->o->name);
		cout << endl;
	}
}

Boundary::Boundary(void){
	begin=end=NULL;
	open=false;
}

void Boundary::init(Edge *ee){
	Vertex *v1=NULL,*v2=NULL,*o1=NULL,*o2=NULL;
	ee->getVertices(v1,v2,o1,o2);
	end=v1;
	begin=v2;
	open=true;
	e.clear();
}

bool Boundary::add(Edge *ee){
	// if edge not found in vector
	if(find(e.begin(),e.end(),ee)==e.end()){
		e.push_back(ee);
		return true;
	} else {
		return false;
	}
}

bool Boundary::edgeExtendsBoundary(Edge *ee){
	Vertex *v1=NULL,*v2=NULL,*o1=NULL,*o2=NULL;
	ee->getVertices(v1,v2,o1,o2);
	if(v1==begin){
		if(add(ee)){ begin=v2; }
		else { cout << "Error: tried to add edge to boundary twice!\n";exit(0);}
		return true;
	} else if(v2==begin) {
		if(add(ee)){ begin=v1; }
		else { cout << "Error: tried to add edge to boundary twice!\n";exit(0);}
		return true;
	} else if(v1==end) {
		if(add(ee)){ end=v2; }
		else { cout << "Error: tried to add edge to boundary twice!\n";exit(0);}
		return true;
	} else if(v2==end) {
		if(add(ee)){ end=v1; }
		else { cout << "Error: tried to add edge to boundary twice!\n";exit(0);}
		return true;
	}
	return false;
}

// #####################################################
// #####################################################

void Contour::print(Controls &cs,std::string name){
	if(strcmp(cs.style,"cp")==false){
		// for each vertex* in contour
		for(v_iterator i=p.begin();i!=p.end();i++){
			cout << (*i)->pN[0] << " "
			<< (*i)->pN[1] << " "
			<< (*i)->pN[2] << " 1 0 0 1\n";
		}
	} else {
		cout << "<Contour name=\"" << name
		<< "\" hidden=\"false\" closed=\"true\" "
		<< "simplified=\"true\" border=\"0 0 1\" "
		<< "fill=\"1 1 0.501961\" mode=\"9\"\n"
		<< " points=\"";
		// for each vertex* in contour
		for(v_iterator i=p.begin();i!=p.end();i++){
			cout << (*i)->pN[0] << " "
			<< (*i)->pN[1] << ",\n";
		}
		cout << "\"/>\n";
	}
}
	
