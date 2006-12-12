// ######################################
// ######################################

class Face {
public:
	int index;		// face index
	Object *o;		// pointer to parent object
	Vertex *v[3];	// pointers to vertices
	Edge   *e[3];	// pointers to edges
	double  n[3];	// outward normal vector
	double  iv[3];	// intersection force vector
	bool intersection;
	int checkEdgeEdgeIntersection(double*[3],int[3][2],double*[3],int,int,bool);
	int checkPolygonEdgeIntersection(double*[3],int[3][2],double[3],double*[3],int,int);
	int* getVertexIndices(int*);
	int* getVertexCoordinates(int*);
	Face(void); 
};

Face::Face(void){
	o=NULL;
	v[0]=v[1]=v[2]=NULL;
	e[0]=e[1]=e[2]=NULL;
	n[0]=n[1]=n[2]=0;
	iv[0]=iv[1]=iv[2]=0;
	intersection=false;
	index=0;
}

// ######################################
// ######################################

class Box {
public:
	vector<Face*> f;			// vector of face pointers
	vector<Face*>::iterator i;	// vector of face pointers
	int x,y,z;					// indices of box in space
	void addFace(Face* fp);
	void clearFaces(void);
	void setIndices(int a,int b,int c);
	Box(void);
};

Box::Box(void){
	x=y=z=-1;
}

Box::addFace(Face *fp){
	f.insert(f.begin(),fp);
}

Box::clearFaces(void){
//	f.erase(f.begin(),f.end());
	f.clear();
}

Box::setIndices(int a,int b,int c){
	x=a;
	y=b;
	z=c;
}

// ######################################
// ######################################

class Edge {
public:
	Object *o;		// pointer to parent object
	Vertex *v1,*v2; // pointers to vertices on edge
	Vertex *o1,*o2; // pointers to vertices on adjacent faces not on edge
	Face *f1,*f2;	// pointers to adjacent faces (i.e. faces that contain edge)
	double angle_cosine; // cosine of angle between adjacent faces
	double l;		// original edge length
	double force[2][3];
	int flip; // indicates direction of force relative to normal 
				// (same,+1 or opposite,-1)
	void initAngles(PolygonClass*);
	void getAngleForceEnergy(PolygonClass*);
	double getLength(void);
	double getEnergy(void);
	double* getForce(double*);
	Edge(void);
};

Edge::Edge(void){
	o=NULL;
	v1=v2=NULL;
	o1=o2=NULL;
	f1=f2=NULL;
	angle_cosine=l=0;
	flip=0;
	for(int i=0;i<3;i++){
		force[0][i]=force[1][i]=0;
	}
}

// ######################################
// ######################################

class Vertex {
	public:
		ofstream Vertexfile;
		int index;
		double p0[3];		// original position coordinates (x,y,z)
		double pN[3];		// current position coordinates (x,y,z)
		double pC[3];		// closest mesh position coordinates (x,y,z)
		bool nice;			// true if vertex is not inside another object, false otherwise
		bool candidate;		// true if vertex is candidate for being moved
		bool frozen;		// if true, then do not move vertex
		vector<Edge*> e;	// pointers to adjacent edges
		vector<Face*> f;	// pointers to adjacent faces
		hash<Vertex*> v;	// pointers to neighborhood vertices
		int getHashVal(void);
		double* getNormal(double*);
		double* getForce(double*); 
		double getEnergy(void); 
		double getDevDist(void); 
		double getPosErr(void); 
		void setNiceness(int);
		void getNewCoords(double); 
		void fileInit (int); 
		void fileOutit (void); 
		void updateLog (int);
		void addVertex (int&); 
		void addEdge (int&); 
		void addFace (int&); 
		void addSortUniqueNeighbor (int&); 
		Vertex(void);
};

Vertex::Vertex(void) {
	index=0;
	p0[0]=p0[1]=p0[2]=0;
	pN[0]=pN[1]=pN[2]=0;
	pC[0]=pC[1]=pC[2]=0;
	nice=candidate=frozen=false;
}

// ######################################
// ######################################

class Object
{
public:
	ofstream Objectfile;
	int index;
	bool frozen;	// if true, object vertices will not move
	char name[SIZE_OF_OBJECT_NAME];	// object name
	double getPosErr(void);	// get cumulative position error of all vertices in object
	double getForce(void);	// get cumulative force of all vertices in object
	double getEnergy(void);	// get cumulative energy of all vertices in object
	int getNonnice(void);	// get cumulative number of all vertices in object
	vector<Vertex*> v;		// container of pointers to all vertices in object
	vector<Face*> f;		// container of pointers to all faces in object
	vector<Edge*> e;		// container of pointers to all edges in object
	int getNumEdges(void);	// return number of edges in object
	int getNumFaces(void);	// return number of faces in object
	int getNumVertices(void);	// return number of vertices in object
	void fileInit (void); 
	void fileOutit (void); 
	void InitVertices (void); 
	void OutitVertices (void); 
	void makeObject(int ,int );
	void addOriginalVertexCoordinates(int ,double *);
	void initializeNewCoordinatesAndDevDist(void);
	void getPolygonCoordinates(int [3],double* [3]);
	void findAdjacentPolygonsAndVertices(PolygonClass*, int, int *);
	void findNeighborhoodVertices(void);
	void computeEdgeLengths(void);
	void computeNewEdgeLengths(void);
	void updateObjectLog(int);
	void updateVertexLog(int);
	void computeGlobalParams(void);
	void initEdges(void);
	void buildEdges(PolygonClass*, int, int *);
	void initEdgeAngles(PolygonClass*);
	void computeEdgeAngleForceEnergy(PolygonClass*);
	void boundObject(double*);
	~ObjectClass(void);
};

Object::~Object(void) {

	delete[] vertices;
	delete[] edges;

}

// ######################################
// ######################################

class ContainerClass
{
	public:
		ofstream Containerfile;
		ofstream ObjectList;
		int object_count;
		int nonnice;
		double position_error;	
		double force;
		double energy;
		ObjectClass *objects;
		void init(void);
		void initObjects(int *,int *);
		void initObjectAndVertexLogFiles(void);
		void OutitObjectAndVertexLogFiles(void);
		void initCoordsAndDist(void);
		void initFindAdj(PolygonClass*,int,int *);
		void initFindNeighbor(void);
		void computeEdges(void);
		void computeNewEdges(void);
		void updateObjectAndVertexLog(int);
		void objectList(void);
		void fileInit(void);
		void fileOutit(void);
		void updateFile(int);
		void computeGlobalParams(void);
		void initObjectEdges(PolygonClass*,int,int *);
		void initObjectEdgeAngles(PolygonClass*);
		void computeObjectEdgeAngleForceEnergy(PolygonClass*);
		void writeVertexClosestDistanceToFile(void);
		~ContainerClass(void);
};

ContainerClass::~ContainerClass(void) {

	delete[] objects;

}

// ######################################
// ######################################

class ScanClass
{
	public:
		ofstream Scanfile;
		int num_files;							// number of input files
		char files[NUMBER_OF_INPUT_FILES][SIZE_OF_OBJECT_NAME];	// array of input file names
		int object_count;						// total number of objects in model
		int vertex_count;						// total number of vertices in model
		int polygon_count;						// total number of polygons in model
		int vertices_array[NUMBER_OF_OBJECTS];	// array of number of vertices per object
		int polygons_array[NUMBER_OF_OBJECTS];	// array of number of polygons per object
		void scanDir(void);
		void scanFiles(ContainerClass&,PolygonClass*, bool);
		void scanFile(ContainerClass&,PolygonClass*, char [32],bool, int);
		void verifyPolygons(ContainerClass&,PolygonClass*);
		void verifyVertices(ContainerClass&);
		void verifyAdjacencies(ContainerClass&);
		void buildMeshAfter(ContainerClass&,PolygonClass*,int);
		void writePolygonsPerBoxHistogram(BoxClass*,int);
		void fileInit(void);
		void fileOutit(void);
		void updateFile(void);
};

// ######################################
// ######################################

// definition of box element
class SpaceClass {
	public:
		ofstream Spacefile;
		int num_space[3];
		int num_boxes;
		BoxClass *boxes;
		double world[6];				// minX, maxX, minY, maxY, minZ, maxZ of world
		void init(void);
		void clear(void);
		void boundWorld(ContainerClass&,int);
		void computeBoxesToCheck(double* [3],int*&,int&);
		void checkVertexOverlap(double [3][3],int,bool*);
		void checkPolygonPlaneCrossesBox(int,double [3][3] ,double*,bool*,bool*);
		void checkPolygonEdgeIntersectBoxFace(int,double[3][3],bool*,bool*);
		void recordPolygon(int*,int,int);
		void checkBoxOverlap(double* [3],int,bool*);
		void fileInit(void);
		void fileOutit(void);
		void updateFile(int);
		~SpaceClass();
};

SpaceClass::~SpaceClass(void) {

	delete[] boxes;
}

// ######################################
// ######################################

class ManipClass
{
	public:
		ofstream Manipfile;
		int object_count;				// total number of objects in model
		int vertex_count;				// total number of vertices in model
		int polygon_count;				// total number of polygons in model
		bool self_intersection;			// flag signifies if self-intersecting object was found
		void fileInit(void);
		void fileOutit(void);
		void computePolygonNormals(ContainerClass&,PolygonClass*);

///////////

		void findNice(ContainerClass&,PolygonClass*,SpaceClass&);
		void findClosestAxis(double[3],double[3]);
		void collectNicePolygons(ContainerClass&,PolygonClass*,SpaceClass& 
								,double [3],int,int*&,int&,int&,int,bool);
		void findIntersectedPolygons(ContainerClass&,PolygonClass*,SpaceClass&,
									int,int,double[3],int*&,int&,
									int*&,int&,int&,
									int*&,int&,int&,bool);
		void findOddMeshes(PolygonClass*,int*,int,int*,int&,bool);

///////////

        int checkPolygonPolygonIntersections(ContainerClass&,PolygonClass*, int,int);

///////////

		void assignPolygonsToBoxes(ContainerClass&,PolygonClass*,SpaceClass&);
		void identifyBoxes(ContainerClass&,PolygonClass*,SpaceClass&,int,int*&,int&,int&);

///////////

		void setDeviationDistance(ContainerClass&,PolygonClass*, SpaceClass&);
		void findClosestVertex(ContainerClass&,PolygonClass*,SpaceClass&,int,int);
		void computeClosest(ContainerClass&,PolygonClass*,SpaceClass&,int,double*,double*);
		void computeVertexEnergyForce(ContainerClass&); 
		void computeNewVertexCoords(ContainerClass&,double); 
		void computePolygonIntersectionForce(ContainerClass&,PolygonClass*,SpaceClass&); 
		void setPolygonIntersection(ContainerClass&,PolygonClass*,SpaceClass&); 
};

//############################################################################
//############################################################################
//############################################################################

