
void Quicksort( int d[], int left, int right);
int Partition( int d[], int left, int right);
void Quicksort_double( double d[], int left, int right);
int Partition_double( double d[], int left, int right);

// ######################################
// ######################################

// definition of box element
class BoxClass {
	public:
		int* polygons;
		int num_polygons;
		int max_polygons;
		double limits[6]; //-x, +x, -y, +y, -z, +z
		void addPolygon(int&);
		~BoxClass(void);
};

// ######################################
// ######################################

// definition of polygon structure
class PolygonClass {
	public:
		int object_index;
		int vertices[3];
		double normal_components[3];
		int* boxes;
		int num_boxes;
		int max_boxes;
		bool intersection;
		int* intersections;
		int num_intersections;
		int max_intersections;
		int checkEdgeEdgeIntersection(double*[3],int[3][2]
										,double*[3],int,int,bool);
		int checkPolygonEdgePolygonIntersection(double*[3],int[3][2],
												double[3],double*[3],int,int);
		void recordBoxes(int*,int);
		~PolygonClass(void);
};

// ######################################
// ######################################

// definition of niceness structure
class NiceClass {
	public:
		int nice;
};

// ######################################
// ######################################

// definition of vertex structure
class EdgeClass {
	public:
		bool valid; // signals that valid data has been loaded into class
		int edge_vertex_indices[2];
		int other_vertex_indices[2];
		double edge_vertex_coords[3];  // coordinates of edge_vertex_indices[0]
		double other_vertex_coords[3];  // coordinates of other_vertex_indices[0]
		int object_index;
		int polygon_indices[2];
		double original_angle_cosine;
		double new_angle_cosine;
		double energy;
		double force[2][3];
		int flip; // indicates direction of force relative to normal 
					// (same,+1 or opposite,-1)
		void initAngles(PolygonClass*);
		void computeAngleForceEnergy(PolygonClass*);
};

// ######################################
// ######################################

// definition of vertex structure
class VertexClass {
	public:
		ofstream Vertexfile;
		int hashval;
		int vertex_index;
		double original_coordinates[3];
		double new_coordinates[3];
		int* adjacent_polygons_index;
		int next_adjacent_polygon;
		int max_adjacent_polygon;
		int* adjacent_vertices_index;
		int next_adjacent_vertex;
		int max_adjacent_vertex;
		int* neighborhood_vertices_index[HASH_TABLE_SIZE];
		int next_neighbor_vertex[HASH_TABLE_SIZE];
		int max_neighborhood_vertices_index[HASH_TABLE_SIZE];
		double* original_edge_lengths;
		int next_original_edge;
		double* new_edge_lengths;
		int next_new_edge;
		int* edges_indices;
		int num_edges_indices;
		int max_edges_indices;
		NiceClass nice_data;
		double outward_normal[3];
		double closest[3]; // coordinates
		double deviation_distance;
		double deviation_distance_old;
		double position_error; // coordinates
		double intersection_force[3];
		double energy;
		double force[3];
		bool candidate;
		void computeMeanOutwardNormal(PolygonClass *, double [3]);
		void computeEnergyForce(double *, double *, double *,int*,int,EdgeClass *); 
		void setOdd(int);
		void computeNewCoords(double); 
		void fileInit (int); 
		void fileOutit (void); 
		void updateLog (int);
		void addNeighbor (int&); 
		void addSortUniqueNeighbor (int&); 
		~VertexClass(void);
};

// ######################################
// ######################################

class ObjectClass
{
	public:
		bool frozen;
		ofstream Objectfile;
		char name[SIZE_OF_OBJECT_NAME];
		int nonnice;
		double position_error;	
		double force;
		double energy;
		VertexClass *vertices;
		EdgeClass *edges;
		int object_index;
		int num_edges;		// number of edges found in object
		int num_vertices; // number of vertices found in object
		int num_polygons; // number of polygons found in object
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
		~ObjectClass(void);
};

// ######################################
// ######################################

class ContainerClass
{
	public:
		ofstream Containerfile;
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

// ######################################
// ######################################

class ScanClass
{
	public:
		ofstream Scanfile;
		int num_files;							// number of input files
		char files[NUMBER_OF_INPUT_FILES][32];	// array of input file names
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
		void buildMeshAfter(ContainerClass&,PolygonClass*);
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

