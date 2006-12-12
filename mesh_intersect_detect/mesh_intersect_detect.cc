#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <string>
#include <dirent.h>

using namespace std;

//####################################################
//#################### parameters ####################
//####################################################

bool optimize_define_info = false;

// set to true to increase accuracy of line/polygon intersections at expense of speed
// set to false to increase speed at expense of accuracy during line/polygon intersections
bool detect_polygon_edge_intersection = false;

//bool print_flag = false;

// for use with "is float close to zero?" 
// in conditional statement
#define DOUBLE_EPSILON	  1E-10

// array sizes
#define NUMBER_OF_INPUT_FILES			600
#define SIZE_OF_DATA_FILE				128
#define SIZE_OF_OBJECT_NAME				32
#define NUMBER_OF_OBJECTS				600

#define NUMBER_OF_POLYGONS_IN_BOXES		7
#define INCREMENT_OF_POLYGONS_IN_BOXES	7

#define NUMBER_OF_BOXES					16
#define INCREMENT_OF_BOXES				8

#define NUMBER_OF_LIST_ENTRIES			8
#define INCREMENT_OF_LIST_ENTRIES		8

// subdivide space 
#define SPACE_LENGTH	40 // nm

//####################################################
//#################### definitions ###################
//####################################################

// definition of box element
class BoxClass {
	public:
		int* polygons;
		int num_polygons;
		int max_polygons;
		int x,y,z;	// indices of box in space
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
		int checkEdgeEdgeIntersection(double*[3],int[3][2],double*[3],int,int,bool);
		int checkPolygonEdgePolygonIntersection(double*[3],int[3][2],double[3],double*[3],int,int);
		void recordBoxes(int*,int);
		~PolygonClass(void);
};

// ######################################
// ######################################

// definition of vertex structure
class VertexClass {
	public:
		int vertex_index;
		double original_coordinates[3];
		double outward_normal[3];
};

// ######################################
// ######################################

class ObjectClass
{
	public:
		VertexClass *vertices;
		int object_index;
		int num_edges;		// number of edges found in object
		int num_vertices; // number of vertices found in object
		int num_polygons; // number of polygons found in object
		void InitVertices (void); 
		void makeObject(int ,int );
		void addOriginalVertexCoordinates(int ,double *);
		void getPolygonCoordinates(int [3],double* [3]);
		void boundObject(double*);
		~ObjectClass(void);
};

// ######################################
// ######################################

class ContainerClass
{
	public:
		int object_count;
		ObjectClass *objects;
		void init(void);
		void initObjects(int *,int *);
		~ContainerClass(void);
};

// ######################################
// ######################################

class ScanClass
{
	public:
		int num_files;							// number of input files
		int object_count;						// total number of objects in model
		int vertex_count;						// total number of vertices in model
		int polygon_count;						// total number of polygons in model
		int vertices_array[NUMBER_OF_OBJECTS];	// array of number of vertices per object
		int polygons_array[NUMBER_OF_OBJECTS];	// array of number of polygons per object
		void scanFiles(ContainerClass&,PolygonClass*, char*,bool);
		void scanFile(ContainerClass&,PolygonClass*, char* ,bool, int);
};

// ######################################
// ######################################

// definition of box element
class SpaceClass {
	public:
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
		~SpaceClass();
};

// ######################################
// ######################################

class ManipClass
{
	public:
		int object_count;				// total number of objects in model
		int vertex_count;				// total number of vertices in model
		int polygon_count;				// total number of polygons in model
		bool self_intersection;			// flag signifies if self-intersecting object was found
		void computePolygonNormals(ContainerClass&,PolygonClass*);
        int checkPolygonPolygonIntersections(ContainerClass&,PolygonClass*, int,int);
		void assignPolygonsToBoxes(ContainerClass&,PolygonClass*,SpaceClass&);
		void identifyBoxes(ContainerClass&,PolygonClass*,SpaceClass&,int,int*&,int&,int&);
		void computePolygonIntersectionForce(ContainerClass&,PolygonClass*,SpaceClass&); 
		void setPolygonIntersection(ContainerClass&,PolygonClass*,SpaceClass&); 
};

//####################################################
//################## subroutines  ####################
//####################################################

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

// #####################################################
// #####################################################

void biggest(double *x, int& big) {

    if ( x[0] < x[1] ) {
        // x[0] < x[1]
        if ( x[0] < x[2] ) {
            // x[0] < x[1]
            // x[0] < x[2]
            if (x[1] < x[2]) {
                // x[0] < x[1]
                // x[1] < x[2]
                // x[0] is smallest
                // x[2] is biggest
                big = 2;
            } else {
                // x[0] < x[1]
                // x[2] <= x[1]
                // x[0] is smallest
                // x[1] is biggest or tied with x[2]
                big = 1;
            }
        } else {
            // x[0] < x[1]
            // x[2] <= x[0]
            // x[2] is smallest or tied with x[0]
            // x[1] is biggest
            big = 1;
        }
    } else {
        // x[1] <= x[0]
        if ( x[0] < x[2] ) {
            // x[1] <= x[0]
            // x[0] < x[2]
            // x1 is smallest or tied with x[0]
            // x2 is biggest
            big = 2;
        } else {
        	big = 0;
        }
    }

}

// #####################################################
// #####################################################

void threeValueSort(double *x, double *biggest, double *smallest) {

    if ( x[0] < x[1] ) {
        // x[0] < x[1]
        if ( x[0] < x[2] ) {
            // x[0] < x[1]
            // x[0] < x[2]
            if (x[1] < x[2]) {
                // x[0] < x[1]
                // x[1] < x[2]
                // x[0] is smallest
                // x[2] is biggest
                *smallest = x[0];
                *biggest = x[2];
            } else {
                // x[0] < x[1]
                // x[2] <= x[1]
                // x[0] is smallest
                // x[1] is biggest or tied with x[2]
                *smallest = x[0];
                *biggest = x[1];
            }
        } else {
            // x[0] < x[1]
            // x[2] <= x[0]
            // x[2] is smallest or tied with x[0]
            // x[1] is biggest
            *smallest = x[2];
            *biggest = x[1];
        }
    } else {
        // x[1] <= x[0]
        if ( x[0] < x[2] ) {
            // x[1] <= x[0]
            // x[0] < x[2]
            // x1 is smallest or tied with x[0]
            // x2 is biggest
            *smallest = x[1];
            *biggest = x[2];
        } else {
            // x[1] <= x[0]
            // x[2] <= x[0]
            if (x[1] < x[2]) {
                // x[1] <= x[0]
                // x[2] <= x[0]
                // x[1] < x[2]
                // x[1] is smallest
                // x[0] is biggest
                *smallest = x[1];
                *biggest = x[0];
            } else {
                // x[1] <= x[0]
                // x[2] <= x[1]
                // x[2] is smallest or tied with x[1]
                // x[0] is biggest or tied with x[1]
                *smallest = x[2];
                *biggest = x[0];
            }
        }
    }

}

// #####################################################
// #####################################################

void checkLinePolygonIntersection(double* pvc[3],double pn[3],double lp[2][3],
									bool *line_flag, bool *poly_flag, bool *poly_edge_flag) {

	//pvc = polygon_vertex_coordinates
	//pn  = polygon_normal
	//lp  = line_points

	int i, j, a, b, tv[2];
    double num, den, pI[3], u , area[3], p, det[3];
	int big;

	*line_flag = false;
	*poly_flag = false;
	*poly_edge_flag = false;

	// do polygons intersect?
	// e.g. use plane-plane intersection and compute overlap
	// e.g. use edge-polygon intersection method (USED THIS HERE)

	// compute polygon area on each of three principal planes
	// xy
	area[0] = fabs(( (pvc[2][0]-pvc[1][0]) * pvc[0][1] +
					(pvc[0][0]-pvc[2][0]) * pvc[1][1] +
					(pvc[1][0]-pvc[0][0]) * pvc[2][1])/2.0);
	// yz
	area[1] = fabs(( (pvc[2][1]-pvc[1][1]) * pvc[0][2] +
					(pvc[0][1]-pvc[2][1]) * pvc[1][2] +
					(pvc[1][1]-pvc[0][1]) * pvc[2][2])/2.0);
	//zx
	area[2] = fabs(( (pvc[2][2]-pvc[1][2]) * pvc[0][0] +
					(pvc[0][2]-pvc[2][2]) * pvc[1][0] +
					(pvc[1][2]-pvc[0][2]) * pvc[2][0])/2.0);

	biggest(&area[0],big);
	if (big == 0) {
		tv[0] = 0;
		tv[1] = 1;
	} else if (big == 1) {
		tv[0] = 1;
		tv[1] = 2;
	} else {
		tv[0] = 0;
		tv[1] = 2;
	}


	// use line points, polygon normal and one polygon vertex
	// to compute where line intersects plane of polygon (if not parallel)

	// denominator of dot products
	den = pn[0]*(lp[1][0]-lp[0][0]) 
		+ pn[1]*(lp[1][1]-lp[0][1]) 
		+ pn[2]*(lp[1][2]-lp[0][2]);
	
	// if line and polygon plane are not parallel
	if  (den) {

		// numerator of dot products
		num = pn[0]*(pvc[0][0]-lp[0][0]) 
			+ pn[1]*(pvc[0][1]-lp[0][1]) 
			+ pn[2]*(pvc[0][2]-lp[0][2]);

		// point of intersection
		u = num/den;

		// if polygon cuts through line
		if (u > 0 && u < 1) {
			*line_flag = true;

			pI[0] = (1-u)*lp[0][0] + u*lp[1][0];
			pI[1] = (1-u)*lp[0][1] + u*lp[1][1];
			pI[2] = (1-u)*lp[0][2] + u*lp[1][2];

			////////// is point of intersection on other polygon? //////////
		
			// does point of intersection lie on polygon edge?
			*poly_edge_flag = false;

			if (detect_polygon_edge_intersection) {
				// if first coordinate of polygon vertex coordinates are distinguishable
				// else second coordinate of polygon vertex coordinates are distinguishable
				if (pvc[0][tv[0]] !=  pvc[1][tv[0]]) {
					p = pvc[0][tv[1]]+(pI[tv[0]]-pvc[0][tv[0]])/ 
						(pvc[1][tv[0]]-pvc[0][tv[0]])* 
						(pvc[1][tv[1]]-pvc[0][tv[1]]);
					if(!distinguishable(p,pI[tv[1]])){*poly_edge_flag = true;}
				} 
				else {
					p = pvc[0][tv[0]]+ (pI[tv[1]] -pvc[0][tv[1]])/ 
						(pvc[1][tv[1]] -pvc[0][tv[1]])* 
						(pvc[1][tv[0]] -pvc[0][tv[0]]);
					if(!distinguishable(p,pI[tv[0]])){*poly_edge_flag = true;}
				}
		
				if (pvc[1][tv[0]] != pvc[2][tv[0]]) {
					p = pvc[1][tv[1]]+ (pI[tv[0]] -pvc[1][tv[0]])/ 
						(pvc[2][tv[0]] -pvc[1][tv[0]])* 
						(pvc[2][tv[1]] -pvc[1][tv[1]]);
					if(!distinguishable(p,pI[tv[1]])){*poly_edge_flag = true;}
				} 
				else {
					p = pvc[1][tv[0]]+ (pI[tv[1]] -pvc[1][tv[1]])/ 
						(pvc[2][tv[1]] -pvc[1][tv[1]])* 
						(pvc[2][tv[0]] -pvc[1][tv[0]]);
					if(!distinguishable(p,pI[tv[0]])){*poly_edge_flag = true;}
				}
		
		
				if (pvc[2][tv[0]] != pvc[0][tv[0]]) {
					p = pvc[2][tv[1]]+ (pI[tv[0]] -pvc[2][tv[0]])/ 
						(pvc[0][tv[0]] -pvc[2][tv[0]])* 
						(pvc[0][tv[1]] -pvc[2][tv[1]]);
					if(!distinguishable(p,pI[tv[1]])){*poly_edge_flag = true;}
				} 
				else {
					p = pvc[2][tv[0]]+ (pI[tv[1]] -pvc[2][tv[1]])/ 
						(pvc[0][tv[1]] -pvc[2][tv[1]])* 
						(pvc[0][tv[0]] -pvc[2][tv[0]]);
					if(!distinguishable(p,pI[tv[0]])){*poly_edge_flag = true;}
				}
			}
		
			// if point of intersection is not on polygon edge
			if (!*poly_edge_flag) {

				// compute three determinants
				det[0] = (pvc[0][tv[0]] -pI[tv[0]]) *(pvc[1][tv[1]] -pI[tv[1]])
						-(pvc[1][tv[0]] -pI[tv[0]]) *(pvc[0][tv[1]] -pI[tv[1]]);
				det[1] = (pvc[1][tv[0]] -pI[tv[0]]) *(pvc[2][tv[1]] -pI[tv[1]])
						-(pvc[2][tv[0]] -pI[tv[0]]) *(pvc[1][tv[1]] -pI[tv[1]]);
	
				if ( ( (det[0] < 0) && (det[1] < 0) ) || ( (det[0] > 0) && (det[1] > 0) )){

					det[2] = (pvc[2][tv[0]] -pI[tv[0]]) *(pvc[0][tv[1]] -pI[tv[1]])
							-(pvc[0][tv[0]] -pI[tv[0]]) *(pvc[2][tv[1]] -pI[tv[1]]);
					if ( ( (det[0] < 0) && (det[1] < 0) && (det[2] < 0) ) || 
						( (det[0] > 0) && (det[1] > 0) && (det[2] > 0) ) ){
						// line intersects polygon plane inside polygon
					    *poly_flag = true;
					}
				}
			}
		}
	}
}

// #####################################################
// #####################################################

void addArrayElement(int*& ptr,int& num,int& max,int incr,int val) {

	int* temp;
	int i, temp_size;

	// before adding element to array
	// check if there is room
	if (num < max) {
		// add element
		ptr[num] = val;
		num++;
	} else {
		if (optimize_define_info) {
			cout << "addArrayElement: Consider making the initial value associated with " 
					<< incr << " bigger.\n";
			cout.flush();
		}
		// first need to make more room in array
		// create temp array
		temp = new int[max];
		if (temp == NULL) { 
			cout << "Not enough memory for temp in addArrayElement\n";
			cout.flush();
			exit(1);
		}
		// store current data in temp array
		temp_size = max;
		for (i=0;i<temp_size;i++) {
			temp[i] = ptr[i];
		}
		// delete polygons array
		delete[] ptr;
		// compute new, bigger array size
		max = max + incr;
		// create bigger array;
		ptr = new int[max];
		if (ptr == NULL) { 
			cout << "Not enough memory for temp in addArrayElement\n";
			cout.flush();
			exit(1);
		}
		// copy temp data back to array
		for (i=0;i<temp_size;i++) {
			ptr[i] = temp[i];
		}
		// delete temp array
		delete[] temp;
		// add new element
		ptr[num] = val;
		num++;
	}
}

// #####################################################
// #####################################################

void BoxClass::addPolygon(int& polygon_index) {

	int* temp;
	int i, temp_size;

	// before adding polygon to polygons array
	// check if there is room
	if (num_polygons < max_polygons) {
		// add polygon
		polygons[num_polygons] = polygon_index;
		num_polygons++;
	} else {
		// first need to make more room in polygons array
		// create temp array
		temp = new int[max_polygons];
		if (temp == NULL) { 
			cout << "Not enough memory for temp in addPolygon\n";
			cout.flush();
			exit(1);
		}
		// store current data in temp array
		temp_size = max_polygons;
		for (i=0;i<temp_size;i++) {
			temp[i] = polygons[i];
		}
		// delete polygons array
		delete[] polygons;
		// compute new, bigger polygons array size
		max_polygons = max_polygons+INCREMENT_OF_POLYGONS_IN_BOXES;
		// create bigger polygons array;
		polygons = new int[max_polygons];
		if (polygons == NULL) { 
			cout << "Not enough memory for polygon in addPolygons\n";
			cout.flush();
			exit(1);
		}
		// copy temp data back to polygons
		for (i=0;i<temp_size;i++) {
			polygons[i] = temp[i];
		}
		// delete temp array
		delete[] temp;
		// add new polygon
		polygons[num_polygons] = polygon_index;
		num_polygons++;
	}
}

// #####################################################
// #####################################################

BoxClass::~BoxClass(void) {

	delete[] polygons;

}

// #####################################################
// #####################################################

int PolygonClass::checkEdgeEdgeIntersection(	double* cpvc[3], int current_pairs[3][2], 
												double* opvc[3],  int current_polygon_index, 
												int other_polygon_index, bool share_edge_flag
											) {

	int i, j, k;
	double cv1[3], cv2[3], ov1[3], ov2[3], x[3], y[3];
	double num, den, u, q, qNum, qDen, uNum, uDen, cos_theta_square;
	double term1;
	bool intersect_flag, parallel_flag;

	// cpvc = current_polygon_vertex_coordinates
	// opvc = other_polygon_vertex_coordinates
	// cv   = current_vertex
	// ov   = other_vertex

	// initialize flag
	intersect_flag = false;

	// for each current polygon edge
	for (i=0;i<3;i++) {
		// for each other polygon edge
		for (j=0;j<3;j++) {
	
			// define two current and two other polygon vertices
			for (k=0;k<3;k++) {
				cv1[k] = cpvc[current_pairs[i][0]][k] ;
				cv2[k] = cpvc[current_pairs[i][1]][k] ;
				ov1[k] = opvc[current_pairs[j][0]][k] ;
				ov2[k] = opvc[current_pairs[j][1]][k] ;
			}
	
			// if the edges do not share a vertex
			if (
				(distinguishable(cv1[0],ov1[0]) &&
				 distinguishable(cv1[1],ov1[1]) &&
				 distinguishable(cv1[2],ov1[2]) ) 
				&&
				(distinguishable(cv1[0],ov2[0]) &&
				 distinguishable(cv1[1],ov2[1]) &&
				 distinguishable(cv1[2],ov2[2]) ) 
				&&
				(distinguishable(cv2[0],ov1[0]) &&
				 distinguishable(cv2[1],ov1[1]) &&
				 distinguishable(cv2[2],ov1[2]) ) 
				&&
				(distinguishable(cv2[0],ov2[0]) &&
				 distinguishable(cv2[1],ov2[1]) &&
				 distinguishable(cv2[2],ov2[2]) ) 
				) {

				// and the edges are not parallel
				parallel_flag = false;
				for (k=0;k<3;k++) {
					y[k] = ov2[k]-ov1[k];
				}

				term1 = (x[0]*y[0]+ x[1]*y[1]+ x[2]*y[2]);
				if ( !distinguishable(term1*term1,
				(x[0]*x[0]+ x[1]*x[1]+ x[2]*x[2])*(y[0]*y[0]+ y[1]*y[1]+ y[2]*y[2])) ) {parallel_flag = true;}
	
				if (!parallel_flag) {
					// compute scalars
					qDen = (cv2[0]-cv1[0])*(ov2[1]-ov1[1])-(cv2[1]-cv1[1])*(ov2[0]-ov1[0]);
					qNum = (cv2[0]-cv1[0])*(cv1[1]-ov1[1])-(cv2[1]-cv1[1])*(cv1[0]-ov1[0]);
					uNum = ov1[0]-cv1[0]+q*(ov2[0]-ov1[0]);
					uDen = cv2[0]-cv1[0];
	
					if(fabs(qDen)>DOUBLE_EPSILON && fabs(uDen)>DOUBLE_EPSILON) {
						q = qNum/qDen;
						u = uNum/uDen;
	
						if ( 
							((u > DOUBLE_EPSILON && u < (1-DOUBLE_EPSILON)) && 
							(q > DOUBLE_EPSILON && q < (1-DOUBLE_EPSILON)) ) ||
							( share_edge_flag && ((u > DOUBLE_EPSILON && u < (1-DOUBLE_EPSILON)) || 
												(q > DOUBLE_EPSILON && q < (1-DOUBLE_EPSILON))) )
							) {

						    return(1);
						}
					}
				}
			}
		}
	}

	return(0);

}

// #####################################################
// #####################################################

int PolygonClass::checkPolygonEdgePolygonIntersection(	double* cpvc[3],
														int cp[3][2], double on[3],
														double* opvc[3],  
														int current_polygon_index,
														int other_polygon_index) {

	// NOTE NEED ONLY CHECK EITHER (1) INTERSECTION OF CURRENT POLYGON EDGES
	// WITH OTHER POLYGON OR (2) VICE VERSA, SINCE THE ROLES OF CURRENT AND OTHER
	// POLYGON WILL GET REVERSED. CHOOSE OPTION 1.
	//
	// FURTHER NOTE THAT POLYGON INTERSECTION CHECK IS ASSYMETRIC
	// I.E. CHECK IS IF CURRENT POLYGON EDGE INTERSECTS OTHER POLYGON ONLY
	// NO CHECKING IS DONE IF CURRENT POLYGON IS INTERSECTED BY OTHER POLYGON EDGE

	//cpvc = current_polygon_vertex_coordinates
	//opvc = other_polygon_vertex_coordinates
	//tv   = target_vertices
	//cv1  = current_vertex 1
	//cv2  = current_vertex 2
	//on   = other_normal
	//cp   = current_pairs

	int i;
	double lp[2][3];
	bool line_flag, poly_flag, poly_edge_flag;

	line_flag = false;
	poly_flag = false;
	// for each current polygon edge
	for (i=0;i<3;i++) {
	
		if(!line_flag && !poly_flag) {

			lp[0][0] = cpvc[cp[i][0]][0];
			lp[0][1] = cpvc[cp[i][0]][1];
			lp[0][2] = cpvc[cp[i][0]][2];
			lp[1][0] = cpvc[cp[i][1]][0];
			lp[1][1] = cpvc[cp[i][1]][1];
			lp[1][2] = cpvc[cp[i][1]][2];

			checkLinePolygonIntersection(opvc,on,lp,&line_flag,&poly_flag,&poly_edge_flag);

		}
	}
	// do polygons intersect?
	if (line_flag && poly_flag) {return(1);}
	else {return(0);}
}

// #####################################################
// #####################################################

void PolygonClass::recordBoxes(int* list,int num_list) {

	int i;

	// record boxes in polygon struct
	if (boxes != NULL) {
		delete [] boxes;
	}
	boxes = new int[NUMBER_OF_BOXES];
	if (boxes == NULL) { 
		cout << "Not enough memory for boxes in PolygonClass::recordBoxes\n";
		cout.flush();
		exit(1);
	}
	num_boxes = 0;
	max_boxes = NUMBER_OF_BOXES;
	for (i=0;i<num_list;i++) {
		addArrayElement(boxes,num_boxes,max_boxes,INCREMENT_OF_BOXES,list[i]);
	}

}

// #####################################################
// #####################################################

PolygonClass::~PolygonClass(void) {

	delete[] boxes;
	delete[] intersections;

}

// #####################################################
// #####################################################

void ObjectClass::InitVertices(void) {

	int i;

	// initialize vertex indices
	for (i=0;i<num_vertices;i++) {
		vertices[i].vertex_index = i;
	}

}

// ######################################
// ######################################

void ObjectClass::makeObject(int number_of_vertices, int number_of_polygons) {

    num_vertices = number_of_vertices;
    num_polygons = number_of_polygons;
    vertices = new VertexClass[number_of_vertices];
	if (vertices == NULL) { 
		cout << "Not enough memory for VertexClass\n";
		cout.flush();
		exit(1);
	}
	
}

// ######################################
// ######################################

void ObjectClass::addOriginalVertexCoordinates(int current_vertex,double *coordinates) {

    if (current_vertex >= num_vertices) {

        cout << "Error: Number of vertices exceeded in object# " << object_index
            << ". Number of vertices = " << num_vertices << "\n";
        cout << "Offending vertex# = " << current_vertex << ", coordinates = "
            << coordinates[0] << " " << coordinates[1] << " " << coordinates[2] << "\n";

    } else {
        vertices[current_vertex].original_coordinates[0] = coordinates[0];
        vertices[current_vertex].original_coordinates[1] = coordinates[1];
        vertices[current_vertex].original_coordinates[2] = coordinates[2];
    }

}

// #####################################################
// #####################################################

void ObjectClass::getPolygonCoordinates(int indices[3],double* vertex[3]) {

	//get coordinates of each vertex index
	vertex[0] = &vertices[indices[0]].original_coordinates[0];
	vertex[1] = &vertices[indices[1]].original_coordinates[0];
	vertex[2] = &vertices[indices[2]].original_coordinates[0];

}

// ######################################
// ######################################

void ObjectClass::boundObject(double* range) {

    int j;
	double xmin,xmax,ymin,ymax,zmin,zmax;
	double x;
	double y;
	double z;

	//initialize mins and maxes
	xmin = vertices[0].original_coordinates[0];
	xmax = vertices[0].original_coordinates[0];
	ymin = vertices[0].original_coordinates[1];
	ymax = vertices[0].original_coordinates[1];
	zmin = vertices[0].original_coordinates[2];
	zmax = vertices[0].original_coordinates[2];


	// for each vertex in object
	for (j=0;j < num_vertices;j++) {

		///////// extract coordinates //////////
		x = vertices[j].original_coordinates[0];
		y = vertices[j].original_coordinates[1];
		z = vertices[j].original_coordinates[2];

		if (x>xmax) {xmax = x;}
		else if (x<xmin) {xmin = x;}
		if (y>ymax) {ymax = y;}
		else if (y<ymin) {ymin = y;}
		if (z>zmax) {zmax = z;}
		else if (z<zmin) {zmin = z;}
	}

	range[0] = xmin;
	range[1] = xmax;
	range[2] = ymin;
	range[3] = ymax;
	range[4] = zmin;
	range[5] = zmax;

}


// #####################################################
// #####################################################

ObjectClass::~ObjectClass(void) {

	delete[] vertices;

}

// #####################################################
// #####################################################

void ContainerClass::init(void) {

	// allocate memory for a specific number (full->object_count)
	// of array elements of type ObjectClass
	objects = new ObjectClass[object_count];
	if (objects == NULL) { 
		cout << "Not enough memory for ObjectClass\n";
		cout.flush();
		exit(1);
	}

}

// #####################################################
// #####################################################

void ContainerClass::initObjects(int *vertices_array, int *polygons_array) {

	int k;

	// for each object in objects array allocate memory 
	// for an array called vertices of size scan.vertices[k]
	// with each array element being a class of type vertex.
	for (k=0;k<object_count;k++) {
		objects[k].makeObject(vertices_array[k],polygons_array[k]);
		objects[k].object_index = k;
	}


}

// #####################################################
// #####################################################

ContainerClass::~ContainerClass(void) {

	delete[] objects;

}

// #####################################################
// #####################################################

void ScanClass::scanFiles (ContainerClass& container,PolygonClass *polygons,char *infile,bool scan_flag) {

	int count;

	// initialize cumulative counts
	vertex_count = 0;
	polygon_count = 0;
	num_files = 1;
	count = 0;
	object_count = 1;

	scanFile(container,polygons,infile,scan_flag,count);

}


// #####################################################
// #####################################################

void ScanClass::scanFile (ContainerClass& container,PolygonClass *polygons, 
							char* filename, bool scan_flag, int count) {

	char line[128];
    int vertex_num, polygon_num;
	int indices[3];
	double coor[3];
    std::string str,coor_str1,coor_str2,coor_str3;
    std::string::size_type pos1,pos2,pos3,pos4,pos5,len;

    std::ifstream inFile(filename);

	polygon_num = 0;
    vertex_num = 0;

    // open input data file
    if (inFile.fail()) // if stream cannot be opened
    {
        cout << "Can't open " << filename ; // display error msg and
        exit(1); // terminate abnormally
    }

	// foreach line in file
    while (inFile.getline(line,1024)) {
        str = line;

        // if vertex
        if (!str.find("Vertex",0)) {

			// count vertex
			vertex_num++;
			vertex_count++;

			if (!scan_flag) {

				// split string into Vertex, index and three coordinates
				len = str.length();
	            pos1 = str.find(" ",0);
	            pos2 = str.find(" ",pos1+1);
	            pos3 = str.find(" ",pos2+1);
	            pos4 = str.find(" ",pos3+1);
	            pos5 = str.find(" ",pos4+1);
				// check if extra space between index and first coordinate
				if (pos3 == pos2+1) {
					coor_str1 = str.substr(pos3+1,pos4-pos3-1);
					coor_str2 = str.substr(pos4+1,pos5-pos4-1);
					coor_str3 = str.substr(pos5+1,len-pos5-1);
				} else {
					coor_str1 = str.substr(pos2+1,pos3-pos2-1);
					coor_str2 = str.substr(pos3+1,pos4-pos3-1);
					coor_str3 = str.substr(pos4+1,len-pos4-1);
				}

				// add vertices to object
				coor[0] = atof(coor_str1.c_str());
				coor[1] = atof(coor_str2.c_str());
				coor[2] = atof(coor_str3.c_str());
           		container.objects[count].addOriginalVertexCoordinates(vertex_num-1,coor);
	
			}
		} else if (!str.find("Face",0)) {

			polygon_num++;
			polygon_count++;

			if (!scan_flag) {

				// split string into Face, index and three indices
				len = str.length();
	            pos1 = str.find(" ",0);
	            pos2 = str.find(" ",pos1+1);
	            pos3 = str.find(" ",pos2+1);
	            pos4 = str.find(" ",pos3+1);
				pos5 = str.find(" ",pos4+1);
				// check if extra space between index and first vertex index
				if (pos3 == pos2+1) {
					coor_str1 = str.substr(pos3+1,pos4-pos3-1);
					coor_str2 = str.substr(pos4+1,pos5-pos4-1);
					coor_str3 = str.substr(pos5+1,len-pos5-1);
				} else {
					coor_str1 = str.substr(pos2+1,pos3-pos2-1);
					coor_str2 = str.substr(pos3+1,pos4-pos3-1);
					coor_str3 = str.substr(pos4+1,len-pos4-1);
				}
	
				// add vertex indices to polygon
				// subtract 1 from each index to account for 0-based array indexing in c++
				indices[0] = atoi(coor_str1.c_str());
				indices[1] = atoi(coor_str2.c_str());
				indices[2] = atoi(coor_str3.c_str());
				polygons[polygon_count-1].object_index = count;
				polygons[polygon_count-1].vertices[0] = indices[0]-1;
				polygons[polygon_count-1].vertices[1] = indices[1]-1;
				polygons[polygon_count-1].vertices[2] = indices[2]-1;

			}
		}
	}
	
	// close file
    inFile.close();

	if (scan_flag) {
		// store number of vertices and polygons found in object
		vertices_array[count] = vertex_num;
		polygons_array[count] = polygon_num;
	}
}

// #####################################################
// #####################################################

void SpaceClass::init(void) {

	int x, y, z, loc, j;

	// subdivide space
	num_space[0] = (int) ceil( (world[1]-world[0])/SPACE_LENGTH );
	num_space[1] = (int) ceil( (world[3]-world[2])/SPACE_LENGTH );
	num_space[2] = (int) ceil( (world[5]-world[4])/SPACE_LENGTH );
	num_boxes = num_space[0]*num_space[1]*num_space[2];

	// allocate memory for boxes array with each 
	// array element of type box (boxes)
	if (boxes != NULL) {
		delete [] boxes;
	}
	boxes = new BoxClass[num_boxes];
	if (boxes == NULL) { 
		cout << "Not enough memory for BoxClass\n";
		cout.flush();
		exit(1);
	}

	// store box limits in boxes class
	// for each box
	for (z =0;z<num_space[2];z++) {
		for (y =0;y<num_space[1];y++) {
			for (x =0;x<num_space[0];x++) {
		
				loc = z*num_space[0]*num_space[1]+y*num_space[0]+x;
				boxes[loc].x = x;
				boxes[loc].y = y;
				boxes[loc].z = z;
			}
		}
	}

	// initialize polygons arrays
	for (j=0;j<num_boxes;j++) {
		boxes[j].polygons = new int[NUMBER_OF_POLYGONS_IN_BOXES];
		if (boxes[j].polygons == NULL) { 
			cout << "Not enough memory for boxes[j].polygons in SpaceClass::init\n";
			cout.flush();
			exit(1);
		}
		boxes[j].num_polygons=0;
		boxes[j].max_polygons=NUMBER_OF_POLYGONS_IN_BOXES;
	}

}

// #####################################################
// #####################################################

SpaceClass::~SpaceClass(void) {

	delete[] boxes;
}

// #####################################################
// #####################################################

void SpaceClass::clear(void) {

	int j;
	// clear boxes polygon list
	for (j=0;j<num_boxes;j++) {
		delete[] boxes[j].polygons;
		boxes[j].num_polygons=0;
		boxes[j].max_polygons=0;
	}

}

// ######################################
// ######################################

void SpaceClass::boundWorld(ContainerClass& container, int object_count) {

    int i, j;
	double xmin,xmax,ymin,ymax,zmin,zmax;
	double x;
	double y;
	double z;
	double range[6];

	//initialize mins and maxes
	xmin = container.objects[0].vertices[0].original_coordinates[0];
	xmax = container.objects[0].vertices[0].original_coordinates[0];
	ymin = container.objects[0].vertices[0].original_coordinates[1];
	ymax = container.objects[0].vertices[0].original_coordinates[1];
	zmin = container.objects[0].vertices[0].original_coordinates[2];
	zmax = container.objects[0].vertices[0].original_coordinates[2];

    ////////// loop through all objects //////////
    // for each object
    for (i=0;i<object_count;i++) {

		// get range of object vertices
		container.objects[i].boundObject(&range[0]);

		if (range[1]>xmax) {xmax = range[1];}
		if (range[0]<xmin) {xmin = range[0];}
		if (range[3]>ymax) {ymax = range[3];}
		if (range[2]<ymin) {ymin = range[2];}
		if (range[5]>zmax) {zmax = range[5];}
		if (range[4]<zmin) {zmin = range[4];}

    }

	if (xmin < 0) {
		world[0] = xmin * 1.01;
	} else {
		world[0] = xmin * 0.99;
	}

	if (xmax < 0) {
		world[1] = xmax * 0.99;
	} else {
		world[1] = xmax * 1.01;
	}

	if (ymin < 0) {
		world[2] = ymin * 1.01;
	} else {
		world[2] = ymin * 0.99;
	}

	if (ymax < 0) {
		world[3] = ymax * 0.99;
	} else {
		world[3] = ymax * 1.01;
	}

	if (zmin < 0) {
		world[4] = zmin * 1.01;
	} else {
		world[4] = zmin * 0.99;
	}

	if (zmax < 0) {
		world[5] = zmax * 0.99;
	} else {
		world[5] = zmax * 1.01;
	}


/*	if (print_flag) {
		cout << "\nworld bounds = [" 
			<< world[0] << " "
			<< world[1] << " "
			<< world[2] << " "
			<< world[3] << " "
			<< world[4] << " "
			<< world[5] << "]\n";
	}
*/
}

// #####################################################
// #####################################################

void SpaceClass::computeBoxesToCheck(double* pvc[3], int* &boxes_to_check,int& num) {

	// pvc = polygon_vertex_coordinates
	
	int x,y,z, box_range[6], loc;
	double poly_limits[6], pvc_x[3], pvc_y[3], pvc_z[3];

	// identify polygon box limits
	pvc_x[0] = pvc[0][0];
	pvc_x[1] = pvc[1][0];
	pvc_x[2] = pvc[2][0];
	pvc_y[0] = pvc[0][1];
	pvc_y[1] = pvc[1][1];
	pvc_y[2] = pvc[2][1];
	pvc_z[0] = pvc[0][2];
	pvc_z[1] = pvc[1][2];
	pvc_z[2] = pvc[2][2];
	threeValueSort(&pvc_x[0],&poly_limits[1], &poly_limits[0]);
	threeValueSort(&pvc_y[0],&poly_limits[3], &poly_limits[2]);
	threeValueSort(&pvc_z[0],&poly_limits[5], &poly_limits[4]);

	// compute box index range that contains polygon box
	// note this range is zero lower-bounded (lowest range is zeroth box)
	// total range is 0..num_space[i]-1
	box_range[0] = (int) floor( (poly_limits[0]-world[0])/SPACE_LENGTH );   // -x
	box_range[1] = (int) floor( (poly_limits[1]-world[0])/SPACE_LENGTH );	//  x
	box_range[2] = (int) floor( (poly_limits[2]-world[2])/SPACE_LENGTH );	// -y
	box_range[3] = (int) floor( (poly_limits[3]-world[2])/SPACE_LENGTH );	//  y
	box_range[4] = (int) floor( (poly_limits[4]-world[4])/SPACE_LENGTH );   // -z
	box_range[5] = (int) floor( (poly_limits[5]-world[4])/SPACE_LENGTH );	//  z

	int count = 0;
	num = (box_range[1]-box_range[0]+1)*(box_range[3]-box_range[2]+1)*(box_range[5]-box_range[4]+1);
	boxes_to_check = new int[num];
	if (boxes_to_check == NULL) { 
		cout << "Not enough memory for boxes_to_check\n";
		cout.flush();
		exit(1);
	}

	for (z = box_range[4];z<box_range[5]+1;z++){
		for (y = box_range[2];y<box_range[3]+1;y++){
			for (x = box_range[0];x<box_range[1]+1;x++){
				loc = z*num_space[0]*num_space[1]+y*num_space[0]+x;
				boxes_to_check[count] = loc;
				count++;
			}
		}
	}
}

// #####################################################
// #####################################################

void SpaceClass::checkVertexOverlap(double pvc[3][3],int box_index,bool *vertex_overlap) {

	int k;

	// pvc = polygon_vertex_coordinates

	// check if the coordinates of any vertex is between 
	// box_limits for all three coordinate axes
	*vertex_overlap = false;

	for (k=0;k<3;k++) {
		if (
			(
			(pvc[k][0] > world[0]+boxes[box_index].x*SPACE_LENGTH) && 
			(pvc[k][0] < world[0]+(boxes[box_index].x+1)*SPACE_LENGTH)  
			) && 
			(
			(pvc[k][1] > world[2]+boxes[box_index].y*SPACE_LENGTH) && 
			(pvc[k][1] < world[2]+(boxes[box_index].y+1)*SPACE_LENGTH)  
			) && 
			(
			(pvc[k][2] > world[4]+boxes[box_index].z*SPACE_LENGTH) && 
			(pvc[k][2] < world[4]+(boxes[box_index].z+1)*SPACE_LENGTH)  
			)  
			) {
			*vertex_overlap = true;
		}
	}
}

// #####################################################
// #####################################################
void SpaceClass::checkPolygonPlaneCrossesBox( int box_index ,double pvc[3][3] ,double *poly_normal
											,bool *pos_dot ,bool *neg_dot) {
	// pvc = polygon_vertex_coordinates

	double dot, A, B, C, limits[6];
	int i, j, k;

	*pos_dot = false;
	*neg_dot = false;

	A = pvc[0][0]*poly_normal[0];
	B = pvc[0][1]*poly_normal[1]; 
	C = pvc[0][2]*poly_normal[2];

		limits[0] = world[0]+boxes[box_index].x*SPACE_LENGTH;
		limits[1] = world[0]+(boxes[box_index].x+1)*SPACE_LENGTH;
		limits[2] = world[2]+boxes[box_index].y*SPACE_LENGTH;
		limits[3] = world[2]+(boxes[box_index].y+1)*SPACE_LENGTH;
		limits[4] = world[4]+boxes[box_index].z*SPACE_LENGTH;
		limits[5] = world[4]+(boxes[box_index].z+1)*SPACE_LENGTH;

	for (i=0;i<2;i++) {
		for (j=3;j>1;j--) {
			for (k=5;k>3;k--) {

				// compute dot product of polygon normal and vertex vectors
				dot = limits[i]*poly_normal[0]-A+limits[j]*poly_normal[1]-B+limits[k]*poly_normal[2]-C;

				// if any two dot products have opposite sign
				// then polygon plane crosses box
				if      (dot > 0 ) { *pos_dot = true;}
				else if (dot < 0 ) { *neg_dot = true;}
			}
		}
	}

}
// #####################################################
// #####################################################

void SpaceClass::checkPolygonEdgeIntersectBoxFace( int box_index ,double pvc[3][3]
											,bool *poly_edge,bool *box_face) {

	// pvc = polygon_vertex_coordinates

	int box_normal[3], m, n, i, val;
	double num, den, u, box_vertex[6][3], pI[3], 
			poly_x3[3], poly_y3[3], poly_z3[3], 
			poly_x4[3], poly_y4[3], poly_z4[3],
			box_limit[6];

	*poly_edge = false;
	*box_face = false;

	poly_x4[0] = pvc[1][0];
	poly_x4[1] = pvc[2][0];
	poly_x4[2] = pvc[0][0];

	poly_y4[0] = pvc[1][1];
	poly_y4[1] = pvc[2][1];
	poly_y4[2] = pvc[0][1];

	poly_z4[0] = pvc[1][2];
	poly_z4[1] = pvc[2][2];
	poly_z4[2] = pvc[0][2];

	box_limit[0] = world[0]+boxes[box_index].x*SPACE_LENGTH;
	box_limit[1] = world[0]+(boxes[box_index].x+1)*SPACE_LENGTH;
	box_limit[2] = world[2]+boxes[box_index].y*SPACE_LENGTH;
	box_limit[3] = world[2]+(boxes[box_index].y+1)*SPACE_LENGTH;
	box_limit[4] = world[4]+boxes[box_index].z*SPACE_LENGTH;
	box_limit[5] = world[4]+(boxes[box_index].z+1)*SPACE_LENGTH;

	box_vertex[0][0] = box_limit[0]; 
	box_vertex[0][1] = box_limit[2]; 
	box_vertex[0][2] = box_limit[4]; 

	box_vertex[1][0] = box_limit[1];
	box_vertex[1][1] = box_limit[3];
	box_vertex[1][2] = box_limit[5];

	box_vertex[2][0] = box_vertex[0][0]; 
	box_vertex[2][1] = box_vertex[0][1];
	box_vertex[2][2] = box_vertex[0][2];

	box_vertex[3][0] = box_vertex[1][0]; 
	box_vertex[3][1] = box_vertex[1][1];
	box_vertex[3][2] = box_vertex[1][2];

	box_vertex[4][0] = box_vertex[0][0]; 
	box_vertex[4][1] = box_vertex[0][1];
	box_vertex[4][2] = box_vertex[0][2];

	box_vertex[5][0] = box_vertex[1][0]; 
	box_vertex[5][1] = box_vertex[1][1];
	box_vertex[5][2] = box_vertex[1][2];

	//for each polygon edge
	for (m=0;m<3;m++) {

		// for each box face
		for (n=0;n<6;n++) {

			if (!*poly_edge  || !*box_face) {

			//use two polygon vertices, box face normal, 
			// and a box face vertex
			// to compute where line of polygon edge intersects
			// plane of box face (if not parallel)
			box_normal[0] = 0;
			box_normal[1] = 0;
			box_normal[2] = 0;
			val = -1;
			for (i=0;i<n;i++) {
				val = val*-1;
			}
			box_normal[n/2] = val;
	
			// dot products
			den = box_normal[0]*(poly_x4[m]-pvc[m][0])
				+ box_normal[1]*(poly_y4[m]-pvc[m][1])
				+ box_normal[2]*(poly_z4[m]-pvc[m][2]);

			// compute point of intersection
			if (den) {
				// point of intersection
				num = box_normal[0]*(box_vertex[n][0]-pvc[m][0])
					+ box_normal[1]*(box_vertex[n][1]-pvc[m][1])
					+ box_normal[2]*(box_vertex[n][2]-pvc[m][2]);
				u = num/den;
				pI[0] = pvc[m][0] + u*(poly_x4[m] -pvc[m][0]);
				pI[1] = pvc[m][1] + u*(poly_y4[m] -pvc[m][1]);
				pI[2] = pvc[m][2] + u*(poly_z4[m] -pvc[m][2]);
	
				// is intersection point between two polygon vertices 
				// and the two coordinate axes ranges?
				// if yes, then polygon is in box

				if ((u > 0) && (u < 1)) {	
					// point of intersection is between polygon vertices
					*poly_edge = true;
				}

				// is point of intersection between box face limits?
				if (n==0 || n==1) { // -x, +x
					if (
						(box_limit[2] < pI[1]) &&
						(box_limit[3] > pI[1]) &&
						(box_limit[4] < pI[2]) &&
						(box_limit[5] > pI[2]) 
						) { 
						*box_face = true;
					}
				} else if (n==2 || n==3) { // -y, y
					if (
						(box_limit[0] < pI[0]) &&
						(box_limit[1] > pI[0]) &&
						(box_limit[4] < pI[2]) &&
						(box_limit[5] > pI[2]) 
						) { 
						*box_face = true;
					}
				} else { // -z, z
					if (
						(box_limit[0] < pI[0]) &&
						(box_limit[1] > pI[0]) &&
						(box_limit[2] < pI[1]) &&
						(box_limit[3] > pI[1]) 
						) { 
						*box_face = true;
					}
				}
			}
			}		
		}	
	}
}

// #####################################################
// #####################################################

void SpaceClass::recordPolygon(int* list, int num_list,int polygon_index) {

	int* temp;
	int i, temp_size;

	// record polygon in boxes struct
	for (i=0;i<num_list;i++) {
		boxes[list[i]].addPolygon(polygon_index);
	}

}

// #####################################################
// #####################################################

void SpaceClass::checkBoxOverlap(double* pvc[3],int box_index,bool *box_overlap) {


    int z, y, x, k;
	double xvec[3], yvec[3], zvec[3];
	double biggestx, biggesty, biggestz, smallestx, smallesty, smallestz;
	double poly_limits[6], box_limits[6];

	xvec[0] = pvc[0][0];
	xvec[1] = pvc[1][0];
	xvec[2] = pvc[2][0];
	yvec[0] = pvc[0][1];
	yvec[1] = pvc[1][1];
	yvec[2] = pvc[2][1];
	zvec[0] = pvc[0][2];
	zvec[1] = pvc[1][2];
	zvec[2] = pvc[2][2];

    // sort object_index_list
	threeValueSort(&xvec[0], &biggestx, &smallestx);
	threeValueSort(&yvec[0], &biggesty, &smallesty);
	threeValueSort(&zvec[0], &biggestz, &smallestz);

    // identify polygon box limits
    poly_limits[0] = smallestx;
    poly_limits[1] = biggestx;
    poly_limits[2] = smallesty;
    poly_limits[3] = biggesty;
    poly_limits[4] = smallestz;
    poly_limits[5] = biggestz;

    // get box limits
	box_limits[0] = world[0]+boxes[box_index].x*SPACE_LENGTH;
	box_limits[1] = world[0]+(boxes[box_index].x+1)*SPACE_LENGTH;
	box_limits[2] = world[2]+boxes[box_index].y*SPACE_LENGTH;
	box_limits[3] = world[2]+(boxes[box_index].y+1)*SPACE_LENGTH;
	box_limits[4] = world[4]+boxes[box_index].z*SPACE_LENGTH;
	box_limits[5] = world[4]+(boxes[box_index].z+1)*SPACE_LENGTH;


    ////////// do boxes overlap //////////
    // yes if box max or min is between poly_limits for all three coordinate axes
    *box_overlap = false;

    if (
        (
        (box_limits[0] <= poly_limits[0]) && (box_limits[1] > poly_limits[0]) ||
        (box_limits[0] > poly_limits[0]) && (box_limits[0] < poly_limits[1])
        ) &&
        (
        (box_limits[2] <= poly_limits[2]) && (box_limits[3] > poly_limits[2]) ||
        (box_limits[2] > poly_limits[2]) && (box_limits[2] < poly_limits[3]) 
        ) &&
        (
        (box_limits[4] <= poly_limits[4]) && (box_limits[5] > poly_limits[4]) ||
        (box_limits[4] > poly_limits[4]) && (box_limits[4] < poly_limits[5]) 
        )
        ) {
        *box_overlap = true;
    }

}

// #####################################################
// #####################################################

void ManipClass::computePolygonNormals(ContainerClass& container ,PolygonClass *polygons) {

	int j, i, object_index, vertices[3];
	double uX, uY, uZ, vX, vY, vZ;
	double* vertex[3];

    // for each polygon structure in master polygons array
	for (j=0;j<polygon_count;j++) {

	    // get vertex_index of each vertex of polygon
		for (i=0;i<3;i++) {
			vertices[i] = polygons[j].vertices[i];
		}

		//get object index of polygon
		object_index = polygons[j].object_index;

		//get coordinates of each vertex index
		container.objects[object_index].getPolygonCoordinates(vertices,vertex);

		// compute vectors 01 and 12
		uX = vertex[1][0]-vertex[0][0];
		uY = vertex[1][1]-vertex[0][1];
		uZ = vertex[1][2]-vertex[0][2];
		vX = vertex[2][0]-vertex[0][0];
		vY = vertex[2][1]-vertex[0][1];
		vZ = vertex[2][2]-vertex[0][2];

		// add cross product (u x v) to normal in polygons
		polygons[j].normal_components[0] = uY*vZ-uZ*vY;
		polygons[j].normal_components[1] = uZ*vX-uX*vZ;
		polygons[j].normal_components[2] = uX*vY-uY*vX;

	}
}

// #####################################################
// #####################################################

int ManipClass::checkPolygonPolygonIntersections(ContainerClass &container ,PolygonClass *polygons,
													int cpi,int opi) {

	int i, j, k, object_index, other_object_index, opvi[3], cpvi[3], current_pairs[3][2], num_unique,
		single_shared_vert;
	double cn[3], on[3], X[3], Y[3], colinear, cos_theta_square;
	double* opvc[3];
	double* cpvc[3];
	double term1;
	bool coplanar_flag, colinear_flag, parallel_flag, share_edge_flag, identical_flag,share_vert_flag;

	// cpvi = current_polygon_vertex_indices
	// cpvc = current_polygon_vertex_coordinates
	// opvc = other_polygon_vertex_coordinates
	// opvi = other_polygon_vertex_indices
	// cpi = current_polygon_index
	// opi = other_polygon_index
	// cn   = current_normal
	// on   = other_normal

	// get current polygon object index
	object_index = polygons[cpi].object_index;

	// get current polygon normal
	cn[0] = polygons[cpi].normal_components[0];
	cn[1] = polygons[cpi].normal_components[1];
	cn[2] = polygons[cpi].normal_components[2];

	// get current polygon vertex indices
	cpvi[0] = polygons[cpi].vertices[0];
	cpvi[1] = polygons[cpi].vertices[1];
	cpvi[2] = polygons[cpi].vertices[2];

	// get current polygon vertex coordinates
	container.objects[object_index].getPolygonCoordinates(polygons[cpi].vertices, cpvc);

	// describe three pairs of current polygon vertices
    current_pairs[0][0] = 0;
    current_pairs[0][1] = 1;
    current_pairs[1][0] = 1;
    current_pairs[1][1] = 2;
    current_pairs[2][0] = 2;
    current_pairs[2][1] = 0;

	// get other_polygon object index
	other_object_index = polygons[opi].object_index;

	// get other_polygon normal
	on[0] = polygons[opi].normal_components[0];
	on[1] = polygons[opi].normal_components[1];
	on[2] = polygons[opi].normal_components[2];

	// get other_polygon vertex indices
	opvi[0] = polygons[opi].vertices[0] ;
	opvi[1] = polygons[opi].vertices[1] ;
	opvi[2] = polygons[opi].vertices[2] ;

	// get other_polygon vertex coordinates
	container.objects[other_object_index].getPolygonCoordinates(polygons[opi].vertices, opvc);

	// initialize flags
	coplanar_flag = false;
	colinear_flag = false;
	parallel_flag = false;
	share_edge_flag = false;
	identical_flag = false;
	share_vert_flag = false;

	// are current polygon and other polygon parallel
	// i.e. is angle between normals equal to zero?
	// i.e. is the square of the cosine of the angle equal to 1?
	term1 = cn[0]*on[0]+ cn[1]*on[1]+ cn[2]*on[2];

	if ( !distinguishable(term1*term1,
	(cn[0]*cn[0]+cn[1]*cn[1]+cn[2]*cn[2])*(on[0]*on[0]+on[1]*on[1]+on[2]*on[2])) ) {parallel_flag = true;}

	// dot product of other polygon normal and 
	// line connecting point on other polygon to point on current polygon
	// try to choose vertices not shared by both polygons

	if ( (opvi[0] != cpvi[0]) && (opvi[0] != cpvi[1]) && (opvi[0] != cpvi[2]) ) { 
		X[0] = opvc[0][0]; 
		X[1] = opvc[0][1]; 
		X[2] = opvc[0][2];
	} else if ( (opvi[1] != cpvi[0]) && (opvi[1] != cpvi[1]) && (opvi[2] != cpvi[2]) ) { 
		X[0] = opvc[1][0]; 
		X[1] = opvc[1][1]; 
		X[2] = opvc[1][2];
	} else { 
		X[0] = opvc[2][0]; 
		X[1] = opvc[2][1]; 
		X[2] = opvc[2][2];
	}

	if ( (cpvi[0] != opvi[0]) && (cpvi[0] != opvi[1]) && (cpvi[0] != opvi[2]) ) {
		Y[0] = cpvc[0][0]; 
		Y[1] = cpvc[0][1]; 
		Y[2] = cpvc[0][2]; 
	} else if ( (cpvi[1] != opvi[0]) && (cpvi[1] != opvi[1]) && (cpvi[1] != opvi[2]) ) {
		Y[0] = cpvc[1][0]; 
		Y[1] = cpvc[1][1]; 
		Y[2] = cpvc[1][2]; 
	} else { 
		Y[0] = cpvc[2][0]; 
		Y[1] = cpvc[2][1]; 
		Y[2] = cpvc[2][2]; 
	}

	colinear = on[0]*(X[0]-Y[0]) + on[1]*(X[1]-Y[1]) + on[2]*(X[2]-Y[2]);

	// if polygons colinear, then normal and line are orthogonal
	// and dot product will be ~zero

	if (fabs(colinear)<DOUBLE_EPSILON) {colinear_flag = true;}
	if (parallel_flag && colinear_flag) {coplanar_flag = true;}
		
	// how many vertices are shared between current and other polygon
	num_unique = 0;
	// for each current polygon vertex
	for (j=0;j<3;j++) {
		if (
			distinguishable(cpvc[j][0], opvc[0][0]) ||
			distinguishable(cpvc[j][1], opvc[0][1]) ||
			distinguishable(cpvc[j][2], opvc[0][2])
			) {num_unique++;} else {single_shared_vert = j;}
		if (
			distinguishable(cpvc[j][0], opvc[1][0]) ||
			distinguishable(cpvc[j][1], opvc[1][1]) ||
			distinguishable(cpvc[j][2], opvc[1][2])
			) {num_unique++;} else {single_shared_vert = j;}
		if (
			distinguishable(cpvc[j][0], opvc[2][0]) ||
			distinguishable(cpvc[j][1], opvc[2][1]) ||
			distinguishable(cpvc[j][2], opvc[2][2])
			) {num_unique++;} else {single_shared_vert = j;}
	}
	if (9-num_unique == 1) {share_vert_flag = true;}
	else if (9-num_unique == 2) {share_edge_flag = true;}
	else if (9-num_unique == 3) {identical_flag = true;}

	////////// begin decision tree //////////

	if (coplanar_flag) {
		if (identical_flag) {
			// polygons are identical
			return(1);
		} else {
			//do polygon edges intersect?	
			// if yes, intersect
			// if no, do not intersect
			if (polygons[cpi].checkEdgeEdgeIntersection(cpvc, current_pairs, opvc,
													cpi, opi, share_edge_flag)) {return(1);}
			else {return(0);}
		}
	} else {
		if (share_vert_flag) {
			// single vertex shared
			double lp[2][3];
			bool line_flag, poly_flag, poly_edge_flag;
			line_flag = false;
			poly_flag = false;
			int m,n;
			if (single_shared_vert == 0){m = 1;n=2;}
			else if (single_shared_vert == 1){m = 0;n=2;}
			else {m = 0;n=1;}
			lp[0][0] = cpvc[m][0];
			lp[0][1] = cpvc[m][1];
			lp[0][2] = cpvc[m][2];
			lp[1][0] = cpvc[n][0];
			lp[1][1] = cpvc[n][1];
			lp[1][2] = cpvc[n][2];
			checkLinePolygonIntersection(opvc,on,lp,&line_flag,&poly_flag,&poly_edge_flag);
			// do polygons intersect?
			if (line_flag && poly_flag) {return(1);}
			else {return(0);}
		} else if (!share_edge_flag) {
			// do polygons intersect?
			if (polygons[cpi].checkPolygonEdgePolygonIntersection( cpvc, current_pairs, on,  
																opvc, cpi, opi)) {return(1);}
			else {return(0);}
		} else { return(0);}
	}
}

// #####################################################
// #####################################################

void ManipClass::assignPolygonsToBoxes(ContainerClass& container ,PolygonClass *polygons, 
											SpaceClass& space) {

	int j;
	int* list;
	int num_list;
	int max_list;

	////////// identify in which boxes each polygon exists ////////

    // for each polygon structure in master polygons array
	for (j=0;j<polygon_count;j++) {

		list = new int[NUMBER_OF_LIST_ENTRIES];
		if (list == NULL) { 
			cout << "Not enough memory for list in ManipClass::assignPolygonsToBoxes\n";
			cout.flush();
			exit(1);
		}
		num_list = 0;
		max_list = NUMBER_OF_LIST_ENTRIES;

		// identify boxes containing part of polygon
		identifyBoxes(container,polygons,space,j,list,num_list,max_list);

		if (num_list == 0) {
			cout << "ERROR: NO BOXES FOR POLYGON INDEX " << j << "\n";
		}

		// record boxes in polygon struct
		polygons[j].recordBoxes(list,num_list);

		// record polygon in boxes struct
		space.recordPolygon(list,num_list,j);

		// free memory
		delete[] list;
	}
}

// #####################################################
// #####################################################

void ManipClass::identifyBoxes(ContainerClass& container ,PolygonClass *polygons,SpaceClass &space
								,int polygon_index,int*& list,int& num_list,int& max_list) {

	int i, num_boxes, j, object_index, polygon_vertex_index, k;
	double poly_normal[3];
	double* pvc[3];
	bool space_overlap, vertex_overlap, pos_dot, neg_dot, poly_edge, box_face, poly_face, box_edge;
	bool box_overlap;

	// pvc = polygon_vertex_coordinates

	// get polygon object index
	object_index = polygons[polygon_index].object_index;

	// get polygon normal
	for (j=0;j<3;j++) {
		poly_normal[j] = polygons[polygon_index].normal_components[j];
	}

	container.objects[object_index].getPolygonCoordinates(polygons[polygon_index].vertices,pvc);

	////////// compute boxes to check//////////
	int* boxes_to_check;
	space.computeBoxesToCheck(pvc,boxes_to_check,num_boxes);

	////////// for each box //////////
	for (k=0;k<num_boxes;k++) {
		i = boxes_to_check[k];
		space.checkBoxOverlap(pvc,i,&box_overlap);
		if (box_overlap) {
			addArrayElement(list,num_list,max_list,INCREMENT_OF_LIST_ENTRIES,i);
		}
	}
	delete[] boxes_to_check;
}

// #####################################################
// #####################################################

void ManipClass::computePolygonIntersectionForce(ContainerClass& container,
												PolygonClass *polygons, SpaceClass& space) {

	int i, v1, v2, v3, object_index;
	double nX, nY, nZ, normal_length, fX, fY, fZ;

	// set intersection bool for each polygon
	setPolygonIntersection(container,polygons,space);

	///// set force /////
	// for each polygon
    for (i=0;i<polygon_count;i++) {

		// if intersecting
		if (polygons[i].intersection) {

			// get vertex indices
			v1 = polygons[i].vertices[0];
			v2 = polygons[i].vertices[1];
			v3 = polygons[i].vertices[2];

			// get object index
			object_index = polygons[i].object_index;

			cout << "Intersected polygon\n"
				<< "Vertex " << v1 << " " 
				<< container.objects[object_index].vertices[v1].original_coordinates[0] << " "
				<< container.objects[object_index].vertices[v1].original_coordinates[1] << " "
				<< container.objects[object_index].vertices[v1].original_coordinates[2] << endl
				<< "Vertex " << v2 << " " 
				<< container.objects[object_index].vertices[v2].original_coordinates[0] << " "
				<< container.objects[object_index].vertices[v2].original_coordinates[1] << " "
				<< container.objects[object_index].vertices[v2].original_coordinates[2] << endl
				<< "Vertex " << v3 << " " 
				<< container.objects[object_index].vertices[v3].original_coordinates[0] << " "
				<< container.objects[object_index].vertices[v3].original_coordinates[1] << " "
				<< container.objects[object_index].vertices[v3].original_coordinates[2] << endl << endl;
		}
	}
}

// #####################################################
// #####################################################

void ManipClass::setPolygonIntersection(ContainerClass& container, PolygonClass *polygons, SpaceClass& space) {

	int i, j, k;
	int* combinations[2];
	int num_combinations;
	int* temp[2];
	int num_temp;
	int next_pair;

	// initialize polygon intersection status
    for (i=0;i<polygon_count;i++) {
		polygons[i].intersection = false;
	}

    // for each subspace
    for (i=0;i<space.num_boxes;i++) {

		// compute number of all possible polygon pairs
		num_combinations = 0;
        for (j=0;j<space.boxes[i].num_polygons;j++) {
			num_combinations += j;
		}

		// allocate memory for combination array
		combinations[0] = new int[num_combinations];
		combinations[1] = new int[num_combinations];

		// load combinations arrays
		num_combinations = 0;
        for (j=0;j<space.boxes[i].num_polygons-1;j++) {
        	for (k=j+1;k<space.boxes[i].num_polygons;k++) {
				combinations[0][num_combinations] = space.boxes[i].polygons[j];
				combinations[1][num_combinations] = space.boxes[i].polygons[k];
				num_combinations++;
			}
		}
		

        // while there are polygon combinations to check 
		next_pair = 0;
		while (next_pair<num_combinations) {

			////// check for intersection of next two polygons //////
			// if polygons intersect
			if (checkPolygonPolygonIntersections(container,polygons,
												combinations[0][next_pair],
												combinations[1][next_pair])) {

				// record intersection
				polygons[combinations[0][next_pair]].intersection = true;
				polygons[combinations[1][next_pair]].intersection = true;
				
				///// remove all instances of each polygon from combinations /////
				// allocate memory for temp array
				temp[0] = new int[num_combinations];
				temp[1] = new int[num_combinations];
				num_temp = 0;

				// for all remaining combinations
				for (k=next_pair+1;k<num_combinations;k++) {
					// if combination does not contain either of the intersection polygons
					if ((combinations[0][k] != combinations[0][next_pair]) &&
						(combinations[1][k] != combinations[0][next_pair]) &&
						(combinations[0][k] != combinations[1][next_pair]) &&
						(combinations[1][k] != combinations[1][next_pair]) ) {
						// add combination to temp
						temp[0][num_temp] = combinations[0][k];
						temp[1][num_temp] = combinations[1][k];
						num_temp++;
					}
				}

				// delete combinations array
				delete[] combinations[0];
				delete[] combinations[1];

				// allocate memory for combination array
				combinations[0] = new int[num_temp];
				combinations[1] = new int[num_temp];
				num_combinations = num_temp;
				next_pair = 0;

				// copy temp to combination
				for (k=0;k<num_temp;k++) {
					combinations[0][k] = temp[0][k];
					combinations[1][k] = temp[1][k];
				}

				// delete temp
				delete[] temp[0];
				delete[] temp[1];

			} else {
				// polygons do not intersect
				next_pair++;
			}

		}

		delete[] combinations[0];
		delete[] combinations[1];

	}
}


//####################################################
//###################### main  #######################
//####################################################

int main(int argc,char *argv[]){

	if (argc != 2)
	{
		printf("\nSyntax: mesh_intersect_detect input_file\n\n");
		printf("Description: Finds intersecting polygons in mesh file,\n");
		printf("		and reports results to stdout.\n");
		return 1;
	}

    fprintf(stderr,"\nInput mesh files is assumed to have vertex and\n");
    fprintf(stderr,"faces with sequentially increasing indices. Run\n");
    fprintf(stderr,"the mesh input file through mesh_renumber first if\n");
    fprintf(stderr,"this criterion is not met.\n\n");

	int i;
	char *infile,*eptr;
	infile = argv[1];

	// declare a pointer (objects) to an array of size zero
	// whose future elements will be of type ObjectClass
	ContainerClass container;

	// declare a pointer (polygons) to an array of size zero
	// whose future elements will be of type polygon
	PolygonClass *polygons; 

	////////// scan input file //////////
	// scan file and count number of objects and 
	// store number of vertices and polygons for each object
	ScanClass scan;
	scan.scanFiles(container,polygons,infile,1);

	////////// transfer scan data to container class //////////
	container.object_count = scan.object_count;

	////////// declare data structure //////////
	container.init();
	container.initObjects(scan.vertices_array,scan.polygons_array);

	// allocate memory for polygons array with each 
	// array element of type polygon (polygons)
	polygons = new PolygonClass[scan.polygon_count];
	if (polygons == NULL) { 
		cout << "Not enough memory for PolygonClass\n";
		cout.flush();
		exit(1);
	}
	for(i=0;i<scan.polygon_count;i++) {
		polygons[i].boxes = NULL;
		polygons[i].intersections = NULL;
	}

	////////// get data //////////
	// load objects with data from file
	scan.scanFiles(container,polygons,infile,0);

	////////// initializing space data structure //////////
	SpaceClass space;
	space.boundWorld(container,scan.object_count);
	space.boxes = NULL;
	space.init();

	////////// initializing manip class //////////
	ManipClass manip;

	////////// transfer scan data to other classes //////////
	manip.object_count = scan.object_count;
	manip.vertex_count = scan.vertex_count;
	manip.polygon_count = scan.polygon_count;

	manip.computePolygonNormals(container,polygons);

	////////// assign polygons to boxes //////////
	manip.assignPolygonsToBoxes(container,polygons,space);

	manip.computePolygonIntersectionForce(container,polygons,space);
//	manip.setPolygonIntersection(container,polygons,space);
}
