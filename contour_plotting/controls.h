// default filename size
#define FILENAME_SIZE 1024

// set to true to increase accuracy of line/polygon intersections at expense of speed
// set to false to increase speed at expense of accuracy during line/polygon intersections
bool detect_polygon_edge_intersection = false;

// small value for tweaking nice ray to avoid
// intersecting polygon edge or
// being parallel with polygon
#define RAY_EPSILON           .00001 //nm

// for use with "is float close to zero?" 
// in conditional statement
#define DOUBLE_EPSILON	  		1E-10
#define DOUBLE_EPSILON_COMP		0.99999999990	// 1-DOUBLE_EPSILON

// subdivide space into cubes with sides of SPACE_LENGTH in size
//#define SPACE_LENGTH	40 // nm
#define SPACE_LENGTH	0.040 // nm

// target mean number faces per space subpartition
//#define FACES_PER_BOX	0.3
#define FACES_PER_BOX	0.1

// When searching for closest surface point to a given vertex,
// the subdivided space box of given vertex will be searched,
// as well as, all boxes within NUK_ADJACENT_BOXES away from 
// given vertex box. If NUM_ADJACENT_BOXES is zero, then only
// the given vertex box will be searched. If NUM_ADJACENT_BOXES
// is one, then all adjacent boxes to given vertex box will also be searched.
//
// choose such that NUM_ADJACENT_BOXES >= SEARCH_RADIUS/SPACE_LENGTH
#define NUM_ADJACENT_BOXES  1

// If a candidate closest point is more than SEARCH_RADIUS 
// distance away from given vertex, then the candidate closest 
// point is not chosen as closest point to given vertex.
// PICK SEARCH_RADIUS < NEIGHBORHOOD_RADIUS else closest point
// will be be just outside neighborhood on same object if no
// other point is found.
//#define SEARCH_RADIUS  70.0 // nm
#define SEARCH_RADIUS  0.070 // nm

// For a given vertex, faces containing vertices within
// the neighborhood of the given vertex are eliminated from
// the search pool for closest point for the given vertex.
// The neighborhood is the collection of all vertices
// on the same object as the given vertex that lie within
// N edges of the given vertex, where N is a function of
// NEIGHBORHOOD_RADIUS and mean edge length for that object.
//#define NEIGHBORHOOD_RADIUS		140 // nm
#define NEIGHBORHOOD_RADIUS		0.140 // nm

#define CLOSEST_POINT_COSINE	0.34202 //  equivalent to a 70 degree sweep around normal

// pi
#define PI 3.14159265358979

// number of bins in histograms
#define NUM_BINS 16
