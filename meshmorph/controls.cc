//####################################################
//#################### parameters ####################
//####################################################
// set to true to write initialization information for each file to stdout, i.e. verbose setting
// set to false for concise initialization informatin, i.e. concise setting
bool write_verbose_init = false;

// set to true to write neighborhood to file specified by NEIGHBORHOOD_FILE
bool write_neighborhood_to_file = false;

// set to true to write closest point distance to .dat files after every iteration
// set to false to only write distances to .dat file after the last iteration has finished
bool write_distances_every_iteration = true;

// set to true to write meshes to .mesh files after every iteration
// set to false to only write meshes to .mesh file after the last iteration has finished
bool write_mesh_every_iteration = true;

// set to true to use iteration number in mesh filename
// set to false to use a generic mesh filename and write over existing mesh files
// e.g. if true, then use 'filename_OUTPUT_SUFFIX_1.mesh'
//		else if false, then use 'filename_OUTPUT_SUFFIX.mesh'
// see below for OUTPUT_SUFFIX details
bool append_iteration_number_to_mesh = false;

// set to true to use iteration number in closest point distances filename
// set to false to use a generic dat filename and write over existing dat file
// e.g. if true, then use 'filename_1.dat'
//		else if false, then use 'filename.dat'
bool append_iteration_number_to_distances = true;

// set to true to increase accuracy of line/polygon intersections at expense of speed
// set to false to increase speed at expense of accuracy during line/polygon intersections
bool detect_polygon_edge_intersection = false;

// set print flag to '1' to enable diagnostic printing
// set print flag to '0' to disable print statements
bool print_flag = false;

// if moving a vertex creates an edge angle
// less than EDGE_ANGLE_THRESHOLD,
// then vertex is not moved
//#define EDGE_ANGLE_THRESHOLD     2e-4 //radians
#define EDGE_ANGLE_THRESHOLD     0.017456 //radians == 1 degree
//#define EDGE_ANGLE_THRESHOLD     0.17453 //radians == 10 degrees
//#define EDGE_ANGLE_THRESHOLD     0.52360 //radians == 30 degrees

// small value for tweaking vertex position to avoid
// identical vertex positions
#define VERTEX_EPSILON           .001 //nm

// small value for tweaking nice ray to avoid
// intersecting polygon edge or
// being parallel with polygon
#define RAY_EPSILON           .00001 //nm

// for use with "is float close to zero?" 
// in conditional statement
#define DOUBLE_EPSILON	  		1E-10
#define DOUBLE_EPSILON_COMP		0.99999999990	// 1-DOUBLE_EPSILON

// subdivide space into cubes with sides of SPACE_LENGTH in size
#define SPACE_LENGTH	80 // nm

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
#define SEARCH_RADIUS  70.0 // nm

// For a given vertex, faces containing vertices within
// the neighborhood of the given vertex are eliminated from
// the search pool for closest point for the given vertex.
// The neighborhood is the collection of all vertices
// on the same object as the given vertex that lie within
// N edges of the given vertex, where N is a function of
// NEIGHBORHOOD_RADIUS and mean edge length for that object.
//#define NEIGHBORHOOD_RADIUS		30 // nm
#define NEIGHBORHOOD_RADIUS		80 // nm

#define CLOSEST_POINT_COSINE	0.34202 //  equivalent to a 70 degree sweep around normal
										//  i.e. a 140 degree window around normal

// desired distance between object surfaces
#define TARGET_SEPARATION		50.0   // nm

// Force Function Weights
// Weights are normalized during the force calculation,
// so all that matters is the ratio of the four weights.
// ANGLE_STRETCH_WEIGHT is 5X to 10X weaker than the other forces

#define INTERSECTION_WEIGHT 	 0.1  // nN
#define EDGE_STRETCH_WEIGHT 	 5.0  // nN/nm
//#define SEPARATION_WEIGHT		65.0  // nN/nm
#define SEPARATION_WEIGHT		75.0  // nN/nm
//#define ANGLE_STRETCH_WEIGHT 	29.9  // nN
#define ANGLE_STRETCH_WEIGHT 	19.9  // nN
// NOTE THESE MUST SUM TO 100!!!!!!!!!!!!!!!

// The displacement of each vertex is equal to
// the product of the calculated force using the relative weights
// from above and the gain, where the gain is equal to TIME_STEP/DAMPING.
#define TIME_STEP				0.6 // s
#define DAMPING					1.0	 // nN/(nm/s)

// After each iteration, the new gain is equal to
// the product of the old gain and GAIN_SCALE.
#define GAIN_SCALE				0.99	// nN/(nm/s)

// Number of primp iterations to execute.
#define ITERATIONS	80

// size of vertex movie refractory period
// i.e. the minimum number of steps between moves of the same vertex
#define REFRACTORY_PERIOD 1000

// max number of vertex moves per iteration
#define MAX_VERTEX_MOVES 1e6

// minimum acceptable virtual displacment of vertex
// if virtual_displacement*virtual_displacement < MIN_DISPLACEMENT_SQ
//  then put vertex into refractory period
#define MIN_DISPLACEMENT_SQ 1E-2

// pi
#define PI 3.14159265358979

//####################################################
//############# files and directories ################
//####################################################

//#define INPUT_DATA_DIR	"/spin/cnl/mcell/justin/mesh_data/full_set/ROUND2/40nm/super_primped_meshes/"
//#define OUTPUT_DATA_DIR	"/spin/cnl/mcell/justin/mesh_data/full_set/LAST_LAP2/"
//#define INPUT_DATA_DIR    "/spin/cnl/mcell/justin/mesh_data/full_set/raw_meshes/"
//#define OUTPUT_DATA_DIR    "/spin/cnl/mcell/justin/mesh_data/full_set/raw_meshes/"
//#define INPUT_DATA_DIR	"/tmp/LAST_LAP/20nm/"
//#define INPUT_DATA_DIR	"/tmp/LAST_LAP/40nmE/"
//#define INPUT_DATA_DIR	"/tmp/LAST_LAP/primped_meshes/"

#define INPUT_DATA_DIR	"/tmp/ROUND2/50nm/super_primped_meshes/"
#define OUTPUT_DATA_DIR	"/tmp/ROUND2/50nm/super_duper_primped_meshes/"
//#define INPUT_DATA_DIR		"/tmp/ROUND2/primped_meshes/"
//#define OUTPUT_DATA_DIR		"/tmp/ROUND2/super_primped_meshes/"

//#define OUTPUT_DATA_DIR	"/tmp/LAST_LAP/40nmF/"
//#define OUTPUT_DATA_DIR	"/tmp/LAST_LAP/"
//#define INPUT_DATA_DIR		"/spin/cnl/mcell/justin/mesh_data/full_set/fixed_meshes/"
//#define OUTPUT_DATA_DIR		"/spin/cnl/mcell/justin/mesh_data/full_set/fixed_meshes/"
//#define OUTPUT_DATA_DIR	"/spin/cnl/mcell/justin/mesh_data/full_set/primped_meshes/primped_meshes/"
#define OUTPUT_SUFFIX	"SMOOTH"

#define MAIN_LOG_FILE		"main.log"
#define SPACE_LOG_FILE		"space.log"
#define CONT_LOG_FILE		"cont.log"
#define OBJECT_LIST_FILE	"object_list.log"
#define CONTROL_FILE		"parameter_values.dat"
#define NEIGHBORHOOD_FILE	"neighborhood.dat"
#define INTERSECTION_FILE	"intersection.dat"
#define VERTEX_SELECTION_FILE	"vertex_select_histo.dat"

//####################################################
//####################################################
//####################################################
