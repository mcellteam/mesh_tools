//####################################################
//#################### parameters ####################
//####################################################
// default filename size
#define FILENAME_SIZE 1024

// set to true to write initialization information
//		for each file to stdout, i.e. verbose setting
// set to false for concise initialization informatin, i.e. concise setting
bool WRITE_VERBOSE_INIT = false;

// set to true to write punished vertices to
//		file specified by PUNISHED_FILE
bool WRITE_PUNISHED_VERTICES_TO_FILE = true;

// set to true to write intersected faces to
//		file specified by INTERSECTED_FILE
bool WRITE_INTERSECTED_FACES_TO_FILE = true;

// choice of output format for intersected faces
//	cp = dreamm custom points format
//	detail = face name and index + vertex details(3X)
// choose one
//#define FORMAT_INTERSECTED_FACES	cp
#define FORMAT_INTERSECTED_FACES	detail

// set to true to write nonnice vertices to
//		file specified by NONNICE_FILE
bool WRITE_NONNICE_VERTICES_TO_FILE = true;

// choice of output format for nonnice vertices
//	cp = dreamm custom points format
//	detail = vertex name,index,x,y,z
// choose one
//#define FORMAT_NONNICE_VERTICES		cp
#define FORMAT_NONNICE_VERTICES	detail

// set to true to write neighborhood to 
//		file specified by NEIGHBORHOOD_FILE
bool WRITE_NEIGHBORHOOD_TO_FILE = false;

// set to true to write closest point distance
//		to .dat files after every iteration
// set to false to only write distances to .dat
//		file after the last iteration has finished
bool WRITE_DISTANCES_EVERY_ITERATION = true;

// set to true to write meshes to .mesh files after every iteration
// set to false to only write meshes to .mesh
//		file after the last iteration has finished
bool WRITE_MESH_EVERY_ITERATION = true;

// set to true to use iteration number in mesh filename
// set to false to use a generic mesh filename and write over existing mesh files
// e.g. if true, then use 'filename_OUTPUT_SUFFIX_1.mesh'
//		else if false, then use 'filename_OUTPUT_SUFFIX.mesh'
// see below for OUTPUT_SUFFIX details
bool APPEND_ITERATION_NUMBER_TO_MESH = false;

// set to true to use iteration number in closest point distances filename
// set to false to use a generic dat filename and write over existing dat file
// e.g. if true, then use 'filename_1.dat'
//		else if false, then use 'filename.dat'
bool APPEND_ITERATION_NUMBER_TO_DISTANCES = true;

// set to true to prevent any new intersections of faces
// set to false to allow new intersecting faces
//	(useful when faces exist with all three vertices nonnice)
bool STRICT_FACE_INTERSECTION_PREVENTION = true;

// set print flag to '1' to enable diagnostic printing
// set print flag to '0' to disable print statements
bool PRINT_FLAG = false;

// set to true to check for complete separation of all faces
// i.e. if Container::ti==0, then set STRICT_FACE_INTERSECTION_PREVENTION=true;
// set to false to do nothing
bool DETECT_COMPLETE_SEPARATION = true;

// set to true to move nonnice vertices
// and intersected face vertices first
// regardless of their virutal displacement ranking
// set to false to strictly follow virtual displacement
// ranking , i.e. topN map
bool MOVE_NONNICE_AND_INTERSECTED_FIRST = false;

// defines the ratio between nanometers and the units of the input data
// 		e.g. in input is in microns, then set SCALE= 1.0/1000.0
#define SCALE	1.0/1.0	// nm/inpu

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
//#define NUM_ADJACENT_BOXES  1
#define NUM_ADJACENT_BOXES  2

#define BIG_BOX 160 // NUM_ADJACENT_BOXES*SPACE_LENGTH nm

// If a candidate closest point is more than SEARCH_RADIUS 
// distance away from given vertex, then the candidate closest 
// point is not chosen as closest point to given vertex.
// PICK SEARCH_RADIUS < NEIGHBORHOOD_RADIUS else closest point
// will be be just outside neighborhood on same object if no
// other point is found.
//#define SEARCH_RADIUS  70.0 // nm
#define SEARCH_RADIUS_SQ  16900.0 // 130 nm * 130 nm

// For a given vertex, faces containing vertices within
// the neighborhood of the given vertex are eliminated from
// the search pool for closest point for the given vertex.
// The neighborhood is the collection of all vertices
// on the same object as the given vertex that lie within
// N edges of the given vertex, where N is a function of
// NEIGHBORHOOD_RADIUS and mean edge length for that object.
//#define NEIGHBORHOOD_RADIUS		80 // nm
#define NEIGHBORHOOD_RADIUS_SQ		16900.0 // 130 nm * 130 nm

//  equivalent to a 45 degree sweep around normal
//  i.e. a 90 degree window around normal
#define CLOSEST_POINT_COSINE	0.70711

// if USE_EDGE_REFERENCE_LENGTH == true
// then compare edge lengths to EDGE_REFERENCE_LENGTH
// else use original length of each edge
bool USE_EDGE_REFERENCE_LENGTH = true;
#define EDGE_REFERENCE_LENGTH	35 // nm

// Force Function Weights
// Weights are normalized during the force calculation,
// so all that matters is the ratio of the four weights.
// ANGLE_STRETCH_WEIGHT is 5X to 10X weaker than the other forces

// original weights
//#define INTERSECTION_WEIGHT 	 0.1  // nN
//#define EDGE_STRETCH_WEIGHT 	29.9  // nN/nm
//#define SEPARATION_WEIGHT		60.0  // nN/nm
//#define ANGLE_STRETCH_WEIGHT 	10.0  // nN
#define INTERSECTION_WEIGHT 	 0.1  // nN
#define EDGE_STRETCH_WEIGHT 	14.9  // nN/nm
#define SEPARATION_WEIGHT	70.0  // nN/nm
#define ANGLE_STRETCH_WEIGHT 	15.0  // nN
//#define INTERSECTION_WEIGHT 	97.0  // nN
//#define EDGE_STRETCH_WEIGHT 	 1.0  // nN/nm
//#define SEPARATION_WEIGHT		 1.0  // nN/nm
//#define ANGLE_STRETCH_WEIGHT 	 1.0  // nN
//#define INTERSECTION_WEIGHT 	 0.1  // nN
//#define EDGE_STRETCH_WEIGHT 	 0.1  // nN/nm
//#define SEPARATION_WEIGHT		99.7  // nN/nm
//#define ANGLE_STRETCH_WEIGHT     0.1  // nN
//#define INTERSECTION_WEIGHT 	 0.1  // nN
//#define EDGE_STRETCH_WEIGHT 	 0.9  // nN/nm
//#define SEPARATION_WEIGHT		 1.0  // nN/nm
//#define ANGLE_STRETCH_WEIGHT    98.0  // nN
// NOTE THESE MUST SUM TO 100!!!!!!!!!!!!!!!

// The displacement of each vertex is equal to
// the product of the calculated force using the relative weights
// from above and the gain, where the gain is equal to TIME_STEP/DAMPING.
#define TIME_STEP				0.3 // s
//#define TIME_STEP				1.0 // s
#define DAMPING					1.0	 // nN/(nm/s)

// After each iteration, the new gain is equal to
// the product of the old gain and GAIN_SCALE.
//#define GAIN_SCALE				0.995	// nN/(nm/s)
//#define GAIN_SCALE				1.0	// nN/(nm/s)

// Number of groups of GROUP_SIZE vertex moves to execute.
#define NUM_GROUPS	30
//#define NUM_GROUPS	1

// size of vertex movie refractory period
// i.e. the minimum number of steps between moves of the same vertex
#define REFRACTORY_PERIOD 1000
//#define REFRACTORY_PERIOD 10

#define MAX_TOUCHES 20
//#define MAX_TOUCHES 1

// energy epsilon used after each group to decide if max_gain is to be lowered
// if ENERGY_WINDOW averaged energy on two consecutive iterations before group
// ends differ by less than ENERGY_EPSILON, then lower max_gain
#define ENERGY_EPSILON	1E-7

// energy averaging window size
#define ENERGY_WINDOW 10000

// upper and lower threshold on allowable
// number ofvertex moves per REFRACTORY_PERIOD
//#define UPPER_MOVES_THRESHOLD 10
//#define LOWER_MOVES_THRESHOLD 2
#define UPPER_MOVES_THRESHOLD 20
#define LOWER_MOVES_THRESHOLD 10

// gain modifcation fraction
// gain change = fraction*overall gain
//#define GAIN_CHANGE_FRACTION 0.03
//#define GAIN_CHANGE_FRACTION 0.002

// number of vertex moves per group
#define GROUP_SIZE 1e6
//#define GROUP_SIZE 1e5
//#define GROUP_SIZE 1e4
//#define GROUP_SIZE 1e3

// minimum acceptable virtual displacment of vertex
// if virtual_displacement*virtual_displacement < MIN_DISPLACEMENT_SQ
//  then put vertex into refractory period
#define MIN_DISPLACEMENT_SQ 1E-2

// pi
#define PI 3.14159265358979

// desired distance between object surfaces
//#define TARGET_SEPARATION			0.5   // nm
//#define TARGET_SEPARATION			9.5   // nm
//#define TARGET_SEPARATION			19.5   // nm
#define TARGET_SEPARATION			30   // nm
#define LOOP_TARGET_SEPARATION	50.0   // nm
//
//#define TARGET_SEPARATION			50.5   // nm
//#define TARGET_SEPARATION			60.0   // nm
//#define LOOP_TARGET_SEPARATION		70.0   // nm
//
//#define TARGET_SEPARATION			80   // nm
//#define LOOP_TARGET_SEPARATION		90.0   // nm

// vertex single-move displacement shall be
// capped at MAX_ACTUAL_DISPL_SQ.
// i.e. if (disp>MAX_ACTUAL_DISPL_SQ)
//      then disp=MAX_ACTUAL_DISPL_SQ
//#define MAX_ACTUAL_DISPL_SQ		9	// nm*nm
#define MAX_ACTUAL_DISPL_SQ		2500	// nm*nm
//#define MAX_ACTUAL_DISPL_SQ		36	// nm*nm

// maximum allowed runtime in seconds
#define MAX_RUNTIME  561600  // 6.5 days

//####################################################
//############# files and directories ################
//####################################################

//#define INPUT_DATA_DIR	"/tmp/ROUND2/20nm/fixed_meshes/"
//#define OUTPUT_DATA_DIR	"/tmp/ROUND2/20nm/morphed_meshes/"
//#define INPUT_DATA_DIR	"/snl/data/mcell/justin/paper2_tortuosity/meshes/60nm/output/"
//#define OUTPUT_DATA_DIR	"/snl/data/mcell/justin/paper2_tortuosity/meshes/60nm/sepdis/"
//#define INPUT_DATA_DIR	"/snl/data/mcell/justin/paper2_tortuosity/meshes/20nm/edge_5_sep_90_angle_5_output/"
//#define OUTPUT_DATA_DIR	"/snl/data/mcell/justin/paper2_tortuosity/meshes/20nm/edge_20_sep_70_angle_20_output/"
//#define INPUT_DATA_DIR	"/tmp/paper2/80nm/input/"
//#define INPUT_DATA_DIR	"/tmp/paper3/20_equal_edges/edge_5_fixed_meshes/"
//#define OUTPUT_DATA_DIR	"/tmp/paper2/80nm/output/"
//#define OUTPUT_DATA_DIR	"/tmp/paper3/20_equal_edges/edge_5_sep_90_angle_5_equal_edges/"

#define OUTPUT_SUFFIX			"MORPHED"
#define MAIN_LOG_FILE			"main.log"
#define SPACE_LOG_FILE			"space.log"
#define CONT_LOG_FILE			"cont.log"
#define SEP_LOG_FILE			"sep.log"
#define OBJECT_LIST_FILE		"object_list.log"
#define CONTROL_FILE			"parameter_values.dat"
#define NEIGHBORHOOD_FILE		"neighborhood.dat"
#define VERTEX_SELECTION_FILE           "vertex_select_histo.dat"
#define PUNISHED_FILE			"punished_vertices.dat"
#define INTERSECTED_FILE		"intersected_faces.dat"
#define NONNICE_FILE			"nonnice_vertices.dat"

//####################################################
//####################################################
//####################################################
