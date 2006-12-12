//
// File name : fea4.c
//
// Author : Justin Kinney
// Date : 12 Dec 2003
// 
// Purpose : FEA of spatial concentration gradients in MCELL model.
// 	     Diffusion across narrow cleft between a release surface
// 	     and a transport surface.
//
// Output : Concentration profiles every usec (each column is profile across cleft)
// 						first column is t=0
// 	    Total Glutamate in cleft every usec (first column is time
// 	    					in usec with 400usec offset.
// 	    					Second column is # GLU molecules)
// 	    					

//####################################################
//#################### parameters ####################
//####################################################

// set to true to write meshes to .mesh files after every iteration
// set to false to only write meshes to .mesh file after the last iteration has finished
bool write_mesh_every_iteration = true;

// set to true to increase accuracy of line/polygon intersections at expense of speed
// set to false to increase speed at expense of accuracy during line/polygon intersections
bool detect_polygon_edge_intersection = false;

// set to true to drive polygon angles to flat
// set to true to drive polygon angles to initial angle
bool init_angle_to_flat = true;

// switch to control printing of each vertex's stats
bool vertex_print = true;

// nonnice vertices are ones that are located inside a neighboring mesh
// the neighboring objects is considered to be violated
// set violated_flag to true to store violated objects for nonnice vertices
// set violated_flag to false to not store violated objects for nonnice vertices
//bool violated_flag = false;

// set print flag to '1' to enable diagnostic printing
// set print flag to '0' to disable print statements
bool print_flag = false;

// small value for tweaking nice ray to avoid
// intersecting polygon edge or
// being parallel with polygon
#define RAY_EPSILON           .001 //nm

// for use with "is float close to zero?" 
// in conditional statement
#define DOUBLE_EPSILON	  1E-10

// array sizes
#define NUMBER_OF_INPUT_FILES			600
#define SIZE_OF_DATA_FILE				128
#define SIZE_OF_OBJECT_NAME				32
#define NUMBER_OF_OBJECTS				600
#define NUMBER_OF_ADJACENT_POLYGONS		10
#define INCREMENT_OF_ADJACENT_POLYGONS	10
#define NUMBER_OF_EDGES_VERTICES		5
#define INCREMENT_OF_EDGES_VERTICES		5
#define NUMBER_OF_ADJACENT_VERTICES		5
#define INCREMENT_OF_ADJACENT_VERTICES	5
#define NUMBER_OF_POLYGONS_IN_BOXES		5
#define INCREMENT_OF_POLYGONS_IN_BOXES	5
#define NUMBER_OF_NEIGHBORS				100
#define INCREMENT_OF_NEIGHBORS			50
#define NUMBER_OF_CROSSED_POLYGONS		50
#define INCREMENT_OF_CROSSED_POLYGONS	25
#define NUMBER_OF_UNIQUE_POLYGONS		300
#define INCREMENT_OF_UNIQUE_POLYGONS	50
#define NUMBER_OF_BOXES					8
#define INCREMENT_OF_BOXES				8
#define NUMBER_OF_INTERSECTIONS			1
#define INCREMENT_OF_INTERSECTIONS		1
#define NUMBER_OF_LIST_ENTRIES			8
#define INCREMENT_OF_LIST_ENTRIES		8

// subdivide space 
#define SPACE_LENGTH	20 // nm

// closest vertex parameters
#define NUM_ADJACENT_BOXES  2
#define SEARCH_RADIUS  20 // nm

// energy function parameters
#define TARGET_SEPARATION		20   // nm
#define SEPARATION_WEIGHT		10 // nN/nm
#define EDGE_STRETCH_WEIGHT 	5  // nN/nm
#define ANGLE_STRETCH_WEIGHT 	10.0  // nN
#define INTERSECTION_WEIGHT 	10.0  // nN
#define TIME_STEP				0.005 // s
#define DAMPING					1	 // nN/(nm/s)
#define GAIN_SCALE				1.1 // nN/(nm/s)
//#define DISPLACEMENT_WEIGHT 	0.01 // nm/nN

// iterations
#define ITERATIONS	30
//#define ITERATIONS	1

// neighborhood radius
#define NEIGHBORHOOD_RADIUS	100 // nm

// edge number per object compensation
// accounts for bubbles
#define EDGE_COMPENSATION	2000

// threshold for computing vertex niceness
#define NICE_THRESHOLD			15 // nm

// set frozen objects flag and 
// specify frozen objects by name
//bool FROZEN_OBJECTS	= true;
//char FROZEN_LIST[SIZE_OF_OBJECT_NAME] = "ER_filtered_50000";

// handle clock values
//#ifndef CLK_TCK
//#define CLK_TCK CLOCKS_PER_SEC
//#endif

//####################################################
//############# files and directories ################
//####################################################

//#define INPUT_DATA_DIR    "/spin/cnl/mcell/justin/mesh_data/single_synapse/filtered_meshes/"
//#define INPUT_DATA_DIR    "/spin/cnl/mcell/justin/mesh_data/full_model/raw_meshes/"
//#define INPUT_DATA_DIR    "/spin/cnl/mcell/justin/mesh_data/single_synapse/temp_meshes/"
//#define OUTPUT_DATA_DIR    "/home/ezhang/SC_CA1/Ea2/mito_test"
#define INPUT_DATA_DIR    "/spin/cnl/mcell/justin/elaine/input/"
#define OUTPUT_DATA_DIR    "/spin/cnl/mcell/justin/elaine/output/"
//#define OUTPUT_DATA_DIR    "/spin/cnl/mcell/justin/mesh_data/single_synapse/manip_meshes/"
//#define OUTPUT_DATA_DIR    "/spin/cnl/mcell/justin/mesh_data/full_model/raw_meshes/"

#define MAIN_LOG_FILE      "main.log"
#define SCAN_LOG_FILE      "scan.log"
#define SPACE_LOG_FILE     "space.log"
#define MANIP_LOG_FILE     "manip.log"
#define CONT_LOG_FILE      "cont.log"
#define OBJECT_LOG_FILE    "object_"
#define VERTEX_LOG_FILE    "vertex_"

//####################################################
//####################################################
//####################################################
