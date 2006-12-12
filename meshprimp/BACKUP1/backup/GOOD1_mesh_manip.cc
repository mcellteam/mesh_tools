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

//#include <stdio.h>
//#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <fstream>
#include <string>
#include <dirent.h>
#include <errno.h>
#include <time.h>
#include <sys/timeb.h>

using namespace std;

//####################################################
//#################### parameters ####################
//####################################################

// set to true to increase accuracy of line/polygon intersections at expense of speed
// set to flase to increase speed at expense of accuracy during line/polygon intersections
bool detect_polygon_edge_intersection = false;

// set to true to drive polygon angles to flat
// set to true to drive polygon angles to initial angle
bool init_angle_to_flat = true;

// switch to control printing of each vertex's stats
bool vertex_print = false;

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
#define SPACE_LENGTH	40 // nm

// closest vertex parameters
#define NUM_ADJACENT_BOXES  1
#define SEARCH_RADIUS  40 // nm

// energy function parameters
#define TARGET_SEPARATION		20   // nm
#define SEPARATION_WEIGHT		10.0  // nN/nm
#define EDGE_STRETCH_WEIGHT 	1  // nN/nm
#define ANGLE_STRETCH_WEIGHT 	10.0  // nN
#define TIME_STEP				0.01 // s
#define DAMPING					1	 // nN/(nm/s)
#define GAIN_SCALE				1.1 // nN/(nm/s)
//#define DISPLACEMENT_WEIGHT 	0.01 // nm/nN

// iterations
#define ITERATIONS	1

// neighborhood radius
#define NEIGHBORHOOD_RADIUS	100 // nm

// edge number per object compensation
// accounts for bubbles
#define EDGE_COMPENSATION	2000

// threshold for computing vertex niceness
#define NICE_THRESHOLD			15 // nm

// handle clock values
//#ifndef CLK_TCK
//#define CLK_TCK CLOCKS_PER_SEC
//#endif

//####################################################
//############# files and directories ################
//####################################################

//#define INPUT_DATA_DIR    "/spin/cnl/mcell/justin/mesh_data/single_synapse/filtered_meshes/"
#define INPUT_DATA_DIR    "/spin/cnl/mcell/justin/mesh_data/full_model/raw_meshes/"
//#define INPUT_DATA_DIR    "/spin/cnl/mcell/justin/mesh_data/single_synapse/temp_meshes/"
//#define OUTPUT_DATA_DIR    "/spin/cnl/mcell/justin/mesh_data/single_synapse/manip_meshes/"
#define OUTPUT_DATA_DIR    "/spin/cnl/mcell/justin/mesh_data/full_model/raw_meshes/"

#define MAIN_LOG_FILE      "main.log"
#define SCAN_LOG_FILE      "scan.log"
#define SPACE_LOG_FILE     "space.log"
#define MANIP_LOG_FILE     "manip.log"
#define CONT_LOG_FILE      "cont.log"
#define OBJECT_LOG_FILE    "object_"
#define VERTEX_LOG_FILE    "vertex_"

//####################################################
//#################### includes ######################
//####################################################

#include </home/jkinney/recon/definitions.cc>
#include </home/jkinney/recon/subroutines.cc>

//####################################################
//###################### main  #######################
//####################################################

main(){

	int k, i, w;
	char log_file[50], str[128];
	time_t currtime;
	double gain;

	int deleteme = 0;

	// initialize time variable
	currtime = time(NULL);
	
	// open log file
	sprintf(log_file,"%s%s",OUTPUT_DATA_DIR,MAIN_LOG_FILE);
	ofstream myfile (log_file);

	// declare a pointer (objects) to an array of size zero
	// whose future elements will be of type ObjectClass
	ContainerClass container;

	// declare a pointer (polygons) to an array of size zero
	// whose future elements will be of type polygon
	PolygonClass *polygons; 

	sleep(deleteme);

	////////// scan input file //////////
	// scan file and count number of objects and 
	// store number of vertices and polygons for each object
	cout << "\n\nscanning file...........................";
	cout.flush();
	ScanClass scan;
	scan.scanDir();
	scan.scanFiles(container,polygons,1);
	cout << "complete.\n";
	cout.flush();
	currtime = recordTime(myfile, currtime,"scanning file:");
	sleep(deleteme);

	////////// transfer scan data to container class //////////
	container.object_count = scan.object_count;

	////////// declare data structure //////////
	cout << "declaring object class data structure...";
	cout.flush();
	container.init();
	container.initObjects(scan.vertices_array,scan.polygons_array);
	cout << "complete.\n";
	cout.flush();
	currtime = recordTime(myfile,currtime,"Declare object class data structure:");
	sleep(deleteme);


	// allocate memory for polygons array with each 
	// array element of type polygon (polygons)
	cout << "declaring polygons data structure.......";
	cout.flush();
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
	cout << "complete.\n";
	cout.flush();
	currtime = recordTime(myfile,currtime,"Declare polygons data structure:");
	sleep(deleteme);


	////////// get data //////////
	// load objects with data from file
	cout << "get data................................";
	cout.flush();
	scan.scanFiles(container,polygons,0);
	cout << "complete.\n";
	cout.flush();
	currtime = recordTime(myfile,currtime,"Get data from files:");
	sleep(deleteme);


	cout << "initialize edges data structure.........";
	cout.flush();
	container.initObjectEdges(polygons,scan.polygon_count,scan.polygons_array);
	cout << "complete.\n";
	currtime = recordTime(myfile,currtime,"Initialize edges data strucutre:");
	sleep(deleteme);


	cout << "initialize new coords and dev_dist......";
	cout.flush();
	container.initCoordsAndDist();
	cout << "complete.\n";
	cout.flush();
	currtime = recordTime(myfile,currtime,"Initialize new coordinates:");
	sleep(deleteme);


	cout << "find adjacent polygons and vertices.....";
	cout.flush();
	container.initFindAdj(polygons,scan.polygon_count,scan.polygons_array);
	cout << "complete.\n";
	currtime = recordTime(myfile,currtime,"Find adjacent polygons and vertices:");
	sleep(deleteme);


	cout << "compute edge lengths....................";
	cout.flush();
	container.computeEdges();
	cout << "complete.\n";
	cout.flush();
	currtime = recordTime(myfile,currtime,"Compute edge lengths:");
	sleep(deleteme);

	
	cout << "find neighborhood vertices..............";
	cout.flush();
	container.initFindNeighbor();
	cout << "complete.\n";
	currtime = recordTime(myfile,currtime,"Find neighborhood vertices:");
	sleep(deleteme);


	////////// initializing space data structure //////////
	cout << "initializing space data structure.......";
	cout.flush();
	SpaceClass space;
	space.boundWorld(container,scan.object_count);
	space.boxes = NULL;
	space.init();
	cout << "complete.\n";
	cout.flush();
	currtime = recordTime(myfile,currtime,"Initialize space data structure:");

	sleep(deleteme);

	////////// initializing manip class //////////
	cout << "instantiating manip class...............";
	cout.flush();
	ManipClass manip;
//	manip.fileInit();
	cout << "complete.\n";
	cout.flush();
	currtime = recordTime(myfile,currtime,"Initialize manip class:");

	sleep(deleteme);

	////////// transfer scan data to other classes //////////
	manip.object_count = scan.object_count;
	manip.vertex_count = scan.vertex_count;
	manip.polygon_count = scan.polygon_count;

	sleep(deleteme);

	cout << "compute polygon normals.................";
	cout.flush();
	manip.computePolygonNormals(container,polygons);
	cout << "complete.\n";
	cout.flush();
	currtime = recordTime(myfile,currtime,"Compute polygon normals:");

	sleep(deleteme);

	cout << "initialize edge angles..................";
	cout.flush();
	container.initObjectEdgeAngles(polygons);
	cout << "complete.\n";
	currtime = recordTime(myfile,currtime,"Initialize edge angles:");

	sleep(deleteme);

	////////// assign polygons to boxes //////////
	cout << "assign polygons to boxes................";
	cout.flush();
	manip.assignPolygonsToBoxes(container,polygons,space);
	cout << "complete.\n";
	cout.flush();
	currtime = recordTime(myfile,currtime,"Assign polygons to boxes:");

	////////// check reconstruction for errors //////////
//	cout << "search for self-intersecting meshes.....";
//	cout.flush();
//	manip.searchForSelfIntersection(container,polygons,space);
//	cout << "complete.\n";
//	cout.flush();
//	currtime = recordTime(myfile,currtime,"Search for self-intersecting meshes:");

	sleep(deleteme);

	////////// initialize and write data to log files //////////
	cout << "writing object data to log files........";
	cout.flush();
	space.fileInit();
	space.updateFile(0);
	scan.fileInit();
	scan.updateFile();
	container.fileInit();
	container.initObjectAndVertexLogFiles();
	cout << "complete.\n\n\n";
	cout.flush();
	currtime = recordTime(myfile,currtime,"Writing data to log files:");

	// initialize gaim
	gain = TIME_STEP/DAMPING;
	sleep(deleteme);

	////////// BEGIN LOOP //////////
	for (w=0;w<ITERATIONS;w++) {

		// if no self-intersecting meshes were found
//		if (!manip.self_intersection) {

			cout << "Iteration " << w << ": ";
			cout << "compute polygon normals...........";
			cout.flush();
			manip.computePolygonNormals(container,polygons);
			cout << "complete.\n";
			cout.flush();
			sprintf(str,"%s%2d%s","Iteration ",w,": Compute polygon normals:");
			currtime = recordTime(myfile,currtime,str);
			sleep(deleteme);

			cout << "Iteration " << w << ": ";
			cout << "compute New edge lengths..........";
			cout.flush();
			container.computeNewEdges();
			cout << "complete.\n";
			cout.flush();
			sprintf(str,"%s%2d%s","Iteration ",w,": Compute new edge lengths:");
			currtime = recordTime(myfile,currtime,str);
			sleep(deleteme);

			cout << "Iteration " << w << ": ";
			cout << "find nice vertices................";
			cout.flush();
			manip.findNice(container,polygons, space);
			cout << "complete.\n";
			cout.flush();
			sprintf(str,"%s%2d%s","Iteration ",w,": Find nice vertices:");
			currtime = recordTime(myfile,currtime,str);
			sleep(deleteme);

			cout << "Iteration " << w << ": ";
			cout << "set deviation distance............";
			cout.flush();
	        manip.setDeviationDistance(container,polygons,space);
			cout << "complete.\n";
			cout.flush();
			sprintf(str,"%s%2d%s","Iteration ",w,": Set deviation distance:");
			currtime = recordTime(myfile,currtime,str);
			sleep(deleteme);
	
			cout << "Iteration " << w << ": ";
			cout << "compute edge energy and force.....";
			cout.flush();
	        container.computeObjectEdgeAngleForceEnergy(polygons);
			cout << "complete.\n";
			cout.flush();
			sprintf(str,"%s%2d%s","Iteration ",w,": Compute edge energy and force:");
			currtime = recordTime(myfile,currtime,str);
			sleep(deleteme);
		
			cout << "Iteration " << w << ": ";
			cout << "compute vertex energy and force...";
			cout.flush();
	        manip.computeVertexEnergyForce(container);
			cout << "complete.\n";
			cout.flush();
			sprintf(str,"%s%2d%s","Iteration ",w,": Compute vertex energy and force:");
			currtime = recordTime(myfile,currtime,str);
			sleep(deleteme);
		
			cout << "Iteration " << w << ": ";
			cout << "compute new vertex coordinates....";
			cout.flush();
	        manip.computeNewVertexCoords(container,gain);
			cout << "complete.\n";
			cout.flush();
			sprintf(str,"%s%2d%s","Iteration ",w,": Compute new vertex coordinates:");
			currtime = recordTime(myfile,currtime,str);
			sleep(deleteme);
	
			cout << "Iteration " << w << ": ";
			cout << "compute world bounds..............";
			cout.flush();
//			space.boundWorld(container,scan.object_count);
//			space.init();
			cout << "complete.\n";
			cout.flush();
			sprintf(str,"%s%2d%s","Iteration ",w,": Compute world bounds:");
			currtime = recordTime(myfile,currtime,str);
			sleep(deleteme);

			////////// assign polygons to boxes //////////
			cout << "Iteration " << w << ": ";
			cout << "assign polygons to boxes..........";
			cout.flush();
//			manip.assignPolygonsToBoxes(container,polygons,space);
			cout << "complete.\n";
			cout.flush();
			sprintf(str,"%s%2d%s","Iteration ",w,": Assign polygons to boxes:");
			currtime = recordTime(myfile,currtime,str);
			sleep(deleteme);

			////////// compute global parameters //////////
			cout << "Iteration " << w << ": ";
			cout << "compute global parameters.........";
			cout.flush();
//			container.computeGlobalParams();
			cout << "complete.\n";
			cout.flush();
			sprintf(str,"%s%2d%s","Iteration ",w,": Compute global parameters:");
			currtime = recordTime(myfile,currtime,str);
			sleep(deleteme);

			////////// write data to log files //////////
			cout << "Iteration " << w << ": ";
			cout << "update log files..................";
			cout.flush();
			space.updateFile(w+1);
			container.updateFile(w+1);
			container.updateObjectAndVertexLog(w+1);
			cout << "complete.\n";
			cout.flush();
			sprintf(str,"%s%2d%s","Iteration ",w,": Update log files:");
			currtime = recordTime(myfile,currtime,str);
			sleep(deleteme);

			// initialize gaim
			gain = gain*GAIN_SCALE;
	
//		}
	
	}
			////////// check reconstruction for errors //////////
//			cout << "\nsearch for self-intersecting meshes...";
//			cout.flush();
//			manip.searchForSelfIntersection(container,polygons,space);
//			cout << "complete.\n";
//			cout.flush();


	////////// close log files //////////
//	manip.fileOutit();
	space.fileOutit();
	scan.fileOutit();
	container.fileOutit();
	container.OutitObjectAndVertexLogFiles();
	myfile.close();

	////////// build 'after picture' as mesh //////////
	cout << "\nbuild Mesh after...";
	cout.flush();
	scan.buildMeshAfter(container,polygons);
	cout << "complete.\n";
	cout.flush();

	////////// write deviation distances to file //////////
	cout << "\nsave deviation distances...";
	cout.flush();
	container.writeVertexClosestDistanceToFile();
	cout << "complete.\n";
	cout.flush();

	delete[] polygons;

}
