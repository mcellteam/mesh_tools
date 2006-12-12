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
//#include <vector>
#include <algorithm>
#include <fstream>
#include <string>
#include <dirent.h>
#include <errno.h>
#include <time.h>
#include <sys/timeb.h>

using namespace std;

//####################################################
//#################### includes ######################
//####################################################

#include </home/jkinney/recon/dkeller/controls.cc>
#include </home/jkinney/recon/dkeller/definitions.cc>
#include </home/jkinney/recon/dkeller/subroutines.cc>

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
	container.objectList();
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
			cout << "compute polygon intersection force.....";
			cout.flush();
	        manip.computePolygonIntersectionForce(container,polygons,space);
			cout << "complete.\n";
			cout.flush();
			sprintf(str,"%s%2d%s","Iteration ",w,": Compute polygon intersection force:");
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
			space.boundWorld(container,scan.object_count);
			space.init();
			cout << "complete.\n";
			cout.flush();
			sprintf(str,"%s%2d%s","Iteration ",w,": Compute world bounds:");
			currtime = recordTime(myfile,currtime,str);
			sleep(deleteme);

			////////// assign polygons to boxes //////////
			cout << "Iteration " << w << ": ";
			cout << "assign polygons to boxes..........";
			cout.flush();
			manip.assignPolygonsToBoxes(container,polygons,space);
			cout << "complete.\n";
			cout.flush();
			sprintf(str,"%s%2d%s","Iteration ",w,": Assign polygons to boxes:");
			currtime = recordTime(myfile,currtime,str);
			sleep(deleteme);

			////////// compute global parameters //////////
			cout << "Iteration " << w << ": ";
			cout << "compute global parameters.........";
			cout.flush();
			container.computeGlobalParams();
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
	
			if (write_mesh_every_iteration) {
				////////// build 'after picture' as mesh //////////
				cout << "\nbuild Mesh after...";
				cout.flush();
				scan.buildMeshAfter(container,polygons);
				cout << "complete.\n";
				cout.flush();
			}

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

	////////// write deviation distances to file //////////
	cout << "\ncopy control settings...";
	cout.flush();
	copyControlFile();
	cout << "complete.\n";
	cout.flush();

	// clean up
	delete[] polygons;

}
