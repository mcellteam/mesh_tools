#include <float.h>
#include <math.h>
#include <iostream>
#include <vector>
#include <set>
#include <list>
#include <algorithm>
#include <fstream>
#include <string>
#include <dirent.h>
#include <ext/hash_map>
#include <ext/hash_set>
#include <numeric>
#include <map>
#include <sys/stat.h>
#include <unistd.h>

using std::map;
using std::vector;
using std::cout;
using std::endl;
using std::left;
using std::set;

//####################################################
//#################### includes ######################
//####################################################

#include "controls.h"
#include "classes.cc"
#include "subroutines.cc"

//####################################################
//###################### main  #######################
//####################################################

int main(int argc,char **argv){
	std::string message = "\n";
	message=message+
	"NAME\n"+
	"       mesh2contours - extract Reconstruct contours from mesh\n"+
	"\nJustin Kinney				2007/10/01\n";

	// instantiate controls class
	Controls cs;

	// parse command line 
	cs.parse(argc,argv,message);

	// create container, objects, vertices, 
	// faces, edges, and find adjacencies
	Container c;

	// if single input file was found
	if(cs.folder==false){
		// save filename
		c.files.push_back(cs.inpath);
		// update index
		c.num_files++;
		// build data structure of mesh
		Object *obj=c.processFile(c.files[0]);
		// if either no vertices or no faces were found
		if (obj!=NULL) {
			// look for contours in object
			obj->findContours(cs.zEM);
			// check mesh integrity
			obj->print(cs);
		}
	} else { 
		// else scan folder
		c.scanDir(cs.inpath.c_str());
		// for each file in folder
		for (int count=0;count<c.num_files;count++) {
			// build data structure of mesh
			Object *obj=c.processFile(cs.inpath+c.files[count]);
			// if neither zero vertices nor zero faces were found
			if (obj!=NULL){
				// look for contours in object
				obj->findContours(cs.zEM);
				// check mesh integrity
				obj->print(cs);
//				c.update(obj);
//				cs.store(obj);
				// print stats
//				obj->print(cs);
				// if no request for inter-object
				// face intersection detection
				if(cs.interf==false){
					// then save space, delete object
					delete obj;
					c.o.clear();
				}
			}
		}
		// print cumulative surface area, volume, 
//		cs.printCumulative(c);
		// detect intersecting faces between objects
	}
}
