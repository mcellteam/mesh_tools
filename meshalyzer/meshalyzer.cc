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
using std::cerr;
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
	"       meshalyzer - mesh quality analyzer\n"+
	"\nSYNOPSIS\n"+
	"       meshalyzer [options] FILE|DIR\n"+
	"\nDESCRIPTION\n"+
	"       Meshalyzer is a general purpose mesh analyzer useful for\n"+
	"       generating a complete summary of the current state of a mesh.\n"+
	"       Meshalyzer assesses mesh integrity (e.g. missing data),\n"+
	"       mesh attributes (e.g. closed, manifold, oriented), and\n"+
	"       mesh characteristics (e.g. number of vertices, faces, edges).\n"+
	"       Batch processing is easy by passing a directory name\n"+
	"       as input on the command line.\n"+
	"\nEXAMPLES\n"+
	"       meshalyzer filename\n"+
	"              Evaluate mesh integrity, attributes,\n"+
	"              and characteristics for the single mesh file.\n\n"+
	"       meshalyzer directoryname\n"+
	"              Evaluate mesh integrity, attributes,\n"+
	"              and characteristics for each single mesh file in directory.\n\n"+
	"       meshalyzer -p filename\n"+
	"              Evaluate mesh integrity, attributes, and characteristics \n"+
	"              for the single mesh file, print the results, and print \n"+
	"              the mesh elements preventing the mesh from being good with\n"+
	"              regards to the mesh characteristics and the attributes, if any.\n\n"+
	"       meshalyzer -a -p filename\n"+
	"              Evaluate the five mesh attributes for the single mesh file,\n"+
	"              print the state of each attribute, and print the mesh elements\n"+
	"              preventing the mesh from being good with regards to the attributes, if any.\n\n"+
	"       meshalyzer -b 10.0 -c 1.0 -p filename\n"+
	"              Evaluate mesh integrity, attributes, and characteristics \n"+
	"              for the single mesh file, print the results, and print \n"+
	"              the mesh elements preventing the mesh from being good with\n"+
	"              regards to the mesh characteristics and the attributes, if any.\n"+
	"              Additionally, screen faces with aspect ratios larger than 10.0 and\n"+
	"              screen edges with lengths larger than 1.0, and print detailed\n"+
	"              information about the offending mesh elements.\n"+
	"\nOPTIONS\n"+
	"       -a\n"+
	"              Evaluate the attributes of the mesh and report the results.\n"+
	"              Skip the evaluation of mesh characteristics.\n\n"+
	"       -b NUM\n"+
	"              Detect edges with length smaller than NUM.\n"+
	"              Units are same as FILE.\n\n"+
	"       -c NUM\n"+
	"              Detect edges with length greater than NUM.\n"+
	"              Units are same as FILE.\n\n"+
	"       -d NUM\n"+
	"              Detect edges with angle between two adjacent faces\n"+
	"              smaller than NUM degrees.\n\n"+
	"       -e NUM\n"+
	"              Detect edges with angle between two adjacent faces\n"+
	"              greater than NUM degrees.\n\n"+
	"       -f NUM\n"+
	"              Detect faces with aspect ratio greater than NUM.\n\n"+
	"       -h\n"+
	"              Print meshalyzer man page.\n\n"+
	"       -i\n"+
	"              Detect intersections between faces from different objects.\n"+
	"              Faceintersection detection is performed once all objects\n"+
	"              are loaded into memory. Single object intersection detection\n"+
	"              is omitted.\n\n"+
	"       -p\n"+
	"              Print detailed information about offending mesh elements \n"+
	"              (i.e. flipped faces, borders, nonmanifold edges,\n"+
	"              nonmanifold vertices, intersecting faces).\n\n"+
	"       -q\n"+
	"              Same as '-p' option, but prints vertex information\n"+
	"              in dreamm custom points format.\n"+
	"       -v\n"+
	"              If folder passed as argument, then only print total\n"+
	"              set volume and nothing else.\n"+
	"\nJustin Kinney				2007/10/01\n";

	// instantiate controls class
	Controls cs;

	// parse command line 
//	char filename[FILENAME_SIZE];
//	cs.parse(argc,argv,message,filename);
	cs.parse(argc,argv,message);

	// create container, objects, vertices, faces, edges, and find adjacencies
	Container c;

	// initialize space data structure
	Space s;

	// if single input file was found
	if(cs.folder==false){
		// save filename
//		c.files.push_back(filename);
		c.files.push_back(cs.inpath);
		// update index
		c.num_files++;
		// build data structure of mesh
		Object *obj=c.processFile(c.files[0]);
		// if either no vertices or no faces were found
		if (obj!=NULL) {
			// check mesh integrity
			obj->checkIntegrity();
			// if mesh file is not fatally flawed
			// i.e. no missing or orphaned vertices
			// and must be contiguous vertex and face indexing.
			// Note score may be set to 1 for other reasons in analyze()
			if(obj->goodIntegrity()==true){
				obj->createEdges();
				obj->findVertexAdjacencies();
				// partition space
				c.boundWorld(s,cs);
				// DEBUG
//				cout << "\nmain: "
//				<< "c.countFace = " << c.countFace()
//				<< ", obj->f.size()=" << obj->f.size() << endl;
				// DEBUG
//				s.initBoxes(c.countFace());
				s.initBoxes(obj->f.size());
				c.assignFacesToBoxes(s);
				// analyze and respond
				obj->analyze(&c,cs,s);
			}
			obj->print(cs);
		}
	} else { 
		// else scan folder
//		c.scanDir(filename);
		c.scanDir(cs.inpath.c_str());
		// for each file in folder
		for (int count=0;count<c.num_files;count++) {
			// build data structure of mesh
			Object *obj=c.processFile(cs.inpath+c.files[count]);
			// if neither zero vertices nor zero faces were found
			if (obj!=NULL){
				obj->checkIntegrity();
				if(obj->goodIntegrity()==true){
					obj->createEdges();
					obj->findVertexAdjacencies();
					// partition space
					c.boundWorld(s,cs);
					s.deleteBoxes();
//					s.initBoxes(c.countFace());
					s.initBoxes(obj->f.size());
					c.assignFacesToBoxes(s);
					// analyze and respond
					obj->analyze(&c,cs,s);
					c.update(obj);
					cs.store(obj);
				} else {
					cs.good_integrity=false;
				}
				// print stats
				if(cs.vol==false){obj->print(cs);}
				// if no request for inter-object
				// face intersection detection
					if(cs.interf==false){
					// then save space, delete object
					delete obj;
						c.o.clear();
				}
			}
		}
		c.boundWorld(s,cs);
		// analyze data as set
//		cs.analyzeCumulative(c);
		cs.analyzeCumulative();
		// print cumulative surface area, volume, 
		if(cs.vol==false){cs.printCumulative(c);}
                else {cout << c.countVol() << endl;}
		// detect intersecting faces between objects
		if(cs.interf==true){
			// partition space
			s.deleteBoxes();
			s.initBoxes(c.countFace());
			c.assignFacesToBoxes(s);
			// analyze and respond
			cs.analyzeBatch(&c,s);
			c.printBatch(cs);
		}
	}
}
