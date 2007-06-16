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
//#include <errno.h>
//#include <time.h>
//#include <sys/timeb.h>
//#include <sys/time.h>
#include <ext/hash_map>
#include <ext/hash_set>
#include <numeric>
#include <map>
#include <sys/stat.h>
//#include <ncurses.h>

//using namespace std;
//using namespace __gnu_cxx;
//using __gnu_cxx::map;
using std::map;
using std::vector;
using std::cout;
using std::endl;
using std::left;
using std::set;

//####################################################
//#################### includes ######################
//####################################################

#include "controls.cc"
#include "classes.cc"
#include "subroutines.cc"

//####################################################
//###################### main  #######################
//####################################################

int main(int argc,char **argv){

	char message[2048];
	sprintf(message,"\nUsage: meshalyzer file|folder ");
	sprintf(message,"%s[-a]|[-p] [-g=t1,t2,t3][-i][-s]\n\n",message);
	sprintf(message,"%sDescription: Analyze and assess the integrity of mesh files.\n\n",message);
	sprintf(message,"%sThe '-a' option evaluates the attributes of the meshes,\n",message);
	sprintf(message,"%s reports the results, and exits.\n\n",message);
	sprintf(message,"%sThe '-p' option prints offending mesh elements \n",message);
	sprintf(message,"%s(i.e. flipped faces, borders, nonmanifold edges,\n",message);
	sprintf(message,"%snonmanifold vertices, intersecting faces).\n\n",message);
	sprintf(message,"%sThe '-g=t1,t2,t3' option explicitly reports whether\n",message);
	sprintf(message,"%sor not mesh is good (i.e. mesh is closed manifold\n",message);
	sprintf(message,"%swith outward normals, no intersecting faces, no orphaned\n",message);
	sprintf(message,"%sor missing vertices, no degenerate faces, contiguous vertex\n",message);
	sprintf(message,"%sand face numbering, maximum face aspect ratio is less than\n",message);
	sprintf(message,"%suser-defined threshold, minimum edge angle is greater than\n",message);
	sprintf(message,"%suser-defined threshold, and minimum edge length is larger\n",message);
	sprintf(message,"%sthan user-defined threshold).\n",message);
	sprintf(message,"%st1=aspect ratio threshold, t2=edge angle threshold(radians)\n",message);
	sprintf(message,"%st3=edge length thresold in same units as input mesh file\n",message);
	sprintf(message,"%sIf no other options used, returns '1' or '0'\n",message);
	sprintf(message,"%sif mesh is good or not good, respectively.\n\n",message);
	sprintf(message,"%sThe '-i' option looks for intersections between\n",message);
	sprintf(message,"%sfaces from different objects.\n\n",message);
	sprintf(message,"%sThe '-s' option computes the seperation distance\n",message);
	sprintf(message,"%sfor every vertex in each object.\n\n",message);

	// instantiate controls class
	Controls cs;

	// parse command line 
	cs.parse(argc,argv,message);

	// adjust input directory
	char filename[1028];
	strcpy(filename,argv[1]);
	if(cs.folder){
    	char *temp=strrchr(filename,'/');
		if(!temp) {strcat(filename,"/");}
		else if(*++temp) {strcat(filename,"/");}
	}

	// create container, objects, vertices, faces, edges, and find adjacencies
	Container c;

	// initialize space data structure
//	Space s(c);
	Space s;

	// if single file
	Object *obj;
	bool badmesh=false;
	if(!cs.folder){
		// save filename
		c.files.push_back(filename);
		// update index
		c.num_files++;

		// build data structure of mesh
		c.scanFiles();
		// if mesh
		obj=c.o.front();
		obj->checkIntegrity(cs);
		// if mesh file is not fatally flawed
		// i.e. no missing or orphaned vertices
		// and must be contiguous vertex and face indexing.
		// Note score may be set to 1 for other reasons in analyze
		if(obj->score!=1){
			obj->createEdges();
//			cout << "pre adjacencies";cout.flush();
			obj->findVertexAdjacencies(cs);
			// partition space
//			cout << "pre boundWorld";cout.flush();
			c.boundWorld();
//			cout << "post boundWorld";cout.flush();
//			s.deleteBoxes();
			s.initBoxes(c);
			c.assignFacesToBoxes(s);
			// analyze and respond
			obj->analyze(&c,cs);
		} else {badmesh=true;}
		// update cumulative score
		cs.num_score+=obj->score;
		obj->print(cs,badmesh);
	} else { 
		// else scan folder
		c.scanDir(cs,filename);
		char file[1024];
		// for each file in folder
		for (int count=0;count<c.num_files;count++) {
			///// create an instance of Object class /////
			// copy char array to string
	        std::string str = c.files[count];
	        // record object name
			std::string::size_type pos1 = str.find(".",0);
			if (!(pos1==std::string::npos)) {
				// ALLOCATE MEMORY FOR NEW OBJECT
				obj = new Object(str.substr(0,pos1));
			} else { cout << "Error! Object name was not found in " << str << "\n";exit(1);}
			// save pointer to object
			c.o.push_back(obj);
			// count object in Controls
			cs.num_obj++;

			// build data structure of mesh
			sprintf(file,"%s%s",filename,c.files[count].c_str());
			c.scanFile(obj,file);
			cs.num_vert+=obj->v.size();
			cs.num_faces+=obj->f.size();
			obj->checkIntegrity(cs);
			if(cs.binaryOutput() && !cs.isGood(obj)){
				break;
			}
			if(obj->score!=1){
				obj->createEdges();
				cs.num_edges+=obj->e.size();
				obj->findVertexAdjacencies(cs);
	
				// partition space
				c.boundWorld();
				s.deleteBoxes();
				s.initBoxes(c);
				c.assignFacesToBoxes(s);
				// analyze and respond
				obj->analyze(&c,cs);
			} else {badmesh=true;}
			// update cumulative score
			cs.num_score+=obj->score;
			if(cs.binaryOutput() && !cs.isGood(obj)){
				break;
			}
			obj->print(cs,badmesh);
			cs.clear();
			// if NOT batch mode
			if(!cs.interf && !cs.sepdist){
 				// then continue as single file
				// delete object from container and clear space
				delete obj;
				s.clearBoxes();
				c.clear();
				badmesh=false;
			}
		}
		// print cumulative surface area, volume, 
		cs.printCumulative(badmesh,obj);
		// if batch mode
		if(cs.interf || cs.sepdist){
			if(cs.num_score == 0){
				// partition space
				c.boundWorld();
				s.deleteBoxes();
				s.initBoxes(c);
				c.assignFacesToBoxes(s);
				// analyze and respond
				obj->analyzeBatch(&c,cs,s);
				c.printBatch(cs);
			}
		}
	}
}
