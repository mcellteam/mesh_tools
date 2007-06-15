#include <math.h>
#include <iostream>
#include <vector>
#include <set>
#include <list>
#include <algorithm>
#include <fstream>
#include <string>
#include <dirent.h>
#include <errno.h>
#include <time.h>
#include <sys/timeb.h>
#include <sys/time.h>
#include <ext/hash_map>
#include <ext/hash_set>
#include <numeric>
#include <map>
#include <ncurses.h>

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

//int main(int argc,char **argv){
int main(){

	copyControlFile();

	///// check that assumption of 32 bit int is correct /////
	if (32!=sizeof(int)*8){cout << "Error. Int is not 32 bits, sizeof(int) " << sizeof(int) << endl;exit(0); }

/*    if (argc != 2 && argc != 3 )
    {
        printf("\nSyntax: meshprimp [-d]\n\n");
        printf("Description: Identifies intersections between faces in each input file.\
n");
        printf("The '-d' option prints a summary report of face intersections to stdout.
\n");
        printf("The '-D' option prints intersecting face vertices in DReAMM custom point
s format to stdout.\n");
        printf("If 4all, then all mesh files in current directory are analyzed.\n\n");
        return 1;
    }*/

    ////////// declare variables /////////
	time_t currtime = time(NULL);
	
	// open log file
	char log_file[50];
	sprintf(log_file,"%s%s",OUTPUT_DATA_DIR,MAIN_LOG_FILE);
	std::ofstream myfile (log_file);

	// instantiate statistics class
	Monitor stats;

	// create container, objects, vertices, faces, edges, and find adjacencies
	Container c;
	currtime = recordTime(myfile, currtime,"scanning file:");

	// initialize space data structure
	Space s(c);
	currtime = recordTime(myfile,currtime,"Initialize space data structure:");

	// assign faces to boxes
	c.assignFacesToBoxes(s);
	currtime = recordTime(myfile,currtime,"Assign faces to boxes:");

	////////// initialize and write data to log files //////////
	cout << "writing object data to log files...............";
	cout.flush();
	c.writeObjectList();
	c.statusFileInit();
	cout << "complete.\n\n";
	cout.flush();
	currtime = recordTime(myfile,currtime,"Writing data to log files:");

	c.clear();

	c.findNice(s);
	currtime = recordTime(myfile,currtime,"Iteration 0: Find nice vertices:");

	c.getSeparationDistances(s,stats);
	currtime = recordTime(myfile,currtime,"Iteration 0: Get separation distances:");

	c.computeFaceIntersectionForce();
	currtime = recordTime(myfile,currtime,"Iteration 0: Compute face intersection force:");
	
	// update container stats
	c.updateStats(0.0);

	c.computeGlobalEnergy();

//	stats.printVerticesWithCP();
	c.checkEdgeAngles();

	// update container class log file
	c.updateFile(0,false,0);
	currtime = recordTime(myfile,currtime,"Iteration 0: Update log files:");

/*
//	int num=0;
//	std::vector<Object*>::iterator x = c.o.begin();
	// for each object
	for(std::vector<Object*>::iterator z = c.o.begin();z!=c.o.end();z++){
		// if object is d000
		if(!strcmp((*z)->name,"d000_FILTERED_SMOOTH")){
			// for each vertex in object
			for(std::vector<Vertex*>::iterator y = (*z)->v.begin();y!=(*z)->v.end();y++){
				printNeighborhood((*y)->o->name,(*y)->index,(*y)->nf,0);
				if(++num>10){break;}
			}
		}
	}
*/
//	exit(0);

	c.writeDistances(0);
//	c.writeDistances50(0);
//	c.writeDistancesNOCP(0);
//	exit(0);
	////////// BEGIN LOOP //////////
//	double min_tim1=1E30,max_tim1=-1E30;
//	double min_tim2=1E30,max_tim2=-1E30;
	// time the loop
	timeval tim;
	int interval =0;
	for (int w=0;w<ITERATIONS;w++) {
		// load topN, save old, and initialize avg
		stats.prep(&c);

		// avoid oscillations
		hashtable_v gate;
		stats.initRefrac();

//		stats.validateTopN("main: top of for");		
		double ref_gain = c.gain;
		stats.clearAvg();
		stats.updateAvg(c.energy);
		c.clear();
		// point multimap iterator to pair with largest separation error
		tv_iterator tvi = stats.topN.begin();
		stats.touchs = 0;
		Vertex* last=NULL,*pen=NULL,*third=NULL;
		// until global energy reaches steady-state
		Vertex *cv=NULL;
		gettimeofday(&tim,NULL);
		double t1=tim.tv_sec+(tim.tv_usec/1000000.0);
					char name[1024];
		int count = 1;
		while ( distinguishable(stats.avg_new,stats.avg_old,1E-12) && count<MAX_VERTEX_MOVES){
//		while ( distinguishable(stats.avg_new,stats.avg_old,1E-12) && count<45000){
//			stats.validateTopN("main: top of while");		
			// if vertex is candidate, i.e. closest point was found for vertex
			// and vertex is different from the last two vertices moved
//			if((*tvi).second->cl!=NULL && (*tvi).second!=last && (*tvi).second!=pen && (*tvi).second!=third){
			if((*tvi).second->cl!=NULL && (gate.find((*tvi).second)==gate.end())
				&& !stats.Refracted((*tvi).second) ){
//			if((*tvi).second->cl!=NULL){
				cv=(*tvi).second;
				// compute new vertex coords
				double pH[3]; // new holding position coordinates (x,y,z)
				cv->computeNewCoords(&c,pH,ref_gain);
				// assign new vertex coords
//				double store = (*tvi).first;
				double store = (pH[0]-cv->pN[0])*(pH[0]-cv->pN[0])+
									(pH[1]-cv->pN[1])*(pH[1]-cv->pN[1])+
									(pH[2]-cv->pN[2])*(pH[2]-cv->pN[2]);
				if(c.assignNewVertexCoords(s,cv,pH,stats)){
//					gettimeofday(&tim,NULL);
//					double t2=tim.tv_sec+(tim.tv_usec/1000000.0);
//					double diff_tim = t2-t1;
//					if(diff_tim<min_tim1){min_tim1=diff_tim;}
//					if(diff_tim>max_tim1){max_tim1=diff_tim;}
					third=pen;
					pen=last;
					last=cv;
					// increment current vertex touch count
					stats.updateTouchMap(cv);
					// update sliding energy average
					stats.updateAvg(c.energy);
					// print
					if(!write_mesh_every_iteration && ++stats.touchs==1000){
						stats.touchs = 0;
						stats.printVertexSelect(c,++interval);
						stats.touch_map.clear();
					}
					// update container stats
					c.updateStats(sqrt(store));
					if(print_flag){
						// update monitor stats
						stats.getBoxesPerFace(c);
						stats.getFacesPerBox(s);
					}
					// update gate
//					updateGate(cv,gate,10);
					stats.updateRefractoryWindow(cv);

					// print to stdout
//					char name[1024];
					if(!(count++%100)){
					sprintf(name,"count %d, ",count);
					cout << name;
					sprintf(name,"iter %d, ",w+1);
					cout << name;
					cout.width(45);
					sprintf(name,"%s->%d, ",cv->o->name.c_str(),cv->index);
					cout << left << name;
					cout.width(20);
					sprintf(name,"vd %.12g, ",sqrt(store));
					cout << left << name;
					cout.width(17);
					sprintf(name,"e %.6g, ",c.energy);
					cout << left << name;
					cout.width(20);
					sprintf(name,"e_avg %.6g, ",stats.avg_new);
					cout << left << name;
//					sprintf(name,"#topN %d, ",stats.topN.size());
//					cout << left << name;
					sprintf(name,"min edge angle %.6g, ",c.min_edge_angle);
					cout << left << name;
					if(print_flag){
						sprintf(name,"bpf (%.6g,%.6g,%.6g), ",stats.bpf_min,stats.bpf_mean,stats.bpf_max);
						cout << left << name;
						sprintf(name,"fpb (%.6g,%.6g,%.6g), ",stats.fpb_min,stats.fpb_mean,stats.fpb_max);
						cout << left << name;
					}
//					sprintf(name,"time1 (%.6g,%.6g,%.6g), ",min_tim1,diff_tim,max_tim1);
//					cout << left << name;
//					gettimeofday(&tim,NULL);
//					double t3=tim.tv_sec+(tim.tv_usec/1000000.0);
//					diff_tim = t3-t2;
//					if(diff_tim<min_tim2){min_tim2=diff_tim;}
//					if(diff_tim>max_tim2){max_tim2=diff_tim;}
//					sprintf(name,"time2 (%.6g,%.6g,%.6g), ",min_tim2,diff_tim,max_tim2);
//					cout << left << name;
//					sprintf(name,"gain %.3g",c.gain);
//					cout << left << name;
					cout << endl;
					}
					// update iterator
					tvi = stats.topN.begin();
//					c.gain=ref_gain;
					ref_gain=c.gain;
				} else {
					if(!(count%100)){
					cout << "vertex " << cv->index << ", vd " << sqrt(store)
					<< ", not moved due to bad geometry, gain " << c.gain << endl;
					}
//					if(c.gain>DOUBLE_EPSILON){
					ref_gain=ref_gain/2.0;
//					if(store>MIN_DISPLACEMENT_SQ){
//						c.gain=c.gain/2.0;
					//	ref_gain=ref_gain/2.0;
//					} else {
					if(store<MIN_DISPLACEMENT_SQ){
						// update gate
						punishGate(cv,gate,static_cast<int>(MAX_VERTEX_MOVES));
//						c.gain=ref_gain;
						ref_gain=c.gain;
						// try next worst vertex in set
						tvi++;
					}
				}
			} else {
				if((*tvi).second->cl==NULL){
					cout << "vertex " << (*tvi).second->index << ", vd " << sqrt((*tvi).first)
					<< ", not moved because there is no closest point.\n";exit(0);
//					cout << "vertex " << (*tvi).second->index << ", vd " << sqrt((*tvi).first)
//					<< ", not moved because there is no closest point.\n";
				} /*else {
					cout << "vertex " << (*tvi).second->index << ", vd " << sqrt((*tvi).first)
					<< ", not moved because it is frozen.\n";
				}*/
				// try next worst vertex in set
				tvi++;
			}
	/*		gettimeofday(&tim,NULL);
			double t3=tim.tv_sec+(tim.tv_usec/1000000.0);
			double diff_tim = t3-t1;
			if(diff_tim<min_tim2){min_tim2=diff_tim;}
			if(diff_tim>max_tim2){max_tim2=diff_tim;}
			sprintf(name,"time2 (%.6g,%.6g,%.6g), ",min_tim2,diff_tim,max_tim2);
			cout << left << name << endl;
			gettimeofday(&tim,NULL);
			t1=tim.tv_sec+(tim.tv_usec/1000000.0);*/
		}

		c.gain *= GAIN_SCALE;
		ref_gain = c.gain;

		gettimeofday(&tim,NULL);
		double t2=tim.tv_sec+(tim.tv_usec/1000000.0);
		c.updateFile(w+1,last!=NULL,t2-t1);
		char str[128];
		sprintf(str,"%s%2d%s","Iteration ",w+1,": Update log files:");
		currtime = recordTime(myfile,currtime,str);

		if (write_mesh_every_iteration) {
			stats.printVertexSelect(c,w+1);
			cout << "Iteration " << w+1 << ": ";
			cout << "build Mesh after..................";
			cout.flush();
			c.buildMeshAfter(w+1);
			cout << "complete.\n";
			cout.flush();
		}
		if (write_distances_every_iteration) {
			cout << "Iteration " << w+1 << ": ";
			cout << "write closest point distances..................";
			cout.flush();
			c.writeDistances(w+1);
			cout << "complete.\n";
			cout.flush();
		}
	}

	stats.freeAvg();

	////////// close log files //////////
	c.fileOutit();
	myfile.close();

	if(!write_mesh_every_iteration){
		stats.printVertexSelect(c,"END");
		cout << "\nbuild Mesh after...";
		cout.flush();
		c.buildMeshAfter(1);
		cout << "complete.\n";
		cout.flush();
	}
	if (!write_distances_every_iteration) {
		cout << "write closest point distances..................";
		cout.flush();
		c.writeDistances(1);
		cout << "complete.\n";
		cout.flush();
	}


	exit(0);
}
