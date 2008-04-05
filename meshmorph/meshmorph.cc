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
#include <cassert>

using std::map;
using std::vector;
using std::cout;
using std::endl;
using std::left;
using std::right;
using std::set;

std::string INPUT_DATA_DIR  = "./";
std::string OUTPUT_DATA_DIR = "./";

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

  if (argc != 3 && argc!=6)
  {
    //printf("\nSyntax: meshmorph input_dir output_dir [frozen_vertices_filename]\n\n");
    printf("\nSyntax: meshmorph input_dir output_dir [target_ecw sep_weight angle_weight]\n\n");
    printf("Default values for optional arguments are 20, 10, and 80, respectively.\n\n");
    //printf("Description: morphs without moving the frozen vertices.\n");
    //printf("In the optional frozen vertex file, specifiy one frozen vertex per\n");
    //printf("line with the syntax 'object_name vertex_index'.\n\n");
    return 1;
  }

  // check that assumption of 32 bit int is correct
  assert(checkIntSize());

  // declare time variable
  time_t currtime = time(NULL);
 // time_t begintime = currtime;

  // save input directory
  char filename[FILENAME_SIZE];
  strcpy(filename,argv[1]);
  char *temp=strrchr(filename,'/');
  if(!temp) {strcat(filename,"/");}
  else if(*++temp) {strcat(filename,"/");}
  INPUT_DATA_DIR = filename;

  // save output directory
  strcpy(filename,argv[2]);
  temp=strrchr(filename,'/');
  if(!temp) {strcat(filename,"/");}
  else if(*++temp) {strcat(filename,"/");}
  OUTPUT_DATA_DIR = filename;

  // open log file
  char log_file[FILENAME_SIZE];
  sprintf(log_file,"%s%s",OUTPUT_DATA_DIR.c_str(),MAIN_LOG_FILE);
  std::ofstream myfile (log_file);
  if(myfile.is_open()==false){cout << "\nFailed to open file " << log_file << ".\n"; }

  // create container, objects, vertices, faces, edges, and find adjacencies
  Container c;
  currtime = recordTime(myfile, currtime,"scanning file:");

  // copy control.cc to OUTPUT_DATA_DIR
  //copyControlFile();

  // read frozen vertices
  //if(argc>3){ c.readFrozen(argv[3]); }

  // read optional arguments
  if(argc==6)
  {
    char *eptr;
    TARGET_SEPARATION     = strtod(argv[3],&eptr);
    SEPARATION_WEIGHT     = strtod(argv[4],&eptr);
    ANGLE_STRETCH_WEIGHT  = strtod(argv[5],&eptr);
  }

  // initialize space data structure
  Space s(c);
  currtime = recordTime(myfile,currtime,"Initialize space data structure:");

  // assign faces to boxes
  c.assignFacesToBoxes(s);
  currtime = recordTime(myfile,currtime,"Assign faces to boxes:");

  // initialize and write data to log files
  c.writeObjectData();
  currtime = recordTime(myfile,currtime,"Writing data to log files:");

  // instantiate statistics class
  Monitor stats;

  // identify vertices that lie inside of another object
  c.findNice(s);
  currtime = recordTime(myfile,currtime,"Iteration 0: Find nice vertices:");

  // identify the closest point on a mesh to each vertex
  c.getSeparationDistances(s,stats);
  currtime = recordTime(myfile,currtime,"Iteration 0: Get separation distances:");

  // compute the vertex forces due to face intersections
  c.computeFaceIntersectionForce(s);
  currtime = recordTime(myfile,currtime,"Iteration 0: Compute face intersection force:");

  // compute the cumulative potential energy of all springs in the model	
  c.computeGlobalEnergy();

  // identify the smallest edge angle
  c.checkEdgeAngles();

  // update container class log file
  c.updateFile(0,false,0);
  currtime = recordTime(myfile,currtime,"Iteration 0: Update log files:");

  // save separation distances to file
  c.writeSeparationDistances();

  // write intersected faces to file
  stats.writeIntersected(0,&c);

  // create instance of moved-vertex tracking structure
  //VTrack pod;

/*  // DEBUG
  cout << "\n\n********** vertices with separation distance < 10 nm **********\n";
  // for each object in container
  for(std::vector<Object*>::iterator i=c.o.begin();i!=c.o.end();i++){
    // for each vertex in object
    for(std::vector<Vertex*>::iterator j=(*i)->v.begin();j!=(*i)->v.end();j++){
      if((*j)->cl!=NULL && (binary_search(c.frozen.begin(),c.frozen.end(),*j)==false)){
        // get closest point
        double pC[3];
        c.computePC((*j)->cl,*j,pC);
        // compute separation vector
        double ss[3];
        for(int ii=0;ii<3;ii++){ ss[ii]=pC[ii]-(*j)->pN[ii]; }
        // compute separation distance //////////
        double sd = sqrt(dot(ss,ss));
        if(sd<1.1){
          //cout << "main: sd = " << sd << endl;
          (*j)->printVertexCP();
          //cout << endl;
          //cout << "main: cp = ["
          //      << pC[0] << " "
          //      << pC[1] << " "
          //      << pC[2] << "]\n";
          //(*j)->cl->printFace((*j)->cl->v[0]->o->name);
          //cout << endl << endl;
        }
      }
    }
  }
  cout << "********** vertices with separation distance < 10 nm **********\n\n";
  exit(0);
*/  // DEBUG


  ////////// BEGIN LOOP //////////
  cout << endl;
  // for each group of iterations
  for (int group=1;group<(NUM_GROUPS+1);group++)
  {
    // enforce maximum runtime policy
    //if( (time(NULL)-begintime) > MAX_RUNTIME){break;}
    // initialize group variables
    stats.groupInit(c);
    // until GROUP_SIZE vertex moves have been made
    while (stats.count<GROUP_SIZE)
    {
      // enforce maximum runtime policy
      //if( (time(NULL)-begintime) > MAX_RUNTIME){break;}
      stats.identifyMeshRegionToUpdate(s,c);
      // for each vertex in vector
      std::vector<Vertex*>::iterator v=stats.vset.begin();
      while(v!=stats.vset.end())
      {
        if(stats.vertexIsMoveCandidate(*v,c))
        {
          // vertex is a move candidate
          stats.computeVertex(*v,&c);
          //pod.clear();
          //pod.addVertex(stats.cv);
          //pod.addOrigP(stats.cv->pN);
          //pod.addOrigSepDis(sqrt(stats.sep_dis));
          /*
          if(true)
          {
            int rrank;
            tv_iterator ttt;
            stats.findTopN(stats.cv,ttt,rrank);
            // get closest point
            double pC[3];
            c.computePC(stats.cv->cl,stats.cv,pC);
            // set
            pod.addOrigCP(pC);
            pod.addOrigVD(sqrt((*ttt).first));
            pod.addOrigTopN(rrank);
          }
          */
          bool int_flag = false,angle_flag=false;
          if(c.assignNewVertexCoords(s,stats.cv,stats.pH,stats,int_flag,angle_flag))
          {
            /*
            pod.addNewP(stats.cv->pN);
            pod.addNewSepDis(sqrt(stats.cv->getSqSepDist(&c)));
            if(true)
            {
              int rrank;
              tv_iterator ttt;
              stats.findTopN(stats.cv,ttt,rrank);
              // get closest point
              double pC[3];
              c.computePC(stats.cv->cl,stats.cv,pC);
              // set
              pod.addNewCP(pC);
              pod.addNewVD(sqrt((*ttt).first));
              pod.addNewTopN(rrank);
            }
            if(pod.isGood()==true){pod.print();}
            else { pod.printBad();exit(0);}
            */
            // vertex was successfully moved
            stats.updateStatsAndPrint(c,group);
            // update iterator
            v++;
            if(v==stats.vset.end()){break;}
          } 
          else
          {
            // vertex was NOT successfully moved
            //std::vector<Vertex*>::iterator q = stats.detectPunishableVertex(v,int_flag,angle_flag,pod,c);
            std::vector<Vertex*>::iterator q = stats.detectPunishableVertex(v,int_flag,angle_flag);
            v=q;
            if(v==stats.vset.end()){break;}
          }
        }
        else
        {
          // vertex is NOT a move candidate
          // so move on to next vertex in set
          v++;
          if(v==stats.vset.end()){break;}
        }
      } // end for each vertex in set
      if(stats.noSetVerticesMoved())
      {
        // NO vertices were moved from last set
        // build set from vertex with next largest virtual displacemnet
        // (assuming it passes test in identifyMeshRegionToUpdate())
        stats.tvi++;
      }
      else
      {
        // vertices were moved from last set, so reset iterator
        stats.tvi = stats.topN.begin();
      }
    } //end while number of moved vertices less than GROUP_SIZE
    // update log file
    char str[128];
    sprintf(str,"%s%2d%s","Iteration ",group,": Update log files:");
    currtime = recordTime(myfile,currtime,str);

    // write output files and update gain
    stats.writeFilesAndUpdateMaxGain(c,group);
  }

  // write output files
  if(WRITE_MESH_EVERY_ITERATION==false)
  {
    cout << "Build Mesh final..................";
    cout.flush();
    c.buildMeshAfter(1);
    cout << "complete.\n";
    cout.flush();
  }
  if(WRITE_DISTANCES_EVERY_ITERATION==false)
  {
    cout << "Write closest point distances final..";
    cout.flush();
    c.writeDistances(1);
    cout << "complete.\n";
    cout.flush();
  }

  // free allocated memory
  stats.freeAvg();

  ////////// close log files //////////
  c.fileOutit();
  myfile.close();
  cout << "meshmorph complete\n\n";
}
