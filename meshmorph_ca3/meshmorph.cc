// Author: Justin Kinney
// Date: Sep 2008

#include "meshmorph.h"

#include <iostream>
#include <cassert>
#include <cmath>

#include "log.h"
#include "nice.h"
#include "energy.h"
#include "controls.h"
#include "container.h"
#include "refractory.h"
#include "virtual_disp.h"
#include "gain_schedule.h"
#include "vertex_schedule.h"
#include "octree_agent_face.h"
#include "intersecting_faces.h"
#include "octree_visitor_face.h"

using std::cout;
using std::endl;

int main (int argc,char **argv)
{
  cout << "\nMeshmorph 4.0\n"
//  << "  Running on azzuri at Wed Sep 24 12:11:51 2008\n\n"
  << "    The Salk Institute for Biological Studies\n\n"
  << "Meshmorph initializing data...\n";

  // record start time
  time_t begintime = time(NULL);

  // check that assumption of 32 bit int is correct
  assert(checkIntSize());

  // create Controls class
  Controls & cs (Controls::instance()); 

  // parse command line
  cs.parseCommandLine(argc,argv,cs.getUsageMessage());

  // add warning
  if (cs.get_use_edge_reference_length()==true)
  {
    cout << "\n\nmain: WARNING. --ref_orig_edge_length option "
         << "may create faces with high aspect\n"
         << "ratios unless a compensatory "
         << "force is added, i.e. force proportional "
         << "to adjacent face aspect ratio.\n\n";
  }

  // create Log class
  Log & log(Log::instance());

  // create container, objects, vertices, faces, edges, and find adjacencies
  Container & c(Container::instance());
  log.recordTime("scanning file:");

  // update number and size of group
  cs.redefineGroup();

  // write command parameter settings
  log.writeCommandSettings();

  // read frozen vertices
  if (cs.get_frozen_vertices_file().empty()==false){ c.readFrozen(cs.get_frozen_vertices_file().c_str()); }

  // build octree
  cout << "Build octree...................................";
  cout.flush();
  Octree<Face> octree( Vector3r(cs.get_octree_min_x(),
                                cs.get_octree_min_y(),
                                cs.get_octree_min_z()),
                       cs.get_octree_width(),
                       cs.get_max_items_per_leaf(),
                       cs.get_max_octree_depth(),
                       cs.get_min_cell_size());
  Octree_Agent_Face agent;
  // for each object in container
  for (o_it i=c.o.begin();i!=c.o.end();++i)
  {
    // for each face in object
    for (f_it j=i->f.begin();j!=i->f.end();++j)
    {
      octree.insert( *j, agent );
    }
  }
  c.octree = &octree; 
  cout << "complete.\n";
  cout.flush();

  // initialize space data structure
  log.recordTime("Initialize space data structure:");

  // initialize and write data to log files
  log.writeObjectData();
  log.recordTime("Writing data to log files:");

  // identify vertices that lie inside of another object
  Nice::instance().findNonniceVertices();
  log.recordTime("Find nice vertices:");
  log.printNumNonnice(cout);

  // find all face intersections
  Intersecting_Faces::instance().findAllFaceIntersections();
  log.recordTime("Find all face intersections:");
  log.printNumInt(cout);

  // identify the closest point on a mesh to each vertex
  c.findClosestFaceToEachVertex();
  log.recordTime("Get extracellular widths:");
  log.printClosestPtStats(cout);

  // measure and write to file important model attributes
  log.writeFiles(0);
  log.recordTime("Update log files:");

  // create instance of Gain_Schedule
  Gain_Schedule & gs(Gain_Schedule::instance());

  // read vertex move sequence
  if (cs.get_vertex_sequence_file().empty()==false){ Vertex_Schedule::instance().readVertexSequence(cs.get_vertex_sequence_file().c_str()); }

  // main loop
  log.recordTime("Begin main loop:");
  // for each group of iterations
  for (int group=1;group<(cs.get_num_groups()+1);group++)
  {
    Vertex_Schedule & vs(Vertex_Schedule::instance());
    // enforce maximum runtime policy
    if (cs.get_max_runtime()>0.0)
    {
      if ((time(NULL)-begintime) > cs.get_max_runtime()) break;
    }
    gs.initGain();
    // initialize count of vertices moved in group
    vs.setNumMovedVertsGroup(0);
    // initialize group variables
    Virtual_Disp::instance().resetForNewGroup(group);
    Refractory::instance().resetForNewGroup();
    log.groupInit();
    // until GROUP_SIZE vertex moves have been made
    while (vs.getNumMovedVertsGroup()<cs.get_group_size())
    {
      // enforce maximum runtime policy
      if (cs.get_max_runtime()>0.0)
      {
        if ((time(NULL)-begintime) > cs.get_max_runtime()) break;
      }
      vs.collectVerticesToMoveNext(group);
      // for each vertex in collection to move 
      vp_cit v=vs.beginVset();
      while (v!=vs.endVset())
      {
        if (Refractory::instance().vertexIsMoveCandidate(*v))
        {
          // vertex is a move candidate
          vs.calculateMoveLocation(*v);
          if (cs.get_enable_vtrack()==true) log.setDetailedInfoPreMove(*v);
          bool int_flag = false,angle_flag=false,outside_octree=false;
          State & s(State::instance());
          vector3 const * const destination = vs.getVertexDestination();
          vector3 origin(*(*v)->getPos());
          if (s.assignNewVertexCoords(*v,destination,origin,int_flag,angle_flag,outside_octree))
          {
            // record seed actual displacement for steady-state analysis
            if (*v==vs.getSeedVertex())
            {
              int iter = vs.getNumMovedVertsGroup()+static_cast<int>(cs.get_group_size()*(group-1.0));
              Virtual_Disp::instance().addAdToSeedAd(iter,sqrt((*destination-origin).dot((*destination-origin))));
            }
            // vertex was successfully moved
            vs.incrementNumMovedVertsGroup();
            Refractory::instance().updateSuccessfulMove(group,*v);
            cs.updatePrintPeriod(vs.getNumMovedVertsGroup());
            if (!(vs.getNumMovedVertsGroup()%cs.get_print_period()))
            {
              log.writeMoveSummary(group,cs.get_group_size());
            }
            if (cs.get_enable_vtrack()==true)
            {
              log.setDetailedInfoPostMove(*v);
              log.writeDetailedMoveInfo();
            }
            log.updateVertDisplStats(*v,origin);
            // automatically strict face intersection prevention
            if (cs.get_strict_face_intersection_prevention()==false)
            {
              if (Intersecting_Faces::instance().getCountOfIntFaces(false)==0)
              {
                cs.set_strict_face_intersection_prevention();
              }
            }
            gs.updateGainPeriod();
            gs.resetGain();
            if (cs.get_write_mesh_now()>0)
            {
              if (vs.getNumMovedVertsGroup()==cs.get_write_mesh_now())
              {
                c.writeMeshData(group);
              }
            }
            // update iterator
            v++;
            if (v==vs.endVset()){break;}
          } 
          else
          {
            // vertex was NOT successfully moved
            vp_cit q = Refractory::instance().getNextVertex(group,v,int_flag,angle_flag,outside_octree);
            v=q;
            if (v==vs.endVset()){break;}
          }
        }
        else
        {
          // vertex is NOT a move candidate
          // so move on to next vertex in set
          v++;
          if (v==vs.endVset()){break;}
        }
      } // end for each vertex in set
      log.updateMovedVertsFromSet();
      if (vs.noSetVerticesMoved())
      {
        // NO vertices were moved from last set
        // build set from vertex with next largest virtual displacemnet
        // (assuming it passes test in collectVerticesToMoveNext())
        Virtual_Disp::instance().advanceSeedToNextLargestVert();
      }
      else
      {
        // vertices were moved from last set, so reset iterator
        Virtual_Disp::instance().resetSeedToLargestVert();
      }
      // break of executing vertex sequence and sequence has ended
      if (cs.get_vertex_sequence_file().empty()==false) break;
    } //end while number of moved vertices less than GROUP_SIZE
    // check lots of stuff
    c.checkFaces("GROUP");
    c.checkFacesInOctree();
    // DEBUG
//    c.checkClosestFace(group,"GROUP");
    // DEBUG
    Virtual_Disp::instance().validateVirtDispMap();
    Virtual_Disp::instance().validateVirtDispMap2();
    Virtual_Disp::instance().validateVirtDispMapComplement();
    // update log file
    log.recordTime(log.format("Iteration %2d: Update log files:",group));
    log.printVertexSchedulingStats(cout);
    log.printClosestPtSearchStats(cout);
    // write output files
    if (cs.get_write_every_group()==true) log.writeFiles(group);
    // update gain
    gs.updateMaxGain();
  } // end of group

  // write output files
  if (cs.get_write_every_group()==false) log.writeFiles(1);
  if (cs.get_write_ecw_to_file()==true)
  {
    log.writeSepDistances(1);
  }

  log.closeFiles();
  cout << "meshmorph complete\n\n";
}


/** Determine if integers are 32 bit.
 * \return True if integers are 32 bit on this machine; false otherwise.
 */

bool checkIntSize (void)
{
  ///// check that assumption of 32 bit int is correct /////
  if (32==sizeof(int)*8){ return true; }
  else
  {
    std::cout << "Error. Int is not 32 bits, sizeof(int) "
          << sizeof(int) << std::endl;
    return false;
  }
}

/** Determine if two floating-point precision numbers
 * are equivalent in value within epsilon.
 * \param[in] a First number.
 * \param[in] b Second number.
 * \param[in] epsilon The difference between the two input values must be
 * greater than the fraction of the largest input value defined by epsilon.
 * \return 1 if Inputs are different; 0 otherwise.
 */

bool distinguishable (double a,double b,double epsilon)
{
  double c;
  c=a-b;
  if (c<0) c=-c;
  if (a<0) a=-a;
  if (a<1) a=1;
  if (b<0) b=-b;
  if (b<a) return (c>a*epsilon);
  else return (c>b*epsilon);
}

/** Determine if two floating-point precision numbers
 * are equivalent in value within MY_DOUBLE_EPSILON. 
 * \param[in] a First number.
 * \param[in] b Second number.
 * \return 1 if Inputs are different; 0 otherwise.
 */

bool distinguishable (double a,double b)
{
  return distinguishable(a, b, Controls::instance().get_my_double_epsilon());
}

