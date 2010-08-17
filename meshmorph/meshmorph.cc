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

  // read cleft vertices
  if (cs.get_cleft_vertices_file().empty()==false)
  {
    c.readCleft(cs.get_cleft_vertices_file().c_str());
    c.findPerisynaptic();
  }

  // build octree
  cout << "Build octree...................................";
  cout.flush();
  /* calculate octree width as max of each axis width */
  double delx = c.getWorld(1)-c.getWorld(0);
  double dely = c.getWorld(3)-c.getWorld(2);
  double delz = c.getWorld(5)-c.getWorld(4);
  double max = 0.0;
  if ((delz >= delx) && (delz >= dely))
  {
    max = delz;
  }
  else if ((dely >= delx) && (dely >= delz))
  {
    max = dely;
  }
  else
  {
    max = delx;
  }

  cs.recordOctreeSize(max);
  Octree<Face> octree( Vector3r(c.getWorld(0),
                                c.getWorld(2),
                                c.getWorld(4)),
                       max,
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
      //boo--;
      //foo++;
      //if (boo==0)
      //{
      //  boo=yah;
      //  cout << "face " << foo << " inserted into octree.\n";cout.flush();
      //}
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

  // DEBUG
  // identify vertices that lie inside of another object
  if (cs.get_assume_nice_vertices()==false)
  {
    Nice::instance().findNonniceVertices();
    log.recordTime("Find nice vertices:");
    log.printNumNonnice(cout);
  }

  // break if executing recon block wrap and sequence has ended
  if (cs.get_recon_block_wrap()=="")
  {
    // find all face intersections
    if (cs.get_measure_ecw_and_exit()==false)
    {
      Intersecting_Faces::instance().findAllFaceIntersections();
      assert(Intersecting_Faces::instance().intFacesAreSymmetric());
      log.recordTime("Find all face intersections:");
      log.printNumInt(cout);
    }
  }
  // DEBUG

  // identify the closest point on a mesh to each vertex
  c.findClosestFaceToEachVertex();
  log.recordTime("Get extracellular widths:");
  log.printClosestPtStats(cout);

  // count vertex neighbors
  if (cs.get_dual_vertex_ecws())
  {
    c.findVertexNeighbors(0);
  }

  // measure and write to file important model attributes
  if (cs.get_measure_ecw_and_exit()==true)
  {
    log.writeSepDistances(1);
    exit(1);
  }
  log.writeFiles(0);
  log.recordTime("Update log files:");

  // create instance of Gain_Schedule
  Gain_Schedule & gs(Gain_Schedule::instance());

  // read vertex move sequence
  if (cs.get_vertex_sequence_file().empty()==false){ Vertex_Schedule::instance().readVertexSequence(cs.get_vertex_sequence_file().c_str()); }

  // compute regions of orthogonality
  if (cs.get_region_orthogonality_threshold() >= 0.0)
  {
    c.findOrthogonalRegions();
    exit(1);
  }

  // main loop
  log.recordTime("Begin main loop:");
  // for each group of iterations
  for (int group=1;group<(cs.get_num_groups()+1);group++)
  {
    // DEBUG
    //int tunnels = 0;
    //int ref_tunnels = 0;
    //int tunnel_wins = 0;
    //int ref_tunnel_wins = 0;
    //int tunnel_fails = 0;
    //int ref_tunnel_fails = 0;
    //int sheets = 0;
    //int ref_sheets = 0;
    //int count_vsets = 1;
    //int ref_vset = 1;
    // DEBUG
    Vertex_Schedule & vs(Vertex_Schedule::instance());
    // enforce maximum runtime policy
    if (cs.get_max_runtime()>0.0)
    {
      if ((time(NULL)-begintime) > cs.get_max_runtime()) break;
    }
    //gs.initGain();
    //einitialize count of vertices moved in group
    vs.setNumMovedVertsGroup(0);
    // initialize group variables
    if (cs.get_recon_block_wrap()=="")
    {
      Virtual_Disp::instance().resetForNewGroup(group);
    }
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
      // DEBUG
      //cout << "meshmorph::main: "
      //      << "set size = " << vs.getVsetSize() << endl;
      //cout.flush();
      //count_vsets++;
      // DEBUG
      // for each vertex in collection to move 
      // DEBUG
      //cout << "GROUP RAN\n";cout.flush();
      //int mycount = 0;
      // DEBUG
      vp_cit v=vs.beginVset();
      while (v!=vs.endVset())
      {
        //if (mycount++>200) exit(1);
        // DEBUG
        //if ((*v)->isSheet()==false) tunnels++;
        // DEBUG
        //assert(Intersecting_Faces::instance().intFacesAreSymmetric());
        // prohibit vertex move if
        //   (1) not a move candidate OR
        //   (2) sheets are frozen and this vertex is a sheet OR
        //   (3) tunnels are frozen and this vertex is a tunnel
        //if ((!cs.get_freeze_sheets() || !((*v)->isSheet()==true)) &&
        //   (!cs.get_freeze_tunnels() || !((*v)->isSheet()==false)) &&
        //   Refractory::instance().vertexIsMoveCandidate(*v)
        // DEBUG
        //cout << "SET RAN\n";cout.flush();
        // DEBUG
        if (Refractory::instance().vertexIsMoveCandidate(*v))
        {
          if (cs.get_freeze_sheets()) assert( (*v)->isSheet()==false );
          if (cs.get_freeze_tunnels()) assert( (*v)->isSheet()==true );
          // vertex is a move candidate
          vs.calculateMoveLocation(*v);
          if (cs.get_enable_vtrack()==true) log.setDetailedInfoPreMove(*v);
          bool int_flag = false,angle_flag=false,outside_octree=false;
          State & s(State::instance());
          vector3 const * const destination = vs.getVertexDestination();
          vector3 origin(*(*v)->getPos());
          if (s.assignNewVertexCoords(*v,destination,origin,int_flag,angle_flag,outside_octree))
          {
            // DEBUG
            //tunnel_wins++;
            // DEBUG
            // record seed actual displacement for steady-state analysis
            if (*v==vs.getSeedVertex())
            {
              int iter = vs.getNumMovedVertsGroup()+static_cast<int>(cs.get_group_size()*(group-1.0));
              Virtual_Disp::instance().addAdToSeedAd(iter,sqrt((*destination-origin).dot((*destination-origin))));
            }
            // vertex was successfully moved
            vs.incrementNumMovedVertsGroup();
            //Refractory::instance().updateSuccessfulMove(group,*v);
            Refractory::instance().updateSuccessfulMove(*v);
            cs.updatePrintPeriod(vs.getNumMovedVertsGroup());
            if (!(vs.getNumMovedVertsGroup()%cs.get_print_period()))
            {
              log.writeMoveSummary(group,cs.get_group_size());
              // DEBUG
              //double num_vsets = count_vsets - ref_vset;
              //double a = (static_cast<double>(sheets)-
              //            static_cast<double>(ref_sheets))/num_vsets;
              //double b = (static_cast<double>(tunnel_fails)-
              //            static_cast<double>(ref_tunnel_fails))/num_vsets;
              //double cc = (static_cast<double>(tunnel_wins)-
              //            static_cast<double>(ref_tunnel_wins))/num_vsets;
              //double d = (static_cast<double>(tunnels)-
              //            static_cast<double>(ref_tunnels))/num_vsets;
              //cout.precision(2);
              //cout << ", num_vsets = " << num_vsets
              //     << ", sheets per vset = " << a
              //      << ", tunnels per vset = " << d
              //      << ", tunnel_fails per vset = " << b
              //      << ", tunnel_wins per vset = " << cc
              //      << endl;
              //cout.precision(12);
              //ref_tunnels = tunnels;
              //ref_sheets = sheets;
              //ref_tunnel_fails = tunnel_fails;
              //ref_tunnel_wins = tunnel_wins;
              //ref_vset = count_vsets;
              // DEBUG
            }
            if (cs.get_enable_vtrack()==true)
            {
              log.setDetailedInfoPostMove(*v);
              log.writeDetailedMoveInfo(*v,group,vs.getNumMovedVertsGroup());
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
            //gs.updateGainPeriod();
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
            //vp_cit q = Refractory::instance().getNextVertex(group,v,int_flag,angle_flag,outside_octree);
            //tunnel_fails++;
            vp_cit q = Refractory::instance().getNextVertex(v,int_flag,angle_flag,outside_octree);
            v=q;
            if (v==vs.endVset()){break;}
          }
        }
        else
        {
          // vertex is NOT a move candidate
          // so move on to next vertex in set
          // DEBUG
          //if (cs.get_freeze_sheets() && (*v)->isSheet()==true) sheets++;
          //else
          //{
          // if (!Refractory::instance().vertexIsMoveCandidate(*v)) tunnel_fails++;
          // else {cout << "WTF!\n"; exit(1);}
          //}
          // DEBUG
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
      // break if executing vertex sequence and sequence has ended
      if (cs.get_vertex_sequence_file().empty()==false) break;
      // break if executing recon block wrap and sequence has ended
      if (cs.get_recon_block_wrap()!="") break;
    } //end while number of moved vertices less than GROUP_SIZE
    if (cs.get_recon_block_wrap()=="")
    {
      // check lots of stuff
      c.checkFaces("GROUP");
      c.checkFacesInOctree();
      Virtual_Disp::instance().validateVirtDispMap();
      Virtual_Disp::instance().validateVirtDispMap2();
      Virtual_Disp::instance().validateVirtDispMapComplement();
    }
    // update log file
    log.recordTime(log.format("Iteration %2d: Update log files:",group));
    log.printVertexSchedulingStats(cout);
    log.printClosestPtSearchStats(cout);
    // write output files
    if (cs.get_write_every_group()==true) log.writeFiles(group);
    // update gain
    gs.updateMaxGain();
    // update vertex neighbor count
    if (cs.get_dual_vertex_ecws() || cs.get_report_vertex_identity())
    {
      if ( group % cs.get_vertex_neighbor_count_period() == 0 )
      {
        c.findVertexNeighbors(group);
      }
    }
  } // end of group

  // write output files
  if (cs.get_write_every_group()==false)
  {
    log.writeFiles(1);
    //if (cs.get_write_ecw_to_file()==true)
    //{
    //  log.writeSepDistances(1);
    //}
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

