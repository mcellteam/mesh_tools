// Author: Justin Kinney
// Date: Sep 2008

#include "log.h"

#include <sys/time.h>

#include <cassert>
#include <cmath>
#include <iostream>

#include "edge.h"
#include "grid.h"
#include "nice.h"
#include "controls.h"
#include "container.h"
#include "refractory.h"
#include "virtual_disp.h"
#include "gain_schedule.h"
#include "vertex_schedule.h"
#include "intersecting_faces.h"
#include "octree_visitor_measure.h"

using std::cout;
using std::cerr;
using std::endl;
using std::left;
using std::right;

Log * Log::only_one = NULL;

Log & Log::instance()
{
  // Not thread-safe.
  // -- lock mutex
  if (only_one == NULL)
    only_one = new Log();
  // -- unlock mutex
  return *only_one;
}

Log::Log (void)
:Cfile(),Mfile(),num_verts(0),leaves_mean(0),leaves_min(1000000),
  leaves_max(-1000000),face_mean(0),face_min(1000000),face_max(-1000000),
  f_check_mean(0),f_check_min(1000000),f_check_max(-1000000),
  num_sets(0),vps_mean(0),vps_min(1000000),vps_max(-1000000),
  moved_mean(0),moved_min(1000000),moved_max(-1000000),
  N(0),d_min(1E30),d_max(-1E30),md(2,0),
  fs_mean(0),fs_min(1000000),fs_max(-1000000),
  ps_mean(0),ps_min(1000000),ps_max(-1000000),
  fs_face_change_mean(0),
  fs_face_change_min(1000000),
  fs_face_change_max(-1000000),
  ps_face_change_mean(0),
  ps_face_change_min(1000000),
  ps_face_change_max(-1000000),
//  fs_mean_face_change_is_adj(0),ps_mean_face_change_is_adj(0),
  topN_val(2,0),vd_val(2,0),sepdis(2,0),p_orig(),
  //topN_val(2,0),vd_val(2,0),sepdis(2,0),p_orig(3,-1.0),
  p_new(NULL),cp_orig(),cp_new(),trouble(0),
  num_faces(0),bpf_mean(0),bpf_min(1000000),bpf_max(-1000000),
  num_boxes(0),num_empty_boxes(0),
  fpb_mean(0),fpb_min(1000000),fpb_max(-1000000),
  t1(0),tim(),currtime()
{
  setTime(time(NULL));
  // track elapsed time
  gettimeofday(&tim,NULL);
  t1=tim.tv_sec+(tim.tv_usec/1000000.0);
  openMainFile();
}

/** Clear members of this class.
*/

void Log::clearVals (void)
{
  N=0;		  // count of moved vertices
  md[0]=md[1]=0;  // mean displacement of N moved vertices
  d_min=1E30;	  // min displacement of N moved vertices
  d_max=-1E30;	  // max displacement of N moved vertices
  fs_mean = 0.0; // mean number of vertices undergoing full closest point search
  fs_min = 1000000; // min number of vertices undergoing full closest point search
  fs_max = -1000000; // max number of vertices undergoing full closest point search
  ps_mean = 0.0; // mean number of vertices undergoing partial closest point search
  ps_min = 1E30; // min number of vertices undergoing partial closest point search
  ps_max = -1E30; // max number of vertices undergoing partial closest point search
  vps_mean = 0.0; // mean number of vertices per set
  vps_min = 1000000; // min number of vertices per set
  vps_max = -1000000; // max number of vertices per set
  moved_mean = 0.0; // mean number of moved vertices per set
  moved_min = 1000000; // min number of moved vertices per set
  moved_max = -1000000; // max number of moved vertices per set
  fs_face_change_mean = 0.0; // fraction of vertices whose closest faces changed
  ps_face_change_mean = 0.0; // fraction of vertices whose closest faces changed
  fs_face_change_min = 1E30; // fraction of vertices whose closest faces changed
  ps_face_change_min = 1E30; // fraction of vertices whose closest faces changed
  fs_face_change_max = -1E30; // fraction of vertices whose closest faces changed
  ps_face_change_max = -1E30; // fraction of vertices whose closest faces changed
//  fs_mean_face_change_is_adj = 0.0; // fraction of new closest faces that are adjacent to moved vertex
//  ps_mean_face_change_is_adj = 0.0; // fraction of new closest faces that are adjacent to moved vertex
}

/** Initialize this class for start of group.
*/

void Log::groupInit (void)
{
  // intialize vertex displacement stats
  clearVals();
  // track elapsed time
  gettimeofday(&tim,NULL);
  t1=tim.tv_sec+(tim.tv_usec/1000000.0);
}

/** Write to file the distribution of number of moves per vertex
 * since start of last initialization. 
 * \param[in] group Current group number (used as filename suffix).
 */

void Log::writeVertMoveDistribution (int const & group)
{
  Controls & cs(Controls::instance());
  int zero=0;
  // open log file
//  char file[cs.get_max_filename_size()];
//  sprintf(file,"%s%s.%d",cs.get_output_data_dir().c_str()
//                        ,cs.get_vertex_selection_file().c_str()
//                        ,group);
  std::string file = format("%s%s.%d",
                            cs.get_output_data_dir().c_str(),
                            cs.get_vertex_selection_file().c_str(),
                            group);
  //std::ofstream this_file (file);
  std::ofstream this_file (file.c_str());
  Container & c(Container::instance());
  // for each object in container
  for (o_cit i=c.o.begin();i!=c.o.end();++i)
  {
    // for each vertex in object
    for (v_cit j=i->v.begin();j!=i->v.end();++j)
    {
      Refractory & r(Refractory::instance());
      // if vertex has moved since the group began 
      Vertex * v = const_cast<Vertex*>(&(*j));
      if (r.vertexMovesAreRecorded(v)==true)
      {
        this_file << r.numRecordedVertexMoves(v) << endl;
      }
      else
      {
        this_file << zero << endl;
      }
    }
  }
  this_file.close();
}

/** Determine if extracellular distance to object is to be stored.
 * \param[in] name Identity of neighboring object.
 * \param[in] dd Extracellular distance to object.
 * \param[out] hausdorff_1 Map of neighboring objects and associated extracellular widths.
 */

void Log::updateHauss_1 (std::string const & name,double const & dd,map_s_d & hausdorff_1)
{
  // if object found in map
  sd_it w = hausdorff_1.find(name);
  if (w!=hausdorff_1.end())
  {
    // if new extracellular width is larger than stored value
    if (dd>(*w).second){hausdorff_1[name]=dd;}
  }
  else
  {
    // add distance to map
    hausdorff_1[name]=dd;
  }
}

/** Determine if extracellular distance to object is to be stored.
 * \param[in] name Identity of neighboring object.
 * \param[in] i Object of interest.
 * \param[in] dd Extracellular distance to object.
 * \param[out] hausdorff_2 Map of neighboring objects and associated extracellular widths.
 */

void Log::updateHauss_2 (std::string const & name,o_it const & i,double const & dd,map_s_d & hausdorff_2)
{
  // if object found in map
  std::string ss = name+"_"+i->getName();
  sd_it w = hausdorff_2.find(ss);
  std::string t = i->getName()+"_"+name;
  sd_it x = hausdorff_2.find(t);
  if (w!=hausdorff_2.end())
  {
    // if new extracellular width is larger than stored value
    if (dd>(*w).second){hausdorff_2[ss]=dd;}
  }
  else if (x!=hausdorff_2.end())
  {
    // if new extracellular width is larger than stored value
    if (dd>(*x).second){hausdorff_2[t]=dd;}
  }
  else
  {
    // add distance to map
    hausdorff_2[ss]=dd;
  }
}

/** Determine if extracellular distance to object is to be stored.
 * \param[in] name Identity of neighboring object.
 * \param[in] dd Extracellular distance to object.
 * \param[out] hausdorff_1_noself Map of neighboring objects and associated extracellular widths.
 */

void Log::updateHauss_1_noself (std::string const & name,double const & dd,map_s_d & hausdorff_1_noself)
{
  // if object found in map
  sd_it w = hausdorff_1_noself.find(name);
  if (w!=hausdorff_1_noself.end())
  {
    // if new extracellular width is larger than stored value
    if (dd>(*w).second){hausdorff_1_noself[name]=dd;}
  }
  else
  {
    // add distance to map
    hausdorff_1_noself[name]=dd;
  }
}

/** Determine if extracellular distance to object is to be stored.
 * \param[in] name Identity of neighboring object.
 * \param[in] i Object of interest.
 * \param[in] dd Extracellular distance to object.
 * \param[out] hausdorff_2_noself Map of neighboring objects and associated extracellular widths.
 */

void Log::updateHauss_2_noself (std::string const & name,o_it const & i,double const & dd,map_s_d & hausdorff_2_noself)
{
  // if object found in map
  std::string ss = name+"_"+i->getName();
  sd_it w = hausdorff_2_noself.find(ss);
  std::string t = i->getName()+"_"+name;
  sd_it x = hausdorff_2_noself.find(t);
  if (w!=hausdorff_2_noself.end())
  {
    // if new extracellular width is larger than stored value
    if (dd>(*w).second){hausdorff_2_noself[ss]=dd;}
  }
  else if (x!=hausdorff_2_noself.end())
  {
    // if new extracellular width is larger than stored value
    if (dd>(*x).second){hausdorff_2_noself[t]=dd;}
  }
  else
  {
    // add distance to map
    hausdorff_2_noself[ss]=dd;
  }
}


void Log::openOrDie (std::ofstream * const handle,std::string str,const int & group)
{
  Controls & cs(Controls::instance());
  std::string file;
  if (cs.get_append_group_number()==true)
  {
    file = format("%s%s.%s.%i",
                  cs.get_output_data_dir().c_str(),
                  cs.get_sep_log_file().c_str(),str.c_str(),group);
  }
  else
  {
    file = format("%s%s.%s",
                  cs.get_output_data_dir().c_str(),
                  cs.get_sep_log_file().c_str(),str.c_str());
  }
  handle->open (file.c_str());
  if (handle->is_open()==false)
  {
    cerr << "\nCould not open output file "
          << file.c_str() << endl;
    assert(handle->is_open()==true);
    exit(1);
  }
  handle->precision(4);
}


/** Calculate and write to file the extracellular width (or some metric)
 * of each vertex in model using the stored closest face.
 * \param[in] group Current group number (used as filename suffix).
 */

void Log::writeSepDistances (int const & group)
{
  Controls & cs(Controls::instance());
  cout << "Writing extracellular widths to file...........";
  cout.flush();	
  // open output files
  std::ofstream out, out_NOCP, out_haus1, out_haus2, out_haus1_noself, out_haus2_noself;
  openOrDie (&out,"",group);
  openOrDie (&out_NOCP,"NOCP",group);
  openOrDie (&out_haus1,"1_sided_hausdorff",group);
  openOrDie (&out_haus2,"2_sided_hausdorff",group);
  openOrDie (&out_haus1_noself,"1_sided_hausdorff_noself",group);
  openOrDie (&out_haus2_noself,"2_sided_hausdorff_noself",group);

  out_NOCP << "x_coordinate y_coordinate z_coordinate state_value "
        << "x_normal y_normal z_normal\n";
  // declare variables
  int total_vertices=0,vertices_used=0,vertices_used_noself=0;
  map_s_d hausdorff_2_noself; // create map: string->double
  map_s_d hausdorff_2;        // create map: string->double
  map_s_d hausdorff_1_noself; // create map: string->double
  map_s_d hausdorff_1;        // create map: string->double
  vec_d storage_1_noself;
  vec_d storage_2_noself;
  vec_d storage_1;
  vec_d storage_2;
  // area of equilateral triangle of side ecw_sampling_length 
  double ecw_sampling_area = cs.get_ecw_sampling_length()*0.86603/2.0;
  // for each object
  Container & c(Container::instance());
  for (o_it i=c.o.begin();i!=c.o.end();++i)
  {
    hausdorff_1.clear();
    hausdorff_1_noself.clear();
    // for each vertex in object
    for (v_it j=i->v.begin();j!=i->v.end();++j)
    {
      total_vertices++;
      // if vertex has a closest face
      if (j->getClosestFace()!=NULL)
      {
        vertices_used++;
        std::string name = j->getClosestFace()->getVertex(0)->getObject()->getName();
        // get closest point
        vector3 closest_point;
        double sqd;
        c.findClosestPtInFaceToLocation(*j->getPos(),j->getClosestFace(),closest_point,sqd);
        double dd = sqrt(sqd);
        // compute extracellular width
        if (Nice::instance().vertexIsNice(&(*j))==false){ dd=-dd;}
        // print extracellular width
        out << dd << endl;
        updateHauss_1(name,dd,hausdorff_1);
        updateHauss_2(name,i,dd,hausdorff_1);
        // if object different from self
        if (i->getName()!=name)
        {
          vertices_used_noself++;
          updateHauss_1_noself(name,dd,hausdorff_1_noself);
          updateHauss_2_noself(name,i,dd,hausdorff_2_noself);
        }
      }
      else
      {
        out_NOCP << j->getCoord(0) << " "
              << j->getCoord(1) << " "
              << j->getCoord(2) << " 1 0 0 1\n";
      }
    }
    // for each face in object
    for(f_it j=i->f.begin();j!=i->f.end();j++)
    {
      // compute face area
      double face_area = j->computeArea();
      // compute number of tiles for face
      double normalized_area = face_area/ecw_sampling_area;
      // need number of tiles on face to be a perfect square
      int n = static_cast<int>(ceil(sqrt( normalized_area )));
      if (n<1) n=1;
      int m=n*n; // number of tiles
      // if number of tiles is greater than 1 then subdivide face
      if(m > 1.0)
      {
        // create grid
        Grid g(&(*j));
        // for each tile on face
        for(int k=0;k<m;k++)
        {
          // compute barycenter
          vector3 pt;
          g.computeBarycenter(pt,k,n,&(*j));
          // find closest point to barycenter
          vector3 p;
          double sqd;
          Face *ncl = NULL;
          if (c.findClosestPtToBarycenter(pt,&(*j),p,sqd,ncl)==true)
          {
            std::string name = ncl->getVertex(0)->getObject()->getName();
            double dd = sqrt(sqd);
            //write to out and process as hausdorffs
            out << dd << endl;
            updateHauss_1(name,dd,hausdorff_1);
            updateHauss_2(name,i,dd,hausdorff_1);
            // if object different from self
            if (j->getVertex(0)->getObject()->getName()!=name)
            {
              updateHauss_1_noself(name,dd,hausdorff_1_noself);
              updateHauss_2_noself(name,i,dd,hausdorff_2_noself);
            }
          }
        }
      }
    }
    for (sd_it w=hausdorff_1.begin();w!=hausdorff_1.end();++w)
    {
      storage_1.push_back((*w).second);
    }
    for (sd_it w=hausdorff_2.begin();w!=hausdorff_2.end();++w)
    {
      storage_2.push_back((*w).second);
    }
    for (sd_it w=hausdorff_1_noself.begin();w!=hausdorff_1_noself.end();++w)
    {
      storage_1_noself.push_back((*w).second);
    }
    for (sd_it w=hausdorff_2_noself.begin();w!=hausdorff_2_noself.end();++w)
    {
      storage_2_noself.push_back((*w).second);
    }
  }
  // haus1
  out_haus1 << "# One-sided Hausdorff distances\n"
        << "# total_vertices " << total_vertices
        << ", vertices_used " << vertices_used << endl;
  for (d_it i=storage_1.begin();i!=storage_1.end();++i)
  {
    out_haus1 << *i << endl;
  }
  // haus1
  // haus2
  out_haus2 << "# Two-sided Hausdorff distances\n"
        << "# total_vertices " << total_vertices
        << ", vertices_used " << vertices_used << endl;
  for (d_it i=storage_2.begin();i!=storage_2.end();++i)
  {
    out_haus2 << *i << endl;
  }
  // haus2
  // haus1_noself
  out_haus1_noself << "# One-sided Hausdorff distances - no self\n"
        << "# total_vertices " << total_vertices
        << ", vertices_used " << vertices_used_noself << endl;
  for (d_it i=storage_1_noself.begin();i!=storage_1_noself.end();++i)
  {
    out_haus1_noself << *i << endl;
  }
  // haus1_noself
  // haus2_noself
  out_haus2_noself << "# Two-sided Hausdorff distances - no self\n"
        << "# total_vertices " << total_vertices
        << ", vertices_used " << vertices_used_noself << endl;
  for (d_it i=storage_2_noself.begin();i!=storage_2_noself.end();++i)
  {
    out_haus2_noself << *i << endl;
  }
  // haus2_noself
  out.close();
  out_NOCP.close();
  out_haus1.close();
  out_haus2.close();
  out_haus1_noself.close();
  out_haus2_noself.close();
  cout << "complete.\n";
  cout.flush();
}

/** Initialize progress file for vertex move campaign.
*/

void Log::statusFileInit (void)
{
  Controls & cs(Controls::instance());
//  char file[cs.get_max_filename_size()];
//  sprintf(file,"%s%s",cs.get_output_data_dir().c_str(),cs.get_cont_log_file().c_str());
  std::string file = format("%s%s",
                            cs.get_output_data_dir().c_str(),
                            cs.get_cont_log_file().c_str());
  Cfile.open(file.c_str());
  if (Cfile.is_open()==false)
  {
    cout << "\nLog::statusFileInit: "
          << "Error. Unable to open file = "
          << file << endl;
    assert(Cfile.is_open()==true);
    exit(1);
  }
  Cfile << "\nLEGEND_______________________\n"
        << "'Iteration'\n"
        << "0th group is the original data.\n"
        << "For 1st group vertices have moved once or not at all.\n"
        << "--\n"
        << "'Nonnice'\n"
        << "Nonnice vertices are inside of another object, "
        << "possibly the parent object.\n"
        << "--\n"
        << "'Self-Nonnice'\n"
        << "Self-nonnice vertices are inside of their parent object.\n"
        << "--\n"
        << "'Int. Faces'\n"
        << "Total number of pairs of intersecting faces in all objects.\n"
        << "--\n"
        << "'Self-Int. Faces'\n"
        << "Number of pairs of self-intersecting faces, "
        << "i.e. faces share a common parent object.\n"
        << "--\n"
        << "'Force'\n"
        << "Cumulative force on all vertices in all objects.\n"
        << "--\n"
        << "'Energy'\n"
        << "Cumulative potential energy in all vertices in all objects.\n"
        << "--\n"
        << "'Mean Disp.'\n"
        << "Mean displacement of all vertices in all objects "
        << "chosen to move during the group.\n"
        << "--\n"
        << "'Min. Edge Angle (rad.)'\n"
        << "The minimum edge angle of all edges in all objects.\n"
        << "Minimum edge angle is the angle between the two faces "
        << "that share the edge.\n"
        << "___________________________________________\n\n\n\n";

  Cfile << "             S\n";
  Cfile << "             .\n";
  Cfile << "             n            S\n";
  Cfile << "             o            .\n";
  Cfile << "G            n            F\n";
  Cfile << "r            n            a\n";
  Cfile << "o            i            c\n";
  Cfile << "u            c            e\n";
  Cfile << "p   Nonnice  e   I.Faces  s   ";
  Cfile.width(10);
  Cfile << left << "N";
  Cfile.width(10);
  Cfile << left << "Max Gain";
  Cfile.width(11);
  Cfile << left << "Disp_Mean";
  Cfile.width(15);
  Cfile << right << "Disp_Min";
  Cfile.width(1);
  Cfile << ":";
  Cfile.width(15);
  Cfile << left << "Disp_Max";
  Cfile.width(19);
  Cfile << left << "MinEdgeAngle(deg)";
  Cfile.width(13);
  Cfile << left << "Duration(sec)" << endl;
  Cfile.flush();
}

/** Close output streams of this class.
*/

void Log::closeFiles (void)
{
  Cfile.close();
  Mfile.close();
}

/** Write to file the progress of vertex move campaign.
 * \param[in] group Current group number.
 * \param[in] not_first_group True if group is NOT the first group;
 * false otherwise.
 * \param[in] elapsed_time Elapsed time since last initialization.
 */

void Log::updateFile (int const & group,bool verts_moved,double const & elapsed_time)
{
//  if (group==0)
//  {
//    cout << "Update log files...............................";
//  }
//  else
//  {
    cout << "Iteration " << group << ": ";
    cout << "Update log files...............................";
//  }
  cout.flush();
  Cfile.precision(10);
  Cfile.width(4);
  Cfile << left << group;
  Cfile.width(9);
  Cfile << left << Nice::instance().getNonniceCount(false); // all nonnice
  Cfile.width(4);
  Cfile << left << Nice::instance().getNonniceCount(true);  // nonnince to self
  Cfile.width(9);
  Cfile << left << Intersecting_Faces::instance().getCountOfIntFaces(false);
  Cfile.width(4);
  Cfile << left << Intersecting_Faces::instance().getCountOfIntFaces(true);
  Cfile.precision(8);
  Cfile.width(10);
  Cfile << left << N;
  //Cfile.width(10);
  //Cfile << left << format("%.3g",Gain_Schedule::instance().getMaxGain());
  if (verts_moved==false)
  {
    Cfile.width(10);
    Cfile << left << "NA";
    Cfile.width(11);
    Cfile << left << "NA";
    Cfile.width(15);
    Cfile << right << "NA";
    Cfile.width(1);
    Cfile << ":";
    Cfile.width(15);
    Cfile << left << "NA";
  }
  else
  {
    Cfile.width(10);
    Cfile << left << format("%.3g",Gain_Schedule::instance().getMaxGain());
    Cfile.width(11);
    Cfile << left << md[1];
    Cfile.width(15);
    Cfile << right << format("%.6g",d_min);
    Cfile.width(1);
    Cfile << ":";
    Cfile.width(15);
    Cfile << left << format("%.6g",d_max);
  }
  Cfile.precision(10);
  Cfile.width(19);
  Cfile << left << Container::instance().getMinEdgeAngle()
                   *180.0/Controls::instance().get_pi();
  Cfile.width(13);
  Cfile.precision(4);
  Cfile << left << elapsed_time << endl;
  Cfile.flush();
  cout << "complete.\n";
  cout.flush();
}

/** Write to file the object list and initialize
 * the vertex move campaing progress file.
 */

void Log::writeObjectData ()
{
  cout << "Writing object data to log files...............";
  cout.flush();
  writeObjectList();
  statusFileInit();
  cout << "complete.\n";
  cout.flush();
}

/** Write to file a list of objects in model.
*/

void Log::writeObjectList (void)
{
  Container & c(Container::instance());
  Controls & cs(Controls::instance());
  //char file[cs.get_max_filename_size()];
  //sprintf(file,"%s%s",cs.get_output_data_dir().c_str(),cs.get_object_list_file().c_str());
  std::string file = format("%s%s",
                            cs.get_output_data_dir().c_str(),
                            cs.get_object_list_file().c_str());
  // open file
  std::ofstream Ofile(file.c_str());
  if (Ofile.is_open()==false)
  {
    cout << "\nLog::writeObjectList: "
          << "Error. Unable to open file = "
          << file << endl;
    assert(Ofile.is_open()==true);
    exit(1);
  }
  // add stuff
  Ofile << "Input data directory = " << cs.get_input_data_dir() << "\n"
        << "Total number of input files = "
        << c.getFileCount() << "\n"
        << "Total number of (objects, vertices, faces, edges) = (" 
        << c.o.size() << ","
        << c.getVertexCount() << ","
        << c.getFaceCount() << ","
        << c.getEdgeCount() << ")\n\n";
  // write world bounds
  Container::instance().writeSummary(Ofile);
  // measure and report octree dimensions
//  int max_depth=0,num_leaves=0;
  Octree_Visitor_Measure visitor;
  Container::instance().octree->visit( visitor );
  double octree_max_x = cs.get_octree_min_x()+cs.get_octree_width();
  double octree_max_y = cs.get_octree_min_y()+cs.get_octree_width();
  double octree_max_z = cs.get_octree_min_z()+cs.get_octree_width();
  Ofile << "Octree size:\n"
        << "   Specified octree envelope ["
        << cs.get_octree_min_x() << " "
        << cs.get_octree_min_y() << " "
        << cs.get_octree_min_z() << "] ["
        << octree_max_x << " "
        << octree_max_y << " "
        << octree_max_z << "]\n"
        << "Octree settings:\n"
        << "   Max faces per leaf = "
        << cs.get_max_items_per_leaf() << endl
        << "   Max allowed octree depth = "
        << cs.get_max_octree_depth() << endl
        << "   Min octree cell size = "
        << cs.get_min_cell_size() << endl
        << "Octree measurements:\n"
        << "   Number of leaves in octree = "
        << visitor.get_num_leaves() << endl
        << "   Octree depth = "
        << visitor.get_max_depth() << endl;
  // report objects
  Ofile << endl << endl;
  Ofile.width(15);
  Ofile << left << "Object name";	
  Ofile.width(15);
  Ofile << left << "#vertices";
  Ofile.width(15);
  Ofile << left << "#faces";
  Ofile.width(15);
  Ofile << left << "#edges" << endl;
  // for each object, write name and index
  for (o_cit i=c.o.begin();i!=c.o.end();++i)
  {
    Ofile.width(15);
    Ofile << left << i->getName();
    Ofile.width(15);
    Ofile << left << i->v.size();
    Ofile.width(15);
    Ofile << left << i->f.size();
    Ofile.width(15);
    Ofile << left << i->e.size() << endl;
  }
  Ofile.close();
}

/** Update statistics on actual vertex displacements.
 * \param[in] vertex_displ Most recent vertex displacement.
 */

//void Log::updateVertDisplStats (double const & vertex_displ)
void Log::updateVertDisplStats (Vertex * const v,vector3 const & old_pos)
{
//  cout << "\nLog::updateVertDisplStats: p_orig\n";
//  cout.flush();
//  old_pos.print(cout);
//  cout << endl;
//  cout.flush();
//  cout << "\nLog::updateVertDisplStats: p_new\n";
//  cout.flush();
//  v->print(cout);
//  cout.flush();
//  //vector3 diff((*p_new)-(*v->getPos());
//  cout << "this ran\n";
//  cout.flush();
  vector3 diff = *v->getPos()-old_pos;
//  cout << "\nLog::updateVertDisplStats: "
//        << "diff ";
//  cout.flush();
//  diff.print(cout);
//  cout.flush();
  double vertex_displ = sqrt(diff.dot(diff));
  N++;
  md[0]=md[1];
  md[1]=md[0]*(N-1)/N+vertex_displ/N;
  if (vertex_displ<d_min){d_min=vertex_displ;}
  if (vertex_displ>d_max){d_max=vertex_displ;}
}

/** Write to file the current position of model and diagnostic data.
 * \param[in] group Current group number.
 */

void Log::writeFiles (int const & group)
{
  Controls & cs(Controls::instance());

  cout << "Iteration " << group << ": ";
  cout << "Write mesh to file.............................";
  cout.flush();
  Container::instance().writeMeshData(group);
  cout << "complete.\n";
  cout.flush();

  if (cs.get_write_refracted_vertices_to_file()==true)
  {
    cout << "Iteration " << group << ": ";
    cout << "Write refracted vertices.......................";
    cout.flush();
    writeRefracted(group);
    cout << "complete.\n";
    cout.flush();
  }
  if (cs.get_write_intersected_faces_to_file()==true)
  {
    cout << "Iteration " << group << ": ";
    cout.flush();
    writeIntersected(group);
  }
  if (cs.get_write_nonnice_vertices_to_file()==true)
  {
    cout << "Iteration " << group << ": ";
    cout << "Write nonnice vertices.........................";
    cout.flush();
    writeNonnice(group);
    cout << "complete.\n";
    cout.flush();
  }
  if (cs.get_write_vertex_move_histogram()==true) 
  {
    writeVertMoveDistribution(group);
  }

  gettimeofday(&tim,NULL);
  double t2=tim.tv_sec+(tim.tv_usec/1000000.0);
  updateFile(group,Vertex_Schedule::instance().getNumMovedVertsGroup()>0,t2-t1);
//  if (cs.get_write_ecw_to_file()==true)
//  {
//    cout << "Iteration " << group << ": ";
//    cout.flush();
//    if (group>0) writeSepDistances(group);
//  }
  if (group>0)
  {
    // write seed virtual displacements to file
    writeVD(group);
    writeAD(group);
  }
}

/** Write to file list of seed vertex virtual displacements.
 * \param[in] group Current group number (used as filename suffix).
 */

void Log::writeAD (int const & group)
{
  //std::string file;
  Controls & cs(Controls::instance());
  //char file[cs.get_max_filename_size()];
  std::string file;
  // create output filename
  if (cs.get_append_group_number()==true)
  {
    //sprintf(file,"%sseed_vd.dat.%i",cs.get_output_data_dir().c_str(),group);
    file = format("%sseed_ad.dat.%i",cs.get_output_data_dir().c_str(),group);
  }
  else
  {
    //sprintf(file,"%sseed_vd.dat",cs.get_output_data_dir().c_str());
    file = format("%sseed_ad.dat",cs.get_output_data_dir().c_str());
  }
  // open output file
  std::ofstream newfile (file.c_str(),std::ios::app);
  if (newfile.is_open()==false)
  {
    cout << "\nLog::writeVD: "
          << "Error. Unable to open file = "
          << file << endl;
    assert(newfile.is_open()==true);
    exit(1);
  }
  newfile.precision(4);
  Virtual_Disp & vd(Virtual_Disp::instance());
  //double scale = 1.0/sqrt(cs.get_search_radius_sq());
  // for each stored seed virtual displacement in history
  for (s_cit j=vd.beginSeedActDisp();j!=vd.endSeedActDisp();++j)
  {
    //double norm = sqrt(*j)*scale;
    newfile << *j << endl;
  }
  newfile.close();
}

/** Write to file list of seed vertex virtual displacements.
 * \param[in] group Current group number (used as filename suffix).
 */

void Log::writeVD (int const & group)
{
  //std::string file;
  Controls & cs(Controls::instance());
  //char file[cs.get_max_filename_size()];
  std::string file;
  // create output filename
  if (cs.get_append_group_number()==true)
  {
    //sprintf(file,"%sseed_vd.dat.%i",cs.get_output_data_dir().c_str(),group);
    file = format("%sseed_vd.dat.%i",cs.get_output_data_dir().c_str(),group);
  }
  else
  {
    //sprintf(file,"%sseed_vd.dat",cs.get_output_data_dir().c_str());
    file = format("%sseed_vd.dat",cs.get_output_data_dir().c_str());
  }
  // open output file
  std::ofstream newfile (file.c_str(),std::ios::app);
  if (newfile.is_open()==false)
  {
    cout << "\nLog::writeVD: "
          << "Error. Unable to open file = "
          << file << endl;
    assert(newfile.is_open()==true);
    exit(1);
  }
  newfile.precision(4);
  Virtual_Disp & vd(Virtual_Disp::instance());
  //double scale = 1.0/sqrt(cs.get_search_radius_sq());
  // for each stored seed virtual displacement in history
  for (s_cit j=vd.beginSeedVirtDisp();j!=vd.endSeedVirtDisp();++j)
  {
    //double norm = sqrt(*j)*scale;
    newfile << *j << endl;
  }
  newfile.close();
}

/** Write to file list of nonnice vertices in model.
 * \param[in] group Current group number (used as filename suffix).
 */

void Log::writeNonnice (int const & group)
{
  Controls & cs(Controls::instance());
  //char file[cs.get_max_filename_size()];
  std::string file;
  // create output filename
  if (cs.get_append_group_number()==true)
  {
    //sprintf(file,"%s%s.%i",cs.get_output_data_dir().c_str(),cs.get_nonnice_file().c_str(),group);
    file = format("%s%s.%i",cs.get_output_data_dir().c_str(),cs.get_nonnice_file().c_str(),group);
  }
  else
  {
    //sprintf(file,"%s%s",cs.get_output_data_dir().c_str(),cs.get_nonnice_file().c_str());
    file = format("%s%s",cs.get_output_data_dir().c_str(),cs.get_nonnice_file().c_str());
  }
  // open output file
  std::ofstream newfile (file.c_str(),std::ios::out);
  if (newfile.is_open()==false)
  {
    cout << "\nLog::writeNonnice: "
          << "Error. Unable to open file = "
          << file << endl;
    assert(newfile.is_open()==true);
    exit(1);
  }
  // create vector set
  v_set vs;
  // for each nonnice hashtable element in object
  // Note sequence of elements in hash_map is not necessarily repeatable.
  for (vhm_cit j=Nice::instance().beginNice();j!=Nice::instance().endNice();++j)
  {
    // add vertex* to set
    vs.insert((*j).first);
  }
  // vector set was collected
  // if cp
  if (cs.get_format_nonnice_vertices()=="dreamm")
  {
    // print header
    newfile << "x_coordinate y_coordinate z_coordinate "
          << "state_value x_normal y_normal z_normal\n";
    // for each vertex in set
    newfile.precision(12);
    for (vs_it i=vs.begin();i!=vs.end();++i)
    {
      // print vertex
      (*i)->printCP(newfile);
    }
  }
  else
  { // else detail
    newfile.precision(12);
    // for each vertex in set
    for (vs_it i=vs.begin();i!=vs.end();++i)
    {
      // print vertex
      (*i)->print(newfile);
    }
  }
  newfile.close();
}

/** Write to file all intersected faces in model.
 * \param[in] group Current group number (used as filename suffix).
 */

void Log::writeIntersected (int const & group)
{
  cout << "Writing intersecting faces to file.............";
  cout.flush();	
  Controls & cs(Controls::instance());
  //char file[cs.get_max_filename_size()];
  std::string file;
  // create output filename
  if (cs.get_append_group_number()==true)
  {
    //sprintf(file,"%s%s.%i",cs.get_output_data_dir().c_str(),cs.get_intersected_file().c_str(),group);
    file = format("%s%s.%i",cs.get_output_data_dir().c_str(),cs.get_intersected_file().c_str(),group);
  }
  else
  {
    //sprintf(file,"%s%s",cs.get_output_data_dir().c_str(),cs.get_intersected_file().c_str());
    file = format("%s%s",cs.get_output_data_dir().c_str(),cs.get_intersected_file().c_str());
  }
  Intersecting_Faces & i_s(Intersecting_Faces::instance());
  // open output file
  std::ofstream newfile (file.c_str(),std::ios::out);
  if (newfile.is_open()==false)
  {
    cout << "\nLog::writeIntersected: "
          << "Error. Unable to open file = "
          << file << endl;
    assert(newfile.is_open()==true);
    exit(1);
  }
  // if cp
  if (cs.get_format_intersected_faces()=="dreamm")
  {
    // create vector set
    v_set vs;
    // for each intersected hashtable element in object
    // Note sequence of elements in hash_map is not necessarily repeatable.
    for (htff_cit j=i_s.begin();j!=i_s.end();++j)
    {
      // grab Face *key
      Face *ff=(*j).first;
      // for each vertex of Face* key
      for (int k=0;k<3;++k)
      {
        // add vertex* to set
        vs.insert(ff->getVertex(k));
      }
      // for each Face* in vector
      const vec_fp *fv=&(*j).second;
      for (fp_cit k=fv->begin();k!=fv->end();++k)
      {
        // for each vertex of Face*
        for (int nn=0;nn<3;++nn)
        {
          // add vertex* to set
          vs.insert((*k)->getVertex(nn));
        }
      }
    }
    // vector set was collected
    // print header
    newfile << "x_coordinate y_coordinate z_coordinate "
          << "state_value x_normal y_normal z_normal\n";
    // for each vertex in set
    newfile.precision(12);
    for (vs_it i=vs.begin();i!=vs.end();++i)
    {
      // print vertex
      (*i)->printCP(newfile);
    }
  }
  else
  { // else detail
    // create face set
    f_set ffs;
    // for each intersected hashtable element
    // Note sequence of elements in hash_map is not necessarily repeatable.
    for (htff_cit j=i_s.begin();j!=i_s.end();++j)
    {
      // add Face* to set
      ffs.insert((*j).first);
      // for each Face* in vector
      const vec_fp *fv=&(*j).second;
      for (fp_cit k=fv->begin();k!=fv->end();++k)
      {
        // add Face* to set
        ffs.insert(*k);
      }
    }
    // face set was collected
    newfile.precision(12);
    // for each face in set
    for (fs_it i=ffs.begin();i!=ffs.end();++i)
    {
      // print face
      (*i)->print(newfile);
    }
  }
  newfile.close();
  cout << "complete.\n";
  cout.flush();
}

/** Write to file all vertices refracted since last initialization.
 * \param[in] group Current group number (used as filename suffix).
 */

void Log::writeRefracted (int const & group)
{
  Controls & cs(Controls::instance());
  //char file[cs.get_max_filename_size()];
  std::string file;
  // create output filename
  if (cs.get_append_group_number())
  {
    //sprintf(file,"%s%s.%i",cs.get_output_data_dir().c_str(),cs.get_refracted_file().c_str(),group);
    file = format("%s%s.%i",cs.get_output_data_dir().c_str(),cs.get_refracted_file().c_str(),group);
  }
  else
  {
    //sprintf(file,"%s%s",cs.get_output_data_dir().c_str(),cs.get_refracted_file().c_str());
    file = format("%s%s",cs.get_output_data_dir().c_str(),cs.get_refracted_file().c_str());
  }
  // open output file
  std::ofstream newfile (file.c_str(),std::ios::out);
  if (newfile.is_open()==false)
  {
    cout << "\nLog::writeRefracted: "
          << "Error. Unable to open file = "
          << file << endl;
    assert(newfile.is_open()==true);
    exit(1);
  }
  Refractory & r(Refractory::instance());
  newfile.precision(12);
  if (r.getNumVertRefractedN() > 0)
  {
    // for each vertex moved too many time
    newfile << "# vertices refracted for moving MAX_TOUCHES="
          << cs.get_max_touches() << "times in group=" << group << endl;
    for (vp_cit i=r.beginN();i!=r.endN();++i)
    {
     // Vertex *v=(*i).first;
      (*i)->print(newfile);
    }
  }
  if (Refractory::instance().getNumVertRefractedInt() > 0)
  {
    // for each refracted vertex
    newfile << "# vertices refracted for causing intersecting faces"
          << " in group=" << group << endl;
    for (vp_cit i=r.beginInt();i!=r.endInt();++i)
    {
      //Vertex *v=(*i).first;
      (*i)->print(newfile);
    }
  }
  if (Refractory::instance().getNumVertRefractedAng() > 0)
  {
    // for each refracted vertex
    newfile << "# vertices refracted for causing "
          << "very small or very large edge angles"
          << " in group=" << group << endl;
    for (vp_cit i=r.beginAng();i!=r.endAng();++i)
    {
      //Vertex *v=(*i).first;
      (*i)->print(newfile);
    }
  }
  if (Refractory::instance().getNumVertRefractedSmallDisp() > 0)
  {
    // for each refracted vertex
    newfile << "# vertices refracted for attempting a small displacement"
          << " in group=" << group << endl;
    for (vp_cit i=r.beginSmallDisp();i!=r.endSmallDisp();++i)
    {
      //Vertex *v=(*i).first;
      (*i)->print(newfile);
    }
  }
  if (Refractory::instance().getNumVertRefractedOctreeVio() > 0)
  {
    // for each refracted vertex
    newfile << "# vertices refracted for breaching octree boundary"
          << " in group=" << group << endl;
    for (vp_cit i=r.beginOctreeVio();i!=r.endOctreeVio();++i)
    {
      //Vertex *v=(*i).first;
      (*i)->print(newfile);
    }
  }
  newfile.close();
}

void Log::writeRefractedNow (int const & group,int code)
{
  Controls & cs(Controls::instance());
  std::string file;
  // create output filename
  if (cs.get_append_group_number())
  {
    file = format("%s%s.%i.%i",cs.get_output_data_dir().c_str(),cs.get_refracted_file().c_str(),group,code);
  }
  else
  {
    file = format("%s%s.%i",cs.get_output_data_dir().c_str(),cs.get_refracted_file().c_str(),code);
  }
  // open output file
  std::ofstream newfile (file.c_str(),std::ios::app);
  if (newfile.is_open()==false)
  {
    cout << "\nLog::writeRefracted: "
          << "Error. Unable to open file = "
          << file << endl;
    assert(newfile.is_open()==true);
    exit(1);
  }
  newfile << N << " 1." << code << "\n";
  newfile.close();
}

/** Write to STDOUT detailed information about last vertex move.
*/

void Log::writeDetailedMoveInfo (void)
{
  cout.precision(12);
  cout <<  "    position   ["
        << p_orig.p[0] << " "
        << p_orig.p[1] << " "
        << p_orig.p[2] << "]->["
        << (*p_new).p[0] << " "
        << (*p_new).p[1] << " "
        << (*p_new).p[2] << "]\n"
        << "closest point  ["
        << cp_orig.p[0] << " "
        << cp_orig.p[1] << " "
        << cp_orig.p[2] << "]->["
        << cp_new.p[0] << " "
        << cp_new.p[1] << " "
        << cp_new.p[2] << "]\n"
        << "         ecw   ["
        << sepdis[0] << "]->["
        << sepdis[1] << "]\n"
        << "virtual disp   ["
        << vd_val[0] << "]->["
        << vd_val[1] << "]\n"
        << "   topN rank   ["
        << topN_val[0] << "]->["
        << topN_val[1] << "]\n\n";
  if (vd_val[0]!=vd_val[0])
  {
    cout << "\n\nLog::printDetailInfo: Error: "
          << "vd_val[0] = " << vd_val[0] << endl << endl;
    assert(vd_val[0]==vd_val[0]);
    exit(1);
  }
}

/** Write to STDOUT a brief summary information about the last vertex move.
 * \param[in] group Current group number (used as filename suffix).
 */

void Log::writeMoveSummary (int const & group, const int & group_size)
{
  Vertex_Schedule & vs(Vertex_Schedule::instance());
  cout.precision(12);
  cout << endl;
//  if (vs.getCurrentVertex()==vs.getSeedVertex())
//  {
//    cout << "SEEDSEEDSEEDSEEDSEEDSEED\n";
//  }
  std::string name = vs.getCurrentVertex()->getObject()->getName();
  cout << "       count   " << vs.getNumMovedVertsGroup()
        << " of "            << group_size << endl
        << "       group   " << group << endl
        << "      object   " << name << endl;
        if (vs.getCurrentVertex()==vs.getSeedVertex())
        {
          cout << " vertex SEED   ";
        }
        else
        {
          cout << "      vertex   ";
        }
        cout << vs.getCurrentVertex()->getIndex() << endl
//        << "vertex         " << vs.getCurrentVertex()->getIndex() << endl
        << "displacement   " << sqrt(State::instance().getVD2()) << endl
        //<< "world energy   " << Container::instance().getEnergy() << endl
        << "        gain   " << Gain_Schedule::instance().getGain() << endl;
}

/** Write to file elapsed time since last call of this function.
 * \param[in] message Message to write to file.
 */

void Log::recordTime (std::string const & message)
{
  time_t after,diff;
  after = time (NULL);
  diff = after-currtime;
  Mfile << message << " " << diff << " seconds\n";
  Mfile.flush();
  currtime = after;
}

/** Initialize time to input value.
 * \param[in] t New time value.
 */

void Log::setTime (time_t t)
{
  currtime = t;
}

/** Initialize output stream to main log file.
*/

void Log::openMainFile (void)
{
  Controls & cs(Controls::instance());
  //char file[cs.get_max_filename_size()];
  //sprintf(file,"%s%s",cs.get_output_data_dir().c_str(),cs.get_main_log_file().c_str());
  std::string file = format("%s%s",cs.get_output_data_dir().c_str(),cs.get_main_log_file().c_str());
  Mfile.open(file.c_str());
  if (Mfile.is_open()==false)
  {
    cout << "\nFailed to open file " << file << ".\n";
    assert(Mfile.is_open()==true);
    exit(1);
  }
}

/** Write to STDOUT debug data for vertex move detailed information.
*/

void Log::printBad (Vertex const * const current_vertex) const
{
  if (current_vertex==NULL){cout << "vertex* is NULL.\n";}
  if (p_orig.p[0]<0){cout << "orig pos not set.\n";}
  if ((*p_new).p[0]<0){cout << "new pos not set.\n";}
  if (cp_orig.p[0]<0){cout << "orig clos pos not set.\n";}
  if (cp_new.p[0]<0){cout << "new clos pos not set.\n";}
  if (topN_val[0]<0){cout << "orig topN_val not set.\n";}
  if (topN_val[1]<0){cout << "new topN_val not set.\n";}
  if (vd_val[0]<0){cout << "orig vd_val not set.\n";}
  if (vd_val[1]<0){cout << "new vd_val not set.\n";}
  if (sepdis[0]<0){cout << "orig sepdis not set.\n";}
  if (sepdis[1]<0){cout << "new sepdis not set.\n";}
  cout << "p_orig[0] = " << p_orig.p[0] << endl;
  cout << "p_new[0] = " << (*p_new).p[0] << endl;
  cout << "cp_orig[0] = " << cp_orig.p[0] << endl;
  cout << "cp_new[0] = " << cp_new.p[0] << endl;
  cout << "topN_val[0] = " << topN_val[0] << endl;
  cout << "topN_val[1] = " << topN_val[1] << endl;
  cout << "vd_val[0] = " << vd_val[0] << endl;
  cout << "vd_val[1] = " << vd_val[1] << endl;
  cout << "sepdis[0] = " << sepdis[0] << endl;
  cout << "sepdis[1] = " << sepdis[1] << endl;
}

/** Assert values are inside valid range
 * for class elements carring detailed vertex move information.
 * \param[in] current_vertex Vertex of interest,
 * likely, the last vertex moved.
 */

bool Log::isGood (Vertex const * const current_vertex) const
{
  return current_vertex!=NULL &&
        (p_orig.p[0]>0) &&
        ((*p_new).p[0]>0) &&
        (cp_orig.p[0]>0) &&
        (cp_new.p[0]>0) &&
        (topN_val[0]>0) &&
        (topN_val[1]>0) &&
        (vd_val[0]>=0) &&
        (vd_val[1]>=0) &&
        (sepdis[0]>=0) &&
        (sepdis[1]>0);
}

/** Record detailed move information of current vertex
 * before the move is attempted.
 * \param[in] v Vertex to be moved.
 */

void Log::setDetailedInfoPreMove (Vertex * const v) 
{
  vector3 const * ptr = v->getPos();
  // WHY MAKE COPY? Store v->getPos isntead.
  p_orig.p[0]=(*ptr).p[0];
  p_orig.p[1]=(*ptr).p[1];
  p_orig.p[2]=(*ptr).p[2];
  sepdis[0] = sqrt(v->getSqSepDist());

  int rank;
  tv_it t;
  if (Virtual_Disp::instance().getVertAndRank(v,t,rank)==true)
  {
    vd_val[0] = sqrt((*t).first);
    topN_val[0] = Virtual_Disp::instance().getNumVertsInVirtDispMap()-rank;
  }
  else
  {
    vd_val[0] = 0;
    topN_val[0] = 0;
  }
  if (v->getClosestFace()!=NULL)
  {
    double dd;
    Container::instance().findClosestPtInFaceToLocation(*v->getPos(),v->getClosestFace(),cp_orig,dd);
  }
}

/** Record detailed move information of current vertex
 * after a successfuly move is made.
 * \param[in] v Last vertex moved.
 */

void Log::setDetailedInfoPostMove (Vertex * v) 
{
  p_new = v->getPos();
  sepdis[1] = sqrt(v->getSqSepDist());
  int rank;
  tv_it t;
  if (Virtual_Disp::instance().getVertAndRank(v,t,rank)==true)
  {
    vd_val[1] = sqrt((*t).first);
    topN_val[1] = Virtual_Disp::instance().getNumVertsInVirtDispMap()-rank;
  }
  else
  {
    vd_val[1] = 0;
    topN_val[1] = 0;
  }
  if (v->getClosestFace()!=NULL)
  {
    double dd;
    Container::instance().findClosestPtInFaceToLocation(*v->getPos(),v->getClosestFace(),cp_new,dd);
  }
}

/** Gather statistics on number of faces and boxes
 * invloved in search for closest point to each vertex.
 * \param[in] face_count Number of faces returned from octree in last search
 * for closest point to current vertex.
 * \param[in] leaves_count Number of leaves checked in last search
 * for closest point to current vertex.
 * \param[in] f_check_count Number of faces checked in last search
 * for closest point to current vertex.
 */

void Log::updateClosestPtStats (int const & face_count,int const & leaves_count,int const & f_check_count)
{
  double previous_f_check_sum = num_verts*f_check_mean;
  double previous_face_sum = num_verts*face_mean;
  double previous_leaves_sum  = num_verts*leaves_mean;

  num_verts++;
  f_check_mean = (previous_f_check_sum+f_check_count)/num_verts; 
  if (f_check_count<f_check_min) f_check_min=f_check_count;
  if (f_check_count>f_check_max) f_check_max=f_check_count;
  face_mean = (previous_face_sum+face_count)/num_verts; 
  if (face_count<face_min) face_min=face_count;
  if (face_count>face_max) face_max=face_count;
  leaves_mean  = (previous_leaves_sum + leaves_count)/num_verts; 
  if (leaves_count<leaves_min) leaves_min=leaves_count;
  if (leaves_count>leaves_max) leaves_max=leaves_count;
}

/** Write to output stream summary statistics
 * for closest point to vertex searches.
 *
 * \param[in] target Pre-initialized output stream.
 */

void Log::printClosestPtStats (std::ostream & target) const
{
  //char str[1024];
  //sprintf(str,"\nLog: number of vertices in closest pt search stats  = ");
  //sprintf(str,"%s%d\n",str,num_verts);
  //target << str;
  //sprintf(str,"Log: number of octree leaves visited in closest pt search\n");
  //sprintf(str,"%s     (min,mean,max) = (%d,%g,%d)\n",str,leaves_min,leaves_mean,leaves_max);
  //target << str;
  //sprintf(str,"Log: number of faces returned from octree in closest pt search\n");
  //sprintf(str,"%s     (min,mean,max) = (%d,%g,%d)\n",str,face_min,face_mean,face_max);
  //target << str;
  //sprintf(str,"Log: number of faces checked in closest pt search\n");
  //sprintf(str,"%s     (min,mean,max) = (%d,%g,%d)\n\n",str,f_check_min,f_check_mean,f_check_max);
  //target << str;
  target << "\nLog: number of vertices in closest pt search stats  = "
         << num_verts << "\n"
         << "Log: number of octree leaves visited in closest pt search\n"
         << "     (min,mean,max) = ("
         << leaves_min  << ","
         << leaves_mean << ","
         << leaves_max  << ")\n"
         << "Log: number of faces returned from octree in closest pt search\n"
         << "     (min,mean,max) = ("
         << face_min  << "," 
         << face_mean << ","
         << face_max  << ")\n"
         << "Log: number of faces checked in closest pt search\n"
         << "     (min,mean,max) = ("
         << f_check_min  << ","
         << f_check_mean << ","
         << f_check_max  << ")\n\n";
}

/** Write to output stream the number of vertices recorded nonnice.
 *
 * \param[in] target Pre-initialized output stream.
 */

void Log::printNumNonnice (std::ostream & target) const
{
  //char str[1024];
  //sprintf(str,"\nLog: number of nonnice vertices = ");
  //sprintf(str,"%s%d\n\n",str,Nice::instance().getNonniceCount(false));
  //target << str;
  //char str[1024];
  target << "\nLog: number of nonnice vertices = "
         << Nice::instance().getNonniceCount(false) << "\n\n";
}

/** Write to output stream the number of faces recorded as intersected.
 *
 * \param[in] target Pre-initialized output stream.
 */

void Log::printNumInt (std::ostream & target) const
{
  //char str[1024];
  //int i = Intersecting_Faces::instance().getCountOfIntFaces(false);
  //sprintf(str,"\nLog: number of face intersections = %d\n\n",i);
  //target << str;
  //char str[1024];
  int i = Intersecting_Faces::instance().getCountOfIntFaces(false);
  target << "\nLog: number of face intersections = " << i << "\n\n";
}

/** Update statistics on the number of boxes assigned to each face.
 * \param[in] box_count Number of space partitions assigned to this face.
 */

void Log::updateBoxesPerFaceStats (int const & box_count)
{
  double previous_box_sum = num_faces*bpf_mean;

  num_faces++;
  bpf_mean = (previous_box_sum+box_count)/num_faces; 
  if (box_count<bpf_min) bpf_min=box_count;
  if (box_count>bpf_max) bpf_max=box_count;
}

/** Write to output stream recorded statistics
 * of the number of faces assigned to each box
 * and the number of boxes assigned to each face.
 * \param[in] target Pre-initialized output stream.
 */

void Log::printPartitioningStats (std::ostream & target) const
{
  //char str[1024];
  //sprintf(str,"\nLog: number of boxes per face\n");
  //sprintf(str,"%s(min,mean,max) = (%d,%g,%d)\n",str,bpf_min,bpf_mean,bpf_max);
  //target << str;
  //sprintf(str,"Log: number of faces per box\n");
  //sprintf(str,"%s(min,mean,max) = (%d,%g,%d)\n",str,fpb_min,fpb_mean,fpb_max);
  //target << str;
  //sprintf(str,"Log: number of boxes with no faces = %d\n\n",num_empty_boxes);
  //target << str;
  target << "\nLog: number of boxes per face\n"
         << "(min,mean,max) = ("
         << bpf_min  << ","
         << bpf_mean << ","
         << bpf_max  << ")\n"
         << "Log: number of faces per box\n"
         << "(min,mean,max) = ("
         << fpb_min  << ","
         << fpb_mean << ","
         << fpb_max  << ")\n"
         << "Log: number of boxes with no faces = " << num_empty_boxes << "\n\n";
}

/** Update statistics on number of vertices per moved set.
 * (noting that not all vertices in a moved set are allowed to move).
 */

void Log::updateMovedVertsFromSet(void)
{
  // decrement num set since Log::updateVertexSchedulingStats
  // already incremented counter
  double previous_moved_sum = (num_sets-1)*moved_mean;

  Vertex_Schedule & vs(Vertex_Schedule::instance());
  int num_moved_verts_last_set = vs.getNumMovedVertsGroup()
        -vs.getReferenceCountMovedVerts(); 
  moved_mean = (previous_moved_sum+num_moved_verts_last_set)/num_sets; 
  if (num_moved_verts_last_set<moved_min) moved_min=num_moved_verts_last_set;
  if (num_moved_verts_last_set>moved_max) moved_max=num_moved_verts_last_set;
}

/** Update statistics on number of vertices per moved set.
 * \param[in] vertex_count Number of vertices in last moved set
 * (noting that not all vertices in a moved set are allowed to move).
 */

void Log::updateVertexSchedulingStats (int const & vertex_count)
{
  double previous_vertex_sum = num_sets*vps_mean;

  num_sets++;
  vps_mean = (previous_vertex_sum+vertex_count)/num_sets; 
  if (vertex_count<vps_min) vps_min=vertex_count;
  if (vertex_count>vps_max) vps_max=vertex_count;
}

/** Write to output streqam statistics on number of vertices per moved set.
 * \param[in] target Pre-initialized output stream.
 */

void Log::printVertexSchedulingStats (std::ostream & target) const
{
  //char str[1024];
  //sprintf(str,"\nLog: ************************************************************\n");
  //sprintf(str,"%sLog: number of vertices per moved set\n",str);
  //sprintf(str,"%s(min,mean,max) = (%d,%g,%d)\n",str,vps_min,vps_mean,vps_max);
  //target << str;
  //sprintf(str,  "Log: (note not all vertices in a \"moved set\" ");
  //sprintf(str,"%sare allowed to move)\n",str);
  //target << str;
  //sprintf(str,"\nLog: number of vertices actually moved per set\n");
  //sprintf(str,"%s(min,mean,max) = (%d,%g,%d)\n\n",str,moved_min,moved_mean,moved_max);
  //target << str;
  target << "\nLog: ************************************************************\n"
         << "Log: number of vertices per moved set\n"
         << "(min,mean,max) = ("
         << vps_min  << ","
         << vps_mean << ","
         << vps_max  << ")\n"
         << "Log: (note not all vertices in a \"moved set\" "
         << "are allowed to move)\n"
         << "\nLog: number of vertices actually moved per set\n"
         << "(min,mean,max) = ("
         << moved_min  << ","
         << moved_mean << ","
         << moved_max  << ")\n\n";
}

/** Gather statistics on number of vertices
 * undergoing full and partial closest point searches.
 * \param[in] full_search Number of vertices given 
 * a full closest point seach.
 * \param[in] partial_search Number of vertices given 
 * a partial closest point seach.
 * \param[in] ss Collection of vertex counts quantifying how often
 * the closest face to vertex changed.
 */

void Log::updateClosestPtSearchStats (int const & full_search,
                                      int const & partial_search,
                                      Search_Stats const & ss)
{
  // Note # moved vertices not updated at this point
  // so N+1 == # vertices moved so far
  double previous_fs_sum = N*fs_mean;
  double previous_ps_sum = N*ps_mean;

  fs_mean = (previous_fs_sum +    full_search)/(N+1.0); 
  ps_mean = (previous_ps_sum + partial_search)/(N+1.0); 
  if    (full_search<fs_min) fs_min=full_search;
  if    (full_search>fs_max) fs_max=full_search;
  if (partial_search<ps_min) ps_min=partial_search;
  if (partial_search>ps_max) ps_max=partial_search;

  double previous_fs_face_change_sum = N*fs_face_change_mean;
  double previous_ps_face_change_sum = N*ps_face_change_mean;
//  double previous_fs_face_adj_sum    = N*fs_mean_face_change_is_adj;
//  double previous_ps_face_adj_sum    = N*ps_mean_face_change_is_adj;

  fs_face_change_mean = (previous_fs_face_change_sum+ss.fs_changed)/(N+1.0); 
  ps_face_change_mean = (previous_ps_face_change_sum+ss.ps_changed)/(N+1.0); 

  if (ss.fs_changed<fs_face_change_min) fs_face_change_min=ss.fs_changed;
  if (ss.fs_changed>fs_face_change_max) fs_face_change_max=ss.fs_changed;
  if (ss.ps_changed<ps_face_change_min) ps_face_change_min=ss.ps_changed;
  if (ss.ps_changed>ps_face_change_max) ps_face_change_max=ss.ps_changed;

//  fs_mean_face_change_is_adj = 0.0;
//  ps_mean_face_change_is_adj = 0.0;
//  if (ss.fs_changed>0) fs_mean_face_change_is_adj = (previous_fs_face_adj_sum+ss.fs_changed_adj)/(N+1); 
//  if (ss.ps_changed>0) ps_mean_face_change_is_adj = (previous_ps_face_adj_sum+ss.ps_changed_adj)/(N+1); 
}

/** Write to output stream summary statistics
 * for closest point to vertex searches.
 *
 * \param[in] target Pre-initialized output stream.
 */

void Log::printClosestPtSearchStats (std::ostream & target) const
{
//  char str[1024];
//  sprintf(str,"Log: ************************************************************\n");
//  sprintf(str,"%sLog: number of vertices undergoing full closest pt search\n",str);
//  sprintf(str,"%s     (min,mean,max) = (%d,%g,%d)\n",str,fs_min,fs_mean,fs_max);
//  target << str;
//  sprintf(str,"\n  Log: ... And recorded a change in closest face to vertex\n");
//  sprintf(str,"%s       (min,mean,max) = (%.15g,%.15g,%.15g)\n\n",str,fs_face_change_min,fs_face_change_mean,fs_face_change_max);
//  target << str;
//  sprintf(str,"Log: number of vertices undergoing partial closest pt search\n");
//  sprintf(str,"%s(min,mean,max) = (%d,%g,%d)\n",str,face_min,face_mean,face_max);
//  target << str;
//  sprintf(str,"\n  Log: ... And recorded a change in closest face to vertex\n");
//  sprintf(str,"%s       (min,mean,max) = (%.15g,%.15g,%.15g)\n\n",str,ps_face_change_min,ps_face_change_mean,ps_face_change_max);
//  target << str;
//  sprintf(str,"Log: number of moved vertices in closest pt search stats  = ");
//  sprintf(str,"%s%d\n\n",str,N);
//  target << str;
  target << "Log: number of vertices undergoing full closest pt search\n"
         << "     (min,mean,max) = ("
         << fs_min  << ","
         << fs_mean << " ,"
         << fs_max  << ")\n"
         << "\n  Log: ... And recorded a change in closest face to vertex\n"
         << "       (min,mean,max) = ("
         << fs_face_change_min  << ","
         << fs_face_change_mean << ","
         << fs_face_change_max  << ")\n\n"
         << "Log: number of vertices undergoing partial closest pt search\n"
         << "(min,mean,max) = ("
         << face_min  << ","
         << face_mean << ","
         << face_max  << ")\n"
         << "\n  Log: ... And recorded a change in closest face to vertex\n"
         << "       (min,mean,max) = ("
         << ps_face_change_min  << ","
         << ps_face_change_mean << ","
         << ps_face_change_max  << ")\n\n"
         << "Log: number of moved vertices in closest pt search stats  = "
         << N << "\n"
         << "Log: ************************************************************\n\n";
}

/** Write command parameter settings to file.
 */

void Log::writeCommandSettings (void)
{
  Controls & cs(Controls::instance());
  // create output filename
  //char file[cs.get_max_filename_size()];
  //sprintf(file,"%s%s",cs.get_output_data_dir().c_str(),cs.get_control_file().c_str());
  std::string file = format("%s%s",cs.get_output_data_dir().c_str(),cs.get_control_file().c_str());
  // open output file
  std::ofstream newfile (file.c_str(),std::ios::out);
  if (newfile.is_open()==false)
  {
    cout << "\nLog::writeCommandSettings: "
          << "Error. Unable to open file = "
          << file << endl;
    assert(newfile.is_open()==true);
    exit(1);
  }
  newfile.precision(12);
  newfile << cs.getCommandSettings();
  newfile.close();
}
