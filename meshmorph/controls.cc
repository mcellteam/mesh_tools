// Author: Justin Kinney
// Date: Sep 2008

#include "controls.h"

#include <iostream>
#include <getopt.h>
#include <cmath>

#include "container.h"

using std::cout;
using std::endl;

Controls * Controls::only_one = NULL;

Controls & Controls::instance(void)
{
  // Not thread-safe.
  // -- lock mutex
  if (only_one == NULL)
    only_one = new Controls();
  // -- unlock mutex
  return *only_one;
}

Controls::Controls (void)
  :MAX_ITEMS_PER_LEAF                  (12),
  MAX_OCTREE_DEPTH                     (12),
  MIN_CELL_SIZE                        (32.0),
  OCTREE_MIN_X                         (700),
  OCTREE_MIN_Y                         (700),
  OCTREE_MIN_Z                         (2900),
  OCTREE_WIDTH                         (9000),
  PI                                   (3.14159265358979),
  MAX_FILENAME_SIZE                    (1024),
  WRITE_MESH_NOW                       (0),
  WRITE_VERBOSE_INIT                   (false),
  WRITE_REFRACTED_VERTICES_TO_FILE     (true),
  WRITE_INTERSECTED_FACES_TO_FILE      (true),
  FORMAT_INTERSECTED_FACES             ("detail"),
  WRITE_NONNICE_VERTICES_TO_FILE       (true),
  FORMAT_NONNICE_VERTICES              ("detail"),
  WRITE_EVERY_GROUP                    (true),
  WRITE_ECW_TO_FILE                    (true),
  ECW_SAMPLING_LENGTH                  (40),
  WRITE_VERTEX_MOVE_HISTOGRAM          (false),
  APPEND_GROUP_NUMBER                  (false),
  STRICT_FACE_INTERSECTION_PREVENTION  (true),
  EDGE_ANGLE_THRESHOLD                 (10*PI/180.0),
  MY_DOUBLE_EPSILON                    (1E-10),
  EPSILON                              (1E-7),
  UPDATE_REGION_SIZE                   (160),
  SEED_REGION_SIZE                     (20),
  SEARCH_RADIUS_SQ                     (16900.0),
  NUMBER_RADIUS_STEPS                  (1),
  MIN_SEARCH_CONE_RADIUS               (60.0),
  CLOSEST_POINT_ANGLE                  (30*PI/180.0),
  CLOSEST_POINT_COSINE                 (0.866025403784439),
  CLOSEST_POINT_SINE                   (0.50),
  USE_EDGE_REFERENCE_LENGTH            (false),
  INTERSECTION_WEIGHT                  (0.001),
  EDGE_LENGTH_WEIGHT                   (0.099),
  ECW_WEIGHT                           (0.45),
  EDGE_ANGLE_WEIGHT                    (0.45),
  ECW_GAIN                             (0.1),
  EDGE_ANGLE_GAIN                      (20.0),
  EDGE_LENGTH_GAIN                     (0.2),
  ASPECT_RATIO_GAIN                    (0.01),
  ASPECT_RATIO_THRESHOLD               (3),
  OVERALL_GAIN                         (1.0),
  GROUP_SIZE                           (10000),
  REDEFINE_GROUP_SIZE                  (1),
  NUM_GROUPS                           (1),
  REDEFINE_NUM_GROUPS                  (1),
  REFRACTORY_PERIOD                    (1000),
  MAX_TOUCHES                          (10),
  ENERGY_WINDOW                        (10000),
  ENERGY_SAMPLE_PERIOD                 (100000),
  MIN_DISPLACEMENT_SQ                  (1E-2),
  TARGET_ECW                           (20.0),
  ECW_THRESHOLD                        (-1.0),
  TARGET_ECW_HIGH                      (300.0),
  TARGET_ECW_LOW                       (15.0),
  GAIN_STEP                            (-1.0),
  MAX_ACTUAL_DISPL_FRACTION            (1.0),
  MAX_RUNTIME                          (0),
  PRINT_PERIOD                         (1000),
  BEGIN_SHORT_PRINT_PERIOD             (0),
  ENABLE_VTRACK                        (false),
  DISABLE_GAIN_SCHEDULING              (false),
  VECTOR_RESERVE                       (1000),
  INPUT_DATA_DIR                       ("./"),
  OUTPUT_DATA_DIR                      ("./"),
  FROZEN_VERTICES_FILE                 (""),
  VERTEX_SEQUENCE_FILE                 (""),
  MESH_OUTPUT_SUFFIX                   (""),
  MAIN_LOG_FILE                        ("main.log"),
  CONT_LOG_FILE                        ("timeline.log"),
  SEP_LOG_FILE                         ("extracellular_widths.dat"),
  OBJECT_LIST_FILE                     ("object_list.log"),
  CONTROL_FILE                         ("command_options.dat"),
  VERTEX_SELECTION_FILE                ("vertex_select_histo.dat"),
  REFRACTED_FILE                       ("refracted_vertices.dat"),
  INTERSECTED_FILE                     ("intersected_faces.dat"),
  NONNICE_FILE                         ("nonnice_vertices.dat"),
  SMALL_ECW_THRESHOLD                  (0.2),
  DISABLE_MESSAGES                     (0),
  MEASURE_ECW_AND_EXIT                 (0)
{
}

void Controls::updatePrintPeriod (int count)
{
  if (BEGIN_SHORT_PRINT_PERIOD > 0)
  {
    if (count>=BEGIN_SHORT_PRINT_PERIOD)
    {
      PRINT_PERIOD = 1;
      BEGIN_SHORT_PRINT_PERIOD=0;
    }
  }
}

void Controls::redefineGroup (void)
{
  // The following settings were experimentally determined to be adequate.
  if (REDEFINE_GROUP_SIZE==true)
  {
    // Let the group size equal the total number of vertices in the model.
    GROUP_SIZE = Container::instance().getVertexCount();
    int p = static_cast<int>(floor(static_cast<float>(GROUP_SIZE)/static_cast<float>(MAX_TOUCHES)));
    if (p<1)
      REFRACTORY_PERIOD = 1;
    else
      REFRACTORY_PERIOD = p;
  }
  if (REDEFINE_NUM_GROUPS==true)
  {
    // Let the number of groups vary linearly with target ecw.
    NUM_GROUPS = 50 + static_cast<int>(fabs(TARGET_ECW-20));
  }
}

void Controls::recordOctreeSize (double size)
{
  Container & c(Container::instance());
  OCTREE_MIN_X = c.getWorld(0);
  OCTREE_MIN_Y = c.getWorld(2);
  OCTREE_MIN_Z = c.getWorld(4);
  OCTREE_WIDTH = size;
}

/** Ensure pathname ends with exactly one '/'.
 * \param[in] ptr Arbitrary pathname.
 * \return Input pathname with '/' added if absent.
 */

std::string Controls::processDir (char * ptr)
{
  std::string pathname = ptr;
  if (pathname.find_last_of('/')==std::string::npos)
  {
    pathname.push_back('/');
    return pathname;
  }
  size_t pos = pathname.find_last_of('/');
  if ((pos+1)==pathname.size())
  {
    return pathname;
  }
  else
  {
    pathname.push_back('/');
    return pathname;
  }
}

/** Create string from floating point number.
 * \param[in] i Number of interest.
 * \return String containing input number.
 */

std::string Controls::d2str (double const & i)
{
  char file[get_max_filename_size()];
  sprintf(file,"%.15g",i);
  return std::string(file);
}

/** Create string from integer number.
 * \param[in] i Number of interest.
 * \return String containing input number.
 */

std::string Controls::i2str (int const & i)
{
  char file[get_max_filename_size()];
  sprintf(file,"%i",i);
  return std::string(file);
}

/** Create meshmorph usage message.
 * \return Meshmorph usage message.
 */

std::string Controls::getUsageMessage (void)
{
 std::string message = "\n";
  message=message+
  "NAME\n"+
  "       meshmorph - move meshes to control extracellular space\n"+
  "\nSYNOPSIS\n"+
  "       meshmorph [options]\n"+
  "\nDESCRIPTION\n"+
  "       Meshmorph moves vertices one at a time to relax\n"+
  "       a spring model of the cell membranes.\n"+
  "\nEXAMPLES\n"+
  "       meshmorph -i input -o output -t 20 -s 0.1 -a 0.8 -v frozen_vertices.dat\n"+
  "              Read meshes from directory 'input' and write new morphed\n"+
  "              meshes to directory 'output'. The target extracellular width\n"+
  "              is 20 data units and the relative weights for extracellular\n"+
  "              width and edge angle are 10 and 80, respectively, out of 100.\n"+
  "              The vertices as specified in the file 'frozen_vertices.dat'\n"+
  "              in the form 'object_name vertex_index' (one vertex per file)\n"+
  "              are not to be moved.\n"+
  "\nOPTIONS\n"+
  "       --verbose_init\n"+
  "              Write initialization information for each file to stdout.\n"+
  "              Default is concise initialization with minimal info writte to stdout.\n\n"+
  "       --no_refracted\n"+
  "              Do not write to file a summary of vertices put into refractory state during morphing.\n"+
  "              Default is to write summary of refacted vertices to file\n"+
  "              specified by --refracted_file in output directory specified by -o\n"+
  "              after completion of every group of moved vertices.\n\n"+
  "       --no_intersected\n"+
  "              Do not write to file a summary of faces found to be intersecting.\n"+
  "              Default is to write summary of intersected faces to file\n"+
  "              specified by --intersected_file in output directory specified by -o\n"+
  "              after completion of every group of moved vertices.\n\n"+
  "       --no_nonnice\n"+
  "              Do not write to file a summary of nonnice vertices, i.e. vertices located inside an object.\n"+
  "              Default is to write summary of nonnice vertices to file\n"+
  "              specified by --nonnice_file in output directory specified by -o\n"+
  "              after completion of every group of moved vertices.\n\n"+
  "       --no_ecw\n"+
  "              Do not measure and do not write to file the extracellular width\n"+
  "              of the model.\n"+
  "              Default is to measure ecw of model with surface sampling density\n"+
  "              defined by --ecw_sampling_length and to write data to file\n"+
  "              specified by --sep_log_file in output directory specified by -o.\n\n"+
  "       --no_every_group\n"+
  "              Do not measure and do not write to file the various diagnostic\n"+
  "              and descriptive output data from the model.\n"+
  "              Default is to measure and write data to file\n"+
  "              in output directory specified by -o\n"+
  "              after completion of every group of moved vertices.\n\n"+
  "       --move_histo\n"+
  "              Write the number of times each vertex in model was moved\n"+
  "              to the file specified by --vertex_selection_file in output\n"+
  "              directory specified by -o after completion of every group of moved vertices.\n"+
  "              Default is to do nothing.\n\n"+
  "       --append_group\n"+
  "              Incorporate group number into output data filenames.\n"+
  "              Default is to not use group number in filenames.\n\n"+
  "       --no_int_prevention\n"+
  "              Allow vertex moves that create intersections between pairs of faces\n"+
  "              which did not intersect each other before the move. This may be necessary\n"+
  "              during the process of separating intersecting objects.\n"+
  "              Default is to prevent vertex moves that create de novo face intersections.\n\n"+
  "       --ref_orig_edge_length\n"+
  "              Use original edge length to measure edge stretch and determine\n"+
  "              force of edge compression or tension. WARNING: this option\n"+
  "              may create faces with high aspect ratios unless a compensatory\n"+
  "              force is added, i.e. force proportional to adjacent face aspect ratio.\n"+
  "              Default is to compare edge length to the instantaneous mean edge\n"+
  "              length of first adjacent face of edge.\n\n"+
  "       --verbose_vertex_move\n"+
  "              Write detailed information about each moved vertex to stdout,\n"+
  "              including position, extracellular width, virtual displacement,\n"+
  "              virtual displacement rank, and closest point, both before and\n"+
  "              after move. Frequency of detailed information output is specified\n"+
  "              by -9. Default is to do nothing.\n\n"+
  "       --disable_gain_scheduling\n"+
  "              Disable gain scheduling so that the proportionality factor between\n"+
  "              force on a vertex and vertex displacement is constant throughout\n"+
  "              execution of the program. The proportionality constant is overall\n"+
  "              gain and is specified by -O. Default is to follow a decaying gain schedule.\n\n"+
  "       --disable_messages\n"+
  "              Disable Refractory and Virtual_Displacement messages for faster program execution.\n"+
  "              Default is to write all messages to STDOUT.\n\n"+
  "       --measure_ecw_and_exit\n"+
  "              Measure the extracellular width in the model at sampling density\n"+
  "              defined by --ecw_sampling_length and exit.\n"+
  "              Default is to do nothing.\n\n"+
  "       -1 NUM, --max_items_per_leaf=NUM\n"+
  "              First argument to construction of octree.\n"+
  "              Probably specifies a criterion for deciding if a cell in tree.\n"+
  "              should be divided, if possible.\n"+
  "              Default is '" + i2str(MAX_ITEMS_PER_LEAF) + "'.\n\n"+
  "       -2 NUM, --max_octree_depth=NUM\n"+
  "              Second argument to construction of octree.\n"+
  "              Probably specifies a maximum number of divisions\n"+
  "              of a cell in tree.\n"+
  "              Default is '" + i2str(MAX_OCTREE_DEPTH) + "'.\n\n"+
  "       -C NUM, --min_cell_size=NUM\n"+
  "              Third and final argument to construction of octree.\n"+
  "              Probably specifies the minimum dimensions of octreen\n"+
  "              cell below which cells will not be further divided.\n"+
  "              Units are same as meshes in input directory.\n"+
  "              Default is '" + d2str(MIN_CELL_SIZE) + "'.\n\n"+
  "       -3 NUM, --max_filename_size=NUM\n"+
  "              Maximum number of characters allowed in input and\n"+
  "              output filenames.\n"+
  "              Default is '" + i2str(MAX_FILENAME_SIZE) + "'.\n\n"+
  "       -4 NUM, --update_region_size=NUM\n"+
  "              The region of space in which the model will be\n"+
  "              updated after a vertex move.\n"+
  "              All vertices of all faces in octree cells that overlap this\n"+
  "              region of space will have their closest point, extracellular\n"+
  "              width, total force, virtual displacement, and virtual\n"+
  "              displacement rank updated.\n"+
  "              Default is '" + d2str(UPDATE_REGION_SIZE) + "'.\n\n"+
  "       -5 NUM, --seed_region_size=NUM\n"+
  "              Vertices are moved in sets consisting of a seed vertex\n"+
  "              (typically the vertex with the largest virtual displacement)\n"+
  "              and vertices of all faces in octree cells intersected by\n"+
  "              cube centered on seed vertex with side length\n"+ 
  "              equal to 2*NUM.\n"+
  "              Default is '" + d2str(SEED_REGION_SIZE) + "'.\n\n"+
  "       -6 NUM, --num_radius_steps=NUM\n"+
  "              The closest point to a vertex is found by first searching\n"+
  "              a region specified by -F. While no closest point is found\n"+
  "              the search region is expanded up to size specified by -G in NUM steps.\n"+
  "              Default is '" + i2str(NUMBER_RADIUS_STEPS) + "'.\n\n"+
  "       -g NUM, --group_size=NUM\n"+
  "              The morphing process is a sequence of vertex moves\n"+
  "              divided into groups of size NUM.\n"+
  "              Default is '" + i2str(GROUP_SIZE) + "'.\n\n"+
  "       -n NUM, --num_groups=NUM\n"+               
  "              The morphing process consists of NUM groups of vertex moves.\n"+
  "              Default is '" + i2str(NUM_GROUPS) + "'.\n\n"+
  "       -7 NUM, --refractory_period=NUM\n"+        
  "              Size of vertex move refractory period, i.e. the minimum number\n"+
  "              of vertex moves between moves of the same vertex.\n"+
  "              Default is '" + i2str(REFRACTORY_PERIOD) + "'.\n\n"+
  "       -8 NUM, --max_touches=NUM\n"+              
  "              Maximum number of moves each vertex may make per group.\n"+
  "              Default is '" + i2str(MAX_TOUCHES) + "'.\n\n"+
  "       -9 NUM, --vertex_move_print_period=NUM\n"+ 
  "              Number of vertex moves in cycle with disabled printing\n"+
  "              alternating with a single vertex move logged to stdout.\n"+
  "              Minimum value == 1 which will print each vertex move.\n"+
  "              Default is '" + i2str(PRINT_PERIOD) + "'.\n\n"+
  "       -Q NUM, --begin_short_print_period=NUM\n"+ 
  "              Number of vertex moves in cycle past which\n"+
  "              every vertex move will be logged to stdout.\n"+
  "              Minimum value == 1 which will print each vertex move.\n"+
  "              Default is '" + i2str(BEGIN_SHORT_PRINT_PERIOD) + "' which does nothing.\n\n"+
  "       -R NUM, --write_mesh_now=NUM\n"+ 
  "              Number of vertex moves in group after which\n"+
  "              the entire mesh model will be written to file once.\n"+
  "              Output file will have the suffix '_iteration_WRITE_MESH_NOW.mesh'.\n"+
  "              Minimum value == 1 which will write mesh after one vertex move.\n"+
  "              Default is '" + i2str(WRITE_MESH_NOW) + "' which does nothing.\n\n"+
  "       -B NUM, --vector_reserve=NUM\n"+           
  "              Initial size of vectors to avoid vector resizing\n"+
  "              and consequent data copying.\n"+
  "              Default is '" + i2str(VECTOR_RESERVE) + "'.\n\n"+
  "       -e NUM, --edge_angle_threshold=NUM\n"+ 
  "              If moving a vertex creates an edge angle less than NUM,\n"+
  "              then vertex is not moved. Units radians.\n"+
  "              Default is '" + d2str(EDGE_ANGLE_THRESHOLD) + "'.\n\n"+
  "       -F NUM, --min_search_cone_radius=NUM\n"+ 
  "              The closest point to a vertex is found by first searching\n"+
  "              a region specified by NUM. While no closest point is found\n"+
  "              the search region is expanded up to size specified by -G\n"+
  "              in a number of steps specified by -6.\n"+
  "              Default is '" + d2str(MIN_SEARCH_CONE_RADIUS) + "'.\n\n"+
  "       -G NUM, --search_radius=NUM\n"+ 
  "              The closest point to a vertex is found by first searching\n"+
  "              a region specified by -F. While no closest point is found\n"+
  "              the search region is expanded up to size NUM\n"+
  "              in a number of steps specified by -6.\n"+
  "              Default is '" + d2str(sqrt(SEARCH_RADIUS_SQ)) + "'.\n\n"+
  "       -H NUM, --search_cone_angle=NUM\n"+ 
  "              Candidate closest points to a vertex must lie\n"+
  "              within NUM degrees of vertex normal vector. Input unit is degrees.\n"+
  "              Default is '" + d2str(CLOSEST_POINT_ANGLE*180.0/PI) + "'.\n\n"+
  "       -f NUM, --face_intersection_weight=NUM\n"+ 
  "              Force on vertices due to intersecting faces\n"+
  "              is scaled by NUM. Range is 0.0 to 1.0.\n"+
  "              Default is '" + d2str(INTERSECTION_WEIGHT) + "'.\n\n"+
  "       -l NUM, --edge_length_weight=NUM\n"+ 
  "              The force generated on each vertex by deviations\n"+
  "              in the measured edge length from the reference\n"+
  "              edge length is scaled by NUM. Range 0.0 to 1.0.\n"+
  "              Default is '" + d2str(EDGE_LENGTH_WEIGHT) + "'.\n\n"+
  "       -s NUM, --ecw_weight=NUM\n"+ 
  "              The force generated on each vertex by deviations\n"+
  "              in the measured extracellular width from the desired\n"+
  "              extracellular width is scaled by NUM. Range is 0.0 to 1.0.\n"+
  "              Default is '" + d2str(ECW_WEIGHT) + "'.\n\n"+
  "       -a NUM, --edge_angle_weight=NUM\n"+        
  "              The force generated on each vertex by deviations\n"+
  "              in the measured angle of nearby edges from 180 degrees\n"+
  "              is scaled by NUM. Range is 0.0 to 1.0.\n"+
  "              Default is '" + d2str(EDGE_ANGLE_WEIGHT) + "'.\n\n"+
  "       -S NUM, --ecw_gain=NUM\n"+ 
  "              The force generated on each vertex by deviations\n"+
  "              in the measured extracellular width from the desired\n"+
  "              extracellular width is scaled by NUM.\n"+
  "              Default is '" + d2str(ECW_GAIN) + "'.\n\n"+
  "       -A NUM, --edge_angle_gain=NUM\n"+ 
  "              The force generated on each vertex by deviations\n"+
  "              in the measured angle of nearby edges from 180 degrees\n"+
  "              is scaled by NUM.\n"+
  "              Default is '" + d2str(EDGE_ANGLE_GAIN) + "'.\n\n"+
  "       -L NUM, --edge_stretch_gain=NUM\n"+ 
  "              The force generated on each vertex by deviations\n"+
  "              in the measured edge length from the reference\n"+
  "              edge length is scaled by NUM.\n"+
  "              Default is '" + d2str(EDGE_LENGTH_GAIN) + "'.\n\n"+
  "       -E NUM, --aspect_ratio_gain=NUM\n"+ 
  "              The force generated on a vertex due to the aspect ratio\n"+
  "              of an adjacent face is scaled by NUM if aspect ratio\n"+
  "              is larger than threshold defined by -c.\n"+
  "              Default is '" + d2str(ASPECT_RATIO_GAIN) + "'.\n\n"+
  "       -c NUM, --aspect_ratio_threshold=NUM\n"+ 
  "              Force is generated on a vertex if adjacent face aspect ratio\n"+
  "              is larger than NUM.\n"+
  "              Default is '" + d2str(ASPECT_RATIO_THRESHOLD) + "'.\n\n"+
  "       -O NUM, --overall_gain=NUM\n"+ 
  "              The cumulative force generated on each vertex\n"+
  "              due to contributions from all adjacent edges\n"+
  "              and the extracellular width is scaled by NUM.\n"+
  "              Default is '" + d2str(OVERALL_GAIN) + "'.\n\n"+
  "       -I NUM, --min_vertex_disp=NUM\n"+ 
  "              The minimum acceptable virtual displacment of any\n"+
  "              vertex is NUM such that if the virtual_displacement\n"+
  "              is less than NUM then the vertex is refracted.\n"+
  "              Default is '" + d2str(sqrt(MIN_DISPLACEMENT_SQ)) + "'.\n\n"+
  "       -t NUM, --target_ecw=NUM\n"+ 
  "              The meshes will try to be morphed so as to have\n"+
  "              an extracellular width of size NUM.\n"+
  "              Units are same as meshes in input directory.\n"+
  "              Default is '" + d2str(TARGET_ECW) + "'.\n\n"+
  "       -d NUM, --ecw_threshold=NUM\n"+ 
  "              If specified, then vertices with an extracellular width\n"+
  "              greater than or equal to NUM will be morphed so as to have\n"+
  "              an extracellular width of size TARGET_ECW_HIGH specified\n"+
  "              by -j option. Vertices with an extracellular width less\n"+
  "              than NUM will be morphed so as to have\n"+
  "              an extracellular width of size TARGET_ECW_LOW specified by -k option.\n"+
  "              Units are same as meshes in input directory.\n"+
  "              Default is to morph all vertices to the same extracellular\n"+
  "              specified by -t option.\n\n"+
  "       -j NUM, --target_ecw_high=NUM\n"+ 
  "              If ECW_THRESHOLD is specified with -d option,\n"+
  "              then vertices with an extracellular width\n"+
  "              greater than or equal to ECW_THRESHOLD will be morphed\n"+
  "              so as to have an extracellular width of size NUM.\n"+
  "              Units are same as meshes in input directory.\n"+
  "              Default is '" + d2str(TARGET_ECW_HIGH) + "'.\n\n"+
  "       -k NUM, --target_ecw_low=NUM\n"+ 
  "              If ECW_THRESHOLD is specified with -d option,\n"+
  "              then vertices with an extracellular width\n"+
  "              less than ECW_THRESHOLD will be morphed\n"+
  "              so as to have an extracellular width of size NUM.\n"+
  "              Units are same as meshes in input directory.\n"+
  "              Default is '" + d2str(TARGET_ECW_LOW) + "'.\n\n"+
  "       -m NUM, --gain_step=NUM\n"+ 
  "              If gain scheduling is not disabled, then after\n"+
  "              each group the maximum allwed gain is decremented\n"+
  "              by NUM*overall_gain/num_groups. Therefore, if NUM is positive\n"+
  "              then the maximum gain decreases. If NUM is negative\n"+
  "              then the maximum gain increases.\n"+
  "              Default is '" + d2str(GAIN_STEP) + "'.\n\n"+
  "       -J NUM, --max_displ_fraction=NUM\n"+       
  "              The maximum ratio of vertex displacement to minimum\n"+
  "              adjacent face edge length is NUM. Vertex displacements\n"+
  "              are scaled down when necessary to conform to NUM.\n"+
  "              Default is '" + d2str(MAX_ACTUAL_DISPL_FRACTION) + "'.\n\n"+
  "       -K NUM, --max_runtime=NUM\n"+              
  "              End program execution when wall clock time exceeds NUM.\n"+
  "              After each vertex move, if elapsed time is greater than\n"+
  "              max runtime then output data files are written and\n"+
  "              the program exists. Default is to do nothing.\n\n"+
  "       -y STRING, --format_intersected_faces=STRING\n"+ 
  "              If STRING equals 'dreamm' then intersecting face\n"+
  "              output data is written in DReAMM custom points format.\n"+
  "              That is to say the location of each vertex of each\n"+
  "              intersecting face is written as 'x y z 1 0 0 1'.\n"+
  "              If STRING is other than 'dreamm',the detailed information about each\n"+
  "              intersecting face and its vertices are written.\n"+
  "              Default is '" + FORMAT_INTERSECTED_FACES + "'.\n\n"+
  "       -M STRING, --format_nonnice_vertices=STRING\n"+  
  "              If STRING equals 'dreamm' then nonnice vertices\n"+
  "              output data is written in DReAMM custom points format.\n"+
  "              That is to say the location of each vertex\n"+
  "              is written as 'x y z 1 0 0 1'.\n"+
  "              If STRING is other than 'dreamm',the detailed\n"+
  "              information about each nonnice vertex is written.\n"+
  "              Default is '" + FORMAT_NONNICE_VERTICES + "'.\n\n"+
  "       -i DIRECTORY, --input_data_dir=DIRECTORY\n"+            
  "              Directory containing input meshes.\n"+
  "              Default is current directory.\n\n"+
  "       -o DIRECTORY, --output_data_dir=DIRECTORY\n"+           
  "              Directory where output datawill be written.\n"+
  "              Default is current directory.\n\n"+
  "       -v FILE, --frozen_vertices_file=FILE\n"+     
  "              The vertices specified in FILE will not be moved and will\n"+
  "              retain their original position from the input meshes.\n"+
  "              The format of FILE must be 'object_name vertex_index'.\n"+
  "              For example, 'd000 134' indicates that vertex number 134\n"+
  "              in object d000 should be frozen.\n"+
  "              In the absence of this option, no vertices are frozen.\n\n"+
  "       -b FILE, --vertex_sequence_file=FILE\n"+     
  "              The vertices specified in FILE will be moved in the sequence given.\n"+
  "              The format of FILE must be 'object_name vertex_index'.\n"+
  "              For example, 'd000 134' indicates that vertex number 134\n"+
  "              in object d000 is to be moved.\n"+
  "              Default behavior in the case that '-f' option is not used\n"+
  "              is to not freeze any vertices.\n\n"+
  "       -N STRING, --mesh_output_suffix=STRING\n"+       
  "              STRING is incorporated into all output mesh data filenames\n"+
  "              and can therefore be used to distinguish input mesh files\n"+
  "              from output mesh files when input and output data directories\n"+
  "              are the same. Default is '" + MESH_OUTPUT_SUFFIX + "'.\n\n"+
  "       -P FILE, --main_log_file=FILE\n"+            
  "              A log of program progress through main function and elapsed\n"+
  "              times are recorded in FILE. Default is '" + MAIN_LOG_FILE + "'.\n\n"+
  "       -T FILE, --cont_log_file=FILE\n"+            
  "              A log of morphing progress is written to FILE\n"+
  "              after completion of every group. Default is '" + CONT_LOG_FILE + "'.\n\n"+
  "       -U FILE, --sep_log_file=FILE\n"+             
  "              A measure of ecw in model with surface sampling density\n"+
  "              defined by --ecw_sampling_length is written to FILE\n"+
  "              in output directory specified by -o unless --no_ecw is set.\n"+
  "              Default is '" + SEP_LOG_FILE + "'.\n\n"+
  "       -V FILE, --object_list_file=FILE\n"+         
  "              A summary of all objects found in data input directory\n"+
  "              are written to FILE. Default is '" + OBJECT_LIST_FILE + "'.\n\n"+
  "       -W FILE, --control_file=FILE\n"+             
  "              A summary of all meshmorph command option settings\n"+
  "              are written to FILE. Default is '" + CONT_LOG_FILE + "'.\n\n"+
  "       -X FILE, --vertex_selection_file=FILE\n"+    
  "              A list of the number of times each vertex in the model\n"+
  "              was moved is written to FILE after completion of every\n"+
  "              group of moved vertices unless --no_every_group is set\n"+
  "              in which case cumulative vertex move data for all groups\n"+
  "              is written once after completion of last group of vertex moves.\n"
  "              The identity of the vertices is not written just the\n"+
  "              number of moves. The ordering of the vertices is not\n"+
  "              necessarily alphabetical by object since it depends\n"+
  "              on dirent.h library. Default is '" + VERTEX_SELECTION_FILE + "'.\n\n"+
  "       -Y FILE, --refracted_file=FILE\n"+           
  "              A list of vertices that were refracted at any time for any reason\n"+
  "              is written to FILE after completion of every\n"+
  "              group of moved vertices unless --no_every_group is set\n"+
  "              in which case cumulative vertex refraction data for all groups\n"+
  "              is written once after completion of last group of vertex moves.\n"
  "              Default is '" + REFRACTED_FILE + "'.\n\n"+
  "       -Z FILE, --intersected_file=FILE\n"+         
  "              A list of faces that are currently intersected\n"+
  "              is written to FILE after completion of every\n"+
  "              group of moved vertices unless --no_every_group is set\n"+
  "              in which case the list of intersected faces\n"+
  "              is written once after completion of last group of vertex moves.\n"
  "              Default is '" + INTERSECTED_FILE + "'.\n\n"+
  "       -z FILE, --nonnice_file=FILE\n"+             
  "              A list of vertices that are currently nonnice, i.e. vertices located inside an object,\n"+
  "              is written to FILE after completion of every\n"+
  "              group of moved vertices unless --no_every_group is set\n"+
  "              in which case the list of nonnice vertices\n"+
  "              is written once after completion of last group of vertex moves.\n"
  "              Default is '" + NONNICE_FILE + "'.\n\n"+
  "       -x NUM, --ecw_sampling_length=NUM\n"+      
  "              Unless --no_ecw is set, meshmorph will measure and write\n"+
  "              to file the extracellular width of the model with surface\n"+
  "              sampling density NUM and to write data to file\n"+
  "              specified by --sep_log_file in output directory specified by -o.\n"+
  "              Default is '" + d2str(ECW_SAMPLING_LENGTH) + "'.\n\n"+
  "       -h, --help\n"+                     
  "              Print meshmorph man page.\n"+
  "\nJustin Kinney				2008/05/01\n";
  return message;
}

/** Parse meshmorph command line.
 * \param[in] argc Argc.
 * \param[in] argv Argv.
 * \param[in] message General error/help message explaining meshmorph.
 * \return True if frozen vertex file name detected; false otherwise.
 */

void Controls::parseCommandLine (int argc,char **argv,std::string const & message)
{

  while (1)
  {
    static struct option long_options[] =
    {
      {"verbose_init"                  , no_argument, & WRITE_VERBOSE_INIT                  , 1}, // default 0
      {"no_refracted"                  , no_argument, & WRITE_REFRACTED_VERTICES_TO_FILE    , 0}, // default 1
      {"no_intersected"                , no_argument, & WRITE_INTERSECTED_FACES_TO_FILE     , 0}, // default 1
      {"no_nonnice"                    , no_argument, & WRITE_NONNICE_VERTICES_TO_FILE      , 0}, // default 1
      {"no_ecw"                        , no_argument, & WRITE_ECW_TO_FILE                   , 0}, // default 1
      {"no_every_group"                , no_argument, & WRITE_EVERY_GROUP                   , 0}, // default 1
      {"move_histo"                    , no_argument, & WRITE_VERTEX_MOVE_HISTOGRAM         , 1}, // default 0
      {"append_group"                  , no_argument, & APPEND_GROUP_NUMBER                 , 1}, // default 0
      {"no_int_prevention"             , no_argument, & STRICT_FACE_INTERSECTION_PREVENTION , 0}, // default 1
      {"ref_orig_edge_length"          , no_argument, & USE_EDGE_REFERENCE_LENGTH           , 1}, // default 0
      {"verbose_vertex_move"           , no_argument, & ENABLE_VTRACK                       , 1}, // default 0
      {"disable_gain_scheduling"       , no_argument, & DISABLE_GAIN_SCHEDULING             , 1}, // default 0
      {"disable_messages"              , no_argument, & DISABLE_MESSAGES                    , 1}, // default 0
      {"measure_ecw_and_exit"          , no_argument, & MEASURE_ECW_AND_EXIT                , 1}, // default 0
      {"max_items_per_leaf"            , required_argument, 0, '1'},
      {"max_octree_depth"              , required_argument, 0, '2'},
      {"max_filename_size"             , required_argument, 0, '3'},
      {"update_region_size"            , required_argument, 0, '4'},
      {"seed_region_size"              , required_argument, 0, '5'},
      {"num_radius_steps"              , required_argument, 0, '6'},
      {"group_size"                    , required_argument, 0, 'g'},
      {"num_groups"                    , required_argument, 0, 'n'},
      {"refractory_period"             , required_argument, 0, '7'},
      {"max_touches"                   , required_argument, 0, '8'},
      {"vertex_move_print_period"      , required_argument, 0, '9'},
      {"begin_short_print_period"      , required_argument, 0, 'Q'},
      {"write_mesh_now"                , required_argument, 0, 'R'},
      {"vector_reserve"                , required_argument, 0, 'B'},
      {"min_cell_size"                 , required_argument, 0, 'C'},
      {"small_ecw_threshold"           , required_argument, 0, 'D'},
      {"edge_angle_threshold"          , required_argument, 0, 'e'},
      {"min_search_cone_radius"        , required_argument, 0, 'F'},
      {"search_radius"                 , required_argument, 0, 'G'},
      {"search_cone_angle"             , required_argument, 0, 'H'},
      {"face_intersection_weight"      , required_argument, 0, 'f'},
      {"edge_length_weight"            , required_argument, 0, 'l'},
      {"ecw_weight"                    , required_argument, 0, 's'},
      {"edge_angle_weight"             , required_argument, 0, 'a'},
      {"ecw_gain"                      , required_argument, 0, 'S'},
      {"edge_angle_gain"               , required_argument, 0, 'A'},
      {"edge_stretch_gain"             , required_argument, 0, 'L'},
      {"aspect_ratio_gain"             , required_argument, 0, 'E'},
      {"aspect_ratio_threshold"        , required_argument, 0, 'c'},
      {"overall_gain"                  , required_argument, 0, 'O'},
      {"min_vertex_disp"               , required_argument, 0, 'I'},
      {"target_ecw"                    , required_argument, 0, 't'},
      {"ecw_threshold"                 , required_argument, 0, 'd'},
      {"target_ecw_high"               , required_argument, 0, 'j'},
      {"target_ecw_low"                , required_argument, 0, 'k'},
      {"gain_step"                     , required_argument, 0, 'm'},
      {"max_displ_fraction"            , required_argument, 0, 'J'},
      {"max_runtime"                   , required_argument, 0, 'K'},
      {"format_intersected_faces"      , required_argument, 0, 'y'},
      {"format_nonnice_vertices"       , required_argument, 0, 'M'},
      {"input_data_dir"                , required_argument, 0, 'i'},
      {"output_data_dir"               , required_argument, 0, 'o'},
      {"frozen_vertices_file"          , required_argument, 0, 'v'},
      {"vertex_sequence_file"          , required_argument, 0, 'b'},
      {"mesh_output_suffix"            , required_argument, 0, 'N'},
      {"main_log_file"                 , required_argument, 0, 'P'},
      {"cont_log_file"                 , required_argument, 0, 'T'},
      {"sep_log_file"                  , required_argument, 0, 'U'},
      {"object_list_file"              , required_argument, 0, 'V'},
      {"control_file"                  , required_argument, 0, 'W'},
      {"vertex_selection_file"         , required_argument, 0, 'X'},
      {"refracted_file"                , required_argument, 0, 'Y'},
      {"intersected_file"              , required_argument, 0, 'Z'},
      {"nonnice_file"                  , required_argument, 0, 'z'},
      {"ecw_sampling_length"           , required_argument, 0, 'x'},
      {"help"                          , no_argument      , 0, 'h'},
      {0, 0, 0, 0}
    };
    /* getopt_long stores the option index here. */
    int option_index = 0;

    int c = getopt_long (argc, argv, "1:2:3:4:5:6:7:8:9:A:B:C:D:E:F:G:H:I:J:K:L:M:N:O:P:Q:R:S:T:U:V:W:X:Y:Z:a:b:c:d:e:f:g:hi:j:k:l:m:n:o:s:t:v:x:y:z:",
                     long_options, &option_index);

    /* Detect the end of the options. */
    if (c == -1)
      break;

    char *eptr=NULL;
    std::string line;
    switch (c)
    {
      case 0:
        /* If this option set a flag, do nothing else now. */
        if (long_options[option_index].flag != 0)
          break;

      case '1':
        MAX_ITEMS_PER_LEAF = static_cast<int>(strtod(optarg,&eptr));
        break;

      case '2':
        MAX_OCTREE_DEPTH = static_cast<int>(strtod(optarg,&eptr));
        break;

      case '3':
        MAX_FILENAME_SIZE = static_cast<int>(strtod(optarg,&eptr));
        break;

      case '4':
        UPDATE_REGION_SIZE = strtod(optarg,&eptr);
        break;

      case '5':
        SEED_REGION_SIZE = strtod(optarg,&eptr);
        break;

      case '6':
        NUMBER_RADIUS_STEPS = static_cast<int>(strtod(optarg,&eptr));
        break;

      case 'g':
        GROUP_SIZE = static_cast<int>(strtod(optarg,&eptr));
        REDEFINE_GROUP_SIZE = 0;
        break;

      case 'n':
        NUM_GROUPS = static_cast<int>(strtod(optarg,&eptr));
        REDEFINE_NUM_GROUPS = 0;
        break;

      case '7':
        REFRACTORY_PERIOD = static_cast<int>(strtod(optarg,&eptr));
        break;

      case '8':
        MAX_TOUCHES = static_cast<int>(strtod(optarg,&eptr));
        break;

      case '9':
        PRINT_PERIOD = static_cast<int>(strtod(optarg,&eptr));
        break;

      case 'Q':
        BEGIN_SHORT_PRINT_PERIOD = static_cast<int>(strtod(optarg,&eptr));
        break;

      case 'R':
        WRITE_MESH_NOW = static_cast<int>(strtod(optarg,&eptr));
        break;

      case 'B':
        VECTOR_RESERVE = static_cast<int>(strtod(optarg,&eptr));
        break;

      case 'C':
        MIN_CELL_SIZE = strtod(optarg,&eptr);
        break;

      case 'D':
        SMALL_ECW_THRESHOLD = strtod(optarg,&eptr);
        break;

      case 'F':
        MIN_SEARCH_CONE_RADIUS = strtod(optarg,&eptr);
        break;

      case 'G':
        SEARCH_RADIUS_SQ = strtod(optarg,&eptr)*strtod(optarg,&eptr);
        break;

      case 'H':
        CLOSEST_POINT_ANGLE  = strtod(optarg,&eptr)*PI/180.0;
        CLOSEST_POINT_COSINE = cos(CLOSEST_POINT_ANGLE);
        CLOSEST_POINT_SINE   = sin(CLOSEST_POINT_ANGLE);
        break;

      case 'I':
        MIN_DISPLACEMENT_SQ = strtod(optarg,&eptr)*strtod(optarg,&eptr);
        break;

      case 'J':
        MAX_ACTUAL_DISPL_FRACTION= strtod(optarg,&eptr);
        break;

      case 'K':
        MAX_RUNTIME = strtod(optarg,&eptr);
        break;

      case 'y':
        FORMAT_INTERSECTED_FACES = optarg;
        break;

      case 'M':
        FORMAT_NONNICE_VERTICES = optarg;
        break;

      case 'N':
        MESH_OUTPUT_SUFFIX = optarg;
        break;

      case 'P':
        MAIN_LOG_FILE = optarg;
        break;

      case 'T':
        CONT_LOG_FILE = optarg;
        break;

      case 'U':
        SEP_LOG_FILE = optarg;
        break;

      case 'V':
        OBJECT_LIST_FILE = optarg;
        break;

      case 'W':
        CONTROL_FILE = optarg;
        break;

      case 'X':
        VERTEX_SELECTION_FILE = optarg;
        break;

      case 'Y':
        REFRACTED_FILE = optarg;
        break;

      case 'Z':
        INTERSECTED_FILE = optarg;
        break;

      case 'z':
        NONNICE_FILE = optarg;
        break;

      case 'e':
        EDGE_ANGLE_THRESHOLD = strtod(optarg,&eptr);
        break;

      case 'h':
        std::cout << message << std::endl;
        exit(1);
        break;

      case 'f':
        INTERSECTION_WEIGHT = strtod(optarg,&eptr);
        break;

      case 'l':
        EDGE_LENGTH_WEIGHT = strtod(optarg,&eptr);
        break;

      case 's':
        ECW_WEIGHT = strtod(optarg,&eptr);
        //EDGE_ANGLE_GAIN = ECW_WEIGHT*100.0;
        break;

      case 'a':
        EDGE_ANGLE_WEIGHT = strtod(optarg,&eptr);
        break;

      case 'S':
        ECW_GAIN = strtod(optarg,&eptr);
        break;

      case 'A':
        EDGE_ANGLE_GAIN = strtod(optarg,&eptr);
        break;

      case 'L':
        EDGE_LENGTH_GAIN = strtod(optarg,&eptr);
        break;

      case 'E':
        ASPECT_RATIO_GAIN = strtod(optarg,&eptr);
        break;

      case 'c':
        ASPECT_RATIO_THRESHOLD = strtod(optarg,&eptr);
        break;

      case 'O':
        OVERALL_GAIN = strtod(optarg,&eptr);
        break;

      case 't':
        TARGET_ECW = strtod(optarg,&eptr);
        break;

      case 'd':
        ECW_THRESHOLD = strtod(optarg,&eptr);
        break;

      case 'j':
        TARGET_ECW_HIGH = strtod(optarg,&eptr);
        break;

      case 'k':
        TARGET_ECW_LOW = strtod(optarg,&eptr);
        break;

      case 'm':
        GAIN_STEP = strtod(optarg,&eptr);
        break;

      case 'i':
        INPUT_DATA_DIR = processDir(optarg);
        break;

      case 'o':
        OUTPUT_DATA_DIR = processDir(optarg);
        break;

      case 'v':
        FROZEN_VERTICES_FILE = optarg;
        break;

      case 'b':
        VERTEX_SEQUENCE_FILE = optarg;
        break;

      case 'x':
        ECW_SAMPLING_LENGTH = strtod(optarg,&eptr);
        break;

      case '?':
        /* getopt_long already printed an error message. */
        exit(1);
        break;

      default:
        std::cout << message << std::endl;
        abort ();
    }
  }

  /* Print any remaining command line arguments (not options). */
  if (optind < argc)
  {
    printf ("non-option ARGV-elements: ");
    while (optind < argc)
      printf ("%s ", argv[optind++]);
    putchar ('\n');
    std::cout << message << std::endl;
    exit(1);
  }
}

/** Get summary of meshmorph command parameter settings.
 * \return Command parameter settings.
 */

std::string Controls::getCommandSettings (void)
{
 std::string message = "";
  message =message +
  "WRITE_VERBOSE_INIT = "                  + i2str(WRITE_VERBOSE_INIT)                  + "\n" +
  "WRITE_REFRACTED_VERTICES_TO_FILE = "    + i2str(WRITE_REFRACTED_VERTICES_TO_FILE)    + "\n" +
  "WRITE_INTERSECTED_FACES_TO_FILE = "     + i2str(WRITE_INTERSECTED_FACES_TO_FILE)     + "\n" + 
  "WRITE_NONNICE_VERTICES_TO_FILE = "      + i2str(WRITE_NONNICE_VERTICES_TO_FILE)      + "\n" + 
  "WRITE_ECW_TO_FILE = "                   + i2str(WRITE_ECW_TO_FILE)                   + "\n" + 
  "WRITE_EVERY_GROUP = "                   + i2str(WRITE_EVERY_GROUP)                   + "\n" + 
  "WRITE_VERTEX_MOVE_HISTOGRAM = "         + i2str(WRITE_VERTEX_MOVE_HISTOGRAM)         + "\n" + 
  "WRITE_MESH_NOW = "                      + i2str(WRITE_MESH_NOW)                      + "\n" + 
  "APPEND_GROUP_NUMBER = "                 + i2str(APPEND_GROUP_NUMBER)                 + "\n" + 
  "STRICT_FACE_INTERSECTION_PREVENTION = " + i2str(STRICT_FACE_INTERSECTION_PREVENTION) + "\n" + 
  "USE_EDGE_REFERENCE_LENGTH = "           + i2str(USE_EDGE_REFERENCE_LENGTH)           + "\n" + 
  "ENABLE_VTRACK = "                       + i2str(ENABLE_VTRACK)                       + "\n" + 
  "DISABLE_GAIN_SCHEDULING = "             + i2str(DISABLE_GAIN_SCHEDULING)             + "\n" + 
  "DISABLE_MESSAGES = "                    + i2str(DISABLE_MESSAGES)                    + "\n" +
  "MEASURE_ECW_AND_EXIT = "                + i2str(MEASURE_ECW_AND_EXIT)                + "\n" +
  "MAX_ITEMS_PER_LEAF = "                  + i2str(MAX_ITEMS_PER_LEAF)                  + "\n" +
  "MAX_OCTREE_DEPTH = "                    + i2str(MAX_OCTREE_DEPTH)                    + "\n" +
  "MAX_FILENAME_SIZE = "                   + i2str(MAX_FILENAME_SIZE)                   + "\n" +
  "OCTREE_MIN_X = "                        + d2str(OCTREE_MIN_X)                        + "\n" +
  "OCTREE_MIN_Y = "                        + d2str(OCTREE_MIN_Y)                        + "\n" +
  "OCTREE_MIN_Z = "                        + d2str(OCTREE_MIN_Z)                        + "\n" +
  "OCTREE_WIDTH = "                        + d2str(OCTREE_WIDTH)                        + "\n" +
  "UPDATE_REGION_SIZE = "                  + d2str(UPDATE_REGION_SIZE)                  + "\n" +
  "SEED_REGION_SIZE = "                    + d2str(SEED_REGION_SIZE)                    + "\n" +
  "NUMBER_RADIUS_STEPS = "                 + i2str(NUMBER_RADIUS_STEPS)                 + "\n" +
  "GROUP_SIZE = "                          + i2str(GROUP_SIZE)                          + "\n" +
  "REDEFINE_GROUP_SIZE = "                 + i2str(REDEFINE_GROUP_SIZE)                 + "\n" +
  "NUM_GROUPS = "                          + i2str(NUM_GROUPS)                          + "\n" +
  "REDEFINE_NUM_GROUPS = "                 + i2str(REDEFINE_NUM_GROUPS)                 + "\n" +
  "REFRACTORY_PERIOD = "                   + i2str(REFRACTORY_PERIOD)                   + "\n" +
  "MAX_TOUCHES = "                         + i2str(MAX_TOUCHES)                         + "\n" +
  "PRINT_PERIOD = "                        + i2str(PRINT_PERIOD)                        + "\n" +
  "BEGIN_SHORT_PRINT_PERIOD = "            + i2str(BEGIN_SHORT_PRINT_PERIOD)            + "\n" +
  "VECTOR_RESERVE = "                      + i2str(VECTOR_RESERVE)                      + "\n" +
  "MIN_CELL_SIZE = "                       + d2str(MIN_CELL_SIZE)                       + "\n" +
  "EDGE_ANGLE_THRESHOLD = "                + d2str(EDGE_ANGLE_THRESHOLD)                + "\n" +
  "MIN_SEARCH_CONE_RADIUS = "              + d2str(MIN_SEARCH_CONE_RADIUS)              + "\n" +
  "SEARCH_RADIUS_SQ = "                    + d2str(SEARCH_RADIUS_SQ)                    + "\n" +
  "CLOSEST_POINT_ANGLE  = "                + d2str(CLOSEST_POINT_ANGLE*180.0/PI)        + "\n" +
  "CLOSEST_POINT_COSINE = "                + d2str(CLOSEST_POINT_COSINE)                + "\n" +
  "CLOSEST_POINT_SINE = "                  + d2str(CLOSEST_POINT_SINE)                  + "\n" +
  "INTERSECTION_WEIGHT = "                 + d2str(INTERSECTION_WEIGHT)                 + "\n" +
  "EDGE_LENGTH_WEIGHT = "                  + d2str(EDGE_LENGTH_WEIGHT)                  + "\n" +
  "ECW_WEIGHT = "                          + d2str(ECW_WEIGHT)                          + "\n" +
  "EDGE_ANGLE_WEIGHT = "                   + d2str(EDGE_ANGLE_WEIGHT)                   + "\n" +
  "ECW_GAIN = "                            + d2str(ECW_GAIN)                            + "\n" +
  "EDGE_ANGLE_GAIN = "                     + d2str(EDGE_ANGLE_GAIN)                     + "\n" +
  "EDGE_LENGTH_GAIN = "                    + d2str(EDGE_LENGTH_GAIN)                    + "\n" +
  "ASPECT_RATIO_GAIN = "                   + d2str(ASPECT_RATIO_GAIN)                   + "\n" +
  "ASPECT_RATIO_THRESHOLD = "              + d2str(ASPECT_RATIO_THRESHOLD)              + "\n" +
  "OVERALL_GAIN = "                        + d2str(OVERALL_GAIN)                        + "\n" +
  "MIN_DISPLACEMENT_SQ = "                 + d2str(MIN_DISPLACEMENT_SQ)                 + "\n" +
  "TARGET_ECW = "                          + d2str(TARGET_ECW)                          + "\n" +
  "ECW_THRESHOLD = "                       + d2str(ECW_THRESHOLD)                       + "\n" +
  "TARGET_ECW_HIGH = "                     + d2str(TARGET_ECW_HIGH)                     + "\n" +
  "TARGET_ECW_LOW = "                      + d2str(TARGET_ECW_LOW)                      + "\n" +
  "GAIN_STEP = "                           + d2str(GAIN_STEP)                           + "\n" +
  "MAX_ACTUAL_DISPL_FRACTION = "           + d2str(MAX_ACTUAL_DISPL_FRACTION)           + "\n" +
  "FORMAT_INTERSECTED_FACES = "            + FORMAT_INTERSECTED_FACES                   + "\n" + 
  "FORMAT_NONNICE_VERTICES = "             + FORMAT_NONNICE_VERTICES                    + "\n" + 
  "ECW_SAMPLING_LENGTH = "                 + d2str(ECW_SAMPLING_LENGTH)                 + "\n" + 
  "MY_DOUBLE_EPSILON = "                   + d2str(MY_DOUBLE_EPSILON)                   + "\n" + 
  "EPSILON = "                             + d2str(EPSILON)                             + "\n" + 
  "ENERGY_WINDOW = "                       + i2str(ENERGY_WINDOW)                       + "\n" + 
  "ENERGY_SAMPLE_PERIOD = "                + i2str(ENERGY_SAMPLE_PERIOD)                + "\n" + 
  "PI = "                                  + d2str(PI)                                  + "\n" + 
  "MAX_RUNTIME = "                         + d2str(MAX_RUNTIME)                         + "\n" + 
  "INPUT_DATA_DIR = "                      + INPUT_DATA_DIR                             + "\n" + 
  "OUTPUT_DATA_DIR = "                     + OUTPUT_DATA_DIR                            + "\n" + 
  "FROZEN_VERTICES_FILE = "                + FROZEN_VERTICES_FILE                       + "\n" + 
  "VERTEX_SEQUENCE_FILE = "                + VERTEX_SEQUENCE_FILE                       + "\n" + 
  "MESH_OUTPUT_SUFFIX = "                  + MESH_OUTPUT_SUFFIX                         + "\n" + 
  "MAIN_LOG_FILE = "                       + MAIN_LOG_FILE                              + "\n" + 
  "CONT_LOG_FILE = "                       + CONT_LOG_FILE                              + "\n" + 
  "SEP_LOG_FILE = "                        + SEP_LOG_FILE                               + "\n" + 
  "OBJECT_LIST_FILE = "                    + OBJECT_LIST_FILE                           + "\n" + 
  "CONTROL_FILE = "                        + CONTROL_FILE                               + "\n" + 
  "VERTEX_SELECTION_FILE = "               + VERTEX_SELECTION_FILE                      + "\n" + 
  "REFRACTED_FILE = "                      + REFRACTED_FILE                             + "\n" + 
  "INTERSECTED_FILE = "                    + INTERSECTED_FILE                           + "\n" + 
  "NONNICE_FILE = "                        + NONNICE_FILE                               + "\n" +
  "SMALL_ECW_THRESHOLD = "                 + d2str(SMALL_ECW_THRESHOLD)                 + "\n"; 
  return message;
}

