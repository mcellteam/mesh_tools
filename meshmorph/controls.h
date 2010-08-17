// Author: Justin Kinney
// Date: Sep 2008


#ifndef CONTROLS_H
#define CONTROLS_H 1

#include <string>

class Controls
{
public:
  static Controls & instance     (void);
  std::string processDir         (char *);
  void        parseCommandLine   (int argc,char **argv,std::string const & message);
  void        redefineGroup      (void);
  void        recordOctreeSize   (double size);
  std::string getUsageMessage    (void);
  std::string getCommandSettings (void);
  std::string d2str              (double const & i);
  std::string i2str              (int const & i);

  void        set_strict_face_intersection_prevention () throw() { STRICT_FACE_INTERSECTION_PREVENTION=true; }

  int    get_freeze_sheets                        () const throw() { return FREEZE_SHEETS; }
  int    get_freeze_tunnels                       () const throw() { return FREEZE_TUNNELS; }
  int    get_dual_vertex_ecws                     () const throw() { return DUAL_TARGET_ECWS; }
  int    get_vertex_neighbor_count_period         () const throw() { return VERTEX_NEIGHBOR_COUNT_PERIOD; }
  int    get_report_vertex_identity               () const throw() { return REPORT_VERTEX_IDENTITY; }
  int    get_report_vertex_area                   () const throw() { return REPORT_VERTEX_AREA; }
  int    get_count_neighbors                      () const throw() { return COUNT_NEIGHBORS; }
  int    get_assume_nice_vertices                 () const throw() { return ASSUME_NICE_VERTICES; }
  int    get_measure_ecw_and_exit                 () const throw() { return MEASURE_ECW_AND_EXIT; }
  int    get_curvature_neighborhood_size          () const throw() { return CURVATURE_NEIGHBORHOOD_SIZE; }
  int    get_max_items_per_leaf                   () const throw() { return MAX_ITEMS_PER_LEAF; }
  int    get_max_octree_depth                     () const throw() { return MAX_OCTREE_DEPTH; }
  int    get_max_filename_size                    () const throw() { return MAX_FILENAME_SIZE; }
  int    get_number_radius_steps                  () const throw() { return NUMBER_RADIUS_STEPS; }
  int    get_group_size                           () const throw() { return GROUP_SIZE; }
  int    get_num_groups                           () const throw() { return NUM_GROUPS; }
  int    get_refractory_period                    () const throw() { return REFRACTORY_PERIOD; }
  int    get_max_touches                          () const throw() { return MAX_TOUCHES; }
  int    get_energy_sample_period                 () const throw() { return ENERGY_SAMPLE_PERIOD; }
  int    get_print_period                         () const throw() { return PRINT_PERIOD; }
  int    get_begin_short_print_period             () const throw() { return BEGIN_SHORT_PRINT_PERIOD; }
  int    get_write_mesh_now                       () const throw() { return WRITE_MESH_NOW; }
  int    get_vector_reserve                       () const throw() { return VECTOR_RESERVE; }
  int    get_write_verbose_init                   () const throw() { return WRITE_VERBOSE_INIT; }
  int    get_write_refracted_vertices_to_file     () const throw() { return WRITE_REFRACTED_VERTICES_TO_FILE; }
  int    get_write_intersected_faces_to_file      () const throw() { return WRITE_INTERSECTED_FACES_TO_FILE; }
  int    get_write_nonnice_vertices_to_file       () const throw() { return WRITE_NONNICE_VERTICES_TO_FILE; }
  int    get_write_every_group                    () const throw() { return WRITE_EVERY_GROUP; }
  int    get_write_ecw_to_file                    () const throw() { return WRITE_ECW_TO_FILE; }
  int    get_write_vertex_move_histogram          () const throw() { return WRITE_VERTEX_MOVE_HISTOGRAM; }
  int    get_write_closest_points                 () const throw() { return WRITE_CLOSEST_POINTS; }
  int    get_assume_max_ecw_error                 () const throw() { return ASSUME_MAX_ECW_ERROR; }
  int    get_append_group_number                  () const throw() { return APPEND_GROUP_NUMBER; }
  int    get_strict_face_intersection_prevention  () const throw() { return STRICT_FACE_INTERSECTION_PREVENTION; }
  int    get_use_edge_reference_length            () const throw() { return USE_EDGE_REFERENCE_LENGTH; }
  int    get_enable_vtrack                        () const throw() { return ENABLE_VTRACK; }
  int    get_disable_gain_scheduling              () const throw() { return DISABLE_GAIN_SCHEDULING; }
  int    get_disable_messages                     () const throw() { return DISABLE_MESSAGES; }
  double get_max_radius_of_curvature              () const throw() { return MAX_RADIUS_OF_CURVATURE; }
  double get_region_orthogonality_threshold       () const throw() { return REGION_ORTHOGONALITY_THRESHOLD; }
  double get_octree_min_x                         () const throw() { return OCTREE_MIN_X; }
  double get_octree_min_y                         () const throw() { return OCTREE_MIN_Y; }
  double get_octree_min_z                         () const throw() { return OCTREE_MIN_Z; }
  double get_octree_width                         () const throw() { return OCTREE_WIDTH; }
  double get_min_cell_size                        () const throw() { return MIN_CELL_SIZE; }
  double get_edge_angle_threshold                 () const throw() { return EDGE_ANGLE_THRESHOLD; }
  double get_my_double_epsilon                    () const throw() { return MY_DOUBLE_EPSILON; }
  double get_epsilon                              () const throw() { return EPSILON; }
  double get_update_region_size                   () const throw() { return UPDATE_REGION_SIZE; }
  double get_seed_region_size                     () const throw() { return SEED_REGION_SIZE; }
  double get_min_search_cone_radius               () const throw() { return MIN_SEARCH_CONE_RADIUS; }
  double get_search_radius_sq                     () const throw() { return SEARCH_RADIUS_SQ; }
  double get_closest_point_angle                  () const throw() { return CLOSEST_POINT_ANGLE; }
  double get_closest_point_cosine                 () const throw() { return CLOSEST_POINT_COSINE; }
  double get_closest_point_sine                   () const throw() { return CLOSEST_POINT_SINE; }
  double get_intersection_weight                  () const throw() { return INTERSECTION_WEIGHT; }
  double get_edge_length_weight                   () const throw() { return EDGE_LENGTH_WEIGHT; }
  double get_ecw_weight                           () const throw() { return ECW_WEIGHT ; }
  double get_edge_angle_weight                    () const throw() { return EDGE_ANGLE_WEIGHT; }
  double get_ecw_gain                             () const throw() { return ECW_GAIN; }
  double get_edge_angle_gain                      () const throw() { return EDGE_ANGLE_GAIN; }
  double get_edge_length_gain                     () const throw() { return EDGE_LENGTH_GAIN; }
  double get_aspect_ratio_gain                    () const throw() { return ASPECT_RATIO_GAIN; }
  double get_aspect_ratio_threshold               () const throw() { return ASPECT_RATIO_THRESHOLD; }
  double get_overall_gain                         () const throw() { return OVERALL_GAIN; }
  double get_gain_step                            () const throw() { return GAIN_STEP; }
  double get_min_displacement_sq                  () const throw() { return MIN_DISPLACEMENT_SQ; }
  double get_pi                                   () const throw() { return PI; }
  double get_target_ecw                           () const throw() { return TARGET_ECW; }
//  double get_ecw_threshold                        () const throw() { return ECW_THRESHOLD; }
  double get_target_ecw_high                      () const throw() { return TARGET_ECW_HIGH; }
  double get_target_ecw_low                       () const throw() { return TARGET_ECW_LOW; }
  double get_max_actual_displ_fraction            () const throw() { return MAX_ACTUAL_DISPL_FRACTION; }
  double get_max_runtime                          () const throw() { return MAX_RUNTIME; }
  double get_ecw_sampling_length                  () const throw() { return ECW_SAMPLING_LENGTH; }
  double get_small_ecw_threshold                  () const throw() { return SMALL_ECW_THRESHOLD; }
  std::string get_format_intersected_faces        () const throw() { return FORMAT_INTERSECTED_FACES; }
  std::string get_format_nonnice_vertices         () const throw() { return FORMAT_NONNICE_VERTICES; }
  std::string get_input_data_dir                  () const throw() { return INPUT_DATA_DIR; }
  std::string get_output_data_dir                 () const throw() { return OUTPUT_DATA_DIR; }
  std::string get_frozen_vertices_file            () const throw() { return FROZEN_VERTICES_FILE; }
  std::string get_cleft_vertices_file             () const throw() { return CLEFT_VERTICES_FILE; }
  std::string get_vertex_sequence_file            () const throw() { return VERTEX_SEQUENCE_FILE; }
  std::string get_mesh_output_suffix              () const throw() { return MESH_OUTPUT_SUFFIX; }
  std::string get_main_log_file                   () const throw() { return MAIN_LOG_FILE; }
  std::string get_cont_log_file                   () const throw() { return CONT_LOG_FILE; }
  std::string get_sep_log_file                    () const throw() { return SEP_LOG_FILE; }
  std::string get_object_list_file                () const throw() { return OBJECT_LIST_FILE; }
  std::string get_control_file                    () const throw() { return CONTROL_FILE; }
  std::string get_vertex_selection_file           () const throw() { return VERTEX_SELECTION_FILE; }
  std::string get_refracted_file                  () const throw() { return REFRACTED_FILE; }
  std::string get_intersected_file                () const throw() { return INTERSECTED_FILE; }
  std::string get_nonnice_file                    () const throw() { return NONNICE_FILE; }
  std::string get_gain_scheduling_file            () const throw() { return GAIN_SCHEDULING_FILE; }
  std::string get_recon_block_wrap                () const throw() { return RECON_BLOCK_WRAP; }
  void updatePrintPeriod (int count);

private:
  static Controls * only_one;
  Controls                (void);
  Controls                (Controls const &);
  Controls & operator =   (Controls const &);

private:
  // octree parameters
  int MAX_ITEMS_PER_LEAF;
  int MAX_OCTREE_DEPTH;
  double MIN_CELL_SIZE;
  double OCTREE_MIN_X;
  double OCTREE_MIN_Y;
  double OCTREE_MIN_Z;
  double OCTREE_WIDTH;

  double PI;
  // maximum file name size
  int MAX_FILENAME_SIZE;

  // Directs the model to be written to file at a
  // user-specified iteration number each group.
  // WRITE_MESH_NOW is the number of vertex moves in group after which
  // the entire mesh model will be written to file once.
  // Minimum value == 1 which will write mesh after one vertex move.
  // Default is '0' which does nothing.
  int WRITE_MESH_NOW;

  // set to true to write initialization information
  //		for each file to stdout, i.e. verbose setting
  // set to false for concise initialization informatin, i.e. concise setting
  int WRITE_VERBOSE_INIT;

  // set to true to write endpoints of extracellular width line segment
  //		when writing extracellular width to file
  // set to false to only write extracellular width to file
  int WRITE_CLOSEST_POINTS;

  // set to true to write refracted vertices to
  //		file specified by REFRACTED_FILE
  int WRITE_REFRACTED_VERTICES_TO_FILE;

  // set to true to write intersected faces to
  //		file specified by INTERSECTED_FILE
  int WRITE_INTERSECTED_FACES_TO_FILE;

  // set to true to assume maximum separation width error
  //            in the event no closest point is found for a vertex
  // set to false to assume zero separation width error
  //            in the event no closest point is found for a vertex
  int ASSUME_MAX_ECW_ERROR;

  // choice of output format for intersected faces
  //	dreamm = dreamm custom points format
  //	detail = face name and index + vertex details(3X)
  // choose one
  std::string FORMAT_INTERSECTED_FACES;

  // set to true to write nonnice vertices to
  //		file specified by NONNICE_FILE
  int WRITE_NONNICE_VERTICES_TO_FILE;

  // choice of output format for nonnice vertices
  //	dreamm = dreamm custom points format
  //	detail = vertex name,index,x,y,z
  // choose one
  std::string FORMAT_NONNICE_VERTICES;

  // set to true to write various output data
  // to file after completion of every group
  // set to false to only write output data
  // to file after the last group has finished
  int WRITE_EVERY_GROUP;

  // set to true to write closest point distance
  //		to file
  // set to false to do nothing
  int WRITE_ECW_TO_FILE;

  // The surface sampling density for extracellular width
  // measurement has an area equal
  // to an equilateral triangle of ECW_SAMPLING_LENGTH.
  double ECW_SAMPLING_LENGTH;

  // set to true to write vertex move histogram to file after each group
  // set to false to do nothing
  int WRITE_VERTEX_MOVE_HISTOGRAM;

  // set to true to use group number in output data filenames
  // set to false to not automatically use group number in output data filenames
  int APPEND_GROUP_NUMBER;

  // set to true to prevent any new intersections of faces
  // set to false to allow new intersecting faces
  //	(useful when faces exist with all three vertices nonnice)
  int STRICT_FACE_INTERSECTION_PREVENTION;

  // if moving a vertex creates an edge angle
  // less than EDGE_ANGLE_THRESHOLD,
  // then vertex is not moved
  double EDGE_ANGLE_THRESHOLD; //radians

  // for use with "is float close to zero?" 
  // in conditional statement
  double MY_DOUBLE_EPSILON;
  double EPSILON;

  // The region of space in which the model will be updated after a vertex move.
  // All vertices of all faces in octree cells that overlap this
  // region of space will have their closest point, extracellular
  // width, total force, virtual displacement, and virtual
  // displacement rank updated.
  double UPDATE_REGION_SIZE;

  // vertices are moved in sets consisting of a seed vertex
  // (typically the vertex with the largest virtual displacement)
  // and vertices of all faces in octree cells intersected by
  // cube centered on seed vertex with side length 
  // equal to 2*SEED_REGION_SIZE
  double SEED_REGION_SIZE;

  // If a candidate closest point is more than SEARCH_RADIUS 
  // distance away from given vertex, then the candidate closest 
  // point is disqualified from being a closest point to given vertex.
  double SEARCH_RADIUS_SQ; // nm * nm

  // The closest point to a vertex is found by first searching
  // a region defined by MIN_SEARCH_CONE_RADIUS. While no
  // closest point is found the search region is expanded up to
  // SEARCH_RADIUS_SQ in NUMBER_RADIUS_STEPS.
  int NUMBER_RADIUS_STEPS;       // must be greater than or equal to zero
  double MIN_SEARCH_CONE_RADIUS; // nm

  //  maximum allowed angle between closest point and vertex normal
  //  precompute cosine and sine of angle
  double CLOSEST_POINT_ANGLE; // radians
  double CLOSEST_POINT_COSINE;
  double CLOSEST_POINT_SINE;

  // if USE_EDGE_REFERENCE_LENGTH == true
  // then compare edge lengths to instantaneous 
  // mean edge length of adjacent faces to edge
  // else use original length of each edge
  int USE_EDGE_REFERENCE_LENGTH;

  // NOTE!! If USE_EDGE_REFERENCE_LENGTH==false
  // need to implemenet aspect_ratio detection
  // using the constants defined below
  // whereby edges associated with adjacent faces
  // with aspect ratios greater than ASPECT_RATIO_THRESHOLD
  // will experience a force to improve the aspect_ratio
  // proportional to ASPECT_GAIN
//  double FACE_ASPECT_RATIO_GAIN;

  // Force Function Weights
  // Weights are normalized during the force calculation,
  // so all that matters is the ratio of the four weights.
  // NOTE THESE MUST SUM TO 100!!!!!!!!!!!!!!!
  double INTERSECTION_WEIGHT; // nN
  double EDGE_LENGTH_WEIGHT;  // nN/nm
  double ECW_WEIGHT;          // nN/nm
  double EDGE_ANGLE_WEIGHT;   // nN

  // Force Function gains
  // The four principal forces are functions of
  // their associated error multiplied by a gain
  double ECW_GAIN;
  double EDGE_ANGLE_GAIN;
  double EDGE_LENGTH_GAIN;
  double ASPECT_RATIO_GAIN;
  double ASPECT_RATIO_THRESHOLD;

  // The displacement of each vertex is equal to
  // the product of the calculated force using the relative weights
  // from above and the overall gain.
  // It could be interesting to think of mapping from force to displacement
  // as motion of a massless particle in a viscous fluid, where
  // the overall gain is the time step (s) divided by damping (nN/(nm/s)).
  double OVERALL_GAIN;

  // number of vertex moves per group
  int GROUP_SIZE;
  
  // honor GROUP_SIZE specified from command line
  int REDEFINE_GROUP_SIZE;

  // Number of groups of GROUP_SIZE vertex moves to execute.
  int NUM_GROUPS;

  // honor NUM_GROUPS specified from command line
  int REDEFINE_NUM_GROUPS;

  // size of vertex move refractory period
  // i.e. the minimum number of vertex moves between moves of the same vertex
  int REFRACTORY_PERIOD;

  // maximum number of moves each vertex may make per group
  int MAX_TOUCHES;

  // energy averaging window size
  int ENERGY_WINDOW;

  // calculate and record total model energy
  // with the following period
  int ENERGY_SAMPLE_PERIOD;

  // minimum acceptable virtual displacment of vertex
  // if virtual_displacement*virtual_displacement < MIN_DISPLACEMENT_SQ
  //  then put vertex into refractory period
  double MIN_DISPLACEMENT_SQ;


  // desired distance between object surfaces
  double TARGET_ECW; // nm

  // If DUAL_TARGET_ECWS is specified, then vertices are classified
  // as either tunnel or sheet. Tunnels are morphed so as to have
  // an extracellular width of size TARGET_ECW_HIGH specified by -j
  // option. Sheet vertices will be morphed so as to have an extracellular
  // width of size TARGET_ECW_LOW specified by -k option.  Units are same
  // as meshes in input directory.  Default is to morph all vertices to
  // the same extracellular specified by -t option.
  int DUAL_TARGET_ECWS;
  // double ECW_THRESHOLD; // nm
  double TARGET_ECW_HIGH; // nm
  double TARGET_ECW_LOW; // nm

  // If gain scheduling is not disabled, then after
  // each group the maximum allwed gain is decremented
  // by NUM*overall_gain/num_groups. If NUM is positive
  // then the maximum gain decreases. If NUM is negative
  // then the maximum gain increases.
  double GAIN_STEP;

  // vertex single-move displacement shall be capped at MAX_ACTUAL_DISPL_SQ.
  // i.e. if (disp>MAX_ACTUAL_DISPL_SQ) then disp=MAX_ACTUAL_DISPL_SQ
  // max ratio of vertex displacement to minimum adjacent face edge length
  // max 5.7392 (asin(0.1)) degree change in edge angle
  // (although vertices on high aspect ratio faces could be moved a lot)
  double MAX_ACTUAL_DISPL_FRACTION;

  // maximum allowed runtime in seconds
  double MAX_RUNTIME;

  // number of vertex moves in cycle with disabled printing
  // alternating with a single vertex move logged to stdout
  // minimum value == 1 which will print each vertex move
  int PRINT_PERIOD;
  int BEGIN_SHORT_PRINT_PERIOD;

  // if true then print detailed information about moved vertex
  // else do nothing
  int ENABLE_VTRACK;

  int DISABLE_GAIN_SCHEDULING;

  // initial size of vectors to avoid vector resizing
  // and consequent data copying
  int VECTOR_RESERVE;

  std::string INPUT_DATA_DIR;
  std::string OUTPUT_DATA_DIR;
  std::string FROZEN_VERTICES_FILE;
  std::string VERTEX_SEQUENCE_FILE;

  std::string MESH_OUTPUT_SUFFIX;
  std::string MAIN_LOG_FILE;
  std::string CONT_LOG_FILE;
  std::string SEP_LOG_FILE;
  std::string OBJECT_LIST_FILE;
  std::string CONTROL_FILE;
  std::string VERTEX_SELECTION_FILE;
  std::string REFRACTED_FILE;
  std::string INTERSECTED_FILE;
  std::string NONNICE_FILE;
  
  // if defined by user, then update gain schedule
  // according to directives specified by user
  // if empty, allow gain scheduling to proceed
  // according to initial gain scheduling values
  std::string GAIN_SCHEDULING_FILE;
  
  // Small extracellular width for which closest point
  // can be outside of search cone
  double SMALL_ECW_THRESHOLD;

  // if true then do not print messages from refractory and virtual_displacement classes
  // else do nothing
  int DISABLE_MESSAGES;

  // if true, then measure and write to file
  // the extracellular width in the model
  // then exit
  int MEASURE_ECW_AND_EXIT;

  // if > 0, then measure and write to file
  // the radius of curvature at each vertex in the model.
  // The mesh subsurface used for sphere fitting is the
  // set of vertices found by a number CURVATURE_NEIGHBORHOOD_SIZE
  // many neighbor-inclusion steps.
  // Minimum valid number of steps is 1,
  // so a value of 0 will do nothing.
  int CURVATURE_NEIGHBORHOOD_SIZE;

  // In case radius of curvature is calculated to be infinite
  // set radius of curvature equal to MAX_RADIUS_OF_CURVATURE.
  double MAX_RADIUS_OF_CURVATURE;

  // if true, then assume all vertices are nice (not inside an object).
  // else explicitely check niceness of each vertex by ray tracing.
  int ASSUME_NICE_VERTICES;

  // if set, then identify and report vertices
  //   whose normal vector lies within threshold of plane of section then exit.
  //   Enter threshold as fraction: 0 <= threshold <= 1.
  //   To convert from degrees (between 0 and 90) to fraction: fraction = degrees/90
  //   If vertex normal is less than or equal to threshold away from plane of section, then report.
  // else do nothing
  double REGION_ORTHOGONALITY_THRESHOLD;

  // if true, then for each vertex count and report number of different neighbor objects
  //   encountered during search for closest point to vertex
  //   Output file is defined by SEP_LOG_FILE.
  // else do nothing
  int COUNT_NEIGHBORS;

  // if true, then for each vertex calculate and report the surace area
  //   vertex as one-third the area of each adjacent face of vertex.
  //   Output file is defined by SEP_LOG_FILE.
  // else do nothing
  int REPORT_VERTEX_AREA;

  // if true, then for each vertex report the index and parent object name.
  //   Output file is defined by SEP_LOG_FILE.
  // else do nothing
  int REPORT_VERTEX_IDENTITY;

  // if true, then prohibit vertices identified as facing tunnel-like
  //   extracellular space from moving. Vertex identity (sheet or tunnel)
  //   is determined by counting the number of neighbor objects. Mutually
  //   exclusive with 'freeze_sheets' option.
  // else, allow all vertex moves.
  int FREEZE_TUNNELS;

  // if true, then prohibit vertices identified as facing sheet-like
  //   extracellular space from moving. Vertex identity (sheet or tunnel)
  //   is determined by counting the number of neighbor objects. Mutually
  //   exclusive with 'freeze_tunnels' option.
  // else, allow all vertex moves.
  int FREEZE_SHEETS;

  // update vertex identity every VERTEX_NEIGHBOR_COUNT_PERIOD groups
  // (that is, of course, if vertex identity is being used)
  int VERTEX_NEIGHBOR_COUNT_PERIOD;

  std::string CLEFT_VERTICES_FILE;

  // if not empty string, interpret RECON_BLOCK_WRAP as name of object
  // wrapping recon to be moved as close as possible to recon.
  // if equal to empty string, then do nothing.
  std::string RECON_BLOCK_WRAP;
};

#endif
