// Author: Justin Kinney
// Date: Feb 2009

#ifndef CONTROLS_H
#define CONTROLS_H 1

#include <string>

class Controls
{
private:
  std::string INPATH;
  int FOLDER;			// true if folder of input mesh files, false if single input mesh file
  int ATTR;				// true if '-a', false otherwise
  // true=eval attributes only
  // false=allow a full analysis, i.e. evaluate characteristics
  int PRINT;				// true if '-p', false otherwise
  // true=print detailed information about offending mesh elements
  // false=print summary information only
  int DREAMM;			// output style for offending mesh elements (-p option)
  // nonzero=dreamm custom points format, i.e. x y z 1 0 0 1
  // zero=detailed face,vertex,edge information
  int INTER;			// true if '-i', false otherwise
  // true=collect all objects,
  //	detect intersections between objects,
  //	in contrast to the intra-object search that is the default
  // false=do not perform inter-object face intersection check
  int PRINT_SET_VOLUME_ONLY;				// true if '-v', false otherwise
  // true=print set volume and nothing else
  // false=do nothing special

public:
  static Controls & instance     (void);
  // command line info
  bool signal[5];			// true if threshold value found on command line
  double thresholds[5];	// user-defined thresholds found on command line
  // [0]=aspect ratio threshold set
  // [1]=min edge angle threshold set
  // [2]=max edge angle threshold set
  // [3]=min edge length threshold set
  // [4]=max edge length threshold set

  void parse(int,char**,std::string);
  std::string getUsageMessage (void);

  int    get_print_set_volume_only            () const throw() { return PRINT_SET_VOLUME_ONLY; }
  int    get_detect_interobject_intersections () const throw() { return INTER; }
  int    get_use_dreamm_output_format         () const throw() { return DREAMM; }
  int    get_print_detailed_info              () const throw() { return PRINT; }
  int    get_compute_attributes_only          () const throw() { return ATTR; }
  int    get_folder_passed                    () const throw() { return FOLDER; }
  int    get_filename_size                    () const throw() { return FILENAME_SIZE; }
  int    get_detect_polygon_edge_intersection () const throw() { return DETECT_POLYGON_EDGE_INTERSECTION;; }
  int    get_num_adjacent_boxes            () const throw() { return NUM_ADJACENT_BOXES; }
  int    get_num_bins                      () const throw() { return NUM_BINS; }
  double get_ray_epsilon                   () const throw() { return RAY_EPSILON; }
  double get_double_epsilon                () const throw() { return DOUBLE_EPSILON; }
  double get_space_length                  () const throw() { return SPACE_LENGTH; }
  double get_faces_per_box                 () const throw() { return FACES_PER_BOX; }
  double get_search_radius                 () const throw() { return SEARCH_RADIUS; }
  double get_neighborhood_radius           () const throw() { return NEIGHBORHOOD_RADIUS; }
  double get_closest_point_cosine          () const throw() { return CLOSEST_POINT_COSINE; }
  double get_pi                            () const throw() { return PI; }
  std::string get_inpath                   () const throw() { return INPATH; }

private:
  static Controls * only_one;
  Controls                (void);
  Controls                (Controls const &);
  Controls & operator =   (Controls const &);

  // default filename size
  int FILENAME_SIZE;

  // set to true to increase accuracy of line/polygon intersections at expense of speed
  // set to false to increase speed at expense of accuracy during line/polygon intersections
  bool DETECT_POLYGON_EDGE_INTERSECTION;

  // small value for tweaking nice ray to avoid
  // intersecting polygon edge or
  // being parallel with polygon
  double RAY_EPSILON; //nm

  // for use with "is float close to zero?" 
  // in conditional statement
  double DOUBLE_EPSILON;

  // subdivide space into cubes with sides of SPACE_LENGTH in size
  //#define SPACE_LENGTH	40 // nm
  double SPACE_LENGTH; // nm

  // target mean number faces per space subpartition
  //#define FACES_PER_BOX	0.3
  double FACES_PER_BOX;

  // When searching for closest surface point to a given vertex,
  // the subdivided space box of given vertex will be searched,
  // as well as, all boxes within NUK_ADJACENT_BOXES away from 
  // given vertex box. If NUM_ADJACENT_BOXES is zero, then only
  // the given vertex box will be searched. If NUM_ADJACENT_BOXES
  // is one, then all adjacent boxes to given vertex box will also be searched.
  //
  // choose such that NUM_ADJACENT_BOXES >= SEARCH_RADIUS/SPACE_LENGTH
  int NUM_ADJACENT_BOXES;

  // If a candidate closest point is more than SEARCH_RADIUS 
  // distance away from given vertex, then the candidate closest 
  // point is not chosen as closest point to given vertex.
  // PICK SEARCH_RADIUS < NEIGHBORHOOD_RADIUS else closest point
  // will be be just outside neighborhood on same object if no
  // other point is found.
  //#define SEARCH_RADIUS  70.0 // nm
  double SEARCH_RADIUS; // nm

  // For a given vertex, faces containing vertices within
  // the neighborhood of the given vertex are eliminated from
  // the search pool for closest point for the given vertex.
  // The neighborhood is the collection of all vertices
  // on the same object as the given vertex that lie within
  // N edges of the given vertex, where N is a function of
  // NEIGHBORHOOD_RADIUS and mean edge length for that object.
  //#define NEIGHBORHOOD_RADIUS		140 // nm
  double NEIGHBORHOOD_RADIUS; // nm

  double CLOSEST_POINT_COSINE; //  equivalent to a 70 degree sweep around normal

  // pi
  double PI;

  // number of bins in histograms
  int NUM_BINS;
};

#endif
