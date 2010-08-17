// Author: Justin Kinney
// Date: Feb 2009

#ifndef CONTROLS_H
#define CONTROLS_H 1

#include <string>

#include "meshheal.h"

class Controls
{
public:
  static Controls & instance  (void);
  void        parse           (int,char**,std::string);
  std::string getUsageMessage (void);

  int         get_max_filename_size   () const throw() { return MAX_FILENAME_SIZE; }
  double      get_double_epsilon      () const throw() { return DOUBLE_EPSILON; }
  double      get_distance_threshold  () const throw() { return DISTANCE_THRESHOLD; }
  vector3     get_random_vector       () const throw() { return RANDOM_VECTOR; }
  std::string get_inpath              () const throw() { return INPATH; }
  std::string d2str                   (double const & i);
  std::string i2str                   (int const & i);

private:
  static Controls * only_one;
  Controls                (void);
  Controls                (Controls const &);
  Controls & operator =   (Controls const &);

private:
  std::string INPATH;

  // default filename size
  int MAX_FILENAME_SIZE;

  // maximum allowed distance between vertices slated for merger
  // Free vertices separated by more than DISTANCE_THRESHOLD
  // will not be merged.
  double DISTANCE_THRESHOLD;

  // for use with "is float close to zero?" 
  // in conditional statement
  double DOUBLE_EPSILON;

  // random vector
  vector3 RANDOM_VECTOR;
};

#endif
