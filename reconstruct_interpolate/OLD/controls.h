// Author: Justin Kinney
// Date: Feb 2009

#ifndef CONTROLS_H
#define CONTROLS_H 1

#include <string.h>

#include <getopt.h>
#include <iostream>
#include <string>
#include <vector>


class Controls
{
public:
  static Controls & instance     (void);
  void   parseCommandLine        (int,char**);
  std::string getUsageMessage    (void);
  int    getNumReserve           () const throw() { return OBJECT_RESERVE_SIZE; }
  int    getMaxContoursPerObject () const throw() { return MAX_CONTOURS_PER_OBJECT; }
  int    getCappingFlag          () const throw() { return CAPPING_FLAG; }
  int    getMaxSection           () const throw() { return MAX_SECTION; }
  int    getMinSection           () const throw() { return MIN_SECTION; }
  int    getSplineSamplesPerPoint                  () const throw() { return SPLINE_SAMPLES_PER_POINT; }
  int    getPrintRawPoints       () const throw() { return PRINT_RAW_POINTS; }
  int    getMaxDeviationAdjustments () const throw() { return MAX_DEVIATION_ADJUSTMENTS; }
  bool   getDiag                 () const throw() { return diag; }
  double getMaxRad               () const throw() { return max_rad; }
  double getT                    () const throw() { return T; }
  double getAmax                 () const throw() { return amax; }
  double getDmax                 () const throw() { return dmax; }
  double getDmin                 () const throw() { return dmin; }
  double getScale                () const throw() { return SCALE; }
  double getDeviationThreshold   () const throw() { return DEVIATION_THRESHOLD; }
  double getLinearThreshold      () const throw() { return LINEAR_THRESHOLD; }
  double getSectionThickness     () const throw() { return SECTION_THICKNESS; }
  char const * getOutputScript   () const throw() { return OUTPUT_SCRIPT.c_str(); }
  char const * getInputDir       () const throw() { return INPUT_DIR.c_str(); }
  char const * getOutputDir      () const throw() { return OUTPUT_DIR.c_str(); }
  char const * getPrefix         () const throw() { return PREFIX.c_str(); }
private:
  static Controls * only_one;
  Controls                (void);
  Controls                (Controls const &);
  Controls & operator =   (Controls const &);
  // thickness of each section (assumed uniform across sections)
  // units [same as input data]
  double SECTION_THICKNESS;
  // scaling factor of output contour position before input to contour_tiler 
  // example, SCALE==1000, convert micron Reconstruct units to nanometer scale
  // before meshing with contour_tiler
  // example, SCALE==1, pass micron Reconstruct units
  // directly to contour_tiler fore meshing
  // units [-]
  double SCALE;
  // Sampled spline points are not allowed to deviated from
  // input contour points by more than DEVIATION_THRESHOLD.
  double DEVIATION_THRESHOLD;
  // All pairs of adjacent raw points on a contour
  // separated by a distance less than LINEAR_THRESHOLD
  // will be interpolated by additional points to satisfy
  // the threshold.
  // units [same as input data]
  double LINEAR_THRESHOLD;
  // if nonzero, then contour_tiler will be instructed to cap the meshes
  // at min and max z sections.
  // if zero, then do nothing. meshes will be left open at min and max z.
  int CAPPING_FLAG;
  int MIN_SECTION;
  int MAX_SECTION;
  std::string INPUT_DIR;
  std::string OUTPUT_DIR;
  std::string PREFIX;
  std::vector<std::string> IGNORED_CONTOURS;
  std::string OUTPUT_SCRIPT;
  int SPLINE_SAMPLES_PER_POINT;		// num is the # samples of the splines between contour points
  // sequential triplets of sampled points are used to
  // compute the radius of curvature , e.g. num/2 between contour points
  bool diag;		// set true to print diagnostic files
  double max_rad;	// calculated radius of curvature of spline will saturate at this value
  double dmin;		// minimum spline sampling distance
  double dmax;		// maximum spline sampling distance
  double T;		// sample period (= time to traverse dmax)
  double amax;		// max radial acceleration
  int OBJECT_RESERVE_SIZE;
  int MAX_CONTOURS_PER_OBJECT;
  int MAX_DEVIATION_ADJUSTMENTS;
  // if nonzero, then print raw_points (possiibly linearly interpolated) to pts files
  // if zero, then do nothing
  int PRINT_RAW_POINTS;
public:

  /** Determine if contour is to be ignored.
   * \param[in] name Name of contour.
   * \return True if contour is to be ignored; false otherwise.
   */

  bool contourIsIgnored (char const * const name) const
  {
    for (std::vector<std::string>::const_iterator i=IGNORED_CONTOURS.begin();i!=IGNORED_CONTOURS.end();i++)
    {
      if (!strcmp((*i).c_str(),name)) return true;
    }
    return false;
  }

  /** Create string from floating point number.
   * \param[in] num Number of interest.
   * \return String containing input number.
   */

  std::string d2str (double const & num)
  {
    char file[1024];
    sprintf(file,"%.15g",num);
    return std::string(file);
  }

  /** Create string from integer number.
   * \param[in] num Number of interest.
   * \return String containing input number.
   */

  std::string i2str (int const & num)
  {
    char file[1024];
    sprintf(file,"%i",num);
    return std::string(file);
  }

  /** Ensure pathname ends with exactly one '/'.
   * \param[in] ptr Arbitrary pathname.
   * \return Input pathname with '/' added if absent.
   */

  std::string processDir (char * ptr)
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
};

#endif
