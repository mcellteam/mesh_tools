// Author: Justin Kinney
// Date: Feb 2009

#ifndef CONTROLS_H
#define CONTROLS_H 1

#include <cassert>
#include <getopt.h>
#include <stdio.h>
#include <string.h>
#include <string>
#include <vector>


class Controls
{
public:
  static Controls & instance             (void);
  void   parseCommandLine                (int,char**);
  std::string getUsageMessage            (void);
  int    getNumReserve                   () const throw() { return OBJECT_RESERVE_SIZE; }
  int    getCappingFlag                  () const throw() { return CAPPING_FLAG; }
  int    getMinSection                   () const throw() { return MIN_SECTION; }
  int    getMaxSection                   () const throw() { return MAX_SECTION; }
  int    getPrintDetailedInfo            () const throw() { return PRINT_DETAILED_INFO; }
  int    getReturnRawContourPoints       () const throw() { return RETURN_RAW_CONTOUR_POINTS; }
  int    getReturnInterpolatedRawPoints  () const throw() { return RETURN_INTERPOLATED_RAW_POINTS; }
  int    getPtPerContourThreshold        () const throw() { return MIN_PT_PER_CONTOUR_THRESHOLD; }
  int    getNumMovesPerIteration         () const throw() { return SIM_ANNEAL_MOVES_PER_ITERATION; }
  int    getMaxNumMoveAttempts           () const throw() { return SIM_ANNEAL_MAX_NUM_MOVE_ATTEMPTS; }
  int    getIntegrationStep              () const throw() { return SIM_ANNEAL_INTEGRATION_STEP; }
  int    getMaxDeviationAdjustments      () const throw() { return MAX_DEVIATION_ADJUSTMENTS; }
  int    getCurvatureEnergyExponent      () const throw() { return SIM_ANNEAL_CURVATURE_ENERGY_EXPONENT; }
  int    getProximityEnergyExponent      () const throw() { return SIM_ANNEAL_PROXIMITY_ENERGY_EXPONENT; }
  double getAdditionalPointsFactor       () const throw() { return ADDITIONAL_POINTS_FACTOR; }
  double getIntegrationStepFactor        () const throw() { return SIM_ANNEAL_INTEGRATION_STEP_FACTOR; }
  double getHighTemp                     () const throw() { return SIM_ANNEAL_HIGH_TEMP; }
  double getTempScale                    () const throw() { return SIM_ANNEAL_TEMP_SCALE; }
  double getBoltzman                     () const throw() { return SIM_ANNEAL_BOLTZMAN; }
  double getCurvatureEnergyGain          () const throw() { return SIM_ANNEAL_CURVATURE_ENERGY_GAIN; }
  double getProximityEnergyGain          () const throw() { return SIM_ANNEAL_PROXIMITY_ENERGY_GAIN; }
  double getMeanAmplitudeNoise           () const throw() { return SIM_ANNEAL_MEAN_AMPLITUDE_NOISE; }
  double getMaxRadiusOfCurvature         () const throw() { return MAXIMUM_RADIUS_OF_CURVATURE; }
  double getSectionThickness             () const throw() { return SECTION_THICKNESS; }
  double getOutputScaleFactor            () const throw() { return OUTPUT_SCALE_FACTOR; }
  double getDeviationThreshold           () const throw() { return DEVIATION_THRESHOLD; }
  double getEpsilon                      () const throw() { return EPSILON; }
  double getMaxSampleInterval            () const throw() { return MAX_SAMPLE_INTERVAL; }
  double getMinSampleInterval            () const throw() { return MIN_SAMPLE_INTERVAL; }
  char const * getOutputSerPrefix        () const throw() { return OUTPUT_SER_PREFIX.c_str(); }
  char const * getInputDir               () const throw() { return INPUT_DIR.c_str(); }
  char const * getOutputDir              () const throw() { return OUTPUT_DIR.c_str(); }
  char const * getPrefix                 () const throw() { return PREFIX.c_str(); }
  char const * getOutputScript           () const throw() { return OUTPUT_SCRIPT.c_str(); }
  char const * getMultiPartSuffix        () const throw() { return MULTI_PART_SUFFIX.c_str(); }

private:
  static Controls * only_one;
  Controls                (void);
  Controls                (Controls const &);
  Controls & operator =   (Controls const &);

  // accelerates program execution by pre-allocating memory
  // set OBJECT_RESERVE_SIZE equal to number of objects to be reconstructed
  // e.g. if reconstructing two dendrites and three axons, then set to five
  // OBJECT_RESERVE_SIZE > 0
  int OBJECT_RESERVE_SIZE;

  // if nonzero, then contour_tiler will be instructed to cap the meshes
  // at min and max z sections.
  // if zero, then do nothing. meshes will be left open at min and max z.
  int CAPPING_FLAG;

  // the range of sections to be analyzed
  int MIN_SECTION;
  int MAX_SECTION;

  // if nonzero, then print raw_points, s parameter values, etc to log files
  // if zero, then do nothing
  int PRINT_DETAILED_INFO;

  // if nonzero, then return input contour points unadulterated
  // if zero then fit splines to contour points and resample
  int RETURN_RAW_CONTOUR_POINTS;

  // if nonzero, then return input contour points
  // linearly interpolated to satisfy minimum and
  // maximum sample interval constraints.
  // if zero then fit splines to contour points and resample
  int RETURN_INTERPOLATED_RAW_POINTS;

  // number of points per contour
  // below which contours will be excluded and not processed
  // Contours with fewer than MIN_PT_PER_CONTOUR_THRESHOLD
  // output sample points will be omitted from output data
  // files and statistics.
  // MIN_PT_PER_CONTOUR_THRESHOLD=0 implies that all contours are included.
  // MIN_PT_PER_CONTOUR_THRESHOLD=3 screens out all contours
  // with fewer than 3 ouput sample points there by insuring that output contours
  // will not generate errors by contour_tiler.
  int MIN_PT_PER_CONTOUR_THRESHOLD;

  // Aside: Contours are fit by cubic splines. The cubic splines are sampled.
  //        The location on the splines of the samples is determined by
  //        simulated annealing.

  // Desired number of successful moves per sample per annealing iteration.
  // An annealing iteration is completed when the number of successful sample moves
  // equals SIM_ANNEAL_MOVES_PER_ITERATION * number of samples in contour.
  // The system is not frozen and annealing continues with another iteration
  // at lower temperature.
  int SIM_ANNEAL_MOVES_PER_ITERATION;

  // Max number of attempted moves per sample per annealing iteration.
  // An annealing iteration is aborted when the number of attempted sample moves
  // equals SIM_ANNEAL_MAX_NUM_MOVE_ATTEMPTS * number of samples in contour.
  // Three aborted iterations in a row indicates that the contour samples
  // are frozen and simulated annealing is complete.
  // SIM_ANNEAL_MAX_NUM_MOVE_ATTEMPTS > SIM_ANNEAL_MOVES_PER_ITERATION
  int SIM_ANNEAL_MAX_NUM_MOVE_ATTEMPTS;

  // Cubic spline length is determined by numerical integration
  // where the function to be integrated is broken into 
  // SIM_ANNEAL_INTEGRATION_STEP-1 pieces.
  // SIM_ANNEAL_INTEGRATION_STEP must be odd and >= 3
  // SIM_ANNEAL_INTEGRATION_STEP_FACTOR is inverse of SIM_ANNEAL_INTEGRATION_STEP-1
  int SIM_ANNEAL_INTEGRATION_STEP;
  double SIM_ANNEAL_INTEGRATION_STEP_FACTOR;

  // If the deviation between sampled cubic spline points and linearly
  // interpolated input contour points is larger than DEVIATION_THRESHOLD
  // then control points are duplicated no more than MAX_DEVIATION_ADJUSTMENTS
  // times to reduce the deviation.
  int MAX_DEVIATION_ADJUSTMENTS;

  // The number N of contour samples will be calculated as follows.
  // Let L = contour length
  //     N_min = L/MAX_SAMPLE_INTERVAL
  //     N_max = L/MIN_SAMPLE_INTERVAL
  // N = N_min + (N_max-N_min)*ADDITIONAL_POINTS_FACTOR
  // where 0 <= ADDITIONAL_POINTS_FACTOR <= 1
  double ADDITIONAL_POINTS_FACTOR;

  // Simulated annealing progresses from an initial high temp
  // (SIM_ANNEAL_HIGH_TEMP) to a very low temperature.
  // SIM_ANNEAL_HIGH_TEMP>0
  double SIM_ANNEAL_HIGH_TEMP; 

  // After each iteration of simulated annealing the temperature
  // of the system is lowered by multiplying the old temperature
  // by SIM_ANNEAL_TEMP_SCALE. 0<SIM_ANNEAL_TEMP_SCALE<1
  double SIM_ANNEAL_TEMP_SCALE; 

  // Inverse of boltzman constant.
  double SIM_ANNEAL_BOLTZMAN; 

  // A sample point contributes to energy in a manner proportional
  // to SIM_ANNEAL_CURVATURE_ENERGY_GAIN and inversely proportional
  // to curvature.
  double SIM_ANNEAL_CURVATURE_ENERGY_GAIN;

  // A sample point contributes to energy in a manner inversely proportional
  // to curvature to the SIM_ANNEAL_CURVATURE_ENERGY_EXPONENT power.
  int SIM_ANNEAL_CURVATURE_ENERGY_EXPONENT;

  // A sample point contributes to energy in a manner proportional
  // to SIM_ANNEAL_PROXIMITY_ENERGY_GAIN and distance between samples.
  double SIM_ANNEAL_PROXIMITY_ENERGY_GAIN;

  // A sample point contributes to energy in a manner proportional to distance
  // between samples to the SIM_ANNEAL_PROXIMITY_ENERGY_EXPONENT power.
  int SIM_ANNEAL_PROXIMITY_ENERGY_EXPONENT;

  // Randomly picked samples are shifted a random amount during simulated
  // annealing. The distribution of shifts is gaussian with a mean
  // of SIM_ANNEAL_MEAN_AMPLITUDE_NOISE. Position along the cubic splines
  // is described a single parameter to which the stochastic perturbation
  // is applied.
  double SIM_ANNEAL_MEAN_AMPLITUDE_NOISE;

  // If radius of curvature of spline is determined to be infinite,
  // i.e. a straight line,  then radius of curvature is set
  // to MAXIMUM_RADIUS_OF_CURVATURE.
  double MAXIMUM_RADIUS_OF_CURVATURE;

  // thickness of each section (assumed uniform across sections)
  // units [same as input data]
  double SECTION_THICKNESS;

  // Scaling factor of output contour positions. Output data
  // will be scaled by OUTPUT_SCALE_FACTOR before being written to file.
  double OUTPUT_SCALE_FACTOR;

  // Sampled spline points are not allowed to deviated from
  // input contour points by more than DEVIATION_THRESHOLD.
  double DEVIATION_THRESHOLD;

  // Value below which a floating point number is considered to be zero.
  double EPSILON;

  // The cubic spline path distance between samples
  // will not exceed MAX_SAMPLE_INTERVAL.
  double MAX_SAMPLE_INTERVAL;

  // The cubic spline path distance between samples
  // will not be less tha MIN_SAMPLE_INTERVAL.
  double MIN_SAMPLE_INTERVAL;

  // Directory containing input contours.
  std::string INPUT_DIR;

  // Directory to which output data will be written.
  std::string OUTPUT_DIR;

  // Whether print to SER or raw points;
  //  0 prints raw points
  //  1 prints to SER files
  int OUTPUT_SER;

  // File name of SER file to output into
  // 'OUTPUT_DIR/OUTPUT_SER_PREFIX.ser'
  // 'OUTPUT_DIR/OUTPUT_SER_PREFIX.MIN_SECTION' to
  // 'OUTPUT_DIR/OUTPUT_SER_PREFIX.MAX_SECTION'
  std::string OUTPUT_SER_PREFIX;

  // The input contours will be read from
  // 'INPUT_DIR/PREFIX.MIN_SECTION' to
  // 'INPUT_DIR/PREFIX.MAX_SECTION'
  std::string PREFIX;

  // Collection of contour names that are to be excluded and not processed.
  std::vector<std::string> EXCLUDED_CONTOURS;

  std::vector<std::string> INCLUDED_CONTOURS;

  // Name of bash script written to output directory
  // which contains commands to tile output contours and convert
  // resulting meshes to Hughes Hoppe mesh format.
  std::string OUTPUT_SCRIPT;

  // Objects with intermediate empty sections (thereby
  // requiring multiple tiling steps) will have name
  // and associated output files appended with
  // MULTI_PART_SUFFIX plus incremented index.
  std::string MULTI_PART_SUFFIX;

public:

  /** Validate the values for control parameters.
   */

  void validate (void)
  {
    // number of integration steps must be even and greater than 1.
    assert(SIM_ANNEAL_INTEGRATION_STEP>1);
    assert(!(SIM_ANNEAL_INTEGRATION_STEP%2));
    // mutually exclusive modes
    assert(!(RETURN_RAW_CONTOUR_POINTS && RETURN_INTERPOLATED_RAW_POINTS));
    // factor between 0 and 1
    assert(ADDITIONAL_POINTS_FACTOR>=0.0);
    assert(ADDITIONAL_POINTS_FACTOR<=1.0);
    // max sample interval greater than or equal to min
    assert(MAX_SAMPLE_INTERVAL>=0.0);
    assert(MIN_SAMPLE_INTERVAL>=0.0);
    assert(MAX_SAMPLE_INTERVAL>=MIN_SAMPLE_INTERVAL);
    // section thickness must be positive and nonzero
    assert(SECTION_THICKNESS>0.0);
    // sections must be positive
    assert(MIN_SECTION>=0);
    assert(MAX_SECTION>=0);
    // simulated annealing must be positive and nonzero
    assert(SIM_ANNEAL_HIGH_TEMP>0.0);
    // simulated annealing temp scaling must be between 0 and 1
    assert(SIM_ANNEAL_TEMP_SCALE>0.0);
    assert(SIM_ANNEAL_TEMP_SCALE<1.0);
    // gains must be positive
    assert(SIM_ANNEAL_CURVATURE_ENERGY_GAIN>=0.0);
    assert(SIM_ANNEAL_PROXIMITY_ENERGY_GAIN>=0.0);
    // exponents must be positive and nonzero
    assert(SIM_ANNEAL_CURVATURE_ENERGY_EXPONENT>0);
    assert(SIM_ANNEAL_PROXIMITY_ENERGY_EXPONENT>0);
    // mean of path parameter displacement distribution must be positive and nonzero
    assert(SIM_ANNEAL_MEAN_AMPLITUDE_NOISE>0.0);
    // maximum radius of curvature must be positive and nonzero
    assert(MAXIMUM_RADIUS_OF_CURVATURE>0.0);
    // output scaling factor must be nonzero
    assert(OUTPUT_SCALE_FACTOR!=0.0);
    // threshold for sample point deviation to original contour must be positive and nonzero
    //assert(DEVIATION_THRESHOLD>0.0);
    // epsilon must be positive and nonzero
    assert(EPSILON>0.0);
  }

  /** Determine if contour is to be excluded.
   * \param[in] name Name of contour.
   * \return True if contour is to be excluded; false otherwise.
   */

  bool contourIsExcluded (char const * const name) const
  {
    for (std::vector<std::string>::const_iterator i=EXCLUDED_CONTOURS.begin();i!=EXCLUDED_CONTOURS.end();i++)
    {
      if (!strcmp((*i).c_str(),name)) return true;
    }
    
    if (INCLUDED_CONTOURS.size() > 0)
    {
      for (std::vector<std::string>::const_iterator i=INCLUDED_CONTOURS.begin();i!=INCLUDED_CONTOURS.end();i++)
      {
        if (!strcmp((*i).c_str(),name)) return false;
      }
      return true;
    }
    
    return false;
  }

  /** Create string from floating point number.
   * \param[in] num Number of interest.
   * \return String containing input number.
   */

  std::string d2str (double const & num) const
  {
    char file[1024];
    sprintf(file,"%.15g",num);
    return std::string(file);
  }

  /** Create string from integer number.
   * \param[in] num Number of interest.
   * \return String containing input number.
   */

  std::string i2str (int const & num) const
  {
    char file[1024];
    sprintf(file,"%i",num);
    return std::string(file);
  }

  /** Ensure pathname ends with exactly one '/'.
   * \param[in] ptr Arbitrary pathname.
   * \return Input pathname with '/' added if absent.
   */

  std::string processDir (char * ptr) const
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
