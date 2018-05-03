// Author: Justin Kinney
// Date: Feb 2009

#ifndef CONTOUR_H
#define CONTOUR_H 1

#include <iostream>
#include <string>
#include <vector>

#include "control_points.h"
#include "histogram.h"
#include "parameter.h"
#include "point.h"
#include "sim_anneal.h"

typedef std::vector<int>                           vec_i;
typedef std::vector<int>::iterator          v_i_iterator;
typedef std::vector<int>::const_iterator  c_v_i_iterator;
typedef std::vector<int>::const_reverse_iterator  c_r_v_i_iterator;
typedef std::vector<double>                        vec_d;
typedef std::vector<double>::iterator         d_iterator;
typedef std::vector<double>::const_iterator c_d_iterator;
typedef std::vector<Point>::iterator          p_iterator;
typedef std::vector<Point>::const_iterator  c_p_iterator;
typedef std::list<int>                            list_i; 
typedef std::list<int>::iterator              i_iterator; 

class Parameter;

class Contour
{
private:
  Parameter                                   s;
  Control_Points                             cp;
  int                                   section;
  double             radius_of_curvature_offset;
  std::string                              name;
  std::string                            header;
  std::vector<Point>                 raw_points;
  std::vector<double>            spline_lengths;
  std::vector<double>       radius_of_curvature;
  std::vector<int>               path_parameter;
  std::vector<double>    uniform_path_parameter;
  std::vector<double> final_radius_of_curvature;
  static const double integration_ratio = 1.0/3.0;
  static const double spline_ratio = 1.0/6.0;
public:
  Contour (char * const str, const char * head, int sec);
  Contour &  operator =  (Contour const &);
  Contour                (Contour const &);
  int processContour                (vec_d & my_path_intervals);
  bool isDegenerate                 (void) const;
  void computeDeviationHistogram    (Histogram & h) const ;
  void computeIntervalHistogram     (Histogram & h) const ;
  void updateDeviations             (Histogram & h) const;
  void calcSampleIntervals          (vec_d & my_path_intervals);
  void updateSampleIntervals        (Histogram & h);
  void setSamplesToRaw              (void);
  void removeDuplicates             (void);
  void writeFilesAfter              (void) const;
  void printPtsFiles                (char const * myname) const;
  void printRawPtsFiles             (char const * myname) const;
  void printScripts                 (void) const;
  void printSimpleScripts           (void) const;
  void printRawPoints               (char const * const filename) const;
  void printSplineSamplesHere       (char const * const filename,
                                     const int & num_splines,
                                     double path_parameter_value) const;
  std::string getSerString          (void) const;
  void duplicateControlPoints       (list_i & control_points) const;
private:
  void filterCurvature              (void);
  bool sampleJumpedNeighbor         (const int & path_parameter_index,
                                     const double & parameter_value_change) const;
  void simAnnealing                 (Sim_Anneal & SA,
                                     const int & num_samples);
  void computeHistogram             (Histogram & h,
                                     const c_d_iterator & data_start,
                                     const c_d_iterator & data_end) const;
  void sampleSplines                (void);
  void uniformSampling              (const int & num_splines,
                                     const int & target_num_samples);
  void calculateRadiusOfCurvature   (const int & target_num_samples);
  void printMaxDeviation            (void) const;
  void printSplineSample            (char const * const filename) const;
  void printArcLengths              (char const * const filename) const;
  void printRadiusOfCurvature       (char const * const filename,
                                     const vec_d & radii) const;
  double getEnergy                  (const int & path_parameter_index) const;
  double getTotalEnergy             (void) const;
  double getSplineLength            (const int & curve_segment_index,
                                     const double & path_parameter_value) const;
  double getInterSampleSplineLength (const int & lhs_parameter_index,
                                     const int & rhs_parameter_index) const;
  double calculateProximityEnergy   (const int & lhs_parameter_index,
                                     const int & rhs_parameter_index) const;
  double calculateSplineLength      (const int & curve_segment_index,
                                     const double & path_parameter_value) const;
  double evalRadiusOfCurvature      (const int & path_parameter_index,
                                     bool return_log) const;
  double getFinalRadiusOfCurvature  (const int & sample_index) const;
  Point sampleSpline                (const double & path_parameter_value,
                                     const int & spline_index) const;
public:
  void printIdentity (void)
  {
    std::cout << name << " " << section << std::endl;
  }
  void addRawPoint (Point p)
  {
    raw_points.push_back(p);
  }
  int getNumSamplePoints (void) const
  {
    return s.getNumSamplePoints();
  }
  int getNumRawPoints (void) const
  {
    return raw_points.size();
  }
  bool pointsAreIdentical (Point const & lhs,Point const & rhs) const
  {
    return lhs.getX()==rhs.getX() && lhs.getY()==rhs.getY();
  }
  char const * getName (void) const
  {
    return name.c_str();
  }
  c_d_iterator firstDeviation (void)
  {
    return s.getFirstDeviationConst();
  }
  c_d_iterator onePastLastDeviation (void)
  {
    return s.getOnePastLastDeviationConst();
  }
  c_p_iterator getFirstRawPoint (void)
  {
    return raw_points.begin();
  }
  c_p_iterator getOnePastLastRawPoint (void)
  {
    return raw_points.end();
  }
  int getSection (void) const
  {
    return section;
  }

  /** Create matrix of control points describing all splines in contour.
   * \param[in] control_points Sequence of raw points to use as control_points.
   */

  void initializeControlPointMatrix (list_i & control_points)
  {
    cp.initialize(control_points);
  }

  /** Create matrix of control points describing all splines in contour.
   * \param[out] control_points Sequence of raw points to use as control_points.
   */

  void initializeControlPointVector (std::list<int> & control_points)
  {
    Controls & cs (Controls::instance()); 
    control_points.clear();
    for (int i=0;i<static_cast<int>(raw_points.size());i++)
    {
      control_points.push_back(i);
      if (cs.getReturnInterpolatedRawPoints())
      {
        control_points.push_back(i);
        control_points.push_back(i);
      }
    }
    control_points.push_back(0);
    if (cs.getReturnInterpolatedRawPoints())
    {
      control_points.push_back(0);
      control_points.push_back(0);
    }
    else
    {
      control_points.push_back(1);
      control_points.push_back(2);
    }
  }

  /** Check if any measured deviation of spline samples from
   *  linearly-connected raw contour points exceeds threshold.
   * \return True if deviation threshold is exceeded; false otherwise.
   */

  bool deviationsLarge (void)
  {
    return s.deviationsLarge();
  }

  /** Compute deviation distance between raw points and splines.
   */

  void computeDeviations (void)
  {
    s.computeDeviations(raw_points);
  }

  /** Check if contour has too few raw points to spline but
   * sufficent number to pass threshold criterion.
   * \return True if keep contour; false otherwise.
   */

  bool fewButKeep (void)
  {
    if (static_cast<int>(raw_points.size())<4 &&
       static_cast<int>(raw_points.size())>=
        Controls::instance().getPtPerContourThreshold())
    {
      std::cout << "Contour '" << getName()
                << "' on section " << getSection()
                << " with " << getNumRawPoints()
                << " raw points passed unaltered since insufficient number "
                << "of points (4) for splining but exceeds user threshold ("
                << Controls::instance().getPtPerContourThreshold()
                << ")." << std::endl;
      return true;
    }
    else
    {
      return false;
    }
  }

  /** Check if contour has too few raw points to spline and
   * an insufficent number to pass threshold criterion.
   * \return True if omit contour; false otherwise.
   */

  bool fewOmit (void)
  {
    if (static_cast<int>(raw_points.size())<4 &&
       static_cast<int>(raw_points.size())<
        Controls::instance().getPtPerContourThreshold())
    {
      std::cout << "Contour '" << getName()
                << "' on section " << getSection()
                << " with " << getNumRawPoints()
                << " raw points omitted since insufficient number "
                << "of points for either splining (4) or to exceed threshold ("
                << Controls::instance().getPtPerContourThreshold()
                << ")." << std::endl;
      return true;
    }
    else
    {
      return false;
    }
  }


  /** Initialize radius of curvature log file to empty.
   */

  void clearRadiusCurvatureLogFiles (void)
  {
    char filename[256];
    sprintf(filename,"%s%s_%i_radius_curvature.log",
            Controls::instance().getOutputDir(),getName(),getSection());
    FILE *F = fopen(filename,"w");
    if (!F) { printf("Could not open output file %s\n",filename); return; }
    fclose(F);
  }

  /** Initialize arc lengths log file to empty.
   */

  void clearArcLengthsLogFiles (void)
  {
    char filename[256];
    sprintf(filename,"%s%s_%i_arclengths.log",
            Controls::instance().getOutputDir(),getName(),getSection());
    FILE *F = fopen(filename,"w");
    if (!F) { printf("Could not open output file %s\n",filename); return; }
    fclose(F);
  }

  /** Initialize log file to empty for spline samples at uniform path parameter.
   */

  void clearUniformSampleLogFiles (void)
  {
    char filename[256];
    sprintf(filename,"%s%s_%i_uniform_sample.log",
            Controls::instance().getOutputDir(),getName(),getSection());
    FILE *F = fopen(filename,"w");
    if (!F) { printf("Could not open output file %s\n",filename); return; }
    fclose(F);
  }

  /** Initialize log file to empty for spline samples.
   */

  void clearSampleLogFiles (void)
  {
    char filename[256];
    sprintf(filename,"%s%s_%i_sample.log",
            Controls::instance().getOutputDir(),getName(),getSection());
    FILE *F = fopen(filename,"w");
    if (!F) { printf("Could not open output file %s\n",filename); return; }
    fclose(F);
  }

  /** Initialize path parameter log file to empty.
   */

  void clearSParamLogFiles (void)
  {
    char filename[256];
    sprintf(filename,"%s%s_%i_param_s.log",
            Controls::instance().getOutputDir(),getName(),getSection());
    FILE *F = fopen(filename,"w");
    if (!F) { printf("Could not open output file %s\n",filename); return; }
    fclose(F);
  }

  /** Initialize log file to empty for indeces of control points.
   */

  void clearControlLogFiles (void)
  {
    char filename[256];
    sprintf(filename,"%s%s_%i_control.log",
            Controls::instance().getOutputDir(),getName(),getSection());
    FILE *F = fopen(filename,"w");
    if (!F) { printf("Could not open output file %s\n",filename); return; }
    fclose(F);
  }
private:
  /** Calculate target number of samples for this contour.
   * \param[in] num_splines Number of splines in contour.
   * \return Target number of contour samples.
     Possible values are >2 and 0.
   */

  int calculateNumSamples (const int & num_splines)
  {
    assert(num_splines>0);
    Controls & cs (Controls::instance()); 
    double total_spline_length = 0.0;
    for (int i=0;i<num_splines;i++)
    {
      total_spline_length += getSplineLength(i,1.0);
    }
    assert(total_spline_length>0.0);
    double N_min = total_spline_length/cs.getMaxSampleInterval();
    double N_max = total_spline_length/cs.getMinSampleInterval();
    int a = N_min + (N_max-N_min)*cs.getAdditionalPointsFactor();
    if (a>2) return a;
    return static_cast<int>(N_max);
  }

  /** Compute offset so radius of curvature values are all positive.
   */

  void setOffset (void)
  {
    double min_radius_of_curvature = 1E30;
    for (c_d_iterator i=final_radius_of_curvature.begin();i!=final_radius_of_curvature.end();i++)
    {
      if (*i<min_radius_of_curvature) min_radius_of_curvature = *i;
    }
    radius_of_curvature_offset = fabs(min_radius_of_curvature)+1.0;
  }

  /** Calculate rate of change of path length with respect
   *  to path parameter evaluated at path parameter.
   * \param[in] t Path parameter. 0<t<1.
   * \param[in] spline_index Index of spline.
   * \return Rate of change of path length.
   */

  double evalSprime (const double & t,const int & spline_index) const
  {
    // evaluate function f(t) at integration steps
    // f(t) = Sqrt[(cx + 2 bx t + 3 ax t^2)^2 + (cy + 2 by t + 3 ay t^2)^2]
    double a = getCx(spline_index)+
           2.0*getBx(spline_index)*t+
           3.0*getAx(spline_index)*t*t;
    double b = getCy(spline_index)+
           2.0*getBy(spline_index)*t+
           3.0*getAy(spline_index)*t*t;
    return sqrt(a*a + b*b);
  }

  /** Calculate cumulative length of all splines fit to contour points. 
   * \return Total spline length.
   */

  double calculateTotalPathLength (void)
  {
    double total_length = 0.0;
    for (c_d_iterator i = spline_lengths.begin();i!=spline_lengths.end();i++)
    {
      total_length += *i;
    }
    return total_length;
  }

  /** Calculate and store the path length of each spline in contour.
   * \param[in] num_splines Number of splines in contour.
   */

  void initCurveLengths (const int & num_splines)
  {
    spline_lengths.clear();
    for (int i=0;i<num_splines;i++)
    {
      spline_lengths.push_back(calculateSplineLength(i,1.0));
    }
  }

  /** Get x coordinate of first control point of spline.
   * \param[in] spline_index Index of spline.
   * \return X coordinate of control point.
   */

  double getControlPointX (const int & spline_index,
                           const int & control_point_number) const
  {
    int raw_points_index = cp.getControlPointIndex(4*spline_index+control_point_number);
#ifndef NOASSERT
    assert(raw_points_index>=0);
    assert(raw_points_index<static_cast<int>(raw_points.size()));
#endif
    return raw_points[raw_points_index].getX();
  }

  /** Get y coordinate of first control point of spline.
   * \param[in] spline_index Index of spline.
   * \return Y coordinate of control point.
   */

  double getControlPointY (const int & spline_index,
                           const int & control_point_number) const
  {
    int raw_points_index = cp.getControlPointIndex(4*spline_index+control_point_number);
#ifndef NOASSERT
    assert(raw_points_index>=0);
    assert(raw_points_index<static_cast<int>(raw_points.size()));
#endif
    return raw_points[raw_points_index].getY();
  }

  /** Calculate coefficient of cubic parameter term in spline x location equation.
   * \param[in] spline_index Index of spline.
   * \return Coefficient of cubic parameter term.
   */

  double getAx (const int & spline_index) const
  {
    return spline_ratio*(-1.0*getControlPointX(spline_index,0)
                         +3.0*getControlPointX(spline_index,1)
                         -3.0*getControlPointX(spline_index,2)
                         +1.0*getControlPointX(spline_index,3));
  }

  /** Calculate coefficient of quadratic parameter term in spline x location equation.
   * \param[in] spline_index Index of spline.
   * \return Coefficient of quadratic parameter term.
   */

  double getBx (const int & spline_index) const
  {
    return spline_ratio*( 3.0*getControlPointX(spline_index,0)
                         -6.0*getControlPointX(spline_index,1)
                         +3.0*getControlPointX(spline_index,2)
                         +0.0*getControlPointX(spline_index,3));
  }

  /** Calculate coefficient of linear parameter term in spline x location equation.
   * \param[in] spline_index Index of spline.
   * \return Coefficient of linear parameter term.
   */

  double getCx (const int & spline_index) const
  {
    return spline_ratio*(-3.0*getControlPointX(spline_index,0)
                         +0.0*getControlPointX(spline_index,1)
                         +3.0*getControlPointX(spline_index,2)
                         +0.0*getControlPointX(spline_index,3));
  }

  /** Calculate constant term in spline x location equation.
   * \param[in] spline_index Index of spline.
   * \return Constant term.
   */

  double getDx (const int & spline_index) const
  {
    return spline_ratio*( 1.0*getControlPointX(spline_index,0)
                         +4.0*getControlPointX(spline_index,1)
                         +1.0*getControlPointX(spline_index,2)
                         +0.0*getControlPointX(spline_index,3));
  }

  /** Calculate coefficient of cubic parameter term in spline y location equation.
   * \param[in] spline_index Index of spline.
   * \return Coefficient of cubic parameter term.
   */

  double getAy (const int & spline_index) const
  {
    return spline_ratio*(-1.0*getControlPointY(spline_index,0)
                         +3.0*getControlPointY(spline_index,1)
                         -3.0*getControlPointY(spline_index,2)
                         +1.0*getControlPointY(spline_index,3));
  }

  /** Calculate coefficient of quadratic parameter term in spline y location equation.
   * \param[in] spline_index Index of spline.
   * \return Coefficient of quadratic parameter term.
   */

  double getBy (const int & spline_index) const
  {
    return spline_ratio*( 3.0*getControlPointY(spline_index,0)
                         -6.0*getControlPointY(spline_index,1)
                         +3.0*getControlPointY(spline_index,2)
                         +0.0*getControlPointY(spline_index,3));
  }

  /** Calculate coefficient of linear parameter term in spline y location equation.
   * \param[in] spline_index Index of spline.
   * \return Coefficient of linear parameter term.
   */

  double getCy (const int & spline_index) const
  {
    return spline_ratio*(-3.0*getControlPointY(spline_index,0)
                         +0.0*getControlPointY(spline_index,1)
                         +3.0*getControlPointY(spline_index,2)
                         +0.0*getControlPointY(spline_index,3));
  }

  /** Calculate constant term in spline y location equation.
   * \param[in] spline_index Index of spline.
   * \return Constant term.
   */

  double getDy (const int & spline_index) const
  {
    return spline_ratio*( 1.0*getControlPointY(spline_index,0)
                         +4.0*getControlPointY(spline_index,1)
                         +1.0*getControlPointY(spline_index,2)
                         +0.0*getControlPointY(spline_index,3));
  }

};

#endif
