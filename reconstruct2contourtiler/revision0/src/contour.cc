#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <stdlib.h>

#include "contour.h"

#include "controls.h"
#include "control_points.h"
#include "histogram.h"
#include "point.h"

extern const double spline_ratio;

Contour::Contour (Contour const & rhs)
:s(rhs.s),cp(rhs.cp),section(rhs.section),
  radius_of_curvature_offset(rhs.radius_of_curvature_offset),
  name(rhs.name),raw_points(rhs.raw_points),
  spline_lengths(rhs.spline_lengths),
  radius_of_curvature(rhs.radius_of_curvature),
  path_parameter(rhs.path_parameter),
  uniform_path_parameter(rhs.uniform_path_parameter),
  final_radius_of_curvature(rhs.radius_of_curvature)
{
}

Contour & Contour::operator = (const Contour& rhs)
{
  std::cout << "Assignment operator prohibited on instances of contour class.\n";
  std::cout << "Object " << rhs.name << std::endl;
  exit(1);
}

Contour::Contour (char * const str, int sec)
:s(),cp(),section(sec),radius_of_curvature_offset(0),
  name(),raw_points(),spline_lengths(),radius_of_curvature(),
  path_parameter(),uniform_path_parameter(),final_radius_of_curvature()
{
  char * ptr = str;
  while(ptr!=NULL){
    ptr=strpbrk(str,"#()");
    if(ptr!=NULL){
      *ptr='_';
    }
  }
  name = str;
  assert(name!="");
}

/** Determine if contour has too few samples to be a valid contour shape.
 *  Specifically, a contour with two or less points has zero circumscribing
 *  area and is an impossibly physiological shape.
 * \return True if contour is degenerate; false otherwise.
 */

bool Contour::isDegenerate (void) const 
{
  if (s.getNumSamplePoints()<3) return true;
  else return false;
}

/** Use the raw contour points as samples.
 */

void Contour::setSamplesToRaw (void)
{
  // for each raw point
  for (c_p_iterator i=raw_points.begin();i!=raw_points.end();i++)
  {
    s.addSplineSample(*i);
    s.addDeviation(0.0);
  }
  if (Controls::instance().getPrintDetailedInfo())
  {
    char filename[256];
    sprintf(filename,"%s%s_%i_sample.log",
            Controls::instance().getOutputDir(),getName(),getSection());
    printSplineSample(filename);
  }
}

/** Retrieve filtered radius of curvature that in some way
  * corresponds to the spline sample of interest.
 * \param[in] sample_index Index of spline sample.
 * \return Radius of curvature of cubic spline.
 */

double Contour::getFinalRadiusOfCurvature (const int & sample_index) const
{
  double t = s.getParameterValue(sample_index);
  t += s.getSplineIndex(sample_index);
  int my_index = floor(t*10.0);
#ifndef NOASSERT
  assert(my_index>=0);
  assert(my_index<static_cast<int>(path_parameter.size()));
#endif
  int ss = path_parameter[my_index];
  // THIS IGNORES DIFFERENCES BETWEEN t
  // AND UNIFORM SAMPLING PATH PARAMETER VALUES.
#ifndef NOASSERT
  assert(ss>=0);
  assert(ss<static_cast<int>(final_radius_of_curvature.size()));
#endif
  return final_radius_of_curvature[ss]+radius_of_curvature_offset;
}

/** Calculate radius of curvature of spline evaluated at path parameter.
 * \param[in] path_parameter_index Index of path parameter.
 * \param[in] add_offset If true offset radius of curvature so all values are positive;
                         else just return radius of curvature.
 * \return Radius of curvature of cubic spline.
 */

double Contour::evalRadiusOfCurvature (const int & path_parameter_index,
                                       bool add_offset) const
{
  double t = s.getParameterValue(path_parameter_index);
  int spline_index = s.getSplineIndex(path_parameter_index);
  double cub_coefX  = getAx(spline_index)*t;
  double cub_coefY  = getAy(spline_index)*t;
  double quad_coefX = 2.0*getBx(spline_index);
  double quad_coefY = 2.0*getBy(spline_index);
  double px = quad_coefX+6.0*cub_coefX;
  double py = quad_coefY+6.0*cub_coefY;
  double qx = getCx(spline_index)+quad_coefX*t+3.0*cub_coefX*t;
  double qy = getCy(spline_index)+quad_coefY*t+3.0*cub_coefY*t;
  double num = qx*qx + qy*qy;
  double den = -(px*qx+py*qy)*(px*qx+py*qy)+(px*px+py*py)*num;
  double r = Controls::instance().getMaxRadiusOfCurvature();
  if (den>0.0)
  {
    r = sqrt((num*num*num)/den);
  }
  if (r>Controls::instance().getMaxRadiusOfCurvature())
  {
    r = Controls::instance().getMaxRadiusOfCurvature();
  }
  double m = log(r);
  // LOG BASE 2 OR 4 OR 8 INSTEAD OF 10
  //int n = 0;
  //frexp(r,&n);
  //int m = n >> 2;
  if (add_offset)
  {
    m+=radius_of_curvature_offset;
  }
  return m;
}

/** Compare sequential contour points and remove duplicate points. 
 */

void Contour::removeDuplicates (void)
{
  if (static_cast<int>(raw_points.size())==1) return;
  p_iterator i=raw_points.begin();
  p_iterator j=i;
  j++;
  while (j!=raw_points.end())
  {
    if (pointsAreIdentical(*i,*j))
    {
      // NOTE this is risky for vector
      // but wanted indexing not available in list
      raw_points.erase(i);
      j=i;
      j++;
    }
    else
    {
      i++;
      j++;
    }
  }
  // compare first and last link in list
  i = raw_points.begin();
  j = raw_points.end();
  j--;
  if (pointsAreIdentical(*i,*j))
  {
    raw_points.erase(i);
  }
}

/** Generate approximately uniform sampling of path length.
 *  Note path length does not vary linearly with path parameter
 *  due to nonuniform spacing of control points.
 *  Hence, approximately uniform.
 * \param[in] num_splines Number of splines in contour.
 * \param[in] target_num_samples Target number samples in contour.
 */

void Contour::uniformSampling (const int & num_splines,
                               const int & target_num_samples)
{
  Controls & cs (Controls::instance()); 
  double total_path_length = calculateTotalPathLength();
  double dL = total_path_length/static_cast<double>(target_num_samples);
  double cumulative_s = 0.0;
  double cumulative_path_length = 0.0;
  // for each spline
  for (int i=0;i<num_splines;i++)
  {
    double a=spline_lengths[i];
    // while path length remainging in current spline
    while (a>cs.getEpsilon())
    {
      double diff = dL-cumulative_path_length;
      // if remainging path length is larger than what we need
      if (diff<a || !distinguishable(diff,a))
      {
        // then finalize sample
        cumulative_s += diff/spline_lengths[i];
        a -= diff;
        double my_s = 1.0-a/spline_lengths[i];
        s.add(my_s,i,cumulative_s);
        cumulative_path_length = 0.0;
        cumulative_s = 0.0;
      }
      // else remainging path length is smaller than or equal to what we need
      else
      {
        // then incorporate remaining spline
        cumulative_s += a/spline_lengths[i];
        cumulative_path_length += a;
        a = 0.0;
      }
    }
  }
  if (s.getNumSamples()!=target_num_samples)
  {
    cout << "Contour::uniformSampling: "
         << "s.getNumSamples = " << s.getNumSamples()
         << ", target_num_samples = " << target_num_samples
         << endl;
    assert(s.getNumSamples()==target_num_samples);
  }
}

/** Calculate length of cubic spline from start of spline
 *  to path parameter value by numerical integration.
 * \param[in] spline_index Index of spline.
 * \param[in] path_parameter_value Path parameter. 0<t<1.
 * \return Spline length.
 */

double Contour::calculateSplineLength (const int & spline_index,
                                       const double & path_parameter_value) const
{
  Controls & cs(Controls::instance());
  if (fabs(path_parameter_value)<cs.getEpsilon()) return 0.0;
  double parameter_step = path_parameter_value*cs.getIntegrationStepFactor();
  double integral = 0.0;
  for (int i=0;i<cs.getIntegrationStep();i+=2)
  {
    double j = static_cast<double>(i);
    double a = evalSprime((j+0.0)*parameter_step,spline_index);
    double b = evalSprime((j+1.0)*parameter_step,spline_index);
    double c = evalSprime((j+2.0)*parameter_step,spline_index);
    integral += a+4.0*b+c;
  }
  return integral*integration_ratio*parameter_step;
}

/** Get length of cubic spline from start of spline
 *  to path parameter value by numerically integration.
 * \param[in] spline_index Index of spline.
 * \param[in] path_parameter_value Path parameter. 0<t<1.
 * \return Spline length.
 */

double Contour::getSplineLength (const int & spline_index,
                                 const double & path_parameter_value) const
{
  assert(spline_index>=0);
  assert(spline_index<static_cast<int>(spline_lengths.size()));
  // if want path length of entire spline
  if (distinguishable(path_parameter_value,1.0)==false)
  {
    assert(spline_lengths[spline_index]>0.0);
    return spline_lengths[spline_index];
  }
  // else want path length of portion of spline
  else
  {
    return calculateSplineLength(spline_index,path_parameter_value);
  }
}

/** Calculate length of cubic splines between two spline sample points.
 * \param[in] lhs_sample_index Index of first sample point.
 * \param[in] rhs_sample_index Index of second sample point.
 * \return Spline length between sample points.
 */

double Contour::getInterSampleSplineLength (const int & lhs_sample_index,
                                            const int & rhs_sample_index) const
{
  Point const * const p1 = s.getPoint(lhs_sample_index);
  Point const * const p2 = s.getPoint(rhs_sample_index);
  return sqrt(p1->distance_squared(*p2));
}

/** Calculate energy of proximity between two spline sample points.
 * \param[in] lhs_sample_index Index of first sample point.
 * \param[in] rhs_sample_index Index of second sample point.
 * \return Energy of proximity between sample points. 
 */

double Contour::calculateProximityEnergy (const int & lhs_sample_index,
                                          const int & rhs_sample_index) const
{
  Controls & cs (Controls::instance()); 
  double arc_length = fabs(getInterSampleSplineLength(lhs_sample_index,
                                                      rhs_sample_index)); 
#ifndef NOASSERT
  assert(arc_length>0.0);
#endif
  double e = cs.getProximityEnergyGain();
  for (int i=0;i<cs.getProximityEnergyExponent();i++)
  {
    e /= arc_length;
    e += 1.0e4*arc_length;
  }
  return e;
}

/** Calculate energy of cubic spline at sample point.
 * \param[in] sample_index Index of sample point.
 * \return Energy of cubic spline sample.
 */

double Contour::getEnergy (const int & sample_index) const
{
  Controls & cs (Controls::instance()); 
  assert(s.getNumSamples()>2);
  int i0 = sample_index-1;
  int i1 = sample_index;
  int i2 = sample_index+1;
  if (i0<0) i0 = s.getNumSamples()-1;
  if (i2==s.getNumSamples()) i2 = 0;
  // E1 = 1/2*k*radius_of_curvature^2
  // add_offset flag is true so we only get positive values
  //double p1 = evalRadiusOfCurvature(i1,true);
  double p1 = getFinalRadiusOfCurvature(i1);
  double E1 = cs.getCurvatureEnergyGain();
  for (int i=0;i<cs.getCurvatureEnergyExponent();i++) E1 *= p1;
  double E2 = calculateProximityEnergy(i0,i1);
  double E3 = calculateProximityEnergy(i1,i2);
  return E1+E2+E3;
}

/** Calculate cumulative energy of cubic spline at all sample points.
 * \return Cumulative energy of all cubic spline samples.
 */

double Contour::getTotalEnergy (void) const
{
  // Note that the proximity energy is included twice
  // but all we care about for simulated annealing is change in energy,
  // so no worries.
  double energy = 0.0;
  // for each spline sample point
  for (int sample_index=0;sample_index<s.getNumSamples();sample_index++)
  {
    energy += getEnergy(sample_index);
  }
  return energy;
}

/** Determine if change in sample point path parameter value will
 *  reorder sample points by path parameter values.
 * \param[in] sample_index Index of sample point.
 * \param[in] parameter_value_change Change in path parameter value.
 * \return True if sample point reorder occurs; false otherwise.
 */

bool Contour::sampleJumpedNeighbor (const int & sample_index,
                                    const double & parameter_value_change) const
{
  if (parameter_value_change < 0.0)
  {
    return s.getLhsS(sample_index)<fabs(parameter_value_change);
  }
  else
  {
    int rhs_sample_index = sample_index+1;
    if (rhs_sample_index==s.getNumSamples()) rhs_sample_index = 0;
    return s.getLhsS(rhs_sample_index)<parameter_value_change;
  }
}

/** Move sample points along splines to minimize
 *  cumulative energy of sample points.
 * \param[in] SA Simulated annealing class.
 * \param[in] num_splines Number of splines in contour.
 */

void Contour::simAnnealing (Sim_Anneal & SA, const int & num_splines)
{
  Controls & cs (Controls::instance()); 
  double total_energy = 0.0;
  if (cs.getPrintDetailedInfo())
  {
    SA.reserveDetailedInfo();
    total_energy = getTotalEnergy();
  }
  SA.setScale(s.getNumSamples());
  while (SA.isFrozen()==false)
  {
    while (SA.moreMoves(s.getNumSamples()))
    {
      // pick spline sample at random
      int sample_index = SA.getRandomSample(s.getNumSamples());
      // save current sample info
      double orig_s = s.getParameterValue(sample_index);
      int orig_spline_index = s.getSplineIndex(sample_index);
      // pick random displacement of sample point along spline
      double s_disp = SA.getRandomPathParameterDisplacement();
      // get new sample info based on displacement
      double new_s = orig_s+s_disp;
      int new_spline_index = orig_spline_index;
      if (new_s<0.0)
      {
        new_s += 1.0;
        new_spline_index = orig_spline_index-1;
        if (new_spline_index<0) new_spline_index = num_splines-1;
      }
      else if (new_s>1.0)
      {
        new_s -= 1.0;
        new_spline_index = orig_spline_index+1;
        if (new_spline_index>=num_splines) new_spline_index = 0;
      }
#ifndef NOASSERT
      assert((new_spline_index+new_s)<cp.getNumSplines());
#endif
      // detect samples crossing
      if (sampleJumpedNeighbor(sample_index,s_disp))
      {
        SA.recordUnsuccessfulMove();
        continue;
      }
      // compute energy at original location
      double energy0 = getEnergy(sample_index);
      // move sample point
      s.set(new_s,new_spline_index,sample_index,s_disp,sampleSpline(new_s,new_spline_index));
      // compute energy at new location
      double energy1 = getEnergy(sample_index);
      // if energy drops
      if (energy1<energy0)
      {
        //then keep move
        SA.recordSuccessfulMove();
        if (cs.getPrintDetailedInfo())
        {
          total_energy += energy1-energy0;
          SA.recordDetailedInfo(total_energy);
        }
      }
      // if energy increases or no change
      else
      {
        // then might keep move
        if (SA.keepMove(energy1-energy0))
        {
          SA.recordSuccessfulMove();
          if (cs.getPrintDetailedInfo())
          {
            total_energy += energy1-energy0;
            SA.recordDetailedInfo(total_energy);
          }
        }
        // or move sample point back
        else
        {
          SA.recordUnsuccessfulMove();
          s.set(orig_s,orig_spline_index,sample_index,-s_disp,sampleSpline(orig_s,orig_spline_index));
        }
      }
    }
    //if (cs.getPrintDetailedInfo())
    //{
    //  SA.debug(name,section,num_samples);
    //}
    // lower temperature and prepare for another set of sample moves
    SA.decreaseTemp();
  }
  if (cs.getPrintDetailedInfo())
  {
    char filename1[256];
    sprintf(filename1,"%s%s_%i_sim_anneal_energy.log",
            cs.getOutputDir(),getName(),getSection());
    SA.printEnergy(filename1);
    sprintf(filename1,"%s_%i_sim_anneal_energy.log",
            getName(),getSection());
    char filename2[256];
    sprintf(filename2,"%s%s_%i_sim_anneal_temp.log",
            cs.getOutputDir(),getName(),getSection());
    SA.printTemperature(filename2);
    sprintf(filename2,"%s_%i_sim_anneal_temp.log",getName(),getSection());
    SA.printScript(filename1,filename2,getName(),getSection());
  }
}

void Contour::filterCurvature (void)
{
  Controls & cs (Controls::instance()); 
  const int order=3;
  const int rotate=0;
  const int multiplicity=3;
  // THIS COULD BE IMPROVED BY NOT COPYING radius_of_curvature
  std::vector<double> input;
  input.reserve(multiplicity*radius_of_curvature.size()+order);
  for (int i=0;i<order;i++) input.push_back(0.0);
  for (int i=0;i<multiplicity;i++)
  {
    input.insert(input.end(),radius_of_curvature.begin(),radius_of_curvature.end());
  }

  std::vector<double> output_first_pass;
  output_first_pass.reserve(multiplicity*radius_of_curvature.size()+order);
  for (int i=0;i<order;i++) output_first_pass.push_back(0);

  // 75 Hz cutoff
  double num[order+1] = {0.031689343849711,  0.095068031549133,  0.095068031549133,  0.031689343849711};
  double den[order+1] = {1.000000000000000, -1.459029062228061,  0.910369000290069, -0.197825187264319};
  // 50 Hz cutoff
  //double num[order+1] = {0.0113,0.0340,0.0340,0.0113};
  //double den[order+1] = {1.0000,-1.9630,1.4000,-0.3464};
  // 15 Hz cutoff
  //double num[order+1] = {0.000416546139076,0.001249638417227,0.001249638417227,0.000416546139076};
  //double den[order+1] = {1.000000000000000,-2.686157396548143,2.419655110966473,-0.730165345305723};
  //double num[7] = {0.0175*1e-5,0.1052*1e-5,0.2630*1e-5,0.3507*1e-5,0.2630*1e-5,0.1052*1e-5,0.0175*1e-5};

  for (int i=order;i<static_cast<int>(input.size());i++)
  {
    double myinput  = num[0]*input[i] +
                      num[1]*input[i-1] +
                      num[2]*input[i-2] +
                      num[3]*input[i-3];
    double myoutput = den[1]*output_first_pass[i-1] +
                      den[2]*output_first_pass[i-2] +
                      den[3]*output_first_pass[i-3];
    //double myinput = num[0]*input[i]   + num[1]*input[i-1] +
    //                 num[2]*input[i-2] + num[3]*input[i-3] +
    //                 num[4]*input[i-4] + num[5]*input[i-5] +
    //                 num[6]*input[i-6];
    //double myoutput = den[1]*output_first_pass[i-1] +
    //                  den[2]*output_first_pass[i-2] +
    //                  den[3]*output_first_pass[i-3] +
    //                  den[4]*output_first_pass[i-4] +
    //                  den[5]*output_first_pass[i-5] +
    //                  den[6]*output_first_pass[i-6];
    output_first_pass.push_back(myinput-myoutput);
    assert((myinput-myoutput)<100.0);
    // DEBUG
    if ((myinput-myoutput)<-100.0)
    {
      cout << "Contour::filterCurvature "
            << "myinput = " << myinput
            << ", myoutput = " << myoutput
            << endl;
    }
    // DEBUG
  }
  for (int i=0;i<order;i++) output_first_pass.push_back(0);

  std::vector<double> output_second_pass;
  output_second_pass.reserve(multiplicity*radius_of_curvature.size()+order);
  for (int i=0;i<order;i++) output_second_pass.push_back(0);
  int i=order;
  double min = 1e30;
  double max = -1e30;
  for (int j=output_first_pass.size()-1-order;j>(order-1);j--)
  {
    double myinput = num[0]*output_first_pass[j] +
                     num[1]*output_first_pass[j+1] +
                     num[2]*output_first_pass[j+2] +
                     num[3]*output_first_pass[j+3];
    double myoutput = den[1]*output_second_pass[i-1] +
                      den[2]*output_second_pass[i-2] +
                      den[3]*output_second_pass[i-3];
    //double myinput = num[0]*output_first_pass[j] +
    //                 num[1]*output_first_pass[j+1] +
    //                 num[2]*output_first_pass[j+2] +
    //                 num[3]*output_first_pass[j+3] +
    //                 num[4]*output_first_pass[j+4] +
    //                 num[5]*output_first_pass[j+5] +
    //                 num[6]*output_first_pass[j+6];
    //double myoutput = den[1]*output_second_pass[i-1] +
    //                  den[2]*output_second_pass[i-2] +
    //                  den[3]*output_second_pass[i-3] +
    //                  den[4]*output_second_pass[i-4] +
    //                  den[5]*output_second_pass[i-5] +
    //                  den[6]*output_second_pass[i-6];
    double diff = myinput-myoutput;
    output_second_pass.push_back(diff);
    if (diff>max) max = diff;
    if (diff<min) min = diff;
    i++;
  }
  double amp = max-min;
  double scale = 1.0;
  if (amp>3.0)
  {
    scale = 3.0/amp;
  }
  const double myscale = scale;
  final_radius_of_curvature.clear();
  final_radius_of_curvature.reserve(radius_of_curvature.size());
  const int shift = multiplicity-1;
  i = shift*radius_of_curvature.size()+shift-radius_of_curvature.size()+rotate;
  for (int j=0;j<rotate;j++)
  {
    final_radius_of_curvature.push_back(myscale*output_second_pass[i]);
    i--;
  }
  i = shift*radius_of_curvature.size()+shift;
  for (int j=0;j<static_cast<int>(radius_of_curvature.size())-rotate;j++)
  {
    assert(i>=0);
    if (i>=static_cast<int>(output_second_pass.size()))
    {
      cout << "Contour::filterCurvature "
            << "i = " << i
            << ", output_second_pass.size = "
            << output_second_pass.size()
            << endl;
    }
    assert(i<static_cast<int>(output_second_pass.size()));
    final_radius_of_curvature.push_back(myscale*output_second_pass[i]);
    i--;
  }
  if (static_cast<int>(final_radius_of_curvature.size())!=s.getNumSamples())
  {
    cout << "Contour::filterCurvature "
          << "final_radius_of_curvature.size = "
          << final_radius_of_curvature.size()
          << ", s.getNumSamples() = "
          << s.getNumSamples() << endl;
  }
  assert(static_cast<int>(final_radius_of_curvature.size())==s.getNumSamples());
  if (cs.getPrintDetailedInfo())
  {
    char filename[256];
    sprintf(filename,"%s%s_%i_radius_curvature_unfiltered.log",cs.getOutputDir(),getName(),getSection());
    printRadiusOfCurvature(filename,radius_of_curvature);
    sprintf(filename,"%s%s_%i_radius_curvature_filtered.log",cs.getOutputDir(),getName(),getSection());
    printRadiusOfCurvature(filename,final_radius_of_curvature);
  }
}

/** Calculate radius of curvature of cubic splines at all path parameter samples.
 * \param[in] target_num_samples Target number samples in contour.
 * \param[out] radius_of_curvature Radius of curvature of cubic splines. 
 */

void Contour::calculateRadiusOfCurvature (const int & target_num_samples)
{
  radius_of_curvature.clear();
  radius_of_curvature.reserve(s.getNumSamples());
  // calculate vector size
  double my_t = s.getParameterValue(s.getNumSamples()-1);
  my_t += s.getSplineIndex(s.getNumSamples()-1);
  int max = floor(my_t*10.0)+4;
  if (max<=target_num_samples) max = target_num_samples+1;
  assert(max>target_num_samples);
  path_parameter.clear();
  uniform_path_parameter.clear();
  uniform_path_parameter.reserve(s.getNumSamples());
  path_parameter.insert(path_parameter.end(),max,-1.0);
  // for each spline sample point
  for (int sample_index=0;sample_index<s.getNumSamples();sample_index++)
  {
    // add_offset flag is false, since we do not know the offset value yet
    double a = evalRadiusOfCurvature(sample_index,false);
    radius_of_curvature.push_back(a);
    double t = s.getParameterValue(sample_index);
    t += s.getSplineIndex(sample_index);
    int my_index = floor(t*10.0);
    // NOTE COULD BE COLLISIONS
#ifndef NOASSERT
    assert(my_index<static_cast<int>(path_parameter.size()));
#endif
    path_parameter[my_index]=sample_index;
    uniform_path_parameter.push_back(t);
  }
  int last_sample_index = -1.0;
  for (c_r_v_i_iterator i = path_parameter.rbegin();i!=path_parameter.rend();i++)
  {
    if (*i>0.0)
    {
      last_sample_index = *i;
      break;
    }
  }
  assert(last_sample_index>0.0);
  // for each path parameter
  for (v_i_iterator i = path_parameter.begin();i!=path_parameter.end();i++)
  {
    while (*i<0.0)
    {
      *i = last_sample_index;
    }
    last_sample_index = *i;
  }
  filterCurvature();
}

/** Identify and print maximum deviation of spline samples
 *  to linearly-interpolated raw contour points.
 */

void Contour::printMaxDeviation (void) const
{
  double deviation_threshold = Controls::instance().getDeviationThreshold();
  int count = 0;
  double max_d = -1.0;
  c_d_iterator i = s.getFirstDeviationConst();
  while (i!=s.getOnePastLastDeviationConst())
  {
    if (*i>max_d) max_d = *i;
    if (*i>deviation_threshold) count++;
    i++;
  }
  cout << "Contour::printMaxDeviation: "
       << "max deviation = " << max_d
       << ", number of hyperthreshold deviations = " << count
       << endl;
}

/** Calculate straight-line distance between successive spline samples.
 * \param[out] my_path_intervals Inter-sample distances.
 */

void Contour::calcSampleIntervals (vec_d & my_path_intervals)
{
  my_path_intervals.clear();
  c_p_iterator i = s.getFirstSplineSampleConst();
  c_p_iterator j = i;
  j++;
  int k = 0;
  while (j!=s.getOnePastLastSplineSampleConst())
  {
    assert(i!=j);
    double a = sqrt((*i).distance_squared(*j));
    my_path_intervals.push_back(a);
    i++;
    j++;
    k++;
  }
  j = s.getFirstSplineSampleConst();
  double a = sqrt((*i).distance_squared(*j));
  my_path_intervals.push_back(a);
}

/** Calculate straight line distance between adjacent sample
 *  points, store, and update sample interval statistics.
 * \param[out] h Instance of Histogram class for sample intervals.
 */

void Contour::updateSampleIntervals (Histogram & h)
{
  c_p_iterator i = s.getFirstSplineSampleConst();
  c_p_iterator j = i;
  j++;
  while (j!=s.getOnePastLastSplineSampleConst())
  {
    assert(i!=j);
    double a = sqrt((*i).distance_squared(*j));
    h.update(a);
    i++;
    j++;
  }
  j = s.getFirstSplineSampleConst();
  double a = sqrt((*i).distance_squared(*j));
  h.update(a);
}

/** Update statistics of deviation of spline samples
 *  to linearly-interpolated raw contour points.
 * \param[out] h Instance of Histogram class for deviations.
 */

void Contour::updateDeviations (Histogram & h) const
{
  c_d_iterator i = s.getFirstDeviationConst();
  while (i!=s.getOnePastLastDeviationConst())
  {
    h.update(*i);
    i++;
  }
}

/** Duplicate control points of splines to reduce deviation
 *  between spline samples and linearly-interpolated
 *  raw contour points.
 * \param[out] control_points Indices of raw contour points
 *  used as control points for splines.
 */

void Contour::duplicateControlPoints (list_i & control_points) const
{
  std::vector<int> raw_points_to_duplicate;
  // identify raw points to duplicate as control points
  s.getRawPointsToDuplicate(raw_points,raw_points_to_duplicate);
  assert(static_cast<int>(raw_points_to_duplicate.size())>0);
  // generate control points vector
  sort(raw_points_to_duplicate.begin(),raw_points_to_duplicate.end());
  std::vector<int>::iterator j = raw_points_to_duplicate.begin();
  i_iterator k,z,i=control_points.begin();
  while (i!=control_points.end())
  {
    k=i;
    k++;k++;
    if (j!=raw_points_to_duplicate.end() && *i==*j)
    {
      if (*k!=*j)
      {
        control_points.insert(i,*j);
        // update tail
        if (*j==0 && raw_points.size()>3)
        {
          z = control_points.end();
          z--;
          if (*z==1 || *z==2) control_points.erase(z);
          if (z==control_points.end()) z--;
          while (*z!=*j) z--;
          control_points.insert(z,*j);
        }
      }
      else
      {
        i++;
        i++;
        i++;
      }
      j++;
    }
    else
    {
      i++;
    }
  }
  if (j!=raw_points_to_duplicate.end())
  {
    cout << "\n\nContour::duplicateControlPoints: "
          << "Error: raw_points_to_duplicate not exhausted.\n"
          << "Contour::duplicateControlPoints: "
          << "raw_points.size = " << raw_points.size()
          << ", raw_points_to_duplicate.size = "
          << raw_points_to_duplicate.size()
          << ", j = " << *j << endl << endl;
    for (j=raw_points_to_duplicate.begin();j!=raw_points_to_duplicate.end();j++)
    {
      cout << "j = " << *j << endl;
    }
    cout << endl;
    exit(1);
    assert(j==raw_points_to_duplicate.end());
  }
}

/** Fit splines to contour points in this contour.
 * \param[out] my_path_intervals Straight-line distance
 * between successive spline samples.
 * \return Expected number of spline samples based on contour
 * spline lengths and minimum and maximum sample intervals.
 */

int Contour::processContour (vec_d & my_path_intervals)
{
  Controls & cs (Controls::instance()); 
  int num_splines = cp.getNumSplines();
  initCurveLengths(num_splines);
  // possible values are >2 and 0.
  int target_num_samples = calculateNumSamples(num_splines);
  // Do not bother annealing contours with less than three points.
  if (target_num_samples<3) return target_num_samples;
  s.reserve(target_num_samples);
  uniformSampling(num_splines,target_num_samples);
  sampleSplines();
  calculateRadiusOfCurvature(target_num_samples);
  setOffset();
  if (cs.getPrintDetailedInfo())
  {
    char filename[256];
    sprintf(filename,"%s%s_%i_param_s.log",cs.getOutputDir(),getName(),getSection());
    s.print(filename);
    sprintf(filename,"%s%s_%i_uniform_sample.log",cs.getOutputDir(),getName(),getSection());
    printSplineSample(filename);
    sprintf(filename,"%s%s_%i_radius_curvature.log",cs.getOutputDir(),getName(),getSection());
    printRadiusOfCurvature(filename,final_radius_of_curvature);
  }
  calcSampleIntervals(my_path_intervals);
  Sim_Anneal SA;
  simAnnealing(SA,num_splines);
  sampleSplines();
  return target_num_samples;
}

/** Write diagnostic information to file.
 */

void Contour::writeFilesAfter (void) const
{
  Controls & cs (Controls::instance()); 
  char filename[256];
  // print spline sample points
  sprintf(filename,"%s%s_%i_sample.log",cs.getOutputDir(),getName(),getSection());
  printSplineSample(filename);
  // print arc lengths
  sprintf(filename,"%s%s_%i_arclengths.log",cs.getOutputDir(),getName(),getSection());
  printArcLengths(filename);
  // print s parameter values
  sprintf(filename,"%s%s_%i_param_s.log",cs.getOutputDir(),getName(),getSection());
  s.print(filename);
  // print starting point of splines
  sprintf(filename,"%s%s_%i_spline_start.log",cs.getOutputDir(),getName(),getSection());
  printSplineSamplesHere(filename,cp.getNumSplines(),0.0);
  // print ending point of splines
  sprintf(filename,"%s%s_%i_spline_end.log",cs.getOutputDir(),getName(),getSection());
  printSplineSamplesHere(filename,cp.getNumSplines(),1.0);
  // print raw contour point indices of spline control points
  sprintf(filename,"%s%s_%i_control.log",cs.getOutputDir(),getName(),getSection());
  cp.print(filename);
}

/** Calculate location on spline.
 * \param[in] path_parameter_value Parameter value describing position along spline. 
 *  Value must be between 0.0 and 1.0, inclusive.
 * \param[in] spline_index Index of spline.
 * \return Location on spline of given index corresponding to path parameter.
 */

Point Contour::sampleSpline (const double & path_parameter_value,
                             const int & spline_index) const
{
    const double s1 = path_parameter_value; 
    const double s2 = s1*s1;
    const double s3 = s2*s1;
    double a = s3*getAx(spline_index)+s2*getBx(spline_index)+
               s1*getCx(spline_index)+   getDx(spline_index);
    double b = s3*getAy(spline_index)+s2*getBy(spline_index)+
               s1*getCy(spline_index)+   getDy(spline_index);
    return Point(a,b);
}

/** Calculate the x and y values of cubic splines at each
 *  path parameter value.
 */

void Contour::sampleSplines (void)
{
  s.clearAndReserveSamples();
  // for each sample point
  for (int sample_index=0;sample_index<s.getNumSamples();sample_index++)
  {
    double s1 = s.getParameterValue(sample_index);
    int i = s.getSplineIndex(sample_index);
    s.addSplineSample(sampleSpline(s1,i));
  }
}

/** Write input contour points to file.
 * \param[in] filename Output file name.
 */

void Contour::printRawPoints (char const * const filename) const
{
  FILE *F = fopen(filename,"a");
  if (!F) { printf("Could not open output file %s\n",filename); return; }
  char line[2048];
  sprintf(line,"# num raw points = %d\n",static_cast<int>(raw_points.size()));
  fputs(line,F);
  // for each point
  for (c_p_iterator i=raw_points.begin();i!=raw_points.end();i++)
  {
    sprintf(line,"%.15g %.15g\n",i->getX(),i->getY());
    fputs(line,F);
  }
  fclose(F);
}

/** Write radius of curvature to file.
 * \param[in] filename Output file name.
 * \param[in] radii Radius of curvature of cubic splines.
 */

void Contour::printRadiusOfCurvature (char const * const filename,
                                      const vec_d & radii) const
{
  FILE *F = fopen(filename,"a");
  if (!F) { printf("Could not open output file %s\n",filename); return; }
  char line[2048];
  sprintf(line,"# num parameter s samples = %d\n",s.getNumSamples());
  fputs(line,F);
  // prepare path parameter
  int index = 0;
  // for each value of parameter s (fraction of cubic curve segment)
  std::vector<double>::const_iterator i=radii.begin();
  for (int sample_index=0;sample_index<s.getNumSamples();sample_index++)
  {
    sprintf(line,"%g %g\n",uniform_path_parameter[sample_index],*i);
    fputs(line,F);
    i++;
    index++;
  }
  fclose(F);
}

/** Write length of cubic splines between pairs of sequential
 *  samples of path parameter.
 * \param[in] filename Output file name.
 */

void Contour::printArcLengths (char const * const filename) const
{
  FILE *F = fopen(filename,"a");
  if (!F) { printf("Could not open output file %s\n",filename); return; }
  char line[2048];
  sprintf(line,"# num spline points = %d\n",s.getNumSamples());
  fputs(line,F);
  // for each value of parameter s (fraction of cubic curve segment)
  for (int sample_index=0;sample_index<s.getNumSamples()-1;sample_index++)
  {
    double arc_length = fabs(getInterSampleSplineLength(sample_index,sample_index+1)); 
    sprintf(line,"%g\n",arc_length);
    fputs(line,F);
  }
  double arc_length = fabs(getInterSampleSplineLength(s.getNumSamples()-1,0)); 
  sprintf(line,"%g\n",arc_length);
  fputs(line,F);
  fclose(F);
}

/** Write spline sample point evaluated at given parameter value.
 * \param[in] filename Output file name.
 * \param[in] num_splines Number of splines in contour.
 * \param[in] path_parameter_value Path parameter value of interest.
 *  Value must be between 0.0 and 1.0, inclusive.
 */

void Contour::printSplineSamplesHere (char const * const filename,
                                      const int & num_splines,
                                      double path_parameter_value) const
{
  FILE *F = fopen(filename,"a");
  if (!F) { printf("Could not open output file %s\n",filename); return; }
  char line[2048];
  sprintf(line,"# num spline points = %d\n",s.getNumSamples());
  fputs(line,F);
  // for each spline
  for (int i=0;i<num_splines;i++)
  {
    Point a = sampleSpline(path_parameter_value,i);
    char line2[2048];
    sprintf(line2,"%g %g\n",a.getX(),a.getY());
    fputs(line2,F);
  }
  fclose(F);
}

/** Write x and y values of spline samples to file.
 * \param[in] filename Output file name.
 */

void Contour::printSplineSample (char const * const filename) const
{
  FILE *F = fopen(filename,"a");
  if (!F) { printf("Could not open output file %s\n",filename); return; }
  char line[2048];
  sprintf(line,"# num spline points = %d\n",s.getNumSamples());
  fputs(line,F);
  // for each point
  c_p_iterator i=s.getFirstSplineSampleConst();
  while (i!=s.getOnePastLastSplineSampleConst())
  {
    i->print(line);
    char line2[2048];
    sprintf(line2,"%s\n",line);
    fputs(line2,F);
    i++;
  }
  fclose(F);
}

/** Write input contour points to file in pts format.
 * \param[in] myname Contour name.
 */

void Contour::printRawPtsFiles (char const * myname) const
{
  Controls & cs (Controls::instance()); 
  char filename[256];
  sprintf(filename,"%sinput_%s%d.pts",cs.getOutputDir(),myname,section);
  FILE * F = fopen(filename,"a");
  if(!F){printf("Couldn't open output file %s\n",filename);exit(1);}
  assert(F);
  char line[2048];
  sprintf(line,"%i\n",static_cast<int>(raw_points.size()));
  fputs(line,F);
  // for each interpolated point in contour
  for (c_p_iterator j=raw_points.begin();j!=raw_points.end();j++)
  {
    sprintf(line,"%.15g %.15g %.15g\n",j->getX()*cs.getOutputScaleFactor(),
            j->getY()*cs.getOutputScaleFactor(),
            section*cs.getSectionThickness()*cs.getOutputScaleFactor());
    fputs(line,F);
  }
  fclose(F);
}

/** Write x and y values of spline samples to file in pts format.
 * \param[in] myname Contour name.
 */

void Contour::printPtsFiles (char const * myname) const
{
  Controls & cs (Controls::instance()); 
  char filename[256];
  sprintf(filename,"%s%s%d.pts",cs.getOutputDir(),myname,getSection());
  FILE * F = fopen(filename,"a");
  if (!F) { printf("Couldn't open output file %s\n",filename);exit(1); }
  assert(F);
  char line[2048];
  sprintf(line,"%i\n",s.getNumSamplePoints());
  fputs(line,F);
  double z = getSection()*cs.getSectionThickness();
  // for each spline sample in contour
  c_p_iterator j=s.getFirstSplineSampleConst();
  while (j!=s.getOnePastLastSplineSampleConst())
  {
    sprintf(line,"%.15g %.15g %.15g\n",j->getX()*cs.getOutputScaleFactor(),
            j->getY()*cs.getOutputScaleFactor(),
            z*cs.getOutputScaleFactor());
    fputs(line,F);
    j++;
  }
  fclose(F);
}

/** Write grace scripts for visualizing contours and spline samples.
 */

void Contour::printSimpleScripts (void) const
{
  Controls & cs (Controls::instance()); 
  char filename[256],line[2048];
  FILE *F;
  // print grace script
  sprintf(filename,"%sbfile_%s_%i",cs.getOutputDir(),getName(),getSection());
  F = fopen(filename,"w");
  if (!F) { printf("Could not open output file %s\n",filename); return; }
  sprintf(line,"#Obligatory descriptive comment\n");
  fputs(line,F);
  sprintf(line,"READ xy \"%s_%i_raw.log\"\n",getName(),getSection());
  fputs(line,F);
  sprintf(line,"READ xy \"%s_%i_sample.log\"\n",getName(),getSection());
  fputs(line,F);
  sprintf(line,"legend on\ns0 legend ");
  sprintf(line,"%s\"Raw\"\n",line);
  sprintf(line,"%ss1 legend \"Samples\"\n",line);
  fputs(line,F);
  sprintf(line,"s0 symbol 1\ns0 symbol color 1\ns0 symbol fill color 1\n");
  sprintf(line,"%ss0 line type 0\ns0 line color 1\ns0 errorbar color 1\n",line);
  fputs(line,F);
  sprintf(line,"s1 symbol 2\ns1 symbol color 2\ns1 symbol fill color 2\n");
  sprintf(line,"%ss1 line type 0\ns1 line color 2\ns1 errorbar color 2\n",line);
  sprintf(line,"%ss1 symbol size 0.330000\n",line);
  fputs(line,F);
  sprintf(line,"legend 0.85, 0.99\n");
  fputs(line,F);
  fclose(F);
  // print shell script
  sprintf(filename,"%s%s_%i.sh",cs.getOutputDir(),getName(),getSection());
  F = fopen(filename,"w");
  if (!F) { printf("Could not open output file %s\n",filename); return; }
  sprintf(line,"#!/bin/bash\n");
  fputs(line,F);
  sprintf(line,"xmgrace -world 0 0 10 10 -autoscale none -noask -nosafe -batch bfile_%s_%i\n",getName(),getSection());
  fputs(line,F);
  fclose(F);
}

/** Write grace scripts for visualizing contours and spline samples.
 */

void Contour::printScripts (void) const
{
  Controls & cs (Controls::instance()); 
  char filename[256],line[2048];
  FILE *F;
  // print grace script
  sprintf(filename,"%sbfile_%s_%i",cs.getOutputDir(),getName(),getSection());
  F = fopen(filename,"w");
  if (!F) { printf("Could not open output file %s\n",filename); return; }
  sprintf(line,"#Obligatory descriptive comment\n");
  fputs(line,F);
  sprintf(line,"view 0.100000, 0.100000, 0.950000, 0.950000\n");
  fputs(line,F);
  //sprintf(line,"world 0, 0, 10, 10\n");
  //sprintf(line,"xmin 0\nymin 0\nxmax 10\nymax 10\n");
  //fputs(line,F);
  sprintf(line,"READ xy \"%s_%i_raw.log\"\n",getName(),getSection());
  fputs(line,F);
  sprintf(line,"READ xy \"%s_%i_uniform_sample.log\"\n",getName(),getSection());
  fputs(line,F);
  sprintf(line,"READ xy \"%s_%i_sample.log\"\n",getName(),getSection());
  fputs(line,F);
  //sprintf(line,"READ xy \"%s_%i_spline_start.log\"\n",getName(),getSection());
  //fputs(line,F);
  //sprintf(line,"READ xy \"%s_%i_spline_end.log\"\n",getName(),getSection());
  //fputs(line,F);
  sprintf(line,"legend on\ns0 legend ");
  sprintf(line,"%s\"Raw\"\n",line);
  sprintf(line,"%ss1 legend \"Uniform Sample\"\n",line);
  fputs(line,F);
  sprintf(line,"s0 symbol 1\ns0 symbol color 1\ns0 symbol fill color 1\n");
  sprintf(line,"%ss0 line type 0\ns0 line color 1\ns0 errorbar color 1\n",line);
  fputs(line,F);
  sprintf(line,"s1 symbol 2\ns1 symbol color 2\ns1 symbol fill color 2\n");
  sprintf(line,"%ss1 line type 0\ns1 line color 2\ns1 errorbar color 2\n",line);
  sprintf(line,"%ss1 symbol size 0.330000\n",line);
  fputs(line,F);
  sprintf(line,"s2 symbol 3\ns2 symbol color 3\ns2 symbol fill color 3\n");
  sprintf(line,"%ss2 line type 0\ns2 line color 3\ns2 errorbar color 3\n",line);
  sprintf(line,"%ss2 symbol size 0.330000\n",line);
  fputs(line,F);
  sprintf(line,"s2 legend \"Anneal Sample\"\n");
  fputs(line,F);
  sprintf(line,"s3 symbol 4\ns3 symbol color 4\ns3 symbol fill color 4\n");
  sprintf(line,"%ss3 line type 0\ns3 line color 4\ns3 errorbar color 4\n",line);
  sprintf(line,"%ss3 symbol size 0.330000\n",line);
  fputs(line,F);
  sprintf(line,"s3 legend \"Spline Start\"\n");
  fputs(line,F);
  sprintf(line,"s4 symbol 5\ns4 symbol color 5\ns4 symbol fill color 5\n");
  sprintf(line,"%ss4 line type 0\ns4 line color 5\ns4 errorbar color 5\n",line);
  sprintf(line,"%ss4 symbol size 0.330000\n",line);
  fputs(line,F);
  sprintf(line,"s4 legend \"Spline End\"\n");
  fputs(line,F);
  sprintf(line,"legend 0.85, 0.99\n");
  fputs(line,F);
  fclose(F);
  // print shell script
  sprintf(filename,"%s%s_%i.sh",cs.getOutputDir(),getName(),getSection());
  F = fopen(filename,"w");
  if (!F) { printf("Could not open output file %s\n",filename); return; }
  sprintf(line,"#!/bin/bash\n");
  fputs(line,F);
  sprintf(line,"xmgrace -world 0 0 10 10 -autoscale none -noask -nosafe -batch bfile_%s_%i\n",getName(),getSection());
  fputs(line,F);
  fclose(F);
}
