// Author: Justin Kinney
// Date: July 2009

#ifndef SIM_ANNEAL_H
#define SIM_ANNEAL_H 1

#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <time.h>

#include "controls.h"

class Sim_Anneal
{
private:
  int num_expirations; // the number of tempertures in succession
                       // for which the target number of successful
                       // sample moves was NOT made.
  int num_moves;       // number of sample points successfully moved at current temperature
  int num_attempts;    // number of attempted sample point moves at current temperature
  int iteration;       // number of different temperatures processed, e.g 0,1,2,... 
  double T;            // current temperature, T >= 0
  double magic_number; // 1/Boltzman constant/T
  static double inverse_rand_max; // inverse of the range of pseudorandom integers 
  static double samples_to_random_range_ratio; // ratio of number of spline samples
                                               // to range of pseudorandom integers
  std::vector<double> energy;      // history of total energy in contour
  std::vector<double> temperature; // history of annealing temperature for contour
public:
  Sim_Anneal             (void);
  int getRandomSample    (const int & num_samples);
  bool moreMoves         (const int & num_raw_points);
  void printEnergy       (char * const filename) const;
  void printTemperature  (char * const filename) const;
  void printScript       (const char * const energy_filename,
                          const char * const temp_filename,
                          const std::string & contour_name,
                          const int & contour_section) const;
  void debug             (const std::string & name,
                          const int & section,
                          const int & num_raw_points);
  double getNormalRandom (void);
public:
  /** Pick random path parameter displacement from normal distribution.
   *  \return Path parameter displacement.
   */

  double getRandomPathParameterDisplacement (void)
  {
    double sign = -1.0;
    if ( (rand()*inverse_rand_max) < 0.5) sign *= -1.0; 
    sign *= getNormalRandom();
    return sign;
  }

  /** Lower temperature and prepare for another set of sample moves.
   */

  void decreaseTemp (void)
  {
    T *= Controls::instance().getTempScale();
    magic_number = 1.0/Controls::instance().getBoltzman()/T;
    num_moves    = 0;
    num_attempts = 0;
    iteration++;
  }

  /** Determine if contour is "frozen" and simulated annealing is complete.
   *  \return True if contour is frozen; false otherwise;
   */

  bool isFrozen (void)
  {
    if (num_expirations>2) return true;
    if (T<Controls::instance().getEpsilon()) return true;
    else return false;
  }

  /** Determine if sample point move that raised energy shoudl be kept.
   *  \param[in] delta_energy Change in energy due to sample point move.
   *  \return True if keep sample point move; false otherwise;
   */

  bool keepMove (double delta_energy)
  {
    // Changes to spline sample positions that raise the energy
    // of the sample system are accepted with some probability.
    // The acceptance probability is determined by the temperature
    // of the system (T), the change in energy due to the sample move (delE),
    // and the boltzman constant (kB): p=exp(-delE/kB/T).
    // This function yields acceptance probabilities near zero for very small
    // values of temperature, i.e. algorithm becomes more greedy as temp drops.
    double a = exp(-delta_energy*magic_number);
    return (rand()*inverse_rand_max) < a;
  }
  void setScale (const int & num_samples)
  {
    samples_to_random_range_ratio = (num_samples-1.0)*inverse_rand_max;
  }
  void reserveDetailedInfo (void)
  {
    energy.reserve(10000);
    temperature.reserve(10000);
  }
  void recordDetailedInfo (const double & current_energy)
  {
    energy.push_back(current_energy);
    temperature.push_back(T);
  }
  void recordSuccessfulMove (void)
  {
    num_moves++;
    num_attempts++;
  }
  void recordUnsuccessfulMove (void)
  {
    num_attempts++;
  }
};

#endif

