// Author: Justin Kinney
// Date: Sep 2008

#include "gain_schedule.h"

#include <iostream>

#include "controls.h"

using std::cout;
using std::cerr;
using std::endl;
using std::left;

Gain_Schedule * Gain_Schedule::only_one = NULL;

Gain_Schedule & Gain_Schedule::instance(void)
{
  // Not thread-safe.
  // -- lock mutex
  if (only_one == NULL)
    only_one = new Gain_Schedule();
  // -- unlock mutex
  return *only_one;
}

Gain_Schedule::Gain_Schedule (void)
  :max_gain(0),ref_gain(0),gain(0),step(0),period(0)
{
  Controls & cs(Controls::instance());
  // inital max gain
  max_gain = cs.get_overall_gain();
  // initialize gain
  gain = max_gain;
  // gain step size
  step = max_gain/cs.get_num_groups();
}

/** Update the maximum allowed value of gain and enforce limit.
 */

void Gain_Schedule::updateMaxGain (void)
{
  if (Controls::instance().get_disable_gain_scheduling()==false)
  {
    if (max_gain>step) max_gain -= step;
    if (gain>max_gain) gain = max_gain;
  }
}

/** Update the gain used to convert force to displacement.
 */

void Gain_Schedule::updateGain (void)
{
  // if gain changeable?
  if (period==0)
  {
    // if change in mean energy < epsilon
    // and gain is larger than gain step size
    //if (Energy::instance().atSteadyState()==true && gain>step)
    //{
    //  // decrease gain
    //  gain -= step;
    //  // reset gain period
    //  period = REFRACTORY_PERIOD;
    //  // update reference gain
    //  ref_gain = gain;
    //} 
  }
}

