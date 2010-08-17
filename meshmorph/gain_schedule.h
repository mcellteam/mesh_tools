// Author: Justin Kinney
// Date: Sep 2008

#ifndef GAIN_SCHEDULE_H
#define GAIN_SCHEDULE_H 1

#include "meshmorph.h"

class Gain_Schedule
{
private:
  float  max_gain; // maximum allowed gain
//  float  ref_gain; // reference gain
  float      gain; // active gain, all the time
  float      step; // gain change increment
//  int      period; // used to schedule gain changes
  static Gain_Schedule * only_one;
  Gain_Schedule                   (void);
  Gain_Schedule                   (Gain_Schedule const &);
  Gain_Schedule & operator =      (Gain_Schedule const &);
  ~Gain_Schedule                  (void);
public:
  static Gain_Schedule & instance (void);
  void   updateMaxGain            (void);
  void   parseGainFile (const char * filename);
//  void   updateGain               (void);

  /** Get max gain.
   * \return Maximum gain.
   */

  double getMaxGain (void) const
  {
    return max_gain;
  }

  /** Get gain.
   * \return Gain.
   */

  double getGain (void) const
  {
    return gain;
  }

//  /** Get reference gain.
//   * \return Reference gain.
//   */
//
//  double getRefGain (void) const
//  {
//    return ref_gain;
//  }

//  /** Reset gain to reference value.
  /** Reset gain to max value.
  */

  void resetGain (void)
  {
    //gain = ref_gain;
    gain = max_gain;
  }

  /** Reduce gain by 50%.
  */

  void halveGain (void)
  {
    gain = gain/2.0;
  }

//  /** Initialize reference gain.
//  */
//
//  void initGain (void)
//  {
//    // initialize reference gain
//    ref_gain = gain;
//  }

//  /** Update decrement time until next gain change.
//  */
//
//  void updateGainPeriod (void)
//  {
//    // update period
//    if (period>0){period--;}
//  }

};

#endif
