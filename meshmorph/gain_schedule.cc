// Author: Justin Kinney
// Date: Sep 2008

#include "gain_schedule.h"

#include <cassert>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <string.h>

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
//  :max_gain(0),ref_gain(0),gain(0),step(0),period(0)
  :max_gain(0),gain(0),step(0)
{
  Controls & cs(Controls::instance());
  // inital max gain
  max_gain = cs.get_overall_gain();
  // initialize gain
  gain = max_gain;
  // gain step size
  step = cs.get_gain_step()*max_gain/cs.get_num_groups();
}

void Gain_Schedule::parseGainFile (const char *filename)
{
  // open file
  std::ifstream f(filename);
  if (f.is_open()==false) { printf("Couldn't open input file %s\n",filename); return;}
  // following line is optional and directs meshmorph
  // to execute the following line exactly once.
  // ONCE max_gain NUM
  // step NUM
  std::string name, value;
  vec_s new_gain_schedule;
  f >> name;
  while(!f.eof())
  {
    //cout << "Gain_Schedule::parseGainFile: "
    //      << "name = " << name << std::endl;
    //cout.flush();
    if (strcmp(name.c_str(),"#")==false)
    {
      f >> name;
      f >> value;
      std::string str = "# " + name + " " + value;
      new_gain_schedule.push_back(str);
      f >> name;
      continue;
    }
    f >> value;
    //cout << "Gain_Schedule::parseGainFile: "
    //      << "value = " << value << std::endl;
    //cout.flush();
    std::string str = "# " + name + " " + value;
    new_gain_schedule.push_back(str);
    if (strcmp(name.c_str(),"max_gain")==false)
    {
      double a = strtod(value.c_str(),NULL);
      cout << "\n\nGain_Schedule::parseGainFile: "
            << "Update max_gain to "
            << a << " (was "
            << max_gain << ").\n\n";
      cout.flush();
      max_gain = a;
    }
    if (strcmp(name.c_str(),"step")==false)
    {
      double b = strtod(value.c_str(),NULL);
      cout << "\n\nGain_Schedule::parseGainFile: "
            << "Update step to "
            << b << " (was "
            << step << ").\n\n";
      cout.flush();
      step = b;
    }
    f >> name;
  }
  // close file
  f.close();
  if (new_gain_schedule.empty()==false)
  {
    std::ofstream newfile;
    newfile.open (filename,std::ios::trunc);
    if (newfile.is_open()==false)
    {
      cout << "\nGain_Schedule::parseGainFile: "
            << "Error. Unable to open file = "
            << newfile << endl;
      assert(newfile.is_open()==true);
      exit(1);
    }
    for (s_cit i=new_gain_schedule.begin();i!=new_gain_schedule.end();++i)
    {
      newfile << *i << std::endl;
    }
    newfile.close();
  }
}

/** Update the maximum allowed value of gain and enforce limit.
 */

void Gain_Schedule::updateMaxGain (void)
{
  if (Controls::instance().get_disable_gain_scheduling()==false)
  {
    // DEBUG
    //cout << "\nGain_Schedule::updateMaxGain: "
    //      << "this ran.\n";
    //cout << "\nGain_Schedule::updateMaxGain: "
    //      << "get_gain_scheduling_file().c_str(),"") = "
    //      << strcmp(Controls::instance().get_gain_scheduling_file().c_str(),"")
    //      << "strcmp(get_gain_scheduling_file().c_str(),"") = "
    //      << strcmp(Controls::instance().get_gain_scheduling_file().c_str(),"")
    //      << std::endl;
    //cout.flush();
    // DEBUG
    max_gain -= step;
    gain = max_gain;
    // if gain schedule file specified by user
    if (strcmp(Controls::instance().get_gain_scheduling_file().c_str(),"")!=0)
    {
      // DEBUG
      //cout << "\nGain_Schedule::updateMaxGain: "
      //      << "this ran2.\n";
      //cout.flush();
      // DEBUG
      // check file and update max_gain and step 
      parseGainFile(Controls::instance().get_gain_scheduling_file().c_str());
    }
  }
}

///** Update the gain used to convert force to displacement.
// */
//
//void Gain_Schedule::updateGain (void)
//{
//  // if gain changeable?
//  if (period==0)
//  {
//    // if change in mean energy < epsilon
//    // and gain is larger than gain step size
//    //if (Energy::instance().atSteadyState()==true && gain>step)
//    //{
//    //  // decrease gain
//    //  gain -= step;
//    //  // reset gain period
//    //  period = REFRACTORY_PERIOD;
//    //  // update reference gain
//    //  ref_gain = gain;
//    //} 
//  }
//}

