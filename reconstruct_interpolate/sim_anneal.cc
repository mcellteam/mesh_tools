#include "sim_anneal.h"

double Sim_Anneal::inverse_rand_max;
double Sim_Anneal::samples_to_random_range_ratio;

Sim_Anneal::Sim_Anneal (void)
  :num_expirations(0),
  num_moves(0),
  num_attempts(0),
  iteration(0),
  T(Controls::instance().getHighTemp()),
  magic_number(1.0/Controls::instance().getBoltzman()/T),
  energy(),temperature()
{
//    srand ( time(NULL) );
    srand ( 1 );
    inverse_rand_max = 1.0/static_cast<double>(RAND_MAX);
}

/** Pick random sample point from uniform distribution.
 *  \param[in] num_samples Number of sample points in contour.
 *  \return Index of sample point.
 */

int Sim_Anneal::getRandomSample (const int & num_samples)
{
  int sample_index = floor(rand()*samples_to_random_range_ratio);
  assert(sample_index>=0);
  assert(sample_index<num_samples);
  return sample_index;
}

/** Get normally distributed random number with equal mean
  *  and standard deviation.
 *  \return Random number.
 */

double Sim_Anneal::getNormalRandom (void)
{
  /* Implements the Polar form of the Box-Muller Transformation
     (c) Copyright 1994, Everett F. Carter Jr.
     Permission is granted by the author to use
     this software for any application provided this
     copyright notice is preserved.
   */
  double x1, x2, w, y1;
  static double y2;
  static int use_last = 0;
  const double mean = Controls::instance().getMeanAmplitudeNoise();
  const double standard_deviation = mean;

  if (use_last) /* use value from previous call */
  {
    y1 = y2;
    use_last = 0;
  }
  else
  {
    do {
      x1 = 2.0 * rand() * inverse_rand_max - 1.0;
      x2 = 2.0 * rand() * inverse_rand_max - 1.0;
      w = x1 * x1 + x2 * x2;
    } while ( w >= 1.0 );

    w = sqrt( (-2.0 * log( w ) ) / w );
    y1 = x1 * w;
    y2 = x2 * w;
    use_last = 1;
  }
  // normal random variate generator
  return mean+y1*standard_deviation;
}

/** Print detailed information for debugging purposes.
 *  \param[in] name Contour name.
 *  \param[in] section Contour section.
 *  \param[in] num_raw_points Number of raw points in contour.
 */

void Sim_Anneal::debug (const std::string & name,
                        const int & section,
                        const int & num_raw_points)
{
  std::cout << "Sim_Annel: "
        << "name = " << name
        << ", section = " << section 
        << ", num_moves " << num_moves
        << ", target_num_successful_moves "
        << Controls::instance().getNumMovesPerIteration()*num_raw_points
        << ", num_attempts = " << num_attempts
        << ", max_num_move_attempts = "
        << Controls::instance().getMaxNumMoveAttempts()*num_raw_points
        << ", num_expirations = " << num_expirations
        << ", T = " << T
        << std::endl;
}

/** Determine if more sample moves should be made at current temperature.
 *  \param[in] num_raw_points Number of raw points in contour.
 *  \return True if move more sampe points; false otherwise;
 */

bool Sim_Anneal::moreMoves (const int & num_raw_points)
{
  Controls & cs (Controls::instance()); 
  bool success = num_moves>=cs.getNumMovesPerIteration()*num_raw_points;
  bool give_up = num_attempts>cs.getMaxNumMoveAttempts()*num_raw_points;
  if (success)
  {
    num_expirations=0;
    return false;
  }
  else if (give_up)
  {
    num_expirations++;
    return false;
  }
  else
  {
    return true;
  }
}

/** Write energy history of simulated annealing to file.
 * \param[in] filename Output file name.
 */

void Sim_Anneal::printEnergy (char * const filename) const
{
  FILE *F = fopen(filename,"w");
  if (!F) { printf("Could not open output file %s\n",filename); return; }
  char line[2048];
  int count = static_cast<int>(energy.size());
  sprintf(line,"# num points = %d\n",count);
  fputs(line,F);
  // for each point
  for (int i=0;i<count;i++)
  {
    sprintf(line,"%d %g\n",i,energy[i]);
    fputs(line,F);
  }
  fclose(F);
}

/** Write temperature history of simulated annealing to file.
 * \param[in] filename Output file name.
 */

void Sim_Anneal::printTemperature (char * const filename) const
{
  FILE *F = fopen(filename,"w");
  if (!F) { printf("Could not open output file %s\n",filename); return; }
  char line[2048];
  int count = static_cast<int>(energy.size());
  sprintf(line,"# num points = %d\n",count);
  fputs(line,F);
  // for each point
  for (int i=0;i<count;i++)
  {
    sprintf(line,"%d %g\n",i,temperature[i]);
    fputs(line,F);
  }
  fclose(F);
}

/** Write grace script for visualizing annealing energy and temperature.
 * \param[in] energy_filename Output file name for energy history.
 * \param[in] temp_filename Output file name for temperature history.
 * \param[in] contour_name Contour name.
 * \param[in] contour_section Contour section.
 */

void Sim_Anneal::printScript (const char * const energy_filename,
                              const char * const temp_filename,
                              const std::string & contour_name,
                              const int & contour_section) const
{
  Controls & cs (Controls::instance()); 
  char filename[256],line[2048];
  FILE *F;
  // print grace script
  sprintf(filename,"%sbfile_anneal_%s_%i",cs.getOutputDir(),contour_name.c_str(),contour_section);
  F = fopen(filename,"w");
  if (!F) { printf("Could not open output file %s\n",filename); return; }
  sprintf(line,"#Obligatory descriptive comment\n");
  fputs(line,F);
  sprintf(line,"arrange ( 1, 2, .11, .3,.5)\n");
  fputs(line,F);
  sprintf(line,"FOCUS G0\n");
  fputs(line,F);
  sprintf(line,"VIEW XMIN 0.15\nVIEW XMAX 1.15\nVIEW YMIN 0.15\nVIEW YMAX 0.85\n");
  fputs(line,F);
  sprintf(line,"READ xy \"%s\"\n",energy_filename);
  fputs(line,F);
  sprintf(line,"s0 line type 1\ns0 line color 1\ns0 errorbar color 1\n");
  fputs(line,F);
  sprintf(line,"yaxis label \"Energy\"\n");
  fputs(line,F);
  sprintf(line,"FOCUS G1\n");
  fputs(line,F);
  sprintf(line,"VIEW XMIN 0.15\nVIEW XMAX 1.15\nVIEW YMIN 0.15\nVIEW YMAX 0.85\n");
  fputs(line,F);
  sprintf(line,"READ xy \"%s\"\n",temp_filename);
  fputs(line,F);
  sprintf(line,"s0 line type 1\ns0 line color 2\ns0 errorbar color 2\n");
  fputs(line,F);
  sprintf(line,"yaxis label \"Temperature\"\n");
  fputs(line,F);
  sprintf(line,"xaxis  off\nyaxis  label place opposite\nyaxis  ticklabel place opposite\n");
  fputs(line,F);
  fclose(F);
  // print shell script
  sprintf(filename,"%s%s_%i_anneal.sh",cs.getOutputDir(),contour_name.c_str(),contour_section);
  F = fopen(filename,"w");
  if (!F) { printf("Could not open output file %s\n",filename); return; }
  sprintf(line,"#!/bin/bash\n");
  fputs(line,F);
  sprintf(line,"xmgrace -noask -nosafe -batch bfile_anneal_%s_%i\n",contour_name.c_str(),contour_section);
  fputs(line,F);
  fclose(F);
}
