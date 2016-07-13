#include <stdlib.h>
#include <string>
#include <cmath>

#include "reconstruct_interpolate.h"

int main (int argc,char *argv[])
{
  Controls & cs (Controls::instance()); 
  cs.parseCommandLine(argc,argv);
  Container c;
  printf("\nReading input files...\n");
  c.getContours();
  // compare sequential points in contours and remove ducplicate points
  c.removeDuplicates();
  printf("\nReading input files complete.\n");
  // remove contours with less than 
  // specified number of points
  cout << endl;
  cout << "Number of input files = " << c.getNumFiles()     << endl;
  cout << "Number of objects     = " << c.getNumObjects()   << endl;
  cout << "Number of contours    = " << c.getNumContours()  << endl;
  cout << "Number of raw points  = " << c.getNumRawPoints() << endl;
  cout << endl;
  if (cs.getPrintDetailedInfo()) c.printRawPoints();
  Histogram h; // sample deviations from linearly-interpolated raw points
  Histogram si; // linear distance between spline samples AFTER annealing
  Histogram si_before; // linear distance between spline samples BEFORE annealing
  c.processContour (h,si,si_before);
  if (!strcmp(cs.getOutputSerPrefix(),"")) {
    c.clearOutputScripts ();
  }
  printf ("Writing output contours.\n");
  c.writeOutputContours ();
  h.printStatistics ("deviation");
  si_before.printStatistics ("sample interval before filtering");
  si.printStatistics ("sample interval after filtering");
  cout << endl;
  cout << endl;
}

/** Determine if two floating-point precision numbers
 * are equivalent in value within epsilon.
 * \param[in] a First number.
 * \param[in] b Second number.
 * \param[in] epsilon The difference between the two input values must be
 * greater than the fraction of the largest input value defined by epsilon.
 * \return 1 if Inputs are different; 0 otherwise.
 */

bool distinguishable (double a,double b,double epsilon)
{
  double c;
  c=a-b;
  if (c<0) c=-c;
  if (a<0) a=-a;
  if (a<1) a=1;
  if (b<0) b=-b;
  if (b<a) return (c>a*epsilon);
  else return (c>b*epsilon);
}

/** Determine if two floating-point precision numbers
 * are equivalent in value within MY_DOUBLE_EPSILON. 
 * \param[in] a First number.
 * \param[in] b Second number.
 * \return 1 if Inputs are different; 0 otherwise.
 */

bool distinguishable (double a,double b)
{
  return distinguishable(a, b, Controls::instance().getEpsilon());
}

