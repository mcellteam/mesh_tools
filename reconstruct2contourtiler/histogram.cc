#include<iostream>

#include "histogram.h"

#include "contour.h"

using std::cout;
using std::endl;
using std::left;
using std::right;

extern const int num_bins;

Histogram::Histogram(void)
  :N(0),min(1e30),max(-1e30),sum(0.0),sum2(0.0),
  count(NULL),bins(NULL),x()
{
  count = new int [num_bins];
  bins = new double[num_bins+1];
  for(int i=0;i<num_bins;i++)
  {
    count[i]=0;
    bins[i]=0.0;
  }
  bins[num_bins] = 0.0;
}

Histogram::~Histogram (void)
{
  delete[] count;
  delete[] bins;
}

/** Record statistics of straight-line distances between sample points.
 * \param[in] my_path_intervals Inter-sample distances.
 */

void Histogram::scan (const vec_d & my_path_intervals)
{
  for (c_d_iterator i=my_path_intervals.begin();i!=my_path_intervals.end();i++)
  {
    update(*i);
  }
}

/** Update statistics of data.
 * \param[in] d New data point.
 */

void Histogram::update (double const & d)
{
  // update min and max deviation distance
  if(d<min) {min=d;}
  if(d>max) {max=d;}
  N++;
  sum+=d;
  sum2+=d*d;
  x.push_back(d);
}

/** Compute mean value.
 *  \return Mean.
 */

double Histogram::mean (void) const
{
  return sum/N;
}

/** Compute variance value.
 *  \return Variance.
 */

double Histogram::variance (void) const
{
  if (N>1) return (sum2-(sum*sum/N))/(N-1);
  else return 0;
}

/** Print summary statistics.
 */

void Histogram::printStats (void) const
{
  cout << "       N        " << N << endl;
  cout << "       min      " << min << endl;
  cout << "       max      " << max << endl;
  cout << "       mean     " << mean() << endl;
  cout << "       variance " << variance() << endl << endl;
}

/** Print histogram.
 */

void Histogram::printHistogram (void)
{
  printf("  %9.4g - %-9.4g    : %9d  | %9.4g - %-9.4g   : %9d\n",
         bins[0], bins[1], count[0], bins[8], bins[9], count[8]);
  printf("  %9.4g - %-9.4g    : %9d  | %9.4g - %-9.4g   : %9d\n",
         bins[1], bins[2], count[1], bins[9], bins[10], count[9]);
  printf("  %9.4g - %-9.4g    : %9d  | %9.4g - %-9.4g   : %9d\n",
         bins[2], bins[3], count[2], bins[10], bins[11], count[10]);
  printf("  %9.4g - %-9.4g    : %9d  | %9.4g - %-9.4g   : %9d\n",
         bins[3], bins[4], count[3], bins[11], bins[12], count[11]);
  printf("  %9.4g - %-9.4g    : %9d  | %9.4g - %-9.4g   : %9d\n",
         bins[4], bins[5], count[4], bins[12], bins[13], count[12]);
  printf("  %9.4g - %-9.4g    : %9d  | %9.4g - %-9.4g   : %9d\n",
         bins[5], bins[6], count[5], bins[13], bins[14], count[13]);
  printf("  %9.4g - %-9.4g    : %9d  | %9.4g - %-9.4g   : %9d\n",
         bins[6], bins[7], count[6], bins[14], bins[15], count[14]);
  printf("  %9.4g - %-9.4g    : %9d  | %9.4g - %-9.4g   : %9d\n",
         bins[7], bins[8], count[7], bins[15], bins[16], count[15]);
}

/** Calculate bin partitions and then bin data.
 */

void Histogram::createHistogram (void)
{
  double bin_width=(4.0*sqrt(variance()))/(static_cast<double>(num_bins)-1.0);

  bool nonneg=false;

  double mm = mean();

  while (nonneg==false)
  {
    // create bins
    bins[0]=0;
    bins[1]=0;

    for (int i=2;i<num_bins;i++)
    {
      bins[i]=mm+(((i-1)-(static_cast<double>(num_bins)/2.0))*bin_width);	
    }
    bins[num_bins]=max;

    nonneg=true;
    //check
    for (int i=2;i<num_bins+1;i++)
    {
      if (bins[i]<0)
      {
        nonneg=false;
      }
    }
    if (nonneg==false)
    {
      // increase mm
      mm = mm*1.1;
    }
  }

  // bin data
  unsigned int start,mid,end;

  // for each element in x
  for (d_iterator i=x.begin();i!=x.end();i++)
  {
    start=0;
    end=num_bins-1;
    while (end-start>1)
    {
      mid=(start+end)/2;
      if (*i >= bins[start] && *i <= bins[mid]) { end=mid;}
      else { start=mid; }
    }
    if (*i >= bins[start] && *i <= bins[end]) { count[start]++;}
    else { count[end]++; }
  }
}

void Histogram::printStatistics (std::string label)
{
  printf("\nSpline %s statistics--------------------------------\n\n",label.c_str());
  printStats();
  createHistogram();
  printHistogram();
}
