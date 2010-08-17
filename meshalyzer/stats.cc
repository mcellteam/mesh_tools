#include "stats.h"

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <iostream>

using std::cout;
using std::endl;

Stats::Stats (void)
  :total(0),min(DBL_MAX),max(-DBL_MAX),sum(0),sum2(0),count(NULL),bins(NULL),x()
{
  int a = Controls::instance().get_num_bins();
  count = new int [a];
  bins  = new double [a+1];
  for (int i=0;i<a;i++)
  {
    count[i]=0;
    bins[i]=0;
  }
  bins[a]=0;
}

Stats::~Stats (void)
{
  delete[] count;
  delete[] bins;
}

void Stats::clear (void)
{
  int a = Controls::instance().get_num_bins();
  sum=sum2=0;
  min=DBL_MAX;
  max=-DBL_MAX;
  for (int i=0;i<a;i++)
  {
    count[i]=0;
    bins[i]=0;
  }
  bins[a]=0;
  x.clear();
}

void Stats::printHistogram (void)
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

void Stats::printAdjacentFaceHistogram (void)
{
  int L[16],U[16];
  for (int i=0;i<16;i++)
  {
    L[i]=1+static_cast<int>(bins[i]);
    U[i]=static_cast<int>(bins[i+1]);
  }
  printf("  %9d - %-9d    : %9d  | %9d - %-9d   : %9d\n",
         0, U[0], count[0], L[8], U[8], count[8]);
  printf("  %9d - %-9d    : %9d  | %9d - %-9d   : %9d\n",
         L[1], U[1], count[1], L[9], U[9], count[9]);
  printf("  %9d - %-9d    : %9d  | %9d - %-9d   : %9d\n",
         L[2], U[2], count[2], L[10], U[10], count[10]);
  printf("  %9d - %-9d    : %9d  | %9d - %-9d   : %9d\n",
         L[3], U[3], count[3], L[11], U[11], count[11]);
  printf("  %9d - %-9d    : %9d  | %9d - %-9d   : %9d\n",
         L[4], U[4], count[4], L[12], U[12], count[12]);
  printf("  %9d - %-9d    : %9d  | %9d - %-9d   : %9d\n",
         L[5], U[5], count[5], L[13], U[13], count[13]);
  printf("  %9d - %-9d    : %9d  | %9d - %-9d   : %9d\n",
         L[6], U[6], count[6], L[14], U[14], count[14]);
  printf("  %9d - %-9d    : %9d  | %9d - %-9d   : %9d\n",
         L[7], U[7], count[7], L[15], U[15], count[15]);
}

void Stats::printAspectRatioHistogram (void)
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
  printf("  %9.4g - %-9.4g    : %9d  | %9.4g - 10000       : %9d\n",
         bins[5], bins[6], count[5], bins[13], count[13]);
  printf("  %9.4g - %-9.4g    : %9d  |     10000 - 100000      : %9d\n",
         bins[6], bins[7], count[6], count[14]);
  printf("  %9.4g - %-9.4g    : %9d  |    100000 -             : %9d\n",
         bins[7], bins[8], count[7], count[15]);
}

void Stats::printStats (void)
{
  cout << "       min      " << min << endl;
  cout << "       max      " << max << endl;
  cout << "       median   " << median() << endl;
  cout << "       mean     " << mean() << endl;
  cout << "       variance " << variance() << endl << endl;
}

double Stats::median (void)
{
  sort(x.begin(),x.end());
  int i = x.size()/2;
  //	cout << "\nStats::median: "
  //	<< "x.size() = " << x.size()
  //	<< ", i = " << i << endl;
  return x[i];
}

void Stats::createHistogram (void)
{
  int a = Controls::instance().get_num_bins();
  double bin_width=(4.0*sqrt(variance()))/(a-1.0);

  bool nonneg=false;

  double mm = mean();

  while (nonneg==false)
  {
    // create bins
    bins[0]=0;
    bins[1]=0;

    for (int i=2;i<a;i++)
    {
      bins[i]=mm+(((i-1)-(static_cast<double>(a)/2.0))*bin_width);	
    }
    bins[Controls::instance().get_num_bins()]=max;

    nonneg=true;
    //check
    for (int i=2;i<a+1;i++)
    {
      if (bins[i]<0)
      {
        nonneg=false;
        //				cout << "Stats::createHistogram: " << bins[i] << " < 0\n";
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
    end=a-1;
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

void Stats::createAdjacentFaceHistogram (void)
{
  // create bins
  bins[0]=0;
  bins[1]=0.99;
  bins[2]=1.99;
  bins[3]=2.99;
  bins[4]=3.99;
  bins[5]=4.99;
  bins[6]=5.99;
  bins[7]=6.99;
  bins[8]=7.99;
  bins[9]=8.99;
  bins[10]=9.99;
  bins[11]=10.99;
  bins[12]=11.99;
  bins[13]=12.99;
  bins[14]=13.99;
  bins[15]=14.99;
  if (max>15.99){bins[16]=max;}
  else{bins[16]=15.99;}
  // bin data
  unsigned int start,mid,end;

  // for each element in x
  int a = Controls::instance().get_num_bins();
  for (d_iterator i=x.begin();i!=x.end();i++)
  {
    start=0;
    end=a-1;
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

void Stats::createAspectRatioHistogram (void)
{
  // create bins
  bins[0]=1.1547;
  bins[1]=1.5;
  bins[2]=2;
  bins[3]=2.5;
  bins[4]=3;
  bins[5]=4;
  bins[6]=6;
  bins[7]=10;
  bins[8]=15;
  bins[9]=25;
  bins[10]=50;
  bins[11]=100;
  bins[12]=300;
  bins[13]=1000;
  bins[14]=10000;
  bins[15]=100000;

  // bin data
  unsigned int start,mid,end;

  // for each element in x
  int a = Controls::instance().get_num_bins();
  for (d_iterator i=x.begin();i!=x.end();i++)
  {
    start=0;
    end=a-1;
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

double Stats::mean (void)
{
  return sum/x.size();
}

double Stats::variance (void)
{
  if (x.size()>1)
  {
    return (sum2-(sum*sum/x.size()))/(x.size()-1);
  } else
  {
    //		cout << "variance not computed since sample size==1.\n";
    //		return -666;
    return 0;
  }
}

