#include "histogram.h"

#include "contour.h"

Histogram::Histogram(void)
  :min(1e30), max(0), mean(0),stddev(0),sum(0),
  num(0)
{
  int i;
  for(i=0;i<16;i++){
    count[i]=0;
  }
}

void Histogram::update (double const & d)
{
  // update min and max deviation distance
  if(d<min) {min=d;}
  if(d>max) {max=d;}
  num++;
  sum+=d;
}

void Histogram::printStatistics (void)
{
  ////////// print deviation statistics //////////
  printf("\n\nSpline deviation statistics:\n\n");
  printf("  Avg deviation = %g +- %g\n",mean,stddev);
  printf("  Smallest deviation:  %10.5g   |  Largest deviation:  %10.5g\n\n",
         min,max);
  printf("  Deviation histogram:\n");
  printf("  %6.5g - %-6.5g       : %9d    | %6.5g - %-6.5g         : %9d\n",
         max/14.0*0.0,max/14.0*1.0,count[0],
         max/14.0*8.0,max/14.0*9.0,count[9]);
  printf("  %6.5g - %-6.5g       : %9d    | %6.5g - %-6.5g         : %9d\n",
         max/14.0*1.0,max/14.0*2.0,count[1],
         max/14.0*9.0,max/14.0*10.0,count[10]);
  printf("  %6.5g - %-6.5g       : %9d    | %6.5g - %-6.5g         : %9d\n",
         max/14.0*2.0,max/14.0*3.0,count[2],
         max/14.0*10.0,max/14.0*11.0,count[11]);
  printf("  %6.5g - %-6.5g       : %9d    | %6.5g - %-6.5g         : %9d\n",
         max/14.0*3.0,max/14.0*4.0,count[3],
         max/14.0*11.0,max/14.0*12.0,count[12]);
  printf("  %6.5g - %-6.5g       : %9d    | %6.5g - %-6.5g         : %9d\n",
         max/14.0*4.0,max/14.0*5.0,count[4],
         max/14.0*12.0,max/14.0*13.0,count[13]);
  printf("  %6.5g - %-6.5g       : %9d    | %6.5g - %-6.5g         : %9d\n",
         max/14.0*5.0,max/14.0*6.0,count[5],
         max/14.0*13.0,max/14.0*14.0,count[14]);
  printf("  %6.5g - %-6.5g       : %9d    | %6.5g - %-6.5g         : %9d\n",
         max/14.0*6.0,max/14.0*7.0,count[6],
         max/14.0*14.0,max/14.0*15.0,count[15]);
  printf("  %6.5g - %-6.5g       : %9d    | %6.5g -                : %9d\n",
         max/14.0*7.0,max/14.0*8.0,count[7],
         max/14.0*15.0,0);
  printf("\n");
}
