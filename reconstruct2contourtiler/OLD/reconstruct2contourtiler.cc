#include <stdlib.h>
#include <string>
#include <cmath>

#include "reconstruct2contourtiler.h"

void createCallingScript (char const * const outdir,char const * script)
{
  char filename[256],line[2048];
  FILE *F;
  sprintf(filename,"%s%s",outdir,script);
  F = fopen(filename,"w");
  if (!F) { printf("Could not open output file %s\n",filename);exit(1);}
  sprintf(line,"#!/bin/bash\n\n/bin/bash mesh.sh\n/bin/bash convert.sh\n");
  fputs(line,F);
  fclose(F);
}

int main (int argc,char *argv[])
{
  Controls & cs (Controls::instance()); 
  cs.parseCommandLine(argc,argv);
  Container c;
  printf("\nReading input files.\n");
  c.getContours();
  c.removeDuplicates();
  c.purgeBadContours();
  if (cs.getLinearThreshold())
  {
    printf("Interpolating contours.\n");
    cout << "Number of raw points: "
          << "begin (" << c.getNumRawPoints() << ")";
    c.interpolateRawPoints();
    cout << ", end (" << c.getNumRawPoints() << ")\n";
    //cout << c.getNumInterpolatedPoints
  }
  else
  {
    cout << "Choice selection: No interpolation of contours.\n";
    cout << "Number of raw points = " << c.getNumRawPoints() << endl;
  }
  //c.countPtsPerContour(); 
  printf("Splining contours.\n");
  Histogram h;
  c.fitSplines(h);
  c.computeHistogram(h);
  c.clearPtsFiles();
  printf("Writing output contours.\n");
  // for each object
  for (o_iterator i = c.firstObjectNonConst();i!=c.onePastLastObjectNonConst();i++)
  {
    //if (!i->degenerateObject() && (cs.getCappingFlag() || i->getMinSection()!=i->getMaxSection()))
    if (!i->degenerateObject())
    {
      ///// print config file /////
      i->printConfigFile();
      /////  append to script files /////
      i->appendScriptFile();
      ///// print pts files of interpolated contour points /////
      i->printPtsFiles();
      if (cs.getPrintRawPoints())
      {
        i->printRawPtsFiles();
      }
      if (cs.getCappingFlag())
      {
        i->printCaps();
      }
    }
    else
    {
      cout << "Warning. Object (" << i->getName() << ") is degenerate.\n";
    }
  }
  createCallingScript(cs.getOutputDir(),cs.getOutputScript());
  if (cs.getDeviationThreshold())
  {
    h.printStatistics();
  }

  return 0;
}
