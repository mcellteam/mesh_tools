#include <string>
#include <cmath>

#include "classes.cc"
#include "functions.cc"

int main(int argc,char *argv[])
{

  Parse p;
  p.parseCommandLine(argc,argv,p.getUsageMessage());

  ////////// declare variables /////////
  char output_script[32] = "mesh_and_convert.csh";
  // spline parameters
  Parameters *pa = new Parameters;
  pa->plot_rad_int=.1;// radius of curvature sampling interval for plotting sampling function
  pa->num=100;		// num is the # samples of the splines between contour points
                        // sequential triplets of sampled points are used to
                        // compute the radius of curvature , e.g. num/2 between contour points
  pa->diag=false;	// set true to print diagnostic files
  pa->max_rad=1E10;	// calculated radius of curvature of spline will saturate at this value
  pa->dmin=.001;	// dmin = minimum spline sampling distance
  pa->dmax=.050;	// dmax = maximum spline sampling distance
  pa->T=1.0;		// sample period (= time to traverse dmax)
  pa->amax=1E-2;	// max radial acceleration

  printf("\nReading input files.\n");
  // create contours
  void_list *ch = getContours(p);

  // remove contours with less than three points
  ch=purgeBadContours(ch);

  // add previous pointers
  for (void_list *q=ch;q!=NULL;q=q->next){((Contour*)q->data)->addPreviousRaw();}

  // check contours for duplicate points
  for (void_list *q=ch;q!=NULL;q=q->next){((Contour*)q->data)->removeDuplicates();}

  // remove contours with less than three points
  ch=purgeBadContours(ch);

  printf("Interpolating contours.\n");
  ///// linearly interpolate raw contour points /////
  if (p.deviation_threshold) {
    for (void_list *q=ch;q!=NULL;q=q->next)
    {
      ((Contour*)q->data)->linearlyInterp(p.deviation_threshold,p.scale);
    }
  }

  ////////// add previous data //////////
  for (void_list *q=ch;q!=NULL;q=q->next){
    ((Contour*)q->data)->addPreviousRaw();
  }

  // NOTE: RAW POINTS IN CONTOURS ARE STORED IN REVERSE ORDER
  // FIRST OFF THE LINKED LIST WAS LAST ADDED TO LIST.

  // create array of number of points in each contour
  int *contour_array = NULL;
  contour_array = createArray(contour_array,ch,countContours(ch));

  printf("Splining contours.\n");
  ////////// fit splines //////////
  Histogram *h = new Histogram();
  fitSplines(ch,h,p.section_thickness,contour_array,p.deviation_threshold,p.scale,p.output_dir.c_str(),pa,countContours(ch));

  ///// compute histogram /////
  computeHistogram(h,ch,contour_array);

  ////////// create objects //////////
  void_list *objectsh=createObjects(ch);

  ////////// clear any existing pts and script files //////////
  clearPtsFiles(p.output_dir.c_str(),objectsh);

  printf("Writing output contours.\n");
  ////////// print each object //////////
  // for each object
  for (void_list *q=objectsh;q!=NULL;q=q->next) {
    Object *o=(Object*)q->data;
    if(!degenerateObject(o) && (p.capping_flag || o->min_section!=o->max_section)){
      ///// print config file /////
      printConfigFile(p.output_dir.c_str(),o,p.capping_flag);
      /////  append to script files /////
      appendScriptFile(p.output_dir.c_str(),o);
      ///// print pts files of interpolated contour points /////
      printPtsFiles(p.output_dir.c_str(),o,p.scale);
      if(p.capping_flag){printCaps(p.output_dir.c_str(),o,p.section_thickness,p.scale);}
    }
  }

  ////////// create calling script //////////
  createCallingScript(p.output_dir.c_str(),output_script);

  ////////// print deviation statistics //////////
  printStatistics(h,p.scale);

  cleanup(objectsh,h,ch,contour_array);
  delete pa;

  return 0;
}
