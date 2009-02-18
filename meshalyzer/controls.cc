#include "controls.h"

#include <cassert>
#include <iostream>
#include <sys/stat.h>

#include "container.h"
#include "object.h"

using std::cerr;
using std::cout;
using std::endl;

Controls * Controls::only_one = NULL;

Controls & Controls::instance(void)
{
  // Not thread-safe.
  // -- lock mutex
  if (only_one == NULL)
    only_one = new Controls();
  // -- unlock mutex
  return *only_one;
}

Controls::Controls(void)
  :inpath(),outpath(),folder(false),attr(false),print(false),interf(false),
  vol(false),good_integrity(true),
  FILENAME_SIZE(1024),
  DETECT_POLYGON_EDGE_INTERSECTION(false),
  RAY_EPSILON                     (.00001),
  DOUBLE_EPSILON                  (1E-10),
  SPACE_LENGTH	                  (0.040),
  FACES_PER_BOX	                  (0.1),
  NUM_ADJACENT_BOXES              (1),
  SEARCH_RADIUS                   (0.070),
  NEIGHBORHOOD_RADIUS		  (0.140),
  CLOSEST_POINT_COSINE	          (0.34202),
  PI                              (3.14159265358979),
  NUM_BINS                        (16)
{
  // command line arguments
  signal[0]=signal[1]=signal[2]=signal[3]=signal[4]=false;
  thresholds[0]=thresholds[1]=thresholds[2]=thresholds[3]=thresholds[4]=0.0;
  // cumulative statistics
  bb[0]=bb[1]=bb[2]=bb[3]=bb[4]=bb[5]=0;
  // cumulative Stats
  wbb[0]=wbb[1]=wbb[2]=1E30;
  wbb[3]=wbb[4]=wbb[5]=-1E30;
}

void Controls::parse(int argc,char **argv,std::string message)
{
  // if no arguments passed
  if(argc==1)
  {
    cout << message << endl;
    exit(1);
  }
  int c;
  opterr=0;
  char *eptr=NULL;
  while((c=getopt(argc,argv,"ab:c:d:e:f:hipqv")) != -1)
    switch(c)
    {
      case 'a':
        // determine attributes only
        attr=true;
        break;
      case 'b':
        // specify max edge length threshold (units of input file)
        signal[4]=true;
        thresholds[4]= strtod(optarg,&eptr);
        break;
      case 'c':
        // specify min edge length threshold (units of input file)
        signal[3]=true;
        thresholds[3]= strtod(optarg,&eptr);
        break;
      case 'd':
        // specify max edge angle threshold (degrees)
        signal[2]=true;
        thresholds[2]= strtod(optarg,&eptr);
        break;
      case 'e':
        // specify min edge angle threshold (degrees)
        signal[1]=true;
        thresholds[1]= strtod(optarg,&eptr);
        break;
      case 'f':
        // specify max aspect ratio threshold
        signal[0]=true;
        thresholds[0]= strtod(optarg,&eptr);
        break;
      case 'h':
        cout << message << endl;
        abort ();
      case 'i':
        interf=true;
        break;
      case 'p':
        print=true;
        sprintf(style,"detail");
        break;
      case 'q':
        print=true;
        sprintf(style,"cp");
        break;
      case 'v':
        vol=true;
        break;
      case '?':
        if (optopt == 'c')
          fprintf (stderr,"Option -%c requires an argument.\n", optopt);
        else if (isprint (optopt))
          fprintf (stderr, "Unknown option '-%c.\n", optopt);
        else
          fprintf (stderr,
                   "Unknown option character '\\x%x'.\n",
                   optopt);
        exit(1);
      default:
        cout << message << endl;
        abort ();
    }

  if(interf==true && folder==false)
  {
    cout << "\n\nError. The '-i' option applies only to folders.\n\n";
    exit(1);
  }

  if(optind<argc)
  {
    for (int index=optind;index<argc;index++)
    {
      // determine if file or folder was passed
      struct stat stat_p;
      stat(argv[index],&stat_p);
      if(stat_p.st_mode & S_IFDIR){ folder=true;}
      else{ folder=false;}

      // adjust input directory
      char filename[FILENAME_SIZE];
      strcpy(filename,argv[index]);
      if(folder)
      {
        char *temp=strrchr(filename,'/');
        if(!temp) {strcat(filename,"/");}
        else if(*++temp) {strcat(filename,"/");}
      }
      if(inpath.empty()==true){ inpath=filename; }
      else if(outpath.empty()==true){ outpath=filename; }
    }
  }
  else 
  {
    fprintf (stderr,"No input file argument found on command line.\n");
    exit(1);
  }
}
