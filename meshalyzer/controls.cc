#include "controls.h"

#include <cassert>
#include <iostream>
#include <stdlib.h>
#include <string.h>
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
  :INPATH(),FOLDER(0),ATTR(0),PRINT(0),DREAMM(0),INTER(0),
  PRINT_SET_VOLUME_ONLY(0),
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
}

std::string Controls::getUsageMessage (void)
{
  std::string message = "\n";
  message=message+
        "NAME\n"+
        "       meshalyzer - mesh quality analyzer\n"+
        "\nSYNOPSIS\n"+
        "       meshalyzer [options] FILE|DIR\n"+
        "\nDESCRIPTION\n"+
        "       Meshalyzer is a general purpose mesh analyzer useful for\n"+
        "       generating a complete summary of the current state of a mesh.\n"+
        "       Meshalyzer assesses mesh integrity (e.g. missing data),\n"+
        "       mesh attributes (e.g. closed, manifold, oriented), and\n"+
        "       mesh characteristics (e.g. number of vertices, faces, edges).\n"+
        "       Batch processing is easy by passing a directory name\n"+
        "       as input on the command line.\n"+
        "\nEXAMPLES\n"+
        "       meshalyzer filename\n"+
        "              Evaluate mesh integrity, attributes,\n"+
        "              and characteristics for the single mesh file.\n\n"+
        "       meshalyzer directoryname\n"+
        "              Evaluate mesh integrity, attributes,\n"+
        "              and characteristics for each single mesh file in directory.\n\n"+
        "       meshalyzer -p filename\n"+
        "              Evaluate mesh integrity, attributes, and characteristics \n"+
        "              for the single mesh file, print the results, and print \n"+
        "              the mesh elements preventing the mesh from being good with\n"+
        "              regards to the mesh characteristics and the attributes, if any.\n\n"+
        "       meshalyzer -a -p filename\n"+
        "              Evaluate the five mesh attributes for the single mesh file,\n"+
        "              print the state of each attribute, and print the mesh elements\n"+
        "              preventing the mesh from being good with regards to the attributes, if any.\n\n"+
        "       meshalyzer -b 10.0 -c 1.0 -p filename\n"+
        "              Evaluate mesh integrity, attributes, and characteristics \n"+
        "              for the single mesh file, print the results, and print \n"+
        "              the mesh elements preventing the mesh from being good with\n"+
        "              regards to the mesh characteristics and the attributes, if any.\n"+
        "              Additionally, screen faces with aspect ratios larger than 10.0 and\n"+
        "              screen edges with lengths larger than 1.0, and print detailed\n"+
        "              information about the offending mesh elements.\n"+
        "\nOPTIONS\n"+
        "       -a\n"+
        "              Evaluate the attributes of the mesh and report the results.\n"+
        "              Skip the evaluation of mesh characteristics.\n\n"+
        "       -b NUM\n"+
        "              Detect edges with length smaller than NUM.\n"+
        "              Units are same as FILE.\n\n"+
        "       -c NUM\n"+
        "              Detect edges with length greater than NUM.\n"+
        "              Units are same as FILE.\n\n"+
        "       -d NUM\n"+
        "              Detect edges with angle between two adjacent faces\n"+
        "              smaller than NUM degrees.\n\n"+
        "       -e NUM\n"+
        "              Detect edges with angle between two adjacent faces\n"+
        "              greater than NUM degrees.\n\n"+
        "       -f NUM\n"+
        "              Detect faces with aspect ratio greater than NUM.\n\n"+
        "       -h\n"+
        "              Print meshalyzer man page.\n\n"+
        "       -i\n"+
        "              Detect intersections between faces from different objects.\n"+
        "              Faceintersection detection is performed once all objects\n"+
        "              are loaded into memory. Single object intersection detection\n"+
        "              is omitted.\n\n"+
        "       -p\n"+
        "              Print detailed information about offending mesh elements \n"+
        "              (i.e. flipped faces, borders, nonmanifold edges,\n"+
        "              nonmanifold vertices, intersecting faces).\n\n"+
        "       -q\n"+
        "              Same as '-p' option, but prints vertex information\n"+
        "              in dreamm custom points format.\n"+
        "       -v\n"+
        "              If folder passed as argument, then only print total\n"+
        "              set volume and nothing else.\n"+
        "\nJustin Kinney				2007/10/01\n";
  return message;
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
        ATTR=1;
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
        INTER=1;
        break;
      case 'p':
        PRINT=1;
        //style = "detail";
        DREAMM = 0;
        break;
      case 'q':
        PRINT=1;
        //style = "cp";
        DREAMM = 1;
        break;
      case 'v':
        PRINT_SET_VOLUME_ONLY=1;
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

  if(optind<argc)
  {
    for (int index=optind;index<argc;index++)
    {
      if(INPATH.empty()==false)
      {
        cout << "Error. Multiple input data paths found on command line.\n";
        exit(1); 
      }
      // determine if file or folder was passed
      struct stat stat_p;
      stat(argv[index],&stat_p);
      if(stat_p.st_mode & S_IFDIR)
      {
        // adjust input directory
        char filename[FILENAME_SIZE];
        strcpy(filename,argv[index]);
        char *temp=strrchr(filename,'/');
        if(!temp) {strcat(filename,"/");}
        else if(*++temp) {strcat(filename,"/");}
        INPATH = filename;
        FOLDER = 1;
      }
      else
      {
        INPATH = argv[index];
        FOLDER = 0;
      }

    }
  }
  else 
  {
    fprintf (stderr,"No input file argument found on command line.\n");
    exit(1);
  }

  if(INTER==1 && FOLDER==0)
  {
    cout << "\n\nError. The '-i' option applies only to folders.\n\n";
    exit(1);
  }

}
