#include "controls.h"

#include <cassert>
#include <iostream>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

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
  :INPATH(),
  MAX_FILENAME_SIZE(1024),
  DISTANCE_THRESHOLD(1E-1),
  DOUBLE_EPSILON(1E-10),
  RANDOM_VECTOR() 
{
  RANDOM_VECTOR =vector3(0.228869482921038,
                         0.897625868343474,
                         0.376678324659227);
}

std::string Controls::getUsageMessage (void)
{
  std::string message = "\n";
  message=message+
        "NAME\n"+
        "       meshheal - merge border vertices separated by small distances\n"+
        "\nSYNOPSIS\n"+
        "       meshheal [options] FILE\n"+
        "\nDESCRIPTION\n"+
        "       Meshheal finds any holes in the input mesh\n"+
        "       and attempts to close the holes by merging vertices\n"+
        "       along the hole. Closest vertices are moved first.\n"+
        "       Merging stops when no holes exist or the distance between\n"+
        "       vertices is larger than threshold.\n"+
        "\nEXAMPLES\n"+
        "       meshheal filename\n"+
        "\nOPTIONS\n"+
        "       -t NUM\n"+
        "              Free vertices separated by more than NUM will not be merged.\n"+
        "              Default is '" + d2str(DISTANCE_THRESHOLD) + "'.\n\n"+
        "       -h, --help\n"+                     
        "              Print meshheal man page.\n"+
        "\nJustin Kinney				2010/06/01\n";
  return message;
}

void Controls::parse (int argc,char **argv,std::string message)
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
  while((c=getopt(argc,argv,"ht:")) != -1)
    switch(c)
    {
      case 'h':
        cout << message << endl;
        abort ();

      case 't':
        DISTANCE_THRESHOLD = strtod(optarg,&eptr);
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
      INPATH = argv[index];
    }
  }
  else 
  {
    fprintf (stderr,"No input file argument found on command line.\n");
    exit(1);
  }
}

/** Create string from floating point number.
 * \param[in] i Number of interest.
 * \return String containing input number.
 */

std::string Controls::d2str (double const & i)
{
  char file[get_max_filename_size()];
  sprintf(file,"%.15g",i);
  return std::string(file);
}

/** Create string from integer number.
 * \param[in] i Number of interest.
 * \return String containing input number.
 */

std::string Controls::i2str (int const & i)
{
  char file[get_max_filename_size()];
  sprintf(file,"%i",i);
  return std::string(file);
}

