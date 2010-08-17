#include <stdlib.h>

#include "controls.h"

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

Controls::Controls (void)
  :SECTION_THICKNESS       (.05),
  SCALE                    (1),
  DEVIATION_THRESHOLD      (0),
  LINEAR_THRESHOLD         (0),
  CAPPING_FLAG             (1),
  MIN_SECTION              (60),
  MAX_SECTION              (160),
  INPUT_DIR                ("./"),
  OUTPUT_DIR               ("./"),
  PREFIX                   ("Volumejosef"),
  IGNORED_CONTOURS         (),
  OUTPUT_SCRIPT            ("mesh_and_convert.sh"),
  SPLINE_SAMPLES_PER_POINT (10),
  diag                     (false),
  max_rad                  (1E10),
  dmin                     (1e-6),
  dmax                     (5e-5),
  T                        (1.0),
  amax                     (1E-5),
  OBJECT_RESERVE_SIZE      (1000),
  MAX_CONTOURS_PER_OBJECT  (1000),
  MAX_DEVIATION_ADJUSTMENTS (10),
  PRINT_RAW_POINTS (0)
{}

std::string Controls::getUsageMessage (void)
{
  std::string message = "\n";
  message=message+
        "NAME\n"+
        "       reconstruct2contourtiler - generate contour_tiler input files from contours\n"+
        "\nSYNOPSIS\n"+
        "       reconstruct2contourtiler [options]\n"+
        "\nDESCRIPTION\n"+
        "       Converts reconstruct contour format to contour_tiler input format.\n"+
        "       All files in input directory are assumed to be\n"+
        "       of the form filename_prefix.section#.\n"+
        "       min_section to max_section is the section range to be converted.\n"+
        "       Section_thickness should be in same scale as x,y contour points.\n"+
        "       x,y and z coordinates of sampled splines will be multipled by scale in output.\n"+
        "       CAPPING_FLAG=1 to attempt end capping.CAPPING_FLAG=0 to leave ends open.\n"+
        "       DEVIATION_THRESHOLD is the maximum allowed deviation of the spline from raw\n"+
        "       contour points in scaled units. Set\n"+
        "       DEVIATION_THRESHOLD to 0 to disable thresholding.\n"+
        "\nEXAMPLES\n"+
        "       reconstruct2contourtiler -i ./contours -f myContours -n 10 -x 100 -t .07 -s 1000 -o ./contour_tiler_output -d 2\n"+
        "              Read contours from directory './contours' and write contour_tiler output\n"+
        "              files to directory 'contour_tiler_output'. The input contour files have the\n"+
        "              name scheme 'myContours.#' where # is the section number which varies from\n"+
        "              10 to 100. The distance between contours in the direction of sectioning\n"+
        "              is .070 microns. The contour_tiler output data will be in nanometers as\n"+
        "              dictated by the 1000 scaling. Capping directives will be included in the\n"+
        "              output data. The interpolated contours will not deviate from the input contours\n"+
        "              by more than 2 (nanometers since scaling is 1000).\n"+
        "\nOPTIONS\n"+
        "       --no_capping\n"+
        "              Do not include directives in the output data to cap the meshes,\n"+
        "              thus creating open surface meshes.\n"+ 
        "              Default is to cap the meshes at the minimum and maximum sections\n"+
        "              thereby creating closed surface meshes.\n\n"+
        "       --print_input_as_pts\n"+
        "              Print input points (possibly linearly interpolated) as .pts\n"+
        "              files in output directory.\n"+ 
        "              Default is to not print input points.\n\n"+
        "       -n NUM, --min_section=NUM\n"+               
        "              The starting section number in the section range.\n"+
        "              Default is '" + i2str(MIN_SECTION) + "'.\n\n"+
        "       -x NUM, --max_section=NUM\n"+        
        "              The ending section number in the section range.\n"+
        "              Default is '" + i2str(MAX_SECTION) + "'.\n\n"+
        "       -t NUM, --section_thickness=NUM\n"+              
        "              Thickness of the section in microns. Each section\n"+
        "              is assumed to be of identical thickness.\n"+
        "              Default is '" + d2str(SECTION_THICKNESS) + "'.\n\n"+
        "       -S NUM, --spline_samples_per_point=NUM\n"+              
        "              Sample each spline NUM times between each\n"+
        "              pair of input points.\n"+
        "              Default is '" + d2str(SPLINE_SAMPLES_PER_POINT) + "'.\n\n"+
        "       -s NUM, --scale=NUM\n"+ 
        "              The output data will be scaled by NUM. If NUM is\n"+
        "              equal to 1 then the output data is in microns.\n"+
        "              Default is '" + d2str(SCALE) + "' so the output\n"+
        "              data is in nanometers.\n\n"+
        "       -d NUM, --deviation_threshold=NUM\n"+ 
        "              The input contours are filtered before output.\n"+
        "              The deviation of the input and output contours\n"+
        "              is constrained to be less than NUM where units\n"+
        "              are input units scaled by --scale value.\n"+
        "              Default is '" + d2str(DEVIATION_THRESHOLD) + "'.\n\n"+
        "       -m NUM, --linear_threshold=NUM\n"+ 
        "              Sequential points on the input contours are constrained\n"+
        "              to have separation distance less than NUM. New points\n"+
        "              are inserted by linear interpolation if required.\n"+
        "              Units are input units scaled by --scale value.\n"+
        "              Default is '" + d2str(LINEAR_THRESHOLD) + "' which means\n"+
        "              no threshold is enforced.\n\n"+
        "       -i DIRECTORY, --input_data_dir=DIRECTORY\n"+            
        "              Directory containing input contours.\n"+
        "              Default is current directory.\n\n"+
        "       -o DIRECTORY, --output_data_dir=DIRECTORY\n"+           
        "              Directory where output data will be written.\n"+
        "              Default is current directory.\n\n"+
        "       -O STRING, --output_script=STRING\n"+     
        "              A bash script named STRING will be written \n"+
        "              to automate contour_tiling process.\n"+
        "              Default script name is '" + OUTPUT_SCRIPT + "'.\n\n"+
        "       -f STRING, --input_filename_prefix=STRING\n"+     
        "              The input contours will be read from\n"+
        "              'input_directory/STRING.min_section' to\n"+
        "              'input_directory/STRING.max_section'.\n"+
        "              Default is '" + PREFIX + "'.\n\n"+
        "       -I STRING, --ignore_contour=STRING\n"+           
        "              Contours with name STRING will not be processed\n"+
        "              and no output data will be written for these contours.\n"+
        "              Default is no ignored contours.\n\n"+
        "       -h, --help\n"+                     
        "              Print reconstruct2contourtiler man page.\n"+
        "\nJustin Kinney				2008/12/16\n";
  return message;
}

void Controls::parseCommandLine (int argc,char **argv)
{
  std::string const message = getUsageMessage();
  while (1)
  {
    static struct option long_options[] =
    {
      {"no_capping"                    , no_argument, & CAPPING_FLAG    , 0}, // default 1
      {"print_input_as_pts"            , no_argument, & PRINT_RAW_POINTS , 1},
      {"min_section"                   , required_argument, 0, 'n'},
      {"max_section"                   , required_argument, 0, 'x'},
      {"section_thickness"             , required_argument, 0, 't'},
      {"scale"                         , required_argument, 0, 's'},
      {"spline_samples_per_point"      , required_argument, 0, 'S'},
      {"deviation_threshold"           , required_argument, 0, 'd'},
      {"linear_threshold"              , required_argument, 0, 'm'},
      {"input_data_dir"                , required_argument, 0, 'i'},
      {"output_data_dir"               , required_argument, 0, 'o'},
      {"output_script"                 , required_argument, 0, 'O'},
      {"input_filename_prefix"         , required_argument, 0, 'f'},
      {"ignore_contour"                , required_argument, 0, 'I'},
      {"help"                          , no_argument      , 0, 'h'},
      {0, 0, 0, 0}
    };
    /* getopt_long stores the option index here. */
    int option_index = 0;

    int c = getopt_long (argc, argv, "n:x:t:s:d:m:hi:o:f:I:O:S:", long_options, &option_index);

    /* Detect the end of the options. */
    if (c == -1)
      break;

    char *eptr=NULL;
    std::string line;
    switch (c)
    {
      case 0:
        /* If this option set a flag, do nothing else now. */
        if (long_options[option_index].flag != 0)
          break;

      case 'n':
        MIN_SECTION = atoi(optarg);
        break;

      case 'x':
        MAX_SECTION = atoi(optarg);
        break;

      case 't':
        SECTION_THICKNESS = strtod(optarg,&eptr);
        break;

      case 'S':
        SPLINE_SAMPLES_PER_POINT = atoi(optarg);
        break;

      case 's':
        SCALE = strtod(optarg,&eptr);
        break;

      case 'd':
        DEVIATION_THRESHOLD = strtod(optarg,&eptr);
        break;

      case 'm':
        LINEAR_THRESHOLD = strtod(optarg,&eptr);
        break;

      case 'f':
        PREFIX = optarg;
        break;

      case 'O':
        OUTPUT_SCRIPT= optarg;
        break;

      case 'I':
        IGNORED_CONTOURS.push_back(optarg);
        break;

      case 'h':
        std::cout << message << std::endl;
        exit(1);
        break;

      case 'i':
        INPUT_DIR = processDir(optarg);
        break;

      case 'o':
        OUTPUT_DIR = processDir(optarg);
        break;

      case '?':
        /* getopt_long already printed an error message. */
        exit(1);
        break;

      default:
        std::cout << message << std::endl;
        abort ();
    }
  }

  /* Print any remaining command line arguments (not options). */
  if (optind < argc)
  {
    printf ("non-option ARGV-elements: ");
    while (optind < argc)
      printf ("%s ", argv[optind++]);
    putchar ('\n');
    std::cout << message << std::endl;
    exit(1);
  }
}

