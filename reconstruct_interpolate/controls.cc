#include <cassert>
#include <stdlib.h>
#include <iostream>

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
  :OBJECT_RESERVE_SIZE                 (10000),
  CAPPING_FLAG                         (1),
  MIN_SECTION                          (0),
  MAX_SECTION                          (0),
  PRINT_DETAILED_INFO                  (0),
  RETURN_RAW_CONTOUR_POINTS            (0),
  RETURN_INTERPOLATED_RAW_POINTS       (0),
  MIN_PT_PER_CONTOUR_THRESHOLD         (4),
  SIM_ANNEAL_MOVES_PER_ITERATION       (3),
  SIM_ANNEAL_MAX_NUM_MOVE_ATTEMPTS     (1000),
  SIM_ANNEAL_INTEGRATION_STEP          (6),
  SIM_ANNEAL_INTEGRATION_STEP_FACTOR   (1.0/static_cast<double>(SIM_ANNEAL_INTEGRATION_STEP)),
  MAX_DEVIATION_ADJUSTMENTS            (10),
  ADDITIONAL_POINTS_FACTOR             (0.5),
  SIM_ANNEAL_HIGH_TEMP                 (100),
  SIM_ANNEAL_TEMP_SCALE                (0.9),
  SIM_ANNEAL_BOLTZMAN                  (1.0/1.380E23),
  SIM_ANNEAL_CURVATURE_ENERGY_GAIN     (1E3),
  SIM_ANNEAL_CURVATURE_ENERGY_EXPONENT (1),
  SIM_ANNEAL_PROXIMITY_ENERGY_GAIN     (1E0),
  SIM_ANNEAL_PROXIMITY_ENERGY_EXPONENT (1),
  SIM_ANNEAL_MEAN_AMPLITUDE_NOISE      (0.1),
  MAXIMUM_RADIUS_OF_CURVATURE          (1E1),
  SECTION_THICKNESS                    (.05),
  OUTPUT_SCALE_FACTOR                  (1),
  DEVIATION_THRESHOLD                  (0),
  EPSILON                              (1E-7),
  MAX_SAMPLE_INTERVAL                  (SECTION_THICKNESS),
  MIN_SAMPLE_INTERVAL                  (SECTION_THICKNESS/5.0),
  INPUT_DIR                            ("./"),
  OUTPUT_DIR                           ("./"),
  OUTPUT_SER_PREFIX                    (""),
  PREFIX                               (""),
  EXCLUDED_CONTOURS                    (),
  INCLUDED_CONTOURS                    (),
  OUTPUT_SCRIPT                        ("mesh_and_convert.sh"),
  MULTI_PART_SUFFIX                    ("_part\%d_")
{}

/** Create usage message.
 * \return usage message.
 */

std::string Controls::getUsageMessage (void)
{
  std::string message = "\n";
  message=message+
        "NAME\n"+
        "       reconstruct_interpolate - generate contour_tiler input files from Reconstruct3D contours\n"+
        "\nSYNOPSIS\n"+
        "       reconstruct_interpolate [options]\n"+
        "\nDESCRIPTION\n"+
        "       Converts Reconstruct3D XML contour format to contour_tiler input format.\n"+
        "       All files in input directory are assumed to be\n"+
        "       of the form filename_prefix.section#.\n"+
        "       min_section to max_section is the section range to be converted.\n"+
        "       Section_thickness should be in same scale as x,y contour points.\n"+
        "       x,y and z coordinates of sampled splines will be multipled by scale in output.\n"+
        "\nEXAMPLES\n"+
        "       reconstruct_interpolate -i ./examples1 -f Volumejosef -n 98 -x 102 -t .05 -s 1000 -o examples1/contour_tiler_output -d 2\n"+
        "              Read contours from directory './examples1' and write contour_tiler output\n"+
        "              files to directory 'examples1/contour_tiler_output'. The input contour files have the\n"+
        "              name scheme 'Volumejosef.#' where # is the section number which varies from\n"+
        "              10 to 100. The distance between contours in the direction of sectioning\n"+
        "              is .050 microns. The contour_tiler output data will be in nanometers as\n"+
        "              dictated by the 1000 scaling. Capping directives will be included in the\n"+
        "              output data. The interpolated contours will not deviate from the input contours\n"+
        "              by more than 2 (nanometers since scaling is 1000).\n"+
        "\nOPTIONS\n"+
        "       --no_capping\n"+
        "              Do not include directives in the output data to cap the meshes,\n"+
        "              thus creating open surface meshes.\n"+ 
        "              Default is to cap the meshes at the minimum and maximum sections\n"+
        "              thereby creating closed surface meshes.\n\n"+
        "       --print_detailed_info\n"+
        "              Print raw contour points, control points, path parameter\n"+
        "              values, etc as .log files in output directory.\n"+ 
        "              Default is to do nothing.\n\n"+
        "       --return_raw_contour_points\n"+
        "              Return input contour points unadulterated.\n"+ 
        "              Default is to fit splines to contour points and resample.\n\n"+
        "       --return_interpolated_raw_points\n"+
        "              Return input contour points linearly interpolated\n"+ 
        "              to satisfy minimum and maximum sample interval constraints.\n"+
        "              Default is to fit splines to contour points and resample.\n\n"+
        "       -n NUM, --min_section=NUM\n"+               
        "              The starting section number in the section range.\n"+
        "              Default is '" + i2str(MIN_SECTION) + "'.\n\n"+
        "       -x NUM, --max_section=NUM\n"+        
        "              The ending section number in the section range.\n"+
        "              Default is '" + i2str(MAX_SECTION) + "'.\n\n"+
        "       -t NUM, --section_thickness=NUM\n"+              
        "              Thickness of the sections. Each section\n"+
        "              is assumed to be of identical thickness.\n"+
        "              Default is '" + d2str(SECTION_THICKNESS) + "'.\n\n"+
        "       -r NUM, --min_point_per_contour=NUM\n"+              
        "              Contours with fewer than NUM output sample points\n"+
        "              will be omitted from output data files and statistics.\n"+
        "              NUM=0 implies that all contours are included.\n"+
        "              Default is '" + d2str(MIN_PT_PER_CONTOUR_THRESHOLD) + "'.\n\n"+
        "       -S NUM, --additional_points_factor=NUM\n"+              
        "              The number of sample points in contour is determined\n"+
        "              by NUM where 0<=NUM<=1. When NUM is 0 the number of contour points\n"+
        "              is equal to the contour length divided by the maximum sample interval.\n"+
        "              When NUM is 1 the number of contour points is equal to the contour length\n"+
        "              divided by the minimum sample interval.\n"+
        "              Default is '" + d2str(ADDITIONAL_POINTS_FACTOR) + "'.\n\n"+
        "       -X NUM, --max_sample_interval=NUM\n"+ 
        "              The linear distance between sampled contour points\n"+
        "              will be less than NUM in most cases.\n"+
        "              Default is '" + d2str(MAX_SAMPLE_INTERVAL) + "'.\n\n"+
        "       -Y NUM, --min_sample_interval=NUM\n"+ 
        "              The linear distance between sampled contour points\n"+
        "              will be greater than NUM in most cases.\n"+
        "              Default is '" + d2str(MIN_SAMPLE_INTERVAL) + "'.\n\n"+
        "       -s NUM, --output_scale=NUM\n"+ 
        "              The output data will be scaled by NUM.\n"+
        "              Default is '" + d2str(OUTPUT_SCALE_FACTOR) + "' which means\n"+
        "              output data has same scale as input data.\n\n"+
        "       -d NUM, --deviation_threshold=NUM\n"+ 
        "              The input contours are filtered before output.\n"+
        "              The deviation of the input and output contours\n"+
        "              is constrained to be less than NUM where units\n"+
        "              are input units scaled by --scale value.\n"+
        "              Default is '" + d2str(DEVIATION_THRESHOLD) + "'\n"+
        "              which means that no threshold is enforced.\n\n"+
        "       -a NUM, --curvature_gain=NUM\n"+ 
        "              A sample point contributes to energy in a manner\n"+
        "              proportional to NUM and inversely proportional to curvature.\n"+
        "              Default is '" + d2str(SIM_ANNEAL_CURVATURE_ENERGY_GAIN) + "'.\n\n"+
        "       -b NUM, --curvature_exponent=NUM\n"+ 
        "              A sample point contributes to energy in a manner\n"+
        "              inversely proportional to curvature to the NUM power.\n"+
        "              Default is '" + i2str(SIM_ANNEAL_CURVATURE_ENERGY_EXPONENT) + "'.\n\n"+
        "       -c NUM, --proximity_gain=NUM\n"+ 
        "              A sample point contributes to energy in a manner\n"+
        "              proportional to NUM and distance between samples.\n"+
        "              Default is '" + d2str(SIM_ANNEAL_PROXIMITY_ENERGY_GAIN) + "'.\n\n"+
        "       -e NUM, --proximity_exponent=NUM\n"+ 
        "              A sample point contributes to energy in a manner\n"+
        "              proportional to distance between samples to the NUM power.\n"+
        "              Default is '" + i2str(SIM_ANNEAL_PROXIMITY_ENERGY_EXPONENT) + "'.\n\n"+
        "       -T NUM, --high_temp=NUM\n"+ 
        "              Starting high temperature of simulated annealing.\n"+
        "              Default is '" + d2str(SIM_ANNEAL_HIGH_TEMP) + "'.\n\n"+
        "       -i DIRECTORY, --input_data_dir=DIRECTORY\n"+            
        "              Directory containing input contours.\n"+
        "              Default is current directory.\n\n"+
        "       -o DIRECTORY, --output_data_dir=DIRECTORY\n"+           
        "              Directory where output data will be written.\n"+
        "              Default is current directory.\n\n"+
        "       -w STRING, --output_ser_prefix=STRING\n" +
        "              Write interpolated traces to SER files using specified prefix. Do not output raw points.\n\n"+
        "       -O STRING, --output_script=STRING\n"+     
        "              A bash script named STRING will be written \n"+
        "              to automate contour_tiling process.\n"+
        "              Default script name is '" + OUTPUT_SCRIPT + "'.\n\n"+
        "       -f STRING, --input_filename_prefix=STRING\n"+     
        "              The input contours will be read from\n"+
        "              'input_directory/STRING.min_section' to\n"+
        "              'input_directory/STRING.max_section'.\n"+
        "              Default is '" + PREFIX + "'.\n\n"+
        "       -E STRING, --exclude_contour=STRING\n"+           
        "              Contours with name STRING will not be processed\n"+
        "              and no output data will be written for these contours.\n"+
        "              Default is no excluded contours.\n\n"+
        "       -I STRING, --include_contour=STRING\n"+           
        "              Contours with name STRING will be included and processed\n"+
        "              and output data will be written for only these contours.\n"+
        "              Default is to include all contours.\n\n"+
        "       -M STRING, --multi_part_suffix=STRING\n"+     
        "              Objects with intermediate empty sections (thereby\n"+
        "              requiring multiple tiling steps) will have name\n"+
        "              and associated output files appended with STRING\n"+
        "              plus incremented index.\n"+
        "              Default script name is '" + MULTI_PART_SUFFIX + "'.\n\n"+
        "       -h, --help\n"+                     
        "              Print reconstruct_interpolate man page.\n"+
        "\nJustin Kinney				Sep 10 2009\n";
  return message;
}

/** Parse command line.
 * \param[in] argc Argc.
 * \param[in] argv Argv.
 */

void Controls::parseCommandLine (int argc,char **argv)
{
  std::string const message = getUsageMessage();
  while (1)
  {
    static struct option long_options[] =
    {
      {"no_capping"                    , no_argument, & CAPPING_FLAG    , 0}, // default 1
      {"print_detailed_info"           , no_argument, & PRINT_DETAILED_INFO , 1},
      {"return_raw_contour_points"     , no_argument, & RETURN_RAW_CONTOUR_POINTS , 1},
      {"return_interpolated_raw_points", no_argument, & RETURN_INTERPOLATED_RAW_POINTS , 1},
      {"min_section"                   , required_argument, 0, 'n'},
      {"max_section"                   , required_argument, 0, 'x'},
      {"min_point_per_contour"         , required_argument, 0, 'r'},
      {"section_thickness"             , required_argument, 0, 't'},
      {"output_scale"                  , required_argument, 0, 's'},
      {"additional_points_factor"      , required_argument, 0, 'S'},
      {"max_sample_interval"           , required_argument, 0, 'X'},
      {"min_sample_interval"           , required_argument, 0, 'Y'},
      {"deviation_threshold"           , required_argument, 0, 'd'},
      {"linear_threshold"              , required_argument, 0, 'm'},
      {"input_data_dir"                , required_argument, 0, 'i'},
      {"output_data_dir"               , required_argument, 0, 'o'},
      {"output_ser_prefix"             , required_argument, 0, 'w'},
      {"output_script"                 , required_argument, 0, 'O'},
      {"multi_part_suffix"             , required_argument, 0, 'M'},
      {"input_filename_prefix"         , required_argument, 0, 'f'},
      {"exclude_contour"               , required_argument, 0, 'E'},
      {"include_contour"               , required_argument, 0, 'I'},
      {"curvature_gain"                , required_argument, 0, 'a'},
      {"curvature_exponent"            , required_argument, 0, 'b'},
      {"proximity_gain"                , required_argument, 0, 'c'},
      {"proximity_exponent"            , required_argument, 0, 'e'},
      {"high_temp"                     , required_argument, 0, 'T'},
      {"help"                          , no_argument      , 0, 'h'},
      {0, 0, 0, 0}
    };
    /* getopt_long stores the option index here. */
    int option_index = 0;

    int c = getopt_long (argc, argv, "a:b:c:e:n:x:t:r:s:d:m:hi:o:f:w:T:I:E:O:S:X:Y:M:", long_options, &option_index);

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

      case 'a':
        SIM_ANNEAL_CURVATURE_ENERGY_GAIN = strtod(optarg,&eptr);
        break;

      case 'b':
        SIM_ANNEAL_CURVATURE_ENERGY_EXPONENT = strtod(optarg,&eptr);
        break;

      case 'c':
        SIM_ANNEAL_PROXIMITY_ENERGY_GAIN = strtod(optarg,&eptr);
        break;

      case 'e':
        SIM_ANNEAL_PROXIMITY_ENERGY_EXPONENT = strtod(optarg,&eptr);
        break;

      case 'T':
        SIM_ANNEAL_HIGH_TEMP = strtod(optarg,&eptr);
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

      case 'r':
        MIN_PT_PER_CONTOUR_THRESHOLD=strtod(optarg,&eptr);
        break;

      case 'S':
        ADDITIONAL_POINTS_FACTOR = strtod(optarg,&eptr);
        break;

      case 'X':
        MAX_SAMPLE_INTERVAL = strtod(optarg,&eptr);
        break;

      case 'Y':
        MIN_SAMPLE_INTERVAL = strtod(optarg,&eptr);
        break;

      case 's':
        OUTPUT_SCALE_FACTOR = strtod(optarg,&eptr);
        break;

      case 'd':
        DEVIATION_THRESHOLD = strtod(optarg,&eptr);
        break;

      case 'f':
        PREFIX = optarg;
        break;

      case 'O':
        OUTPUT_SCRIPT = optarg;
        break;

      case 'M':
        MULTI_PART_SUFFIX = optarg;
        break;

      case 'E':
        EXCLUDED_CONTOURS.push_back(optarg);
        break;

      case 'I':
        INCLUDED_CONTOURS.push_back(optarg);
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

      case 'w':
        OUTPUT_SER_PREFIX = optarg;
        printf("Writing output to SER...\n");
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
  // validate control values
  validate();
}

