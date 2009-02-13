#include <getopt.h>
#include <iostream>
#include <string>
#include <vector>

struct Parse
{
  double section_thickness;
  double scale;
  double deviation_threshold;
  int capping_flag;
  int min_section;
  int max_section;
  std::string input_dir;
  std::string output_dir;
  std::string prefix;
  std::vector<std::string> ignored_contours;
  Parse (void)
        : section_thickness(.05),scale(1000),deviation_threshold(2),
        capping_flag(1),min_section(60),max_section(160),input_dir("./"),
        output_dir("./"),prefix("Volumejosef"),ignored_contours()
  {}
  void parseCommandLine (int,char**,std::string const &);
  std::string getUsageMessage (void);

  bool contourIsIgnored (char const * const name) const
  {
    for (std::vector<std::string>::const_iterator i=ignored_contours.begin();i!=ignored_contours.end();i++)
    {
      if (!strcmp((*i).c_str(),name)) return true;
    }
    return false;
  }

  /** Create string from floating point number.
   * \param[in] i Number of interest.
   * \return String containing input number.
   */

  std::string d2str (double const & i)
  {
    char file[1024];
    sprintf(file,"%.15g",i);
    return std::string(file);
  }

  /** Create string from integer number.
   * \param[in] i Number of interest.
   * \return String containing input number.
   */

  std::string i2str (int const & i)
  {
    char file[1024];
    sprintf(file,"%i",i);
    return std::string(file);
  }

  /** Ensure pathname ends with exactly one '/'.
   * \param[in] ptr Arbitrary pathname.
   * \return Input pathname with '/' added if absent.
   */

  std::string processDir (char * ptr)
  {
    std::string pathname = ptr;
    if (pathname.find_last_of('/')==std::string::npos)
    {
      pathname.push_back('/');
      return pathname;
    }
    size_t pos = pathname.find_last_of('/');
    if ((pos+1)==pathname.size())
    {
      return pathname;
    }
    else
    {
      pathname.push_back('/');
      return pathname;
    }
  }

};

std::string Parse::getUsageMessage (void)
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
        "       capping_flag=1 to attempt end capping.capping_flag=0 to leave ends open.\n"+
        "       deviation_threshold is the maximum allowed deviation of the spline from raw\n"+
        "       contour points in scaled units. Set\n"+
        "       deviation_threshold to 0 to disable thresholding.\n"+
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
        "       -n NUM, --min_section=NUM\n"+               
        "              The starting section number in the section range.\n"+
        "              Default is '" + i2str(min_section) + "'.\n\n"+
        "       -x NUM, --max_section=NUM\n"+        
        "              The ending section number in the section range.\n"+
        "              Default is '" + i2str(max_section) + "'.\n\n"+
        "       -t NUM, --section_thickness=NUM\n"+              
        "              Thickness of the section in microns. Each section\n"+
        "              is assumed to be of identical thickness.\n"+
        "              Default is '" + d2str(section_thickness) + "'.\n\n"+
        "       -s NUM, --scale=NUM\n"+ 
        "              The output data will be scaled by NUM. If NUM is\n"+
        "              equal to 1 then the output data is in microns.\n"+
        "              Default is '" + d2str(scale) + "' so the output\n"+
        "              data is in nanometers.\n\n"+
        "       -d NUM, --deviation_threshold=NUM\n"+ 
        "              The input contours are filtered before output.\n"+
        "              The deviation of the input and output contours\n"+
        "              is constrained to be less than NUM where units\n"+
        "              are input units scaled by --scale value.\n"+
        "              Default is '" + d2str(deviation_threshold) + "'.\n\n"+
        "       -i DIRECTORY, --input_data_dir=DIRECTORY\n"+            
        "              Directory containing input contours.\n"+
        "              Default is current directory.\n\n"+
        "       -o DIRECTORY, --output_data_dir=DIRECTORY\n"+           
        "              Directory where output data will be written.\n"+
        "              Default is current directory.\n\n"+
        "       -f STRING, --input_filename_prefix=STRING\n"+     
        "              The input contours will be read from\n"+
        "              'input_directory/STRING.min_section' to\n"+
        "              'input_directory/STRING.max_section'.\n"+
        "              Default is '" + prefix + "'.\n\n"+
        "       -I STRING, --ignore_contour=STRING\n"+           
        "              Contours with name STRING will not be processed\n"+
        "              and no output data will be written for these contours.\n"+
        "              Default is no ignored contours.\n\n"+
        "       -h, --help\n"+                     
        "              Print reconstruct2contourtiler man page.\n"+
        "\nJustin Kinney				2008/12/16\n";
  return message;
}

void Parse::parseCommandLine (int argc,char **argv,std::string const & message)
{

  while (1)
  {
    static struct option long_options[] =
    {
      {"no_capping"                    , no_argument, & capping_flag    , 0}, // default 1
      {"min_section"                   , required_argument, 0, 'n'},
      {"max_section"                   , required_argument, 0, 'x'},
      {"section_thickness"             , required_argument, 0, 't'},
      {"scale"                         , required_argument, 0, 's'},
      {"deviation_threshold"           , required_argument, 0, 'd'},
      {"input_data_dir"                , required_argument, 0, 'i'},
      {"output_data_dir"               , required_argument, 0, 'o'},
      {"input_filename_prefix"         , required_argument, 0, 'f'},
      {"ignore_contour"                , required_argument, 0, 'I'},
      {"help"                          , no_argument      , 0, 'h'},
      {0, 0, 0, 0}
    };
    /* getopt_long stores the option index here. */
    int option_index = 0;

    int c = getopt_long (argc, argv, "n:x:t:s:d:hi:o:f:I:", long_options, &option_index);

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
        min_section = atoi(optarg);
        break;

      case 'x':
        max_section = atoi(optarg);
        break;

      case 't':
        section_thickness = strtod(optarg,&eptr);
        break;

      case 's':
        scale = strtod(optarg,&eptr);
        break;

      case 'd':
        deviation_threshold = strtod(optarg,&eptr);
        break;

      case 'f':
        prefix = optarg;
        break;

      case 'I':
        ignored_contours.push_back(optarg);
        break;

      case 'h':
        std::cout << message << std::endl;
        exit(1);
        break;

      case 'i':
        input_dir = processDir(optarg);
        break;

      case 'o':
        output_dir = processDir(optarg);
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

class void_list
{
public:
  void_list *previous;
  void_list *next;
  void *data;
};

class Parameters
{
public:
  double plot_rad_int;	// radius of curvature sampling interval for plotting sampling function
  double dmin;			// minimum spline sampling distance
  double dmax;			// maximum spline sampling distance
  double max_rad;		// calculated radius of curvature of spline will saturate at this value
  bool diag;			// set true to print diagnostic files
  int num;				// num is the # samples of the splines between contour points
  // sequential triplets of sampled points are used to
  // compute the radius of curvature , e.g. num/2 between contour points
  double T;				// sample period (= time to traverse dmax)
  double amax;			// max radial acceleration
};

class Object
{
public:
  char name[128];
  int min_section,max_section;
  void_list *contours;
  Object(char *str,int sec);
  ~Object(void);
};

Object::Object(char *str,int sec)
{
  strcpy(name,str);
  contours = NULL;
  min_section = sec;
  max_section = sec;
}

Object::~Object(void){
  void_list *p,*q;
  p=contours;
  while (p!=NULL) {
    q=p->next;
    delete p;
    p=q;
  }
}

class Point
{
public:
  double x,y,z;
  Point(char *str,int section,double thickness,double *transform);
  Point(double xval,double yval,double zval);
};

Point::Point(double xval, double yval, double zval)
{
  x = xval;
  y = yval;
  z = zval;
}

Point::Point(char *str, int section, double thickness,double *t)
{
  char val[80];
  char *eptr;
  int i;
  double xval,yval;

  // set z coordinate
  z = section*thickness;

  // get past 'points'
  while (strchr(" points=\"\t",*str)!=NULL) {str++;}

  // grab x coordinate
  i=0;
  while (strchr("0123456789+-eE.",*str)!=NULL) { val[i++] = *str++; }
  val[i]=0;
  xval = strtod(val,&eptr);
  if (val==eptr) {
    x=y=z=0;
    printf("Error in reading x coordinate\n");
    printf("str =%s\n",str);
    return;
  }

  // grab y coordinate
  while (strchr(" \t,",*str)!=NULL) { str++; }
  i=0;
  while (strchr("0123456789+-eE.",*str)!=NULL) { val[i++] = *str++; }
  val[i]=0;
  yval = strtod(val,&eptr);
  if (val==eptr) {
    x=y=z=0;
    printf("Error in reading y coordinate\n");
    return;
  }
  //	if(section==64){
  //	printf("xcoef = %g %g %g %g %g %g\n",t[0],t[1],t[2],t[3],t[4],t[5]);
  //	printf("ycoef = %g %g %g %g %g %g\n\n",t[6],t[7],t[8],t[9],t[10],t[11]);
  //	}
  x = t[0]+t[1]*xval+t[2]*yval+t[3]*xval*yval+t[4]*xval*xval+t[5]*yval*yval;
  y = t[6]+t[7]*xval+t[8]*yval+t[9]*xval*yval+t[10]*xval*xval+t[11]*yval*yval;
}

class Contour
{
public:
  char name[128];
  int section,num_interp_points;
  void_list *raw_points,*rawend;
  void_list *interp_points,*interpend;
  double *deviations;
  Contour(char* str, int sec);
  ~Contour(void);
  void removeDuplicates(void);
  void clearSpline(void);
  void linearlyInterp(double,double);
  void addPreviousRaw(void);
  void addPreviousInterp(void);
  int getNumRawPoints(void)
  {
    int i=0;
    void_list *p=raw_points;
    while (p!=NULL)
    {
      void_list *q=p->next;
      p=q;
      i++;
    }
    return i;
  }
  void checkRawPoints(void)
  {
    void_list *p=raw_points;
    while (p!=NULL)
    {
      void_list *q=p->next;
      p=q;
    }
  }
};

Contour::~Contour(void)
{
  void_list *p,*q;
  p=raw_points;
  while (p!=NULL) {
    q=p->next;
    delete (Point*)p->data;
    delete p;
    p=q;
  }
  p=interp_points;
  while (p!=NULL) {
    q=p->next;
    delete (Point*)p->data;
    delete p;
    p=q;
  }
  if (deviations!=NULL)
    delete[] deviations;
}

Contour::Contour(char *str, int sec)
:deviations(NULL)
{
  char *ptr = str;
  while(ptr!=NULL){
    //      ptr=strchr(val,'#');
    ptr=strpbrk(str,"#()");
    if(ptr!=NULL){
      *ptr='_';
    }
  }

  strcpy(name,str);
  section = sec;
  raw_points = NULL;
  rawend=NULL;
  interp_points = NULL;
}

void Contour::clearSpline(void)
{
  void_list *p,*q;
  p=interp_points;
  while(p!=NULL) {
    delete (Point*)p->data;
    q=p->next;
    delete p;
    p=q;
  }
  interp_points = NULL;
  interpend = NULL;
}

void Contour::addPreviousRaw(void){
  void_list *q,*prev;
  prev=NULL;
  rawend=NULL;
  // for each point
  for (q=raw_points;q!=NULL;q=q->next) {
    q->previous = prev;
    prev = q;
    rawend = q;
  }
}

void Contour::addPreviousInterp(void){
  void_list *q,*prev;
  prev=NULL;
  interpend=NULL;
  // for each point
  for (q=interp_points;q!=NULL;q=q->next) {
    q->previous = prev;
    prev = q;
    interpend = q;
  }
}

void_list * removeLink(void_list* L) {
  void_list *q;
  // and remove face from candidate face list
  if (L->previous!=NULL) {
    if (L->next!=NULL) {
      // if both previous and next exist
      (L->previous)->next = L->next;
      (L->next)->previous = L->previous;
    } else {
      // if previous exists and next does not
      (L->previous)->next = NULL;
    }
  } else {
    if (L->next!=NULL) {
      // if previous does not exist and next does
      (L->next)->previous = NULL;
    } // else { // if neither previous nor next exists }
  }
  // update pointer
  q=L->next;
  delete L;
  return q;
}

void_list * deletePoint(void_list *q,void_list *p,void_list *&ptr)
{
  Point *pt1,*pt2;
  pt1=(Point*)q->data;
  pt2=(Point*)p->data;
  // if points are identical
  if ((pt1->x==pt2->x)&&(pt1->y==pt2->y))
  {
    // delete point
    delete pt1;
    // adjust list pointer
    if (q==ptr) {ptr = q->next; }
    // remove current point from list
    q=removeLink(q);
  }
  return p;
}

void Contour::removeDuplicates(void)
{
  void_list *q,*ptr;
  //Point *pt1,*pt2;
  ptr = raw_points;
  q=raw_points;
  while (q->next!=NULL) {q=deletePoint(q,q->next,ptr);}
  // adjust pointer
  raw_points = ptr;

  // compare first and last link in list
  q=deletePoint(raw_points,rawend,ptr);
  // adjust pointer
  raw_points = ptr;
}

void_list * interpPoints(void_list *q,void_list *p,void_list *&ptr,
                         double maxdev,double scale,int flag){
  void_list *pp;
  double dist,distx,disty,x,y,count;
  int num;
  Point *v,*pt1,*pt2;
  pt1=(Point*)q->data;
  pt2=(Point*)p->data;
  // compute distance between points
  distx = pt2->x-pt1->x;
  disty = pt2->y-pt1->y;
  dist = sqrt(distx*distx+disty*disty);
  num = (int)(dist/(maxdev/scale)/3);
  if (num)
  {
    // linearly interpolate num evenly spaced points
    count = num+1;
    ptr = q;
    while(num)
    {
      // insert point
      x = pt1->x+(count-(double)num)/count*distx;
      y = pt1->y+(count-(double)num)/count*disty;
      v = new Point(x,y,pt1->z);
      pp = new void_list();
      pp->next = p;
      pp->previous = ptr;
      pp->data = (void*)v;
      ptr->next = pp;
      p->previous = pp;
      // decrement num
      num--;
      ptr = pp;
    }
    if(flag){ptr->next=NULL;}
  }
  return p;
}

void Contour::linearlyInterp(double maxdev,double scale)
{
  void_list *q,*ptr=NULL;
  q=raw_points;

  if (q==NULL)
  {
    printf("WTF? raw_points==NULL.\n");
    exit(1);
  }

  while (q->next!=NULL)
  {
    q=interpPoints(q,q->next,ptr,maxdev,scale,0);
  }
  // compare first and last link in list
  q=interpPoints(rawend,raw_points,ptr,maxdev,scale,1);
  // adjust pointer
  rawend=ptr;
}

class Histogram
{
public:
  double min,max,mean,stddev,sum;
  int count[16],num;
  Histogram(void);
  void update(void_list*,int);
};

Histogram::Histogram(void)
{
  int i;
  for(i=0;i<16;i++){
    count[i]=0;
  }
  min=1e30;
  max=0;
  mean=stddev=sum=0;
  num=0;
}

void Histogram::update(void_list* q,int total) {
  ///// update deviation distance statistics /////
  int i;
  double d;
  Contour* c=(Contour*)q->data;
  for (i=0;i<total;i++) {
    d=c->deviations[i];
    // update min and max deviation distance
    if(d<min) {min=d;}
    if(d>max) {max=d;}
    num++;
    sum+=d;
  }
}

class SplinePoint
{
public:
  double t,x,y,r,intfac;
};

class Weights
{
public:
  double *bx,*by;
};

