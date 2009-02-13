//####################################################
//####################################################
class Object
{
public:
  std::string name;   // object name
  std::vector<int> s; // sequence numbers of mesh files
  Object(std::string const &str)
  {
    name = str;
  }
};

struct lts
{
  bool operator()(const std::string s1, const std::string s2) const
  {
    return s1 < s2;
  }
};

typedef std::map<std::string,Object*,lts,std::allocator<Object*> > map_so;
typedef std::map<std::string,Object*,lts,std::allocator<Object*> >::iterator mso_iterator;

class Cabinet{
public:
  std::string indir,outdir;
  std::ofstream F,M;
  std::vector<std::string> files;	// array of input file names
  int num_files;	// number of input files found
  int range[2]; // min,max
  map_so so;
  Cabinet(void)
  {
    num_files=0;
  }
  void scanDir(void);
  void parse(int,char**,std::string);
  bool endswith(std::string const &, std::string const &);
  std::string strip_tail(std::string const &,uint);
  void saveObject(std::string const&);
  std::string getName(std::string const&);
  Object * findName(std::string const &);
  Object * createObject(std::string const &);
  int getSeq(std::string const&);
  void validateObjects(void);
  void writeMain(void);
  void writeScript(void);
};

void Cabinet::writeScript(void)
{
  // open script
  std::string str = outdir + SCRIPT;
  F.open(str.c_str());
  if(F.is_open()==false){
    fprintf(stderr,"\nCouldn't open output file %s\n",str.c_str());
    exit(0);
  }
  // for each object
  for (mso_iterator i=so.begin();i!=so.end();i++)
  {
    Object *p=(*i).second;
    // for each sequence
    for (std::vector<int>::iterator j=p->s.begin();j!=p->s.end();j++)
    {
      F << "echo \"converting "
            << indir << p->name << "_" << *j << ".mdl to mdl format\"\n"
            << "mesh2mcell " << p->name << " " << indir << p->name << "_" 
            << *j << ".mesh > " << outdir << p->name << "_" << *j << ".mdl\n";
    }
  }
  F << "cd " << outdir << "\n" <<
        "echo \"initializing sequence.mdl\"\n" <<
        "echo \"i=" << range[0] << "\" > sequence.mdl\n" <<
        "a=" << range[0] << "\n" <<
        "b=" << range[1] << "\n" <<
        "while [ \"$a\" -lt \"$b\" ]\n" <<
        "do\n" <<
        "  echo \"running mcell $a\"\n" <<
        "  ~/bin/mcell3.1-opteron main.mdl\n" <<
        "  let \"a+=1\"\n" <<
        "  done\n" <<
        "echo \"cleaning up\"\n" <<
        "rm *.mdl\n" <<
        "rm chkpt\n" <<
        "rm " << SCRIPT <<
        "\n";
  F.close();
}

void Cabinet::writeMain(void)
{
  // open main.mdl
  std::string str = outdir + MAIN;
  M.open(str.c_str());
  if(M.is_open()==false){
    fprintf(stderr,"\nCouldn't open output file %s\n",str.c_str());
    exit(0);
  }
  // prep main
  M << "ITERATIONS = " << range[1] << "\n" << 
        "TIME_STEP = 1e-6\n" << 
        "EFFECTOR_GRID_DENSITY = 15000  \n" << 
        "\n" << 
        "PARTITION_X = [[0.00001 TO 10000.00001 STEP 2000.00000001]]\n" << 
        "PARTITION_Y = [[0.00001 TO 10000.00001 STEP 2000.00000001]]\n" << 
        "PARTITION_Z = [[0.00001 TO 10000.00001 STEP 2000.00000001]]\n" << 
        "\n" << 
        "CHECKPOINT_INFILE = \"chkpt\"\n" << 
        "CHECKPOINT_OUTFILE = \"chkpt\"\n" << 
        "CHECKPOINT_ITERATIONS = 1\n" << 
        "\n" << 
        "/* defines i*/\n" << 
        "INCLUDE_FILE=\"sequence.mdl\"\n" << 
        "\n" <<
        "/* read geometry using i*/\n";
  // for each object 
  for (mso_iterator i=so.begin();i!=so.end();i++)
  {
    M << "sprintf(geo_f,\"" << (*i).first << "_%g.mdl\",i)\n" << 
          "INCLUDE_FILE=geo_f\n"; 
  }
  M << "\n" << 
        "/* write incremented value of i*/\n" << 
        "fname = fopen(\"sequence.mdl\",\"w\")\n" << 
        "fprintf(fname,\"i=%g\\n\",i+1)\n" << 
        "fclose(fname)\n" << 
        "\n" << 
        "INSTANTIATE a OBJECT {\n";
  // for each object 
  for (mso_iterator i=so.begin();i!=so.end();i++)
  {
    M << (*i).first << " OBJECT " 
          << (*i).first << " {}\n";
  }
  M << "}\n" << "\n" << 
        "VIZ_OUTPUT\n" <<
        "{\n" << 
        "  MODE=DREAMM_V3\n" << 
        "  FILENAME = \"./a\"\n" << 
        "  MESHES\n" << 
        "  {\n" << 
        "   NAME_LIST {ALL_MESHES }\n" << 
        "   ITERATION_NUMBERS { ALL_DATA @ ALL_ITERATIONS}\n" << 
        "  }\n" << 
        "}\n";
  M.close();
}

std::string Cabinet::strip_tail(std::string const &str, uint len)
{
  if (str.size() < len)
    return str;

  return str.substr(0, str.size() - len);
}

bool Cabinet::endswith(std::string const &str, std::string const &suffix)
{
  if (str.size() < suffix.size())
    return false;

  return str.substr(str.size() - suffix.size()) == suffix;
}

void Cabinet::parse(int argc,char **argv,std::string message){

  // if no arguments passed
  if(argc!=3){
    cout << message << endl;
    exit(0);
  }

  int c;
  opterr=0;
  while((c=getopt(argc,argv,"h")) != -1){
    switch(c){ 
      case 'h':
        // print help message
        cout << message << endl;
        exit(0);
      case '?':
        if (optopt == 'c')
          fprintf (stderr,"Option -%c requires an argument.\n", optopt);
        else if (isprint (optopt))
          fprintf (stderr, "Unknown option '-%c.\n", optopt);
        else
          fprintf (stderr,
                   "Unknown option character '\\x%x'.\n",
                   optopt);
        exit(0);
      default:
        cout << message << endl;
        abort ();
    }
  }

  if(optind<argc){
    bool found=false;
    indir=outdir="";
    for (int index=optind;index<argc;index++){
      // determine if file or folder was passed
      struct stat stat_p;
      stat(argv[index],&stat_p);
      if((stat_p.st_mode & S_IFDIR)==false){
        cout << "\n\nvoid parse: Error."
              << " argument #" << index
              << "=" << argv[index]
              << " was not classified as a directory.\n";
        cout << "NAME\n"
              << "       mesh2groupdx - multiple meshes into single dx file\n"
              << "\nSYNOPSIS\n       mesh2groupdx INPUTDIR OUTPUTDIR\n";
        exit(0);
      }

      // adjust input directory
      char filename[1028];
      strcpy(filename,argv[index]);
      char *temp=strrchr(filename,'/');
      if(!temp) {strcat(filename,"/");}
      else if(*++temp) {strcat(filename,"/");}
      if(found==false){indir=filename;found=true;}
      else {outdir=filename;}
    }
    if(indir=="" || outdir==""){
      fprintf (stderr,"No output directory argument found on command line.\n");
      exit(0);
    }
  } else {
    fprintf (stderr,"No input file argument found on command line.\n");
    exit(0);
  }
}

void Cabinet::validateObjects(void)
{
  int min = -1,max = -1;
  // for each object
  for (mso_iterator i=so.begin();i!=so.end();i++)
  {
    Object *p = (*i).second;
    // sort sequences
    sort(p->s.begin(),p->s.end());
    // assert that # sequences equals sequence range
    uint range_count = p->s.back()-p->s.front()+1;
    if(range_count!=p->s.size())
    {
      cout << "Cabinet::validateObjects: Error: "
            << "sequence numbers not continuous for object "
            << p->name << endl;
      for (std::vector<int>::iterator j=p->s.begin();j!=p->s.end();j++)
      {
        cout << "sequence # " << *j << endl;
      }
      exit(0);
    }
    // check sequence overlap for objects
    if (min<0){min = p->s.front();}
    if (max<0){max = p->s.back();}
    if (p->s.front()!=min || p->s.back()!=max)
    {
      cout << "Cabinet::validateObjects: Error: "
            << "sequence range does not perfectly overlap previous object for object "
            << p->name << endl;
      for (std::vector<int>::iterator j=p->s.begin();j!=p->s.end();j++)
      {
        cout << "sequence # " << *j << endl;
      }
      exit(0);
    }
  }
  range[0]=min;
  range[1]=max;
}

std::string Cabinet::getName(std::string const &str)
{
  // extract name (== everything before last '_')
  if (str.find_last_of('_') == std::string::npos)
  {
    cout << "std::string Cabinet: Error:"
          << " '_' not found in file name==\""
          << str << "\"\n";
    exit(0);
  }

  return str.substr(0, str.find_last_of('_'));
}

Object * Cabinet::findName(std::string const &str)
{
  // if element exists in map with matching key, then get Object pointer
  if (so.count(str)>0){ return so[str]; }
  return NULL;
}

Object * Cabinet::createObject(std::string const &str)
{
  Object *p = new Object(str);
  so[str]=p;
  return p;
}

int Cabinet::getSeq(std::string const &str)
{
  // extract sequence# ( == number after last '_' and before '.mesh')
  std::string s = strip_tail(str,SUFFIX.size());
  std::string seq = s.substr(s.find_last_of('_')+1,s.size());
  char *eptr;
  return static_cast<int>(strtod(seq.c_str(),&eptr));
}

void Cabinet::saveObject(std::string const &s)
{
  // extract object name from file name
  std::string name = getName(s);
  // look for object with the same name
  Object *p = findName(name);
  // if no object found with same name, then create new object
  if (p==NULL)
  {
    p = createObject(name);
  }
  // extract sequence# ( == number after last '_' and before '.mesh')
  int seq = getSeq(s);
  // save sequence# to Object
  p->s.push_back(seq);
}

void Cabinet::scanDir(void)
{
  // pointer to dirent structure
  struct dirent *pent;
  char filename[1024];
  strcpy(filename,indir.c_str());
  // pointer to a directory data structure
  DIR *pdir = opendir(filename);
  if (!pdir) {printf("Error. Could not open %s.\n",filename);exit(1);}
  cout << "\nInput directory found: " << filename << endl << endl;
  while ((pent=readdir(pdir)))
  {
    // copy char array to string
    std::string str = pent->d_name;
    // if found
    if (endswith(str,".mesh")==true)
    {
      // save filename
      files.push_back(str);
      // update index
      num_files++;
      // save object name and mesh file sequence number
      saveObject(str);
      // print file found to screen
      cout << "file found: " << str << "\n"; cout.flush();
    }
  }
  cout << endl;
  closedir(pdir);
}

//####################################################
//####################################################
