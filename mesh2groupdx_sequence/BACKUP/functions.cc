//####################################################
//####################################################

class Cabinet{
public:
	std::string indir,outdir;
        std::ofstream F,M,I;
	std::vector<std::string> files;	// array of input file names
	int num_files;	// number of input files found
	Cabinet(void);
	void scanDir(void);
	void parse(int,char**,std::string);
	void mesh2mcell(int);
	bool endswith(std::string const &, std::string const &);
	std::string strip_tail(std::string const &,uint);
	void updateMain(int);
	void updateInc(int);
	void finalizeMain(void);
	void finalizeScript(void);
	void close(void);
	void openFiles(void);
        char *FILENAME_FORMAT;
        char *STARTING_ITERATION;
        char *ENDING_ITERATION;
};


void Cabinet::finalizeScript(void){
        std::string START(STARTING_ITERATION);
        std::string END(ENDING_ITERATION);
	// prep main
	std::string message = "";
	message = message +
			"cd "+ outdir+ "\n"+
//			"echo \"running mcell\"\n"+
//			"~/bin/mcell3.1-amd64 main.mdl\n"+
//			"~/bin/mcell3 main.mdl\n"+
                        "echo \"initializing sequence.mdl\"\n"+
                        "echo \"i="+START+"\" > sequence.mdl\n"+
                        "a="+START+"\n"+
                        "b="+END+"\n"+
                        "while [ \"$a\" -lt \"$b\" ]\n"+
                        "do\n"+
                        "  echo \"running mcell $a\"\n"+
                        "  ~/bin/mcell3 main.mdl\n"+
                        "  let \"a+=1\"\n"+
                        "  done\n"+
			"echo \"cleaning up\"\n"+
			"#rm *.mdl"+
//			SCRIPT+
//			" n_ME-n_SSV.dat\n";
			"\n";
	F << message;
}

void Cabinet::finalizeMain(void){
	// prep main
	std::string message = "";
	message = message +
                        "A OBJECT A {}\n"+
			"}\n"+
			"\n"+
/*			"VIZ_DATA_OUTPUT {\n"+
			"  MODE=DX\n"+
			"  ITERATION_FRAME_DATA {\n"+
			"    SURFACE_POSITIONS = [0]\n"+
			"    SURFACE_STATES = [0]\n"+
			"  }\n"+
			"  STATE_VALUES {\n"+
			"      a = 1\n"+
			"  }\n"+
			"  OBJECT_FILE_DESIGNATORS {\n"+
			"     a = \"a\"\n"+
			"  }\n"+
			"}\n";
*/
			"VIZ_OUTPUT\n"+
			"{\n"+
			"  MODE=DREAMM_V3\n"+
			"  FILENAME = \"./a\"\n"+
			"  MESHES\n"+
			"  {\n"+
			"   NAME_LIST {ALL_MESHES }\n"+
			"   ITERATION_NUMBERS { ALL_DATA @ ALL_ITERATIONS}\n"+
			"  }\n"+
			"}\n";
	M << message;
}

void Cabinet::mesh2mcell(int count){
	// create paths
	std::string name = strip_tail(files[count],SUFFIX.size());
	std::string in = indir + name + SUFFIX;
	std::string out = outdir + name + ".mdl";
	// clear file
//    std::ofstream T(out.c_str());
//	if(T.is_open()==false){
//    	fprintf(stderr,"\nCouldn't open output file %s\n",out.c_str());
//		exit(0);
//	}
	// prep main
	F << "echo \"converting "
	<< in << " to mdl format\"\n";
	F << "mesh2mcell A "
	<< in << " > "
	<< out << "\n";
}

//void Cabinet::updateMain(int count){
//	std::string name = strip_tail(files[count],SUFFIX.size());
//	M << name << " OBJECT " << name << " {}\n";
//}

//void Cabinet::updateInc(int count){
//	std::string name = strip_tail(files[count],SUFFIX.size());
////	I << "INCLUDE_FILE=\"" << outdir << name << ".mdl\"\n";
//	I << "INCLUDE_FILE=\"" << name << ".mdl\"\n";
//}

Cabinet::Cabinet(void){
	num_files=0;
}

void Cabinet::openFiles(void){
	// open script
	std::string str = outdir + SCRIPT;
	F.open(str.c_str());
	if(F.is_open()==false){
    	fprintf(stderr,"\nCouldn't open output file %s\n",str.c_str());
		exit(0);
	}
	// open main.mdl
	str = outdir + MAIN;
	M.open(str.c_str());
	if(M.is_open()==false){
    	fprintf(stderr,"\nCouldn't open output file %s\n",str.c_str());
		exit(0);
	}
	// prep main
        std::string END(ENDING_ITERATION);
        std::string FILENAME(FILENAME_FORMAT);
	std::string message = "";
	message = message +
			"ITERATIONS = "+END+"\n"+
			"TIME_STEP = 1e-6\n"+
			"EFFECTOR_GRID_DENSITY = 15000  \n"+
			"\n"+
			"PARTITION_X = [[0.00001 TO 10000.00001 STEP 2000.00000001]]\n"+
			"PARTITION_Y = [[0.00001 TO 10000.00001 STEP 2000.00000001]]\n"+
			"PARTITION_Z = [[0.00001 TO 10000.00001 STEP 2000.00000001]]\n"+
			"\n"+
			"CHECKPOINT_INFILE = \"chkpt\"\n"+
			"CHECKPOINT_OUTFILE = \"chkpt\"\n"+
			"CHECKPOINT_ITERATIONS = 1\n"+
			"\n"+
			"/* defines i*/\n"+
			"INCLUDE_FILE=\"sequence.mdl\"\n"+
			"\n"+
			"/* read geometry using i*/\n"+
			"sprintf(geo_f,\""+FILENAME+"%g.mdl\",i)\n"+
			"INCLUDE_FILE=geo_f\n"+
			"\n"+
			"/* write incremented value of i*/\n"+
			"fname = fopen(\"sequence.mdl\",\"w\")\n"+
			"fprintf(fname,\"i=%g\\n\",i+1)\n"+
			"fclose(fname)\n"+
			"\n"+
			"INSTANTIATE a OBJECT {\n";
	M << message;
	// open include.mdl
//	str = outdir + INC;
//	I.open(str.c_str());
//	if(I.is_open()==false){
//    	fprintf(stderr,"\nCouldn't open output file %s\n",str.c_str());
//		exit(0);
//	}
}

std::string Cabinet::strip_tail(std::string const &str, uint len)
{
  if (str.size() < len)
    return str;

  return str.substr(0, str.size() - len);
}

void Cabinet::close(void){
    F.close();
    M.close();
    I.close();
}

void Cabinet::scanDir(void) {
  struct dirent *pent;			// pointer to dirent structure

  char filename[1024];
  strcpy(filename,indir.c_str());
  DIR *pdir = opendir(filename);	// pointer to a directory data structure
  if (!pdir) {printf("Error. Could not open %s.\n",filename);exit(1);}
  cout << "\nInput directory found: " << filename << endl << endl;
  while ((pent=readdir(pdir))){
    // copy char array to string
    std::string str = pent->d_name;
    // if found
    if (endswith(str,".mesh")==true) {
      // save filename
      files.push_back(str);
      // update index
      num_files++;
      // print file found to screen
      cout << "file found: " << str << "\n"; cout.flush();
    }
  }
  cout << endl;
  closedir(pdir);
}

bool Cabinet::endswith(std::string const &str, std::string const &suffix)
{
  if (str.size() < suffix.size())
    return false;

  return str.substr(str.size() - suffix.size()) == suffix;
}

void Cabinet::parse(int argc,char **argv,std::string message){

  // if no arguments passed
  if(argc!=9){
    cout << message << endl;
    exit(0);
  }

  int c;
  opterr=0;
  while((c=getopt(argc,argv,"hf:s:e:")) != -1){
    switch(c){ 
      case 'h':
        // print help message
        cout << message << endl;
        exit(0);
      case 'f':
        // mesh filename format
        FILENAME_FORMAT = optarg;
        break;
      case 's':
        STARTING_ITERATION = optarg;
        break;
      case 'e':
        ENDING_ITERATION = optarg;
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

//####################################################
//####################################################
