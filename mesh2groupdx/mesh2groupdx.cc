#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <dirent.h>
#include <sys/stat.h>
#include <unistd.h>

using std::cout;
using std::endl;

//####################################################
//####################################################

std::string BASH("/bin/bash");
std::string SCRIPT("mesh2groupdx.sh");
std::string MAIN("main.mdl");
std::string INC("include.mdl");
std::string SUFFIX(".mesh");

extern char **environ;

//####################################################
//####################################################

#include "functions.cc"

//####################################################
//###################### main  #######################
//####################################################

int main(int argc,char **argv){
	std::string message = "\n";
	message=message+
	"NAME\n"+
	"       mesh2groupdx - save multiple meshes as a single dx file\n"+
	"\nSYNOPSIS\n"+
//	"       mesh2groupdx DIRECTORY|FILE...\n"+
	"       mesh2groupdx INPUTDIR OUTPUTDIR\n"+
	"\nDESCRIPTION\n"+
	"       mesh2groupdx includes all input meshes into\n"+
	"       a single dx file for group visualization.\n\n"+
	"       Uses mesh2mcell and mcell3.\n\n"+
	"       Files created in OUTPUTDIR directory:\n"+
	"          1. An mcell3 version (.mdl) of each input\n"+
	"             mesh file is written to the current\n"+
	"             working directory.\n"+
	"          2. mcell3 files are created: {main.mdl,\n"+
	"             include.mdl}\n"+
	"          3. A bash script 'mesh2groupdx.sh' is\n"+
	"             written to orchestrate the process.\n"+
	"          4. A single dx file and multiple .bin\n"+
	"             files are written by mcell3.\n"+
	"       The mdl files created by mesh2groupdx\n"+
	"       as well as the script are deleted.\n"+
	"       Files in INPUTDIR directory ending in '.mesh'\n"+
	"       are gathered for inclusion.\n"+
	"\nOPTIONS\n"+
	"       None.\n"+
	"\nEXAMPLES\n"+
	"       mesh2groupdx ./input_meshes ./output_meshes\n"+
	"\nTODO\n"+
	"       Add sequence of filenames on command line as input.\n"+
	"       Add verbose option for file-by-file updates to STDOUT.\n"+
	"       Add option for specifying path to mcell executable.\n"+
	"\nJustin Kinney				2007/10/01\n";

	// instantiate cabinet
	Cabinet c;
	// parse command line 
	c.parse(argc,argv,message);
	// open output files
	c.openFiles();
	// scan folder
	c.scanDir();
	// for each file in folder
	for (int count=0;count<c.num_files;count++) {
		// wrtie mesh2mcell command
		c.mesh2mcell(count);
		// update main
		c.updateMain(count);
		// update include
		c.updateInc(count);
	}
	c.finalizeMain();
	c.finalizeScript();
	c.close();
	std::string cmd = BASH + " " + c.outdir + SCRIPT;
	system(cmd.c_str());
}
