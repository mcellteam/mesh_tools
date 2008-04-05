#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <iostream>
#include <fstream>

#include <vector>
#include <algorithm>
#include <ext/hash_map>
#include <ext/hash_set>
#include <map>
#include <float.h>

using std::cout;
using std::endl;

#include "classes.cc"
#include "functions.cc"

int main(int argc,char *argv[]){

	if (argc != 3)
	{
        printf("\nSyntax: meshrefine input_file [threshold]\n\n");
        printf("Detail: Polygons are subdivided so that");
        printf(" no edge is longer than threshold.\n");
        printf("If threshold is negative each edge is bisected");
        printf(" once if edge length\n is longer than");
        printf(" absolute value of threshold.\n");
        printf("New mesh is written to stdout.\n\n");
        printf("Input mesh must be fully closed, consistently ");
		printf("oriented with vertex and face indices sequentially numbered.\n\n");
		return 1;
	}
	char *eptr;
	double threshold = strtod(argv[2],&eptr);

	// build object
	Object o(argv[1]);

	int j=0;
	do {
		fprintf(stderr,"iteration %i: ",j++);
		fprintf(stderr,"max vertex %i max face %i\n",o.v.size(),o.f.size());
		fprintf(stderr,"Building edges..................");
		fflush(stderr);
		o.createEdges();
		fprintf(stderr,"complete.\n");
		fflush(stderr);
	} while(refineMesh(o,threshold)==true);
	// print new mesh to stdout
	//printMesh(vlh,flh);
	printMesh(o);
	// clean up faces and vertices
//	cleanup2(flh,vlh);
	return 0;
}
