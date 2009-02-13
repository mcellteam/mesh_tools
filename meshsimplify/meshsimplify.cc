#include <stdio.h>
#include <stdlib.h>
#include <string>

#include <vector>
#include <algorithm>
#include <ext/hash_map>

#define MAX_ITER 20

#include "classes.cc"
#include "functions.cc"

using std::vector;


int main(int argc,char *argv[]){

	if (argc != 3)
	{
        printf("\nSyntax: meshsimplify input_file threshold\n\n");
        printf("Detail: Vertices are merged so that no edge is shorter than threshold.\n");
        printf("        New mesh is written to stdout.\n\n");
		return 1;
	}

	////////// declare variables /////////
	char *eptr;
	double threshold = strtod(argv[2],&eptr)*strtod(argv[2],&eptr);
	hashtable_t hme;
	hashtable_v hmv;
	std::vector<Edge*> e;
//	vector<Vertex*> vert_array;
	Edge *ee;
	bool flag=true;
	int max_verts;

	// vertex and face linked lists
	void_list *vlh,*flh;
	vlh=flh=NULL;

	////////// get data /////////
	getData(argv[1],flh,vlh);
	flh=addPrevious(flh);
	vlh=addPrevious(vlh);
	max_verts=maxVert(vlh);
	fprintf(stderr,"main: original max_verts %i\n",max_verts);

	// map of pointers to vertices
	buildVertMap(vlh,hmv);
	addPointersToFaces(flh,hmv);

	while (flag) {
		///// build edge list ////
		clearEdges(e,hme);
		clearFaceEdges(flh);
		fprintf(stderr,"main: getting edges...");
		getEdges(flh,hme,e);
		fprintf(stderr,"complete.\n");
		// find adjacencies
		clearVertexAdjacencies(vlh);
		findVertexAdjacencies(flh,e);
		// find bad edge
		fprintf(stderr,"main: find bad edge...");
		ee=findBadEdge(e,threshold);
		fprintf(stderr,"complete.\n");
		if(ee!=NULL) {simplifyMesh(ee,flh,vlh,max_verts);flag=true;}
		else {flag=false;}
		fprintf(stderr,"main: max_verts %i\n",max_verts);

//		fprintf(stderr,"main: postcheck vertices...");
//		check(vlh);
//		fprintf(stderr,"complete.\n");
	}

	printMesh(vlh,flh);

	// clean up faces and vertices
	cleanup2(flh,vlh);
	clearEdges(e,hme);

	return 0;
}
