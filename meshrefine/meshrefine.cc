#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <vector>
#include <algorithm>
#include <ext/hash_map>

#include "classes.cc"
#include "functions.cc"

using std::vector;

int main(int argc,char *argv[]){

	if (argc != 3)
	{
        printf("\nSyntax: meshrefine input_file threshold\n\n");
        printf("Detail: Polygons are subdivided so that no edge is longer than threshold.\n");
        printf("        New mesh is written to stdout.\n\n");
		return 1;
	}

	////////// declare variables /////////
	char *eptr;
	bool flag=true;
	double threshold = strtod(argv[2],&eptr)*strtod(argv[2],&eptr);
    hashtable_t hme;
    hashtable_v hmv;
    std::vector<Edge*> e;
	int max_verts,max_faces;

	// vertex and face linked lists
	void_list *vlh,*flh;
	vlh=flh=NULL;

	////////// get data /////////
    getData(argv[1],flh,vlh);
    flh=addPrevious(flh);
    vlh=addPrevious(vlh);
    max_verts=maxVert(vlh);
    max_faces=maxFace(flh);
    fprintf(stderr,"main: original max_verts %i\n",max_verts);

    // map of pointers to vertices
    buildVertMap(vlh,hmv);
    addPointersToFaces(flh,hmv);

	int j=0;
	while (flag) {
        fprintf(stderr,"iteration %i: max vertex %i max face %i\n",j++,max_verts,max_faces);
        ///// build edge list ////
        clearEdges(e,hme);
        clearFaceEdges(flh);
        getEdges(flh,hme,e);
        // find adjacencies
        clearVertexAdjacencies(vlh);
        findVertexAdjacencies(flh,e);
		// process mesh
		flag = refineMesh(e,flh,vlh,threshold,max_faces,max_verts);
	}
	// print new mesh to stdout
	printMesh(vlh,flh);
	// clean up faces and vertices
	cleanup2(flh,vlh);
	return 0;
}
