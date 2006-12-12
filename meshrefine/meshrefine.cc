#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "classes.cc"
#include "functions.cc"

int main(int argc,char *argv[]){

	if (argc != 3)
	{
        printf("\nSyntax: meshrefine input_file threshold\n\n");
        printf("Detail: Polygons are subdivided so that no edge is longer than threshold.\n");
        printf("        New mesh is written to stdout.\n\n");
		return 1;
	}

	////////// declare variables /////////
	void_list *q,*qq,*p,*pp,*prev;
	char *eptr;
	int num_verts,num_faces,num_edges;
	EdgeBlock *eb;
	bool flag=true;
	double threshold = strtod(argv[2],&eptr)*strtod(argv[2],&eptr);

	// vertex and face linked lists
	void_list *vlh,*flh;
	vlh=flh=NULL;

	////////// get data /////////
	getData(argv[1],flh,vlh);

	flh=addPrevious(flh);
	vlh=addPrevious(vlh);

	int j=0;
	while (flag) {
		// process data
		num_verts=maxVert(vlh);
		num_faces=maxFace(flh);

		fprintf(stderr,"iteration %i: max vertex %i max face %i\n",j++,num_verts,num_faces);

		// array of pointers to vertices
		Vertex **vert_array;
		buildVertArray(vlh,vert_array,num_verts);

		///// build edge list ////
		num_edges=num_faces+num_verts-2+1000;
		eb = new EdgeBlock(num_edges);
		eb->ht = new HashTable(eb->n);
		getEdges(flh,eb,num_faces,vert_array);

		// process mesh
		flag = refineMesh(eb,flh,vlh,threshold,vert_array,num_faces,num_verts);

		// clean up edges
		cleanup1(eb,vert_array);

		// reset face flags and pointers
		resetFaces(flh);
	}

	// print new mesh to stdout
	printMesh(vlh,flh);

	// clean up faces and vertices
	cleanup2(flh,vlh);

	return 0;
}
