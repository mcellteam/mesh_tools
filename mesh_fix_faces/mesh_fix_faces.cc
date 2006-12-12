#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "classes.cc"
#include "functions.cc"

int main(int argc,char *argv[]){

    if (argc != 3)
    {
        printf("\nSyntax: mesh_fix_faces input_file print_flag\n\n");
        printf("Description: Identifies and fixes polygons with\n");
        printf("             reversed orientation relative to other\n");
        printf("             polygons, i.e. not all ccw or cw.\n");
        printf("             Writes output to stdout. print_flag = 1 to print\n");
        printf("             updated mesh to stdout. print_flag = 0 to not print.\n\n");
        return 1;
    }

	////////// declare variables /////////
	int num_verts,num_faces,num_edges;
	EdgeBlock *eb;
	char *eptr;

	// vertex and face linked lists
	void_list *vlh,*flh;
	vlh=flh=NULL;

	////////// get data /////////
	int print_flag = (int)strtod(argv[2],&eptr);

	getData(argv[1],flh,vlh);
	flh=addPrevious(flh);
	vlh=addPrevious(vlh);
	num_verts=maxVert(vlh);
	num_faces=maxFace(flh);

	///// build edge list ////
	num_edges=num_faces+num_verts-2+1000;
	eb = new EdgeBlock(num_edges);
	eb->ht = new HashTable(eb->n);
	getEdges(flh,eb,num_faces);

	// group faces
	groupFaces(eb,print_flag,flh);
	printFlipActivity(flh);
	countFaceGroups(flh,print_flag);

	// print mesh to stdout
	if (print_flag){ printMesh(vlh,flh); }

	cleanup(eb,flh,vlh);

	return 0;
}
