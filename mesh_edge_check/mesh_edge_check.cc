#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include</home/jkinney/mesh_tools/mesh_edge_check/classes.cc>
#include</home/jkinney/mesh_tools/mesh_edge_check/functions.cc>

int main(int argc,char *argv[]){

	if (argc != 3)
	{
		printf("\nSyntax: mesh_edge_check input_file print_flag\n\n");
		printf("Description: Reports all polygon edges ");
		printf("that are not shared by exactly two polygons.\n");
		printf("If print_flag=1, progress is reported to stdout.\n");
		printf("If print_flag=0, only a summary is reported to stdout.\n\n");
		return 1;
	}

	////////// declare variables /////////
	void_list *q,*qq,*p,*pp,*prev;
	int num_verts,num_faces,num_edges;
	EdgeBlock *eb;
	char *eptr;
	int print_flag = (int)strtod(argv[2],&eptr);

	// vertex and face linked lists
	void_list *vlh,*flh;
	vlh=flh=NULL;

	// pointers to start of vertex and face linked lists 
	void_list *vend,*fend;

	////////// get data /////////
	getData(argv[1],flh,vlh);
	flh=addPrevious(flh);
	vlh=addPrevious(vlh);
	num_verts=maxVert(vlh);
	num_faces=maxFace(flh);

	///// build edge list ////
	num_edges=num_faces+num_verts-2+1000;
	eb = new EdgeBlock(num_edges);
	eb->ht = new HashTable(eb->n);
	getEdges(flh,eb,num_faces,print_flag);

	// identify bad faces
	identifyBadFaces(eb,flh,vlh);

	cleanup(eb,flh,vlh);

	return 0;
}
