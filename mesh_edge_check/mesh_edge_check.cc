#include <stdio.h>
#include <stdlib.h>
#include <string>

#include <dirent.h>
#include <errno.h>

#include "classes.cc"
#include "functions.cc"

using namespace std;

int main(int argc,char *argv[]){

	if (argc != 3)
	{
		printf("\nSyntax: mesh_edge_check input_file|4all print_flag\n\n");
		printf("Description: Reports all polygon edges ");
		printf("that are not shared by exactly two polygons.\n");
		printf("If 4all, then all mesh files in current directory are analyzed.\n");
		printf("If print_flag=1, progress is reported to stdout.\n");
		printf("If print_flag=0, only a summary is reported to stdout.\n\n");
		return 1;
	}

	////////// declare variables /////////
	void_list *q,*files;
	files=NULL;
	int num_verts,num_faces,num_edges;
	EdgeBlock *eb;
	char *eptr;
	int print_flag = (int)strtod(argv[2],&eptr);

	// vertex and face linked lists
	void_list *vlh,*flh;
	vlh=flh=NULL;

	// pointers to start of vertex and face linked lists 
	void_list *vend,*fend;

	// get files to analyze
	if(!strcmp(argv[1],"4all")){
		files=scanDir(files);
	} else {
		files = new void_list();
		files->next=NULL;
		files->data=(void*)argv[1];
	}

	// for all mesh files in current directory
	for(q=files;q!=NULL;q=q->next){
		printf("\nanalyzing %s...\n",(char*)q->data);
		////////// get data /////////
		getData((char*)q->data,flh,vlh);
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
		flh=vlh=NULL;
		eb=NULL;
		printf("complete\n\n");
	}
	return 0;
}
