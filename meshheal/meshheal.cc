#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <string>

#include "classes.cc"
#include "functions.cc"

int main(int argc,char *argv[]){

	if (argc != 4)
	{
		printf("\nSyntax: meshheal input_file epsilon print_flag\n\n");
		printf("Description: Fills holes in meshes caused by vertices separated\n");
		printf("             by small a  distance, epsilon, and therefore not being treated\n");
		printf("             as the same vertex. This script takes a user-defined epsilon\n");
		printf("             and identifies and removes duplicate vertices,\n");
		printf("             updates polygon vertex info, and writes\n");
		printf("             output to stdout. print_flag = 1 to print\n");
		printf("             updated mesh to stdout. print_flag = 0 to not print.\n\n");
		return 1;
	}

/*
This program works as follows. For each vertex on a border, i.e. open edge, 
compute the squared distances to neighboring vertices in 3D to identify the closest
vertex that is also on a border. Declare these two vertices to be duplicates.
*/

	////////// declare variables /////////

	char *eptr;
	int i,num_vfree,d_count,num_verts,num_faces,num_edges;
    EdgeBlock *eb;
	void_list** dist_array;
    Vertex **vert_array;

	int delay = 0;

	// linked lists
	void_list *flh,*vlh,*vfree,*vbad;
	flh=vlh=vfree=vbad=NULL;

	////////// get data /////////
	double epsilon = strtod(argv[2],&eptr);
	int print_flag = (int)strtod(argv[3],&eptr);
	getData(argv[1],flh,vlh);
    flh=addPrevious(flh);
    vlh=addPrevious(vlh);
    num_verts=maxVert(vlh);
	num_faces=maxFace(flh);

    // create array of pointers to vertex list
	initVertArray(vlh,vert_array,num_verts);

    ///// build edge list ////
    num_edges=num_faces+num_verts-2+10000;
    eb = new EdgeBlock(num_edges);
    eb->ht = new HashTable(eb->n);
    getEdges(flh,eb,num_faces,print_flag);

	// get max edge length squared
	double threshold = computeLongestEdge(eb,vert_array,epsilon);

	///// gather free vertices /////
	vfree = gatherFreeVertices(eb,num_verts,vert_array,print_flag);

	///// compute distances /////
	dist_array = computeDistances(eb,dist_array,vfree,flh,d_count,epsilon,print_flag,threshold);

	///// count number free vertices /////
	num_vfree = countFreeVertices(vfree,print_flag);

	///// merge vertices /////
	bool flag=true;
	while(flag){
		///// identify free vertex pair separated by the smallest distance /////
		i=findVerticesToMerge(dist_array,d_count,vfree,print_flag,epsilon);
		if (i==d_count){flag=false;}
		else if (sqrt(((Distance*)dist_array[i]->data)->d*epsilon*epsilon)>1e-9){flag=false;}
		else {
			///// remove second vertex in pair from vertex list /////
			vlh = fixVertices(vlh,(Distance*)dist_array[i]->data,print_flag,vbad);

			///// replace all instances of second vertex in face list with first vertex /////
			flh = fixFaces(flh,(Distance*)dist_array[i]->data,print_flag);
		
			///// delete edge list /////
			deleteEdges(eb);
		
		    ///// build edge list ////
		    num_edges=num_faces+num_verts-2+10000;
		    eb = new EdgeBlock(num_edges);
		    eb->ht = new HashTable(eb->n);
		    getEdges(flh,eb,num_faces,print_flag);
		
			///// delete free vertices /////
			deleteFreeVertices(vfree);
			vfree=NULL;
		
			///// gather free vertices /////
			vfree = gatherFreeVertices(eb,num_verts,vert_array,print_flag);
		
			///// count number free vertices /////
			num_vfree = countFreeVertices(vfree,print_flag);
			if(!num_vfree){flag=false;}
		}
	}
	// print to stdout
	if(print_flag){printVerticesFaces(vlh,flh);}

	///// cleanup /////
	cleanup(eb,vfree,vlh,flh,dist_array,d_count,vert_array,vbad);

	return 0;
}
