#include <stdio.h>
#include <stdlib.h>
#include <string.h>

class void_list
{
public:
  void_list *previous;
  void_list *next;
  void *data;
};

class Vertex
{
public:
  double x,y,z;
  int index;
  Vertex(char *triplet);
};

Vertex::Vertex(char *triplet)
{
	char val[80];
	char *eptr;
	int i;

	// get past 'Vertex'
  	while (strchr("Vertx",*triplet)!=NULL) {triplet++;}

	// grab vertex index
	while (strchr(" \t,",*triplet)!=NULL) { triplet++; }
	i=0;
	while (strchr("0123456789+-eE.",*triplet)!=NULL)
	{
		val[i++] = *triplet++;
	}
	val[i]=0;
	index = (int) strtod(val,&eptr);
	if (val==eptr)
	{
		index=0;
		x=y=z=0;
		printf("Error in reading vertex index\n");
		return;
	}

	// grab x coord
	while (strchr(" \t,",*triplet)!=NULL) { triplet++; }
	i=0;
	while (strchr("0123456789+-eE.",*triplet)!=NULL)
	{
		val[i++] = *triplet++;
	}
	val[i]=0;
	x = strtod(val,&eptr);
	if (val==eptr)
	{
		x=y=z=0;
		printf("Error in reading vertex\n");
		return;
	}

	// grab y coord
	while (strchr(" \t,",*triplet)!=NULL) triplet++;
	i=0;
	while (strchr("0123456789+-eE.",*triplet))
	{
		val[i++] = *triplet++;
	}
	val[i]=0;
	y = strtod(val,&eptr);
	if (val==eptr)
	{
		x=y=z=0;
		printf("Error in reading vertex\n");
		return;
	}

	// grab z coord
	while (strchr(" \t,",*triplet)!=NULL) triplet++;
	i=0;
	while (strchr("0123456789+-eE.",*triplet))
	{
		val[i++] = *triplet++;
	}
	val[i]=0;
	z = strtod(val,&eptr);
	if (val==eptr)
	{
		x=y=z=0;
		printf("Error in reading vertex\n");
		return;
	}
}

class Face
{
public:
	int index;	// Face index
	int v1,v2,v3;	// vertex indices
	Face(char *triplet);
};

Face::Face(char *triplet)
{
	char val[80];
	char *eptr;
	int i;

	// get past 'Face'
  	while (strchr("Face",*triplet)!=NULL) {triplet++;}

	// grab Face index
	while (strchr(" \t,",*triplet)!=NULL) { triplet++; }
	i=0;
	while (strchr("0123456789+-eE.",*triplet)!=NULL)
	{
		val[i++] = *triplet++;
	}
	val[i]=0;
	index = (int) strtod(val,&eptr);
	if (val==eptr)
	{
		index=0;
		v1=v2=v3=0;
		printf("Error in reading face index\n");
		return;
	}

	// grab first vertex index
	while (strchr(" \t,",*triplet)!=NULL) { triplet++; }
	i=0;
	while (strchr("0123456789+-eE.",*triplet)!=NULL)
	{
		val[i++] = *triplet++;
	}
	val[i]=0;
	v1 = (int) strtod(val,&eptr);
	if (val==eptr)
	{
		v1=v2=v3=0;
		printf("Error in reading vertex index\n");
		return;
	}

	// grab second vertex index
	while (strchr(" \t,",*triplet)!=NULL) triplet++;
	i=0;
	while (strchr("0123456789+-eE.",*triplet))
	{
		val[i++] = *triplet++;
	}
	val[i]=0;
	v2 = (int) strtod(val,&eptr);
	if (val==eptr)
	{
		v1=v2=v3=0;
		printf("Error in reading vertex index\n");
		return;
	}

	// grab third vertex index
	while (strchr(" \t,",*triplet)!=NULL) triplet++;
	i=0;
	while (strchr("0123456789+-eE.",*triplet))
	{
		val[i++] = *triplet++;
	}
	val[i]=0;
	v3 = (int) strtod(val,&eptr);
	if (val==eptr)
	{
		v1=v2=v3=0;
		printf("Error in reading vertex index\n");
		return;
	}
}

int compare (const void* a, const void* b ) {

  return ( *((int*)b+1) - *((int*)a+1) );

}

int main(int argc,char *argv[]){

	if (argc != 2)
	{
		printf("\nSyntax: mesh_separate input_file\n\n");
		printf("Detail: The input mesh file is scanned for separate closed objects.\n");
		printf("Output: One file per separate object with the form input_file.#.mesh.\n\n");
		return 1;
	}

	fprintf(stderr,"\nInput mesh file is assumed to have vertex and\n");
	fprintf(stderr,"faces with sequentially increasing indeces. Run\n");
	fprintf(stderr,"the mesh input file through mesh_renumber first if\n");
	fprintf(stderr,"this criterion is not met.\n");

	////////// declare variables /////////

	char *infile;
	char outfile[128];
	char line[2048];
	char *str;
	FILE *F;
	void_list *q,*p,*prev;
	bool flag;
	int num_vertices, num_faces, i, j;

	// vertex linked lists for first and second file
	void_list *vlh,*vl;
	Vertex *v;
	vlh = NULL;

	// face linked lists for first and second file
	void_list *flh,*fl;
	Face *f;
	flh = NULL;

	// pointers to start of vertex and face linked lists for first and second file
	void_list *vend;
	void_list *fend;
	vend = NULL;
	fend = NULL;

	////////// get data /////////

	// open first file
	infile = argv[1];
	F = fopen(infile,"r");
	if (!F) {
	printf("Couldn't open input file %s\n",infile);
	return 1;
	}

	// for every line in first file
	num_vertices = 0;
	num_faces = 0;
	for (str=fgets(line,2048,F) ; str!=NULL ; str=fgets(line,2048,F)) {

		// skip leading whitespace
		while (strchr(" \t,",*str)!=NULL) { str++;}

		// if first character is V for Vertex, add new linked list class instance
		if (strchr("V",*str)!=NULL){
			vl = new void_list();
			vl->next = vlh;
			v = new Vertex(str);
			vl->data = (void*)v;
			vlh = vl;
			num_vertices++;
		} 
		// if first character is F for Face, add new linked list class instance
		else if (strchr("F",*str)!=NULL){
			fl = new void_list();
			fl->next = flh;
			f = new Face(str);
			fl->data = (void*)f;
			flh = fl;
			num_faces++;
		}
	}
	fclose(F);

	// go through linked list backwards and add previous pointers
	// for each vertex
	prev = NULL;
	for (p=vlh;p!=NULL;p=p->next) {
		p->previous = prev;
		prev = p;
		if (p->next==NULL) break;
	}

	// for each face
	prev = NULL;
	for (p=flh;p!=NULL;p=p->next) {
		p->previous = prev;
		prev = p;
		if (p->next==NULL) break;
	}

    // set end pointers
    for (q=vlh;q!=NULL;q=q->next) { vend=q; if (q->next==NULL) break; }
    for (q=flh;q!=NULL;q=q->next) { fend=q; if (q->next==NULL) break; }


	///// create and initialize faces array /////
	int faces[num_faces];
	for (i=1;i<num_faces+1;i++) {
		faces[i] = 0;
	}

	///// create and initialize vertex array /////
	int vertex[num_vertices];
	for (i=1;i<num_vertices+1;i++) {
		vertex[i] = 0;
	}

	///// create vertices array /////
	void_list* vertices[num_vertices];
	void_list *index;
	for (i=1;i<num_vertices+1;i++) {
		index = new void_list();
		index->next = NULL;
		index->data = (void*)NULL;
		vertices[i] = index;
	}

	///// load vertices array /////
	// for each face
	for (q=flh;q!=NULL;q=q->next) {

		// load face index into v1
		index = new void_list();
		index->next = vertices[(int)((Face*)q->data)->v1];
		index->data = (void*)((Face*)q->data)->index;
		vertices[(int)((Face*)q->data)->v1] = index;


		// load face index into v2
		index = new void_list();
		index->next = vertices[(int)((Face*)q->data)->v2];
		index->data = (void*)((Face*)q->data)->index;
		vertices[(int)((Face*)q->data)->v2] = index;

		// load face index into v3
		index = new void_list();
		index->next = vertices[(int)((Face*)q->data)->v3];
		index->data = (void*)((Face*)q->data)->index;
		vertices[(int)((Face*)q->data)->v3] = index;

		if (q->next==NULL) break;
	}


	///// coarsely separate closed surfaces /////
	int next_group = 1;
	int current_group;
	// for each vertex
	for (i=1;i<num_vertices+1;i++) {

		current_group = 0;
		// for each associated face
		for (q=vertices[i];q!=NULL;q=q->next) {
			// look up face group
			// if not last link in list
			if (q->next!=NULL) {
				if (faces[(int)q->data]){
					current_group = faces[(int)q->data];
				}
			}
			if (q->next==NULL) break;
		}

		// if a group was found
		if (current_group) {
			// set all faces to current_group
			// for each associated face
			for (q=vertices[i];q!=NULL;q=q->next) {
				// if not last link in list
				if (q->next!=NULL) {
					faces[(int)q->data] = current_group;
				}
				if (q->next==NULL) break;
			}
		} else {
		// no group was found
			// set all faces to next_group
			// for each associated face
			for (q=vertices[i];q!=NULL;q=q->next) {
				// if not last link in list
				if (q->next!=NULL) {
					faces[(int)q->data] = next_group;
				}
				if (q->next==NULL) break;
			}
			// increment next_group
			next_group++;
		}

	}

	///// amalgamate groups /////
	bool amalgamate = true;
	int smallest_group;
	bool swap;
	while(amalgamate) {
		amalgamate = false;
		// for each vertex
		for (i=1;i<num_vertices+1;i++) {
	
			// for each associated face
			swap = false;
			q=vertices[i];
			smallest_group = faces[(int)q->data];
			for (q=vertices[i];q!=NULL;q=q->next) {
				// look up face group
				// if not last link in list
				if (q->next!=NULL) {
					if (faces[(int)q->data]!=smallest_group){ swap=true;amalgamate = true;}
					if (faces[(int)q->data]<smallest_group){smallest_group=faces[(int)q->data];}
				}
				if (q->next==NULL) break;
			}
	
			// if two or more groups were found
			if (swap) {
				// set all faces to smallest_group
				// for each associated face
				for (q=vertices[i];q!=NULL;q=q->next) {
					// if not last link in list
					if (q->next!=NULL) {
						// if not smallest_group
						if (faces[(int)q->data]!=smallest_group){
							//swap all occurences of group in faces array with smallest_group
							for (j=1;j<num_faces+1;j++) {
								if(faces[j]==faces[(int)q->data]) { faces[j]=smallest_group; }
							}
						}
					}
					if (q->next==NULL) break;
				}
			}
	
		}
	}



	// for each vertex
	for (i=1;i<num_vertices+1;i++) {
		// for first associated face
		q=vertices[i];
		// set vertex group to face group
		vertex[i] = faces[(int)q->data];
	}

	///// count number of faces in each group /////
	int max_groups = 65536;
	int face_count[max_groups][2];
	int num_groups = 0;
	int loc;
	//for each face
	for (i=1;i<num_faces+1;i++) {
		//search for group in face_count
		flag = false;
		for (j=0;j<num_groups;j++) {
			// if face group found in face_count record location
			if(face_count[j][0]==faces[i]){loc = j;flag = true;}
		}
		// if face group found in face_count, increment count
		if (flag) {face_count[loc][1]++;}
		else {
			//add face group to face_count
			face_count[num_groups][0] = faces[i];
			face_count[num_groups][1] = 1;
			//check if room for more groups
			if (num_groups+1 < max_groups) {num_groups++;}
			else {printf("Increase size of max_groups (currently = %i)\n",max_groups);}
		}
	}

	///// sort face_count /////
    qsort(face_count,num_groups,2*sizeof(int),compare);

	///// output groups /////
	// for each group
	for (j=0;j<num_groups;j++) {

	    // open file
            char prefix[1024];
            int L = strlen(infile)-5;
            strncpy (prefix,infile,L);
            prefix[L]='\0';
	    //sprintf(outfile,"%s.%i",infile,j);
	    sprintf(outfile,"%s_%i.mesh",prefix,j);
	    printf("\n%s written (%i faces).\n\n",outfile,face_count[j][1]);
	    F = fopen(outfile,"w");
	    if (!F) {
	        printf("Couldn't open output file %s\n",outfile);
	        return 1;
	    }

		// for each vertex
	    for (p=vend;p!=NULL;p=p->previous) {
			// if vertex is in current group
			if(vertex[((Vertex*)p->data)->index]==face_count[j][0]) {
		        sprintf(line,"Vertex %i  %.15g %.15g %.15g\n",
		                ((Vertex*)p->data)->index,
		                ((Vertex*)p->data)->x,
		                ((Vertex*)p->data)->y,
		                ((Vertex*)p->data)->z);
		        fputs(line,F);
			}
			if (p->previous==NULL) break;
	    }

		// for each face
	    for (p=fend;p!=NULL;p=p->previous) {
			// if face is in current group
			if(faces[((Face*)p->data)->index]==face_count[j][0]) {
		        sprintf(line,"Face %i  %i %i %i\n",
		                ((Face*)p->data)->index,
		                ((Face*)p->data)->v1,
		                ((Face*)p->data)->v2,
		                ((Face*)p->data)->v3);
		        fputs(line,F);
			}
			if (p->previous==NULL) break;
	    }
		fclose(F);
		
	}

	return 0;
}
