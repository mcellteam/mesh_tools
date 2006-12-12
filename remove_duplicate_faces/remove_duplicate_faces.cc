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

	double i;
	void_list *j,*k;

	j = *(void_list**)a;
	k = *(void_list**)b;

	i = ((Face*)j->data)->v1 - ((Face*)k->data)->v1;
	if (i<0) {	return -1;}
	else if (i>0) {	return 1;}
	else {
		i = ((Face*)j->data)->v2 - ((Face*)k->data)->v2;
		if (i<0) {	return -1;}
		else if (i>0) {	return 1;}
		else {
			i = ((Face*)j->data)->v3 - ((Face*)k->data)->v3;
			if (i<0) {	return -1;}
			else if (i>0) {	return 1;}
			else { return (0);}
		}
	}
}

void removeLink(void_list* L) {
    // and remove face from candidate face list
    if (L->previous!=NULL) {
        if (L->next!=NULL) {
            // if both previous and next exist
            (L->previous)->next = L->next;
            (L->next)->previous = L->previous;
        } else {
            // if previous exists and next does not
            (L->previous)->next = NULL;
        }
    } else {
        if (L->next!=NULL) {
            // if previous does not exist and next does
            (L->next)->previous = NULL;
        } // else { // if neither previous nor next exists }
    }
}


int main(int argc,char *argv[]){

	if (argc != 3)
	{
		printf("\nSyntax: remove_duplicate_faces input_file print_flag\n\n");
		printf("Description: Identifies and removes duplicate faces\n");
		printf("		and writes output to stdout. print_flag = 1 to print\n");
		printf("		updated mesh to stdout. print_flag = 0 to not print.\n\n");
		return 1;
	}

	////////// declare variables /////////

	char *infile;
	char line[2048];
	char *str;
	char *eptr;
	FILE *F;
	void_list *q,*qq,*p,*pp,*prev;

	// vertex and face linked lists
	void_list *vlh,*vl;
	Vertex *v;
	vlh = NULL;
	void_list *flh,*fl;
	Face *f;
	flh = NULL;

	// pointers to start of vertex and face linked lists 
	void_list *vend;
	void_list *fend;

	////////// get data /////////

	int print_flag = (int)strtod(argv[2],&eptr);

	// open first file
	infile = argv[1];
	F = fopen(infile,"r");
	if (!F) {
	printf("Couldn't open input file %s\n",infile);
	return 1;
	}

	// for every line in first file
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
		} 
		// if first character is F for Face, add new linked list class instance
		else if (strchr("F",*str)!=NULL){
			fl = new void_list();
			fl->next = flh;
			f = new Face(str);
			fl->data = (void*)f;
			flh = fl;
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

	// find max vertex index
    int max_vertex = 0;
    for (p=vlh;p!=NULL;p=p->next) {
        if (((Vertex*)p->data)->index > max_vertex) {max_vertex = ((Vertex*)p->data)->index;}
        if (p->next==NULL) break;
    }

	// create array of pointers to void list
	int num = 3*max_vertex;
	void_list* face_array[num];

	// null face_array pointers;
	int i;
	for (i=0;i<num;i++) {
		face_array[i] = NULL;
	}

	///// group faces by sum of vertex indices /////
	// for each face
	int sum;
	void_list *n,*nh;
	for (p=flh;p!=NULL;p=p->next) {
		// sum of vertex indices
		sum = ((Face*)p->data)->v1 + ((Face*)p->data)->v2 + ((Face*)p->data)->v3;
		// create new link in appropriate list
		n = new void_list;
		n->next = face_array[sum];
		n->data = (void*)p;
		face_array[sum] = n;

		if (p->next==NULL) break;
	}
	
	///// turn array into linked list /////
//	void_list *sort_h,*sort;
//	sort_h = NULL;
//	int i;
//	for (i=0;i<num_face;i++) {
//		sort = new void_list();
//		sort->next = sort_h;
//		sort->data = (void*)face_array[i];
//		sort_h = sort;
//	}

	// go through linked list backwards and add previous pointers
	// for each vertex
//	prev = NULL;
//	for (p=sort_h;p!=NULL;p=p->next) {
//		p->previous = prev;
//		prev = p;
//		if (p->next==NULL) break;
//	}

	///// remove duplicate faces from list////
	void_list *ptr,*ptr2;
	ptr = flh;
	// for each group in face_array
	for (i=0;i<num;i++) {

		// if pointer is not null
		if (face_array[i]!=NULL) {
			ptr2 = face_array[i];
			// for all combinations of links in list
			for (p=face_array[i];p!=NULL;p=p->next) {
				for (q=ptr2;q!=NULL;q=q->next) {

					if ( ((Face*)((void_list*)p->data)->data)->index != ((Face*)((void_list*)q->data)->data)->index ) {

						// if two faces are identical
						if (
							(
							(((Face*)((void_list*)p->data)->data)->v1 == ((Face*)((void_list*)q->data)->data)->v1) ||
							(((Face*)((void_list*)p->data)->data)->v1 == ((Face*)((void_list*)q->data)->data)->v2) ||
							(((Face*)((void_list*)p->data)->data)->v1 == ((Face*)((void_list*)q->data)->data)->v3)
							) &&
							(
							(((Face*)((void_list*)p->data)->data)->v2 == ((Face*)((void_list*)q->data)->data)->v1) ||
							(((Face*)((void_list*)p->data)->data)->v2 == ((Face*)((void_list*)q->data)->data)->v2) ||
							(((Face*)((void_list*)p->data)->data)->v2 == ((Face*)((void_list*)q->data)->data)->v3)
							) &&
							(
							(((Face*)((void_list*)p->data)->data)->v3 == ((Face*)((void_list*)q->data)->data)->v1) ||
							(((Face*)((void_list*)p->data)->data)->v3 == ((Face*)((void_list*)q->data)->data)->v2) ||
							(((Face*)((void_list*)p->data)->data)->v3 == ((Face*)((void_list*)q->data)->data)->v3)
							) 
							) {
	
							if (!print_flag) {
								printf("Face %i %i %i %i equal to Face %i %i %i %i\n",
									((Face*)((void_list*)p->data)->data)->index,
									((Face*)((void_list*)p->data)->data)->v1,
									((Face*)((void_list*)p->data)->data)->v2,
									((Face*)((void_list*)p->data)->data)->v3,
									((Face*)((void_list*)q->data)->data)->index,
									((Face*)((void_list*)q->data)->data)->v1,
									((Face*)((void_list*)q->data)->data)->v2,
									((Face*)((void_list*)q->data)->data)->v3
									);
							}
	
							// remove face link from list
					        removeLink((void_list*)q->data);        
							// adjust list pointer
					        if ((void_list*)q->data==ptr) {ptr = ((void_list*)q->data)->next; }
						}
					}
				}
				// remove link from p
		        removeLink(p);        
				// adjust list pointer
		        if (p==ptr2) {ptr2 = p->next;}
			}
		}
	}
	// adjust pointers
	flh = ptr;

	////////// add pointers to ends of linked lists /////////
	// read out vertex linked list
	for (q=vlh;q!=NULL;q=q->next) {
		vend = q;
		if (q->next==NULL) break;
	}

	// read out face linked list
	for (q=flh;q!=NULL;q=q->next) {
		fend = q;
		if (q->next==NULL) break;
	}

if(print_flag) {
	////////// write data to stdout /////////
	// write out vertex linked list
	for (q=vend;q!=NULL;q=q->previous) {
      	printf("Vertex %i  %.15g %.15g %.15g\n",((Vertex*)q->data)->index
								,((Vertex*)q->data)->x
								,((Vertex*)q->data)->y
								,((Vertex*)q->data)->z);
		if (q->previous==NULL) break;
	}

	// write out face linked list
	for (q=fend;q!=NULL;q=q->previous) {
      	printf("Face %i  %i %i %i\n",((Face*)q->data)->index
								,((Face*)q->data)->v1
								,((Face*)q->data)->v2
								,((Face*)q->data)->v3);
		if (q->previous==NULL) break;
	}
}
	return 0;
}
