#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int distinguishable(double a,double b,double eps)
{
  double c;

  c=a-b;

  if (c<0) c=-c;
  if (a<0) a=-a;
  if (a<1) a=1;
  if (b<0) b=-b;

  if (b<a) eps*=a;
  else eps*=b;
  return (c>eps);
}

class void_list
{
public:
  void_list *previous;
  void_list *next;
  void *data;
};

class Swap
{
public:
  int index_old;
  int index_new;
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

	i = ((Vertex*)j->data)->x - ((Vertex*)k->data)->x;
	if (i<0) {	return -1;}
	else if (i>0) {	return 1;}
	else {
		i = ((Vertex*)j->data)->y - ((Vertex*)k->data)->y;
		if (i<0) {	return -1;}
		else if (i>0) {	return 1;}
		else {
			i = ((Vertex*)j->data)->z - ((Vertex*)k->data)->z;
			if (i<0) {	return -1;}
			else if (i>0) {	return 1;}
			else { return (0);}
		}
	}
}


int main(int argc,char *argv[]){

	if (argc != 2)
	{
		printf("\nSyntax: remove_duplicates input_file\n\n");
		printf("Description: Identifies and removes duplicate vertices,\n");
		printf("		updates polygon vertex info, and writes\n");
		printf("		output to stdout.\n\n");
		return 1;
	}

	////////// declare variables /////////

	char *infile;
	char line[2048];
	char *str;
	char *eptr;
	FILE *F;
	void_list *q,*qq,*p,*pp,*prev;
	double epsilon = 1e-4;


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

	// array of pointers to vertex linked list
	int num_vert=0;

	////////// get data /////////

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
			num_vert++;
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

	// create array of pointers to vertex list
	void_list* vert_array[num_vert];
	num_vert=0;

	// go through linked list backwards and add previous pointers
	// for each vertex
	prev = NULL;
	for (p=vlh;p!=NULL;p=p->next) {
		p->previous = prev;
		prev = p;
		vert_array[num_vert++]=p;
		if (p->next==NULL) break;
	}

	// for each face
	prev = NULL;
	for (p=flh;p!=NULL;p=p->next) {
		p->previous = prev;
		prev = p;
		if (p->next==NULL) break;
	}

	////////// sort vertices by coordinates /////////
	qsort(vert_array,num_vert,sizeof(void_list*),compare);

	///// turn array into linked list /////
	void_list *sort_h,*sort;
	sort_h = NULL;
	int i;
	for (i=0;i<num_vert;i++) {
		sort = new void_list();
		sort->next = sort_h;
		sort->data = (void*)vert_array[i];
		sort_h = sort;
	}

	// go through linked list backwards and add previous pointers
	// for each vertex
	prev = NULL;
	for (p=sort_h;p!=NULL;p=p->next) {
		p->previous = prev;
		prev = p;
		if (p->next==NULL) break;
	}

	///// remove duplicate vertices from list////
	void_list *swap_h,*swap;
	void_list* j,*jp1;
	swap_h = NULL;
	Swap* s;
	// for each vertex in sorted linked list
	q=sort_h;
	while (q!=NULL) {

		j = (void_list*)q->data;
		jp1 = (void_list*)((q->next)->data);


		// if current vertex is identical to next vertex
		if (
			( !distinguishable(((Vertex*)j->data)->x,((Vertex*)jp1->data)->x,epsilon) &&
			!distinguishable(((Vertex*)j->data)->y,((Vertex*)jp1->data)->y,epsilon) &&
			!distinguishable(((Vertex*)j->data)->z,((Vertex*)jp1->data)->z,epsilon)
			) ){

			printf("%i %g %g %g indistinguishables from %i %g %g %g\n",
				((Vertex*)j->data)->index,
				((Vertex*)j->data)->x,
				((Vertex*)j->data)->y,
				((Vertex*)j->data)->z,
				((Vertex*)jp1->data)->index,
				((Vertex*)jp1->data)->x,
				((Vertex*)jp1->data)->y,
				((Vertex*)jp1->data)->z
					);
		}

		// if current vertex is identical to next vertex
		if (
			(((Vertex*)j->data)->x == ((Vertex*)jp1->data)->x) &&
			(((Vertex*)j->data)->y == ((Vertex*)jp1->data)->y) &&
			(((Vertex*)j->data)->z == ((Vertex*)jp1->data)->z)
			) {

//			printf("%i %g %g %g equal to %i %g %g %g\n",
//				((Vertex*)j->data)->index,
//				((Vertex*)j->data)->x,
//				((Vertex*)j->data)->y,
//				((Vertex*)j->data)->z,
//				((Vertex*)jp1->data)->index,
//				((Vertex*)jp1->data)->x,
//				((Vertex*)jp1->data)->y,
//				((Vertex*)jp1->data)->z
//					);

			// save vertex indices
			swap = new void_list();
			swap->next = swap_h;
			s = new Swap;
			s->index_old = ((Vertex*)j->data)->index;
			s->index_new = ((Vertex*)jp1->data)->index;
			swap->data = (void*)s;
			swap_h = swap;

			// remove link from list
			// choose first identical vertex to be deleted
			// use index of second identical vertex to replace instances of first
			if (j->previous!=NULL) {
				if (j->next!=NULL) {
					// if both previous and next exist
					(j->previous)->next = j->next;
					(j->next)->previous = j->previous;
				} else {
					// if previous exists and next does not
					(j->previous)->next = NULL;
				}
			} else {
				if (j->next!=NULL) {
					// if previous does not exist and next does
					(j->next)->previous = NULL;
				} // else { // if neither previous nor next exists }
			}
		}

		q=q->next;	// supposed to be before following if statement.
		if (q->next==NULL) break; // do not change!!
	}

if (0) {
		// search swap linked list
		for (p=swap_h;p!=NULL;p=p->next) {
			printf("%i is changed to %i\n",((Swap*)p->data)->index_old,((Swap*)p->data)->index_new);
			if (p->next==NULL) break;
		}
}

	///// fix faces /////
	// for each face in linked list
	for (q=flh;q!=NULL;q=q->next) {

		// search swap linked list
		for (p=swap_h;p!=NULL;p=p->next) {

      		if ( ((Face*)q->data)->v1 == ((Swap*)p->data)->index_old ) {
				((Face*)q->data)->v1 = ((Swap*)p->data)->index_new;
			}
      		if ( ((Face*)q->data)->v2 == ((Swap*)p->data)->index_old ) {
				((Face*)q->data)->v2 = ((Swap*)p->data)->index_new;
			}
      		if ( ((Face*)q->data)->v3 == ((Swap*)p->data)->index_old ) {
				((Face*)q->data)->v3 = ((Swap*)p->data)->index_new;
			}
		
			if (p->next==NULL) break;
		}
		if (q->next==NULL) break;
	}
	

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

if(0) {
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
