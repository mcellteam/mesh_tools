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

class Edge
{
public:
  int v1,v2; // vertex indices
  int c12,c21;  // count of times edge was traversed from v1 to v2 and v2 to v1
  void_list *f;
  bool crumb;
};

class Vertex
{
public:
  double x,y,z;
  int index;
  bool checked;
  Vertex(char *triplet);
};

Vertex::Vertex(char *triplet)
{
	char val[80];
	char *eptr;
	int i;

	checked = false;

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

int main(int argc,char *argv[]){

	if (argc != 2)
	{
		printf("\nSyntax: mesh_renumber input_file\n\n");
		printf("Description: Renumbers vertex and face indices of mesh file in sequentially increasing fashion. Writes output to stdout.\n\n");
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

	////////// create arrays //////////

	// find max verterx
	int max_vertex=0;
	for (p=vlh;p!=NULL;p=p->next) {
		if(((Vertex*)p->data)->index > max_vertex){max_vertex = ((Vertex*)p->data)->index;}
		if (p->next==NULL) break;
	}
	max_vertex++;

	// find max face
	int max_face=0;
	for (p=flh;p!=NULL;p=p->next) {
		if(((Face*)p->data)->index > max_face){max_face = ((Face*)p->data)->index;}
		if (p->next==NULL) break;
	}
	max_face++;

	// arrays
	int *verts,*faces;
	verts = new int[max_vertex];
	faces = new int[max_face];

	// load arrays
	int n = 1;
	for (p=vend;p!=NULL;p=p->previous) {
		verts[((Vertex*)p->data)->index] = n++;
		if (p->previous==NULL) break;
	}
	n = 1;
	for (p=fend;p!=NULL;p=p->previous) {
		faces[((Face*)p->data)->index] = n++;
		if (p->previous==NULL) break;
	}

	////////// fix linked lists //////////
	for (p=vlh;p!=NULL;p=p->next) {
		((Vertex*)p->data)->index = verts[((Vertex*)p->data)->index];
		if (p->next==NULL) break;
	}
	
	for (p=flh;p!=NULL;p=p->next) {
		((Face*)p->data)->index = faces[((Face*)p->data)->index];
		((Face*)p->data)->v1 = verts[((Face*)p->data)->v1];
		((Face*)p->data)->v2 = verts[((Face*)p->data)->v2];
		((Face*)p->data)->v3 = verts[((Face*)p->data)->v3];
		if (p->next==NULL) break;
	}
	
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

	return 0;
}
