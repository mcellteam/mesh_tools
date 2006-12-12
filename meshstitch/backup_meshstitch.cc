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

class Vertex
{
public:
  double x,y,z;
  int index;
  Vertex(char *triplet);
};

class ExtraVertex
{
public:
  Vertex *v;
  Vertex *lateral1,*lateral2;
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

void_list* addPrevious(void_list* L) {
	// go through linked list backwards and add previous pointers
	void_list *p,*prev;
	prev = NULL;
	for (p=L;p!=NULL;p=p->next) {
		p->previous = prev;
		prev = p;
		if (p->next==NULL) break;
	}
	return L;
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

int vertexInFace(void_list *F,void_list *V){
	return ( ( ((Face*)F->data)->v1 == ((Vertex*)V->data)->index ) ||
			( ((Face*)F->data)->v2 == ((Vertex*)V->data)->index ) ||
			( ((Face*)F->data)->v3 == ((Vertex*)V->data)->index ) );
}

void_list * removeBadFaces(void_list *F,void_list *B) {
	void_list *p,*q,*ptr1,*ptr2;
	bool flag;
	ptr1 = F;
	ptr2 = B;
	// for each face
	for (p=F;p!=NULL;p=p->next) {
		// search bad face list
		flag = true;
		q=ptr2;
		while(q!=NULL && flag) {
			// if face is bad
    	  	if ( ((Face*)p->data)->index == ((Face*)q->data)->index ) {
				// remove face from face list
				removeLink(p);
				// adjust list pointer
				if (p==ptr1) { ptr1 = p->next; }
				// remove bad face from bad face list
				removeLink(q);
				// adjust list pointer
				if (q==ptr2) { ptr2 = q->next; }
				flag = false;
			}
			if (q->next==NULL) break;
			q=q->next;
		}
		if (p->next==NULL) break;
	}
	return ptr1;
}



void_list * removeBadVertices(void_list *V,void_list *B) {
	void_list *p,*q, *ptr1,*ptr2;
	bool flag;
	ptr1 = V;
	ptr2 = B;
	// for each vertex
	for (p=V;p!=NULL;p=p->next) {
		// search bad vertex list
		flag = true;
		q=ptr2;
		while(q!=NULL && flag) {
			// if vertex is bad
    	  	if ( ((Vertex*)p->data)->index == ((Vertex*)q->data)->index ) {
				// remove vertex from vertex list
				removeLink(p);
				// adjust list pointer
				if (p==ptr1) { ptr1 = p->next; }
				// remove bad vertex from bad vertex list
				removeLink(q);
				// adjust list pointer
				if (q==ptr2) { ptr2 = q->next; }
				flag = false;
			}
			if (q->next==NULL) break;
			q=q->next;
		}
		if (p->next==NULL) break;
	}
	return ptr1;
}


void_list* matchZ(void_list *L,int z_value){
	// build linked list of first vertex elements with z_value matching target
	void_list *q,*t,*th;
	th = NULL;
	for (q=L;q!=NULL;q=q->next) {
      	if (((Vertex*)q->data)->z == z_value) {
			t = new void_list();
			t->next = th;
			t->data = (Vertex*)q->data;
			th = t;
		}
		if (q->next==NULL) break;
	}
	return th;
}



void_list* candidateBadVertices(void_list *V, int val){
	// gather candidate bad vertices from first vertex linked list
	void_list *q,*t,*th;
	th = NULL;
	// for each vertex
	for (q=V;q!=NULL;q=q->next) {
		// if z value matches target
		if (((Vertex*)q->data)->z == val) {
			// store vertex 
			t = new void_list();
			t->next = th;
			t->data = (Vertex*)q->data;
			th = t;
		}
		if (q->next==NULL) break;
	}
	return th;
}

void_list* candidateBadFaces(void_list *F,void_list *V){
	///// find faces that have at least one vertex from cbv lists /////
	void_list *q,*p,*t,*th;
	th = NULL;
	bool flag;
	// for each face from first face linked list
	for (q=F;q!=NULL;q=q->next) {
		// search cbv1 list 
		flag = true;
		p=V;
		while(p!=NULL && flag) {
			// if candidate bad vertex is at least one of three face vertices
    	  	if (vertexInFace(q,p)) {
				// add face to candidate bad face list
				t = new void_list();
				t->next = th;
				t->data = (Face*)q->data;
				th = t;
				flag = false;
			}
			if (p->next==NULL) break;
			p=p->next;
		}
		if (q->next==NULL) break;
	}
	return th;
}

void_list* findBadFaces(void_list* CBF,void_list* V,void_list* BF) {
	// scan candidate bad faces searching for known bad vertices
	void_list *q,*p,*t,*th;
	th = BF;
	bool flag;
	// for each candidate bad face
	for (p=CBF;p!=NULL;p=p->next) {
		// search bad vertex vertex list
		flag = true;
		q=V;
		while(q!=NULL && flag) {
			// if bad vertex is at least one of three face vertices
    	  	if (vertexInFace(p,q)) {
				// add face to bad face list
				t = new void_list();
				t->next = th;
				t->data = (Face*)p->data;
				th = t;
				flag = false;
				// and remove face from candidate face list
				removeLink(p);
			}
			if (q->next==NULL) break;
			q=q->next;
		}
		if (p->next==NULL) break;
	}
	return th;
}

void_list * extraVertices(void_list *M,void_list *S,void_list *F) {
	// find extra vertices
	void_list *p,*q,*t,*th;
	ExtraVertex *EV;
	th=NULL;
	bool found;
	// for each target z vertices
	for (p=M;p!=NULL;p=p->next) {
		// check if vertex is shared
		found = false;
		q = S;
		while(q!=NULL && !found) {
			if ( ((Vertex*)p->data)->index == ((Vertex*)q->data)->index) { found = true;}
			if (q->next==NULL) break;
			q=q->next;
		}
		// if not shared
		if (!found) {
			// is vertex part of a bad face
			found = false;
			q=F;
			while(q!=NULL && !found) {
				if (vertexInFace(q,p)) { found = true;}
				if (q->next==NULL) break;
				q=q->next;
			}

			// if vertex has target z, is not shared, and is part of bad face
			// then it is an extra vertex
			if (found) {
				t = new void_list();
				EV = new ExtraVertex;
				EV->v = (Vertex*)p->data;
				EV->lateral1 = NULL;
				EV->lateral2 = NULL;
				t->next = th;
				t->data = (void*)EV;
				th = t;
			}
		}
		if (p->next==NULL) break;
	}
	return th;
}

void lateralVertices(void_list *EV,void_list *BF,void_list *S) {
	// find lateral vertices
	void_list *p,*q,*pp;
	bool found;
	// for each extra vertex
	for (p=EV;p!=NULL;p=p->next) {
	
		// for each bad face
		for (q=BF;q!=NULL;q=q->next) {
		
			// is extra vertex a member of bad face?
			if ( 
				((((ExtraVertex*)p->data)->v)->index == ((Face*)q->data)->v1) || 
				((((ExtraVertex*)p->data)->v)->index == ((Face*)q->data)->v2) || 
				((((ExtraVertex*)p->data)->v)->index == ((Face*)q->data)->v3) ){
				// find shared lateral vertex in bad face
				found = false;
				pp=S;
				while(pp!=NULL && !found) {
				if (vertexInFace(q,pp)) {found = true;}
					if (pp->next==NULL) break;
					if(!found) {pp=pp->next;}
				}
				// add shared lateral vertex to extra vertex
				if (found) {
					if( !((ExtraVertex*)p->data)->lateral1 ) 
						{((ExtraVertex*)p->data)->lateral1 = (Vertex*)pp->data;}
					else if( !((ExtraVertex*)p->data)->lateral2 ) 
						{((ExtraVertex*)p->data)->lateral2 = (Vertex*)pp->data;}
				}
			}
			if (q->next==NULL) break;
		}
		if (p->next==NULL) break;
	}
}


void convertLateral(void_list *V,void_list *EV,double epsilon) {
	// convert lateral vertices to other object
	void_list *p,*q;
	for (p=V;p!=NULL;p=p->next) {
		for (q=EV;q!=NULL;q=q->next) {

			if ( !distinguishable(((Vertex*)p->data)->x,(((ExtraVertex*)q->data)->lateral1)->x,epsilon) &&
				!distinguishable(((Vertex*)p->data)->y,(((ExtraVertex*)q->data)->lateral1)->y,epsilon) &&
				!distinguishable(((Vertex*)p->data)->z,(((ExtraVertex*)q->data)->lateral1)->z,epsilon)
				){((ExtraVertex*)q->data)->lateral1 = (Vertex*)p->data;}

			if ( !distinguishable(((Vertex*)p->data)->x,(((ExtraVertex*)q->data)->lateral2)->x,epsilon) &&
				!distinguishable(((Vertex*)p->data)->y,(((ExtraVertex*)q->data)->lateral2)->y,epsilon) &&
				!distinguishable(((Vertex*)p->data)->z,(((ExtraVertex*)q->data)->lateral2)->z,epsilon)
				){((ExtraVertex*)q->data)->lateral2 = (Vertex*)p->data;}
			if (q->next==NULL) break;
		}
		if (p->next==NULL) break;
	}
}


int faceContainsLaterals(void_list *p,void_list *q){
	return ( ( ((((ExtraVertex*)p->data)->lateral1)->index == ((Face*)q->data)->v1) ||
			((((ExtraVertex*)p->data)->lateral1)->index == ((Face*)q->data)->v2) ||
			((((ExtraVertex*)p->data)->lateral1)->index == ((Face*)q->data)->v3) ) && (
			((((ExtraVertex*)p->data)->lateral2)->index == ((Face*)q->data)->v1) ||
			((((ExtraVertex*)p->data)->lateral2)->index == ((Face*)q->data)->v2) ||
			((((ExtraVertex*)p->data)->lateral2)->index == ((Face*)q->data)->v3)));
}



void_list * addFaces(void_list *F,void_list *EV,int MF, int MV,int O,int I) {
	// add faces
	void_list *t,*th;
	Face *f;
	th = F;
	char buffer[128];
	sprintf (buffer, "Face 0 0 0 0");

	t = new void_list();
	t->next = th;
	f = new Face(buffer);
	t->data = (void*)f;
	th = t;
	((Face*)th->data)->index = MF;
	if(O) { ((Face*)th->data)->v1 = (((ExtraVertex*)EV->data)->lateral1)->index;}
	else { ((Face*)th->data)->v1 = (((ExtraVertex*)EV->data)->lateral2)->index;}
	((Face*)th->data)->v2 = MV;
	((Face*)th->data)->v3 = I;
	MF++;
	t = new void_list();
	t->next = th;
	f = new Face(buffer);
	t->data = (void*)f;
	th = t;
	((Face*)th->data)->index = MF;
	((Face*)th->data)->v1 = I;
	((Face*)th->data)->v2 = MV;
	if(O) { ((Face*)th->data)->v3 = (((ExtraVertex*)EV->data)->lateral2)->index; }
	else { ((Face*)th->data)->v3 = (((ExtraVertex*)EV->data)->lateral1)->index; }
	MF++;
	return th;
}



void_list * addVertex(void_list *V,void_list *EV,int MV) {
	// add extra vertex
	void_list *t,*th;
	Vertex *v;
	th = V;
	char buffer[128];
	sprintf (buffer, "Vertex 0 0 0 0");

	t = new void_list();
	t->next = th;
	v = new Vertex(buffer);
	t->data = (void*)v;
	th = t;
	((Vertex*)th->data)->index = MV;
	((Vertex*)th->data)->x = (((ExtraVertex*)EV->data)->v)->x;
	((Vertex*)th->data)->y = (((ExtraVertex*)EV->data)->v)->y;
	((Vertex*)th->data)->z = (((ExtraVertex*)EV->data)->v)->z;
	return th;
}

int main(int argc,char *argv[]){

	if (argc != 5)
	{
		printf("\nSyntax: meshstitch input_file1 input_file2 z_value orientation\n\n");
		printf("Detail: The z_value refers to the z axis coordinate of\n");
		printf("        the common contour between the two input meshes. \n");
		printf("        Orientation is '1' if the polygons to be removed\n");
		printf("        from input_file2 lie on the +z side of the z_value, and\n");
		printf("        and orientation is '-1' if the polygons to be removed\n");
		printf("        from input_file2 lie on the -z side of the z_value.\n");
		printf("        The polygons to be removed from input_file1 are assumed\n");
		printf("        to lie on the opposite of the z_value specified by orientation.\n");
		printf("Output: input_file1.clipped and input_file2.clipped.\n\n");
		return 1;
	}

	////////// declare variables /////////

	char *infile1;
	char *infile2;
	char outfile1[128];
	char outfile2[128];
	char line[2048];
	char *str;
	char *eptr;
	FILE *F;
	void_list *q,*qq,*p,*pp,*prev;
	bool flag;
	bool new_flag;
	double epsilon = 1e-6;


	// vertex linked lists for first and second file
	void_list *vlh1,*vl1;
	Vertex *v1;
	vlh1 = NULL;
	void_list *vlh2,*vl2;
	Vertex *v2;
	vlh2 = NULL;

	// face linked lists for first and second file
	void_list *flh1,*fl1;
	Face *f1;
	flh1 = NULL;
	void_list *flh2,*fl2;
	Face *f2;
	flh2 = NULL;

	// vertex with target z_value linked lists for first and second file
	void_list *v1_match_h,*v2_match_h;
	v1_match_h = NULL;
	v2_match_h = NULL;

	// shared vertex linked lists for first and second file
	void_list *v1_shared,*v1_shared_h;
	v1_shared_h = NULL;
	void_list *v2_shared,*v2_shared_h;
	v2_shared_h = NULL;

	// candidate bad vertex linked lists for first and second file
	void_list *cbv1h,*cbv2h;
	cbv1h = NULL;
	cbv2h = NULL;

	// candidate bad face linked lists for first and second file
	void_list *cbf1h,*cbf2h;
	cbf1h = NULL;
	cbf2h = NULL;

	// bad vertex linked lists for first and second file
	void_list *bv1,*bv1h;
	bv1h = NULL;
	void_list *bv2,*bv2h;
	bv2h = NULL;

	// bad face linked lists for first and second file
	void_list *bf1,*bf1h;
	bf1h = NULL;
	void_list *bf2,*bf2h;
	bf2h = NULL;

	// pointers to start of vertex and face linked lists for first and second file
	void_list *v1end;
	void_list *v2end;
	void_list *f1end;
	void_list *f2end;
	v1end = NULL;
	v2end = NULL;
	f1end = NULL;
	f2end = NULL;

	// variables used to find extra vertices
	void_list *ev1h;
	void_list *ev2h;
	ev1h = NULL;
	ev2h = NULL;
	int max_vertex1;
	int max_vertex2;
	int max_face1;
	int max_face2;
	int third_index;
	int orient;

	// get z_value
	int z_value;
	z_value = (int) (strtod(argv[3],&eptr));

	// get target z_value
	int orientation;
	orientation = (int) strtod(argv[4],&eptr);
	int val;

	////////// get data /////////

	// open first file
	infile1 = argv[1];
	F = fopen(infile1,"r");
	if (!F) {
	printf("Couldn't open input file %s\n",infile1);
	return 1;
	}

	// for every line in first file
	for (str=fgets(line,2048,F) ; str!=NULL ; str=fgets(line,2048,F)) {

		// skip leading whitespace
		while (strchr(" \t,",*str)!=NULL) { str++;}

		// if first character is V for Vertex, add new linked list class instance
		if (strchr("V",*str)!=NULL){
			vl1 = new void_list();
			vl1->next = vlh1;
			v1 = new Vertex(str);
			vl1->data = (void*)v1;
			vlh1 = vl1;
		} 
		// if first character is F for Face, add new linked list class instance
		else if (strchr("F",*str)!=NULL){
			fl1 = new void_list();
			fl1->next = flh1;
			f1 = new Face(str);
			fl1->data = (void*)f1;
			flh1 = fl1;
		}
	}
	fclose(F);

	// opend second file
	infile2 = argv[2];
	F = fopen(infile2,"r");
	if (!F) {
	printf("Couldn't open input file %s\n",infile2);
	return 1;
	}

	// for every line in second file
	for (str=fgets(line,2048,F) ; str!=NULL ; str=fgets(line,2048,F)) {

		// skip leading whitespace
		while (strchr(" \t,",*str)!=NULL) { str++;}

		// if first character is V for Vertex, add new linked list class instance
		if (strchr("V",*str)!=NULL){
			vl2 = new void_list();
			vl2->next = vlh2;
			v2 = new Vertex(str);
			vl2->data = (void*)v2;
			vlh2 = vl2;
		} 
		// if first character is F for Face, add new linked list class instance
		else if (strchr("F",*str)!=NULL){
			fl2 = new void_list();
			fl2->next = flh2;
			f2 = new Face(str);
			fl2->data = (void*)f2;
			flh2 = fl2;
		}
	}
	fclose(F);

	// maked linked lists bidirectional
	vlh1 = addPrevious(vlh1);
	flh1 = addPrevious(flh1);
	vlh2 = addPrevious(vlh2);
	flh2 = addPrevious(flh2);

	////////// find shared vertices /////////

	// build linked list of first and second vertex elements with z_value matching target
	v1_match_h = matchZ(vlh1,z_value);
	v2_match_h = matchZ(vlh2,z_value);

	// find shared vertices
	for (q=v1_match_h;q!=NULL;q=q->next) {
		for (qq=v2_match_h;qq!=NULL;qq=qq->next) {
			if ( !distinguishable(((Vertex*)q->data)->x,((Vertex*)qq->data)->x,epsilon) &&
				!distinguishable(((Vertex*)q->data)->y,((Vertex*)qq->data)->y,epsilon) &&
				!distinguishable(((Vertex*)q->data)->z,((Vertex*)qq->data)->z,epsilon)){
				v1_shared = new void_list();
				v1_shared->next = v1_shared_h;
				v1_shared->data = (Vertex*)q->data;
				v1_shared_h = v1_shared;
				v2_shared = new void_list();
				v2_shared->next = v2_shared_h;
				v2_shared->data = (Vertex*)qq->data;
				v2_shared_h = v2_shared;
			}
			if (qq->next==NULL) break;
		}
		if (q->next==NULL) break;
	}

	///// STRATEGY /////
	// (1) collect all vertices with z value half section
	//		thickness away from target z_value in orientation direction
	// (2) collect all polygons connected to vertices from above
	// (3) for each collected polygon, if at least one vertex is a shared vertex from above
	// 		then any vertex on that polygon with z value half section
	//		thickness away from target z_value in orientation direction is bad
	// (4) for each collected polygon, if at least one vertex is a bad vertex from above
	// 		then polygon is bad

	// gather candidate bad vertices from first and second vertex linked list
	val = z_value+orientation*25*(-1);
	cbv1h = candidateBadVertices(vlh1,val);
	val = z_value+orientation*25;
	cbv2h = candidateBadVertices(vlh2,val);

	///// find faces that have at least one vertex from cbv lists /////
	cbf1h = candidateBadFaces(flh1,cbv1h);
	cbf2h = candidateBadFaces(flh2,cbv2h);

	// maked linked lists bidirectional
	cbv1h = addPrevious(cbv1h);
	cbv2h = addPrevious(cbv2h);
	cbf1h = addPrevious(cbf1h);
	cbf2h = addPrevious(cbf2h);

	///// find bad vertices and faces from first vertex and first face linked lists
	void_list *ptr1,*ptr2;
	ptr1=cbf1h;
	ptr2=cbv1h;
	// for each candidate bad face
	for (p=cbf1h;p!=NULL;p=p->next) {

		// search first shared vertex list, v1_shared_h
		flag = true;
		q=v1_shared_h;
		while(q!=NULL && flag) {
			// if shared vertex is at least one of three face vertices
    	  	if (vertexInFace(p,q)) {
				// add face to bad face list
				bf1 = new void_list();
				bf1->next = bf1h;
				bf1->data = (Face*)p->data;
				bf1h = bf1;
				flag = false;
				// and remove face from candidate face list
				removeLink(p);
				// adjust list pointer
				if (p==ptr1) { ptr1 = p->next; }
				// search candidate bad vertex list for each face vertex
				pp=ptr2;
				while(pp!=NULL) {
					// if bad vertex is at least one of three face vertices
		    	  	if (vertexInFace(p,pp)) {
						// add vertex to bad vertex list
						bv1 = new void_list();
						bv1->next = bv1h;
						bv1->data = (Vertex*)pp->data;
						bv1h = bv1;
						// and remove vertex from candidate vertex list
						removeLink(pp);
						// adjust list pointer
						if (pp==ptr2) { ptr2 = pp->next; }
					}
					if (pp->next==NULL) break;
					pp=pp->next;
				}
			}
			if (q->next==NULL) break;
			q=q->next;
		}
		if (p->next==NULL) break;
	}
	// adjust pointers
	cbf1h=ptr1;
	cbv1h=ptr2;

	// rescan candidate bad faces searching for known bad vertices
	bf1h = findBadFaces(cbf1h,bv1h,bf1h);

	///// find bad vertices and faces from second vertex and second face linked lists
	ptr1=cbf2h;
	ptr2=cbv2h;
	// for each candidate bad face
	for (p=cbf2h;p!=NULL;p=p->next) {

		// search second shared vertex list, v2_shared_h
		flag = true;
		q=v2_shared_h;
		while(q!=NULL && flag) {
			// if shared vertex is at least one of three face vertices
    	  	if (vertexInFace(p,q)) {
				// add face to bad face list
				bf2 = new void_list();
				bf2->next = bf2h;
				bf2->data = (Face*)p->data;
				bf2h = bf2;
				flag = false;
				// and remove face from candidate face list
				removeLink(p);
				// adjust list pointer
				if (p==ptr1) { ptr1 = p->next; }
				// search candidate bad vertex list for each face vertex
				pp=ptr2;
				while(pp!=NULL) {
					// if bad vertex is at least one of three face vertices
		    	  	if (vertexInFace(p,pp)) {
						// add vertex to bad vertex list
						bv2 = new void_list();
						bv2->next = bv2h;
						bv2->data = (Vertex*)pp->data;
						bv2h = bv2;
						// and remove vertex from candidate vertex list
						removeLink(pp);
						// adjust list pointer
						if (pp==ptr2) { ptr2 = pp->next; }
					}
					if (pp->next==NULL) break;
					pp=pp->next;
				}
			}
			if (q->next==NULL) break;
			q=q->next;
		}
		if (p->next==NULL) break;
	}
	// adjust pointers
	cbf2h=ptr1;
	cbv2h=ptr2;

	// rescan candidate bad faces searching for known bad vertices
	bf2h = findBadFaces(cbf2h,bv2h,bf2h);

	///// remove bad vertices from first and second vertex list /////
	vlh1 = removeBadVertices(vlh1,bv1h);
	vlh2 = removeBadVertices(vlh2,bv2h);

if(0){
    for (p=vlh1;p!=NULL;p=p->next) {
        printf("vertices1 %i %.15g %.15g %.15g\n",
                ((Vertex*)p->data)->index,
                ((Vertex*)p->data)->x,
                ((Vertex*)p->data)->y,
                ((Vertex*)p->data)->z);
        if (p->next==NULL) break;
    }
    for (p=vlh2;p!=NULL;p=p->next) {
        printf("vertices2 %i %.15g %.15g %.15g\n",
                ((Vertex*)p->data)->index,
                ((Vertex*)p->data)->x,
                ((Vertex*)p->data)->y,
                ((Vertex*)p->data)->z);
        if (p->next==NULL) break;
    }
}
	///// remove bad faces from first and second face list /////
	flh1 = removeBadFaces(flh1,bf1h);
	flh2 = removeBadFaces(flh2,bf2h);

if(0){
    for (p=flh1;p!=NULL;p=p->next) {
        printf("faces1 %i %i %i %i\n",
                ((Face*)p->data)->index,
                ((Face*)p->data)->v1,
                ((Face*)p->data)->v2,
                ((Face*)p->data)->v3);
        if (p->next==NULL) break;
    }
    for (p=flh2;p!=NULL;p=p->next) {
        printf("faces2 %i %i %i %i\n",
                ((Face*)p->data)->index,
                ((Face*)p->data)->v1,
                ((Face*)p->data)->v2,
                ((Face*)p->data)->v3);
        if (p->next==NULL) break;
    }
}
	////////// find extra vertices //////////

if(0){
    for (p=v2_match_h;p!=NULL;p=p->next) {
        printf("match vertices2 %i %.15g %.15g %.15g\n",
                ((Vertex*)p->data)->index,
                ((Vertex*)p->data)->x,
                ((Vertex*)p->data)->y,
                ((Vertex*)p->data)->z);
        if (p->next==NULL) break;
    }
    for (p=v2_shared;p!=NULL;p=p->next) {
        printf("shared vertices2 %i %.15g %.15g %.15g\n",
                ((Vertex*)p->data)->index,
                ((Vertex*)p->data)->x,
                ((Vertex*)p->data)->y,
                ((Vertex*)p->data)->z);
        if (p->next==NULL) break;
    }
    for (p=bf2h;p!=NULL;p=p->next) {
        printf("bad faces2 %i %i %i %i\n",
                ((Face*)p->data)->index,
                ((Face*)p->data)->v1,
                ((Face*)p->data)->v2,
                ((Face*)p->data)->v3);
        if (p->next==NULL) break;
    }
}

	ev1h = extraVertices(v1_match_h,v1_shared_h,bf1h);
	ev2h = extraVertices(v2_match_h,v2_shared_h,bf2h);
	lateralVertices(ev1h,bf1h,v1_shared);
	lateralVertices(ev2h,bf2h,v2_shared);
	convertLateral(vlh2,ev1h,epsilon);
	convertLateral(vlh1,ev2h,epsilon);

if(0){
    for (p=ev1h;p!=NULL;p=p->next) {
        printf("extra vertices1 %i %.15g %.15g %.15g\n",
                (((ExtraVertex*)p->data)->v)->index,
                (((ExtraVertex*)p->data)->v)->x,
                (((ExtraVertex*)p->data)->v)->y,
                (((ExtraVertex*)p->data)->v)->z);
        if (p->next==NULL) break;
    }
    for (p=ev2h;p!=NULL;p=p->next) {
        printf("extra vertices2 %i %.15g %.15g %.15g\n",
                (((ExtraVertex*)p->data)->v)->index,
                (((ExtraVertex*)p->data)->v)->x,
                (((ExtraVertex*)p->data)->v)->y,
                (((ExtraVertex*)p->data)->v)->z);
        printf("lateral vertices2 %i %.15g %.15g %.15g\n",
				(((ExtraVertex*)p->data)->lateral1)->index,
				(((ExtraVertex*)p->data)->lateral1)->x,
				(((ExtraVertex*)p->data)->lateral1)->y,
				(((ExtraVertex*)p->data)->lateral1)->z);
        printf("lateral vertices2 %i %.15g %.15g %.15g\n",
				(((ExtraVertex*)p->data)->lateral2)->index,
				(((ExtraVertex*)p->data)->lateral2)->x,
				(((ExtraVertex*)p->data)->lateral2)->y,
				(((ExtraVertex*)p->data)->lateral2)->z);
        if (p->next==NULL) break;
    }
}
	// find max vertex indices
	max_vertex1 = 0;
	for (p=vlh1;p!=NULL;p=p->next) {
		if (((Vertex*)p->data)->index > max_vertex1) {max_vertex1 = ((Vertex*)p->data)->index;}
		if (p->next==NULL) break;
	}
	max_vertex1++;
	max_vertex2 = 0;
	for (p=vlh2;p!=NULL;p=p->next) {
		if (((Vertex*)p->data)->index > max_vertex1) {max_vertex2 = ((Vertex*)p->data)->index;}
		if (p->next==NULL) break;
	}
	max_vertex2++;

	// find max face index
	max_face1 = 0;
	for (p=flh1;p!=NULL;p=p->next) {
		if (((Face*)p->data)->index > max_face1) {max_face1 = ((Face*)p->data)->index;}
		if (p->next==NULL) break;
	}
	max_face1++;
	max_face2 = 0;
	for (p=flh2;p!=NULL;p=p->next) {
		if (((Face*)p->data)->index > max_face2) {max_face2 = ((Face*)p->data)->index;}
		if (p->next==NULL) break;
	}
	max_face2++;

	// for each extra vertex
	ptr2 = flh2;
	for (p=ev1h;p!=NULL;p=p->next) {

		// for each face in object1
		for (q=flh2;q!=NULL;q=q->next) {

			// if face contains both lateral vertices
			if(faceContainsLaterals(p,q)){

				// identify third vertex
				if( ((((ExtraVertex*)p->data)->lateral1)->index != ((Face*)q->data)->v1) &&
					((((ExtraVertex*)p->data)->lateral2)->index != ((Face*)q->data)->v1)
					) {third_index = ((Face*)q->data)->v1;}
				else if( 
					((((ExtraVertex*)p->data)->lateral1)->index != ((Face*)q->data)->v2) &&
					((((ExtraVertex*)p->data)->lateral2)->index != ((Face*)q->data)->v2)
					) {third_index = ((Face*)q->data)->v2;}
				else {third_index = ((Face*)q->data)->v3;}

				// identify orientation
				if( 
					( ((((ExtraVertex*)p->data)->lateral1)->index == ((Face*)q->data)->v1) &&
					((((ExtraVertex*)p->data)->lateral2)->index == ((Face*)q->data)->v2)) ||
					( ((((ExtraVertex*)p->data)->lateral1)->index == ((Face*)q->data)->v3) &&
					((((ExtraVertex*)p->data)->lateral2)->index == ((Face*)q->data)->v1)) ||
					( ((((ExtraVertex*)p->data)->lateral1)->index == ((Face*)q->data)->v2) &&
					((((ExtraVertex*)p->data)->lateral2)->index == ((Face*)q->data)->v3))
					){ orient = 1;}
				else {orient = 0;}

				// remove face
				removeLink(q);
				// adjust list pointer
				if (q==ptr2) { ptr2 = q->next; }
			}
			if (q->next==NULL) break;
		}
		// adjust list pointer
		flh2=ptr2;

		// add faces
		flh2 = addFaces(flh2,p,max_face2,max_vertex2,orient,third_index);
		max_face2++;
		max_face2++;
		ptr2 = flh2;

		// add extra vertex
		vlh2 = addVertex(vlh2,p,max_vertex2);
		max_vertex2++;

		if (p->next==NULL) break;
	}

	// for each extra vertex
	ptr1 = flh1;
	for (p=ev2h;p!=NULL;p=p->next) {

		// for each face in object1
		for (q=flh1;q!=NULL;q=q->next) {

			// if face contains both lateral vertices
			if(faceContainsLaterals(p,q)){

				// identify third vertex
				if( ((((ExtraVertex*)p->data)->lateral1)->index != ((Face*)q->data)->v1) &&
					((((ExtraVertex*)p->data)->lateral2)->index != ((Face*)q->data)->v1)
					) {third_index = ((Face*)q->data)->v1;}
				else if( 
					((((ExtraVertex*)p->data)->lateral1)->index != ((Face*)q->data)->v2) &&
					((((ExtraVertex*)p->data)->lateral2)->index != ((Face*)q->data)->v2)
					) {third_index = ((Face*)q->data)->v2;}
				else {third_index = ((Face*)q->data)->v3;}

				// identify orientation
				if( 
					( ((((ExtraVertex*)p->data)->lateral1)->index == ((Face*)q->data)->v1) &&
					((((ExtraVertex*)p->data)->lateral2)->index == ((Face*)q->data)->v2)) ||
					( ((((ExtraVertex*)p->data)->lateral1)->index == ((Face*)q->data)->v3) &&
					((((ExtraVertex*)p->data)->lateral2)->index == ((Face*)q->data)->v1)) ||
					( ((((ExtraVertex*)p->data)->lateral1)->index == ((Face*)q->data)->v2) &&
					((((ExtraVertex*)p->data)->lateral2)->index == ((Face*)q->data)->v3))
					){ orient = 1;}
				else {orient = 0;}

				// remove face
				removeLink(q);
				// adjust list pointer
				if (q==ptr1) { ptr1 = q->next; }
			}
			if (q->next==NULL) break;
		}
		// adjust list pointer
		flh1=ptr1;

		// add faces
		flh1 = addFaces(flh1,p,max_face1,max_vertex1,orient,third_index);
		max_face1++;
		max_face1++;
		ptr1 = flh1;

		// add extra vertex
		vlh1 = addVertex(vlh1,p,max_vertex1);
		max_vertex1++;

		if(1){
		        printf("\n\nextra vertices2 %i %.15g %.15g %.15g\n",
		                (((ExtraVertex*)p->data)->v)->index,
		                (((ExtraVertex*)p->data)->v)->x,
		                (((ExtraVertex*)p->data)->v)->y,
		                (((ExtraVertex*)p->data)->v)->z);
				fflush(stdout);
		        printf("lateral vertices2 %i %.15g %.15g %.15g\n",
						(((ExtraVertex*)p->data)->lateral1)->index,
						(((ExtraVertex*)p->data)->lateral1)->x,
						(((ExtraVertex*)p->data)->lateral1)->y,
						(((ExtraVertex*)p->data)->lateral1)->z);
				fflush(stdout);
		        printf("lateral vertices2 %i %.15g %.15g %.15g\n",
						(((ExtraVertex*)p->data)->lateral2)->index,
						(((ExtraVertex*)p->data)->lateral2)->x,
						(((ExtraVertex*)p->data)->lateral2)->y,
						(((ExtraVertex*)p->data)->lateral2)->z);
				fflush(stdout);
	        	printf("added faces2 %i %i %i %i\n",
		                ((Face*)flh1->data)->index,
		                ((Face*)flh1->data)->v1,
		                ((Face*)flh1->data)->v2,
		                ((Face*)flh1->data)->v3);
				fflush(stdout);
	        	printf("added faces2 %i %i %i %i\n",
		                ((Face*)(flh1->next)->data)->index,
		                ((Face*)(flh1->next)->data)->v1,
		                ((Face*)(flh1->next)->data)->v2,
		                ((Face*)(flh1->next)->data)->v3);
				fflush(stdout);
        		printf("added vertex2 %i %.15g %.15g %.15g\n",
		                ((Vertex*)vlh1->data)->index,
		                ((Vertex*)vlh1->data)->x,
		                ((Vertex*)vlh1->data)->y,
		                ((Vertex*)vlh1->data)->z);
        		printf("third_index2 = %i\n",third_index);
				fflush(stdout);
		}


		if (p->next==NULL) break;
	}

	// go through linked list backwards and add previous pointers
	vlh1 = addPrevious(vlh1);
	flh1 = addPrevious(flh1);
	vlh2 = addPrevious(vlh2);
	flh2 = addPrevious(flh2);

	////////// write output to files //////////
	// set end pointers
    for (q=vlh1;q!=NULL;q=q->next) { v1end=q; if (q->next==NULL) break; }
    for (q=flh1;q!=NULL;q=q->next) { f1end=q; if (q->next==NULL) break; }
    for (q=vlh2;q!=NULL;q=q->next) { v2end=q; if (q->next==NULL) break; }
    for (q=flh2;q!=NULL;q=q->next) { f2end=q; if (q->next==NULL) break; }

	// open first file
	sprintf(outfile1,"%s.clipped",infile1);
	printf("\n%s written.\n",outfile1);
	F = fopen(outfile1,"w");
	if (!F) {
		printf("Couldn't open output file %s\n",outfile1);
		return 1;
	}

	// for every vertex in first file
	for (p=v1end;p!=NULL;p=p->previous) {
		sprintf(line,"Vertex %i  %.15g %.15g %.15g\n",
				((Vertex*)p->data)->index,
				((Vertex*)p->data)->x,
				((Vertex*)p->data)->y,
				((Vertex*)p->data)->z);
		fputs(line,F);
		if (p->previous==NULL) break;
	}

	// for every face in first file
	for (p=f1end;p!=NULL;p=p->previous) {
		sprintf(line,"Face %i  %i %i %i\n",
				((Face*)p->data)->index,
				((Face*)p->data)->v1,
				((Face*)p->data)->v2,
				((Face*)p->data)->v3);
		fputs(line,F);
		if (p->previous==NULL) break;
	}
	fclose(F);
	
	// open second file
	sprintf(outfile2,"%s.clipped",infile2);
	printf("%s written.\n\n",outfile2);
	F = fopen(outfile2,"w");
	if (!F) {
		printf("Couldn't open output file %s\n",outfile2);
		return 1;
	}

	// for every vertex in second file
	for (p=v2end;p!=NULL;p=p->previous) {
		sprintf(line,"Vertex %i  %.15g %.15g %.15g\n",
				((Vertex*)p->data)->index,
				((Vertex*)p->data)->x,
				((Vertex*)p->data)->y,
				((Vertex*)p->data)->z);
		fputs(line,F);
		if (p->previous==NULL) break;
	}

	// for every face in second file
	for (p=f2end;p!=NULL;p=p->previous) {
		sprintf(line,"Face %i  %i %i %i\n",
				((Face*)p->data)->index,
				((Face*)p->data)->v1,
				((Face*)p->data)->v2,
				((Face*)p->data)->v3);
		fputs(line,F);
		if (p->previous==NULL) break;
	}
	fclose(F);

	return 0;
}
