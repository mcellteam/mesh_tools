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

class Distance
{
public:
  double d; // distance between vertices with indices vA and vB
  Vertex *vA,*vB;
  Distance(double,void_list*,void_list*);
};

Distance::Distance(double dist,void_list *p,void_list *q)
{
	d = dist;
	vA = (Vertex*)p->data;
	vB = (Vertex*)q->data;
}

class Edge
{
public:
  Vertex *va,*vb; // pointers to vertex links
  int c12,c21;  // count of times edge was traversed from va to vb and vb to va
  void_list *f;
  Edge(Vertex*,Vertex*);
  void addFace(Face*);
  ~Edge(void);
};

Edge::Edge(Vertex *v1,Vertex *v2)
{
	va = v1;
	vb = v2;
	c12 = 1;
	c21 = 0;
	f=NULL;
}

Edge::~Edge(void)
{
	void_list *q,*p;
	q=f;
	while (q!=NULL) {
		p=q->next;
		delete q;
		q=p;
	}
}

void Edge::addFace(Face *fc) {
	//add face to edge
	void_list *qq;
	qq = new void_list();
	qq->next = f;
	qq->data = (void*)fc;
	f = qq;
}

class Data
{
public:
  Edge *ptr; // pointer to edge link
  int vert; // vertex of pair 
  Data(Edge*,int);
};

Data::Data(Edge *e,int v)
{
	ptr = e;
	vert = v;
}

int compare (const void* a, const void* b ) {

	double i;
	void_list *j,*k;

	j = *(void_list**)a;
	k = *(void_list**)b;

	i = ((Distance*)j->data)->d - ((Distance*)k->data)->d;
	if (i<0) {	return -1;}
	else if (i>0) {	return 1;}
	else { return (0);}
}

int compare_v (const void* a, const void* b ) {

	double i;
	void_list *j,*k;

	j = *(void_list**)a;
	k = *(void_list**)b;

	i = ((Vertex*)j->data)->index - ((Vertex*)k->data)->index;
	if (i<0) {	return -1;}
	else if (i>0) {	return 1;}
	else { return (0);}
}

void_list* removeLink(void_list* L) {
	void_list *q;
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
	// update pointer
	q=L->next;
	delete L;
	return q;
}

void_list * findEdge(void_list *F,void_list *E,int v1,int v2,void_list** vert_array,void_list** pairs) {
    ///// find corresponding edge /////
    void_list *edge,*p,*qq;
    Edge* e;
	Data* d;
    bool flag=false;

	// look up first vertex
    p=pairs[v1];
    while (p!=NULL && !flag) {
		if((int)((Data*)p->data)->vert==v2){flag=true;}
		else{p=p->next;}
	}
    if (flag) {
		d=(Data*)p->data;
		e=(Edge*)d->ptr;
		// update traversal
		e->c12 += 1;
		//add face to edge
  		e->addFace((Face*)F->data);
    } 
	// if not found under first vertex,
	// look up second vertex
	if(!flag){
	    p=pairs[v2];
	    while (p!=NULL && !flag) {
			if((int)((Data*)p->data)->vert==v1){flag=true;}
			else{p=p->next;}
		}
	    if (flag) {
			d=(Data*)p->data;
			e=(Edge*)d->ptr;
			// update traversal
			e->c21 += 1;
			//add face to edge
  			e->addFace((Face*)F->data);
	    } 
	}

    // if corresponding edge was not found
	if (!flag){
        // add edge
        edge = new void_list();
        edge->next = E;
        e = new Edge((Vertex*)vert_array[v1]->data,(Vertex*)vert_array[v2]->data);
        // add face to edge
  		e->addFace((Face*)F->data);
        edge->data = (void*)e;
        E = edge;
		// update pairs array
		qq = new void_list();
		qq->next = pairs[v1];
		d = new Data(e,v2);
		qq->data = (void*)d;
		pairs[v1] = qq;
    }
    return E;
}

void deletePairs(void_list **pairs,int max_vert) {
	///// delete pairs /////
	int i;
	void_list *p,*q;
	for(i=0;i<max_vert+1;i++){
		// for each vertex pair
		p=pairs[i];
		while (p!=NULL) {
			q=p->next;
			delete (Data*)p->data;
			delete p;
			p=q;
		}
	}
	delete[] pairs;
}

void_list * gatherFreeEdges(void_list *F,void_list *E,void_list** vert_array,int print_flag,int max_vert,void_list ** &pairs) {
	if (!print_flag) {printf("Building edge list...");fflush(stdout);}
    ///// build edge list ////
	void_list *q,*ptr,*prev,*p;
	// initialize vertex pairs data structure
	int i;
//	void_list** pairs;
	pairs = new void_list*[max_vert+1];
	for(i=0;i<max_vert+1;i++){
		pairs[i]=NULL;
	}

    // for each face
	i = 0;
    for (q=F;q!=NULL;q=q->next) {
        E = findEdge(q,E,((Face*)q->data)->v1,((Face*)q->data)->v2,vert_array,pairs);
        E = findEdge(q,E,((Face*)q->data)->v2,((Face*)q->data)->v3,vert_array,pairs);
        E = findEdge(q,E,((Face*)q->data)->v3,((Face*)q->data)->v1,vert_array,pairs);
    }

	// go through linked list backwards and add previous pointers
	// for each free face
	prev = NULL;
	for (p=E;p!=NULL;p=p->next) {
		p->previous = prev;
		prev = p;
		if (p->next==NULL) break;
	}

	///// gather free edges /////
	ptr=E;
	q=E;
    while (q!=NULL) {
        // if edge is not free
		if ((((Edge*)q->data)->c12+((Edge*)q->data)->c21) != 1) {
			// delete edge
			delete (Edge*)q->data;
			// adjust list pointer
			if (q==ptr) { ptr=q->next;}
			// remove link from list
			q=removeLink(q);
		} else {q=q->next;}
    }
	// adjust pointer
	E=ptr;

	if(0){
	    for (q=E;q!=NULL;q=q->next) {
			printf("Free edge vertices %i %i\n",
				((Vertex*)((Edge*)q->data)->va)->index,
				((Vertex*)((Edge*)q->data)->vb)->index);
	    }
	}

	///// delete pairs /////
//	deletePairs(pairs,max_vert);

	if (!print_flag) {printf("complete.\n");fflush(stdout);}
	return E;
}

void_list * deleteFreeEdges(void_list *E) {
    ///// delete edge list ////
	void_list *p,*q,*pp,*qq;
	// for each free edge
	p=E;
	while (p!=NULL) {
		q=p->next;
		delete (Edge*)p->data;
		delete p;
		p=q;
	}
	delete p;
	E=NULL;
	return E;
}

void_list * gatherFreeFaces(void_list *E,int print_flag) {
	if (!print_flag) {printf("Gathering free faces...");fflush(stdout);}
	///// gather free faces /////
	void_list *q,*p,*qq,*ptr,*ffree,*ffreeh,*prev;
	ffreeh = NULL;
	int max_face=0,i;
	// for each free edge
    for (q=E;q!=NULL;q=q->next) {
		//  for each associated face
    	for (qq=((Edge*)q->data)->f;qq!=NULL;qq=qq->next) {
			ffree = new void_list();
			ffree->next = ffreeh;
			ffree->data = (void*)(qq->data);
			ffreeh = ffree;
			if(max_face<((Face*)ffreeh->data)->index)
				{max_face=((Face*)ffreeh->data)->index;}
	    }
    }

	// go through linked list backwards and add previous pointers
	// for each free face
	prev = NULL;
	for (p=ffreeh;p!=NULL;p=p->next) {
		p->previous = prev;
		prev = p;
		if (p->next==NULL) break;
	}

	///// gather unique free faces /////
	int f_array[max_face+1];
    for (i=0;i<max_face+1;i++) { f_array[i]=0; }
	ptr=ffreeh;
	q=ffreeh;
    while (q!=NULL) {
		if(!f_array[((Face*)q->data)->index])
			{f_array[((Face*)q->data)->index]=1;q=q->next;}
		else {
			// adjust list pointer
			if (q==ptr) { ptr=q->next;}
			// remove free face link
			q=removeLink(q);
		}
    }
	// adjust pointer
	ffreeh=ptr;

	if(0){
	    for (q=ffreeh;q!=NULL;q=q->next) {
			printf("Free face %i %i %i %i\n",
				((Face*)q->data)->index,
				((Face*)q->data)->v1,
				((Face*)q->data)->v2,
				((Face*)q->data)->v3);
	    }
	}

	if (!print_flag) {printf("complete.\n");fflush(stdout);}
	return ffreeh;
}


void_list * gatherFreeVertices(void_list *E,int max_vert,int print_flag) {
	if (!print_flag) {printf("Gathering free vertices...");fflush(stdout);}
	///// gather free vertices /////
	void_list *q,*p,*qq,*ptr,*vfree,*vfreeh,*prev;
	vfreeh = NULL;
	int i;
	// for each edge
    for (q=E;q!=NULL;q=q->next) {
		vfree = new void_list();
		vfree->next = vfreeh;
		vfree->data = (void*)((Edge*)q->data)->va;
		vfreeh = vfree;
		vfree = new void_list();
		vfree->next = vfreeh;
		vfree->data = (void*)((Edge*)q->data)->vb;
		vfreeh = vfree;
	}

	// go through linked list backwards and add previous pointers
	// for each free vertex
	prev = NULL;
	for (p=vfreeh;p!=NULL;p=p->next) {
		p->previous = prev;
		prev = p;
		if (p->next==NULL) break;
	}

	///// gather unique free vertices /////
	int v_array[max_vert+1];
    for (i=0;i<max_vert+1;i++) { v_array[i]=0; }
	ptr=vfreeh;
    q=vfreeh;
    while(q!=NULL) {
		if(!v_array[((Vertex*)q->data)->index])
			{v_array[((Vertex*)q->data)->index]=1;q=q->next;}
		else {
			// adjust list pointer
			if (q==ptr) { ptr=q->next;}
			// remove free face link
			q = removeLink(q);
		}
    }
	// adjust pointer
	vfreeh=ptr;

	if(0){
	    for (q=vfreeh;q!=NULL;q=q->next) {
			printf("Free vertex %i %.15g %.15g %.15g\n",
				((Vertex*)q->data)->index,
				((Vertex*)q->data)->x,
				((Vertex*)q->data)->y,
				((Vertex*)q->data)->z);
	    }
	}

	if (!print_flag) {printf("complete.\n");fflush(stdout);}
	return vfreeh;
}

void deleteFreeVertices(void_list *v) {
    ///// delete edge list ////
	void_list *p,*q;
	// for each free edge
	p=v;
	while (p!=NULL) {
		q=p->next;
		delete p;
		p=q;
	}
}

int countFreeVertices(void_list *v,int print_flag) {
	///// count number free vertices /////
	void_list *q;
	int i=0;
    for (q=v;q!=NULL;q=q->next) { i++; }
	if(!print_flag){printf("num free vertices = %i\n",i);fflush(stdout);}
	return i;
}

void_list ** computeDistances(void_list **dist_array,void_list *v,void_list *F,int& d_count,double epsilon,int print_flag,void_list ** pairs) {
	///// compute distances /////
	if (!print_flag) {printf("Computing vertex distances...");fflush(stdout);}
	void_list *q,*qq,*pp,*dl,*dlh,*p;
	double diffx,diffy,diffz,squared_dist;
	bool face_found;
	int i;
	Distance *d;
	dlh = NULL;
	d_count = 0;
	// for each free vertex
    for (q=v;q!=NULL;q=q->next) {
		// for each free vertex
	    for (qq=v;qq!=NULL;qq=qq->next) {
			// if combination is unique
			if (((Vertex*)q->data)->index>((Vertex*)qq->data)->index){

				// are vertices members of the same face?
				// look up first vertex
				face_found=false;
			    p=pairs[((Vertex*)q->data)->index];
			    while (p!=NULL && !face_found) {
					if((int)((Data*)p->data)->vert==((Vertex*)qq->data)->index){face_found=true;}
					else{p=p->next;}
				}
			    if (!face_found) {
					// look up second vertex
				    p=pairs[((Vertex*)qq->data)->index];
				    while (p!=NULL && !face_found) {
						if((int)((Data*)p->data)->vert==((Vertex*)q->data)->index){face_found=true;}
						else{p=p->next;}
					}
				}
				// if not
				if (!face_found) {
					diffx = ((Vertex*)q->data)->x/epsilon-((Vertex*)qq->data)->x/epsilon;
					diffy = ((Vertex*)q->data)->y/epsilon-((Vertex*)qq->data)->y/epsilon;
					diffz = ((Vertex*)q->data)->z/epsilon-((Vertex*)qq->data)->z/epsilon;
					squared_dist = diffx*diffx+diffy*diffy+diffz*diffz;
					// create distance
					dl = new void_list();
					dl->next = dlh;
					d = new Distance(squared_dist,q,qq);
					dl->data = (void*)d;
					dlh = dl;
					d_count++;
				}
			}
		}
	}
	
	if(0){
	    for (q=dlh;q!=NULL;q=q->next) {
			printf("Free vertices %i and %i, dist = %.15g\n",
				((Vertex*)((Distance*)q->data)->vA)->index,
				((Vertex*)((Distance*)q->data)->vB)->index,
				((Distance*)q->data)->d);
	    }
	}

	///// sort distances /////
	// create array of distances
	dist_array = new void_list*[d_count];
	i=0;
    for (q=dlh;q!=NULL;q=q->next) {
		dist_array[i++]=q;
    }
	qsort(dist_array,d_count,sizeof(int),compare);

	if(0){ for (i=0;i<d_count;i++){ printf("Free vertices %i and %i, dist = %.15g\n",
				((Vertex*)((Distance*)(dist_array[i])->data)->vA)->index,
				((Vertex*)((Distance*)(dist_array[i])->data)->vB)->index,
				((Distance*)(dist_array[i])->data)->d); } }

	if (!print_flag) {printf("complete.\n");fflush(stdout);}
	return dist_array;
}

int findVerticesToMerge(void_list **array,int count,void_list *v,int print_flag) {
	if (!print_flag) {printf("Identify target vertex pair...");fflush(stdout);}
	int i=0;
	void_list *q;
	bool flag=false,vflag1,vflag2;
	while (i<count && !flag) {
		if(0){
			printf("candidate vertices to merge: %i %i, dist = %.15g\n",
				((Vertex*)((Distance*)(array[i])->data)->vA)->index,
				((Vertex*)((Distance*)(array[i])->data)->vB)->index,
				((Distance*)(array[i])->data)->d
				);
		}
		// is ith vertex pair free?
		vflag1=false;vflag2=false;
		q=v;
		// for each free vertex
		while(q!=NULL && (!vflag1 || !vflag2)){
			if ( (((Vertex*)((Distance*)(array[i])->data)->vA)->index == ((Vertex*)q->data)->index) 
				){vflag1=true;}
			if ( (((Vertex*)((Distance*)(array[i])->data)->vB)->index == ((Vertex*)q->data)->index) 
				){vflag2=true;}
			q=q->next;
		}
		if (vflag1&&vflag2){flag=true;}
		else{i++;}
	}
	if (!print_flag) {printf("complete.\n");fflush(stdout);}
	if (!print_flag){ printf("vertices to merge: %i %i, dist = %.15g\n",
				((Vertex*)((Distance*)(array[i])->data)->vA)->index,
				((Vertex*)((Distance*)(array[i])->data)->vB)->index,
				((Distance*)(array[i])->data)->d); }
	if (i==count){printf("Error! No eligible vertices to merge were found.\n");exit(0);}
	else {return i;}
}

void_list * fixVertices(void_list *v,Distance *d,int print_flag) {
	///// remove second vertex in pair from vertex list /////
	void_list *ptr,*p;
	if (!print_flag) {printf("Remove vertex from vertices...");fflush(stdout);}
	ptr=v;
	p=v;
	while (p!=NULL) {
		if(((Vertex*)d->vB)->index == ((Vertex*)p->data)->index){
			// delete vertex
			delete (Vertex*)p->data;
			// adjust list pointer
			if (p==ptr) { ptr=p->next;}
			// remove link from list
			p=removeLink(p);
		}else{p=p->next;}
	}
	// adjust pointer
	v=ptr;
	if (!print_flag) {printf("complete.\n");fflush(stdout);}
	return v;
}

void_list * fixFaces(void_list *f,Distance *d,int print_flag) {
	///// replace all instances of second vertex in face list with first vertex /////
	void_list *p;
	if (!print_flag) {printf("Replace vertex in faces...");fflush(stdout);}
	for (p=f;p!=NULL;p=p->next) {
		if(((Vertex*)d->vB)->index == ((Face*)p->data)->v1){((Face*)p->data)->v1 = ((Vertex*)d->vA)->index;}
		if(((Vertex*)d->vB)->index == ((Face*)p->data)->v2){((Face*)p->data)->v2 = ((Vertex*)d->vA)->index;}
		if(((Vertex*)d->vB)->index == ((Face*)p->data)->v3){((Face*)p->data)->v3 = ((Vertex*)d->vA)->index;}
	}
	if (!print_flag) {printf("complete.\n");fflush(stdout);}
	return f;
}

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

	char *infile;
	char line[2048];
	char *str;
	char *eptr;
	FILE *F;
	void_list *q,*qq,*p,*pp,*prev,*ptr;
	int max_vert,i,j,num_vfree,d_count;
	void_list** dist_array;
	void_list** pairs;

	// linked lists
	void_list *vlh,*vl,*vfreeh,*vfree;
	Vertex *v;
	vlh = NULL;
	vfreeh = NULL;
	void_list *flh,*fl,*ffreeh,*ffree;
	Face *f;
	flh = NULL;
	ffreeh = NULL;
    void_list *edgeh;
    edgeh = NULL;
	void_list *dlh;
	dlh=NULL;

	// pointers to start of vertex and face linked lists 
	void_list *vend;
	void_list *fend;
	void_list *swapend;

	////////// get data /////////
	double epsilon = strtod(argv[2],&eptr);
	int print_flag = (int)strtod(argv[3],&eptr);

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

	// find maximum vertex index
	max_vert = ((Vertex*)vlh->data)->index;
	for (p=vlh;p!=NULL;p=p->next) {
		if (max_vert<((Vertex*)p->data)->index) {max_vert=((Vertex*)p->data)->index;}
	}

    // create array of pointers to vertex list
    void_list* vert_array[max_vert+1];

	// go through linked list backwards and add previous pointers to each vertex
	prev = NULL;
	for (p=vlh;p!=NULL;p=p->next) {
		p->previous = prev;
		prev = p;
		vert_array[((Vertex*)p->data)->index]=p;
		if (p->next==NULL) break;
	}

	// go through linked list backwards and add previous pointers to each face
	prev = NULL;
	for (p=flh;p!=NULL;p=p->next) {
		p->previous = prev;
		prev = p;
		if (p->next==NULL) break;
	}

    ///// build edge list ////
	edgeh = gatherFreeEdges(flh,edgeh,vert_array,print_flag,max_vert,pairs);

	///// gather free vertices /////
	vfreeh = gatherFreeVertices(edgeh,max_vert,print_flag);

	///// compute distances /////
	dist_array = computeDistances(dist_array,vfreeh,flh,d_count,epsilon,print_flag,pairs);

	///// delete pairs /////
	deletePairs(pairs,max_vert);
	
	///// count number free vertices /////
	num_vfree = countFreeVertices(vfreeh,print_flag);

	///// merge vertices /////
	while(num_vfree){

		///// identify free vertex pair separated by the smallest distance /////
		i=findVerticesToMerge(dist_array,d_count,vfreeh,print_flag);

		///// remove second vertex in pair from vertex list /////
		vlh = fixVertices(vlh,(Distance*)dist_array[i]->data,print_flag);

		///// replace all instances of second vertex in face list with first vertex /////
		flh = fixFaces(flh,(Distance*)dist_array[i]->data,print_flag);

		///// delete edge list /////
		edgeh = deleteFreeEdges(edgeh);

	    ///// build edge list ////
		edgeh = gatherFreeEdges(flh,edgeh,vert_array,print_flag,max_vert,pairs);

		///// delete pairs /////
		deletePairs(pairs,max_vert);

		///// delete free vertices /////
		deleteFreeVertices(vfreeh);
		vfreeh=NULL;

		///// gather free vertices /////
		vfreeh = gatherFreeVertices(edgeh,max_vert,print_flag);

		///// count number free vertices /////
		num_vfree = countFreeVertices(vfreeh,print_flag);
	}

	///// cleanup /////
	edgeh = deleteFreeEdges(edgeh);
	deleteFreeVertices(vfreeh);

	////////// add pointers to ends of linked lists /////////
	// read out vertex linked list
	for (q=vlh;q!=NULL;q=q->next) { vend = q; if (q->next==NULL) break; }

	// read out face linked list
	for (q=flh;q!=NULL;q=q->next) { fend = q; if (q->next==NULL) break; }

	////////// write data to stdout /////////
	if(print_flag) {
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


	////////// delete vertices and faces /////////
	p=vlh;
	while (p!=NULL) {
		q=p->next;
		delete (Vertex*)p->data;
		delete p;
		p=q;
	}
	p=flh;
	while (p!=NULL) {
		q=p->next;
		delete (Face*)p->data;
		delete p;
		p=q;
	}
	////////// delete dist_array /////////
/*	void_list** dist_array;
	dist_array = new void_list*[d_count];
	i=0;
    for (q=dlh;q!=NULL;q=q->next) {
		dist_array[i++]=q;
    }
					// create distance
					dl = new void_list();
					dl->next = dlh;
					d = new Distance(squared_dist,q,qq);
					dl->data = (void*)d;
					dlh = dl;
					d_count++;
	p=vlh;
	while (p!=NULL) {
		q=p->next;
		delete (Vertex*)p->data;
		delete p;
		p=q;
	}
*/
	for(i=0;i<d_count;i++){
		delete (Distance*)dist_array[i]->data;
		delete dist_array[i];
	}
	delete[] dist_array;

	return 0;
}
