
typedef  unsigned long int  u4;   /* unsigned 4-byte type */
typedef  unsigned     char  u1;   /* unsigned 1-byte type */

// ######################################
// ######################################
// declarations - so compiler doesn't complain
class Vertex;
class Face;
class Edge;

u4 computeHashValue(int,int);
u4 keyPair(int,int);

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
  std::vector<Vertex*> v;	// adjacent vertices
  std::vector<Edge*> e;		// adjacent edges
  std::vector<Face*> f;		// adjacent faces
  void_list *p; // points to void_list link holding this vertex
  Vertex(char *,void_list*);
  Vertex(int,double,double,double,void_list*);
};

Vertex::Vertex(int i,double xval,double yval,double zval,void_list *q)
{
    index=i;
    x=xval;
    y=yval;
    z=zval;
	p=q;
}


Vertex::Vertex(char *triplet,void_list *q)
{
	char val[80];
	char *eptr;
	int i;

	p=q;

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
	Vertex *v1,*v2,*v3;	// vertex indices
	Edge *e1,*e2,*e3;	// edge pointers
	void_list *p; // points to void_list link holding this face
	Face(char *,void_list*);
    Face(int,Vertex*,Vertex*,Vertex*,void_list*);
	void clearEdges(void);
	bool subdivided;
	bool bisubdivided;
};

void Face::clearEdges(void){
	e1=e2=e3=NULL;
	subdivided=bisubdivided=false;
}

Face::Face(int i,Vertex *va,Vertex *vb,Vertex *vc,void_list *q)
{
    subdivided = false;
    bisubdivided = false;
    e1=e2=e3=NULL;
    index=i;
    v1=va;
    v2=vb;
    v3=vc;
	p=q;
}

Face::Face(char *triplet,void_list *q)
{
	v1=v2=v3=NULL;
	e1=e2=e3=NULL;
    subdivided = false;
    bisubdivided = false;

	char val[80];
	char *eptr;
	int i;

	p=q;

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
		v1=v2=v3=NULL;
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
	v1 = (Vertex*)(int) strtod(val,&eptr);
	if (val==eptr)
	{
		v1=v2=v3=NULL;
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
	v2 = (Vertex*)(int) strtod(val,&eptr);
	if (val==eptr)
	{
		v1=v2=v3=NULL;
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
	v3 = (Vertex*)(int) strtod(val,&eptr);
	if (val==eptr)
	{
		v1=v2=v3=NULL;
		printf("Error in reading vertex index\n");
		return;
	}
}

class Edge
{
public:
  Vertex *v1,*v2; // pointers to edge vertices
  Vertex *va,*vb; // pointers to vertices in adjacent faces not on edge
  Face *f1,*f2;	// pointers to adjacent faces (i.e. faces that contain edge)
  Vertex* v; // pointer to new vertex if edge is bisected
  double l; // edge length squared
  bool bisected;
  Edge(Face*,Vertex*,Vertex*);
  void update(Face*);
  void updateLength(void);
};

void Edge::updateLength(void){
	l=(v1->x-v2->x)*(v1->x-v2->x)+(v1->y-v2->y)*
		(v1->y-v2->y)+(v1->z-v2->z)*(v1->z-v2->z);
}

Edge::Edge(Face *f,Vertex *a,Vertex *b){
	v1=a;
	v2=b;
	f1=f;
	f2=NULL;
	v=NULL;
	bisected=false;
	// compute original edge length
	l=(v1->x-v2->x)*(v1->x-v2->x)+
		(v1->y-v2->y)*(v1->y-v2->y)+
		(v1->z-v2->z)*(v1->z-v2->z);
	// assign va
	if (f->v1!=v1 && f->v1!=v2){va=f->v1;}
	else if (f->v2!=v1 && f->v2!=v2){va=f->v2;}
	else if (f->v3!=v1 && f->v3!=v2){va=f->v3;}
	else {fprintf(stderr,"Error assigning va in edge.\n");exit(1);}
	vb=NULL;
	// add edge* to face*
	if(f->e1==NULL&&((a==f->v1&&b==f->v2)||(a==f->v2&&b==f->v1))){f->e1=this;}
	else if(f->e2==NULL&&((a==f->v2&&b==f->v3)||(a==f->v3&&b==f->v2))){f->e2=this;}
	else if(f->e3==NULL&&((a==f->v3&&b==f->v1)||(a==f->v1&&b==f->v3))){f->e3=this;}
	else {fprintf(stderr,"Error. No available edge* in face.\n");exit(0);}
}

void Edge::update(Face *f){
    //add face to edge
	if(f1==NULL) {f1=f;}
	else if (f2==NULL) {f2=f;}
	else {
		fprintf(stderr,"Error. Tried to add third face to edge.\n");
		fprintf(stderr,"Face %i %i ",f->index,f->v1->index); 
		fprintf(stderr,"%i %i\n",f->v2->index,f->v3->index);
		fprintf(stderr,"Existing Face %i %i ",f1->index,f1->v1->index); 
		fprintf(stderr,"%i %i\n",f1->v2->index,f1->v3->index);
		fprintf(stderr,"Existing Face %i %i ",f2->index,f2->v1->index); 
		fprintf(stderr,"%i %i\n",f2->v2->index,f2->v3->index);
		fprintf(stderr,"Edge vertices: v1 %i, v2 %i\n",v1->index,v2->index);
		char s[32];
		fprintf(stderr,"keyPair %u\n",keyPair(v1->index,v2->index)); 
		exit(1);
	}
	// assign vb
	if (f->v1!=v1 && f->v1!=v2){vb=f->v1;}
	else if (f->v2!=v1 && f->v2!=v2){vb=f->v2;}
	else if (f->v3!=v1 && f->v3!=v2){vb=f->v3;}
	else {fprintf(stderr,"Error assigning vb in edge.\n");exit(1);}
	// add edge* to face*
	if(f->e1==NULL){f->e1=this;}
	else if(f->e2==NULL){f->e2=this;}
	else if(f->e3==NULL){f->e3=this;}
	else {fprintf(stderr,"Error. No available edge* in face.\n");exit(0);}
}

struct equ4
{
  bool operator()(const u4 s1, const u4 s2) const
  {
    return s1==s2;
  }
};

struct pair_hash
{
    u4 operator()(u4 i) const { return i; }
};

typedef __gnu_cxx::hash_map<u4,void_list*,pair_hash,equ4,std::allocator<void_list*> > hashtable_t;
typedef __gnu_cxx::hash_map<u4,void_list*,pair_hash,equ4,std::allocator<void_list*> >::iterator ht_iterator;
typedef __gnu_cxx::hash_map<int,Vertex*> hashtable_v;

