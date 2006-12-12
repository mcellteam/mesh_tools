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
  Vertex(int i,double xval,double yval,double zval);
};

Vertex::Vertex(int i,double xval,double yval,double zval)
{
	index=i;
	x=xval;
	y=yval;
	z=zval;
}

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

class Edge
{
public:
  int hashval;
  int v1,v2; // vertex indices
  int c12,c21;	// count of times edge was traversed from v1 to v2 and v2 to v1
  void_list *f;	// list of pointers to faces
  Vertex* v; // pointer to new vertex if edge is bisected
  double l; // edge length squared
  bool bisected;
  Edge(void);
  ~Edge(void);
};

Edge::Edge(void){
	hashval=-1;
	l=v1=v2=c12=c21=0;
	f=NULL;
	v=NULL;
	bisected = false;
}

Edge::~Edge(void){
	void_list *p,*q;
    p=f;
    while (p!=NULL) {
        q=p->next;
        delete p;
        p=q;
    }
}

class HashTable{
public:
	void_list **t;	// pointer to void_list *
					// use for array of void_list *
	int s;			// size of array
	int n;			// edges per hash bin
	HashTable(int);
	~HashTable(void);
};

HashTable::HashTable(int a){
	n=4;
	// a = number of edges in block
	int bins = a/n+1;
	int i=0;
	while(bins=bins/2){i++;}
	s=1;
	while(i){s=s*2;i--;}
	t=new void_list*[s];
	for(i=0;i<s;i++){
		t[i]=NULL;
	}
}

HashTable::~HashTable(void){
	void_list *p,*q;
	int i;
	for(i=0;i<s;i++){
	    p=t[i];
	    while (p!=NULL) {
	        q=p->next;
	        delete p;
	        p=q;
	    }
	}
	delete[] t;
}

class EdgeBlock{
public:
	EdgeBlock *next;	// pointer to next block of edges
	Edge *e;			// pointer to edge
	int n;				// number of edges in block
	int p;				// next open edge in block
	HashTable *ht;
	EdgeBlock(int m);
	~EdgeBlock(void);
};

EdgeBlock::EdgeBlock(int m){
	n=m;
	e=new Edge[m];
	next=NULL;
	ht=NULL;
	p=0;
}

EdgeBlock::~EdgeBlock(void){
	delete[] e;
	delete ht;
}

class Face
{
public:
	int index;	// Face index
	int v1,v2,v3;	// vertex indices
	Face(char *triplet);
	Face(int i,int va,int vb,int vc);
	~Face(void);
	void_list *e;	// list of pointers to edges
	bool subdivided;
	bool bisubdivided;
};

Face::~Face(void)
{
    void_list *q,*p;
    q=e;
    while (q!=NULL) {
        p=q->next;
        delete q;
        q=p;
    }
}

Face::Face(int i,int va,int vb,int vc)
{
	subdivided = false;
	bisubdivided = false;
	e=NULL;
	index=i;
	v1=va;
	v2=vb;
	v3=vc;
}

Face::Face(char *triplet)
{
	subdivided = false;
	bisubdivided = false;
	e=NULL;

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
