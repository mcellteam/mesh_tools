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
  ~Distance(void);
  bool deleteme;
};

Distance::Distance(double dist,void_list *p,void_list *q)
{
	d = dist;
	vA = (Vertex*)p->data;
	vB = (Vertex*)q->data;
	deleteme = false;
}

Distance::~Distance(void)
{
	if(deleteme){delete vB;}
}

class Edge
{
public:
  int hashval;
//  Vertex *va,*vb; // pointers to vertex links
  int v1,v2; // vertex indices
  int c12,c21;  // count of times edge was traversed from va to vb and vb to va
  void_list *f;
//  Edge(Vertex*,Vertex*);
  Edge(void);
  void addFace(Face*);
  ~Edge(void);
};

Edge::Edge(void){
    hashval=-1;
    v1=v2=c12=c21=0;
    f=NULL;
}

/*Edge::Edge(Vertex *v1,Vertex *v2)
{
	va = v1;
	vb = v2;
	c12 = 1;
	c21 = 0;
	f=NULL;
}
*/
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

class HashTable{
public:
    void_list **t;  // pointer to void_list *
                    // use for array of void_list *
    int s;          // size of array
    int n;          // edges per hash bin
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
    EdgeBlock *next;    // pointer to next block of edges
    Edge *e;            // pointer to edge
    int n;              // number of edges in block
    int p;              // next open edge in block
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

class Data
{
public:
  Edge *ptr; // pointer to edge
  int vert; // vertex of pair 
  Data(Edge*,int);
};

Data::Data(Edge *e,int v)
{
	ptr = e;
	vert = v;
}
