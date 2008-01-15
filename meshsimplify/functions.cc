// #####################################################
// #####################################################

typedef  unsigned long int  u4;   /* unsigned 4-byte type */
typedef  unsigned     char  u1;   /* unsigned 1-byte type */

/* The mixing step */
#define mix(a,b,c) \
{ \
  a=a-b;  a=a-c;  a=a^(c>>13); \
  b=b-c;  b=b-a;  b=b^(a<<8);  \
  c=c-a;  c=c-b;  c=c^(b>>13); \
  a=a-b;  a=a-c;  a=a^(c>>12); \
  b=b-c;  b=b-a;  b=b^(a<<16); \
  c=c-a;  c=c-b;  c=c^(b>>5);  \
  a=a-b;  a=a-c;  a=a^(c>>3);  \
  b=b-c;  b=b-a;  b=b^(a<<10); \
  c=c-a;  c=c-b;  c=c^(b>>15); \
}

/* The whole new hash function */
u4 hash( register u1* k, u4 length, u4 initval)
//register u1 *k;        /* the key */
//u4           length;   /* the length of the key in bytes */
//u4           initval;  /* the previous hash, or an arbitrary value */
{
   register u4 a,b,c;  /* the internal state */
   u4          len;    /* how many key bytes still need mixing */

   /* Set up the internal state */
   len = length;
   a = b = 0x9e3779b9;  /* the golden ratio; an arbitrary value */
   c = initval;         /* variable initialization of internal state */

   /*---------------------------------------- handle most of the key */
   while (len >= 12)
   {
      a=a+(k[0]+((u4)k[1]<<8)+((u4)k[2]<<16) +((u4)k[3]<<24));
      b=b+(k[4]+((u4)k[5]<<8)+((u4)k[6]<<16) +((u4)k[7]<<24));
      c=c+(k[8]+((u4)k[9]<<8)+((u4)k[10]<<16)+((u4)k[11]<<24));
      mix(a,b,c);
      k = k+12; len = len-12;
   }

   /*------------------------------------- handle the last 11 bytes */
   c = c+length;
   switch(len)              /* all the case statements fall through */
   {
   case 11: c=c+((u4)k[10]<<24);
   case 10: c=c+((u4)k[9]<<16);
   case 9 : c=c+((u4)k[8]<<8);
      /* the first byte of c is reserved for the length */
   case 8 : b=b+((u4)k[7]<<24);
   case 7 : b=b+((u4)k[6]<<16);
   case 6 : b=b+((u4)k[5]<<8);
   case 5 : b=b+k[4];
   case 4 : a=a+((u4)k[3]<<24);
   case 3 : a=a+((u4)k[2]<<16);
   case 2 : a=a+((u4)k[1]<<8);
   case 1 : a=a+k[0];
     /* case 0: nothing left to add */
   }
   mix(a,b,c);
   /*-------------------------------------------- report the result */
   return c;
}

// #####################################################
// #####################################################

u4 computeHashValue(int key1, int key2) {
	char cval[10];
	u4 result;
	sprintf(cval,"%i",key1);
	result = hash((u1*)cval,(u4)strlen(cval),(u4)key2);
	return result;
}

// #####################################################
// #####################################################

void getData(char *infile,void_list *&flh,void_list *&vlh){

	char line[2048],*str,*eptr;
	FILE *F;
	void_list *vl;
	Vertex *v;
	void_list *fl;
	Face *f;

	// open first file
	F = fopen(infile,"r");
	if (!F) {
	printf("Couldn't open input file %s\n",infile);
	return;
	}

	// for every line in file
	for (str=fgets(line,2048,F) ; str!=NULL ; str=fgets(line,2048,F)) {

		// skip leading whitespace
		while (strchr(" \t,",*str)!=NULL) { str++;}

		// if first character is V for Vertex, add new linked list class instance
		if (strchr("V",*str)!=NULL){
			vl = new void_list();
			vl->next = vlh;
			v = new Vertex(str,vl);
			vl->data = (void*)v;
			vlh = vl;
		} 
		// if first character is F for Face, add new linked list class instance
		else if (strchr("F",*str)!=NULL){
			fl = new void_list();
			fl->next = flh;
			f = new Face(str,fl);
			fl->data = (void*)f;
			flh = fl;
		}
	}
	fclose(F);
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

// #####################################################
// #####################################################

int maxVert(void_list *L){
    void_list *p;
    int max=0;
    for (p=L;p!=NULL;p=p->next) {
        if (max<((Vertex*)p->data)->index) {max=((Vertex*)p->data)->index;}
    }
    return max;
}

int maxFace(void_list *L){
    void_list *p;
    int max=0;
    for (p=L;p!=NULL;p=p->next) {
        if (max<((Face*)p->data)->index) {max=((Face*)p->data)->index;}
    }
    return max;
}

// #####################################################
// #####################################################

void addPointersToFaces(void_list *flh,hashtable_v &hmv){
	void_list *q;
	Face *f;
	int i=0;
	// for each face
	for (q=flh;q!=NULL;q=q->next) {
		f=(Face*)q->data;
		f->v1=hmv[(int)f->v1];
		f->v2=hmv[(int)f->v2];
		f->v3=hmv[(int)f->v3];
	}
}

void buildVertMap(void_list *vlh,hashtable_v &hmv){
	void_list *p;
	Vertex *v;
	// for each vertex
	for (p=vlh;p!=NULL;p=p->next){
		v=(Vertex*)p->data;
	    hmv[v->index]=v;
	}
}

// #####################################################
// #####################################################

u4 keyPair(int a,int b){
	if (a<b){ return computeHashValue(a,b);}
	else {return computeHashValue(b,a);}
}

bool edgeMatch(Edge *e,int va,int vb) {
	if ( (e->v1->index==va && e->v2->index==vb) ||
		(e->v1->index==vb && e->v2->index==va) ){return true;}
	else {return false;}
}

void checkEdge(hashtable_t &hme,Face *f,Vertex *va,Vertex *vb,std::vector<Edge*> &e,int j) {
	char s[32];
	Edge *ee=NULL,*et=NULL,*en=NULL;
	void_list *p;
	///// find edge /////
	// if element exists given key
	if (hme.count(keyPair(va->index,vb->index))>0){
		///// get Edge pointer /////
		// for each link in list
		for(p=hme[keyPair(va->index,vb->index)];p!=NULL;p=p->next){
			// does edge match
			et=(Edge*)p->data;
			if (edgeMatch(et,va->index,vb->index)){ee=et;}
		}
	}
	///// if exists, then update /////
	if(ee!=NULL){ ee->update(f); }
	else { ///// edge pointer not found /////
		// new edge
		en = new Edge(f,va,vb);
		// add edge to vector
		e.push_back(en);
		// create new link in void_list
		p = new void_list();
		p->data=(void*)en;
		p->previous=NULL;
		// hash_map empty
		if (et==NULL){
			p->next=NULL;
		} else { // hash_map not empty 
			p->next=hme[keyPair(va->index,vb->index)];
			hme[keyPair(va->index,vb->index)]->previous=p;
		}
	    // store edge pointer in hash table
	    hme[keyPair(va->index,vb->index)]=p;
	}
}

void getEdges(void_list *flh,hashtable_t &hme,std::vector<Edge*> &e,int j){
	void_list *q;
	Face *f;
	int i=0;
	// for each face
	for (q=flh;q!=NULL;q=q->next) {
		f=(Face*)q->data;
		checkEdge(hme,f,f->v1,f->v2,e,j);
		checkEdge(hme,f,f->v2,f->v3,e,j);
		checkEdge(hme,f,f->v3,f->v1,e,j);
	}
}

// #####################################################
// #####################################################

void clearVertexAdjacencies(void_list *vlh){
	Vertex *v;
	// for each vertex
	for (void_list *p=vlh;p!=NULL;p=p->next) {
		v=(Vertex*)p->data;
		v->f.clear();
		v->v.clear();
		v->e.clear();
	}
}

void findVertexAdjacencies(void_list *flh,int j,std::vector<Edge*> &e){
	Face *ff;
	std::vector<Edge*>::iterator i;
	// for each face, add face* to each face vertex
	for (void_list *p=flh;p!=NULL;p=p->next) {
		ff=(Face*)p->data;
		ff->v1->f.push_back(ff);
		ff->v2->f.push_back(ff);
		ff->v3->f.push_back(ff);
	}
	// for each edge
	for(i=e.begin();i!=e.end();i++){
		// add each vertex to other as adjacent
		(*i)->v1->v.push_back((*i)->v2);
		(*i)->v2->v.push_back((*i)->v1);
		// add edge* to each vertex
		(*i)->v1->e.push_back(*i);
		(*i)->v2->e.push_back(*i);
	}
}

// #####################################################
// #####################################################

int threshold(Edge *e,double t){return(e->l<t);}

Edge* findBadEdge(std::vector<Edge*> &e,double t){
	std::vector<Edge*>::iterator i;
	// for each edge
	i=e.begin();
	while (i!=e.end()) {
		if (threshold(*i,t)) {return (*i);}
		i++;
	}
	return NULL;
}

// #####################################################
// #####################################################

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

void deleteVertices(Edge *e,void_list *&vlh){
	fprintf(stderr,"deleteVertices: v1 %i, v2 %i\n",e->v1->index,e->v2->index);
	void_list *q=e->v1->p;
	// adjust list pointer
	if (q==vlh) {vlh=q->next;}
	// remove vertex from linked list
	q=removeLink(q);
	// assign pointer
	q=e->v2->p;
	// adjust list pointer
	if (q==vlh) {vlh=q->next;}
	// remove vertex from linked list
	q=removeLink(q);
}


void deleteFaces(Edge *e,void_list *&flh){
	void_list *p;
	if (e->f1==NULL){fprintf(stderr,"deleteFaces: f1 is NULL\n");exit(0);}
	if (e->f2==NULL){fprintf(stderr,"deleteFaces: f2 is NULL\n");exit(0);}
	if(e->f1!=NULL){
		fprintf(stderr,"deleteFace: Face %i %i %i %i was deleted.\n",
				e->f1->index,e->f1->v1->index,e->f1->v2->index,e->f1->v3->index);
		///// delete first face /////
		// assign pointer
		p=e->f1->p;
		// adjust list pointer
		if (p==flh) {flh=p->next;}
		// remove face from linked list
		p=removeLink(p);
	}
	else {fprintf(stderr,"deleteFaces: f1 is NULL\n");exit(0);}
	if (e->f2!=NULL){
		fprintf(stderr,"deleteFace: Face %i %i %i %i was deleted.\n",
				e->f2->index,e->f2->v1->index,e->f2->v2->index,e->f2->v3->index);
		///// delete second face /////
		// assign pointer
		p=e->f2->p;
		// adjust list pointer
		if (p==flh) {flh=p->next;}
		// remove face from linked list
		p=removeLink(p);
	}
	else {fprintf(stderr,"deleteFaces: f2 is NULL\n");exit(0);}
}

void fixFaces(Vertex *nv,Vertex *ov){
	///// fix all adjacent faces of old vertex //////
	std::vector<Face*>::iterator i;
	// for each face of edge vertex v
	for(i=ov->f.begin();i!=ov->f.end();i++){
		// replace reference to either old edge vertex
		// with reference to new vertex
		if((*i)->v1->index==ov->index||(*i)->v1->index==ov->index){(*i)->v1=nv; 
          fprintf(stderr,"fixFaces: Face %i %i %i %i fixed.\n",
                  (*i)->index,(*i)->v1->index,(*i)->v2->index,(*i)->v3->index);
		}
		else if((*i)->v2->index==ov->index||(*i)->v2->index==ov->index){(*i)->v2=nv; 
          fprintf(stderr,"fixFaces: Face %i %i %i %i fixed.\n",
                  (*i)->index,(*i)->v1->index,(*i)->v2->index,(*i)->v3->index);
		}
		else if((*i)->v3->index==ov->index||(*i)->v3->index==ov->index){(*i)->v3=nv; 
          fprintf(stderr,"fixFaces: Face %i %i %i %i fixed.\n",
                  (*i)->index,(*i)->v1->index,(*i)->v2->index,(*i)->v3->index);
		}
	}
}

Vertex* addVertex(Edge *e,void_list *&vlh,int &max_verts){
	// compute new vertex
	double x=(e->v1->x+e->v2->x)/2.0;
	double y=(e->v1->y+e->v2->y)/2.0;
	double z=(e->v1->z+e->v2->z)/2.0;
	// create new vertex
	void_list *p = new void_list();
	vlh->previous=p;
	p->next = vlh;
	p->previous = NULL;
	Vertex *v = new Vertex(++max_verts,x,y,z,p);
	p->data = (void*)v;
	vlh = p;
	fprintf(stderr,"addVertex: %i %.15g %.15g %.15g\n",max_verts,x,y,z);
	return v;
}

void check(void_list *vlh){
	void_list *q;
	Face *f;
	Vertex *v;
	// for each vertex
	for (q=vlh;q!=NULL;q=q->next) {
		v=(Vertex*)q->data;
		if(v->index<0){fprintf(stderr,"index is negative\n");exit(1);}
		if(v->p==NULL){fprintf(stderr,"p is NULL\n");exit(1);}
	}
}

bool equalSums(Face *lhs,Face *rhs){
	int L1=lhs->v1->index,L2=lhs->v2->index,L3=lhs->v3->index;
	int R1=rhs->v1->index,R2=rhs->v2->index,R3=rhs->v3->index;
	return ((L1+L2+L3)==(R1+R2+R3));
}

bool identicalFaces(Face *lhs,Face *rhs){
	int L1=lhs->v1->index,L2=lhs->v2->index,L3=lhs->v3->index;
	int R1=rhs->v1->index,R2=rhs->v2->index,R3=rhs->v3->index;
	return ( (L1==R1||L1==R2||L1==R3)&&
			(L2==R1||L2==R2||L2==R3)&&
			(L3==R1||L3==R2||L3==R3));
}

bool sharedVertex(Vertex *n,Edge *e,Vertex *v){
	// if vertex is not on edge and not new26 vertex
	if(v->index!=e->v1->index && v->index!=e->v2->index && v->index!=n->index){
		// if vertex is on an edge face
		if( (v->index==e->f1->v1->index) ||
			(v->index==e->f1->v2->index) ||
			(v->index==e->f1->v3->index) ||
			(v->index==e->f2->v1->index) ||
			(v->index==e->f2->v2->index) ||
			(v->index==e->f2->v3->index)){return true;}
		else {return false;}
	} else {return false;}
}

void collectEdgePtrs(Vertex *t,Edge *e,std::vector<Edge*> &ep){
	std::vector<Edge*>::iterator i;
	// for each edge* of third vertex
	for (i=t->e.begin();i!=t->e.end();i++){
		// if edge includes either of target edge vertices, then keep
		if ( ((*i)->v1->index==e->v1->index) ||
			((*i)->v1->index==e->v2->index) ||
			((*i)->v2->index==e->v1->index) ||
			((*i)->v2->index==e->v2->index)
			){ep.push_back(*i); }
	}
	// check
	if (ep.size()!=2){fprintf(stderr,"collectEdgePtrs: ERROR. Wrong number of edges - %i\n",ep.size());exit(0);}
}

void addFace(Edge *e,int v1,int v2,std::vector<Face*> &up,std::vector<Face*> &dn){
	int va,vb,vc;
	Face *f=e->f1;
	va=f->v1->index;
	vb=f->v2->index;
	vc=f->v3->index;
	if ( (va==v1 && vb==v2) ||
		(vb==v1 && vc==v2) ||
		(vc==v1 && va==v2)
		){up.push_back(f);dn.push_back(e->f2);}
	else{dn.push_back(f);up.push_back(e->f2);}
}

void addFacesToUpDn(Edge *e,std::vector<Face*> &up,std::vector<Face*> &dn,std::vector<Edge*> &ep){
	std::vector<Edge*>::iterator i;
	int v1,v2;
	// process e->v1,e->v2 edge
	addFace(e,e->v1->index,e->v2->index,up,dn);
	// for each edge of bounding triangle
	for (i=ep.begin();i!=ep.end();i++){
		// if edge is e->v2,t edge
		if ((*i)->v1->index==e->v2->index || (*i)->v2->index==e->v2->index) {
			// assign edge vertices
			if ((*i)->v1->index==e->v2->index){v1=(*i)->v1->index;v2=(*i)->v2->index;}
			else {v2=(*i)->v1->index;v1=(*i)->v2->index;}
			addFace(*i,v1,v2,up,dn);
		} else {
			// edge is t,e->v1 edge
			// assign edge vertices
			if ((*i)->v1->index==e->v1->index){v2=(*i)->v1->index;v1=(*i)->v2->index;}
			else {v1=(*i)->v1->index;v2=(*i)->v2->index;}
			addFace(*i,v1,v2,up,dn);
		}
	}
}

void centroid(Vertex *v1,Vertex *v2,Vertex *v3,std::vector<double> &c){
	c.push_back((v1->x+v2->x+v3->x)/3.0);
	c.push_back((v1->y+v2->y+v3->y)/3.0);
	c.push_back((v1->z+v2->z+v3->z)/3.0);
}

double computeSqDist(std::vector<Face*> f,std::vector<double> c){
	double a=0;
	std::vector<double> o;
	std::vector<Face*>::iterator i;
	// for each face in set (either up or dn)
	for (i=f.begin();i!=f.end();i++){
		// compute centroid of each face
		centroid((*i)->v1,(*i)->v2,(*i)->v3,o);
		// compute sum of squared distance between face centroid and vertices' centroid
		a += (c[0]-o[0])*(c[0]-o[0])+(c[1]-o[1])*(c[1]-o[1])+(c[2]-o[2])*(c[2]-o[2]);
	}
	return a;
}

bool vNotInSet(Vertex *v,Vertex *a,Vertex *b,Vertex *c){
	return v->index!=a->index && v->index!=b->index && v->index!=c->index;
}

struct eSortComp
{
	bool operator()(const Edge *lhs,const Edge *rhs)
	{
/*		if(lhs->v1!=rhs->v1){return (lhs->v1->index>rhs->v1->index);}
		// then v1s are equal
		else if(lhs->v2!=rhs->v2){return (lhs->v2->index>rhs->v2->index);}
		// then v1s and v2s are equal, so return false
		else {return false;}
*/
		return (lhs<rhs);
	}
};

struct vSortComp
{
	bool operator()(const Vertex *lhs,const Vertex *rhs)
	{
		//return (lhs->index < rhs->index);
		return (lhs<rhs);
	}
};

struct fSortComp
{
	bool operator()(const Face *lhs,const Face *rhs)
	{
		//return (lhs->index>rhs->index);
		return (lhs<rhs);
	}
};

inline bool vUniqComp(Vertex *a,Vertex* b) { return a->index == b->index; }
inline bool fUniqComp(Face *a,Face* b) { return a->index == b->index; }
inline bool eUniqComp(Edge *a,Edge* b) { return ((a->v1==b->v1)&&(a->v2==b->v2)); }

void addVertices(std::vector<Vertex*> &v,std::vector<Face*> &f,Edge *e,Vertex *t){
	std::vector<Face*>::iterator i;
	std::vector<Vertex*>::iterator j;
	// for each inside face
	for(i=f.begin();i!=f.end();i++){
		// add face vertex that is not part of bounding triangle to interior vertices
		if (vNotInSet((*i)->v1,e->v1,e->v2,t)){v.push_back((*i)->v1);}
		else if (vNotInSet((*i)->v2,e->v1,e->v2,t)){v.push_back((*i)->v2);}
		else if (vNotInSet((*i)->v3,e->v1,e->v2,t)){v.push_back((*i)->v3);}
		else {fprintf(stderr,"addVertices: ERROR. No face vertex is not in set.\n");exit(0);}
	}
	// sort vertices
	sort(v.begin(),v.end(),vSortComp());
	// keep unique set
	j = unique(v.begin(),v.end(),vUniqComp);
	v.assign(v.begin(),j);
}

bool getInteriorVerts(std::vector<Vertex*> &iv,Edge *e,Vertex *t){
	std::vector<Vertex*> av;
	std::vector<Vertex*>::iterator i;
	int count=0;
	do{
		if(av.size()>0){
			// save as interior vertex set
			iv.assign(av.begin(),av.end());
			av.clear();
		}
		// for each interior vertex, collect adjacent vertices
		for(i=iv.begin();i!=iv.end();i++){
			av.push_back(*i);
			av.insert(av.end(),(*i)->v.begin(),(*i)->v.end());
		}
		// sort vertices
		sort(av.begin(),av.end(),vSortComp());
		// keep unique set
		i = unique(av.begin(),av.end(),vUniqComp);
		av.assign(av.begin(),i);
		//remove bounding triangle vertices
		i=av.begin();
		while(i!=av.end()){
			if ((*i)->index==e->v1->index || (*i)->index==e->v2->index || (*i)->index==t->index ){i=av.erase(i);}
			else {i++;}
		}
//		fprintf(stderr,"getInteriorVerts: iv.size() %i, av.size() %i\n",iv.size(),av.size());
	} while (iv.size()!=av.size()&&count++<MAX_ITER);
//	fprintf(stderr,"getInteriorVerts: count %i.\n",count);
	if ((count-1)<MAX_ITER){return false;}
	else {return true;}
}

void deleteVertsFaces(std::vector<Vertex*> &iv,void_list *&vlh,void_list *&flh,Edge *e){
	std::vector<Vertex*>::iterator i;
	std::vector<Face*>::iterator j;
	std::vector<Face*> f;
	void_list *p;
	///// collect associated faces and delete vertices /////
	// for each interior vertex
	for(i=iv.begin();i!=iv.end();i++){
		// add adjacent faces to vector
		f.insert(f.end(),(*i)->f.begin(),(*i)->f.end());
		///// delete vertex /////
		fprintf(stderr,"deleteVertsFaces: removing vertex %i\n",(*i)->index);
	}
	///// delete associated faces /////
	// sort faces
	sort(f.begin(),f.end(),fSortComp());
	// keep unique
	j = unique(f.begin(),f.end(),fUniqComp);
	f.assign(f.begin(),j);
	// for each adjacent face
//	for(j=f.begin();j!=f.end();j++){
//		fprintf(stderr,"deleteVertsFaces: face %i\n",(*j)->index);
//	}
	// for each adjacent face
	for(j=f.begin();j!=f.end();j++){
		// if face is not associated with edge
		if((*j)!=e->f1&&(*j)!=e->f2){
			///// delete face /////
			fprintf(stderr,"deleteVertsFaces: removing face %i\n",(*j)->index);
			// assign pointer
			p=(*j)->p;
			// adjust list pointer
			if (p==flh) {flh=p->next;}
			// remove face from linked list
			p=removeLink(p);
		}
	}
}

void gatherNewFaces(std::vector<Edge*> &ec,std::vector<Face*> &nf){
	std::vector<Edge*>::iterator i;
	std::vector<Face*>::iterator j;
	// for each edge
	for(i=ec.begin();i!=ec.end();i++){
		// add edge-associated faces to face collection
		nf.push_back((*i)->f1);
		nf.push_back((*i)->f2);
	}
	// sort faces and keep unique
	sort(nf.begin(),nf.end(),fSortComp());
	j = unique(nf.begin(),nf.end(),fUniqComp);
	nf.assign(nf.begin(),j);
}

void gatherEdges(std::vector<Face*> &nf,std::vector<Edge*> &ne){
	std::vector<Edge*>::iterator i;
	std::vector<Face*>::iterator j;
	// for each face
	for(j=nf.begin();j!=nf.end();j++){
		// add face-associated edges
		ne.push_back((*j)->e1);
		ne.push_back((*j)->e2);
		ne.push_back((*j)->e3);
	}
	// sort edges and keep unique
	sort(ne.begin(),ne.end(),eSortComp());
	///// remove any edges that occur twice /////
	// for each edge
	i=ne.begin();
	while((i+1)!=ne.end()){
		// if sequential pair is identical, then delete them
		if(*i==*(i+1)){
//			fprintf(stderr,"gatherEdges: identical edge pair to be deleted - v1 %i, v2 %i\n",(*i)->v1->index,(*i)->v2->index);
			i=ne.erase(i);
//			fprintf(stderr,"gatherEdges: identical edge pair to be deleted - v1 %i, v2 %i\n",(*i)->v1->index,(*i)->v2->index);
			i=ne.erase(i);
		}else{i++;}
		if(i==ne.end()){break;}
	}
}

void keepNewEdges(std::vector<Edge*> &ec,std::vector<Edge*> &ne){
	std::vector<Edge*>::iterator i;
	std::pair<std::vector<Edge*>::iterator,std::vector<Edge*>::iterator> p;

//	for(i=ec.begin();i!=ec.end();i++){
//		fprintf(stderr,"keepNewEdges: collection edge- v1 %i, v2 %i\n",(*i)->v1->index,(*i)->v2->index);
//	}

//	for(i=ne.begin();i!=ne.end();i++){
//		fprintf(stderr,"keepNewEdges: BEFORE new edge- v1 %i, v2 %i\n",(*i)->v1->index,(*i)->v2->index);
//	}

	// for each edge in collection
	for(i=ec.begin();i!=ec.end();i++){
		// search for edge in new edge* set
		p=equal_range(ne.begin(),ne.end(),*i);
		// if found
		if(p.first!=p.second){
			// remove from new edge* set
			ne.erase(p.first);
		} 
		//else {fprintf(stderr,"growToThird: Error. Edge from collection not found in new edge set.\n");exit(0);}
	}

//	for(i=ne.begin();i!=ne.end();i++){
//		fprintf(stderr,"keepNewEdges: AFTER new edge- v1 %i, v2 %i\n",(*i)->v1->index,(*i)->v2->index);
//	}
}

void removeThirdEdges(Triangle &tri,std::vector<Edge*> &ne,Vertex *&v){
	std::vector<Edge*>::iterator i;
	std::pair<std::vector<Edge*>::iterator,std::vector<Edge*>::iterator> p;

//	std::vector<Vertex*>::iterator j;
//	for(j=tri.t.begin();j!=tri.t.end();j++){
//		fprintf(stderr,"removeThirdEdges: third vertex %i\n",(*j)->index);
//	}

	// for each new edge
	i=ne.begin();
	while(i!=ne.end()){
		if ( // if either new edge vertex is in third vertex vector
			(count(tri.t.begin(),tri.t.end(),(*i)->v1)!=0 ||
			count(tri.t.begin(),tri.t.end(),(*i)->v2)!=0) &&
			// and either new edge vertex is on Edge *e
			((*i)->v1==tri.e->v1 || (*i)->v1==tri.e->v2 || (*i)->v2==tri.e->v1 || (*i)->v2==tri.e->v2)
			){
			// save third vertex
			if((*i)->v1==tri.e->v1 || (*i)->v1==tri.e->v2){v=(*i)->v2;}
			else if((*i)->v2==tri.e->v1 || (*i)->v2==tri.e->v2){v=(*i)->v1;}
			else {fprintf(stderr,"removeThirdEdges: Error. No third vertex found.\n");exit(0);}
			// then remove new edge from vector
			i=ne.erase(i);
		} else {i++;}
	}
	// sort new edges
	sort(ne.begin(),ne.end(),eSortComp());
	///// remove target edge, e /////
	p=equal_range(ne.begin(),ne.end(),tri.e);
	// if found
	if(p.first!=p.second){
		// remove from new edge* set
		ne.erase(p.first);
	}
}

void keepNewFaces(std::vector<Face*> &nf,std::vector<Face*> &fc){
	std::vector<Face*>::iterator i;
	std::pair<std::vector<Face*>::iterator,std::vector<Face*>::iterator> p;

//	for(i=fc.begin();i!=fc.end();i++){
//		fprintf(stderr,"keepNewFaces: collection face %i, ptr %i\n",(*i)->index,*i);
//	}

//	for(i=nf.begin();i!=nf.end();i++){
//		fprintf(stderr,"keepNewFaces: BEFORE new face %i, ptr %i\n",(*i)->index,*i);
//	}

	// sort new edges
	sort(nf.begin(),nf.end(),fSortComp());
	// for each face in collection
	for(i=fc.begin();i!=fc.end();i++){
		// search for face in new face* set
		p=equal_range(nf.begin(),nf.end(),*i);
		// if found
		if(p.first!=p.second){
//			fprintf(stderr,"keepNewFaces: erasing face %i\n",(*p.first)->index);
			// remove from new face* set
			nf.erase(p.first);
		}  
		//else {fprintf(stderr,"growToThird: Error. Edge from collection not found in new edge set.\n");exit(0);}
//		else {fprintf(stderr,"keepNewFaces: Insertion point for face %i in new face set %i.\n",(*i)->index,(*p.first)->index);exit(0);}
	}

//	for(i=nf.begin();i!=nf.end();i++){
//		fprintf(stderr,"keepNewFaces: AFTER new face %i\n",(*i)->index);
//	}
}

bool growToThird(Wave &w,Triangle &tri){
	std::vector<Face*> nf;	// new face* collection
	std::vector<Edge*> ne;	// new edge* collection
	std::vector<Edge*>::iterator i;
	std::vector<Face*>::iterator j;
	Vertex *v=NULL;
	int z=0;
	// while new edges (i.e. third vertex not reached)
	// and iterations < MAX_ITER ()
	while(!w.ec.empty()&&(z++<MAX_ITER)){
		// clear vectors
		nf.clear();
		ne.clear();
		// gather new faces from edge collection
		gatherNewFaces(w.ec,nf);
		// remove new faces found in face collection
		keepNewFaces(nf,w.fc);

//		for(j=nf.begin();j!=nf.end();j++){
//			fprintf(stderr,"growToThird: new face %i\n",(*j)->index);
//		}
//		fprintf(stderr,"growToThird: iter %i, # new faces %i\n",z,nf.size());

		// gather edges from new faces
		gatherEdges(nf,ne);

//		for(i=ne.begin();i!=ne.end();i++){
//			fprintf(stderr,"growToThird: new edge- v1 %i, v2 %i\n",(*i)->v1->index,(*i)->v2->index);
//		}
//		fprintf(stderr,"growToThird: iter %i, # new edges %i\n",z,ne.size());

		// keep edges not in edge collection
		keepNewEdges(w.ec,ne);
//		fprintf(stderr,"growToThird: iter %i, # kept edges %i\n",z,ne.size());
		// remove e->v1/3rd or e->v2/3rd, if any, and target edge, e
		removeThirdEdges(tri,ne,v);
//		fprintf(stderr,"growToThird: iter %i,after removing third edges, # edges %i\n",z,ne.size());
		// assign new edges to edge collection
		w.ec.assign(ne.begin(),ne.end());
		// save faces
		w.fc.insert(w.fc.end(),nf.begin(),nf.end());
	}
	// if third vertex was found
	if(w.ec.empty()&&z<MAX_ITER){
		if(v!=NULL){
			fprintf(stderr,"growToThird: iter %i,third vertex was found, v %i\n",z,v->index);
		}else{
			fprintf(stderr,"growToThird: ERROR. iter %i,third vertex was found, but v=NULL\n",z);
			exit(0);
		}
		// remove third vertex from list
		tri.t.erase(find(tri.t.begin(),tri.t.end(),v));
		///// load edge collection with third vertex edges /////
		// for each adjacent edge of v
		for(i=v->e.begin();i!=v->e.end();i++){
			// if edge includes e->v1 or e->v2
			if ((*i)->v1==tri.e->v1||(*i)->v1==tri.e->v2||(*i)->v2==tri.e->v1||(*i)->v2==tri.e->v2){
				// add edge to ec
				w.ec.push_back(*i);
			}
		}
		// save third vertex
		if(!tri.flag){tri.v1=v;}
		else{tri.v2=v;}
		// save faces to collection
//		for(j=w.fc.begin();j!=w.fc.end();j++){
//			fprintf(stderr,"growToThird: added to collection- face %i\n",(*j)->index);
//		}
		tri.f.insert(tri.f.end(),w.fc.begin(),w.fc.end());
		//return
		return true;
	}
	fprintf(stderr,"growToThird: iter %i,NO third vertex was found\n",z);
	return false;
}

//void initializeEdgesFaces(Edge *e,Vertex *v,std::vector<Edge*> &ec,std::vector<Face*> &fc){
void initializeEdgesFaces(Vertex *v,Wave &w,Triangle &tri){
	Face *f;
	// identify correct edge face
	if(v==tri.e->va){f=tri.e->f1;}
	else if(v==tri.e->vb){f=tri.e->f2;}
	else {fprintf(stderr,"Error. Vertex matches neither e->va nor e->vb\n");exit(0);}
	// add edges e->v1/v and e->v2/v to edge collection
	if(f->e1!=tri.e){w.ec.push_back(f->e1);}
	if(f->e2!=tri.e){w.ec.push_back(f->e2);}
	if(f->e3!=tri.e){w.ec.push_back(f->e3);}
	// add face to new face
	w.fc.push_back(f);
}

void growSide(Vertex *v,Triangle &tri){
	Wave w;
	// add edges e->v1/v and e->v2/v to edge collection where v=e->va or v=e->vb
	initializeEdgesFaces(v,w,tri);

//	std::vector<Edge*>::iterator j;
//	for(j=w.ec.begin();j!=w.ec.end();j++){
//		fprintf(stderr,"growSide: initialized edge- v1 %i, v2 %i\n",(*j)->v1->index,(*j)->v2->index);
//	}

	// while growing finds third vertices
	while(growToThird(w,tri)){if(tri.t.empty())break;}
}

void findTriangles(Triangle &tri){
	// search e->va side
	fprintf(stderr,"findTriangles: growSide e->va\n");
	growSide(tri.e->va,tri);
	// if more third vertices to find, then search e->vb side
	if (!tri.t.empty()) {
		fprintf(stderr,"findTriangles: growSide e->vb\n");
		tri.flag=true;
		growSide(tri.e->vb,tri);
	}
	// WHEN finished all third vertices should have been found
	if (!tri.t.empty()) {fprintf(stderr,"findTriangles: Error. Not all third vertices found\n");exit(1);}
}

void Triangle::gatherTriangleVerts(void){
	std::vector<Face*>::iterator i;
	std::vector<Vertex*>::iterator j;
	std::pair<std::vector<Vertex*>::iterator,std::vector<Vertex*>::iterator> p;
	// for each face
	for(i=f.begin();i!=f.end();i++){
		// gather face verts
		v.push_back((*i)->v1);
		v.push_back((*i)->v2);
		v.push_back((*i)->v3);
	}
	// sort vertices and keep unique
	sort(v.begin(),v.end(),vSortComp());
	j = unique(v.begin(),v.end(),vUniqComp);
	v.assign(v.begin(),j);
	///// delete edge vertices /////
	p=equal_range(v.begin(),v.end(),e->v1);
	// if found then remove
	if(p.first!=p.second){v.erase(p.first);}
	p=equal_range(v.begin(),v.end(),e->v2);
	// if found then remove
	if(p.first!=p.second){v.erase(p.first);}
	///// delete third vertices /////
	p=equal_range(v.begin(),v.end(),v1);
	// if found then remove
	if(p.first!=p.second){v.erase(p.first);}
	p=equal_range(v.begin(),v.end(),v2);
	// if found then remove
	if(p.first!=p.second){v.erase(p.first);}
}

void Triangle::deleteTriangleVerts(void_list *&vlh){
	std::vector<Vertex*>::iterator i;
	void_list *q;
	// for each vertex
	for(i=v.begin();i!=v.end();i++){
		fprintf(stderr,"deleteTriangleVerts: deleting vertex %i\n",(*i)->index);
		// assign pointer
		q=(*i)->p;
		// adjust list pointer
		if (q==vlh) {vlh=q->next;}
		// remove vertex from linked list
		q=removeLink(q);
	}
}

void Triangle::deleteTriangleFaces(void_list *&flh){
	std::vector<Face*>::iterator i;
	void_list *q;
	std::pair<std::vector<Face*>::iterator,std::vector<Face*>::iterator> p;
	// sort vertices and keep unique
	sort(f.begin(),f.end(),fSortComp());
	i = unique(f.begin(),f.end(),fUniqComp);
	f.assign(f.begin(),i);
	///// remove edge faces if exist /////
	p=equal_range(f.begin(),f.end(),e->f1);
	// if found then remove
	if(p.first!=p.second){f.erase(p.first);}
	p=equal_range(f.begin(),f.end(),e->f2);
	// if found then remove
	if(p.first!=p.second){f.erase(p.first);}
	///// delete faces /////
	// for each face
	for(i=f.begin();i!=f.end();i++){
		fprintf(stderr,"deleteTriangleFaces: deleting Face %i\n",(*i)->index);
		// assign pointer
		q=(*i)->p;
		// adjust list pointer
		if (q==flh) {flh=q->next;}
		// remove face from linked list
		q=removeLink(q);
	}
}

void Triangle::findThirdVertices(void){
	std::vector<Vertex*> av;
	std::vector<Vertex*>::iterator i;
	// for both edge vertices, collect adjacent vertices
	av.insert(av.end(),e->v1->v.begin(),e->v1->v.end());
	av.insert(av.end(),e->v2->v.begin(),e->v2->v.end());
	// sort adjacent vertices
	sort(av.begin(),av.end());
	// for each pair of sequential vertices in sorted list
	for(i=av.begin();(i+1)!=av.end();i++){
		// if vertex pair is identical
		if ((*i)->index==(*(i+1))->index){
			// if pair are not part of associated vertices of edge
			if ((*i)->index!=e->v1->index && (*i)->index!=e->v2->index ){
				// if pair are not part of associated faces of edge
				if (*i!=e->va && *i!=e->vb){
					// add member of pair to list of third vertices
					t.push_back(*i);
					fprintf(stderr,"findThirdVertices: third vertex found %i\n",(*i)->index);
				}
			}
		}
	}
}

void collapseEdge(Edge *e,void_list *&flh,void_list *&vlh){
	Triangle tri(e);
	tri.findThirdVertices();
	if (!tri.t.empty()) {
		findTriangles(tri);
		tri.gatherTriangleVerts();
		tri.deleteTriangleVerts(vlh);
		tri.deleteTriangleFaces(flh);
	}
}

void simplifyMesh(Edge *e,void_list *&flh,void_list *&vlh,int &max_verts){
	fprintf(stderr,"removing edge v1 %i, v2 %i,length %.15g\n",e->v1->index,e->v2->index,sqrt(e->l));
	Vertex *v;
	///// delete redundant vertices and faces /////
	collapseEdge(e,flh,vlh);
	// create new vertex at edge midpoint
	v = addVertex(e,vlh,max_verts);
	///// fix all adjacent faces of two vertices //////
	fixFaces(v,e->v1);
	fixFaces(v,e->v2);
	///// delete faces of current edge//////
	deleteFaces(e,flh);
	///// delete current edge vertices //////
	deleteVertices(e,vlh);
}

// #####################################################
// #####################################################

void cleanup2(void_list *f,void_list *v){
	void_list *p,*q;
	// faces
	p=f;
	while (p!=NULL) {
		q=p->next;
		delete (Face*)p->data;
		delete p;
		p=q;
	}
	// vertices
	p=v;
	while (p!=NULL) {
		q=p->next;
		delete (Vertex*)p->data;
		delete p;
		p=q;
	}
}

void printMesh(void_list *vlh,void_list *flh){
	void_list *vend,*fend,*q;
	Vertex *v;
	Face *f;

	flh=addPrevious(flh);
	vlh=addPrevious(vlh);

    ////////// add pointers to ends of linked lists /////////
    // read out vertex linked list
    for (q=vlh;q!=NULL;q=q->next) { vend = q; if (q->next==NULL) break; }

    // read out face linked list
    for (q=flh;q!=NULL;q=q->next) { fend = q; if (q->next==NULL) break; }

    ////////// write data to stdout /////////
	// write out vertex linked list
	for (q=vend;q!=NULL;q=q->previous) {
		v = (Vertex*)q->data;
		printf("Vertex %i  %.15g %.15g %.15g\n",v->index,v->x,v->y,v->z);
	}
	// write out face linked list
	for (q=fend;q!=NULL;q=q->previous) {
		f=(Face*)q->data;
		printf("Face %i  %i %i %i\n",f->index,f->v1->index,f->v2->index,f->v3->index);
	}
}

void clearEdges(std::vector<Edge*> &e,hashtable_t &hme){
	// for each edge, delete pointer
	std::vector<Edge*>::iterator i;
	for(i=e.begin();i!=e.end();i++){delete (*i);}
	e.clear();
	// for each void_list, delete links
	void_list *p,*q;
	__gnu_cxx::hash_map<u4, void_list*,pair_hash,equ4,std::allocator<void_list*> >::iterator j;
//	ht_iterator j;
	for(j=hme.begin();j!=hme.end();j++){
	    p=(*j).second;
	    while (p!=NULL) {
	        q=p->next;
	        delete p;
	        p=q;
	    }
	}
	hme.clear();
}

void clearFaceEdges(void_list *flh){
	void_list *q;
    // read out face linked list
    for (q=flh;q!=NULL;q=q->next) {((Face*)q->data)->clearEdges();}
}
