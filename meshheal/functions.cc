void getData(char *infile,void_list *&flh,void_list *&vlh){
    char line[2048],*str,*eptr;
    FILE *F;
    void_list *vl;
    Vertex *v;
    void_list *fl;
    Face *f;

	// open first file
	F = fopen(infile,"r");
	if (!F) { printf("Couldn't open input file %s\n",infile);return;}

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

u4 computeHashValue(int key) {
    char cval[10];
    u4 result;
    sprintf(cval,"%i",key);
    result = hash((u1*)cval,(u4)strlen(cval),0xa205b064);
    return result;
}

int getEdgeHashVal(int va,int vb,int mask){
    u4 hashval1 = computeHashValue(va);
    u4 hashval2 = computeHashValue(vb);
    u4 result = hashval1+hashval2;
    return (int)(result&(mask-1));
}

// #####################################################
// #####################################################

Edge* findEdge(EdgeBlock *eb,int va,int vb){
    EdgeBlock *q;
    void_list *p;
    bool found=false;
    int hashval;
	Edge *e;

    // for each EdgeBlock
    q=eb;
    while(q!=NULL && !found){
        // compute hashval of va,vb
        hashval = getEdgeHashVal(va,vb,q->ht->s);
		if(hashval>q->ht->s-1){printf("hashval %i, size of hashtable %i\n",hashval,q->ht->s);}
        // for each edge in hashtable element pointed to by hashval
        p = eb->ht->t[hashval];
        while (p!=NULL && !found){
			e=(Edge*)p->data;
            // if v1,v2
            if  ((e->v1==va||e->v2==va)&&(e->v1==vb||e->v2==vb)){found=true;}
            else {p=p->next;}
        }if(!found){q=q->next;}
    }
    if(found){return e;}
    else{return NULL;}
}

void updateEdge(Face *F,Edge *e,int va,int vb){
    void_list *qq;
    //add face to edge
    qq = new void_list();
    qq->next = e->f;
    qq->data = (void*)F;
    e->f = qq;
    if  (e->v1==va && e->v2==vb) {e->c12+=1;}
    else {e->c21+=1;}
}


EdgeBlock* createEdge(Face *F,EdgeBlock *eb,int va,int vb){
    EdgeBlock *q;
    Edge *e;
    bool found=false;
    void_list *p;

    // for each EdgeBlock
    q=eb;
    while(q!=NULL && !found){
        // is there an open Edge slot
        if (q->p<(q->n-1)){
            // if yes grab Edge pointer
            e=&q->e[q->p];
            q->p++;
            found = true;
        } else {q=q->next;}
    }

    // if no open Edge slot found in any EdgeBlock
    if (!found) {
		printf("Increase number of edges in EdgeBlock!\n");exit(1);
/*
		printf("new edge block created!\n");
        // create a new EdgeBlock
        q = new EdgeBlock(10000);
        // store pointer to last EdgeBlock
        q->next=eb;
        // create HashTable
        q->ht = new HashTable(q->n);
        // grab Edge pointer
        e=&q->e[q->p];
        q->p++;
		eb=q;
*/
    }

    // load Edge
    ////////////////////
    e->v1 = va;
    e->v2 = vb;
    e->c12 = 1;
    e->c21 = 0;
    // add face to edge
    // set link data to point to Face class instance
    e->f = new void_list();
    (e->f)->next = NULL;
    (e->f)->data = (void*)F;
    // store edge pointer in hash table
    // compute hashval of va,vb
    int hashval = getEdgeHashVal(va,vb,q->ht->s);
    // load hash table
    p=new void_list();
    p->next = q->ht->t[hashval];
    p->data= (void*)e;
    q->ht->t[hashval] = p;
	return eb;
}


void checkEdge(Face *F,EdgeBlock *&eb,int va,int vb) {
    Edge *e;
    if(e=findEdge(eb,va,vb)){ updateEdge(F,e,va,vb);}
    else {eb=createEdge(F,eb,va,vb);}
}

void getEdges(void_list *flh,EdgeBlock *&eb,int num,int print_flag){
    void_list *q;
    Face *f;
    int i=0;
    // for each face
    for (q=flh;q!=NULL;q=q->next) {
        f=(Face*)q->data;
        checkEdge(f,eb,f->v1,f->v2);
        checkEdge(f,eb,f->v2,f->v3);
        checkEdge(f,eb,f->v3,f->v1);
    }
}

// #####################################################
// #####################################################

void initVertArray(void_list *vlh,Vertex **&vert_array,int max_vert){
	void_list *p;
	vert_array = new Vertex*[max_vert+1];
	for (p=vlh;p!=NULL;p=p->next){
		vert_array[((Vertex*)p->data)->index]=(Vertex*)p->data;
	}
if(0){
	for (p=vlh;p!=NULL;p=p->next){
		printf("initVertArray index %i\n",vert_array[((Vertex*)p->data)->index]->index);
	}
}
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

void_list* addLink(void_list* L, Vertex *v) {
	void_list *q;
	q = new void_list();
	q->next = L;
	q->data = (void*)v;
	return q;
}

double computeLongestEdge(EdgeBlock *eb,Vertex **vert_array,double epsilon) {
	EdgeBlock *q;
	Edge *e;
	Vertex *va,*vb;
	double d,max=0.0;

    // for each EdgeBlock
    for(q=eb;q!=NULL;q=q->next){
        // for each edge in block
        for(int i=0;i<q->n;i++){
            e=&q->e[i];
			// if edge is valid
			if(e->v1!=0 &&e->v2!=0){
				va=vert_array[e->v1];
				vb=vert_array[e->v2];
				// compute squared edge length
				d=(va->x/epsilon-vb->x/epsilon)*(va->x/epsilon-vb->x/epsilon)+(va->y/epsilon-vb->y/epsilon)*(va->y/epsilon-vb->y/epsilon)+(va->z/epsilon-vb->z/epsilon)*(va->z/epsilon-vb->z/epsilon);
				// if edge is longer than max, then save d
	            if (d>max){max=d;}
			}
		}
	}
	return max; // return threshold squared
}

void_list * gatherFreeVertices(EdgeBlock *eb,int max_vert,Vertex **vert_array,int print_flag) {
	///// gather free vertices /////
	EdgeBlock *q;
	Edge *e;
	void_list *p,*ptr,*vfree;
	vfree = NULL;
	Vertex *v;

    // for each EdgeBlock
    for(q=eb;q!=NULL;q=q->next){
        // for each edge in block
        for(int i=0;i<q->n;i++){
            e=&q->e[i];
			// if edge is free
            if ((e->c12+e->c21)==1) {
				p = new void_list();
				p->next = vfree;
				p->data = (void*)vert_array[e->v1];
				vfree = p;
				p = new void_list();
				p->next = vfree;
				p->data = (void*)vert_array[e->v2];
				vfree = p;
			}
		}
	}

	// go through linked list backwards and add previous pointers
	vfree=addPrevious(vfree);

	///// gather unique free vertices /////
	int v_array[max_vert+1];
    for (int i=0;i<max_vert+1;i++) {v_array[i]=0;}
	ptr=vfree;
    p=vfree;
    while(p!=NULL) {
		v=(Vertex*)p->data;
		if(!v_array[v->index]){v_array[v->index]=1;p=p->next;}
		else {
			// adjust list pointer
			if (p==ptr) { ptr=p->next;}
			// remove free face link
			p = removeLink(p);
		}
    }
	// adjust pointer
	vfree=ptr;

	if(0){
	    for (p=vfree;p!=NULL;p=p->next) {
			v=(Vertex*)p->data;
			printf("Free vertex %i %.15g %.15g %.15g\n",v->index,v->x,v->y,v->z);
	    }
	}

	return vfree;
}

// #####################################################
// #####################################################

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

void_list ** computeDistances(EdgeBlock *eb,void_list **dist_array,void_list *v,void_list *F,int& d_count,double epsilon,int print_flag,double threshold) {
	///// compute distances /////
	fprintf(stderr,"Computing vertex distances...");fflush(stderr);
	void_list *q,*qq,*dl,*dlh;
	double diffx,diffy,diffz,squared_dist;
	int i=0;
	Distance *d;
	dlh = NULL;
	d_count = 0;
	Vertex *vi,*vo;
	Edge *e;
	// for each free vertex
    for (q=v;q!=NULL;q=q->next) {
		vo=(Vertex*)q->data;
		// for each free vertex
	    for (qq=v;qq!=NULL;qq=qq->next) {
			vi=(Vertex*)qq->data;
			// if combination is unique
			if (vo->index > vi->index){
				// are vertices members of the same face?
    			if(e=findEdge(eb,vo->index,vi->index)){}
				else { // if not
					// do not factor out epsilon, leave as is
					diffx = vo->x/epsilon-vi->x/epsilon;
					diffy = vo->y/epsilon-vi->y/epsilon;
					diffz = vo->z/epsilon-vi->z/epsilon;
					squared_dist = diffx*diffx+diffy*diffy+diffz*diffz;
					// if distance less than longest edge, then save
					if (squared_dist<threshold){
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
    for (q=dlh;q!=NULL;q=q->next) {dist_array[i++]=q;}
	qsort(dist_array,d_count,sizeof(int),compare);

	if(0){ for (i=0;i<d_count;i++){ printf("Free vertices %i and %i, dist = %.15g\n",
				((Vertex*)((Distance*)(dist_array[i])->data)->vA)->index,
				((Vertex*)((Distance*)(dist_array[i])->data)->vB)->index,
				((Distance*)(dist_array[i])->data)->d); } }

	if (!print_flag) {printf("complete.\n");fflush(stdout);}
	return dist_array;
}

// #####################################################
// #####################################################

int countFreeVertices(void_list *v,int print_flag) {
	///// count number free vertices /////
	void_list *q;
	int i=0;
    for (q=v;q!=NULL;q=q->next) { i++; }
	fprintf(stderr,"num free vertices = %i\n",i);fflush(stderr);
	return i;
}

int findVerticesToMerge(void_list **array,int count,void_list *vfree,
						int print_flag,double epsilon) {
	int i=0;
	void_list *q;
	Vertex *v;
	Distance *d;
	bool flag=false,vflag1,vflag2;
	while (i<count && !flag) {
		d=(Distance*)array[i]->data;
		// is ith vertex pair free?
		vflag1=false;vflag2=false;
		q=vfree;
		// for each free vertex
		while(q!=NULL && (!vflag1||!vflag2)){
			v=(Vertex*)q->data;
			if (((Vertex*)d->vA)->index==v->index){vflag1=true;}
			if (((Vertex*)d->vB)->index==v->index){vflag2=true;}
			q=q->next;
		}
		if (vflag1&&vflag2){flag=true;}
		else{i++;}
	}
	if (i==count){
		fprintf(stderr,"ERROR! No eligible vertices to merge were found.\n");
	} else {
		fprintf(stderr,"vertices to merge: %i %i, dist = %.15g\n",
				((Vertex*)d->vA)->index,((Vertex*)d->vB)->index,sqrt(d->d*epsilon*epsilon));
	}
	return i;
}

void_list * fixVertices(void_list *v,Distance *d,int print_flag,void_list *&vbad) {
	///// remove second vertex in pair from vertex list /////
	void_list *ptr,*p;
	ptr=v;
	p=v;
	bool found = false;
	while (p!=NULL && !found) {
		if(((Vertex*)d->vB)->index == ((Vertex*)p->data)->index){
			// add pointer to bad vertex list
			vbad = addLink(vbad,(Vertex*)p->data);
			// adjust list pointer
			if (p==ptr) { ptr=p->next;}
			// remove link from list
			p=removeLink(p);
			found = true;
		}else{p=p->next;}
	}
	// adjust pointer
	v=ptr;
	return v;
}

void_list * fixFaces(void_list *flh,Distance *d,int print_flag) {
	///// replace all instances of second vertex in face list with first vertex /////
	void_list *p;
	Face *f;
	Vertex *va,*vb;
	va=(Vertex*)d->vA;
	vb=(Vertex*)d->vB;
	for (p=flh;p!=NULL;p=p->next) {
		f=(Face*)p->data;
		if(vb->index==f->v1){f->v1=va->index;}
		if(vb->index==f->v2){f->v2=va->index;}
		if(vb->index==f->v3){f->v3=va->index;}
	}

	// tag distance class for future vertex b deletion
	d->deleteme=true;

	return flh;
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

void_list * printVerticesFaces(void_list *vlh,void_list *flh){
	void_list *q,*vend,*fend;
	Vertex *v;
	Face *f;
	for (q=vlh;q!=NULL;q=q->next) {vend=q;}
	for (q=flh;q!=NULL;q=q->next) {fend=q;}
	// write out vertex linked list
	for (q=vend;q!=NULL;q=q->previous) {
		v=(Vertex*)q->data;		
      	printf("Vertex %i  %.15g %.15g %.15g\n",v->index,v->x,v->y,v->z);
	}
	// write out face linked list
	for (q=fend;q!=NULL;q=q->previous) {
		f=(Face*)q->data;
      	printf("Face %i  %i %i %i\n",f->index,f->v1,f->v2,f->v3);
	}
}

void deleteEdges(EdgeBlock *eb) {
    EdgeBlock *pp,*qq;
    // for each EdgeBlock
    pp=eb;
    while (pp!=NULL) {
        qq=pp->next;
        delete pp;
        pp=qq;
    }
}

void cleanup(EdgeBlock *eb,void_list *vf,void_list *v,
			void_list *f,void_list **d,int count,Vertex **vert_array,void_list *vbad){
	void_list *p,*q;
	deleteEdges(eb);
	deleteFreeVertices(vf);
	////////// delete vertices and faces /////////
	p=v;
	while (p!=NULL) {
		q=p->next;
		delete (Vertex*)p->data;
		delete p;
		p=q;
	}
	p=f;
	while (p!=NULL) {
		q=p->next;
		delete (Face*)p->data;
		delete p;
		p=q;
	}
	////////// delete dist_array /////////
	for(int i=0;i<count;i++){
		delete (Distance*)d[i]->data;
		delete d[i];
	}
	delete[] d;
	///////////////////////
	delete[] vert_array;
	////////// delete bad vertices /////////
	p=vbad;
	while (p!=NULL) {
		q=p->next;
		delete (Vertex*)p->data;
		delete p;
		p=q;
	}
}
