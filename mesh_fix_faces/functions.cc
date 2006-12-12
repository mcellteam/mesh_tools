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
	Edge *e;
	void_list *p;
	bool found=false;
	int hashval;

	// for each EdgeBlock
	q=eb;
	while(q!=NULL && !found){
		// compute hashval of va,vb
		hashval = getEdgeHashVal(va,vb,q->ht->s);
		// for each edge in hashtable element pointed to by hashval
		p = eb->ht->t[hashval];
		while (p!=NULL && !found){
			e=(Edge*)p->data;
			// if v1,v2
			if  ((e->v1==va||e->v2==va)&&(e->v1==vb||e->v2==vb)){found=true;}
			else {p=p->next;}
		}if(!found){q=q->next;}
	}
	if(found){return (Edge*)p->data;}
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

void createEdge(Face *F,EdgeBlock *eb,int va,int vb){
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
		// create a new EdgeBlock
		q = new EdgeBlock(10000);
		// store pointer to last EdgeBlock
		q->next=eb;
		// create HashTable
		q->ht = new HashTable(q->n);
		// grab Edge pointer
		e=&q->e[q->p];
		q->p++;
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
}

void checkEdge(Face *F,EdgeBlock *&eb,int va,int vb) {
	Edge *e;
	if(e=findEdge(eb,va,vb)){updateEdge(F,e,va,vb);}
	else {createEdge(F,eb,va,vb);}
}

void getEdges(void_list *flh,EdgeBlock *eb,int num){
	void_list *q;
	Face *f;
	int i=0;
	// for each face
	for (q=flh;q!=NULL;q=q->next) {
		f=(Face*)q->data;
//		fprintf(stderr,"face %i of %i\n",i++,num-1);
		checkEdge(f,eb,f->v1,f->v2);
		checkEdge(f,eb,f->v2,f->v3);
		checkEdge(f,eb,f->v3,f->v1);
	}
}

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

void cleanup(EdgeBlock *eb,void_list *f,void_list *v){
	void_list *p,*q;
	EdgeBlock *pp,*qq;
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
	// EdgeBlocks
	// for each EdgeBlock
	pp=eb;
	while (pp!=NULL) {
		qq=pp->next;
		delete pp;
		pp=qq;
	}
}

int evaluateEdge(Edge *e) {
    // evaluate the orientation of adjacent faces to given edge
    int a;
    void_list *p;
	Face *f;
    //reset traversal counts
    e->c12=e->c21=0;
    // for each associated face
     for (p=e->f;p!=NULL;p=p->next) {
		f=(Face*)p->data;
        // if face traverses edge ccw
        if  ((f->v1==e->v1)&&(f->v2==e->v2)||
            (f->v2==e->v1)&&(f->v3==e->v2) ||
            (f->v3==e->v1)&&(f->v1==e->v2))
            {e->c12++;}
        // if face traverses edge cw
        else {e->c21++;}
    }

    if ((e->c21==1)&&(e->c12==1)){a=0;}   // perfect
    else if (((e->c21==0)&&(e->c12==2))||
            ((e->c21==2)&&(e->c12==0))){a=1;} // needs simple flip
    else if ((e->c21+e->c12)==1){a=2;} // free
    else {a=-1;} // weird
    return a;
}

void flipFace(Face *f){
	printf("face was flipped\n");
    // flip the face
    int temp;
    temp = f->v3;
    f->v3 = f->v1;
    f->v1 = temp;
    f->flips++;
}

void switchFlipFace(Face *o,Face *n,void_list *flh,int val){
    // flip all faces with old group # and set group # to new one
    void_list *p;
	Face *f;
    // for every face in mesh file
    for (p=flh;p!=NULL;p=p->next) {
		f=(Face*)p->data;
        // if the face has group # of second face, then assign group # of first face
        if (f->g==o->g) {f->g=n->g;if(val){flipFace(f);}}
    }
}

void badEdge(Edge *e){
	Face *f;
	void_list *p;
	fprintf(stderr,"\nAbnormal edge skipped.\n");
	fprintf(stderr,"Edge traversal numbers are %i and %i.\n",e->c12,e->c21);
	// print out each associated face
	for (p=e->f;p!=NULL;p=p->next) {
		f=(Face*)p->data;
	    fprintf(stderr,"Face %i %i %i %i\n",f->index,f->v1,f->v2,f->v3);
	}
}

void getEdgeGroupsFaces(Edge *e,int &g1,int &g2,Face *&F1,Face *&F2){
	void_list *p;
	p=e->f;
	F1=(Face*)p->data;
	g1=F1->g;
	p=p->next;
	F2=(Face*)p->data;
	g2=F2->g;
}

void newGroup(Face *F1,Face *F2,int g,int val){
	F1->g=g;
	F2->g=g;
	if(val){flipFace(F2);}
}

void assignGroup(Face *F1,Face *F2,int g1,int g2,int val){
	if(g1&&!g2){F2->g=g1;if(val){flipFace(F2);}}
	else if(!g1&&g2){F1->g=g2;if(val){flipFace(F1);}}
}

void switchGroup(Face *F1,Face *F2,void_list *flh,int g1,int g2,int val){
	// if group # 1 is smaller, then switch and flip all faces with second group #
	if (g1<g2){switchFlipFace(F2,F1,flh,val);}
	else {switchFlipFace(F1,F2,flh,val);}
}

void groupFaces(EdgeBlock *eb,int print_flag,void_list *flh){
	EdgeBlock *q;
	Edge *e;
	Face *F1,*F2;
	int i,count=0,val,next_group=1,g1,g2;
	bool amalgamate=true,swap;
    while(amalgamate) {
        amalgamate = false;
		// for each EdgeBlock
		for(q=eb;q!=NULL;q=q->next){
			// for each edge in block
			for(i=0;i<q->n;i++){
				e=&q->e[i];
				if(e->hashval!=-1){
					val = evaluateEdge(e);
					if(val<0){badEdge(e);}
					else if(val<2){
						getEdgeGroupsFaces(e,g1,g2,F1,F2);
						if (!g1&&!g2) {newGroup(F1,F2,next_group++,val);amalgamate=true;}
		                // else if only first face has a group
		                else if ((g1&&!g2)||(!g1&&g2)){assignGroup(F1,F2,g1,g2,val);amalgamate=true;}
						// else if both faces have a group,
		                else if (g1&&g2) {
		                    // if group #s are different
		                    if (g1!=g2) {switchGroup(F1,F2,flh,g1,g2,val);amalgamate=true;}
	    	                else {
			                    // if group #s are same
		                        // if edge has two faces oriented differently
		                        if (val) {
		                            // problem!
		                            fprintf(stderr,"ERROR! Face %i and %i are both ",F1->index,F2->index);
		                            fprintf(stderr,"in group %ibut with different orientations.\n",F1->g);
		                            exit(1);
								} //else { // if oriented same // do nothing //}
		                    }
		                }
					}
				}
			}
		}
	}
}

int maxGroup(void_list *L){
    void_list *p;
    int max=0;
    for (p=L;p!=NULL;p=p->next) {
        if (max<((Face*)p->data)->g) {max=((Face*)p->data)->g;}
    }
    return max;
}

int compare (const void* a, const void* b ) {
  return ( *(int*)b - *(int*)a );
}

void countFaceGroups(void_list *flh,int print_flag){
	void_list *p;
	Face *f;
	int max_group,i;
	// find maxgroup#
	max_group = maxGroup(flh);
	// declare array of size maxgroup#+1
	int count[max_group+1];
	for(i=0;i<max_group+1;i++){
		count[i]=0;
	}
	// for each face
	for(p=flh;p!=NULL;p=p->next){
		count[((Face*)p->data)->g]++;
	}

	qsort(count,max_group+1,sizeof(int),compare);

	///// output groups /////
    fprintf(stderr,"\n***SUMMARY***\n\n");
    // for each group
    for (i=0;i<max_group+1;i++) {
		if(count[i]){
	        fprintf(stderr," group %i: %i faces\n\n",i,count[i]);
		}
    }
}

void printFlipActivity(void_list *flh){
	void_list *p;
	Face *f;
	// for each face
	for(p=flh;p!=NULL;p=p->next){
		f=(Face*)p->data;
		if(f->flips){ printf("Face %i, group %i was flipped %i times\n",f->index,f->g,f->flips);}
	}
}

void printMesh(void_list *vlh,void_list *flh){
	void_list *vend,*fend,*q;
	Vertex *v;
	Face *f;
    ////////// add pointers to ends of linked lists /////////
    // read out vertex linked list
    for (q=vlh;q!=NULL;q=q->next) {vend=q;}

    // read out face linked list
    for (q=flh;q!=NULL;q=q->next) {fend=q;}

	////////// write data to stdout /////////
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

