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

	// for each EdgeBlock
	q=eb;
	while(q!=NULL && !found){
		// compute hashval of va,vb
		hashval = getEdgeHashVal(va,vb,q->ht->s);
//		printf("EdgeBlock: va %i,vb %i,hashval %i\n",va,vb,hashval);
		// for each edge in hashtable element pointed to by hashval
		p = eb->ht->t[hashval];
		while (p!=NULL && !found){
			// if v1,v2
			if  ((((Edge*)p->data)->v1==va||((Edge*)p->data)->v2==va)&& 
				(((Edge*)p->data)->v1==vb||((Edge*)p->data)->v2==vb)){found=true;}
			else	{p=p->next;}
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
//	printf("checkEdge: findEdge...");
	if(e=findEdge(eb,va,vb)){	
//		printf("checkEdge: updateEdge...");
		updateEdge(F,e,va,vb);
//		printf("complete.\n");
	}
	else {
//		printf("checkEdge: createEdge...");
		createEdge(F,eb,va,vb);
//		printf("complete.\n");
	}
}

void getEdges(void_list *flh,EdgeBlock *eb,int num,int print_flag){
	void_list *q;
	Face *f;
	int i=0;
	// for each face
	for (q=flh;q!=NULL;q=q->next) {
		f=(Face*)q->data;
		if (print_flag) {printf("face %i of %i\n",i++,num-1);}
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

void_list * addFace(void_list *as,Face *F){
	void_list *p;
	p=new void_list();
	p->next=as;
	p->data=(void*)F;
	as=p;
	return as;
}

void_list * getAssociatedFaces(Edge *e,void_list *flh){
	void_list *p,*as;
	as=NULL;
	Face* f;
	printf("Associated Faces and Vertices\n");
	for(p=flh;p!=NULL;p=p->next){
		f=(Face*)p->data;
		if( (f->v1==e->v1||f->v1==e->v2)||
			(f->v2==e->v1||f->v2==e->v2)||
			(f->v3==e->v1||f->v3==e->v2)
			){
			printf("Face %i %i %i %i\n",f->index,f->v1,f->v2,f->v3);
			as=addFace(as,f);
		}
	}
	return as;
}

int compare (const void* a, const void* b ) {
  return ( *(int*)a - *(int*)b );
}

void unique(int *array,int count){
	int i;
	for(i=0;i<count-1;i++){
		if(array[i]==array[i+1]){array[i]=0;}
	}
}

void getAssociatedVertices(void_list *af,void_list *vlh){

	// gather unique list of vertices in faces
	// count number of vertices in faces
	int count=0;
	void_list *p;
	Face *f;
	for(p=af;p!=NULL;p=p->next){count+=3;}
	// declare array
	int array[count];
	int i=0;
	for(p=af;p!=NULL;p=p->next){
		f=(Face*)p->data;
		array[i++]=f->v1;
		array[i++]=f->v2;
		array[i++]=f->v3;
	}
	// sort array
    qsort(array,count,sizeof(int),compare);
	unique(array,count);
	bool flag;
	Vertex *v;
	for(i=0;i<count;i++){
		if(array[i]){
			flag=false;
			p=vlh;
			while(p!=NULL && !flag){
				v=(Vertex*)p->data;
				if(v->index==array[i]){
					printf("Vertex %i %.15g %.15g %.15g\n",v->index,v->x,v->y,v->z);
					flag=true;
				}else{p=p->next;}
			}
		}
	}
}

void identifyBadFaces(EdgeBlock *eb,void_list *flh,void_list *vlh){
	EdgeBlock *q;
	Edge *e;
	void_list *p,*af;
	Face *f;
	int i,count=0;
	
	// for each EdgeBlock
	for(q=eb;q!=NULL;q=q->next){
		// for each edge in block
		for(i=0;i<q->n;i++){
			e=&q->e[i];
			// if edge traversal is not exactly once in each direction
			if ((e->c12+e->c21)==1) {
				printf("\nError: Edge is open.");
				printf(" va %i vb %i,",e->v1,e->v2);
				f = (Face*)(e->f)->data;
				printf(" Face %i %i %i %i\n",f->index,f->v1,f->v2,f->v3);
				count++;
			}
			else if (((e->c12==2)&&(e->c21==0))||((e->c12==0)&&(e->c21==2))){
				printf("\nError: Flip candidates.");
				printf(" va %i vb %i, c12 %i, c21 %i\n",e->v1,e->v2,e->c12,e->c21);
				for (p=e->f;p!=NULL;p=p->next) {
					f=(Face*)p->data;
					printf("     Face %i %i %i %i\n",f->index,f->v1,f->v2,f->v3);
					count++;
				}
			}
			else if ((e->c12&&(e->c21>1))||((e->c12>1)&&e->c21)) {
				printf("\nError: Topological Problem.");
				printf(" va %i vb %i, c12 %i, c21 %i\n",e->v1,e->v2,e->c12,e->c21);
				for (p=e->f;p!=NULL;p=p->next) {
					f=(Face*)p->data;
					printf("     Face %i %i %i %i\n",f->index,f->v1,f->v2,f->v3);
					count++;
				}
				af=getAssociatedFaces(e,flh);
				getAssociatedVertices(af,vlh);
			}
		}
	}
	printf("\n\n # BAD FACES = %i\n\n",count);
}

